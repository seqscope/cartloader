import sys, os, gzip, logging, argparse, inspect, subprocess
import json
import math
import pandas as pd

from collections import Counter
from cartloader.utils.utils import create_custom_logger


SQRT_3 = math.sqrt(3)


def _hexagon_vertices(center_x, center_y, flat_width, orientation):
    """Return a closed ring describing a hexagon centred on (x, y).

    The generated coordinates stay in the same projected coordinate system as the
    input data (typically Web Mercator metres). Tippecanoe receives the CRS via
    ``-s EPSG:3857`` so it will reproject these planar coordinates as needed.
    """

    if flat_width <= 0:
        raise ValueError("flat_width must be positive")

    if orientation == 'pointy':
        # flat_width equals the left-to-right distance between opposite sides
        edge_length = flat_width / SQRT_3
        angle_offset = 30
    else:
        # flat_width equals the top-to-bottom distance between opposite sides
        edge_length = flat_width / 2.0
        angle_offset = 0
    coords = []
    for i in range(6):
        angle = math.radians(angle_offset + (i * 60))
        vertex_x = center_x + edge_length * math.cos(angle)
        vertex_y = center_y + edge_length * math.sin(angle)
        coords.append([vertex_x, vertex_y])

    coords.append(coords[0])
    return coords


def _build_parser(repo_dir, prog_name):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {prog_name}",
        description="Convert a generic TSV file into PMTiles format that can be converted to pmtiles"
    )

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-tsv', type=str, required=True, help='Generic TSV file to convert to PMTiles. Typically named as *.tsv.gz, and must include spatial coordinates.')
    inout_params.add_argument('--out-prefix', type=str, required=True, help='Prefix of output files. New directory will be created if needed')

    geom_params = parser.add_argument_group("Geometry Parameters", "Parameters controlling feature geometry.")
    geom_params.add_argument('--geometry-type', choices=['point', 'hexagon'], default='hexagon', help='Geometry to generate for each row. Use point for centroid output or hexagon for polygonal tiles (default: hexagon)')
    geom_params.add_argument('--hex-width', type=float, help='Flat-to-flat width of the hexagon in the input coordinate system (typically Web Mercator metres when used with tippecanoe) (used if --geometry-type is hexagon)')
    geom_params.add_argument('--hex-orientation', choices=['pointy', 'flat'], default='pointy', help='Orientation of generated hexagons in Web Mercator (pointy = vertex at north/south; flat = flat edge at top/bottom)')
    geom_params.add_argument('--lon-column', type=str, default='lon', help='Name of longitude column in the processed data (used for geometry generation)')
    geom_params.add_argument('--lat-column', type=str, default='lat', help='Name of latitude column in the processed data (used for geometry generation)')

    iocol_params = parser.add_argument_group("Input/Output Columns Parameters", "Input/output column parameters.")
    iocol_params.add_argument('--remove-column', type=str, nargs='+', help='List of column names to remove')
    iocol_params.add_argument('--rename-column', type=str, nargs='+', help='List of columns to rename in the format of [old_name1:new_name1] [old_name2:new_name2] .... Note that lon/lat must exist in the output columns')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters frequently used by users")
    aux_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    aux_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')
    aux_params.add_argument('--skip-pmtiles', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    aux_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')
    aux_params.add_argument('--in-delim', type=str, default='\t', help="Delimiter of the input file")
    aux_params.add_argument('--out-delim', type=str, default=',', help="Delimiter of the output file")
    aux_params.add_argument('--out-csv-suffix', type=str, default='.csv', help="Suffix of output CSV file")
    aux_params.add_argument('--out-geojson-suffix', type=str, default='.geojson', help="Suffix of output GeoJSON file when generating polygons")
    aux_params.add_argument('--out-pmtiles-suffix', type=str, default='.pmtiles', help="Suffix of output PMTiles file")
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level')
    aux_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level')
    aux_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles')
    aux_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles')
    aux_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory to be used (default: out-dir/tmp; specify /tmp if needed)')

    return parser

def _compute_output_columns(hdr_cols_input, args, logger):
    if args.remove_column:
        set_cols_to_remove = set(args.remove_column)
        if len(set_cols_to_remove) != len(args.remove_column):
            duplicates = list(filter(lambda x: Counter(args.remove_column)[x] > 1, set_cols_to_remove))
            logger.error(f"Duplicate column names {duplicates} in --remove-column")
            sys.exit(1)
    else:
        set_cols_to_remove = set()

    if args.rename_column:
        dict_cols_to_rename = dict([rename.split(':') for rename in args.rename_column])
        if len(dict_cols_to_rename) != len(args.rename_column):
            rename_keys = [rename.split(':')[0] for rename in args.rename_column]
            duplicates = list(filter(lambda x: Counter(rename_keys)[x] > 1, dict_cols_to_rename.keys()))
            logger.error(f"Duplicate column names {duplicates} in --rename-column")
            sys.exit(1)
    else:
        dict_cols_to_rename = {}

    if set_cols_to_remove & set(dict_cols_to_rename.keys()):
        overlapping_cols = set_cols_to_remove & set(dict_cols_to_rename.keys())
        logger.error(f"Overlapping column names {overlapping_cols} between --remove-column and --rename-column")
        sys.exit(1)

    hdr_cols_output = []
    hdr_index_output = []
    n_removed = 0
    n_renamed = 0
    for i, colname in enumerate(hdr_cols_input):
        if colname not in set_cols_to_remove:
            if colname in dict_cols_to_rename:
                hdr_cols_output.append(dict_cols_to_rename[colname])
                n_renamed += 1
            else:
                hdr_cols_output.append(colname)
            hdr_index_output.append(i)
        else:
            n_removed += 1

    if len(hdr_cols_output) != len(set(hdr_cols_output)):
        duplicates = list(filter(lambda x: Counter(hdr_cols_output)[x] > 1, hdr_cols_output))
        logger.error(f"Duplicate column names {duplicates} in the output columns")
        sys.exit(1)

    if args.remove_column and len(args.remove_column) != n_removed:
        missing_cols = set(set_cols_to_remove) - set(hdr_cols_input)
        logger.error(f"Error in removing columns. The following columns are missing {missing_cols}")
        sys.exit(1)

    if args.rename_column and len(args.rename_column) != n_renamed:
        missing_cols = set(dict_cols_to_rename.keys()) - set(hdr_cols_input)
        logger.error(f"Error in renaming columns. The following columns are missing {missing_cols}")
        sys.exit(1)

    return hdr_cols_output, hdr_index_output


def _validate_geometry_requirements(args, hdr_cols_output, logger):
    if args.geometry_type != 'hexagon':
        return

    if args.hex_width is None:
        logger.error("--hex-width must be provided when --geometry-type is hexagon")
        sys.exit(1)
    if args.hex_width <= 0:
        logger.error("--hex-width must be greater than 0")
        sys.exit(1)

    missing_geom_cols = [col for col in (args.lon_column, args.lat_column) if col not in hdr_cols_output]
    if missing_geom_cols:
        logger.error(f"Missing required coordinate columns {missing_geom_cols} in the processed header")
        sys.exit(1)

def _write_hexagon_features(chunk, geojson_handle, args, feature_state):
    skipped_geometry_rows = 0
    records = chunk.to_dict(orient='records')
    for record in records:
        lon_val = record.get(args.lon_column)
        lat_val = record.get(args.lat_column)

        if lon_val is None or lat_val is None or pd.isna(lon_val) or pd.isna(lat_val):
            skipped_geometry_rows += 1
            continue

        try:
            lon = float(lon_val)
            lat = float(lat_val)
        except (TypeError, ValueError):
            skipped_geometry_rows += 1
            continue

        properties = {}
        for key, value in record.items():
            if key == args.lon_column:
                properties[key] = lon
            elif key == args.lat_column:
                properties[key] = lat
            else:
                properties[key] = None if pd.isna(value) else value

        vertices = _hexagon_vertices(lon, lat, args.hex_width, args.hex_orientation)
        feature = {
            "type": "Feature",
            "properties": properties,
            "geometry": {"type": "Polygon", "coordinates": [vertices]}
        }

        if feature_state["first"]:
            geojson_handle.write('\n')
            feature_state["first"] = False
        else:
            geojson_handle.write(',\n')

        geojson_handle.write(json.dumps(feature, separators=(',', ':')))
        feature_state["count"] += 1

    return skipped_geometry_rows


def _process_input_chunks(args, hdr_cols_output, hdr_index_output, data_path, logger):
    logger.info(f"Converting the input data to {'GeoJSON' if args.geometry_type == 'hexagon' else 'CSV/TSV'} format")
    n_chunks = 0
    skipped_geometry_rows = 0
    geojson_handle = None
    feature_state = None

    try:
        if args.geometry_type == 'hexagon':
            geojson_handle = open(data_path, 'w', encoding='utf-8')
            geojson_handle.write('{"type":"FeatureCollection","features":[')
            feature_state = {"first": True, "count": 0}

        for chunk in pd.read_csv(args.in_tsv, sep=args.in_delim, chunksize=args.chunk_size, usecols=hdr_index_output):
            chunk.columns = hdr_cols_output

            if args.geometry_type == 'point':
                chunk.to_csv(
                    data_path,
                    sep=args.out_delim,
                    index=False,
                    header=(n_chunks == 0),
                    mode='w' if n_chunks == 0 else 'a'
                )
            else:
                skipped_geometry_rows += _write_hexagon_features(chunk, geojson_handle, args, feature_state)

            n_chunks += 1
            logger.info(f"Finished processing chunk {n_chunks} of size {args.chunk_size}...")
    finally:
        if geojson_handle is not None:
            if feature_state["first"]:
                geojson_handle.write(']}')
            else:
                geojson_handle.write('\n]}')
            geojson_handle.close()

    n_features = 0 if feature_state is None else feature_state["count"]
    return n_chunks, skipped_geometry_rows, n_features

def _run_tippecanoe(data_path, args, logger):
    logger.info(f"Converting the processed {'GeoJSON' if args.geometry_type == 'hexagon' else 'CSV/TSV'} data to PMTiles format")
    my_env = os.environ.copy()
    my_env["TIPPECANOE_MAX_THREADS"] = str(args.threads)
    cmd_parts = [
        f"'{args.tippecanoe}'",
        f"-t {args.tmp_dir}",
        f"-o {args.out_prefix}{args.out_pmtiles_suffix}",
        f"-Z {args.min_zoom}",
        f"-z {args.max_zoom}",
        "--force",
        "-s EPSG:3857",
        f"-M {args.max_tile_bytes}",
        f"-O {args.max_feature_counts}",
        "--drop-densest-as-needed",
        "--extend-zooms-if-still-dropping",
        f"'--preserve-point-density-threshold={args.preserve_point_density_thres}'",
    ]

    if args.geometry_type == 'hexagon':
        cmd_parts.append("--no-line-simplification")
        cmd_parts.append("--no-simplification-of-shared-nodes")
        cmd_parts.append("--generate-ids")
    else:
        cmd_parts.append("--no-duplication")
        cmd_parts.append("--no-clipping")
        cmd_parts.append("--buffer 0")
    cmd_parts.append(data_path)

    cmd = " ".join(cmd_parts)
    print("Command to run:" + cmd)
    result = subprocess.run(cmd, shell=True, env=my_env)
    if result.returncode != 0:
        logger.error("Error while converting processed data into PMTiles")
        sys.exit(1)

def render_hexagons(_args):
    """Convert a generic TSV file into PMTiles format."""
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    prog_name = inspect.getframeinfo(inspect.currentframe()).function
    parser = _build_parser(repo_dir, prog_name)
    args = parser.parse_args(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_render_hexagons" + args.log_suffix if args.log else None)

    # prepare dirs
    out_dir = os.path.dirname(args.out_prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(out_dir, "tmp")
    if not os.path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir, exist_ok=True)

    # metadata
    logger.info("Reading the metadata")
    with gzip.open(args.in_tsv, 'rt', encoding='utf-8') as rf:
        hdr_cols_input =  rf.readline().rstrip().split(args.in_delim)

    hdr_cols_output, hdr_index_output = _compute_output_columns(hdr_cols_input, args, logger)
    _validate_geometry_requirements(args, hdr_cols_output, logger)

    data_suffix = args.out_geojson_suffix if args.geometry_type == 'hexagon' else args.out_csv_suffix
    data_path = f"{args.out_prefix}{data_suffix}"

    _, skipped_geometry_rows, n_features = _process_input_chunks(args, hdr_cols_output, hdr_index_output, data_path, logger)
    if args.geometry_type == 'hexagon':
        if n_features == 0:
            logger.warning("Generated GeoJSON contains no polygon features. Check coordinate columns and width settings.")
        if skipped_geometry_rows:
            logger.warning(f"Skipped {skipped_geometry_rows} rows with missing or invalid coordinates while generating hexagons")

    if not args.skip_pmtiles:
        _run_tippecanoe(data_path, args, logger)

    # cleanup intermediate
    if not args.keep_intermediate_files and os.path.exists(data_path):
        logger.info(f"Removing the intermediate file {data_path}...")
        os.remove(data_path)

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

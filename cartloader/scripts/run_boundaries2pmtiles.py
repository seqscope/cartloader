import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, inspect, csv
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, execute_makefile, flexopen, unquote_str

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Split and convert transcripts TSV file into pmtiles")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-csv', type=str, help='Input CSV file to be containing cell boundaries')
#   inout_params.add_argument('--in-parquet', type=str, help='Input Parquet file to be containing cell boundaries')
    inout_params.add_argument('--in-clust', type=str, required=True, help='Input CSV/TSV file containing cluster information per cell')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='The output prefix. New directory will be created if needed')  
    inout_params.add_argument('--out-colname-clust', type=str, default='topK', help='Column name to be used for cluster information in the output pmtiles. Default: topK')
    inout_params.add_argument('--in-format', type=str, choices=['10x_xenium', 'cosmx_smi', 'vizgen_merscope', 'generic'], default='generic', help='Input format type. Default: generic')
    inout_params.add_argument('--colname-cell-id', type=str, default='cell_id', help='Column name for cell ID in the input CSV/Parquet file. Default: cell_id')
    inout_params.add_argument('--colname-x', type=str, default='vertex_x', help='Column name for vertex X coordinate in the input CSV/Parquet file. Default: vertex_x')
    inout_params.add_argument('--colname-y', type=str, default='vertex_y', help='Column name for vertex Y coordinate in the input CSV/Parquet file. Default: vertex_y')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    conv_params = parser.add_argument_group("Parameters for pmtiles conversion")
    conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level')
    conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level')
    conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles')
    conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles')
    conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory to be used (default: out-dir/tmp; specify /tmp if needed)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_boundaries2pmtiles(_args):
    """
    Convert cell boundary file to pmtiles
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_boundaries2pmtiles" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    # start mm
    mm = minimake()

    # create output directory if needed
    out_dir = os.path.dirname(args.out_prefix)
    out_base = os.path.basename(args.out_prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(out_dir, "tmp")
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir, exist_ok=True)

    ## read the cluster information
    logger.info(f"  * Reading cluster information from: {args.in_clust}")
    bcd2clust = {}
    with flexopen(args.in_clust, 'rt') as rf:
        sep = None
        for line in rf:
            if sep is None:
                if ',' in line:
                    sep = ','
                elif '\t' in line:
                    sep = '\t'
                elif ' ' in line:
                    sep = ' '
                else:
                    raise ValueError("Input cluster file must be CSV/TSV/SPACE delimited")
            toks = line.strip().split(sep)
            if len(toks) < 2:
                raise ValueError("Input cluster file must have at least two columns: bcd and cluster")
            bcd = toks[0]
            clust = toks[1]
            bcd2clust[bcd] = clust

    logger.info("  * Converting input CSV into GEOJSON format")
    if args.in_format == '10x_xenium':
        logger.info(f"Processing cell boundary information from {args.in_csv} in 10x_xenium format")

        ## read the cell CSV files
        bound_out=f"{args.out_prefix}.geojson"
        # col_id = "cell_id"
        # col_x = "vertex_x"
        # col_y = "vertex_y"
        col_id = args.colname_cell_id
        col_x = args.colname_x
        col_y = args.colname_y
        with flexopen(args.in_csv, "rt") as f, flexopen(bound_out, "wt") as wf:
            reader = csv.DictReader(f)
            hdrs = reader.fieldnames
            assert hdrs[0] == col_id and hdrs[1] == col_x and hdrs[2] == col_y
            current_cell_id = None
            current_vertices = []
            for row in reader:
                cell_id = unquote_str(row[col_id])
                x = float(row[col_x])
                y = float(row[col_y])
                if current_cell_id is None:
                    current_cell_id = cell_id
                if cell_id != current_cell_id:
                    clusteridx = bcd2clust.get(current_cell_id, "NA")
                    wf.write(f'{{"type": "Feature", "geometry": {{"type": "Polygon", "coordinates": [[{",".join(current_vertices)}]]}}, "properties": {{"cell_id": "{current_cell_id}", "topK": "{clusteridx}"}}}}\n')
                    current_cell_id = cell_id
                    current_vertices = []
                current_vertices.append(f'[{x},{y}]')
            if current_cell_id is not None:
                clusteridx = bcd2clust.get(current_cell_id, "NA")
                wf.write(f'{{"type": "Feature", "geometry": {{"type": "Polygon", "coordinates": [[{",".join(current_vertices)}]]}}, "properties": {{"cell_id": "{current_cell_id}", "topK": "{clusteridx}"}}}}\n')
    else:
        raise NotImplementedError(f"Input format {args.in_format} is not implemented yet")

    ## convert CSV into PMtiles
    tippecanoe_cmd = " ".join([
        f"TIPPECANOE_MAX_THREADS={args.threads}",
        f"'{args.tippecanoe}'",
        f"-t {args.tmp_dir}",
        f"-o {args.out_prefix}.pmtiles",
        f"-Z {args.min_zoom} -z {args.max_zoom} --force",
        f"-s EPSG:3857 -M {args.max_tile_bytes} -O {args.max_feature_counts}",
        f"--drop-densest-as-needed",
        f"--extend-zooms-if-still-dropping",
        f"'--preserve-point-density-threshold={args.preserve_point_density_thres}'",
        f"--no-clipping",
        f"--buffer 0",
        f"{args.out_prefix}.geojson"
    ])

    logger.info(f"  * Running tippecanoe command: {tippecanoe_cmd}")
    result = subprocess.run(tippecanoe_cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        logger.error(f"Command {tippecanoe_cmd}\nfailed with error: {result.stderr.decode()}")
        sys.exit(1)
    else:
        logger.info("   * PMTiles creation command completed successfully")
    
    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

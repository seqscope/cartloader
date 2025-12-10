import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, csv, yaml, inspect, json
import math
import pandas as pd
import numpy as np
from scipy.io import mmread
import subprocess

from shapely.geometry import shape, mapping
from shapely.affinity import scale as shapely_scale

from cartloader.utils.utils import create_custom_logger, flexopen, unquote_str, smartsort, write_dict_to_file, load_file_to_dict, scheck_app
from cartloader.scripts.import_xenium_cell import process_cluster_csv, read_de_csv, write_de_tsv, write_cmap_tsv, tile_csv_into_pmtiles, make_factor_dict, write_umap_tsv, umap_tsv2pmtiles, umap_tsv2png, umap_tsv2indpng
from cartloader.scripts.sge_convert import extract_unit2px_from_json

def _rescale_geometry(geom, units_per_um):
    """Rescale geometry coordinates into microns when needed."""
    if units_per_um is None or math.isclose(units_per_um, 1.0):
        return geom
    if units_per_um == 0:
        raise ValueError("units_per_um must be non-zero")
    scale_factor = 1.0 / units_per_um
    return shapely_scale(geom, xfact=scale_factor, yfact=scale_factor, origin=(0, 0))

def process_cell_geojson_w_mtx(cells_geojson, cell_ftr_mex, cells_out, bcd2clusteridx, units_per_um):
    # Load barcodes
    bcd_path = os.path.join(cell_ftr_mex, "barcodes.tsv.gz")
    df_bcd_idx = pd.read_csv(bcd_path, header=None, names=["cell_id"])
    df_bcd_idx["barcode_index"] = df_bcd_idx.index + 1  # 1-based indexing (MatrixMarket style)

    # Map barcode to cluster idx
    df_bcd_idx["clusteridx"] = df_bcd_idx["cell_id"].map(bcd2clusteridx)

    # Load mtx
    mtx_path = os.path.join(cell_ftr_mex, "matrix.mtx.gz")
    with gzip.open(mtx_path, "rt") as f:
        matrix = mmread(f).tocsr()

    # Compute UMI counts per barcode (i.e., column-wise sum)
    barcode_sums = matrix.sum(axis=0).A1
    df_bcd_ct = pd.DataFrame({
        "barcode_index": range(1, len(barcode_sums) + 1),
        "count": barcode_sums
    })

    # Load convert cells_geojson into csv
    with open(cells_geojson, "r") as f:
        cell_data = json.load(f)
    
    def _iter_cell_centroids():
        for feature in cell_data["features"]:
            formatted_id = f"cellid_{feature['properties']['cell_id']:09d}-1"
            geom = shape(feature["geometry"])
            scaled_geom = _rescale_geometry(geom, units_per_um)
            centroid = scaled_geom.centroid
            yield {
                "cell_id": formatted_id,
                "lon": centroid.x,
                "lat": centroid.y,
            }

    df_bcd_xy = pd.DataFrame(_iter_cell_centroids())

    # Merge and select final columns
    df_bcd = (
        df_bcd_ct
        .merge(df_bcd_idx, on="barcode_index", how="left")
        .merge(df_bcd_xy, on="cell_id", how="left")
        [["lon", "lat", "cell_id", "count", "clusteridx"]]
        .rename(columns={"clusteridx": "topK"})
    )
    df_bcd["topK"] = df_bcd["topK"].apply(lambda x: str(x) if pd.notna(x) else "NA")
    df_bcd.to_csv(cells_out)

def process_boundaries_geojson(input_geojson, output_geojson, bcd2clusteridx, units_per_um):
    with open(input_geojson) as f:
        geojson = json.load(f)

    for feature in geojson["features"]:
        raw_id = feature["properties"]["cell_id"]
        formatted_id = f"cellid_{raw_id:09d}-1"
        clusteridx = bcd2clusteridx.get(formatted_id, "NA")
        geom = _rescale_geometry(shape(feature["geometry"]), units_per_um)
        feature["properties"] = {
            "cell_id": formatted_id,
            "topK": str(clusteridx)
        }
        feature["geometry"] = mapping(geom)

    with open(output_geojson, "w") as f:
        json.dump({"type": "FeatureCollection", "features": geojson["features"]}, f, indent=2)

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Import cell segmentation results from Visium HD Space Ranger outputs and package cells/boundaries into PMTiles (optional catalog update)"
    )
    run_params = parser.add_argument_group("Run Options", "Execution controls for generating and running the Makefile")
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe) (default: 4)')
    run_params.add_argument('--log', action='store_true', default=False, help='Write logs to a file alongside outputs')
    run_params.add_argument('--log-suffix', type=str, default=".log", help='Log filename suffix; final path is <outprefix>_importxenium<suffix> (default: .log)')

    cmd_params = parser.add_argument_group("Commands", "Select one or more actions to perform")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Enable all actions: --cells, --boundaries')
    cmd_params.add_argument('--cells', action='store_true', default=False, help='Import segmented cells and generate PMTiles')
    cmd_params.add_argument('--boundaries', action='store_true', default=False, help='Import segmented cell boundaries and generate GeoJSON')
    cmd_params.add_argument('--umap', action='store_true', default=False, help='Import UMAP projection and generate PMTiles for UMAP points')
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update an existing catalog.yaml generated by run_cartload2 (not included in --all)')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'Two input modes: 1) JSON — set --in-json with keys (CELL_FEATURE_MEX, CELL_GEOJSON, CLUSTER, DE); 2) Manual — set --in-dir and provide locations, see "Manual Input Parameters"')
    inout_params.add_argument('--in-json', type=str, help='Path to input JSON with paths for cells, boundaries, clusters, and differential expression.')
    inout_params.add_argument('--outprefix', type=str, required=True, help='Output prefix')
    inout_params.add_argument('--id', type=str, help='Identifier for the cell factor; if omitted, uses basename of --outprefix')
    inout_params.add_argument('--name', type=str, help='Display name for the cell factor; if omitted, uses basename of --outprefix')

    aux_inout_params = parser.add_argument_group("Manual Input Files Parameters", "Manually specify input directory (--in-dir) and input file locations under --in-dir")
    aux_inout_params.add_argument('--in-dir', type=str, help='Input directory containing the Space Ranger output files.')
    aux_inout_params.add_argument('--geojson-cells', type=str, default="segmented_outputs/cell_segmentations.geojson", help='Location of GEOJSON with cell locations under --in-dir (default: segmented_outputs/cell_segmentations.geojson)')
    aux_inout_params.add_argument('--mtx-cells', type=str, default="segmented_outputs/filtered_feature_cell_matrix", help='Directory location of the cell feature MatrixMarket files under --in-dir (default: segmented_outputs/filtered_feature_cell_matrix)')
    aux_inout_params.add_argument('--csv-clust', type=str, default="analysis/clustering/gene_expression_graphclust/clusters.csv", help='Location of CSV with cell cluster assignments under --in-dir (default: analysis/clustering/gene_expression_graphclust/clusters.csv)')
    aux_inout_params.add_argument('--csv-diffexp', type=str, default="analysis/diffexp/gene_expression_graphclust/differential_expression.csv", help='Location of CSV with differential expression results under --in-dir (default: analysis/diffexp/gene_expression_graphclust/differential_expression.csv)')
    aux_inout_params.add_argument('--csv-umap', type=str, default="analysis/pca/gene_expression_10_components/projection.csv", help='Location of CSV with UMAP results under --in-dir (default: analysis/pca/gene_expression_10_components/projection.csv')

    # - scaling
    aux_inout_params.add_argument('--scale-json', type=str, default=None, help=f'Location of scale JSON under --in-dir. If set, defaults --units-per-um from microns_per_pixel in this JSON file (No default value applied. Typical locations: square_002um/spatial/scalefactors_json.json; binned_outputs/square_002um/spatial/scalefactors_json.json")')
    aux_inout_params.add_argument('--units-per-um', type=float, default=1, help='Coordinate units per µm in inputs (default: 1).')

    aux_colnames_params = parser.add_argument_group("Auxiliary Colname Parameters", "Override column names for input files")
    aux_colnames_params.add_argument('--clust-colname-barcode', type=str, default='Barcode', help='Column name for cell barcode in --csv-clust (default: Barcode)')
    aux_colnames_params.add_argument('--clust-colname-cluster', type=str, default='Cluster', help='Column name for cluster label in --csv-clust (default: Cluster)')
    aux_colnames_params.add_argument('--umap-colname-barcode', type=str, default='Barcode', help='Column name for cell barcode in --csv-umap (default: Barcode)')
    aux_colnames_params.add_argument('--umap-colname-x', type=str, default='UMAP-1', help='Column name for UMAP X coordinate in --csv-umap (default: UMAP-1)')
    aux_colnames_params.add_argument('--umap-colname-y', type=str, default='UMAP-2', help='Column name for UMAP Y coordinate in --csv-umap (default: UMAP-2)')

    aux_conv_params = parser.add_argument_group("Auxiliary PMTiles Conversion Parameters")
    aux_conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level (default: 10)')
    aux_conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level (default: 18)')
    aux_conv_params.add_argument('--umap-min-zoom', type=int, default=0, help='Minimum zoom level for UMAP (default: 0)')
    aux_conv_params.add_argument('--umap-max-zoom', type=int, default=18, help='Maximum zoom level for UMAP (default: 18)')
    aux_conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles (default: 5000000)')
    aux_conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles (default: 500000)')
    aux_conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles (default: 1024)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Advanced settings; defaults work for most cases")    
    aux_params.add_argument('--tsv-cmap', type=str, default=f"{repo_dir}/assets/fixed_color_map_60.tsv", help=f'Location of TSV with color mappings for clusters (default: {repo_dir}/assets/fixed_color_map_60.tsv)')
    aux_params.add_argument('--de-max-pval', type=float, default=0.01, help='Maximum p-value for differential expression (default: 0.01)')
    aux_params.add_argument('--de-min-fc', type=float, default=1.2, help='Minimum fold change for differential expression (default: 1.2)')
    aux_params.add_argument('--catalog-yaml', type=str, help='Path to catalog.yaml to update (used with --update-catalog; default: <in-dir>/catalog.yaml)')
    aux_params.add_argument('--out-catalog-yaml', type=str, help='Path to save the updated catalog.yaml as a new file instead of overwriting the input (--catalog-yaml). Defaults to the same path as --catalog-yaml (used with --update-catalog)')
    #aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory for intermediate files (default: out-dir/tmp or /tmp if specified)')

    env_params = parser.add_argument_group("Env Parameters", "Tool paths (override defaults if needed)")
    env_params.add_argument('--R', type=str, default="R", help='Path to R binary (default: R).')
    env_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # Sanity check: JSON mode and manual mode are mutually exclusive
    if args.in_json and args.in_dir:
        parser.error("Cannot enable both JSON mode (--in-json) and manual input mode (--in-dir and/or --geojson-cells/--mtx-cells/--csv-clust/--csv-diffexp). Choose one.")

    return args

def import_visiumhd_cell(_args):
    """
    Import cell segmentation results from Space Ranger output
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.outprefix + "_import_visium_cell" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    if args.all:
        args.cells = True
        args.boundaries = True
        args.umap = True
    
    if not args.cells and not args.boundaries and not args.update_catalog and not args.umap:
        raise ValueError("At least one action should be activated.")

    # create output directory if needed
    out_dir = os.path.dirname(args.outprefix)
    out_base = os.path.basename(args.outprefix)
    if not os.path.exists(out_dir) and out_dir != "":
        os.makedirs(out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(out_dir, "tmp")
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir, exist_ok=True)

    #temp_fs=[]

    # read in_json if provided 
    scale_json = None

    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"File not found: {args.in_json} (--in-json)"
        raw_data = load_file_to_dict(args.in_json)
        cell_data = raw_data.get("CELLS", raw_data)  # support flat dicts from older scripts
        scale_json = raw_data.get("SGE", {}).get("SCALE", None)
    else:
        if args.in_dir is None:
            raise ValueError("--in-dir is required when --in-json is not provided")
        scale_json =  os.path.join(args.in_dir, args.scale_json) if args.scale_json else None
        cell_data = {
            "CELL_FEATURE_MEX": f"{args.in_dir}/{args.mtx_cells}",
            "CELL_GEOJSON": f"{args.in_dir}/{args.geojson_cells}",
            "CLUSTER": f"{args.in_dir}/{args.csv_clust}",
            "DE": f"{args.in_dir}/{args.csv_diffexp}",
            "UMAP_PROJ": f"{args.in_dir}/{args.csv_umap}"
        }

    # set units_per_um if scale_json is provided
    if scale_json is not None:
        assert os.path.exists(scale_json), f"File not found: {scale_json} (--scale-json)"
        args.units_per_um = extract_unit2px_from_json(scale_json)
        logger.info(f"Setting --units-per-um = {args.units_per_um} from {scale_json}")
    else:
        logger.warning(f"No scale JSON provided; assuming --units-per-um = {args.units_per_um}")
    
    # Cluster/DE
    if args.cells or args.boundaries or args.umap:
        clust_in = cell_data.get("CLUSTER", None)
        # presence and existence for cluster file
        assert clust_in is not None, ('Path not provided: "CLUSTER" in --in-json' if args.in_json is not None else 'Path not provided: --csv-clust')
        assert os.path.exists(clust_in), (f'File not found: {clust_in} ("CLUSTER" in --in-json)' if args.in_json is not None else f'File not found: {clust_in} (--csv-clust)')

        logger.info(f"Loading cell cluster data from {clust_in}")
        sorted_clusters, cluster2idx, bcd2clusteridx = process_cluster_csv(
            clust_in,
            barcode_col=args.clust_colname_barcode,
            cluster_col=args.clust_colname_cluster,
        )
        logger.info(f"  * Loaded {len(bcd2clusteridx)} cells")

        ## write the color map
        cmap_out=f"{args.outprefix}-rgb.tsv"
        assert os.path.exists(args.tsv_cmap), f"File not found: {args.tsv_cmap} (--tsv-cmap)"        
        
        logger.info(f"  * Writing color map from {args.tsv_cmap} to {cmap_out}")
        write_cmap_tsv(cmap_out, args.tsv_cmap, sorted_clusters)

    if args.cells or args.boundaries:
        ## read/write DE results
        de_in = cell_data.get("DE", None)
        assert de_in is not None, ('Path not provided: "DE" in --in-json' if args.in_json is not None else 'Path not provided: --csv-diffexp')
        assert os.path.exists(de_in), (f'File not found: {de_in} ("DE" in --in-json)' if args.in_json is not None else f'File not found: {de_in} (--csv-diffexp)')        

        logger.info(f"  * Reading DE results from {de_in}")
        clust2genes=read_de_csv(de_in, cluster2idx, args.de_min_fc, args.de_max_pval)

        ## write DE results
        de_out=f"{args.outprefix}-cells-bulk-de.tsv"
        logger.info(f"  * Writing DE results for {len(clust2genes)} clusters) to {de_out}")
        write_de_tsv(clust2genes, de_out, sorted_clusters)

    # Process segmented calls 
    if args.cells or args.boundaries:
        cells_json=cell_data.get("CELL_GEOJSON", None)
        assert cells_json is not None, (f'Path not provided: CELL_GEOJSON in --in-json' if args.in_json is not None else f'Path not provided: --geojson-cells')
        assert os.path.exists(cells_json), (f'File not found: {cells_json} (CELL_GEOJSON in --in-json)' if args.in_json is not None else f'File not found: {cells_json} (--geojson-cells)')
    
    if args.cells:
        cell_ftr_mex = cell_data.get("CELL_FEATURE_MEX", None)
        assert cell_ftr_mex is not None, ('Path not provided: CELL_FEATURE_MEX in --in-json' if args.in_json is not None else 'Path not provided: --mtx-cells')
        assert os.path.isdir(cell_ftr_mex), (f'Directory not found: {cell_ftr_mex} (CELL_FEATURE_MEX in --in-json)' if args.in_json is not None else f'Directory not found: {cell_ftr_mex} (--mtx-cells)')
        # required files exist
        bcd_path = os.path.join(cell_ftr_mex, "barcodes.tsv.gz")
        mtx_path = os.path.join(cell_ftr_mex, "matrix.mtx.gz")
        assert os.path.exists(bcd_path), (f"File not found: {bcd_path} (barcodes.tsv.gz in CELL_FEATURE_MEX from --in-json)" if args.in_json is not None else f"File not found: {bcd_path} (barcodes.tsv.gz from --mtx-cells)")
        assert os.path.exists(mtx_path), (f"File not found: {mtx_path} (matrix.mtx.gz in CELL_FEATURE_MEX from --in-json)" if args.in_json is not None else f"File not found: {mtx_path} (matrix.mtx.gz from --mtx-cells)")

        logger.info(f"Processing cell information from {cells_json}; {cell_ftr_mex}")

        # convert cell geojson into csv format with cell_id, x, y, count
        cells_out=f"{args.outprefix}-cells.csv"
        logger.info(f"  * Reading cell data from {cells_json}; {cell_ftr_mex} and extracting geometry to {cells_out}")
        process_cell_geojson_w_mtx(cells_json, cell_ftr_mex, cells_out, bcd2clusteridx, args.units_per_um)

        cells_pmtiles = f"{args.outprefix}-cells.pmtiles"
        logger.info(f"  * Generating PMTiles from cell geometry data into cell pmtiles: {cells_pmtiles}")
        tile_csv_into_pmtiles(cells_out, cells_pmtiles, args, logger, no_dup=True)
        #temp_fs.append(cells_out)

    # Process cell boundaries
    if args.boundaries:
        logger.info(f"Processing cell boundary information from {cells_json}")
        ## read the cell CSV files
        bound_out=f"{args.outprefix}-boundaries.geojson"
        logger.info(f"  * Reading cell boundary data from {cells_json} and extracting geometry to {bound_out}")
        process_boundaries_geojson(cells_json, bound_out, bcd2clusteridx, args.units_per_um)

        ## Run the tippecanoe command
        bound_pmtiles = f"{args.outprefix}-boundaries.pmtiles"
        logger.info(f"  * Generating PMTiles from boundary geometry data into boundary pmtiles: {bound_pmtiles}")
        tile_csv_into_pmtiles(bound_out, bound_pmtiles, args, logger, no_dup=False)
        #temp_fs.append(bound_out)

    # UMAP
    if args.umap:
        scheck_app(args.R)

        umap_in = cell_data.get("UMAP_PROJ", None)
        umap_tsv_out = f"{args.outprefix}-umap.tsv.gz"
        umap_pmtiles = f"{args.outprefix}-umap.pmtiles"

        assert umap_in is not None, ('Path not provided: "UMAP_PROJ" in --in-json' if args.in_json is not None else 'Path not provided: --csv-umap')
        assert os.path.exists(umap_in), (f'File not found: {umap_in} ("UMAP_PROJ" in --in-json)' if args.in_json is not None else f'File not found: {umap_in} (--csv-umap)')

        logger.info(f"Processing UMAP projection from {umap_in}")
        write_umap_tsv(umap_in, umap_tsv_out, bcd2clusteridx, args)

        logger.info(f"  * Generated PMTiles for UMAP projection:{umap_pmtiles}")
        umap_tsv2pmtiles(umap_tsv_out, umap_pmtiles, args)

        logger.info(f"  * UMAP Visualization for all factors...")
        umap_tsv2png(umap_tsv_out, args.outprefix, cmap_out)

        logger.info(f"  * UMAP Visualization (plot for individual factors; colorized by cluster)...")
        umap_tsv2indpng(umap_tsv_out, args.outprefix, cmap_out)

    # JSON/YAML (always summary)
    factor_id = out_base if args.id is None else args.id
    factor_name = out_base if args.name is None else args.name
    pmtiles_keys=[]
    if os.path.exists(f"{args.outprefix}-cells.pmtiles"):
        pmtiles_keys.append("cells")
    if os.path.exists(f"{args.outprefix}-boundaries.pmtiles"):
        pmtiles_keys.append("boundaries")

    if os.path.exists(f"{args.outprefix}-umap.pmtiles") and os.path.exists(f"{args.outprefix}-umap.tsv.gz") and os.path.exists(f"{args.outprefix}.umap.png") and os.path.exists(f"{args.outprefix}.umap.single.binary.png"):
        umap_src = True
    else:
        umap_src = False

    out_assets_f=f"{args.outprefix}_assets.json"
    logger.info(f"Summarizing assets information into {out_assets_f}")
    new_factor = make_factor_dict(factor_id, factor_name, args.outprefix, factor_type="cell", pmtiles_keys=pmtiles_keys, umap_src=umap_src)
    write_dict_to_file(new_factor, out_assets_f, check_equal=True)

    if args.update_catalog:
        logger.info(f"Updating catalog YAML file: {args.catalog_yaml}")

        ## read the input catalog.yaml file
        if args.catalog_yaml is None:
            out_dir = os.path.dirname(args.outprefix)
            args.catalog_yaml = os.path.join(out_dir, "catalog.yaml")
        
        with open(args.catalog_yaml, 'r') as f:
            catalog = yaml.load(f, Loader=yaml.FullLoader)  # Preserves order

        ## add files to the catalog
        new_factor = make_factor_dict(factor_id, factor_name, out_base, factor_type="cell", pmtiles_keys=pmtiles_keys, umap_src=umap_src)
        print(new_factor)
        if "factors" not in catalog["assets"]:
            catalog["assets"]["factors"] = [new_factor]
        else:
            factors_list = catalog["assets"]["factors"]
            factors_list = [item for item in factors_list if item.get('id') != factor_id]
            factors_list.append(new_factor)
            catalog["assets"]["factors"] = factors_list

        ## write the updated catalog.yaml file
        out_yaml = args.out_catalog_yaml if args.out_catalog_yaml is not None else args.catalog_yaml
        with open(out_yaml, 'w') as f:
            yaml.dump(catalog, f, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)
        logger.info(f"Successfully wrote the catalog.yaml file: {out_yaml}")

    ## clean the temp files
    # if not args.keep_intermediate_files:
    #     logger.info(f"Cleaning intermediate files")
    #     if len(temp_fs) >0:
    #         for temp_f in temp_fs:
    #             if os.path.exists(temp_f):
    #                 os.remove(temp_f)


    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

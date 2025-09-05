import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, csv, yaml, inspect
import pandas as pd
import numpy as np
from scipy.stats import chi2

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, flexopen, unquote_str, smartsort, write_dict_to_file, load_file_to_dict


def process_cluster_csv(clust_csv):
    bcd2cluster = {}
    cluster2cnt = {}
    with flexopen(clust_csv, "rt") as f:
        for line in f:
            if line.startswith("Barcode"):
                continue
            (bcd, cluster) = line.strip().split(",")
            bcd = unquote_str(bcd)
            bcd2cluster[bcd] = cluster
            cluster2cnt[cluster] = cluster2cnt.get(cluster, 0) + 1
    #sorted_clusters, cluster2idx = smartsort(clusters)
    ## sort clusters from largest to smallest
    sorted_clusters = sorted(cluster2cnt.keys(), key=lambda x: cluster2cnt[x], reverse=True)
    cluster2idx = {cluster: idx for idx, cluster in enumerate(sorted_clusters)}
    bcd2clusteridx = {bcd: cluster2idx[clust] for bcd, clust in bcd2cluster.items() if clust in cluster2idx}
    return sorted_clusters, cluster2idx, bcd2cluster, bcd2clusteridx

def process_cells_csv(cells_csv, out_csv, bcd2clusteridx):
    with flexopen(cells_csv, "rt") as f, flexopen(out_csv, "wt") as wf:
        reader = csv.DictReader(f)
        required_cols = {"cell_id", "x_centroid", "y_centroid", "transcript_counts"}
        if not required_cols.issubset(reader.fieldnames):
            raise ValueError(f"Missing expected columns: {required_cols - set(reader.fieldnames)}")
        wf.write("lon,lat,cell_id,count,topK\n")
        for row in reader:
            cell_id = row["cell_id"]
            x = row["x_centroid"]
            y = row["y_centroid"]
            count = row["transcript_counts"]
            clusteridx = bcd2clusteridx.get(cell_id, "NA")
            wf.write(",".join([x, y, cell_id, count, str(clusteridx)]) + "\n")

def read_de_csv(de_path, cluster2idx, fc_threshold, pval_threshold):
    clust2genes = {}
    with flexopen(de_path, "rt") as f:
        hdrs = f.readline().strip().split(",")
        assert hdrs[1] == "Feature Name"
        for line in f:
            toks = line.strip().split(",")
            gene_name = toks[1]
            for i in range(2, len(toks), 3):
                cluster = str((i - 2) // 3 + 1)
                clusteridx = cluster2idx.get(cluster, "NA")
                if clusteridx == "NA":
                    continue
                pval = float(toks[i+2]) 
                fc = 2 ** float(toks[i+1])
                avgcount = float(toks[i])
                if fc > fc_threshold and pval < pval_threshold:
                    clust2genes.setdefault(clusteridx, []).append([gene_name, avgcount, pval, fc])
    return clust2genes

def write_de_tsv(clust2genes, out_path, sorted_clusters):
    with flexopen(out_path, "wt") as wf:
        wf.write("\t".join(["gene","factor","Chi2","pval","FoldChange","gene_total","log10pval"]) + "\n")
        for clusteridx in range(len(sorted_clusters)):
            if clusteridx not in clust2genes:
                continue
            for gname, avgcount, pval, fc in sorted(clust2genes[clusteridx], key=lambda x: x[3], reverse=True):
                if pval == 0:
                    log10pval, chisq = 320, 1465.911
                else:
                    log10pval = -np.log10(pval)
                    chisq = chi2.isf(pval, 1)
                wf.write("\t".join([gname, str(clusteridx), f"{chisq:.4f}", f"{pval:.4e}",f"{fc:.4f}", f"{avgcount:.4f}", f"{log10pval:.4f}"]) + "\n")

def write_cmap_tsv(out_cmap, tsv_cmap, sorted_clusters):
    with flexopen(out_cmap, "wt") as wf:
        with flexopen(tsv_cmap, "rt") as f:
            # Write header line
            wf.write(f.readline())

            # Write color lines for each cluster
            for i in range(len(sorted_clusters)):
                line = f.readline()
                if not line:
                    raise ValueError(
                        f"Not enough colors in the color map file {tsv_cmap}"
                    )
                wf.write(line)

def tile_csv_into_pmtiles(file_in, file_out, args, logger, no_dup=True):
    tippecanoe_cmd = " ".join([
        f"TIPPECANOE_MAX_THREADS={args.threads}",
        f"'{args.tippecanoe}'",
        f"-t {args.tmp_dir}",
        f"-o {file_out}",
        f"-Z {args.min_zoom} -z {args.max_zoom} --force",
        f"-s EPSG:3857 -M {args.max_tile_bytes} -O {args.max_feature_counts}",
        f"--drop-densest-as-needed",
        f"--extend-zooms-if-still-dropping",
        f"'--preserve-point-density-threshold={args.preserve_point_density_thres}'",
        f"--no-duplication" if no_dup else "",
        f"--no-clipping",
        f"--buffer 0",
        f"{file_in}"
    ])

    logger.info(f"  * Running tippecanoe command: {tippecanoe_cmd}")
    result = subprocess.run(tippecanoe_cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        logger.error(f"Command {tippecanoe_cmd}\nfailed with error: {result.stderr.decode()}")
        sys.exit(1)
    else:
        logger.info("   * PMTiles creation command completed successfully")

def make_factor_dict(factor_id, factor_name, outprefix, pmtiles_keys=[]):
    pmtiles={}
    for key in pmtiles_keys:
        pmtiles[key]=f"{outprefix}-{key}.pmtiles"
    return {
        "id": factor_id,
        "name": factor_name,
        "cells_id": factor_id,
        "rgb": f"{outprefix}-rgb.tsv",
        "de": f"{outprefix}-cells-bulk-de.tsv",
        "raw_pixel_col": None,
        "pmtiles": pmtiles
    }


def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Import cell segmentation results from Xenium Ranger output")

    run_params = parser.add_argument_group("Run Options", "Run options for GNU Make")
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')
    run_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    run_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Enable all actions: --cells, --boundaries, and --summary')
    cmd_params.add_argument('--cells', action='store_true', default=False, help='Import segmented cells and generate PMTiles')
    cmd_params.add_argument('--boundaries', action='store_true', default=False, help='Import segmented cell boundaries and generate GeoJSON')
    cmd_params.add_argument('--summary', action='store_true', default=False, help='Write a JSON summary of parameters and output paths')
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update an existing catalog.yaml generated by run_cartload2 (not included in --all)')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'Two input modes: 1) JSON — set --in-json with keys (CELL, BOUNDARY, CLUSTER, DE); 2) Manual — set --indir and provide locations, see "Manual Input Parameters"')
    inout_params.add_argument('--in-json', type=str, help="Path to input JSON specifying paths for cells, boundaries, clusters, and differential expression.")
    inout_params.add_argument('--outprefix', type=str, required=True, help='Output prefix')
    inout_params.add_argument('--id', type=str, help='Identifier for the cell factor; if omitted, uses basename of --outprefix')
    inout_params.add_argument('--name', type=str, help='Display name for the cell factor; if omitted, uses basename of --outprefix')

    aux_inout_params = parser.add_argument_group("Manual Input Parameters", "Manually specify input directory (--indir) and input file locations under --indir")
    aux_inout_params.add_argument('--indir', type=str, help='Input directory containing the Xenium Ranger output files.')
    aux_inout_params.add_argument('--csv-cells', type=str, default="cells.csv.gz", help='Location of CSV or Parquet containing cell locations under --indir (default: cells.csv.gz)')
    aux_inout_params.add_argument('--csv-boundaries', type=str, default="cell_boundaries.csv.gz", help='Location of CSV containing cell boundary coordinates under --indir (default: cell_boundaries.csv.gz)')
    aux_inout_params.add_argument('--csv-clust', type=str, default="analysis/clustering/gene_expression_graphclust/clusters.csv", help='Location of CSV with cell cluster assignments under --indir (default: analysis/clustering/gene_expression_graphclust/clusters.csv)')
    aux_inout_params.add_argument('--csv-diffexp', type=str, default="analysis/diffexp/gene_expression_graphclust/differential_expression.csv", help='Location of CSV with differential expression results under --indir (default: analysis/diffexp/gene_expression_graphclust/differential_expression.csv)')

    aux_conv_params = parser.add_argument_group("Auxiliary PMTiles Conversion Parameters")
    aux_conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level (default: 10)')
    aux_conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level (default: 18)')
    aux_conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles (default: 5000000)')
    aux_conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles (default: 500000)')
    aux_conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles (default: 1024)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--tsv-cmap', type=str, default=f"{repo_dir}/assets/fixed_color_map_60.tsv", help=f'Location of TSV with color mappings for clusters under --indir (default: {repo_dir}/assets/fixed_color_map_60.tsv)')
    aux_params.add_argument('--de-max-pval', type=float, default=0.01, help='Maximum p-value for differential expression (default: 0.01)')
    aux_params.add_argument('--de-min-fc', type=float, default=1.2, help='Minimum fold change for differential expression (default: 1.2)')
    aux_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')
    aux_params.add_argument('--catalog-yaml', type=str, help='Path to catalog.yaml to update (used with --update-catalog; default: <out_dir>/catalog.yaml)')
    #aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory for intermediate files (default: out-dir/tmp or /tmp if specified)')

    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools binary. Required if a Parquet file is provided to --csv-cells (default: parquet-tools)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # Sanity check: JSON mode and manual mode are mutually exclusive
    if args.in_json and args.indir:
        parser.error("Cannot enable both JSON mode (--in-json) and manual input mode (--indir and/or --csv-cells/--csv-boundaries/--csv-clust/--csv-diffexp). Choose one.")

    return args

def import_xenium_cell(_args):
    """
    Import cell segmentation results from Xenium Ranger output
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.outprefix + "_importxenium" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    if args.all:
        args.cells = True
        args.boundaries = True
        args.summary = True
        # args.update_catalog = True
    
    if not args.cells and not args.boundaries and not args.summary and not args.update_catalog:
        raise ValueError("At least one action must be enabled.")

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
    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"File not found: {args.in_json} (--in-json)"
        raw_data=load_file_to_dict(args.in_json)
        cell_data=raw_data.get("CELLS", raw_data) # # use raw_data as default to support the flat dict build in the old scripts
    else:
        cell_data={
            "CELL": f"{args.indir}/{args.csv_cells}",
            "BOUNDARY": f"{args.indir}/{args.csv_boundaries}",
            "CLUSTER": f"{args.indir}/{args.csv_clust}",
            "DE": f"{args.indir}/{args.csv_diffexp}",
        }

    # Cluster/DE
    if args.cells or args.boundaries:
        clust_in = cell_data.get("CLUSTER", None)
        # presence and existence for cluster file
        assert clust_in is not None, ('Path not provided: "CLUSTER" in --in-json' if args.in_json is not None else 'Path not provided: --csv-clust')
        assert os.path.exists(clust_in), (f'File not found: {clust_in} ("CLUSTER" in --in-json)' if args.in_json is not None else f'File not found: {clust_in} (--csv-clust)')

        logger.info(f"Loading cell cluster data from {clust_in}")
        sorted_clusters, cluster2idx, bcd2cluster, bcd2clusteridx=process_cluster_csv(clust_in)
        logger.info(f"  * Loaded {len(bcd2cluster)} cells")
        
        ## read/write DE results
        de_in = cell_data.get("DE", None)
        # presence and existence for differential expression file
        assert de_in is not None, ('Path not provided: "DE" in --in-json' if args.in_json is not None else 'Path not provided: --csv-diffexp')
        assert os.path.exists(de_in), (f'File not found: {de_in} ("DE" in --in-json)' if args.in_json is not None else f'File not found: {de_in} (--csv-diffexp)')        

        logger.info(f"  * Reading DE results from {de_in}")
        clust2genes=read_de_csv(de_in, cluster2idx, args.de_min_fc, args.de_max_pval)

        ## write DE results
        de_out=f"{args.outprefix}-cells-bulk-de.tsv"
        logger.info(f"  * Writing DE results for {len(clust2genes)} clusters) to {de_out}")
        write_de_tsv(clust2genes, de_out, sorted_clusters)

        ## write the color map
        cmap_out=f"{args.outprefix}-rgb.tsv"
        assert os.path.exists(args.tsv_cmap), f"File not found: {args.tsv_cmap} (--tsv-cmap)"        

        logger.info(f"  * Writing color map from {args.tsv_cmap} to {cmap_out}")
        write_cmap_tsv(cmap_out, args.tsv_cmap, sorted_clusters)

    # Process segmented calls 
    if args.cells:
        cells_in = cell_data.get("CELL", None)
        # presence and existence for cells file
        assert cells_in is not None, ('Path not provided: "CELL" in --in-json' if args.in_json is not None else 'Path not provided: --csv-cells')
        assert os.path.exists(cells_in), (f'File not found: {cells_in} ("CELL" in --in-json)' if args.in_json is not None else f'File not found: {cells_in} (--csv-cells)')
        
        logger.info(f"Processing cell information from {cells_in}")

        # convert cell parquet into csv format
        if cells_in.endswith(".parquet"):
            cells_parquet=cells_in
            cells_csv = os.path.join(out_dir, "cells.csv.gz")
            par2csv_cmd=f"{args.parquet_tools} csv {cells_parquet} |  gzip -c > {cells_csv}"
            
            logger.info(f"  * Converting {cells_in} from parquet to CSV: {par2csv_cmd}")
            result = subprocess.run(par2csv_cmd, shell=True, capture_output=True)
            if result.returncode != 0:
                logger.error(f"Command {par2csv_cmd}\nfailed with error: {result.stderr.decode()}")
                sys.exit(1)
            # else:
            #     temp_fs.append(cells_csv)
        else:
            cells_csv = cells_in
        
        ## create a cell output file with cell_id, x, y, count
        cells_out=f"{args.outprefix}-cells.csv"
        logger.info(f"  * Reading cell data from {cells_csv} and extracting geometry to {cells_out}")
        process_cells_csv(cells_csv, cells_out, bcd2clusteridx)

        cells_pmtiles = f"{args.outprefix}-cells.pmtiles"
        logger.info(f"  * Generating PMTiles from cell geometry data into cell pmtiles: {cells_pmtiles}")
        tile_csv_into_pmtiles(cells_out, cells_pmtiles, args, logger, no_dup=True)
        #temp_fs.append(cells_out)
        

    # Process cell boundaries
    if args.boundaries:
        bound_in = cell_data.get("BOUNDARY", None)
        # presence and existence for boundary file
        assert bound_in is not None, ('Path not provided: "BOUNDARY" in --in-json' if args.in_json is not None else 'Path not provided: --csv-boundaries')
        assert os.path.exists(bound_in), (f'File not found: {bound_in} ("BOUNDARY" in --in-json)' if args.in_json is not None else f'File not found: {bound_in} (--csv-boundaries)')

        logger.info(f"Processing cell boundary information from {bound_in}")

        ## read the cell CSV files
        bound_out=f"{args.outprefix}-boundaries.geojson"
        logger.info(f"  * Reading cell boundary data from {bound_in} and extracting geometry to {bound_out}")
        with flexopen(bound_in, "rt") as f:
            with flexopen(bound_out, "wt") as wf:
                current_cell_id = None
                current_vertices = []
                hdrs = f.readline().rstrip().split(",")
                assert unquote_str(hdrs[0]) == "cell_id"
                assert unquote_str(hdrs[1]) == "vertex_x"
                assert unquote_str(hdrs[2]) == "vertex_y"
                for line in f:
                    toks = line.strip().split(",")
                    cell_id = unquote_str(toks[0])
                    x = float(toks[1])
                    y = float(toks[2])
                    if current_cell_id is None:
                        current_cell_id = cell_id
                    if cell_id != current_cell_id:
                        cluster = bcd2cluster.get(current_cell_id, "NA")
                        clusteridx = cluster2idx.get(cluster, "NA") if cluster != "NA" else "NA"
                        wf.write(f'{{"type": "Feature", "geometry": {{"type": "Polygon", "coordinates": [[{",".join(current_vertices)}]]}}, "properties": {{"cell_id": "{current_cell_id}", "topK": "{clusteridx}"}}}}\n')
                        current_cell_id = cell_id
                        current_vertices = []
                    current_vertices.append(f'[{x},{y}]')
                if current_cell_id is not None:
                    cluster = bcd2cluster.get(current_cell_id, "NA")
                    clusteridx = cluster2idx.get(cluster, "NA") if cluster != "NA" else "NA"
                    wf.write(f'{{"type": "Feature", "geometry": {{"type": "Polygon", "coordinates": [[{",".join(current_vertices)}]]}}, "properties": {{"cell_id": "{current_cell_id}", "topK": "{clusteridx}"}}}}\n')

        ## Run the tippecanoe command
        bound_pmtiles = f"{args.outprefix}-boundaries.pmtiles"
        logger.info(f"  * Generating PMTiles from boundary geometry data into boundary pmtiles: {bound_pmtiles}")
        tile_csv_into_pmtiles(bound_out, bound_pmtiles, args, logger, no_dup=False)
        #temp_fs.append(bound_out)
    
    # JSON/YAML
    if args.summary or args.update_catalog:
        ## add the new factor to the catalog
        factor_id = out_base if args.id is None else args.id
        factor_name = out_base if args.name is None else args.name
        pmtiles_keys=[]
        if os.path.exists(f"{args.outprefix}-cells.pmtiles"):
            pmtiles_keys.append("cells")
        if os.path.exists(f"{args.outprefix}-boundaries.pmtiles"):
            pmtiles_keys.append("boundaries")

    if args.summary:
        out_assets_f=f"{args.outprefix}_assets.json"
        logger.info(f"Summarizing assets information into {out_assets_f}")
        new_factor = make_factor_dict(factor_id, factor_name, args.outprefix, pmtiles_keys)
        write_dict_to_file(new_factor, out_assets_f, check_equal=True)

    if args.update_catalog:
        ## read the input catalog.yaml file
        if args.catalog_yaml is None:
            args.catalog_yaml = os.path.join(args.indir, "catalog.yaml")
        
        logger.info(f"Updating catalog YAML file: {args.catalog_yaml}")
        ##TO-DO: catalog_yaml should be within the same dir with the output assets

        ## load the YAML file
        with open(args.catalog_yaml, 'r') as f:
            catalog = yaml.load(f, Loader=yaml.FullLoader)  # Preserves order

        ## add files to the catalog
        new_factor = make_factor_dict(factor_id, factor_name, out_base, pmtiles_keys)
        if "factors" not in catalog["assets"]:
            raise ValueError("No factors found in the catalog.yaml file. Check if the file is correct.")
        catalog["assets"]["factors"].append(new_factor)

        ## write the updated catalog.yaml file
        out_yaml = f"{args.outprefix}-catalog.yaml"
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

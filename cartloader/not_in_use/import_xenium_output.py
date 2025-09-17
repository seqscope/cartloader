import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, csv, yaml, inspect
import pandas as pd
import numpy as np
from scipy.stats import chi2

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, flexopen, unquote_str, smartsort, write_dict_to_file, load_file_to_dict

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Import cell segmentation results from Xenium Ranger output")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Run all commands (cells, boundaries, summary)')
    cmd_params.add_argument('--cells', action='store_true', default=False, help='Add segmented cells to PMTiles output as factors')
    cmd_params.add_argument('--boundaries', action='store_true', default=False, help='Add segmented cell bounaries to PMTiles output')
    cmd_params.add_argument('--summary', action='store_true', default=False, help='Generate a JSON file summarizing parameters and output paths.')
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update the YAML files (NOT included in --all)')

    inout_params = parser.add_argument_group("Input/Output Parameters", """
                                             Input/output directory/files. There are two ways to specify the input files. 
                                             One is to use --in-json to provide the paths. 
                                             The other is to use --indir with auxiliary input parameters (--cells, --csv-boundaries, --csv-clust, --csv-diffexp) to specify input paths """)
    inout_params.add_argument('--in-json', type=str, help="Provide a json file specified the path for the cells, boundaries, cluster, diffexp")
    inout_params.add_argument('--indir', type=str, help='Input directory containing the Xenium Ranger output files. If --in-json is specified, skip this argument')
    inout_params.add_argument('--outprefix', type=str, help='Prefix of output files')
    inout_params.add_argument('--id', type=str, help='Identifier of the factor')
    inout_params.add_argument('--name', type=str, help='Name of the factor')

    conv_params = parser.add_argument_group("Parameters for pmtiles conversion")
    conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level (default: 10)')
    conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level (default: 18)')
    conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles (default: 5000000)')
    conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles (default: 500000)')
    conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles (default: 1024)')

    run_params = parser.add_argument_group("Run Options", "Run options for GNU Make")
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')
    run_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    run_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--catalog-yaml', type=str, help='YAML file to be updated when --update-catalog is enabled')
    aux_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')
    aux_params.add_argument('--csv-cells', type=str, default="cells.csv.gz", help='Location of the CSV or Parquet file containing cell locations in the --indir (default: cells.csv.gz)')
    aux_params.add_argument('--csv-boundaries', type=str, default="cell_boundaries.csv.gz", help='Location of the CSV file containing cell boundary coordinates in the --indir (default: cell_boundaries.csv.gz)')
    aux_params.add_argument('--csv-clust', type=str, default="analysis/clustering/gene_expression_graphclust/clusters.csv", help='Location of the CSV file containing cell cluster assignments in the --indir (default: analysis/clustering/gene_expression_graphclust/clusters.csv)')
    aux_params.add_argument('--csv-diffexp', type=str, default="analysis/diffexp/gene_expression_graphclust/differential_expression.csv", help='Location of the CSV file with differential expression results in the --indir (default: analysis/diffexp/gene_expression_graphclust/differential_expression.csv)')
    aux_params.add_argument('--tsv-cmap', type=str, default=f"{repo_dir}/assets/fixed_color_map_60.tsv", help=f'Location of the TSV file with color mappings for clusters in the --indir (default: {repo_dir}/assets/fixed_color_map_60.tsv)')
    aux_params.add_argument('--de-max-pval', type=float, default=0.01, help='Maximum p-value threshold for differential expression (default: 0.01)')
    aux_params.add_argument('--de-min-fc', type=float, default=1.2, help='Minimum fold change threshold for differential expression (default: 1.2)')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files (default: False)')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory for intermediate files (default: out-dir/tmp or /tmp if specified)')

    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools binary. Required if a parquet file is provided to --csv-cells (defaults: parquet-tools)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def import_xenium_output(_args):
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

    # create output directory if needed
    out_dir = os.path.dirname(args.outprefix)
    out_base = os.path.basename(args.outprefix)
    if not os.path.exists(out_dir) and out_dir != "":
        os.makedirs(out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(out_dir, "tmp")
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir, exist_ok=True)

    # read in_json if provided 
    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"The input json file doesn't exist: {args.in_json}"
        cellbounds=load_file_to_dict(args.in_json)
    else:
        cellbounds={
            "CELL": f"{args.indir}/{args.csv_cells}",
            "BOUNDARY": f"{args.indir}/{args.csv_boundaries}",
            "CLUSTER": f"{args.indir}/{args.csv_clust}",
            "DE": f"{args.indir}/{args.csv_diffexp}",
        }

    clust_f=cellbounds.get("CLUSTER", None)
    assert os.path.exists(clust_f), f"The input cluster file doesn't exist {clust_f}"

    ## read cluster information 
    logger.info(f"Reading the cell cluster information from {clust_f}")
    bcd2cluster = {}
    cluster2cnt = {}
    with flexopen(clust_f, "rt") as f:
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

    logger.info(f"Read {len(bcd2cluster)} cellss")
    
    # Process segmented calls 
    if args.cells:
        cells_f=cellbounds.get("CELL", None)
        assert os.path.exists(cells_f), f"The input cells file doesn't exist {cells_f}"

        if cells_f.endswith(".parquet"):
            cells_parquet=cells_f
            cells_csv = os.path.join(out_dir, "cells.csv.gz")
            par2csv_cmd=f"{args.parquet_tools} csv {cells_parquet} |  gzip -c > {cells_csv}"
            logger.info(f"Converting {cells_f} from parquet to CSV: {par2csv_cmd}")
            result = subprocess.run(par2csv_cmd, shell=True, capture_output=True)
            if result.returncode != 0:
                logger.error(f"Command {par2csv_cmd}\nfailed with error: {result.stderr.decode()}")
                sys.exit(1)
            else:
                logger.info("Cells format conversion command completed successfully")
        else:
            cells_csv = cells_f
        
        ## read the cell boundary CSV files
        logger.info(f"Reading the cell CSV file {cells_csv}")
        with flexopen(cells_csv, "rt") as f:
            with flexopen(f"{args.outprefix}-cells.csv", "wt") as wf:
                wf.write(",".join(["lon","lat","cell_id","count","topK"]) + "\n")
                reader = csv.DictReader(f)
                for row in reader:
                    cell_id = row["cell_id"]
                    x = row["x_centroid"]
                    y = row["y_centroid"]
                    transcript_count = row["transcript_counts"]
                    cluster = bcd2cluster.get(cell_id, "NA")
                    clusteridx = cluster2idx.get(cluster, "NA") if cluster != "NA" else "NA"
                    wf.write(",".join([x, y, cell_id, transcript_count, str(clusteridx)]) + "\n")

        ## Run the tippecanoe command
        pmtiles_path = f"{args.outprefix}-cells.pmtiles"
        tippecanoe_cmd = f"TIPPECANOE_MAX_THREADS={args.threads} '{args.tippecanoe}' -t {args.tmp_dir} -o {pmtiles_path} -Z {args.min_zoom} -z {args.max_zoom} --force -s EPSG:3857 -M {args.max_tile_bytes} -O {args.max_feature_counts} --drop-densest-as-needed --extend-zooms-if-still-dropping '--preserve-point-density-threshold={args.preserve_point_density_thres}' --no-duplication --no-clipping --buffer 0 {args.outprefix}-cells.csv"
        logger.info(f"Running tippecanoe command: {tippecanoe_cmd}")
        result = subprocess.run(tippecanoe_cmd, shell=True, capture_output=True)
        if result.returncode != 0:
            logger.error(f"Command {tippecanoe_cmd}\nfailed with error: {result.stderr.decode()}")
            sys.exit(1)
        else:
            logger.info("PMTiles creation command completed successfully")

        ## read DE results
        de_f=cellbounds.get("DE", None)
        assert os.path.exists(de_f), f"The input differentially expressed profile file doesn't exist {de_f}"

        logger.info(f"Reading the DE results from {de_f}")
        clust2genes = {}
        with flexopen(de_f, "rt") as f:
            hdrs = f.readline().strip().split(",")
            assert hdrs[1] == "Feature Name" ## make sure the header is correct
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
                    if fc > args.de_min_fc and pval < args.de_max_pval:
                        if clusteridx not in clust2genes:
                            clust2genes[clusteridx] = []
                        clust2genes[clusteridx].append([gene_name, avgcount, pval, fc])
        ## write DE results
        logger.info(f"Writing the DE results to {args.outprefix}-cells-bulk-de.tsv")
        with flexopen(f"{args.outprefix}-cells-bulk-de.tsv", "wt") as wf:
            wf.write("\t".join(["gene","factor","Chi2","pval","FoldChange","gene_total","log10pval"]) + "\n")
            for clusteridx in range(len(sorted_clusters)):
                if clusteridx in clust2genes:
                    genes = clust2genes[clusteridx]
                    sorted_genes = sorted(genes, key=lambda x: x[3], reverse=True)
                    for gene in sorted_genes:
                        (gname, avgcount, pval, fc) = gene
                        if pval == 0:
                            log10pval = 320
                            chisq = 1465.911
                        else:
                            log10pval = -np.log10(pval)
                            chisq = chi2.isf(pval, 1)
                        wf.write("\t".join([gname, str(clusteridx), f"{chisq:.4f}", f"{pval:.4e}", f"{fc:.4f}", f"{avgcount:.4f}", f"{log10pval:.4f}"]) + "\n")
        logger.info(f"Successfully written DE results")

        ## write the color map
        logger.info(f"Writing the color map to {args.outprefix}-rgb.tsv")
        with flexopen(f"{args.outprefix}-rgb.tsv", "wt") as wf:
            with flexopen(args.tsv_cmap, "rt") as f:
                wf.write(f.readline())
                for i in range(len(sorted_clusters)):
                    line = f.readline()
                    if line == "" or line is None:
                        raise ValueError(f"Not enough colors in the color map file {args.tsv_cmap}")
                    wf.write(line)
        logger.info("Color map written successfully")

    # Process cell boundaries
    if args.boundaries:
        bound_f=cellbounds.get("BOUNDARY", None)
        assert os.path.exists(bound_f), f"The input cell boundary files doesn't exist {bound_f}"

        ## read the cell CSV files
        logger.info(f"Reading the cell CSV file {bound_f}")
        with flexopen(bound_f, "rt") as f:
            with flexopen(f"{args.outprefix}-boundaries.geojson", "wt") as wf:
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
        logger.info(f"Successfully written cell boundaries")

        ## Run the tippecanoe command
        pmtiles_path = f"{args.outprefix}-boundaries.pmtiles"
        tippecanoe_cmd = f"TIPPECANOE_MAX_THREADS={args.threads} '{args.tippecanoe}' -t {args.tmp_dir} -o {pmtiles_path} -Z {args.min_zoom} -z {args.max_zoom} --force -s EPSG:3857 -M {args.max_tile_bytes} -O {args.max_feature_counts} --drop-densest-as-needed --extend-zooms-if-still-dropping '--preserve-point-density-threshold={args.preserve_point_density_thres}' --no-clipping --buffer 0 {args.outprefix}-boundaries.geojson"
        logger.info(f"Running tippecanoe command: {tippecanoe_cmd}")
        result = subprocess.run(tippecanoe_cmd, shell=True, capture_output=True)
        if result.returncode != 0:
            logger.error(f"Command {tippecanoe_cmd}\nfailed with error: {result.stderr.decode()}")
            sys.exit(1)
        else:
            logger.info("PMTiles creation command completed successfully")

    # JSON/YAML
    if args.summary or args.update_catalog:
        ## add the new factor to the catalog
        factor_id = out_base if args.id is None else args.id
        factor_name = out_base if args.name is None else args.name

    if args.summary:
        out_assets_f=f"{args.outprefix}_assets.json"
        new_factor = {
            "id": factor_id,
            "name": factor_name,
            "cells_id": factor_id,
            "rgb": f"{args.outprefix}-rgb.tsv",
            "de": f"{args.outprefix}-cells-bulk-de.tsv",
            "pmtiles": {
                "cells": f"{args.outprefix}-cells.pmtiles",
                "boundaries": f"{args.outprefix}-boundaries.pmtiles"
            }
        }
        write_dict_to_file(new_factor, out_assets_f, check_equal=True)

    if args.update_catalog:
        logger.info(f"Updating the YAML files with the results")

        ## read the input catalog.yaml file
        if args.catalog_yaml is None:
            args.catalog_yaml = os.path.join(args.indir, "catalog.yaml")

        ## load the YAML file
        with open(args.catalog_yaml, 'r') as f:
            catalog = yaml.load(f, Loader=yaml.FullLoader)  # Preserves order

        ## add files to the catalog
        new_factor = {
            "id": factor_id,
            "name": factor_name,
            "cells_id": factor_id,
            "rgb": f"{out_base}-rgb.tsv",
            "de": f"{out_base}-cells-bulk-de.tsv",
            "pmtiles": {
                "cells": f"{out_base}-cells.pmtiles",
                "boundaries": f"{out_base}-boundaries.pmtiles"
            }
        }

        if "factors" not in catalog["assets"]:
            raise ValueError("No factors found in the catalog.yaml file. Check if the file is correct.")
        catalog["assets"]["factors"].append(new_factor)

        ## write the updated catalog.yaml file
        out_yaml = f"{args.outprefix}-catalog.yaml"
        with open(out_yaml, 'w') as f:
            yaml.dump(catalog, f, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)
        logger.info(f"Successfully wrote the catalog.yaml file: {out_yaml}")

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

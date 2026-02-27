import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, csv, yaml, inspect
import pandas as pd
import numpy as np
from scipy.stats import chi2
import subprocess

from cartloader.utils.utils import scheck_app, create_custom_logger, flexopen, unquote_str, smartsort, write_dict_to_file, load_file_to_dict, run_command

from cartloader.utils.color_helper import normalize_rgb, rgb_to_hex

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def process_cluster_csv(clust_csv, barcode_col="Barcode", cluster_col="Cluster", output_filename=None):
    bcd2cluster = {}
    cluster2cnt = {}
    with flexopen(clust_csv, "rt") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bcd = unquote_str(row[barcode_col])
            clust = row[cluster_col]
            bcd2cluster[bcd] = clust
            cluster2cnt[clust] = cluster2cnt.get(clust, 0) + 1
    
    sorted_clusters = sorted(cluster2cnt.keys(), key=lambda x: cluster2cnt[x], reverse=True)
    cluster2idx = {cluster: idx for idx, cluster in enumerate(sorted_clusters)}
    bcd2clusteridx = {bcd: cluster2idx[clust] for bcd, clust in bcd2cluster.items()}

    if output_filename is not None:
        with flexopen(output_filename, "wt") as wf:
            wf.write("cell_id\tcluster\n")
            for bcd, clustidx in bcd2clusteridx.items():
                wf.write(f"{bcd}\t{clustidx}\n")
                
    return sorted_clusters, cluster2idx, bcd2clusteridx

def process_cells_csv(cells_csv, out_csv, bcd2clusteridx, cell_id_col="cell_id", x_col="x_centroid", y_col="y_centroid", count_col="transcript_counts"):
    with flexopen(cells_csv, "rt") as f, flexopen(out_csv, "wt") as wf:
        reader = csv.DictReader(f)
        required_cols = {cell_id_col, x_col, y_col, count_col}
        if not required_cols.issubset(reader.fieldnames):
            raise ValueError(f"Missing expected columns: {required_cols - set(reader.fieldnames)}")
        wf.write("lon,lat,cell_id,count,topK\n")
        for row in reader:
            cell_id = row[cell_id_col]
            x = row[x_col]
            y = row[y_col]
            count = row[count_col]
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

# def write_cmap_tsv(out_cmap, tsv_cmap, sorted_clusters):
#     with flexopen(out_cmap, "wt") as wf:
#         with flexopen(tsv_cmap, "rt") as f:
#             # Write header line
#             wf.write(f.readline())

#             # Write color lines for each cluster
#             for i in range(len(sorted_clusters)):
#                 line = f.readline()
#                 if not line:
#                     raise ValueError(
#                         f"Not enough colors in the color map file {tsv_cmap}"
#                     )
#                 wf.write(line)

# this rgb will be 0-1
def write_cmap_tsv(out_cmap, tsv_cmap, sorted_clusters):
    with flexopen(tsv_cmap, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    required = ["R", "G", "B", "Name"]
    for col in required:
        assert col in fieldnames, f"Color map file {tsv_cmap} is missing required column: {col}"

    assert len(rows) >= len(sorted_clusters), f"Not enough colors in {tsv_cmap}: needed {len(sorted_clusters)}, found {len(rows)}"
    has_color_hex = "Color_hex" in fieldnames
    if not has_color_hex:
        fieldnames.append("Color_hex")

    for row in rows:
        r = float(row["R"])
        g = float(row["G"])
        b = float(row["B"])
        r01, g01, b01 = normalize_rgb(r, g, b)
        row["R"] = f"{r01:.6f}"
        row["G"] = f"{g01:.6f}"
        row["B"] = f"{b01:.6f}"

        if not has_color_hex:
            row["Color_hex"] = rgb_to_hex(r01, g01, b01)

    with flexopen(out_cmap, "wt") as wf:
        writer = csv.DictWriter(wf, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for i in range(len(sorted_clusters)):
            writer.writerow(rows[i])

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

def write_umap_tsv(umap_in, umap_tsv_out, bcd2clusteridx, args):
    umap_df = pd.read_csv(umap_in)
    for col in [args.umap_colname_barcode, args.umap_colname_x, args.umap_colname_y]:
        if col not in umap_df.columns:
            raise ValueError(f"Column '{col}' not found in UMAP projection file: {umap_in}")

    # here use idx to ensure the cluster is the reclustered one (after ranking by size)
    umap_df['topK'] = umap_df[args.umap_colname_barcode].map(bcd2clusteridx).fillna("NA").astype(str)
    if args.umap_colname_barcode == "Barcode":
        umap_df.rename(columns={args.umap_colname_x: "UMAP1", args.umap_colname_y: "UMAP2"}, inplace=True)
    else:
        umap_df.rename(columns={args.umap_colname_barcode: "Barcode", args.umap_colname_x: "UMAP1", args.umap_colname_y: "UMAP2"}, inplace=True)
    
    umap_df = umap_df[['Barcode', 'topK', 'UMAP1', 'UMAP2']]
    umap_df.to_csv(umap_tsv_out, index=False, compression='gzip', sep="\t")

def umap_tsv2pmtiles(umap_tsv_out, umap_pmtiles, args):
    # 1) Use NDJSON (Newline-Delimited JSON) than JSON: 1) this file is light; 2) tippecanoe can process each line without reading all info
    umap_ndjson=umap_pmtiles.replace(".pmtiles", ".ndjson")
    convert_cmd = " ".join([
        "cartloader", "render_umap",
        f"--input {umap_tsv_out}",
        f"--out {umap_ndjson}",
        f"--colname-factor topK",
        f"--colname-x UMAP1",
        f"--colname-y UMAP2"
    ])
    run_command(convert_cmd)
    
    # 2) ndjson to pmtiles
    tippecanoe_cmd = " ".join([
        f"TIPPECANOE_MAX_THREADS={args.threads}",
        f"'{args.tippecanoe}'",
        f"-t '{args.tmp_dir}'",
        f"-o '{umap_pmtiles}'",
        "-Z", str(args.umap_min_zoom),
        "-z", str(args.umap_max_zoom),
        "-l", "umap",
        "--force",
        "--drop-densest-as-needed",
        "--extend-zooms-if-still-dropping",
        "--no-duplication",
        f"--preserve-point-density-threshold={args.preserve_point_density_thres}",
        umap_ndjson
    ])
    run_command(tippecanoe_cmd)

def umap_tsv2png(umap_tsv, model_prefix, color_map, title="Cell Segmentation"):
    draw_umap_rscript=f"{repo_dir}/cartloader/r/draw_umap.r"

    plot_cmd = " ".join([
        f"Rscript '{draw_umap_rscript}'",
        f"--input '{umap_tsv}'",
        f"--out-prefix '{model_prefix}'",
        f"--cmap '{color_map}'",
        f'--subtitle \"{title}\"',
        ])
    run_command(plot_cmd)

def umap_tsv2indpng(umap_tsv, model_prefix, color_map, mode="binary", title="Cell Segmentation"):
    draw_umap_single_rscript=f"{repo_dir}/cartloader/r/draw_umap_single.r"

    plot_cmd = " ".join([
        f"Rscript '{draw_umap_single_rscript}'",
        f"--input '{umap_tsv}'",
        f"--out-prefix '{model_prefix}'",
        f"--cmap '{color_map}'",
        f'--subtitle  \"{title}\"',
        f"--mode {mode}"
        ])
    run_command(plot_cmd)

def make_factor_dict(factor_id, factor_name, outprefix, factor_type, pmtiles_keys=[], umap_src=False, pseudobulk_src=False):
    assert factor_type in ["cells", "square"], "Currently only support cells and square type..."
    pmtiles={}
    for key in pmtiles_keys:
        if factor_type=="cells":
            pmtiles[key]=f"{outprefix}-{key}.pmtiles"
        elif factor_type=="square":
            pmtiles["boundaries"]=f"{outprefix}.pmtiles"
    model_label=f"{factor_type}_id"
    factor_dict ={
        "id": factor_id,
        "name": factor_name,
        model_label: factor_id,
        "rgb": f"{outprefix}-rgb.tsv",
        "de": f"{outprefix}-cells-bulk-de.tsv" if factor_type=="cells" else f"{outprefix}-bulk-de.tsv",
        "raw_pixel_col": False,
        "pmtiles": pmtiles,
    }
    if pseudobulk_src:
        factor_dict["post"] = f"{outprefix}-pseudobulk.tsv.gz"
        factor_dict["info"] = f"{outprefix}-info.tsv"
    if umap_src:
        factor_dict["umap"] = {
            "tsv": f"{outprefix}-umap.tsv.gz",
            "png": f"{outprefix}.umap.png",
            "ind_png": f"{outprefix}.umap.single.binary.png",
            "pmtiles": f"{outprefix}-umap.pmtiles"
        }
    return factor_dict

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
    cmd_params.add_argument('--all', action='store_true', default=False, help='Enable all actions: --cells and --boundaries')
    cmd_params.add_argument('--cells', action='store_true', default=False, help='Import segmented cells and generate PMTiles')
    cmd_params.add_argument('--boundaries', action='store_true', default=False, help='Import segmented cell boundaries and generate GeoJSON')
    cmd_params.add_argument('--umap', action='store_true', default=False, help='Import UMAP projection and generate PMTiles for UMAP points')
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update an existing catalog.yaml generated by run_cartload2 (not included in --all)')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'Two input modes: 1) JSON — set --in-json with keys (CELL, BOUNDARY, CLUSTER, DE); 2) Manual — set --in-dir and provide locations, see "Manual Input Parameters"')
    inout_params.add_argument('--in-json', type=str, help="Path to input JSON specifying paths for cells, boundaries, clusters, and differential expression.")
    inout_params.add_argument('--outprefix', type=str, required=True, help='Output prefix')
    inout_params.add_argument('--id', type=str, help='Identifier for the cell factor; if omitted, uses basename of --outprefix')
    inout_params.add_argument('--name', type=str, help='Display name for the cell factor; if omitted, uses basename of --outprefix')

    aux_inout_params = parser.add_argument_group("Manual Input Parameters", "Manually specify input directory (--in-dir) and input file locations under --in-dir")
    aux_inout_params.add_argument('--in-dir', type=str, help='Input directory containing the Xenium Ranger output files.')
    aux_inout_params.add_argument('--csv-cells', type=str, default="cells.csv.gz", help='Location of CSV or Parquet containing cell locations under --in-dir (default: cells.csv.gz)')
    aux_inout_params.add_argument('--csv-boundaries', type=str, default="cell_boundaries.csv.gz", help='Location of CSV containing cell boundary coordinates under --in-dir (default: cell_boundaries.csv.gz)')
    aux_inout_params.add_argument('--csv-clust', type=str, default="analysis/clustering/gene_expression_graphclust/clusters.csv", help='Location of CSV with cell cluster assignments under --in-dir (default: analysis/clustering/gene_expression_graphclust/clusters.csv)')
    aux_inout_params.add_argument('--csv-diffexp', type=str, default="analysis/diffexp/gene_expression_graphclust/differential_expression.csv", help='Location of CSV with differential expression results under --in-dir (default: analysis/diffexp/gene_expression_graphclust/differential_expression.csv)')
    ## cell-level MEX format files
    aux_inout_params.add_argument('--mex-dir', type=str, default="cell_feature_matrix", help='Directory location of 10x Genomic MatrixMarket files under --in-dir (default: cell_feature_matrix)')
    aux_inout_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='Filename for barcodes in the MatrixMarket directory (default: barcodes.tsv.gz)')
    aux_inout_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='Filename for features in the MatrixMarket directory (default: features.tsv.gz)')
    aux_inout_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='Filename for matrix in the MatrixMarket directory (default: matrix.mtx.gz)')
    ## pixel-level TSV file with cell_id
    aux_inout_params.add_argument('--pixel', type=str, help='Location of pixel-level TSV file with cell_id under --in-dir')
    aux_inout_params.add_argument('--pixel-colname-cell-id', type=str, default='cell_id', help='Column name for cell ID in --pixel (default: cell_id)')
    aux_inout_params.add_argument('--pixel-colname-gene', type=str, default='gene', help='Column name for gene in --pixel (default: gene)')
    aux_inout_params.add_argument('--pixel-colname-count', type=str, default='count', help='Column name for count in --pixel (default: count)')
    ## UMAP
    aux_inout_params.add_argument('--csv-umap', type=str, default="analysis/umap/gene_expression_2_components/projection.csv", help='Location of CSV with UMAP results under --in-dir (default: analysis/umap/gene_expression_2_components/projection.csv')

    aux_conv_params = parser.add_argument_group("Auxiliary PMTiles Conversion Parameters")
    aux_conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level for cells and boundaries (default: 10)')
    aux_conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level for cells and boundaries (default: 18)')
    aux_conv_params.add_argument('--umap-min-zoom', type=int, default=0, help='Minimum zoom level for UMAP (default: 0)')
    aux_conv_params.add_argument('--umap-max-zoom', type=int, default=18, help='Maximum zoom level for UMAP (default: 18)')
    aux_conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles (default: 5000000)')
    aux_conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles (default: 500000)')
    aux_conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles (default: 1024)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--skip-redo-diffexp', action='store_true', default=False, help='Skip computing differential expression from the MatrixMarket files. Use files from --csv-diffexp instead.') 
    aux_params.add_argument('--skip-redo-pseudobulk', action='store_true', default=False, help='Skip generating pseudobulk files from the MatrixMarket files.')
    aux_params.add_argument('--tsv-cmap', type=str, default=f"{repo_dir}/assets/fixed_color_map_256.tsv", help=f'Location of TSV with color mappings for clusters under --in-dir (default: {repo_dir}/assets/fixed_color_map_60.tsv)')
    aux_params.add_argument('--de-max-pval', type=float, default=0.01, help='Maximum p-value for differential expression (default: 0.01)')
    aux_params.add_argument('--de-min-fc', type=float, default=1.2, help='Minimum fold change for differential expression (default: 1.2)')
    aux_params.add_argument('--catalog-yaml', type=str, help='Path to catalog.yaml to update (used with --update-catalog; default: <out_dir>/catalog.yaml)')
    aux_params.add_argument('--out-catalog-yaml', type=str, help='Path to save the updated catalog.yaml as a new file instead of overwriting the input (--catalog-yaml). Defaults to the same path as --catalog-yaml (used with --update-catalog)')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory for intermediate files (default: out-dir/tmp or /tmp if specified)')

    aux_colnames_params = parser.add_argument_group("Auxiliary Colname Parameters", "Override column names for input files")
    aux_colnames_params.add_argument('--cells-colname-cell-id', type=str, default='cell_id', help='Column name for cell ID in --csv-cells (default: cell_id)')
    aux_colnames_params.add_argument('--cells-colname-x', type=str, default='x_centroid', help='Column name for X coordinate in --csv-cells (default: x_centroid)')
    aux_colnames_params.add_argument('--cells-colname-y', type=str, default='y_centroid', help='Column name for Y coordinate in --csv-cells (default: y_centroid)')
    aux_colnames_params.add_argument('--cells-colname-count', type=str, default='transcript_counts', help='Column name for transcript count in --csv-cells (default: transcript_counts)')
    aux_colnames_params.add_argument('--clust-colname-barcode', type=str, default='Barcode', help='Column name for cell barcode in --csv-clust (default: Barcode)')
    aux_colnames_params.add_argument('--clust-colname-cluster', type=str, default='Cluster', help='Column name for cluster label in --csv-clust (default: Cluster)')
    aux_colnames_params.add_argument('--boundaries-colname-cell-id', type=str, default='cell_id', help='Column name for cell ID in --csv-boundaries (default: cell_id)')
    aux_colnames_params.add_argument('--boundaries-colname-x', type=str, default='vertex_x', help='Column name for X coordinate in --csv-boundaries (default: vertex_x)')
    aux_colnames_params.add_argument('--boundaries-colname-y', type=str, default='vertex_y', help='Column name for Y coordinate in --csv-boundaries (default: vertex_y)')
    aux_colnames_params.add_argument('--umap-colname-barcode', type=str, default='Barcode', help='Column name for barcode in --csv-umap (default: Barcode)')
    aux_colnames_params.add_argument('--umap-colname-x', type=str, default='UMAP-1', help='Column name for UMAP X coordinate in --csv-umap (default: UMAP-1)')
    aux_colnames_params.add_argument('--umap-colname-y', type=str, default='UMAP-2', help='Column name for UMAP Y coordinate in --csv-umap (default: UMAP-2)')

    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')
    env_params.add_argument('--spatula', type=str, default=f"{repo_dir}/submodules/spatula/bin/spatula", help='Path to spatula binary (default: <cartloader_dir>/submodules/spatula/bin/spatula)')
    env_params.add_argument('--python', type=str, default="python3",  help='Python3 binary')
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"), help='Path to punkst (ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools binary. Required if a Parquet file is provided to --csv-cells (default: parquet-tools)')
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary (default: gzip)')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--R', type=str, default="Rscript", help='Path to Rscript binary (default: Rscript)')
    
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # Sanity check: JSON mode and manual mode are mutually exclusive
    if args.in_json and args.in_dir:
        parser.error("Cannot enable both JSON mode (--in-json) and manual input mode (--in-dir and/or --csv-cells/--csv-boundaries/--csv-clust/--csv-diffexp). Choose one.")

    return args

def import_xenium_cell(_args):
    """
    Import cell segmentation results from Xenium Ranger output
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.outprefix + "_import_xenium_cell" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    if args.all:
        args.cells = True
        args.boundaries = True
        args.umap = True
    
    if not args.cells and not args.boundaries and not args.update_catalog and not args.umap:
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

    temp_fs=[]

    # read in_json if provided 
    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"File not found: {args.in_json} (--in-json)"
        raw_data=load_file_to_dict(args.in_json)
        cell_data=raw_data.get("CELLS", raw_data) # # use raw_data as default to support the flat dict build in the old scripts
    else:
        cell_data={
            "CELL": f"{args.in_dir}/{args.csv_cells}",
            "BOUNDARY": f"{args.in_dir}/{args.csv_boundaries}",
            "CLUSTER": f"{args.in_dir}/{args.csv_clust}",
            "DE": f"{args.in_dir}/{args.csv_diffexp}",
            "UMAP_PROJ": f"{args.in_dir}/{args.csv_umap}",
            "CELL_FEATURE_MEX": f"{args.in_dir}/{args.mex_dir}",
            # "MEX_BCD": os.path.join(args.in_dir, args.mex_dir, args.mex_bcd),
            # "MEX_FTR": os.path.join(args.in_dir, args.mex_dir, args.mex_ftr),
            # "MEX_MTX": os.path.join(args.in_dir, args.mex_dir, args.mex_mtx),
        }
        if args.pixel is not None:
            cell_data["PIXEL"] = args.pixel

    if cell_data.get("CELL_FEATURE_MEX") is not None:
        mex_ftr_dir = cell_data["CELL_FEATURE_MEX"]
        cell_data["MEX_BCD"] = os.path.join(mex_ftr_dir, args.mex_bcd)
        cell_data["MEX_FTR"] = os.path.join(mex_ftr_dir, args.mex_ftr)
        cell_data["MEX_MTX"] = os.path.join(mex_ftr_dir, args.mex_mtx)

    # Convert pixel-level TSV data into sptsv format
    if not args.skip_redo_pseudobulk or not args.skip_redo_diffexp:
        ## if pixel file is provided, convert to sptsv
        if args.pixel is not None:
            logger.info(f"  * Generating spTSV from pixel TSV file")
            pixelf = cell_data.get("PIXEL", None)
            cmd = f"{args.spatula} pixel2sptsv --pixel {pixelf} --out {args.outprefix}.sptsv --in-col-id {args.pixel_colname_cell_id} --in-col-gene {args.pixel_colname_gene} --in-col-count {args.pixel_colname_count} --gzip {args.gzip} --sort {args.sort}"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            if result.returncode != 0:
                logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                sys.exit(1)
            temp_fs.append(f"{args.outprefix}.sptsv.tsv")
            temp_fs.append(f"{args.outprefix}.sptsv.json")
            temp_fs.append(f"{args.outprefix}.sptsv.feature.counts.tsv")
        else:
            mex_bcd = cell_data.get("MEX_BCD", None)
            mex_ftr = cell_data.get("MEX_FTR", None)
            mex_mtx = cell_data.get("MEX_MTX", None)

            assert mex_bcd is not None and os.path.exists(mex_bcd), (f'Path not provided or file not found: "MEX_BCD" in --in-json' if args.in_json is not None else f'Path not provided or file not found: --mex-bcd')
            assert mex_ftr is not None and os.path.exists(mex_ftr), (f'Path not provided or file not found: "MEX_FTR" in --in-json' if args.in_json is not None else f'Path not provided or file not found: --mex-ftr')
            assert mex_mtx is not None and os.path.exists(mex_mtx), (f'Path not provided or file not found: "MEX_MTX" in --in-json' if args.in_json is not None else f'Path not provided or file not found: --mex-mtx')

            logger.info(f"  * Generating spTSV from MEX files")
            cmd = f"'{args.spatula}' mex2sptsv --bcd '{mex_bcd}' --ftr '{mex_ftr}' --mtx '{mex_mtx}' --out '{args.outprefix}.sptsv'"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            if result.returncode != 0:
                logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                sys.exit(1)
            temp_fs.append(f"{args.outprefix}.sptsv.tsv")
            temp_fs.append(f"{args.outprefix}.sptsv.json")
            temp_fs.append(f"{args.outprefix}.sptsv.feature.counts.tsv")

    # Cluster/DE
    if args.cells or args.boundaries or args.umap:
        ## import cluster
        clust_in = cell_data.get("CLUSTER", None)
        # presence and existence for cluster file
        assert clust_in is not None, ('Path not provided: "CLUSTER" in --in-json' if args.in_json is not None else 'Path not provided: --csv-clust')
        assert os.path.exists(clust_in), (f'File not found: {clust_in} ("CLUSTER" in --in-json)' if args.in_json is not None else f'File not found: {clust_in} (--csv-clust)')

        logger.info(f"Loading cell cluster data from {clust_in}")
        sorted_clusters, cluster2idx, bcd2clusteridx = process_cluster_csv(
            clust_in,
            barcode_col=args.clust_colname_barcode,
            cluster_col=args.clust_colname_cluster,
            output_filename=f"{args.outprefix}-clust.tsv.gz"   
        )
        logger.info(f"  * Loaded {len(bcd2clusteridx)} cells")
        temp_fs.append(f"{args.outprefix}-clust.tsv.gz")

        ## Create pseudobulk DE results from MEX
        if not args.skip_redo_pseudobulk:
            logger.info(f"  * Generating pseudobulk expression matrix for {len(sorted_clusters)} clusters from MEX files")
            pseudobulk_out = f"{args.outprefix}-pseudobulk.tsv.gz"

            cmd = f"'{args.spatula}' sptsv2model --tsv '{args.outprefix}.sptsv.tsv' --clust '{args.outprefix}-clust.tsv.gz' --json '{args.outprefix}.sptsv.json' --out '{pseudobulk_out}'"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            if result.returncode != 0:
                logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                sys.exit(1)

        ## write the color map
        cmap_out=f"{args.outprefix}-rgb.tsv"
        assert os.path.exists(args.tsv_cmap), f"File not found: {args.tsv_cmap} (--tsv-cmap)"        
        
        logger.info(f"  * Writing color map from {args.tsv_cmap} to {cmap_out}")
        write_cmap_tsv(cmap_out, args.tsv_cmap, sorted_clusters)

        if args.cells or args.boundaries:
            if args.skip_redo_diffexp:
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
            else:
                logger.info(f"  * Generating DE results from Pseudobulk files")
                if args.skip_redo_pseudobulk:
                    raise ValueError("Cannot perform MEX DE generation when pseudobulk generation from MEX is skipped (--skip-mex-pseudobulk)")
                ## Generate DE from pseudobulk
                cmd = f"'{args.spatula}' diffexp-model-matrix --tsv1 '{args.outprefix}-pseudobulk.tsv.gz' --out '{args.outprefix}-pseudobulk' --min-fc {args.de_min_fc} --max-pval {args.de_max_pval}"
                result = subprocess.run(cmd, shell=True, capture_output=True)
                if result.returncode != 0:
                    logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                    sys.exit(1)
                pseudobulk_prefix = f"{args.outprefix}-pseudobulk"
                de_out=f"{args.outprefix}-cells-bulk-de.tsv"
                cmd = f"('{args.gzip}' -cd '{pseudobulk_prefix}.de.marginal.tsv.gz' | head -1 | sed 's/^Feature/gene/'; '{args.gzip}' -cd '{pseudobulk_prefix}.de.marginal.tsv.gz' | tail -n +2 | '{args.sort}' -k 2,2n -k 3,3gr;) > '{de_out}'"
                result = subprocess.run(cmd, shell=True, capture_output=True)
                if result.returncode != 0:
                    logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                    sys.exit(1)
                logger.info(f"  * Wrote differential expression results to {de_out}")

                ## create factor reporting files
                ficture2report = os.path.join(args.ficture2, "ext/py/factor_report.py")                        
                assert os.path.exists(ficture2report), f"File not found: {ficture2report}. Please check the path to punkst/ficture2 provided by --ficture2"

                cmd = f"'{args.gzip}' -dc '{pseudobulk_prefix}.tsv.gz' > '{pseudobulk_prefix}.tsv'"
                result = subprocess.run(cmd, shell=True, capture_output=True)
                if result.returncode != 0:
                    logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                    sys.exit(1)                

                cmd = f"head -n $(head -1 '{pseudobulk_prefix}.tsv' | wc -w) '{args.tsv_cmap}' > '{pseudobulk_prefix}.cmap.tsv'"
                result = subprocess.run(cmd, shell=True, capture_output=True)
                if result.returncode != 0:
                    logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                    sys.exit(1)

                cmd = " ".join([
                    args.python,
                    ficture2report,
                    f"--de '{de_out}'",
                    f"--pseudobulk '{pseudobulk_prefix}.tsv'",
                    f"--factor_label factor",
                    f"--feature_label Feature",
                    f"--color_table '{pseudobulk_prefix}.cmap.tsv'",
                    f"--output_pref '{pseudobulk_prefix}'",
                    ])
                result = subprocess.run(cmd, shell=True, capture_output=True)
                if result.returncode != 0:
                    logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                    sys.exit(1)

                cmd = f"cp '{pseudobulk_prefix}.info.tsv' '{args.outprefix}-info.tsv'"
                result = subprocess.run(cmd, shell=True, capture_output=True)
                if result.returncode != 0:
                    logger.error(f"Command {cmd}\nfailed with error: {result.stderr.decode()}")
                    sys.exit(1)
                logger.info(f"  * Wrote factor info to {args.outprefix}-info.tsv")

                temp_fs.append(f"{args.outprefix}-pseudobulk.de.marginal.tsv.gz")
                temp_fs.append(f"{args.outprefix}-pseudobulk.tsv")
                temp_fs.append(f"{args.outprefix}-pseudobulk.cmap.tsv")
                temp_fs.append(f"{args.outprefix}-pseudobulk.info.tsv")

                
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
        process_cells_csv(
            cells_csv,
            cells_out,
            bcd2clusteridx,
            cell_id_col=args.cells_colname_cell_id,
            x_col=args.cells_colname_x,
            y_col=args.cells_colname_y,
            count_col=args.cells_colname_count
        )

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
        col_id = args.boundaries_colname_cell_id
        col_x = args.boundaries_colname_x
        col_y = args.boundaries_colname_y
        with flexopen(bound_in, "rt") as f, flexopen(bound_out, "wt") as wf:
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
                    # cluster = bcd2cluster.get(current_cell_id, "NA")
                    # clusteridx = cluster2idx.get(cluster, "NA") if cluster != "NA" else "NA"
                    clusteridx = bcd2clusteridx.get(current_cell_id, "NA")
                    wf.write(f'{{"type": "Feature", "geometry": {{"type": "Polygon", "coordinates": [[{",".join(current_vertices)}]]}}, "properties": {{"cell_id": "{current_cell_id}", "topK": "{clusteridx}"}}}}\n')
                    current_cell_id = cell_id
                    current_vertices = []
                current_vertices.append(f'[{x},{y}]')
            if current_cell_id is not None:
                # cluster = bcd2cluster.get(current_cell_id, "NA")
                # clusteridx = cluster2idx.get(cluster, "NA") if cluster != "NA" else "NA"
                clusteridx = bcd2clusteridx.get(current_cell_id, "NA")
                wf.write(f'{{"type": "Feature", "geometry": {{"type": "Polygon", "coordinates": [[{",".join(current_vertices)}]]}}, "properties": {{"cell_id": "{current_cell_id}", "topK": "{clusteridx}"}}}}\n')

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
    ## add the new factor to the catalog
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
    new_factor = make_factor_dict(factor_id, factor_name, args.outprefix, factor_type="cells", pmtiles_keys=pmtiles_keys, umap_src=umap_src, pseudobulk_src=args.skip_redo_pseudobulk==False)
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
        new_factor = make_factor_dict(factor_id, factor_name, out_base, factor_type="cells", pmtiles_keys=pmtiles_keys, umap_src=umap_src, pseudobulk_src=args.skip_redo_pseudobulk==False)
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
    if not args.keep_intermediate_files:
        logger.info(f"Cleaning intermediate files")
        if len(temp_fs) >0:
            for temp_f in temp_fs:
                if os.path.exists(temp_f):
                    os.remove(temp_f)


    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

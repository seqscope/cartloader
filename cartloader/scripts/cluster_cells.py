import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, csv, yaml, inspect
import pandas as pd
import numpy as np
from scipy.stats import chi2
import subprocess

from cartloader.utils.utils import scheck_app, create_custom_logger, flexopen, unquote_str, smartsort, write_dict_to_file, load_file_to_dict, run_command

from cartloader.utils.color_helper import normalize_rgb, rgb_to_hex

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Perform LDA-leiden cell clustering")

    run_params = parser.add_argument_group("Run Options", "Run options for GNU Make")
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')
    run_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    run_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Enable all actions: --cells and --boundaries')
    cmd_params.add_argument('--sptsv', action='store_true', default=False, help='Create SPTSV files for LDA clustering')
    cmd_params.add_argument('--lda', action='store_true', default=False, help='Perform LDA factorization')
    cmd_params.add_argument('--leiden', action='store_true', default=False, help='Generate Leiden clusters based on LDA factorization')
    cmd_params.add_argument('--tsne', action='store_true', default=False, help='Generate TSNE manifolds based on LDA factorization')
    cmd_params.add_argument('--umap', action='store_true', default=False, help='Generate UMAP manifolds based on LDA factorization')
    cmd_params.add_argument('--pseudobulk', action='store_true', default=False, help='Generate pseudobulk files based on Leiden clusters')
    cmd_params.add_argument('--heatmap', action='store_true', default=False, help='Generate heamap between LDA factors and Leiden clusters')

    inout_params = parser.add_argument_group("Input output parameters", "Input and output parameters")
    inout_params.add_argument('--mex-dir', type=str, help='Directory containing MEX files')
    inout_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='Barcode files in MEX format')
    inout_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='Feature files in MEX format')
    inout_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='Matrix files in MEX format')
    inout_params.add_argument('--pixel', type=str, help='Pixel-level transcript file TSV format')
    inout_params.add_argument('--out', type=str, help='Output prefix')

    lda_params = parser.add_argument_group("LDA Parameters", "Parameters for LDA factorization")
    lda_params.add_argument('--n-topics', type=int, help='Number of topics for LDA factorization')
    lda_params.add_argument('--seed', type=int, default=42, help='Seed for LDA factorization')
    lda_params.add_argument('--model', type=str, help='LDA model file to use for projection')

    aux_inout_params = parser.add_argument_group("Auxilary Input Parameters that can be specified manually", "Parameters to specify manually")
    aux_inout_params.add_argument('--min-feature-count', type=int, default=1, help='Minimum feature count for LDA factorization')
    aux_inout_params.add_argument('--min-cell-count', type=int, default=1, help='Minimum cell count for LDA factorization')
    aux_inout_params.add_argument('--sptsv-prefix', type=str, help='Prefix for SPTSV files, when --sptsv is skipped')
    aux_inout_params.add_argument('--fit', type=str, help='Results from LDA fit/projection')
    aux_inout_params.add_argument('--clust', type=str, help='Results from Leiden clustering')
    aux_inout_params.add_argument('--in-col-id', type=str, default="cell_id", help='Column name for cell ID in pixel file')
    aux_inout_params.add_argument('--in-col-ftr', type=str, default="gene", help='Column name for feature ID in pixel file')
    aux_inout_params.add_argument('--in-col-cnt', type=str, default="count", help='Column name for count in pixel file')
    aux_inout_params.add_argument('--in-col-x', type=str, default="X", help='Column name for x coordinate in pixel file')
    aux_inout_params.add_argument('--in-col-y', type=str, default="Y", help='Column name for y coordinate in pixel file')
    aux_inout_params.add_argument('--ignore-ids', type=str, default="UNASSIGNED,NA,0,-1", help='IDs to ignore in pixel file')
    aux_inout_params.add_argument('--resolution', type=float, default=1.0, help='Resolution for Leiden clustering')

    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary (default: spatulk)')
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"),  help='Path to punkst (ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--R', type=str, default="Rscript", help='Path to Rscript binary (default: Rscript)')
    
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    return args

def cluster_cells(_args):
    """
    Import cell segmentation results from Xenium Ranger output
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.outprefix + "_import_xenium_cell" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    scheck_app(args.spatula)
    scheck_app(args.R)

    if args.all:
        args.sptsv = True
        args.lda = True
        args.leiden = True
        args.pseudobulk = True
        args.tsne = True
        args.umap = True
        args.heatmap = True
    
    if not args.sptsv and not args.lda and not args.leiden and not args.pseudobulk and not args.tsne and not args.umap and not args.heatmap:
        raise ValueError("At least one action must be enabled.")

    # create output directory if needed
    out_dir = os.path.dirname(args.out)
    out_base = os.path.basename(args.out)
    if not os.path.exists(out_dir) and out_dir != "":
        os.makedirs(out_dir, exist_ok=True)

    # create sptsv
    if args.sptsv:
        ## check if mex or pixel is specified:
        if not args.mex_dir and not args.pixel:
            raise ValueError("At least one of --mex-dir or --pixel-dir must be specified.")
        if args.mex_dir is not None and args.pixel is not None:
            raise ValueError("Cannot specify both --mex-dir and --pixel-dir.") 
        ## perform mex2sptsv
        if args.mex_dir is not None:
            logger.info(f"Creating sptsv from {args.mex_dir}")
            cmd = f"{args.spatula} mex2sptsv --in-dir {args.mex_dir} --bcd {args.mex_bcd} --ftr {args.mex_ftr} --mtx {args.mex_mtx} --out {args.out}.sptsv --min-feature-count {args.min_feature_count}"
        else:
            logger.info(f"Creating sptsv from {args.pixel}")
            cmd = f"{args.spatula} pixel2sptsv --pixel {args.pixel} --in-col-id {args.in_col_id} --in-col-ftr {args.in_col_ftr} --in-col-cnt {args.in_col_cnt} --in-col-x {args.in_col_x} --in-col-y {args.in_col_y} --ignore-ids {args.ignore_ids} --out {args.out}.sptsv"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to create sptsv: {cmd}")
        sptsv_prefix = f"{args.out}.sptsv"

        ## randomize sptsv
        logger.info(f"Randomizing sptsv")
        cmd = f"sort -k 1,1 {sptsv_prefix}.tsv > {sptsv_prefix}.randomized.tsv"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to create sptsv: {cmd}")
    else:
        if args.sptsv_prefix is not None:
            sptsv_prefix = args.sptsv_prefix
        else:
            sptsv_prefix = f"{args.out}.sptsv"

    # run LDA
    if args.lda:
        ficture2bin = os.path.join(args.ficture2, "bin/punkst")
        assert os.path.exists(ficture2bin), f"File not found: {ficture2bin}. FICTURE2 Directory should include bin/punkst (--ficture2)"

        if args.model is None:  ## run LDA to generate model
            if args.n_topics is None:
                raise ValueError("--n-topics must be specified when --model is not specified with --lda ON.")
            logger.info(f"Running LDA")
            cmd = f"{ficture2bin} lda4hex --in-data {sptsv_prefix}.randomized.tsv --in-meta {sptsv_prefix}.json --out-prefix {args.out}.lda --sort-topics --n-topics {args.n_topics} --transform --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            result = subprocess.run(cmd, shell=True, check=True)
            if result.returncode != 0:
                raise ValueError(f"Failed to create sptsv: {cmd}")
        else:  ## use existing model
            logger.info(f"Projecting existing LDA model: {args.model}")
            cmd = f"{ficture2bin} lda4hex --model {args.model} --projection-only --in-data {sptsv_prefix}.randomized.tsv --in-meta {sptsv_prefix}.json --out-prefix {args.out}.lda --transform --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            result = subprocess.run(cmd, shell=True, check=True)
            if result.returncode != 0:
                raise ValueError(f"Failed to create sptsv: {cmd}")

    # run Leiden clustering, TSNE, UMAP generation
    if args.leiden:
        ## perform leiden clustering
        logger.info(f"Running Leiden clustering")
        cmd = f"cartloader lda_leiden_cluster --tsv {args.out}.lda.results.tsv --out {args.out}.leiden.tsv.gz --resolution {args.resolution}"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to run Leiden clustering: {cmd}")

    if args.tsne:
        ## generate TSNE manifolds
        logger.info(f"Generating TSNE manifolds")
        cmd = f"cartloader lda_tsne --tsv {args.out}.lda.results.tsv --out {args.out}.tsne.tsv.gz"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to generate TSNE manifolds: {cmd}")

        ## draw TSNE manifolds
        logger.info(f"Drawing TSNE manifolds")
        draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
        cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold {args.out}.tsne.tsv.gz --tsv-clust {args.out}.leiden.tsv.gz --tsv-colname-x TSNE1 --tsv-colname-y TSNE2 --out {args.out}.tsne.png --tsv-colname-clust cluster"
        print(cmd, file=sys.stderr)
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to draw TSNE manifolds: {cmd}")

    if args.umap:
        ## generate UMAP manifolds
        scheck_app(args.R)
        logger.info(f"Generating UMAP manifolds")
        create_umap_rscript=f"{repo_dir}/cartloader/r/create_umap.r"
        cmd = f"{args.R} '{create_umap_rscript}' --input {args.out}.lda.results.tsv --out-prefix {args.out} --tsv-colname-meta random_key {args.in_col_id}"
        print(cmd, file=sys.stderr)
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to generate UMAP manifolds: {cmd}")

        ## draw UMAP manifolds
        logger.info(f"Drawing UMAP manifolds")
        draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
        cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold {args.out}.umap.tsv.gz --tsv-clust {args.out}.leiden.tsv.gz --tsv-colname-x UMAP1 --tsv-colname-y UMAP2 --out {args.out}.umap.png --tsv-colname-clust cluster"
        print(cmd, file=sys.stderr)
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to draw UMAP manifolds: {cmd}")

    # Produce pseodobulk from leiden clustering results
    if args.pseudobulk:
        logger.info(f"Producing pseudobulk from leiden clustering results")
        cmd = f"{args.spatula} sptsv2model --tsv {sptsv_prefix}.randomized.tsv --clust {args.out}.leiden.tsv.gz --json {sptsv_prefix}.json --out {args.out}.leiden.pseudobulk.tsv.gz"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to produce pseudobulk: {cmd}")

    if args.heatmap:
        ## perform DE test on the model matrix
        logger.info(f"Performing DE test on the model matrix")
        model_tsv = args.model if args.model is not None else f"{args.out}.lda.model.tsv" 
        cmd = f"{args.spatula} diffexp-model-matrix --tsv1 {model_tsv} --out {args.out}.lda.model"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to perform DE test on the model matrix: {cmd}")

        ## perform DE test on the cell pseudobulk matrix
        logger.info(f"Performing DE test on the cell pseudobulk matrix")
        cmd = f"{args.spatula} diffexp-model-matrix --tsv1 {args.out}.leiden.pseudobulk.tsv.gz --out {args.out}.leiden.pseudobulk"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to perform DE test on the cell pseudobulk matrix: {cmd}")

        ## create heatmap
        logger.info(f"Creating heatmap")
        heatmap_rscript=f"{repo_dir}/cartloader/r/create_heatmap.r"
        cmd = f"{args.R} '{heatmap_rscript}' --results {args.out}.lda.results.tsv --clust {args.out}.leiden.tsv.gz --de-results {args.out}.lda.model.de.marginal.tsv.gz --de-clust {args.out}.leiden.pseudobulk.de.marginal.tsv.gz --offset-data 2 --out {args.out}.heatmap --colname-clust cluster --draw"
        result = subprocess.run(cmd, shell=True, check=True)
        if result.returncode != 0:
            raise ValueError(f"Failed to create heatmap: {cmd}")
    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

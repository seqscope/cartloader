import sys, os, gzip, argparse, logging, shutil, subprocess, inspect, json
from venv import logger
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, flexopen, execute_makefile

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Run FICTURE2")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile, and print commands without executing them')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and start from the beginning')
    run_params.add_argument('--threads', type=int, default=8, help='Maximum number of threads per job (default: 8)')
    run_params.add_argument('--n-jobs', type=int, default=2, help='Number of parallel jobs to run (default: 2)')
    run_params.add_argument('--makefn', type=str, default="run_ficture2_multi_cells.mk", help='File name of Makefile to write (default: run_ficture2_multi.mk)')

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Enable all actions: --cells and --boundaries')
    cmd_params.add_argument('--sptsv', action='store_true', default=False, help='Create SPTSV files for LDA clustering')
    cmd_params.add_argument('--lda', action='store_true', default=False, help='Perform LDA factorization')
    cmd_params.add_argument('--leiden', action='store_true', default=False, help='Generate Leiden clusters based on LDA factorization')
    cmd_params.add_argument('--tsne', action='store_true', default=False, help='Generate TSNE manifolds based on LDA factorization')
    cmd_params.add_argument('--umap', action='store_true', default=False, help='Generate UMAP manifolds based on LDA factorization')
    cmd_params.add_argument('--pseudobulk', action='store_true', default=False, help='Generate pseudobulk files based on Leiden clusters')
    cmd_params.add_argument('--heatmap', action='store_true', default=False, help='Generate heamap between LDA factors and Leiden clusters')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='Perform pixel-level decoding based on cell clusters. Onlt available when --pixel is specified.')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output parameters for FICTURE")
    inout_params.add_argument('--out-dir', required=True, type=str, help='Output directory')
    inout_params.add_argument('--out-prefix', type=str, default="cells", help='Prefix for output files (default: cells)')
    inout_params.add_argument('--out-json', type=str, default=None, help="Path to output JSON file to store analysis parameters (default: <out-dir>/ficture.cells.params.json)")
    inout_params.add_argument('--in-dir', type=str, default=None, help='Path to input directory containing FICTURE output files. Same to out_dit if not specified')
    inout_params.add_argument('--mex-dir', type=str, help='Directory containing MEX files')
    inout_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='Barcode files in MEX format')
    inout_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='Feature files in MEX format')
    inout_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='Matrix files in MEX format')
    inout_params.add_argument('--mex-list', type=str, help='TSV file containing sample IDs and paths to MEX files')
    inout_params.add_argument('--sptsv-prefix', type=str, help='Prefix for SPTSV files')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--n-factor', type=int, help='Number of factors for LDA training.')
    key_params.add_argument('--resolution', type=float, default=1.0, help='Resolution for Leiden clustering (default: 1.0)')
    key_params.add_argument('--anchor-res', type=int, default=6, help='Anchor resolution for decoding (default: 6)')
    key_params.add_argument('--cmap-file', type=str, default=os.path.join(repo_dir, "assets", "fixed_color_map_256.tsv"), help='Path to fixed color map TSV (default: <cartloader_dir>/assets/fixed_color_map_256.tsv)')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--colidx-feature',  type=int, default=3, help='Column index for Y-axis in the --in-transcript (default: 3)')
    aux_params.add_argument('--colidx-count',  type=int, default=4, help='Column index for intensity in the --in-transcript (default: 4)')
    aux_params.add_argument('--colidx-cell-id', type=int, default=5, help='Column index for cell ID in the --in-transcript (default: 5)')
    aux_params.add_argument('--ignore-ids', type=str, default="UNASSIGNED,NA,0,-1", help='IDs to ignore in pixel file')
    # segmentation - ficture
    aux_params.add_argument('--min-count-per-sample', type=int, default=50, help='Minimum count per sample in the tiled SGE (default: 50)')
    aux_params.add_argument('--min-ct-per-cell', type=int, default=50, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format (default: 50)')
    # train
    aux_params.add_argument('--train-epoch', type=int, default=2, help='Training epoch for LDA model (default: 2)')
    aux_params.add_argument('--skip-umap', action='store_true', default=False, help='Skip creating umap')
    aux_params.add_argument('--decode-scale', type=int, default=1, help='Decode scale (default: 1)')
    aux_params.add_argument('--seed', type=int, default=1, help='Random seed for random number generation (default: 0)')

    # others parameters shared across steps
    aux_params.add_argument('--min-feature-count', type=int, default=20, help='Minimum feature count for LDA factorization')
    aux_params.add_argument('--min-cell-count', type=int, default=50, help='Minimum cell count for LDA factorization')
    aux_params.add_argument('--de-min-ct-per-feature', type=int, default=20, help='Minimum count per feature for differential expression (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='P-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')
    aux_params.add_argument('--decode-fit-width', type=int, default=18, help='Fitting width (in microns) for decoding (default: 10)')
    # project from external model
    aux_params.add_argument('--pretrained-model', type=str, help='Path to a pre-trained model to use for projection. If provided, LDA training will be skipped, and the provided model will be used for projection.')
    aux_params.add_argument('--list-cluster', type=str, help='Path to a existing cluster files to create pseudobulk matrix in the format of [SAMPLE_ID] [CLUSTER_FILE]. If provided, Leiden clustering will be skipped, and the provided cluster files will be used for pseudobulk generation.')
    aux_params.add_argument('--list-xy', type=str, help='Path to a existing file containing X/Y locations of each cell. If provided, X/Y locations will be read from the provided file instead of computing from pixel file.')
    aux_params.add_argument('--list-boundaries', type=str, help='Path to a existing file containing cell boundaries. This file will simply be stored in the output JSON for future use.')
    aux_params.add_argument('--xy-colname-cell-id', type=str, default="cell_id", help='Column name for cell IDs in the metadata file (default: cell_id)')
    aux_params.add_argument('--xy-colname-x', type=str, default="X", help='Column name for X coordinates in the metadata file (default: X)')
    aux_params.add_argument('--xy-colname-y', type=str, default="Y", help='Column name for Y coordinates in the metadata file (default: Y)')

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Customizing Auxiliary Parameters", "Customize features (typically genes) used by FICTURE without altering the original feature TSV") # This ensures the original feature TSV file is retained in the output JSON file for downstream processing 
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='Regex of feature names to include')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default="1G", help='Memory size for each process (default: 1G)')
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary (default: "spatula" in the system PATH)') # default=f"{repo_dir}/submodules/spatula/bin/spatula",
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"), help='Path to punkst (ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--python', type=str, default="python3",  help='Python3 binary')
    env_params.add_argument('--R', type=str, default="Rscript", help='Path to R binary for UMAP generation (default: R)')


    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_ficture2_multi_cells(_args):
    """Run all functions in FICTURE2 cell clustering with multi-sample pipeline
    This function is meant to be used in a local environment that has sufficient resources to run all functions in FICTURE at once.
    This function performs the following tasks:
    (1) Take the input parameters relevant to the FICTURE runs
    (2) Identify the sequence of commands to run FICTURE
    (3) Create a GNU makefile to run the commands in parallel
    (4) Run the GNU makefile
    """
    # args
    args=parse_arguments(_args)

    # input/output/other files
    # dirs
    os.makedirs(args.out_dir, exist_ok=True)

    # ficture2
    ficture2bin = os.path.join(args.ficture2, "bin/punkst")
    assert os.path.exists(ficture2bin), f"File not found: {ficture2bin}. FICTURE2 Directory should include bin/punkst (--ficture2)"
    
    ficture2report = args.python + " " + os.path.join(args.ficture2, "ext/py/factor_report.py")
    
    # out files
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.{args.out_prefix}.params.json")

    if args.in_dir is None:
        args.in_dir = args.out_dir

    ## parse the input list file
    in_samples = []
    ## list directories in args.in_dir/samples/
    samples_dir = os.path.join(args.in_dir, "samples")
    for entry in os.listdir(samples_dir):
        entry_path = os.path.join(samples_dir, entry)
        if os.path.isdir(entry_path):
            in_samples.append(entry)
    logger.info(f"Found {len(in_samples)} samples in {samples_dir}: {in_samples}")
#   in_tsvs = []
#    with flexopen(args.in_list, "rt") as f:
#        for line in f:
#            toks = line.strip().split("\t")
#            in_samples.append(toks[0])
#           in_tsvs.append(toks[1])
    
    # cmap
    assert os.path.exists(args.cmap_file), f"File not found: {args.cmap_file} (--cmap-file)"
    
    # start mm
    mm = minimake()

    # assume that multi-sample tiling already exists
    # check the existence of tiled files
    for sample_id in in_samples:
        sample_prefix = f"{args.in_dir}/samples/{sample_id}/{sample_id}.tiled"
        if not os.path.exists(f"{sample_prefix}.tsv"):
            raise FileNotFoundError(f"File not found: {sample_prefix}.tsv. Please run run_ficture2_multi first.")
        if not os.path.exists(f"{sample_prefix}.index"):
            raise FileNotFoundError(f"File not found: {sample_prefix}.index. Please run run_ficture2_multi first.")

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
        args.decode = True

    ## create cell-based SPTSV files
    if args.sptsv:
        if args.sptsv_prefix is not None:
            raise ValueError("When --sptsv is ON, --sptsv-prefix should not be provided.")
        sptsv_prefix = os.path.join(args.out_dir, args.out_prefix) + ".sptsv"
        cmds = cmd_separator([], f"Creating cell-based SPTSV files...")
        samp2sptsv = {} ## sample ID to SPTSV file mapping
        ## if mex_dir is provided, create SPTSV from MEX files, assuming that MEX files contain cell-level data across all samples
        if args.mex_dir is not None:
            if args.mex_list is not None:
                raise ValueError("When --mex-dir is provided, --mex-list should not be provided.")
            cmd = f"{args.spatula} mex2sptsv --in-dir {args.mex_dir} --bcd {args.mex_bcd} --ftr {args.mex_ftr} --mtx {args.mex_mtx} --out {sptsv_prefix} --min-feature-count {args.min_feature_count}"
            cmds.append(cmd)
            ## in this case, we assume that merged SPTSV file already exists
        elif args.mex_list is not None:
            with flexopen(args.mex_list, 'rt') as rf:
                for line in rf:
                    toks = line.strip().split("\t")
                    sample_id = toks[0]
                    if len(toks) == 2:
                        mex_dir = toks[1]
                        mex_bcd = os.path.join(mex_dir, args.mex_bcd)
                        mex_ftr = os.path.join(mex_dir, args.mex_ftr)
                        mex_mtx = os.path.join(mex_dir, args.mex_mtx)
                    elif len(toks) == 4:
                        mex_bcd = toks[1]
                        mex_ftr = toks[2]
                        mex_mtx = toks[3]
                    else:
                        raise ValueError(f"Each line in --mex-list must have 2 or 4 columns. Found {len(toks)} columns in line: {line}")
                    sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
                    cmd = f"{args.spatula} mex2sptsv --bcd {mex_bcd} --ftr {mex_ftr} --mtx {mex_mtx} --out {sample_sptsv_prefix} --min-feature-count {args.min_feature_count}"
                    cmds.append(cmd)
                    samp2sptsv[sample_id] = sample_sptsv_prefix
        else:
            for sample_id in in_samples:
                #sample_id = in_samples[0]  ## use the first sample's tiled file to create SPTSV
                pixelf = f"{args.in_dir}/samples/{sample_id}/{sample_id}.tiled"
                sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
                cmd = f"{args.spatula} pixel2sptsv --pixel {pixelf}.tsv --no-header --idx-col-x {args.colidx_x} --idx-col-y {args.colidx_y} --idx-col-ftr {args.colidx_feature} --idx-col-cnt {args.colidx_count} --idx-col-id {args.colidx_cell_id} --ignore-ids {args.ignore_ids} --out {sample_sptsv_prefix} --min-feature-count {args.min_feature_count}"
                cmds.append(cmd)
                samp2sptsv[sample_id] = sample_sptsv_prefix
        
        ## merge SPTSV files if needed
        if len(samp2sptsv) > 0:
            ## create a list file
            samp_listf = f"{sptsv_prefix}.list.tsv"
            with flexopen(samp_listf, "wt") as wf:
                for sample_id in samp2sptsv:
                    sample_sptsv_prefix = samp2sptsv[sample_id]
                    #wf.write(f"{sample_id}\t{sample_sptsv_prefix}.feature.counts.tsv\t{sample_sptsv_prefix}.tsv\t{sample_sptsv_prefix}.json\t-2\n")
                    wf.write(f"{sample_id}\t{sample_sptsv_prefix}.feature.counts.tsv\t{sample_sptsv_prefix}.tsv\t{sample_sptsv_prefix}.json\n")
            #cmd = f"{ficture2bin} merge-units --in-list {samp_listf} --out-pref {sptsv_prefix} --temp-dir {sptsv_prefix}.tmp --threads {args.threads}"
            cmd = f"{args.spatula} merge-sptsv --list {samp_listf} --out {sptsv_prefix}"
            cmds.append(cmd)
        ## randomize SPTSV file
        cmd = f"sort -k 1,1 {sptsv_prefix}.tsv > {sptsv_prefix}.randomized.tsv"
        cmds.append(cmd)
        cmds.append(f"[ -f {sptsv_prefix}.randomized.tsv ] && touch {sptsv_prefix}.done" )
        mm.add_target(f"{sptsv_prefix}.done", [f"{args.mex_dir}/{args.mex_bcd}", f"{args.mex_dir}/{args.mex_ftr}", f"{args.mex_dir}/{args.mex_mtx}"] if args.mex_dir is not None else [], cmds)
    elif args.sptsv_prefix is not None:
        sptsv_prefix = args.sptsv_prefix
    else:
        raise ValueError("When --sptsv is not provided, --sptsv-prefix must be provided.")

    ## perform LDA training
    if args.lda:
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        if args.pretrained_model is None:  ## run LDA to generate model
            cmds = cmd_separator([], f"Performing LDA training/projection...")
            if args.n_factor is None:
                raise ValueError("--n-factor must be specified when --model is not specified with --lda ON.")
            cmd = f"{ficture2bin} lda4hex --in-data {sptsv_prefix}.randomized.tsv --in-meta {sptsv_prefix}.json --out-prefix {lda_prefix} --sort-topics --n-topics {args.n_factor} --transform --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            cmds.append(cmd)
            cmds.append(f"[ -f {lda_prefix}.model.tsv ] && [ -f {lda_prefix}.results.tsv ] && touch {lda_prefix}.multi.done" )
            mm.add_target(f"{lda_prefix}.multi.done", [f"{sptsv_prefix}.done"], cmds)
        else:  ## use existing model
            cmds = cmd_separator([], f"Projecting existing LDA model...")
            cmd = f"{ficture2bin} lda4hex --model-prior {args.pretrained_model} --projection-only --in-data {sptsv_prefix}.randomized.tsv --in-meta {sptsv_prefix}.json --out-prefix {lda_prefix} --transform --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            cmds.append(cmd)
            cmds.append(f"[ -f '{lda_prefix}.results.tsv' ] && touch '{lda_prefix}.multi.done'" )
            mm.add_target(f"{lda_prefix}.multi.done", [f"{sptsv_prefix}.done"], cmds)
        ## project the LDA model to each sample separately
        deps = [f"{lda_prefix}.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Performing LDA projection for {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
            cmd = f"{ficture2bin} lda4hex --model-prior {lda_prefix}.model.tsv --projection-only --in-data {sample_sptsv_prefix}.tsv --in-meta {sample_sptsv_prefix}.json --out-prefix {sample_lda_prefix} --transform --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            cmds.append(cmd)
            cmds.append(f"[ -f {sample_lda_prefix}.results.tsv ] && touch {sample_lda_prefix}.done" )
            mm.add_target(f"{sample_lda_prefix}.done", [f"{lda_prefix}.multi.done"], cmds)
            deps.append(f"{sample_lda_prefix}.done")
        ## final target
        cmds = cmd_separator([], f"Finalizing LDA projection for all samples...")
        cmds.append(f"touch {lda_prefix}.done")
        mm.add_target(f"{lda_prefix}.done", deps, cmds);

    if args.leiden:
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        cmds = cmd_separator([], f"Generating Leiden clusters...")
        if args.list_cluster is None:
            cmd = f"cartloader lda_leiden_cluster --offset-data 4 --tsv '{lda_prefix}.results.tsv' --out '{leiden_prefix}.tsv.gz' --resolution {args.resolution} --colname-cluster topK --key-ids sample_id cell_id"
            cmds.append(cmd)
            for sample_id in in_samples:
                sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
                sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
                sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
                cmd = f"cartloader lda_leiden_cluster --offset-data 3 --tsv '{sample_lda_prefix}.results.tsv' --out '{sample_leiden_prefix}.tsv.gz' --resolution {args.resolution} --colname-cluster topK --key-ids cell_id"
                cmds.append(cmd)
        else:
            samp2clust = {}
            with flexopen(args.list_cluster, "rt") as rf:
                for line in rf:
                    toks = line.strip().split("\t")
                    if len(toks) != 2:
                        raise ValueError(f"Each line in --list-cluster must have exactly 3 columns containing [SAMPLE_ID] [CLUSTER_FILE] [METADTA_FILE]")
                    sample_id = toks[0]
                    cluster_file = toks[1]
                    if not os.path.exists(cluster_file):
                        raise FileNotFoundError(f"File not found: {cluster_file} (from --list-cluster)")
                    samp2clust[sample_id] = cluster_file
            with flexopen(f"{leiden_prefix}.tsv.gz", "wt") as wf:
                wf.write("sample_id\tcell_id\ttopK\n")
                for sample_id in in_samples: 
                    sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
                    clustf = samp2clust[sample_id]
                    logger.info(f"Reformatting existing cluster file {clustf} for sample {sample_id}...")
                    with flexopen(clustf, "rt") as rf, flexopen(f"{sample_leiden_prefix}.tsv.gz", "wt") as wf_sample:
                        delim = None
                        nlines = 0
                        wf_sample.write("cell_id\ttopK\n")
                        for line in rf:
                            if delim is None:
                                if line.find("\t") != -1:
                                    delim = "\t"
                                elif line.find(",") != -1:
                                    delim = ","
                                elif line.find(" ") != -1:
                                    delim = " "
                                else:
                                    raise ValueError(f"Cannot determine delimiter in existing cluster file based on the first line {line}.")
                            toks = line.strip().split(delim)
                            if len(toks) != 2:
                                raise ValueError(f"Each line in existing cluster file must have exactly 2 columns. Found {len(toks)} columns in line: {line}")
                            cell_id = toks[0].replace('"', '')
                            cluster_id = toks[1].replace('"', '')
                            if nlines > 0 or cluster_id.isdigit():
                                int_cluster_id = int(cluster_id)-1 ## convert to 0-based
                                wf.write(f"{sample_id}\t{cell_id}\t{int_cluster_id}\n")
                                wf_sample.write(f"{cell_id}\t{int_cluster_id}\n")
                            nlines += 1
            ## write per-sample metadata file

        ## spatial visualization of leiden clusters
        samp2xy = {}
        if args.list_xy is not None:
            with flexopen(args.list_xy, "rt") as rf:
                for line in rf:
                    toks = line.strip().split("\t")
                    if len(toks) != 2:
                        raise ValueError(f"Each line in --list-xy must have exactly 2 columns containing [SAMPLE_ID] [XY_FILE]")
                    sample_id = toks[0]
                    xy_file = toks[1]
                    if not os.path.exists(xy_file):
                        raise FileNotFoundError(f"File not found: {xy_file} (from --list-xy)")
                    samp2xy[sample_id] = xy_file
        for sample_id in in_samples:
            metaf = f"{sample_sptsv_prefix}.cell.metadata.tsv"
            if sample_id in samp2xy:
                xyf = samp2xy[sample_id]
                sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
                metaf = f"{sample_sptsv_prefix}.cell.xy.tsv"
                with flexopen(xyf, "rt") as rf, flexopen(metaf, "wt") as wf_sample:
                    delim = None
                    header = None
                    nlines = 0
                    for line in rf:
                        if delim is None:
                            if line.find("\t") != -1:
                                delim = "\t"
                            elif line.find(",") != -1:
                                delim = ","
                            elif line.find(" ") != -1:
                                delim = " "
                            else:
                                raise ValueError(f"Cannot determine delimiter in existing metadata file based on the first line {line}.")
                        toks = line.strip().split(delim)
                        if nlines == 0:
                            col2idx = {colname.replace('"', ''): idx for idx, colname in enumerate(toks)}
                            idx_cell_id = col2idx.get(args.xy_colname_cell_id, None)
                            idx_x = col2idx.get(args.xy_colname_x, None)
                            idx_y = col2idx.get(args.xy_colname_y, None)
                            if idx_cell_id is None or idx_x is None or idx_y is None:
                                raise ValueError(f"Column names {args.xy_colname_cell_id}, {args.xy_colname_x}, {args.xy_colname_y} not found in existing metadata file {metaf}.")
                            wf_sample.write("cell_id\tX\tY\n")
                        else:
                            cell_id = toks[idx_cell_id].replace('"', '')
                            x = toks[idx_x]
                            y = toks[idx_y]
                            wf_sample.write(f"{cell_id}\t{x}\t{y}\n")
                        nlines += 1
            draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
            cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{metaf}' --tsv-clust '{sample_leiden_prefix}.tsv.gz' --tsv-colname-x X --tsv-colname-y Y --out '{sample_leiden_prefix}.xy.png' --out-tsv '{sample_leiden_prefix}.xy.tsv.gz' --tsv-colname-clust topK"
            cmds.append(cmd)
        cmds.append(f"[ -f '{leiden_prefix}.tsv.gz' ] && touch '{leiden_prefix}.done'" )
        mm.add_target(f"{leiden_prefix}.done", [f"{lda_prefix}.done"], cmds)

        ## spatial visualization of leiden clusters
        samp2boundaries = {}
        if args.list_boundaries is not None:
            with flexopen(args.list_boundaries, "rt") as rf:
                for line in rf:
                    toks = line.strip().split("\t")
                    if len(toks) != 2:
                        raise ValueError(f"Each line in --list-boundaries must have exactly 2 columns containing [SAMPLE_ID] [BOUNDARIES_FILE]")
                    sample_id = toks[0]
                    boundaries_file = toks[1]
                    if not os.path.exists(boundaries_file):
                        raise FileNotFoundError(f"File not found: {boundaries_file} (from --list-boundaries)")
                    samp2boundaries[sample_id] = boundaries_file
    if args.tsne:
        ## generate TSNE manifolds
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        tsne_prefix = os.path.join(args.out_dir, args.out_prefix) 
        cmds = cmd_separator([], f"Generating TSNE manifolds...")
        cmd = f"cartloader lda_tsne --offset-data 4 --tsv '{lda_prefix}.results.tsv' --out '{tsne_prefix}.tsne.tsv.gz'"
        cmds.append(cmd)

        draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
        cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{tsne_prefix}.tsne.tsv.gz' --tsv-clust '{leiden_prefix}.tsv.gz' --tsv-colname-x TSNE1 --tsv-colname-y TSNE2 --out '{tsne_prefix}.tsne.png' --out-tsv '{tsne_prefix}.tsne.leiden.tsv.gz' --tsv-colname-clust topK"
        cmds.append(cmd)
        cmds.append(f"[ -f '{tsne_prefix}.tsne.leiden.tsv.gz' ] && [ -f '{tsne_prefix}.tsne.png' ] && touch '{tsne_prefix}.tsne.multi.done'" )
        mm.add_target(f"{tsne_prefix}.tsne.multi.done", [f"{lda_prefix}.done", f"{leiden_prefix}.done"], cmds)

        deps = [f"{tsne_prefix}.tsne.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Generating TSNE manifolds for sample {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_tsne_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}"
            cmd = f"cartloader lda_tsne --offset-data 3 --tsv '{sample_lda_prefix}.results.tsv' --out '{sample_tsne_prefix}.tsne.tsv.gz'"
            cmds.append(cmd)

            draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
            cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{sample_tsne_prefix}.tsne.tsv.gz' --tsv-clust '{sample_leiden_prefix}.tsv.gz' --tsv-colname-x TSNE1 --tsv-colname-y TSNE2 --out '{sample_tsne_prefix}.tsne.png' --out-tsv '{sample_tsne_prefix}.tsne.leiden.tsv.gz' --tsv-colname-clust topK"
            cmds.append(cmd)
            cmds.append(f"[ -f '{sample_tsne_prefix}.tsne.leiden.tsv.gz' ] && [ -f '{sample_tsne_prefix}.tsne.png' ] && touch '{sample_tsne_prefix}.tsne.done'" )
            mm.add_target(f"{sample_tsne_prefix}.tsne.done", [f"{sample_lda_prefix}.done", f"{leiden_prefix}.done"], cmds)
            deps.append(f"{sample_tsne_prefix}.tsne.done")
        ## final target 
        cmds = cmd_separator([], f"Finalizing TSNE manifolds for all samples...")
        cmds.append(f"touch {tsne_prefix}.tsne.done")
        mm.add_target(f"{tsne_prefix}.tsne.done", deps, cmds);

    if args.umap:
        scheck_app(args.R)
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        umap_prefix = os.path.join(args.out_dir, args.out_prefix) 
        cmds = cmd_separator([], f"Generating UMAP manifolds...")
        create_umap_rscript=f"{repo_dir}/cartloader/r/create_umap.r"
        cmd = f"{args.R} '{create_umap_rscript}' --input '{lda_prefix}.results.tsv' --out-prefix '{umap_prefix}' --tsv-colname-meta random_key sample_id cell_id"
        cmds.append(cmd)

        ## draw UMAP manifolds
        draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
        cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{umap_prefix}.umap.tsv.gz' --tsv-clust '{leiden_prefix}.tsv.gz' --tsv-colname-x UMAP1 --tsv-colname-y UMAP2 --out '{umap_prefix}.umap.png' --out-tsv '{umap_prefix}.umap.leiden.tsv.gz' --tsv-colname-clust topK"
        cmds.append(cmd)
        cmds.append(f"[ -f '{umap_prefix}.umap.leiden.tsv.gz' ] && [ -f '{umap_prefix}.umap.png' ] && touch '{umap_prefix}.umap.multi.done'" )
        mm.add_target(f"{umap_prefix}.umap.multi.done", [f"{lda_prefix}.done", f"{leiden_prefix}.done"], cmds)

        deps = [f"{umap_prefix}.umap.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Generating UMAP manifolds for sample {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_umap_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}"
            cmd = f"{args.R} '{create_umap_rscript}' --input '{sample_lda_prefix}.results.tsv' --out-prefix '{sample_umap_prefix}' --tsv-colname-meta random_key cell_id"
            cmds.append(cmd)

            draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
            cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{sample_umap_prefix}.umap.tsv.gz' --tsv-clust '{sample_leiden_prefix}.tsv.gz' --tsv-colname-x UMAP1 --tsv-colname-y UMAP2 --out '{sample_umap_prefix}.umap.png' --out-tsv '{sample_umap_prefix}.umap.leiden.tsv.gz' --tsv-colname-clust topK"
            cmds.append(cmd)
            cmds.append(f"[ -f '{sample_umap_prefix}.umap.leiden.tsv.gz' ] && [ -f '{sample_umap_prefix}.umap.png' ] && touch '{sample_umap_prefix}.umap.done'" )
            mm.add_target(f"{sample_umap_prefix}.umap.done", [f"{sample_lda_prefix}.done", f"{leiden_prefix}.done"], cmds)
            deps.append(f"{sample_umap_prefix}.umap.done")
        ## final target
        cmds = cmd_separator([], f"Finalizing UMAP manifolds for all samples...")
        cmds.append(f"touch {umap_prefix}.umap.done")
        mm.add_target(f"{umap_prefix}.umap.done", deps, cmds);

    if args.pseudobulk:
        sptsv_prefix = os.path.join(args.out_dir, args.out_prefix) + ".sptsv"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        pseudobulk_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden.pseudobulk"

        cmds = cmd_separator([], f"Generating pseudobulk matrix...")
        cmd = f"{args.spatula} sptsv2model --min-count {args.min_feature_count} --tsv '{sptsv_prefix}.randomized.tsv' --clust '{leiden_prefix}.tsv.gz' --features '{sptsv_prefix}.feature.counts.tsv' --json '{sptsv_prefix}.json' --out '{pseudobulk_prefix}.tsv'"
        cmds.append(cmd)
        cmds.append(f"[ -f '{pseudobulk_prefix}.tsv' ] && touch '{pseudobulk_prefix}.multi.done'" )
        mm.add_target(f"{pseudobulk_prefix}.multi.done", [f"{sptsv_prefix}.done", f"{leiden_prefix}.done"], cmds)

        deps = [f"{pseudobulk_prefix}.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Generating pseudobulk matrix for sample {sample_id}...")
            sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_pseudobulk_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden.pseudobulk"
            cmd = f"{args.spatula} sptsv2model --min-count {args.min_feature_count} --tsv '{sample_sptsv_prefix}.tsv' --clust '{sample_leiden_prefix}.tsv.gz' --features '{sample_sptsv_prefix}.feature.counts.tsv' --json '{sample_sptsv_prefix}.json' --out '{sample_pseudobulk_prefix}.tsv'"
            cmds.append(cmd)

            cmd = f"head -n $(head -1 '{sample_pseudobulk_prefix}.tsv' | wc -w) '{args.cmap_file}' > '{sample_pseudobulk_prefix}.cmap.tsv'"
            cmds.append(cmd)

            cmds.append(f"[ -f '{sample_pseudobulk_prefix}.tsv' ] && touch '{sample_pseudobulk_prefix}.done'" )
            mm.add_target(f"{sample_pseudobulk_prefix}.done", [f"{sptsv_prefix}.done", f"{leiden_prefix}.done"], cmds)
            deps.append(f"{sample_pseudobulk_prefix}.done")
        ## final target
        cmds = cmd_separator([], f"Finalizing pseudobulk matrix for all samples...")
        cmds.append(f"touch {pseudobulk_prefix}.done")
        mm.add_target(f"{pseudobulk_prefix}.done", deps, cmds);

    if args.heatmap:
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        pseudobulk_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden.pseudobulk"
        heatmap_prefix = os.path.join(args.out_dir, args.out_prefix) + ".heatmap"

        cmds = cmd_separator([], f"Generating heatmap between LDA factors and Leiden clusters...")
        model_tsv = args.pretrained_model if args.pretrained_model is not None else f"{lda_prefix}.model.tsv"
        cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{model_tsv}' --out '{lda_prefix}.model'"
        cmds.append(cmd)
        cmds.append(f"({args.gzip} -cd {lda_prefix}.model.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {lda_prefix}.model.de.marginal.tsv.gz | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > {lda_prefix}.model.de.tsv")
        cmds.append(f"rm -f '{lda_prefix}.model.de.marginal.tsv.gz'")

        ## perform DE test on the cell pseudobulk matrix
        cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{pseudobulk_prefix}.tsv' --out '{pseudobulk_prefix}'"
        cmds.append(cmd)
        cmds.append(f"({args.gzip} -cd {pseudobulk_prefix}.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {pseudobulk_prefix}.de.marginal.tsv.gz | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > {pseudobulk_prefix}.de.tsv")
        cmds.append(f"rm -f '{pseudobulk_prefix}.de.marginal.tsv.gz'")

        ## create heatmap
        heatmap_rscript=f"{repo_dir}/cartloader/r/create_heatmap.r"
        cmd = f"{args.R} '{heatmap_rscript}' --results '{lda_prefix}.results.tsv' --clust '{leiden_prefix}.tsv.gz' --de-results '{lda_prefix}.model.de.marginal.tsv.gz' --de-clust '{pseudobulk_prefix}.de.marginal.tsv.gz' --offset-data 3 --out '{heatmap_prefix}' --colname-clust topK --cell-clust sample_id cell_id --draw"
        cmds.append(cmd)
        cmds.append(f"[ -f '{heatmap_prefix}.pdf' ] && touch '{heatmap_prefix}.multi.done'" )
        mm.add_target(f"{heatmap_prefix}.multi.done", [f"{lda_prefix}.done", f"{leiden_prefix}.done", f"{pseudobulk_prefix}.done"], cmds)
        deps = [f"{heatmap_prefix}.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Generating heatmap between LDA factors and Leiden clusters for sample {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_pseudobulk_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden.pseudobulk"
            sample_heatmap_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.heatmap"
            model_tsv = args.pretrained_model if args.pretrained_model is not None else f"{sample_lda_prefix}.model.tsv"
            #cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{model_tsv}' --out '{sample_lda_prefix}.model'"
            #cmds.append(cmd)

            ## perform DE test on the cell pseudobulk matrix
            cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{sample_pseudobulk_prefix}.tsv' --out '{sample_pseudobulk_prefix}'"
            cmds.append(cmd)

            cmds.append(f"({args.gzip} -cd {sample_pseudobulk_prefix}.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {sample_pseudobulk_prefix}.de.marginal.tsv.gz | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > {sample_pseudobulk_prefix}.de.tsv")
            cmds.append(f"rm -f '{sample_pseudobulk_prefix}.de.marginal.tsv.gz'")

            ## create heatmap
            heatmap_rscript=f"{repo_dir}/cartloader/r/create_heatmap.r"
            cmd = f"{args.R} '{heatmap_rscript}' --results '{sample_lda_prefix}.results.tsv' --clust '{sample_leiden_prefix}.tsv.gz' --de-results '{lda_prefix}.model.de.tsv' --de-clust '{sample_pseudobulk_prefix}.de.tsv' --offset-data 2 --out '{sample_heatmap_prefix}' --colname-clust topK --cell-clust cell_id --draw"
            cmds.append(cmd)
            cmds.append(f"[ -f '{sample_heatmap_prefix}.pdf' ] && touch '{sample_heatmap_prefix}.done'" )
            mm.add_target(f"{sample_heatmap_prefix}.done", [f"{sample_lda_prefix}.done", f"{leiden_prefix}.done", f"{sample_pseudobulk_prefix}.done", f"{heatmap_prefix}.multi.done"], cmds)
            deps.append(f"{sample_heatmap_prefix}.done")
        ## final target
        cmds = cmd_separator([], f"Finalizing heatmap for all samples...")
        cmds.append(f"touch {heatmap_prefix}.done")
        mm.add_target(f"{heatmap_prefix}.done", deps, cmds);

    if args.decode:
        pseudobulk_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden.pseudobulk"
        decode_all_prefix = os.path.join(args.out_dir, args.out_prefix) + ".decode"
        modelf = f"{pseudobulk_prefix}.tsv"
        modelf_done = f"{pseudobulk_prefix}.done"
        sample_decode_done_files = []
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Performing pixel-level decoding for sample {sample_id}...")
            sample_prefix = f"{args.in_dir}/samples/{sample_id}/{sample_id}.tiled"
            decode_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.pixel"
            fit_width = args.decode_fit_width  ## e.g., 10um
            fit_n_move = fit_width // 2 + 1
            decode_id = f"p{fit_width}_a{args.anchor_res}"
            model_path= modelf
            cmd = " ".join([
                ficture2bin, "pixel-decode",
                f"--model '{model_path}'",
                f"--in-tsv '{sample_prefix}.tsv'",
                f"--in-index '{sample_prefix}.index'",
                f"--temp-dir '{decode_prefix}.tmp'",
                f"--out '{decode_prefix}.tsv'",
                f"--icol-x {args.colidx_x-1}",
                f"--icol-y {args.colidx_y-1}",
                f"--icol-feature 2",
                f"--icol-val 3",
                f"--hex-grid-dist {fit_width}",
                f"--n-moves {fit_n_move}",
                f"--pixel-res 0.5",
                f"--threads {args.threads}",
                f"--seed {args.seed}",
                f"--output-original"
                ])
            cmds.append(cmd)
            cmd = " ".join([
                ficture2bin, "draw-pixel-factors",
                f"--in-tsv '{decode_prefix}.tsv'",
                f"--header-json '{decode_prefix}.json'",
                f"--in-color '{args.cmap_file}'",
                f"--out '{decode_prefix}.png'",
                f"--scale {args.decode_scale}",
                f"--range '{sample_prefix}.coord_range.tsv'"
                ])
            cmds.append(cmd)
            cmd = f"gzip -f '{decode_prefix}.tsv'"
            cmds.append(cmd)
            cmd = f"[ -f '{decode_prefix}.tsv.gz' ] && touch '{decode_prefix}.done'"
            cmds.append(cmd)
            mm.add_target(f"{decode_prefix}.done", [modelf_done], cmds)
            sample_decode_done_files.append(f"{decode_prefix}.done")
        cmds = cmd_separator([], f"Finishing pixel-level decoding for all samples...")
        cmd = ""
        for sample_donef in sample_decode_done_files:
            cmd += f"[ -f '{sample_donef}' ] && "
        cmd += f"touch {decode_all_prefix}.done"
        cmds.append(cmd)
        mm.add_target(f"{decode_all_prefix}.done", sample_decode_done_files, cmds)

    ## step 4. write the output JSON file for each sample
    json_each_targets = []
    for i in range(len(in_samples)):
        sample = in_samples[i]
        #cmds = cmd_separator([], f"Writing output JSON file for sample {sample}...")
        sample_out_dir = os.path.join(args.out_dir, "samples", sample)
        sample_out_json = os.path.join(sample_out_dir, f"ficture.{args.out_prefix}.params.json")
        
        summary_aux_args = []
        prerequisities = [f"{args.out_dir}/multi.done"]

        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        sample_prefix = f"{sample_out_dir}/{sample}.{args.out_prefix}"

        out_manifolds = {
            "tsne": {
                "tsv": f"{sample_prefix}.tsne.leiden.tsv.gz",
                "png": f"{sample_prefix}.tsne.png"
            },
            "umap": {
                "tsv": f"{sample_prefix}.umap.leiden.tsv.gz",
                "png": f"{sample_prefix}.umap.png"
            }
        }
        out_cell_params = { 
            "model_type": "lda",
            "model_id": args.out_prefix,
            "cmap": f"{sample_prefix}.leiden.pseudobulk.cmap.tsv",
            "model_path": f"{lda_prefix}.model.tsv",
            "fit_path": f"{sample_prefix}.lda.results.tsv",
            "sptsv_prefix": f"{sample_prefix}.sptsv",
            "cell_xy_path": f"{sample_prefix}.leiden.xy.tsv.gz",
            "cluster_path": f"{sample_prefix}.leiden.tsv.gz",
            "cluster_pseudobulk": f"{sample_prefix}.leiden.pseudobulk.tsv",
            "cluster_de": f"{sample_prefix}.leiden.pseudobulk.de.tsv",
            "cluster_model_heatmap_pdf": f"{sample_prefix}.heatmap.pdf",
            "cluster_model_heatmap_tsv": f"{sample_prefix}.heatmap.counts.tsv",
            "pixel_png_path": f"{sample_prefix}.pixel.png",
            "pixel_tsv_path": f"{sample_prefix}.pixel.tsv.gz",
            "manifolds": out_manifolds
        }
        if sample in samp2boundaries:
            out_cell_params["cell_boundaries_path"] = samp2boundaries[sample]
        out_json = { "cell_params": out_cell_params }

        out_json_path = sample_out_json
        json.dump(out_json, flexopen(out_json_path, "wt"), indent=4)

    ## write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

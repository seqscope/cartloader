import sys, os, gzip, argparse, logging, shutil, subprocess, inspect, json
from venv import logger
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, flexopen, execute_makefile
from cartloader.utils.geometry_helper import iter_geojson_cell_centroids, CENTROID_SUPPORTED_FORMATS

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Run FICTURE2")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile, and print commands without executing them')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and start from the beginning')
    run_params.add_argument('--threads', type=int, default=8, help='Maximum number of threads per job (default: 8)')
    run_params.add_argument('--n-jobs', type=int, default=2, help='Number of parallel jobs to run (default: 2)')
    run_params.add_argument('--makefn', type=str, help='File name of Makefile to write (default: run_ficture2_multi.mk)')

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Enable all actions: --cells and --boundaries')
    cmd_params.add_argument('--sptsv', action='store_true', default=False, help='Create SPTSV files for LDA clustering')
    cmd_params.add_argument('--lda', action='store_true', default=False, help='Perform LDA factorization')
    cmd_params.add_argument('--leiden', action='store_true', default=False, help='Generate Leiden clusters based on LDA factorization')
    cmd_params.add_argument('--tsne', action='store_true', default=False, help='Generate TSNE manifolds based on LDA factorization')
    cmd_params.add_argument('--umap', action='store_true', default=False, help='Generate UMAP manifolds based on LDA factorization')
    cmd_params.add_argument('--pseudobulk', action='store_true', default=False, help='Generate pseudobulk files based on Leiden clusters')
    cmd_params.add_argument('--heatmap', action='store_true', default=False, help='Generate heamap between LDA factors and Leiden clusters')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='Perform pixel-level decoding based on cell clusters.')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output parameters for FICTURE")
    inout_params.add_argument('--out-dir', required=True, type=str, help='Output directory')
    inout_params.add_argument('--out-prefix', type=str, default="cells", help='Prefix for output files (default: cells)')
    inout_params.add_argument('--out-json', type=str, default=None, help="Path to output JSON file to store analysis parameters (default: <out-dir>/ficture.cells.params.json)")
    inout_params.add_argument('--in-dir', type=str, default=None, help='Path to input directory containing FICTURE output files. Same to out_dit if not specified')
    # inout_params.add_argument('--mex-dir', type=str, help='Directory containing MEX files')
    inout_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='Barcode files in MEX format')
    inout_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='Feature files in MEX format')
    inout_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='Matrix files in MEX format')
    inout_params.add_argument('--mex-list', type=str, help='TSV file containing sample IDs and paths to MEX files')
    inout_params.add_argument('--tsv-list', type=str, help='TSV file of [SAMPLE_ID] [PIXEL_TSV] naming an external headerless pixel TSV that carries a cell-id column. Use when the cell assignment lives in a separate file that cannot be mapped onto the tiled transcript (e.g. a Stereo-seq cell-bin GEM). The file is read with the same --colidx-* columns as the tiled transcript (X, Y, gene, count, cell_id).')
    inout_params.add_argument('--sptsv-prefix', type=str, help='Prefix for SPTSV files')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--n-factor', type=int, help='Number of factors for LDA training.')
    key_params.add_argument('--leiden-resolution', type=float, default=1.0, help='Resolution for Leiden clustering (default: 1.0)')
    key_params.add_argument('--anchor-resolution', type=int, default=6, help='Anchor resolution for decoding (default: 6)')
    key_params.add_argument('--cmap-file', type=str, default=os.path.join(repo_dir, "assets", "default_color_map.tsv"), help='Path to fixed color map TSV (default: <cartloader_dir>/assets/default_color_map.tsv)')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--colidx-feature',  type=int, default=3, help='Column index for feature in the --in-transcript (default: 3)')
    aux_params.add_argument('--colidx-count',  type=int, default=4, help='Column index for intensity in the --in-transcript (default: 4)')
    aux_params.add_argument('--colidx-cell-id', type=int, default=5, help='Column index for cell ID in the --in-transcript (default: 5)')
    aux_params.add_argument('--ignore-ids', type=str, default="UNASSIGNED,NA,0,-1", help='IDs to ignore in pixel file')
    # train
    aux_params.add_argument('--train-epoch', type=int, default=2, help='Training epoch for LDA model (default: 2)')
    #aux_params.add_argument('--skip-umap', action='store_true', default=False, help='Skip creating umap')
    aux_params.add_argument('--decode-scale', type=int, default=1, help='Decode scale (default: 1)')
    aux_params.add_argument('--seed', type=int, default=1, help='Random seed for random number generation (default: 1)')
    aux_params.add_argument('--single-molecule', action='store_true', default=False, help='Turn on single-molecule mode for pixel decode')
    aux_params.add_argument('--decode-pixel-res', type=float, default=0.5, help='Decode resolution (default: 0.5)')

    # others parameters shared across steps
    aux_params.add_argument('--min-feature-count', type=int, default=20, help='Minimum feature count for LDA factorization')
    aux_params.add_argument('--min-cell-count', type=int, default=50, help='Minimum cell count for LDA factorization (not in effect with --mex-dir option)')
    aux_params.add_argument('--de-min-ct-per-feature', type=int, default=20, help='Minimum count per feature for differential expression (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='P-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')
    aux_params.add_argument('--decode-fit-width', type=int, default=18, help='Fitting width (in microns) for decoding (default: 18)')
    # project from external model
    aux_params.add_argument('--pretrained-model', type=str, help='Path to a pre-trained model to use for projection. If provided, LDA training will be skipped, and the provided model will be used for projection.')
    aux_params.add_argument('--list-samples', type=str, help='Path to a TSV file containing sample IDs and paths to their transcript TSV files for multi-sample analysis. If provided, the samples listed in the file will be used for analysis.')
    aux_params.add_argument('--list-cluster', type=str, help='Path to a existing cluster files to create pseudobulk matrix in the format of [SAMPLE_ID] [CLUSTER_FILE]. If provided, Leiden clustering will be skipped, and the provided cluster files will be used for pseudobulk generation.')
    aux_params.add_argument('--list-xy', type=str, help='Path to a existing file containing X/Y locations of each cell. If provided, X/Y locations will be read from the provided file instead of computing from pixel file.')
    aux_params.add_argument('--list-boundaries', type=str, help='Path to an existing file containing cell boundaries, as [SAMPLE_ID] [BOUNDARIES_FILE] or [SAMPLE_ID] [BOUNDARIES_FILE] [SCALE_JSON]. The path is stored in the output JSON; when a sample has boundaries but no --list-xy entry and --boundaries-format supports it, per-cell centroids are derived from the polygons (rescaled into microns using microns_per_pixel from the optional third-column SCALE_JSON).')
    aux_params.add_argument('--xy-colname-cell-id', type=str, default="cell_id", help='Column name for cell IDs in the metadata file (default: cell_id)')
    aux_params.add_argument('--xy-colname-x', type=str, default="X", help='Column name for X coordinates in the metadata file (default: X)')
    aux_params.add_argument('--xy-colname-y', type=str, default="Y", help='Column name for Y coordinates in the metadata file (default: Y)')
    aux_params.add_argument('--boundaries-format', type=str, default=None, choices=[None, "geojson"], help='Format of the --list-boundaries files. Set (e.g. "geojson") to enable deriving per-cell centroids from the boundary polygons for samples that have boundaries but no --list-xy entry (default: None, i.e. boundaries are pass-through only).')
    aux_params.add_argument('--boundaries-cell-id-format', type=str, default="cellid_{:09d}-1", help='Python format string applied to each boundary feature cell id so the derived centroid id matches the cell barcode convention used for clustering (default: "cellid_{:09d}-1", the Visium HD segmented convention). Use "{}" to keep the raw id.')
    aux_params.add_argument('--boundaries-cell-id-prop', type=str, default="cell_id", help='GeoJSON feature property holding the cell id (default: cell_id)')
    aux_params.add_argument('--boundaries-units-key', type=str, default="microns_per_pixel", help='Key in the third-column SCALE_JSON giving microns per pixel; centroids are rescaled by this value into microns (default: microns_per_pixel)')
    aux_params.add_argument('--zero-based-clust-id', action='store_true', default=False, help='Whether the cluster IDs in the existing cluster files provided by --list-cluster are zero-based. By default, it is assumed that the cluster IDs are one-based and will be converted to zero-based by subtracting 1. If the cluster IDs are already zero-based, please turn on this option to avoid incorrect cluster ID conversion.')

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Customizing Auxiliary Parameters", "Customize features (typically genes) used by FICTURE without altering the original feature TSV") # This ensures the original feature TSV file is retained in the output JSON file for downstream processing 
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='Regex of feature names to include')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    #env_params.add_argument('--sort-mem', type=str, default="1G", help='Memory size for each process (default: 1G)')
    env_params.add_argument('--spatula', type=str, default=f"{repo_dir}/submodules/spatula/bin/spatula",  help='Path to spatula binary (default: "spatula" in the system PATH)') # default=f"{repo_dir}/submodules/spatula/bin/spatula",
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"), help='Path to punkst (ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--python', type=str, default="python3",  help='Python3 binary')
    env_params.add_argument('--R', type=str, default="Rscript", help='Path to R binary for UMAP generation (default: Rscript)')


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

    if args.makefn is None:
        args.makefn = f"run_ficture2_multi_cells.{args.out_prefix}.mk"

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
    if args.list_samples is not None:
        assert os.path.exists(args.list_samples), f"File not found: {args.list_samples} (--list-samples)"
        with flexopen(args.list_samples, "rt") as f:
            for line in f:
                toks = line.strip().split("\t")
                sample_id = toks[0]
                ## make sure that the sample directory exists
                sample_dir = os.path.join(args.in_dir, "samples", sample_id)
                if not os.path.exists(sample_dir):
                    raise FileNotFoundError(f"Sample directory not found: {sample_dir}. Please make sure that the sample directory exists in --in-samples")
                in_samples.append(sample_id)
        logger.info(f"Found {len(in_samples)} samples: {in_samples}") 

    else:     ## list directories in args.in_dir/samples/
        samples_dir = os.path.join(args.in_dir, "samples")
        for entry in os.listdir(samples_dir):
            entry_path = os.path.join(samples_dir, entry)
            if os.path.isdir(entry_path):
                in_samples.append(entry)
        logger.info(f"Found {len(in_samples)} samples in {samples_dir}: {in_samples}") 

    n_samples = len(in_samples)

    ## parse the MEX list (cell x gene matrices) into samp2mex: sample_id -> (bcd, ftr, mtx).
    ## A sample present here is clustered from its MEX counts (mex2sptsv); a sample absent
    ## here is clustered from the tiled transcript's cell_id column (pixel2sptsv). This is
    ## resolved per sample so a mixed run (some samples MEX, some transcript-based) works.
    samp2mex = {}
    if args.mex_list is not None:
        with flexopen(args.mex_list, 'rt') as rf:
            for line in rf:
                toks = line.strip().split("\t")
                if len(toks) == 0 or toks[0] == "":
                    continue
                sample_id = toks[0]
                if len(toks) == 2:
                    mex_dir = toks[1]
                    samp2mex[sample_id] = (os.path.join(mex_dir, args.mex_bcd),
                                           os.path.join(mex_dir, args.mex_ftr),
                                           os.path.join(mex_dir, args.mex_mtx))
                elif len(toks) == 4:
                    samp2mex[sample_id] = (toks[1], toks[2], toks[3])
                else:
                    raise ValueError(f"Each line in --mex-list must have 2 or 4 columns. Found {len(toks)} columns in line: {line}")

    ## parse the external pixel-TSV list into samp2tsv: sample_id -> pixel TSV path.
    ## Such a sample is clustered from that file's cell-id column instead of the tiled
    ## transcript's, for platforms whose cell assignment cannot be mapped back onto the
    ## pixel-level data (Stereo-seq cell bins). --mex-list wins if a sample is in both.
    samp2tsv = {}
    if args.tsv_list is not None:
        with flexopen(args.tsv_list, 'rt') as rf:
            for line in rf:
                toks = line.strip().split("\t")
                if len(toks) == 0 or toks[0] == "":
                    continue
                if len(toks) != 2:
                    raise ValueError(f"Each line in --tsv-list must have exactly 2 columns containing [SAMPLE_ID] [PIXEL_TSV]. Found {len(toks)} columns in line: {line}")
                if not os.path.exists(toks[1]):
                    raise FileNotFoundError(f"File not found: {toks[1]} (from --tsv-list)")
                samp2tsv[toks[0]] = toks[1]

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

    cmd_ftr_include_exclude = ""
    if args.include_feature_regex is not None:
        cmd_ftr_include_exclude += f" --include-feature-regex '{args.include_feature_regex}'"
    if args.exclude_feature_regex is not None:
        cmd_ftr_include_exclude += f" --exclude-feature-regex '{args.exclude_feature_regex}'"

    ## create cell-based SPTSV files
    if args.sptsv:
        if args.sptsv_prefix is not None:
            raise ValueError("When --sptsv is ON, --sptsv-prefix should not be provided.")
        sptsv_prefix = os.path.join(args.out_dir, args.out_prefix) + ".sptsv"
        cmds = cmd_separator([], f"Creating cell-based SPTSV files...")
        cmds.append(f"touch '{sptsv_prefix}.begin'")
        samp2sptsv = {} ## sample ID to SPTSV file mapping
        deps = []
        ## Resolve the cell-count source per sample: a sample listed in --mex-list is
        ## clustered from its MEX matrix (mex2sptsv); one listed in --tsv-list from that
        ## external pixel TSV's cell_id column; any other from the tiled transcript's
        ## cell_id column (both via pixel2sptsv). This mix lets a joint run combine
        ## MEX-based samples (e.g. MERSCOPE cell_by_gene without boundaries) with
        ## transcript/boundary-based samples in a single decode. Every branch writes the
        ## same per-sample sptsv prefix, so the steps below are source-agnostic.
        def cmd_pixel2sptsv(pixel_tsv, out_prefix):
            return (f"{args.spatula} pixel2sptsv --min-cell-count {args.min_cell_count} --pixel {pixel_tsv} "
                    f"--no-header --idx-col-x {args.colidx_x} --idx-col-y {args.colidx_y} "
                    f"--idx-col-ftr {args.colidx_feature} --idx-col-cnt {args.colidx_count} "
                    f"--idx-col-id {args.colidx_cell_id} --ignore-ids {args.ignore_ids} "
                    f"--out {out_prefix} --min-feature-count {args.min_feature_count} {cmd_ftr_include_exclude}")

        for sample_id in in_samples:
            sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
            if sample_id in samp2mex:
                mex_bcd, mex_ftr, mex_mtx = samp2mex[sample_id]
                cmd = f"{args.spatula} mex2sptsv --bcd {mex_bcd} --ftr {mex_ftr} --mtx {mex_mtx} --out {sample_sptsv_prefix} --min-feature-count {args.min_feature_count} {cmd_ftr_include_exclude}"
                deps.extend([mex_bcd, mex_ftr, mex_mtx])
            elif sample_id in samp2tsv:
                cmd = cmd_pixel2sptsv(samp2tsv[sample_id], sample_sptsv_prefix)
                deps.append(samp2tsv[sample_id])
            else:
                pixelf = f"{args.in_dir}/samples/{sample_id}/{sample_id}.tiled"
                cmd = cmd_pixel2sptsv(f"{pixelf}.tsv", sample_sptsv_prefix)
                deps.append(f"{pixelf}.tsv")
            cmds.append(cmd)
            samp2sptsv[sample_id] = sample_sptsv_prefix
        
        ## merge SPTSV files if needed
        if len(samp2sptsv) > 0:
            ## create a list file
            samp_listf = f"{sptsv_prefix}.list.tsv"
            with flexopen(samp_listf, "wt") as wf:
                for sample_id in samp2sptsv:
                    sample_sptsv_prefix = samp2sptsv[sample_id]
                    wf.write(f"{sample_id}\t{sample_sptsv_prefix}.feature.counts.tsv\t{sample_sptsv_prefix}.tsv\t{sample_sptsv_prefix}.json\n")
            cmd = f"{args.spatula} merge-sptsv --list {samp_listf} --out {sptsv_prefix}"
            cmds.append(cmd)
        ## randomize SPTSV file
        cmd = f"sort -k 1,1 {sptsv_prefix}.tsv > {sptsv_prefix}.randomized.tsv"
        cmds.append(cmd)
        cmds.append(f"[ -f {sptsv_prefix}.randomized.tsv ] && touch {sptsv_prefix}.done" )
        mm.add_target(f"{sptsv_prefix}.done", deps, cmds)
    elif args.sptsv_prefix is not None:
        sptsv_prefix = args.sptsv_prefix
    else:
        raise ValueError("When --sptsv is not provided, --sptsv-prefix must be provided.")

    ## perform LDA training
    if args.lda:
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        if args.pretrained_model is None:  ## run LDA to generate model
            cmds = cmd_separator([], f"Performing LDA training/projection...")
            cmds.append(f"touch {lda_prefix}.multi.begin")
            if args.n_factor is None:
                raise ValueError("--n-factor must be specified when --model is not specified with --lda ON.")
            cmd = f"{ficture2bin} lda4hex --in-data {sptsv_prefix}.randomized.tsv --in-meta {sptsv_prefix}.json --out-prefix {lda_prefix} --sort-topics --n-topics {args.n_factor} --transform --residuals--minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            cmds.append(cmd)
            cmds.append(f"[ -f {lda_prefix}.model.tsv ] && [ -f {lda_prefix}.results.tsv ] && touch {lda_prefix}.multi.done" )
            mm.add_target(f"{lda_prefix}.multi.done", [f"{sptsv_prefix}.done"], cmds)
        else:  ## use existing model
            cmds = cmd_separator([], f"Projecting existing LDA model...")
            cmds.append(f"touch {lda_prefix}.multi.begin")
            ## copy the pretrained model to lda_prefix
            if args.pretrained_model.endswith(".gz"):
                cmd = f"{args.gzip} -dc {args.pretrained_model} > {lda_prefix}.model.tsv"
            else:
                cmd = f"cp {args.pretrained_model} {lda_prefix}.model.tsv"
            cmds.append(cmd)
            cmd = f"{ficture2bin} lda4hex --model-prior {lda_prefix}.model.tsv --projection-only --in-data {sptsv_prefix}.randomized.tsv --in-meta {sptsv_prefix}.json --out-prefix {lda_prefix} --transform --residuals --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            cmds.append(cmd)
            cmds.append(f"[ -f '{lda_prefix}.results.tsv' ] && touch '{lda_prefix}.multi.done'" )
            mm.add_target(f"{lda_prefix}.multi.done", [f"{sptsv_prefix}.done"], cmds)
        ## project the LDA model to each sample separately
        deps = [f"{lda_prefix}.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Performing LDA projection for {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
            cmds.append(f"touch {sample_lda_prefix}.begin")
            cmd = f"{ficture2bin} lda4hex --model-prior {lda_prefix}.model.tsv --projection-only --in-data {sample_sptsv_prefix}.tsv --in-meta {sample_sptsv_prefix}.json --out-prefix {sample_lda_prefix} --transform --residuals --minibatch-size 500 --seed {args.seed} --n-epochs 2 --threads {args.threads}"
            cmds.append(cmd)
            cmds.append(f"[ -f {sample_lda_prefix}.results.tsv ] && touch {sample_lda_prefix}.done" )
            mm.add_target(f"{sample_lda_prefix}.done", [f"{lda_prefix}.multi.done"], cmds)
            deps.append(f"{sample_lda_prefix}.done")
        ## final target
        cmds = cmd_separator([], f"Finalizing LDA projection for all samples...")
        cmds.append(f"touch {lda_prefix}.done")
        mm.add_target(f"{lda_prefix}.done", deps, cmds);

    samp2boundaries = {}
    xy_samples = set()   # sample ids for which a per-cell scatter (cell.xy) was produced
    if args.leiden:
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        cmds = cmd_separator([], f"Generating Leiden clusters...")
        cmds.append(f"touch '{leiden_prefix}.begin'")
        if args.list_cluster is None:
            cmd = f"cartloader lda_leiden_cluster_fast --offset-data 4 --tsv '{lda_prefix}.results.tsv' --out '{leiden_prefix}.tsv.gz' --resolution {args.leiden_resolution} --colname-cluster topK --key-ids sample_id cell_id"
            cmds.append(cmd)
            for sample_id in in_samples:
                sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
                sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
                sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
                #cmd = f"cartloader lda_leiden_cluster --offset-data 3 --tsv '{sample_lda_prefix}.results.tsv' --out '{sample_leiden_prefix}.tsv.gz' --resolution {args.leiden_resolution} --colname-cluster topK --key-ids cell_id"
                cmd = f"({args.gzip} -cd '{leiden_prefix}.tsv.gz' | head -1 | cut -f 2-; {args.gzip} -cd '{leiden_prefix}.tsv.gz' | grep -w ^{sample_id} | cut -f 2- ;) | {args.gzip} -c > '{sample_leiden_prefix}.tsv.gz'"
                cmds.append(cmd)
        else:
            samp2clust = {}
            with flexopen(args.list_cluster, "rt") as rf:
                for line in rf:
                    toks = line.strip().split("\t")
                    if len(toks) != 2:
                        raise ValueError(f"Each line in --list-cluster must have exactly 2 columns containing [SAMPLE_ID] [CLUSTER_FILE] [METADTA_FILE]")
                    sample_id = toks[0]
                    cluster_file = toks[1]
                    if not os.path.exists(cluster_file):
                        raise FileNotFoundError(f"File not found: {cluster_file} (from --list-cluster)")
                    samp2clust[sample_id] = cluster_file
            with flexopen(f"{leiden_prefix}.tsv.gz", "wt") as wf:
                wf.write("sample_id\tcell_id\ttopK\n")
                for sample_id in in_samples: 
                    sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
                    if sample_id not in samp2clust:
                        logger.warning(f"Sample {sample_id} not found in --list-cluster. Writing non-informative file")
                        with flexopen(f"{sample_leiden_prefix}.tsv.gz", "wt") as wf_sample:
                            dummy_cell_id = f"{sample_id}_dummy_cell_id"
                            wf_sample.write("cell_id\ttopK\n")
                            wf_sample.write(f"{dummy_cell_id}\tNA\n")
                            wf.write(f"{sample_id}\t{dummy_cell_id}\tNA\n")
                    else:
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
                                    if args.zero_based_clust_id:
                                        int_cluster_id = int(cluster_id)
                                    else:
                                        int_cluster_id = int(cluster_id)-1 ## convert to 0-based
                                    if int_cluster_id < 0:
                                        raise ValueError(f"Cluster ID must be >= {0 if args.zero_based_clust_id else 1} in existing cluster file. Found {cluster_id} in line: {line}")
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
        # Boundary files (optional 3rd column = per-sample scale JSON for unit rescaling).
        # Parsed before the scatter loop so that centroids can be derived for samples that
        # have boundaries but no --list-xy entry (e.g. default Visium HD segmentation).
        samp2boundaries_scale = {}
        if args.list_boundaries is not None:
            with flexopen(args.list_boundaries, "rt") as rf:
                for line in rf:
                    toks = line.strip().split("\t")
                    if len(toks) not in (2, 3):
                        raise ValueError("Each line in --list-boundaries must have 2 or 3 columns: [SAMPLE_ID] [BOUNDARIES_FILE] [SCALE_JSON (optional)]")
                    if not os.path.exists(toks[1]):
                        raise FileNotFoundError(f"File not found: {toks[1]} (from --list-boundaries)")
                    samp2boundaries[toks[0]] = toks[1]
                    if len(toks) == 3 and toks[2]:
                        samp2boundaries_scale[toks[0]] = toks[2]

        def _boundary_units_per_um(sid):
            # units_per_um = coordinate units per micron = 1 / microns_per_pixel, so the
            # shared helper rescales polygon coordinates (pixels) into microns. No scale
            # JSON => assume coordinates are already in microns (units_per_um = 1).
            sj = samp2boundaries_scale.get(sid)
            if not sj:
                return 1.0
            with open(sj) as jf:
                mpp = json.load(jf).get(args.boundaries_units_key)
            if not mpp:
                raise ValueError(f"'{args.boundaries_units_key}' missing or zero in scale JSON {sj} (from --list-boundaries)")
            return 1.0 / float(mpp)

        merge_cmd = ""
        for sample_id in in_samples:
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_sptsv_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.sptsv"
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
            elif sample_id in samp2boundaries and args.boundaries_format in CENTROID_SUPPORTED_FORMATS:
                # No pre-computed cell XY, but boundary polygons are available (e.g. the
                # default Visium HD segmentation, which ships cell polygons but no
                # centroids): derive per-cell centroids from the polygons so the
                # leiden-cluster scatter and cell-point PMTiles can still be produced.
                metaf = f"{sample_sptsv_prefix}.cell.xy.tsv"
                upp = _boundary_units_per_um(sample_id)
                derived_ids = []
                with flexopen(metaf, "wt") as wf_sample:
                    wf_sample.write("cell_id\tX\tY\n")
                    for cid, cx, cy in iter_geojson_cell_centroids(
                            samp2boundaries[sample_id], upp,
                            args.boundaries_cell_id_format, args.boundaries_cell_id_prop):
                        wf_sample.write(f"{cid}\t{cx}\t{cy}\n")
                        derived_ids.append(cid)
                # Sanity check (printed so it is visible in the run log): how many derived
                # centroid ids match the cell barcodes that drive clustering? Leiden ids
                # are a subset of these, so a low/zero match ratio means the cell_id
                # conventions disagree and the cell-point layer will be (near) empty --
                # verify --boundaries-cell-id-format / --boundaries-cell-id-prop.
                msg = (f"[boundaries->centroids] sample {sample_id}: derived "
                       f"{len(derived_ids)} centroids (units_per_um={upp:g})")
                if sample_id in samp2mex:
                    with flexopen(samp2mex[sample_id][0], "rt") as bf:
                        barcodes = {ln.split("\t")[0].strip().strip('"') for ln in bf if ln.strip()}
                    matched = len(set(derived_ids) & barcodes)
                    pct = (100.0 * matched / len(barcodes)) if barcodes else 0.0
                    msg += f"; matched {matched}/{len(barcodes)} MEX barcodes ({pct:.1f}%)"
                    if barcodes and matched == 0:
                        msg += "  <-- WARNING: ZERO matches; cell-point layer will be empty"
                print(msg, file=sys.stderr, flush=True)
            elif sample_id in samp2mex:
                # MEX-based clustering carries no cell coordinates (mex2sptsv writes no
                # per-cell metadata), so a MEX sample without an xy file or boundaries has
                # nothing to place spatially; skip its per-cell leiden-cluster scatter
                # rather than failing on a missing metadata file. A transcript/boundary-
                # based sample (not in samp2mex) still has pixel2sptsv metadata, drawn below.
                continue
            xy_samples.add(sample_id)   # a per-cell scatter (cell.xy) is produced below
            draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
            cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{metaf}' --tsv-clust '{sample_leiden_prefix}.tsv.gz' --tsv-colname-x X --tsv-colname-y Y --out '{sample_leiden_prefix}.xy.png' --out-tsv '{sample_leiden_prefix}.xy.tsv.gz' --tsv-colname-clust topK"
            cmds.append(cmd)
            merge_cmd += f"[ -f '{sample_leiden_prefix}.xy.done' ] && "
            cmds.append(f"[ -f '{sample_leiden_prefix}.xy.tsv.gz' ] && [ -f '{sample_leiden_prefix}.xy.png' ] && touch '{sample_leiden_prefix}.xy.done'" )
        merge_cmd += f"[ -f '{leiden_prefix}.tsv.gz' ] && touch '{leiden_prefix}.done'"
        cmds.append(merge_cmd)
        mm.add_target(f"{leiden_prefix}.done", [f"{lda_prefix}.done"], cmds)
    if args.tsne:
        ## generate TSNE manifolds
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        leiden_prefix = os.path.join(args.out_dir, args.out_prefix) + ".leiden"
        tsne_prefix = os.path.join(args.out_dir, args.out_prefix) 
        cmds = cmd_separator([], f"Generating TSNE manifolds...")
        cmds.append(f"touch {tsne_prefix}.tsne.begin")
        cmd = f"cartloader lda_tsne --offset-data 4 --tsv '{lda_prefix}.results.tsv' --out '{tsne_prefix}.tsne.tsv.gz' --key-ids sample_id cell_id"
        cmds.append(cmd)

        draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
        cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{tsne_prefix}.tsne.tsv.gz' --tsv-clust '{leiden_prefix}.tsv.gz' --tsv-colname-x TSNE1 --tsv-colname-y TSNE2 --out '{tsne_prefix}.tsne.png' --tsv-colname-ids sample_id cell_id --out-tsv '{tsne_prefix}.tsne.leiden.tsv.gz' --tsv-colname-clust topK"
        cmds.append(cmd)
        cmds.append(f"[ -f '{tsne_prefix}.tsne.leiden.tsv.gz' ] && [ -f '{tsne_prefix}.tsne.png' ] && touch '{tsne_prefix}.tsne.multi.done'" )
        mm.add_target(f"{tsne_prefix}.tsne.multi.done", [f"{lda_prefix}.done", f"{leiden_prefix}.done"], cmds)

        deps = [f"{tsne_prefix}.tsne.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Generating TSNE manifolds for sample {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_tsne_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}"
            if n_samples == 1: ## if sample size is 1, do not rerun TSNE
                cmd = f"{args.gzip} -cd {tsne_prefix}.tsne.tsv.gz | cut -f 2- | {args.gzip} -c > '{sample_tsne_prefix}.tsne.tsv.gz'"
            else:
                cmd = f"cartloader lda_tsne --offset-data 3 --tsv '{sample_lda_prefix}.results.tsv' --out '{sample_tsne_prefix}.tsne.tsv.gz'"
            cmds.append(cmd)

            draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
            cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{sample_tsne_prefix}.tsne.tsv.gz' --tsv-clust '{sample_leiden_prefix}.tsv.gz' --tsv-colname-x TSNE1 --tsv-colname-y TSNE2 --out '{sample_tsne_prefix}.tsne.png' --out-tsv '{sample_tsne_prefix}.tsne.leiden.tsv.gz' --tsv-colname-clust topK"
            cmds.append(cmd)
            cmds.append(f"[ -f '{sample_tsne_prefix}.tsne.leiden.tsv.gz' ] && [ -f '{sample_tsne_prefix}.tsne.png' ] && touch '{sample_tsne_prefix}.tsne.done'" )
            mm.add_target(f"{sample_tsne_prefix}.tsne.done", [f"{sample_lda_prefix}.done", f"{leiden_prefix}.done", f"{tsne_prefix}.tsne.multi.done"], cmds)
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
        cmds.append(f"touch {umap_prefix}.umap.begin")
        create_umap_rscript=f"{repo_dir}/cartloader/r/create_umap.r"
        cmd = f"{args.R} '{create_umap_rscript}' --input '{lda_prefix}.results.tsv' --out-prefix '{umap_prefix}' --tsv-colname-meta random_key sample_id cell_id"
        cmds.append(cmd)

        ## draw UMAP manifolds
        draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
        cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{umap_prefix}.umap.tsv.gz' --tsv-clust '{leiden_prefix}.tsv.gz' --tsv-colname-x UMAP1 --tsv-colname-y UMAP2 --out '{umap_prefix}.umap.png' --out-tsv '{umap_prefix}.umap.leiden.tsv.gz' --tsv-colname-ids sample_id cell_id --tsv-colname-clust topK"
        cmds.append(cmd)
        cmds.append(f"[ -f '{umap_prefix}.umap.leiden.tsv.gz' ] && [ -f '{umap_prefix}.umap.png' ] && touch '{umap_prefix}.umap.multi.done'" )
        mm.add_target(f"{umap_prefix}.umap.multi.done", [f"{lda_prefix}.done", f"{leiden_prefix}.done"], cmds)

        deps = [f"{umap_prefix}.umap.multi.done"]
        for sample_id in in_samples:
            cmds = cmd_separator([], f"Generating UMAP manifolds for sample {sample_id}...")
            sample_lda_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.lda"
            sample_leiden_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}.leiden"
            sample_umap_prefix = f"{args.out_dir}/samples/{sample_id}/{sample_id}.{args.out_prefix}"
            if n_samples == 1: ## if sample size is 1, do not rerun TSNE
                cmd = f"{args.gzip} -cd {umap_prefix}.umap.tsv.gz | cut -f 1,3- | {args.gzip} -c > '{sample_umap_prefix}.umap.tsv.gz'"
            else:
                cmd = f"{args.R} '{create_umap_rscript}' --input '{sample_lda_prefix}.results.tsv' --out-prefix '{sample_umap_prefix}' --tsv-colname-meta random_key cell_id"
            cmds.append(cmd)

            draw_manifold_rscript=f"{repo_dir}/cartloader/r/draw_manifold_clust.r"
            cmd = f"{args.R} '{draw_manifold_rscript}' --tsv-manifold '{sample_umap_prefix}.umap.tsv.gz' --tsv-clust '{sample_leiden_prefix}.tsv.gz' --tsv-colname-x UMAP1 --tsv-colname-y UMAP2 --out '{sample_umap_prefix}.umap.png' --out-tsv '{sample_umap_prefix}.umap.leiden.tsv.gz' --tsv-colname-clust topK"
            cmds.append(cmd)
            cmds.append(f"[ -f '{sample_umap_prefix}.umap.leiden.tsv.gz' ] && [ -f '{sample_umap_prefix}.umap.png' ] && touch '{sample_umap_prefix}.umap.done'" )
            mm.add_target(f"{sample_umap_prefix}.umap.done", [f"{sample_lda_prefix}.done", f"{leiden_prefix}.done", f"{umap_prefix}.umap.multi.done"], cmds)
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
        cmds.append(f"touch {pseudobulk_prefix}.begin")
        cmd = f"{args.spatula} sptsv2model --min-count {args.min_feature_count} --tsv '{sptsv_prefix}.randomized.tsv' --clust '{leiden_prefix}.tsv.gz' --features '{sptsv_prefix}.feature.counts.tsv' --json '{sptsv_prefix}.json' --out '{pseudobulk_prefix}.tsv'"
        cmds.append(cmd)

        cmd = f"head -n $(head -1 '{pseudobulk_prefix}.tsv' | wc -w) '{args.cmap_file}' > '{pseudobulk_prefix}.cmap.tsv'"
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
        cmds.append(f"touch {heatmap_prefix}.begin")
        #model_tsv = args.pretrained_model if args.pretrained_model is not None else f"{lda_prefix}.model.tsv"
        model_tsv = f"{lda_prefix}.model.tsv"
        cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{model_tsv}' --out '{lda_prefix}.model' --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}"
        cmds.append(cmd)
        cmds.append(f"({args.gzip} -cd {lda_prefix}.model.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {lda_prefix}.model.de.marginal.tsv.gz | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > {lda_prefix}.model.de.tsv")
        cmds.append(f"rm -f '{lda_prefix}.model.de.marginal.tsv.gz'")

        ## perform DE test on the cell pseudobulk matrix
        cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{pseudobulk_prefix}.tsv' --out '{pseudobulk_prefix}' --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}"
        cmds.append(cmd)
        cmds.append(f"({args.gzip} -cd {pseudobulk_prefix}.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {pseudobulk_prefix}.de.marginal.tsv.gz | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > {pseudobulk_prefix}.de.tsv")
        cmds.append(f"rm -f '{pseudobulk_prefix}.de.marginal.tsv.gz'")

        ## create factor report for cluster annotation
        cmd = " ".join([
            ficture2report,
            f"--factor_label factor",
            f"--de '{pseudobulk_prefix}.de.tsv'",
            f"--pseudobulk '{pseudobulk_prefix}.tsv'",
            f"--feature_label Feature",
            f"--color_table '{pseudobulk_prefix}.cmap.tsv'",
            f"--output_pref '{pseudobulk_prefix}.factor'",
            ])
        cmds.append(cmd)

        ## create heatmap
        heatmap_rscript=f"{repo_dir}/cartloader/r/create_heatmap.r"
        cmd = f"{args.R} '{heatmap_rscript}' --results '{lda_prefix}.results.tsv' --clust '{leiden_prefix}.tsv.gz' --de-results '{lda_prefix}.model.de.tsv' --de-clust '{pseudobulk_prefix}.de.tsv' --offset-data 3 --out '{heatmap_prefix}' --colname-clust topK --cell-clust sample_id cell_id --draw"
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

            ## perform DE test on the cell pseudobulk matrix
            cmd = f"{args.spatula} diffexp-model-matrix --tsv1 '{sample_pseudobulk_prefix}.tsv' --out '{sample_pseudobulk_prefix}' --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}"
            cmds.append(cmd)

            cmds.append(f"({args.gzip} -cd {sample_pseudobulk_prefix}.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {sample_pseudobulk_prefix}.de.marginal.tsv.gz | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > {sample_pseudobulk_prefix}.de.tsv")
            cmds.append(f"rm -f '{sample_pseudobulk_prefix}.de.marginal.tsv.gz'")

            cmd = " ".join([
                ficture2report,
                f"--factor_label factor",
                f"--de '{sample_pseudobulk_prefix}.de.tsv'",
                f"--pseudobulk '{sample_pseudobulk_prefix}.tsv'",
                f"--feature_label Feature",
                f"--color_table '{sample_pseudobulk_prefix}.cmap.tsv'",
                f"--output_pref '{sample_pseudobulk_prefix}.factor'",
                ])
            cmds.append(cmd)

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
            cmds.append(f"touch '{decode_prefix}.begin'")
            fit_width = args.decode_fit_width  ## e.g., 18um
            fit_n_move = fit_width // args.anchor_resolution + 1
            decode_id = f"p{fit_width}_a{args.anchor_resolution}"
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
                f"--single-molecule" if args.single_molecule else f"--pixel-res {args.decode_pixel_res}",
                f"--output-binary",
                f"--threads {args.threads}",
                f"--seed {args.seed}"
                #f"--output-original"
            ])
            cmds.append(cmd)
            cmd = f"'{ficture2bin}' draw-pixel-factors --in '{decode_prefix}' --binary --in-color '{args.cmap_file}' --out '{decode_prefix}.png' --scale {args.decode_scale} --range '{sample_prefix}.coord_range.tsv'"
            # cmd = " ".join([
            #     ficture2bin, "draw-pixel-factors",
            #     f"--in-tsv '{decode_prefix}.tsv'",
            #     f"--header-json '{decode_prefix}.json'",
            #     f"--in-color '{args.cmap_file}'",
            #     f"--out '{decode_prefix}.png'",
            #     f"--scale {args.decode_scale}",
            #     f"--range '{sample_prefix}.coord_range.tsv'"
            #     ])
            cmds.append(cmd)
            # cmd = f"{args.gzip} -f '{decode_prefix}.tsv'"
            # cmds.append(cmd)
            #cmd = f"[ -f '{decode_prefix}.tsv.gz' ] && touch '{decode_prefix}.done'"
            cmd = f"[ -f '{decode_prefix}.bin' ] && touch '{decode_prefix}.done'"
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

        shared_prefix = os.path.join(args.out_dir, args.out_prefix)
        lda_prefix = os.path.join(args.out_dir, args.out_prefix) + ".lda"
        sample_prefix = f"{sample_out_dir}/{sample}.{args.out_prefix}"

        out_manifolds = {
            "shared": {},
            "sample": {}
        }
        if args.tsne:
            out_manifolds["shared"]["tsne"] = {
                "tsv": f"{shared_prefix}.tsne.leiden.tsv.gz",
                "png": f"{shared_prefix}.tsne.png"
            }
            out_manifolds["sample"]["tsne"] = {
                "tsv": f"{sample_prefix}.tsne.leiden.tsv.gz",
                "png": f"{sample_prefix}.tsne.png"
            }
        
        if args.umap:
            out_manifolds["shared"]["umap"] = {
                "tsv": f"{shared_prefix}.umap.leiden.tsv.gz",
                "png": f"{shared_prefix}.umap.png"
            }
            out_manifolds["sample"]["umap"] = {
                "tsv": f"{sample_prefix}.umap.leiden.tsv.gz",
                "png": f"{sample_prefix}.umap.png"
            }

        out_cell_params = { 
            "model_type": "lda",
            "model_id": args.out_prefix,
        }
        
        if args.pseudobulk:
            out_cell_params["cmap"] = f"{sample_prefix}.leiden.pseudobulk.cmap.tsv"
            out_cell_params["cluster_pseudobulk"] = f"{sample_prefix}.leiden.pseudobulk.tsv"
            out_cell_params["cluster_de"] = f"{sample_prefix}.leiden.pseudobulk.de.tsv"
            out_cell_params["cluster_info"] = f"{sample_prefix}.leiden.pseudobulk.factor.info.tsv"
            out_cell_params["shared_cluster_pseudobulk"] = f"{shared_prefix}.leiden.pseudobulk.tsv"
            out_cell_params["shared_cluster_de"] = f"{shared_prefix}.leiden.pseudobulk.de.tsv"
            out_cell_params["shared_cluster_info"] = f"{shared_prefix}.leiden.pseudobulk.factor.info.tsv"
            
        if args.lda:
            out_cell_params["model_path"] = f"{lda_prefix}.model.tsv"
            out_cell_params["fit_path"] = f"{sample_prefix}.lda.results.tsv"
            
        if args.sptsv:
            out_cell_params["sptsv_prefix"] = f"{sample_prefix}.sptsv"
            
        if args.leiden:
            # cell_xy_path is only produced when the per-cell scatter actually ran — i.e.
            # the sample had cell XY (--list-xy) or centroids derived from its boundary
            # polygons. Coordinate-less MEX clustering produces none, so run_cartload2
            # skips cell-point PMTiles instead of failing on a missing file.
            if sample in xy_samples:
                out_cell_params["cell_xy_path"] = f"{sample_prefix}.leiden.xy.tsv.gz"
            out_cell_params["cluster_path"] = f"{sample_prefix}.leiden.tsv.gz"
            
        if args.heatmap:
            out_cell_params["cluster_model_heatmap_pdf"] = f"{sample_prefix}.heatmap.pdf"
            out_cell_params["cluster_model_heatmap_tsv"] = f"{sample_prefix}.heatmap.normfrac.tsv"
            out_cell_params["shared_cluster_model_heatmap_pdf"] = f"{shared_prefix}.heatmap.pdf"
            out_cell_params["shared_cluster_model_heatmap_tsv"] = f"{shared_prefix}.heatmap.normfrac.tsv"
            
        if args.decode:
            out_cell_params["pixel_png_path"] = f"{sample_prefix}.pixel.png"
            out_cell_params["pixel_bin_prefix"] = f"{sample_prefix}.pixel"
            out_cell_params["pixel_res"] = "0" if args.single_molecule else str(args.decode_pixel_res) ## use 0 to indicate single-molecule decoding
            out_cell_params["decode_scale"] = str(args.decode_scale)
            #out_cell_params["pixel_tsv_path"] = f"{sample_prefix}.pixel.tsv.gz"

        if args.tsne or args.umap:
            out_cell_params["manifolds"] = out_manifolds
        if n_samples > 1:
            out_cell_params["analysis_type"] = "multi-sample"
        # Advertise the boundary path for the cell-boundaries PMTiles layer only for
        # formats run_cartload2 can tile (vertex CSV). A geojson consumed here is used
        # solely to derive centroids (its boundaries layer is produced separately, e.g.
        # by import_visiumhd_cell), so it must not be passed on as cell_boundaries_path.
        if sample in samp2boundaries and args.boundaries_format not in CENTROID_SUPPORTED_FORMATS:
            out_cell_params["cell_boundaries_path"] = samp2boundaries[sample]
        out_json = { "cell_params": out_cell_params }

        out_json_path = sample_out_json
        json.dump(out_json, flexopen(out_json_path, "wt"), indent=4)

    ## write the shared multi-sample cell manifest (points to per-sample JSONs + shared cell components)
    sp = args.out_prefix   # paths relative to --out-dir
    shared = {"model_type": "lda", "model_id": args.out_prefix}
    if args.lda:
        shared["shared_model_path"] = f"{sp}.lda.model.tsv"
    if args.pseudobulk:
        shared["shared_cmap"] = f"{sp}.leiden.pseudobulk.cmap.tsv"
        shared["shared_cluster_pseudobulk"] = f"{sp}.leiden.pseudobulk.tsv"
        shared["shared_cluster_de"] = f"{sp}.leiden.pseudobulk.de.tsv"
        shared["shared_cluster_info"] = f"{sp}.leiden.pseudobulk.factor.info.tsv"
    if args.heatmap:
        shared["shared_cluster_model_heatmap_pdf"] = f"{sp}.heatmap.pdf"
        shared["shared_cluster_model_heatmap_tsv"] = f"{sp}.heatmap.normfrac.tsv"
    shared_manifolds = {}
    if args.tsne:
        shared_manifolds["tsne"] = {"tsv": f"{sp}.tsne.leiden.tsv.gz", "png": f"{sp}.tsne.png"}
    if args.umap:
        shared_manifolds["umap"] = {"tsv": f"{sp}.umap.leiden.tsv.gz", "png": f"{sp}.umap.png"}
    if shared_manifolds:
        shared["manifolds"] = shared_manifolds

    multi_manifest = {
        "analysis_type": "multi-sample",
        "n_samples": n_samples,
        "out_prefix": args.out_prefix,
        "samples": {s: os.path.join("samples", s, f"ficture.{args.out_prefix}.params.json") for s in in_samples},
        "shared": shared,
    }
    multi_json_path = os.path.join(args.out_dir, f"ficture.multi.{args.out_prefix}.params.json")
    json.dump(multi_manifest, flexopen(multi_json_path, "wt"), indent=4)

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

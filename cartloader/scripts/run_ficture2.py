import sys, os, gzip, argparse, logging, shutil, subprocess, inspect
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, flexopen, write_dict_to_file, load_file_to_dict, execute_makefile

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Run FICTURE2 tiling, segmentation, LDA training, decoding;  writes a Makefile and executes steps."
    )
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    run_params = parser.add_argument_group("Run Options", "Run options.")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile and print commands without executing them')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and start from the beginning')
    run_params.add_argument('--makefn', type=str, default="run_ficture2.mk", help='File name of Makefile to write (default: run_ficture2.mk)')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (default: 4)')

    cmd_params = parser.add_argument_group("Commands", "FICTURE steps to run")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run all main functions, including tile, segment, init-lda, decode') #, and summary')
    cmd_params.add_argument('--tile', action='store_true', default=False, help='(Main function) Perform tiling step')
    cmd_params.add_argument('--segment', action='store_true', default=False, help='(Main function) Hexagon segmentation into FICTURE-compatible format')
    cmd_params.add_argument('--init-lda', action='store_true', default=False, help='(Main function) Train LDA model(s)')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='(Main function) Pixel-level decoding')
    #cmd_params.add_argument('--summary', action='store_true', default=False, help='(Main function) Write JSON summarizing FICTURE parameters and outputs')
    cmd_params.add_argument('--segment-10x', action='store_true', default=False, help='Hexagon segmentation into 10x MEX format')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output paths")
    inout_params.add_argument('--out-dir', required=True, type=str, help='Output directory')
    inout_params.add_argument('--out-json', type=str, default=None, help='Path to output JSON manifest summarizing FICTURE parameters (default: <out-dir>/ficture.params.json)')
    inout_params.add_argument('--in-json', type=str, default=None, help='Path to input manifest JSON with --in-transcript/--in-minmax/--in-feature; if set, omit --in-transcript/--in-minmax/--in-feature (typical: sge_assets.json from cartloader sge_convert)')
    inout_params.add_argument('--in-transcript', type=str, default=None, help='Path to input unsorted transcript-indexed SGE TSV (default: <out-dir>/transcripts.unsorted.tsv.gz)')
    inout_params.add_argument('--in-minmax', type=str, default=None, help='Path to input coordinate min/max TSV (defaults to --out-dir/transcripts.tiled.coord_range.tsv after tiling)')
    inout_params.add_argument('--in-feature', type=str, default=None,  help='Path to input per-feature UMI count TSV')
    inout_params.add_argument('--in-feature-ficture', type=str, default=None, help='(Optional) Path to custom feature list for FICTURE; restricts analysis to listed features')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires users' attention")
    key_params.add_argument('--width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (µm) for LDA training')
    key_params.add_argument('--n-factor', type=str, default=None, help='Comma-separated factor counts for LDA training')
    key_params.add_argument('--anchor-res', type=int, default=6, help='Anchor resolution for decoding (default: 6)')
    # key_params.add_argument('--radius-buffer', type=int, default=1, help='Buffer to radius(=anchor_res + radius_buffer) for pixel-level decoding (default: 1)')
    key_params.add_argument('--cmap-file', type=str, default=os.path.join(repo_dir, "assets", "fixed_color_map_256.tsv"), help='Path to fixed color map TSV (default: <cartloader_dir>/assets/fixed_color_map_256.tsv)')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Additional parameters (defaults recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='1-based column index of X in --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='1-based column index of Y in --in-transcript (default: 2)')
    aux_params.add_argument('--colname-count', type=str, default="count", help='Column name for count in inputs (default: count)')
    aux_params.add_argument('--colname-feature', type=str, default="gene", help='Column name for feature in inputs (default: gene)')
    aux_params.add_argument('--tile-size', type=int, default=500, help='Tile size for tiling (default: 500)')
    aux_params.add_argument('--tile-buffer', type=int, default=1000, help='Tile buffer for tiling (default: 1000)')
    aux_params.add_argument('--seed', type=int, default=1, help='Random seed for random number generation (default: 1)')
    # segmentation - ficture
    aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=50, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format (default: 50)')
    # minibatch
    aux_params.add_argument('--minibatch-size', type=int, default=500, help='Batch size used in minibatch processing (default: 500)')
    # train
    aux_params.add_argument('--min-ct-per-unit-train', type=int, default=50, help='Minimum count for training (default: 50)')
    aux_params.add_argument('--train-epoch', type=int, default=2, help='Training epoch for LDA model (default: 2)')
    # fit
    aux_params.add_argument('--fit-width',  type=str, default=None, help='Hexagon flat-to-flat width (µm) during model fitting (defaults to --width)')
    # aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=20, help='Minimum count per hexagon unit during model fitting (default: 20)')
    # aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot (default: 1)')   # in Scopeflow, this is set to 2
    # 10x segment:
    key_params.add_argument('--segment-width-10x', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (µm) in 10x format (required if --segment-10x)')
    # decode
    aux_params.add_argument('--decode-scale', type=int, default=1, help='Scale factor from input coordinates to output image pixels (default: 1)')
    # others parameters shared across steps
    aux_params.add_argument('--min-ct-per-feature', type=int, default=20, help='Minimum count per feature during LDA training, transform and decoding (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='P-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')

    # AUX feacture-filtering params
    aux_ftrfilter_params = parser.add_argument_group("Feature Customizing Auxiliary Parameters", "Customize features (typically genes) used by FICTURE without altering the original feature TSV")
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='Regex of feature names to include')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude')
    
    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Paths to external tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary (default: gzip). For speed, consider "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary (default: sort). For faster processing, you may include flags like "sort -T /tmp --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default="1G", help='Sort memory limit per process (default: 1G)')
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary (default: spatula)')
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"),  help='Path to punkst (ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--python', type=str, default="python3",  help='Path to Python (default: python3)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args(_args)

def define_lda_runs(args):
    assert args.init_lda, "--init-lda must be ON when running define_lda_runs()"
    assert args.width is not None, "When --init-lda is ON, provide at least one train width for LDA training using --width"

    if args.n_factor is None:
        print("Warning: --n-factor is not provided for LDA. Using default values: 12,24")
        args.n_factor = "12,24"

    train_widths = [int(x) for x in args.width.split(",")] if args.width else [] #and not (args.use_external_model and args.external_model_type == "custom") else []
    n_factors = [int(x) for x in args.n_factor.split(",")] if args.n_factor else []  #and not (args.use_external_model and args.external_model_type == "custom") else []

    train_params= [
        {
         "model_type": "lda",
         "train_width": train_width,
         "n_factor": n_factor,
         "model_id":f"t{train_width}_f{n_factor}",
        }
        for train_width in train_widths
        for n_factor in n_factors
    ]
    return train_params

def define_training_runs(args):
    train_params = define_lda_runs(args)
    return train_params

def define_decode_runs(args):
    decode_runs = []
    train_params = define_training_runs(args)
    for train_param in train_params:
        model_type = train_param["model_type"]
        train_width = train_param["train_width"]
        n_factor = train_param["n_factor"]
        model_id = train_param["model_id"]
        model_prefix = os.path.join(args.out_dir, model_id)
        model_path= f"{model_prefix}.model.tsv"
        fit_widths = [train_width] if args.fit_width is None else [int(x) for x in args.fit_width.split(",")]
        for fit_width in fit_widths:
            decode_id = f"{model_id}_p{fit_width}_a{args.anchor_res}"
            cmap_path = f"{model_prefix}.cmap.tsv" #define_cmap(args, model_id)
            decode_runs.append({
                "model_type": model_type,
                "model_id": model_id,
                "model_path": model_path,
                "decode_id": decode_id,
                "n_factor": n_factor,
                "fit_width": fit_width,
                "cmap_path": cmap_path,
                "prerequisite_path": f"{model_prefix}.done"
            })
    return decode_runs

def run_ficture2(_args):
    """Run all functions in FICTURE by using GNU Makefile
    This function is meant to be used in a local environment that has sufficient resources to run all functions in FICTURE at once.
    This function performs the following tasks:
    (1) Take the input parameters relevant to the FICTURE runs
    (2) Identify the sequence of commands to run FICTURE
    (3) Create a GNU makefile to run the commands in parallel
    (4) Run the GNU makefile
    """
    # args
    args=parse_arguments(_args)

    if args.main:
        args.tile = True
        args.segment = True
        args.init_lda = True
        args.decode = True
        # args.summary = True

    #assert not (args.init_lda and args.init_ext), "Cannot choose both --init-lda and --init-ext"

    # input/output/other files
    # out dirs
    os.makedirs(args.out_dir, exist_ok=True)

    # out files
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.params.json")

    # in files
    ## read in_json if provided 
    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"File not found: {args.in_json} (--in-json)"

        print(f"Reading inputs from manifest: {args.in_json}")
        sge_data=load_file_to_dict(args.in_json)
        args.in_transcript = sge_data.get("transcript", None)
        args.in_minmax = sge_data.get("minmax", None)
        args.in_feature = sge_data.get("feature", None)
        print(f"Resolved input: --in-transcript {args.in_transcript}")
        print(f"Resolved input: --in-feature {args.in_feature}")
        print(f"Resolved input: --in-minmax {args.in_minmax}")
    
    if args.in_transcript is None:
        args.in_transcript = os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")
    assert os.path.exists(args.in_transcript), f"File not found: {args.in_transcript} (--in-transcript)"

    if args.in_minmax is not None:
        assert os.path.exists(args.in_minmax), f"File not found: {args.in_minmax} (--in-minmax)"
    
    if args.in_feature is not None:
        assert os.path.exists(args.in_feature), f"File not found: {args.in_feature} (--in-feature)"

    if args.in_feature_ficture is not None: 
        assert os.path.exists(args.in_feature_ficture), f"File not found: {args.in_feature_ficture} (--in-feature-ficture)"
    
    assert os.path.exists(args.cmap_file), f"File not found: {args.cmap_file} (--cmap-file)"

    # FICTURE installations
    ficture2bin = os.path.join(args.ficture2, "bin/punkst")
    assert os.path.exists(ficture2bin), f"File not found: {ficture2bin}. FICTURE2 Directory should include bin/punkst (--ficture2)"
    
    ficture2de = args.python + " " + os.path.join(args.ficture2, "ext/py/de_bulk.py")
    ficture2report = args.python + " " + os.path.join(args.ficture2, "ext/py/factor_report.py")

    # start mm
    mm = minimake()

    # 1. tiling :
    tile_flag = f"{args.out_dir}/transcripts.tiled.done"
    if args.tile:
        scheck_app(args.gzip)
        cmds = cmd_separator([], f"Creating tiled tsv from {os.path.basename(args.in_transcript)}...")
        cmd = " ".join([
            ficture2bin, "pts2tiles",
            f"--in-tsv {args.in_transcript}",
            f"--out-prefix {args.out_dir}/transcripts.tiled",
            f"--icol-x {args.colidx_x-1}",
            f"--icol-y {args.colidx_y-1}",
            f"--icol-feature 2",
            f"--icol-int 3",
            f"--skip 1",
            f"--temp-dir {args.out_dir}/tmp",
            f"--tile-size {args.tile_size}",
            f"--tile-buffer {args.tile_buffer}",
            f"--threads {args.threads}"
        ])
        cmds.append(cmd)

        cmds.append(f"[ -f {args.out_dir}/transcripts.tiled.tsv ] && [ -f {args.out_dir}/transcripts.tiled.index ] && touch {args.out_dir}/transcripts.tiled.done" )
        mm.add_target(f"{args.out_dir}/transcripts.tiled.done", [args.in_transcript], cmds)

    # - assign minmax
    # Note: for minmax/ftr, defaults to the minmax/ftr from the tiled step. Use the corresponding prerequisites to check if it is feasible.
    if args.in_minmax is None:
        args.in_minmax = os.path.join(args.out_dir, "transcripts.tiled.coord_range.tsv")
        if args.tile:
            minmax_prereq = f"{args.out_dir}/transcripts.tiled.done"
        else:
            minmax_prereq = args.in_minmax
    else:
        minmax_prereq = args.in_minmax
    
    # - Prepare feature files
    if args.in_feature_ficture is not None:
        in_feature_ficture = args.in_feature_ficture
    else:
        in_feature_ficture = args.in_feature
    
    if in_feature_ficture is not None:
        ## if a specific feature file is provided, use it as "selected features" for FICTURE analysis            
        assert os.path.exists(in_feature_ficture), f"File not found: {in_feature_ficture}" + (f"(--in-feature-ficture)" if args.in_feature_ficture is not None else "--in-feature")
        feature_nohdr = f"{args.out_dir}/transcripts.tiled.selected_features.tsv"       
        feature_plain = f"{args.out_dir}/transcripts.tiled.selected_features.hdr.tsv"

        cmds = cmd_separator([], f"Preparing a feature file for FICTURE from a specified feature file: {in_feature_ficture}")

        # create a feature without header 
        if in_feature_ficture.lower().endswith(".gz"):
            cmds.append(f"""gzip -dc {in_feature_ficture} | tail -n +2 | awk -F'\\t' '{{print $1"\\t"$(NF);}}' > {feature_nohdr}""")
        else:
            cmds.append(f"""tail -n +2 {in_feature_ficture} | awk -F'\\t' '{{print $1"\\t"$(NF);}}' > {feature_nohdr}""")
        
        # plain file 
        cmds.append(f'echo -e "{args.colname_feature}\\t{args.colname_count}" > {feature_plain}')
        cmds.append(f"cat {feature_nohdr} >> {feature_plain}")
        
        mm.add_target(feature_plain, [in_feature_ficture], cmds)
    else:
        ## if feature file is not provided, create transcripts.tiled.features.tsv from pts2tiles
        feature_nohdr = f"{args.out_dir}/transcripts.tiled.features.tsv"
        feature_plain = f"{args.out_dir}/transcripts.tiled.features.hdr.tsv"
        if args.tile:
            prereq=[tile_flag]
        else:
            prereq=[feature_nohdr]

        cmds = cmd_separator([], f"Preparing a feature file for FICTURE using the feature file from pts2tiles: {feature_nohdr}")
        # plain file 
        cmds.append(f"(echo {args.colname_feature} {args.colname_count} | tr ' ' '\\t'; cat {feature_nohdr};) > {feature_plain}")

        mm.add_target(feature_plain, prereq, cmds)

        # update the args.in_feature
        args.in_feature = feature_plain

    # - segment
    if args.segment:
        scheck_app(args.gzip)
        assert args.width is not None, "When --segment, provide at least one hexagon width for segmentation in FICTURE-compatible format using --width"
        hexagon_widths=[int(x) for x in args.width.split(",")]

        for hexagon_width in hexagon_widths:
            hexagon_prefix=f"{args.out_dir}/hexagon.d_{hexagon_width}"
            cmds = cmd_separator([], f"Creating hexagon-indexed SGE in FICTURE-compatible format for {hexagon_width}um...")
            cmd = " ".join([
                ficture2bin, "tiles2hex",
                f"--in-tsv {args.out_dir}/transcripts.tiled.tsv",
                f"--in-index {args.out_dir}/transcripts.tiled.index",  
                f"--feature-dict {feature_nohdr}",
                f"--icol-x {args.colidx_x-1}",
                f"--icol-y {args.colidx_y-1}",
                f"--icol-feature 2",
                f"--icol-int 3",
                f"--out {hexagon_prefix}.tsv",
                f"--temp-dir {args.out_dir}/tmp/{hexagon_width}",
                f"--threads {args.threads}",
                f"--hex-grid-dist {hexagon_width}",
                f"--min-count {args.min_ct_per_unit_hexagon}"
                ])
            cmds.append(cmd)                
            cmds.append(f"{args.sort} -S {args.sort_mem} -k 1,1 {hexagon_prefix}.tsv > {hexagon_prefix}.randomized.tsv")
            cmds.append(f"rm -f {hexagon_prefix}.tsv")
            cmds.append(f"[ -f {hexagon_prefix}.randomized.tsv ] && [ -f {hexagon_prefix}.json ] && touch {hexagon_prefix}.done" )
            mm.add_target(f"{hexagon_prefix}.done", [f"{args.out_dir}/transcripts.tiled.done", feature_plain], cmds)

    # - segment 10x 
    if args.segment_10x:
        scheck_app(args.gzip)
        assert args.segment_width_10x is not None, "When --segment-10x, provide at least one hexagon width for segmentation in 10X MEX format using --segment-width-10x"
        hexagon_widths=[int(x) for x in args.segment_width_10x.split(",")]

        for hexagon_width in hexagon_widths:
            hexagon_10x_prefix=f"{args.out_dir}/hexagon.d_{hexagon_width}.min1"
            cmds = cmd_separator([], f"Creating hexagon-indexed SGE in 10X MEX format for {hexagon_width}um...")
            cmd = " ".join([
                ficture2bin, "tiles2hex",
                f"--in-tsv {args.out_dir}/transcripts.tiled.tsv",
                f"--in-index {args.out_dir}/transcripts.tiled.index",  
                f"--feature-dict {feature_nohdr}",
                f"--icol-x {args.colidx_x-1}",
                f"--icol-y {args.colidx_y-1}",
                f"--icol-feature 2",
                f"--icol-int 3",
                f"--out {hexagon_10x_prefix}.tsv",
                f"--temp-dir {args.out_dir}/tmp/{hexagon_width}",
                f"--threads {args.threads}",
                f"--hex-grid-dist {hexagon_width}",
                f"--min-count 1"
                ])
            cmds.append(cmd)                
            cmds.append(f"{args.sort} -S {args.sort_mem} -k 1,1 {hexagon_10x_prefix}.tsv > {hexagon_10x_prefix}.randomized.tsv")
            
            hexagon_10x_dir=f"{args.out_dir}/hexagon.d_{hexagon_width}.10x"
            os.makedirs(hexagon_10x_dir, exist_ok=True)
            
            cmd = " ".join([
                f"{args.spatula} sptsv2mex",
                f"--tsv {hexagon_10x_prefix}.randomized.tsv",
                f"--json {hexagon_10x_prefix}.json",
                f"--out-dir {hexagon_10x_dir}"])
            cmds.append(cmd)
            cmds.append(f'if [ -f {hexagon_10x_dir}/barcodes.tsv.gz ] && [ -f {hexagon_10x_dir}/features.tsv.gz ]  && [ -f {hexagon_10x_dir}/matrix.mtx.gz ] && [ -f {hexagon_10x_prefix}.randomized.tsv ]; then rm {hexagon_10x_prefix}.randomized.tsv {hexagon_10x_prefix}.tsv; fi')
            cmds.append(f"[ -f {hexagon_10x_dir}/barcodes.tsv.gz ] && [ -f {hexagon_10x_dir}/features.tsv.gz ]  && [ -f {hexagon_10x_dir}/matrix.mtx.gz ] && touch {hexagon_10x_dir}/hexagon_10x.done" )
            mm.add_target(f"{hexagon_10x_dir}/hexagon_10x.done", [f"{args.out_dir}/transcripts.tiled.done", feature_plain], cmds)

    # - lda
    if args.init_lda:
        lda_runs = define_lda_runs(args)
        for lda_params in lda_runs:
            # params & prefix
            train_width = lda_params["train_width"]
            n_factor = lda_params["n_factor"]

            model_id = lda_params["model_id"]
            model_prefix = os.path.join(args.out_dir, model_id)
            #lda_fillr = int(train_width // 2 + 1)

            # files
            hexagon = f"{args.out_dir}/hexagon.d_{train_width}.randomized.tsv"
            meta = f"{args.out_dir}/hexagon.d_{train_width}.json"
            lda_model_matrix = f"{model_prefix}.model.tsv"
            lda_fit_tsv = f"{model_prefix}.results.tsv"
            #lda_postcount_tsv = f"{model_prefix}.pseudobulk.tsv"
            lda_de = f"{model_prefix}.bulk_chisq.tsv"

            # 1) fit model
            cmds = cmd_separator([], f"LDA training for {train_width}um and {n_factor} factors...")
            cmd = " ".join([
                ficture2bin, "lda4hex",
                f"--in-data {hexagon}",
                f"--in-meta {meta}",
                f"--out-prefix {model_prefix}.unsorted",
                f"--n-topics {n_factor}",
                f"--transform",
                f"--min-count-train {args.min_ct_per_unit_train}",
                f"--min-count-per-feature {args.min_ct_per_feature}",
                f"--features {feature_nohdr}",
                f"--include-feature-regex '{args.include_feature_regex}'" if args.include_feature_regex is not None else "",
                f"--exclude-feature-regex '{args.exclude_feature_regex}'" if args.exclude_feature_regex is not None else "",
                f"--minibatch-size {args.minibatch_size}",
                f"--seed {args.seed}",
                f"--n-epochs {args.train_epoch}",
                f"--threads {args.threads}",
                ])
            cmds.append(cmd)
            #cmd = f"cut -f 2- {model_prefix}.unsorted.results.tsv > {model_prefix}.unsorted.results.nohex.tsv"
            #cmds.append(cmd)
            cmd = " ".join([
                args.spatula, "append-topk-tsv",
                f"--in-model {model_prefix}.unsorted.model.tsv",
                f"--in-json {meta}",
                f"--out-model {model_prefix}.model.tsv",
                f"--reorder",
                f"--in-tsv {model_prefix}.unsorted.results.tsv",
                f"--out-tsv {model_prefix}.results.tsv.gz",
                f"--offset-model 1"
            ])
            cmds.append(cmd)
            #cmds.append(f"rm -f {model_prefix}.unsorted.model.tsv {model_prefix}.unsorted.results.tsv {model_prefix}.unsorted.results.nohex.tsv")
            cmds.append(f"rm -f {model_prefix}.unsorted.model.tsv {model_prefix}.unsorted.results.tsv")
            cmds.append(f"[ -f {lda_fit_tsv}.gz ] && [ -f {lda_model_matrix} ] && touch {model_prefix}.done" )
            mm.add_target(f"{model_prefix}.done", [f"{args.out_dir}/transcripts.tiled.done", f"{args.out_dir}/hexagon.d_{train_width}.done", feature_plain], cmds)

            # create color table
            cmds = cmd_separator([], f"Generate the color map ")
            color_map=f"{model_prefix}.cmap.tsv"
            cmds.append(f'head -n $(({n_factor} + 1)) "{args.cmap_file}" > "{color_map}"')
            mm.add_target(color_map, [args.cmap_file], cmds)

            # 2) DE
            cmds = cmd_separator([], f" LDA DE/report for {train_width}um and {n_factor} factors...")
            cmds.append(f"{ficture2de} --input {lda_model_matrix} --output {lda_de} --feature_label Feature --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold}")
            cmd = " ".join([
                ficture2report,
                f"--de {lda_de}",
                f"--pseudobulk {lda_model_matrix}",
                f"--feature_label Feature",
                f"--color_table {color_map}",
                f"--output_pref {model_prefix}"
                ])
            cmds.append(cmd)
            cmds.append(f"[ -f {lda_de} ] && [ -f {model_prefix}.factor.info.html ] && touch {model_prefix}_summary.done")
            mm.add_target(f"{model_prefix}_summary.done", [f"{model_prefix}.done", color_map], cmds)

    # - decode
    if args.decode:
        scheck_app(args.sort)
        scheck_app(args.gzip)

        decode_runs = define_decode_runs(args)
        for decode_params in decode_runs:
            # input
            model_prefix = os.path.join(args.out_dir, decode_params["model_id"])
            model_path = decode_params["model_path"]
            color_map = decode_params["cmap_path"]
            
            # prerequisities
            #fit_prereq = decode_params["prerequisite_path"]
            
            # params & prefix
            fit_width = decode_params["fit_width"]
            decode_id = decode_params["decode_id"]
            decode_prefix = os.path.join(args.out_dir, decode_id)

            fit_n_move = int(fit_width / args.anchor_res)
            decode_postcount = f"{decode_prefix}.pseudobulk.tsv"
            decode_fit_tsv = f"{decode_prefix}.tsv"
            decode_flag = f"{decode_prefix}.done"

            decode_de = f"{decode_prefix}.bulk_chisq.tsv"
            decode_report = f"{decode_prefix}.factor.info.html"
            decode_summary_flag=f"{decode_prefix}_summary.done"

            #1) transform/fit
            cmds=cmd_separator([], f"Creating decode, ID: {decode_id}")
            cmd = " ".join([
                ficture2bin, "pixel-decode",
                f"--model {model_path}",
                f"--in-tsv {args.out_dir}/transcripts.tiled.tsv",
                f"--in-index {args.out_dir}/transcripts.tiled.index",
                f"--temp-dir {args.out_dir}/tmp/{decode_id}",
                f"--out {decode_prefix}.tsv",
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
            # compress the decode tsv file
            cmds.append(f"{args.gzip} -f {decode_fit_tsv}")
            cmds.append(f"{args.gzip} -f {decode_postcount}")

            cmds.append(f"[ -f {decode_fit_tsv}.gz ] && [ -f {decode_postcount}.gz ] && touch {decode_flag}" )
            mm.add_target(decode_flag, [f"{args.out_dir}/transcripts.tiled.done", f"{model_prefix}.done"], cmds)

            # 3) DE/report
            cmds=cmd_separator([], f"Decode DE and report, ID: {decode_id}")
            # - transform-DE
            cmd = " ".join([
                ficture2de,
                f"--input {decode_postcount}.gz", ## TBC can this use .gz file?
                f"--output {decode_de}",
                f"--min_ct_per_feature {args.min_ct_per_feature}",
                f"--max_pval_output {args.de_max_pval}",
                f"--min_fold_output {args.de_min_fold}",
                f"--feature_label Feature"
                ])
            cmds.append(cmd)
            # - transform-report
            cmd = " ".join([
                ficture2report,
                f"--de {decode_de}",
                f"--pseudobulk {decode_postcount}.gz",
                f"--feature_label Feature",
                f"--color_table {color_map}",
                f"--output_pref {decode_prefix}"
                ])
            cmds.append(cmd)
            # - done & target
            cmds.append(f"[ -f {decode_de} ] && [ -f {decode_report} ] && touch {decode_summary_flag}")
            mm.add_target(decode_summary_flag, [f"{decode_prefix}.done", color_map], cmds)

            # 7) visualization
            cmds=cmd_separator([], f"Decode visualization, ID: {decode_id}")
            cmd = " ".join([
                f"{args.gzip} -dc {decode_fit_tsv}.gz |",
                ficture2bin, "draw-pixel-factors",
                f"--in-tsv /dev/stdin",
                f"--header-json {decode_prefix}.json",
                f"--in-color {color_map}",
                f"--out {decode_prefix}.png",
                f"--scale {args.decode_scale}",
                f"--range {args.in_minmax}"
                ])
            cmds.append(cmd)
            mm.add_target(f"{decode_prefix}.png", [decode_summary_flag, color_map, minmax_prereq], cmds)

    # - summary (update: always run)
    # if args.summary:
    prerequisities=[feature_plain] 
    summary_aux_args=[]
    # lda or external model
    if args.init_lda:
        summary_aux_args_models = ["--lda-model"]
        train_params = define_lda_runs(args)
        for train_param in train_params:
            train_width = train_param["train_width"]
            n_factor = train_param["n_factor"]
            model_prefix = os.path.join(args.out_dir, train_param["model_id"])
            # prerequisities
            prerequisities.append(f"{model_prefix}.done")
            # args
            summary_cmap = f"{model_prefix}.cmap.tsv"
            summary_aux_args_models.append(f"lda,{model_prefix}.model.tsv,{train_param['model_id']},{train_width},{n_factor},{summary_cmap}")
        summary_aux_args.append(" ".join(summary_aux_args_models))
    # projection & decode
    if args.decode:
        summary_aux_args_decode = ["--decode"]
        decode_runs = define_decode_runs(args)
        for decode_params in decode_runs:
            model_type = decode_params["model_type"]
            model_id = decode_params["model_id"]
            fit_width = decode_params["fit_width"]
            decode_id = decode_params["decode_id"]
            # prerequisities
            if args.decode:
                prerequisities.append(f"{args.out_dir}/{decode_id}.done")
            # args
            if args.decode:
                summary_aux_args_decode.append(f"{model_type},{model_id},{decode_id},{fit_width},{args.anchor_res}")
        if args.decode and len(summary_aux_args_decode) > 1:
            summary_aux_args.append(" ".join(summary_aux_args_decode))
    # summary
    cmds = cmd_separator([], f"Summarizing output into to the <out_json> files...")
    cmd = " ".join([
        "cartloader", "write_json_for_ficture2",
            "--merge --merge-override",
            f"--in-transcript {args.in_transcript}",
            f"--in-feature {args.in_feature}", # use the original feature file for SGE
            f"--in-feature-ficture {feature_plain}",
            f"--in-minmax {args.in_minmax}",
            f"--out-json {args.out_json}",
            " ".join(summary_aux_args)
        ])
    cmds.append(cmd)
    prerequisities.append(minmax_prereq)
    mm.add_target(args.out_json, prerequisities, cmds)

    ## write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that ast least one run option was turned on")
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

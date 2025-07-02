import sys, os, gzip, argparse, logging, shutil, subprocess
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, flexopen

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader run_ficture2", description="Run FICTURE2")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--threads', type=int, default=8, help='Maximum number of threads to use in each process')
    run_params.add_argument('--n-jobs', type=int, default=2, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default="run_ficture2.mk", help='The name of the Makefile to generate (default: run_ficture2.mk)')

    cmd_params = parser.add_argument_group("Commands", "FICTURE commands to run together")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run the main functions (tile, segment, init-lda, decode, summary)')
    cmd_params.add_argument('--tile', action='store_true', default=False, help='(Main function) Perform tiling step')
    cmd_params.add_argument('--segment', action='store_true', default=False, help='(Main function) Perform hexagon segmentation into FICTURE-compatible format')
    cmd_params.add_argument('--init-lda', action='store_true', default=False, help='(Main function) Initialize model with LDA model training')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='(Main function) Perform pixel-level decoding')
    cmd_params.add_argument('--summary', action='store_true', default=False, help='(Main function) Generate a JSON file summarizing all fixture parameters for which outputs are available in the <out-dir>.')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output parameters for FICTURE")
    inout_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    inout_params.add_argument('--out-json', type=str, default=None, help="Output JSON file for summarizing the ficture parameters (default: <out-dir>/ficture.params.json)")
    inout_params.add_argument('--in-transcript', type=str, default=None, help='Path to the input unsorted transcript-indexed SGE file in TSV format (default: <out-dir>/transcripts.unsorted.tsv.gz)')
    inout_params.add_argument('--in-minmax', type=str, default=None, help='Path to the input coordinate minmax TSV file. (default: <out-dir>/coordinate_minmax.tsv)')
    inout_params.add_argument('--in-feature', type=str, default=None,  help='Path to the input UMI count per gene TSV file.(default: feature.clean.tsv.gz).')
    inout_params.add_argument('--in-feature-ficture', type=str, default=None, help='(Optional) Use --in-feature-ficture to provide a feature file for FICTURE analysis if such feature file exists. Alternatively, see "Feature Customizing Auxiliary Parameters" to customize the features in FICTURE analysis.')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in um) for LDA training (default: None)')
    key_params.add_argument('--n-factor', type=str, default=None, help='Comma-separated list of factor counts for LDA training.')
    key_params.add_argument('--anchor-res', type=int, default=6, help='Anchor resolution for decoding (default: 6)')
    key_params.add_argument('--radius-buffer', type=int, default=1, help='Buffer to radius(=anchor_res + radius_buffer) for pixel-level decoding (default: 1)')
    key_params.add_argument('--decode-scale', type=int, default=1, help='scales input coordinates to pixels in the output image (default: 1)')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default="1G", help='Memory size for each process')
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",
    env_params.add_argument('--ficture2', type=str, default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "submodules", "punkst"),  help='Path to punkst(ficture2) repositor (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--python', type=str, default="python3",  help='Python3 binary')

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Customizing Auxiliary Parameters",
                                                      "Auxiliary parameters for customizing features in FICTURE analysis without modifying the original input feature TSV file. This ensures the original feature TSV file is retained in the output JSON file for downstream processing .")
    # given the input sge should be standardized, the csv-delim, csv-colname-feature-name, ftr-delim, ftr-colname-feature-name are not necessary
    aux_ftrfilter_params.add_argument('--filter-by-overlapping-features', action='store_true', default=False, help='When the input SGE is stitched SGE, it is optional to filter the features in FICTURE analysis by only shared features')
    aux_ftrfilter_params.add_argument('--in-feature-dist', type=str, default=None, help='Path to the input feature distribution file. This file is used to identify overlapping features for FICTURE analysis. (default: None)')
    aux_ftrfilter_params.add_argument('--min-ct-per-ftr-tile', type=int, default=0, help='Apply a minimum count to filter overlapping feature. Filtering process will be applied if --min-ct-per-overlapftr > 0. (default: 0)')
    # aux_ftrfilter_params.add_argument('--out-feature-ficture', type=str, default="features.ficture.tsv.gz", help='File name for the output TSV file of feature used in FICTURE analysis (default: None)')
    # aux_ftrfilter_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    # aux_ftrfilter_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')
    # type regex
    # aux_ftrfilter_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). When --include-feature-type-regex, use --colname-feature-type or --feature-type-ref to provide gene type information.') # (e.g. protein_coding|lncRNA)
    # aux_ftrfilter_params.add_argument('--colname-feature-type', type=str, default=None, help='Column name in the --in-transcript that has gene type information (default: None). ')
    # aux_ftrfilter_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated reference file to provide gene type information for each each per row (default: None)')
    # aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='Column index for gene name in the reference file (default: None).')
    # aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='Column index for gene type in the reference file (default: None).')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--colname-count', type=str, default="count", help='Columns from the input transcript file to be used as key for count')
    aux_params.add_argument('--colname-feature', type=str, default="gene", help='Columns from the input transcript file to be used as key for feature name')
    aux_params.add_argument('--tile-size', type=int, default=500, help='Tile size for tiling (default: 500)')
    aux_params.add_argument('--tile-buffer', type=int, default=1000, help='Tile buffer for tiling (default: 1000)')
    aux_params.add_argument('--seed', type=int, default=1, help='Random seed for random number generation (default: 0)')
    # segmentation - ficture
    aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=50, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format (default: 50)')
    # minibatch
    aux_params.add_argument('--minibatch-size', type=int, default=500, help='Batch size used in minibatch processing (default: 500)')
    #aux_params.add_argument('--minibatch-buffer', type=int, default=30, help='Batch buffer used in minibatch processing (default: 30)')
    # train
    aux_params.add_argument('--min-ct-per-unit-train', type=int, default=50, help='Minimum count for training (default: 50)')
    aux_params.add_argument('--train-epoch', type=int, default=2, help='Training epoch for LDA model (default: 2)')
    #aux_params.add_argument('--train-epoch-id-len', type=int, default=2, help='Training epoch ID length (default: 2)')
    #aux_params.add_argument('--lda-rand-init', type=int, default=10, help='Number of random initialization during model training (default: 10)')
    #aux_params.add_argument('--lda-plot-um-per-pixel', type=float, default=1, help='Image resolution for LDA plot (default: 1)')
    # fit
    aux_params.add_argument('--fit-width',  type=str, default=None, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    #aux_params.add_argument('--fit-precision', type=float, default=2, help='Output precision of model fitting (default: 2)')
    aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=50, help='Minimum count per hexagon unit during model fitting (default: 20)')
    aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot (default: 1)')   # in Scopeflow, this is set to 2
    # decode
    #aux_params.add_argument('--decode-top-k', type=int, default=3, help='Top K columns to output in pixel-level decoding results (default: 3)')
    #aux_params.add_argument('--decode-block-size', type=int, default=100, help='Block size for pixel decoding output (default: 100)')
    #aux_params.add_argument('--decode-scale', type=int, default=100, help='Scale parameters for pixel decoding output (default: 100)')
    #aux_params.add_argument('--decode-precision', type=float, default=0.01, help='Precision of pixel level decoding (default: 0.01)')
    #aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=0.5, help='Image resolution for pixel decoding plot (default: 0.5)')
    # merge_by_pixel
    # aux_params.add_argument('--merge-max-dist-um', type=float, default=0.1, help='Maximum distance in um for merging pixel-level decoding results (default: 0.1)')
    # aux_params.add_argument('--merge-max-k', type=int, default=1, help='Maximum number of K columns to output in merged pixel-level decoding results (default: 1)')
    # aux_params.add_argument('--merge-max-p', type=int, default=1, help='Maximum number of P columns to output in merged pixel-level decoding results (default: 1)')
    # color map
    aux_params.add_argument('--cmap-file', type=str, default=None, help='Define the path to the fixed color map (default: <cartloader_dir>/assets/fixed_color_map_60.tsv)')
    # others parameters shared across steps
    aux_params.add_argument('--min-ct-per-feature', type=int, default=20, help='Minimum count per feature during LDA training, transform and decoding (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='p-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def define_lda_runs(args):
    assert args.init_lda, "--init-lda must be ON when running define_lda_runs()"
    assert args.width is not None, "When --init-lda is ON, provide at least one train width for LDA training using --train-width"

    if args.n_factor is None:
        print("Warning: --n-factor is not provided for LDA. Using default values: 12,24")
        args.n_factor = "12,24"
    #assert args.n_factor is not None, "When --init-lda is ON, provide at least one n.factor for LDA training using --n-factor"

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
        args.summary = True

    #assert not (args.init_lda and args.init_ext), "Cannot choose both --init-lda and --init-ext"

    # input/output/other files
    # dirs
    os.makedirs(args.out_dir, exist_ok=True)
    #os.makedirs(os.path.join(args.out_dir, "tmp"), exist_ok=True)
    
    # start mm
    mm = minimake()

    # in files
    if args.in_transcript is None:
        args.in_transcript = os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")

    # do not require to have the feature and minmax files        
    # if args.in_minmax is None:
    #     args.in_minmax = os.path.join(args.out_dir, "coordinate_minmax.tsv")
    
    # if args.in_feature is None:
    #     args.in_feature = os.path.join(args.out_dir, "feature.clean.tsv.gz")
    
    assert os.path.exists(args.in_transcript), "Provide at least one valid input transcript-indexed SGE file by --in-transcript or --in-cstranscript"
    if args.in_minmax is not None:
        assert os.path.exists(args.in_minmax), "Provide a valid input coordinate minmax file by --in-minmax, or skip specifying it"
    if args.in_feature is not None:
        assert os.path.exists(args.in_feature), "Provide a valid input feature file by --in-feature, or skip specifying it"
    
    if args.cmap_file is None:
        script_path = os.path.abspath(__file__)
        cartloader_root = os.path.abspath(os.path.join(script_path, "..", "..", ".."))
        args.cmap_file=os.path.join(cartloader_root, "assets", "fixed_color_map_256.tsv")
    
    print(f"cmap path: {args.cmap_file}")
    if not os.path.exists(args.cmap_file):
        raise FileNotFoundError(f"Color map not found at: {args.cmap_file}")

    assert os.path.exists(args.ficture2), f"Provide a valid path to the FICTURE2 repository by --ficture2, or set it to the default path: {args.ficture2}"
    ficture2bin = os.path.join(args.ficture2, "bin/punkst")
    ficture2de = args.python + " " + os.path.join(args.ficture2, "ext/py/de_bulk.py")
    ficture2report = args.python + " " + os.path.join(args.ficture2, "ext/py/factor_report.py")
    
    # out files
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.params.json")

    # 1. tiling :
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

        if args.in_minmax is None: ## specify minmax file if not provided
            args.in_minmax = os.path.join(args.out_dir, "transcripts.tiled.coord_range.tsv")
        
        # feature file
        if args.in_feature_ficture is None:
            in_feature_ficture = args.in_feature_ficture
        else:
            in_feature_ficture = args.in_feature
        
        ## Features step 1. write the feature file from the tiled SGE
        if in_feature_ficture is not None:
            ## if a specific feature file is provided, use it as "selected features" for FICTURE analysis
            in_feature_ficture_flag = in_feature_ficture
            feature_plain = f"{args.out_dir}/transcripts.tiled.selected_features.hdr.tsv"
            feature_nohdr = f"{args.out_dir}/transcripts.tiled.selected_features.tsv"
            with flexopen(in_feature_ficture, "rt") as f:
                with open(feature_plain, "wt") as wf:
                    with open(feature_nohdr, "wt") as wf2:
                        wf.write(f.readline())
                        for line in f:
                            toks = line.strip().split("\t")
                            count = int(toks[-1])
                            if count >= args.min_ct_per_feature:
                                wf.write(line)
                                wf2.write(line)
            feature_nohdr_flag = feature_nohdr
        else: 
            ## if feature file is not provided, create one from the input transcript file
            in_feature_ficture_flag = f"{args.out_dir}/transcripts.tiled.done"
            feature_plain = f"{args.out_dir}/transcripts.tiled.features.hdr.tsv"
            feature_nohdr = f"{args.out_dir}/transcripts.tiled.features.tsv"
            feature_nohdr_flag = f"{args.out_dir}/transcripts.tiled.done"
            in_feature_ficture = feature_plain
            args.in_feature = feature_plain
            cmds.append(f"(echo {args.colname_feature} {args.colname_count} | tr ' ' '\\t'; cat {feature_nohdr};) > {feature_plain}")

        cmds.append(f"[ -f {args.out_dir}/transcripts.tiled.tsv ] && [ -f {args.out_dir}/transcripts.tiled.index ] && touch {args.out_dir}/transcripts.tiled.done" )
        mm.add_target(f"{args.out_dir}/transcripts.tiled.done", [args.in_transcript], cmds)

        #  Features step 2. generate overlapping features for FICTURE analysis (optional only if the input SGE is stitched SGE)
        if args.filter_by_overlapping_features:
            cmds = cmd_separator([], f"Customizing features for FICTURE analysis: limited to shared features and features with a minimal count in the stitched SGE...")
            feature_overlapping_plain = os.path.join(args.out_dir, "transcripts.tiled.overlapping_features.hdr.tsv") if args.min_ct_per_ftr_tile == 0 else os.path.join(args.out_dir, f"transcripts.tiled.overlapping_features.min{args.min_ct_per_ftr_tile}.hdr.tsv")
            cmd = " ".join(["cartloader feature_overlapping",
                                        f"--in-dist {args.in_feature_dist}", 
                                        f"--in-feature {in_feature_ficture}",
                                        f"--output {feature_overlapping_plain}", 
                                        f"--min-ct-per-ftr-tile {args.min_ct_per_ftr_tile}",
                                        f"--colname-count {args.colname_count}",
                                        f"--colname-feature-name {args.colname_feature}",
                                        f"--log"
                                        ])
            cmds.append(cmd)
            mm.add_target(feature_overlapping_plain, [args.in_feature_dist, in_feature_ficture_flag], cmds)
            
            feature_overlapping_nohdr = os.path.join(args.out_dir, "transcripts.tiled.overlapping_features.tsv") if args.min_ct_per_ftr_tile == 0 else os.path.join(args.out_dir, f"transcripts.tiled.overlapping_features.min{args.min_ct_per_ftr_tile}.tsv")
            cmds = cmd_separator([], f"Generating feature without header for overlapping features...")
            cmds.append(f"cat {feature_overlapping_plain} | tail +2 > {feature_overlapping_nohdr}")
            mm.add_target(feature_overlapping_nohdr, [feature_overlapping_plain], cmds)
            
            feature_plain = feature_overlapping_plain
            feature_nohdr = feature_overlapping_nohdr
            feature_nohdr_flag = feature_overlapping_nohdr

    # 2. segment
    if args.segment:
        scheck_app(args.gzip)
        assert args.width is not None, "When --segment, provide at least one hexagon width for segmentation in FICTURE-compatible format using --hexagon-width"
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
            mm.add_target(f"{hexagon_prefix}.done", [f"{args.out_dir}/transcripts.tiled.done", feature_nohdr_flag], cmds)

    # 3. lda
    if args.init_lda:
        lda_runs = define_lda_runs(args)
        for lda_params in lda_runs:
            # params & prefix
            train_width = lda_params["train_width"]
            n_factor = lda_params["n_factor"]

            model_id = lda_params["model_id"]
            model_prefix = os.path.join(args.out_dir, model_id)

            lda_fillr = int(train_width // 2 + 1)
            # files
            hexagon = f"{args.out_dir}/hexagon.d_{train_width}.randomized.tsv"
            meta = f"{args.out_dir}/hexagon.d_{train_width}.json"
            lda_model_matrix = f"{model_prefix}.model.tsv"
            lda_fit_tsv = f"{model_prefix}.results.tsv"
#           lda_postcount_tsv = f"{model_prefix}.pseudobulk.tsv"
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
            mm.add_target(f"{model_prefix}.done", [f"{args.out_dir}/transcripts.tiled.done", f"{args.out_dir}/hexagon.d_{train_width}.done", feature_nohdr_flag], cmds)

            # create color table
            out_cmap = f"{model_prefix}.cmap.tsv"
            with open(args.cmap_file, "r") as f:
                with open(out_cmap, "w") as f2:
                    for i in range(n_factor+1):
                        f2.write(f.readline())

            # 2) DE
            cmds = cmd_separator([], f" LDA DE/report for {train_width}um and {n_factor} factors...")
            cmds.append(f"{ficture2de} --input {lda_model_matrix} --output {lda_de} --feature_label Feature --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold}")
            cmd = " ".join([
                ficture2report,
                f"--de {lda_de}",
                f"--pseudobulk {lda_model_matrix}",
                f"--feature_label Feature",
                f"--color_table {out_cmap}",
                f"--output_pref {model_prefix}"
                ])
            cmds.append(cmd)
            cmds.append(f"[ -f {lda_de} ] && [ -f {model_prefix}.cmap.tsv ] && [ -f {model_prefix}.factor.info.html ] && touch {model_prefix}_summary.done")
            #cmds.append(f"[ -f {lda_de} ] && touch {model_prefix}_summary.done")
            mm.add_target(f"{model_prefix}_summary.done", [f"{model_prefix}.done"], cmds)

    if args.decode:
        scheck_app(args.sort)
        tiled_tsv = f"{args.out_dir}/transcript.tiled.tsv"
        decode_runs = define_decode_runs(args)
        for decode_params in decode_runs:
            # input
            model_prefix = os.path.join(args.out_dir, decode_params["model_id"])
            model_path = decode_params["model_path"]
            cmap_path = decode_params["cmap_path"]
            # prerequisities
            fit_prereq = decode_params["prerequisite_path"]
            # params & prefix
            fit_width = decode_params["fit_width"]
            decode_id = decode_params["decode_id"]
            decode_prefix = os.path.join(args.out_dir, decode_id)

            fit_n_move = int(fit_width / args.anchor_res)
            decode_postcount = f"{decode_prefix}.pseudobulk.tsv"
            decode_fit_tsv = f"{decode_prefix}.tsv"
            decode_de = f"{decode_prefix}.bulk_chisq.tsv"
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
            cmds.append(f"[ -f {decode_fit_tsv} ] && [ -f {decode_postcount} ] && touch {decode_prefix}.done" )
            mm.add_target(f"{decode_prefix}.done", [f"{args.out_dir}/transcripts.tiled.done", f"{model_prefix}.done"], cmds)

            # 3) DE/report
            cmds=cmd_separator([], f"Decode DE and report, ID: {decode_id}")
            # - transform-DE
            cmd = " ".join([
                ficture2de,
                f"--input {decode_postcount}",
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
                f"--pseudobulk {decode_postcount}",
                f"--feature_label Feature",
                f"--color_table {cmap_path}",
                f"--output_pref {decode_prefix}"
                ])
            cmds.append(cmd)
            # compress the decode tsv file
            cmd = " ".join([
                args.gzip, "-f", decode_fit_tsv
            ])
            cmds.append(cmd)
            # - done & target
            cmds.append(f"[ -f {decode_de} ] && [ -f {decode_prefix}.factor.info.html ] && [ -f {decode_fit_tsv}.gz ] && touch {decode_prefix}_summary.done")
            mm.add_target(f"{decode_prefix}_summary.done", [f"{decode_prefix}.done", cmap_path], cmds)

            # 4) visualization
            cmds=cmd_separator([], f"Decode visualization, ID: {decode_id}")
            cmd = " ".join([
                f"{args.gzip} -dc {decode_fit_tsv}.gz |",
                ficture2bin, "draw-pixel-factors",
                f"--in-tsv /dev/stdin",
                f"--header-json {decode_prefix}.json",
                f"--in-color {cmap_path}",
                f"--out {decode_prefix}.png",
                f"--scale {args.decode_scale}",
                f"--range {args.in_minmax}"
                ])
            cmds.append(cmd)
            mm.add_target(f"{decode_prefix}.png", [f"{decode_prefix}_summary.done", cmap_path], cmds)

    if args.summary:
        prerequisities=[feature_nohdr_flag] # since feature_nohdr and feature_plain are generated at the same step, use feature_nohdr_flag as the prerequisite for feature_plain
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
                "--merge",
                f"--in-transcript {args.in_transcript}",
                f"--in-feature {args.in_feature}", # use the original feature file for SGE
                f"--in-feature-ficture {feature_plain}",
                f"--in-minmax {args.in_minmax}",
                f"--out-dir {args.out_dir}",
                f"--out-json {args.out_json}",
                " ".join(summary_aux_args)
            ])
        cmds.append(cmd)
        mm.add_target(args.out_json, prerequisities, cmds)

    ## write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that ast least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    if args.dry_run:
        dry_cmd=f"make -f {make_f} -n {'-B' if args.restart else ''} "
        os.system(dry_cmd)
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}")
    else:
        exe_cmd=f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}"
        result = subprocess.run(exe_cmd, shell=True)
        if result.returncode != 0:
            print(f"Error in executing: {exe_cmd}")
            sys.exit(1)

if __name__ == "__main__":
    # Get the path to the cartloader repository
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

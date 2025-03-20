import sys, os, gzip, argparse, logging, shutil, subprocess
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, find_major_axis, add_param_to_cmd


def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader run_ficture", description="Run FICTURE")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default="run_ficture.mk", help='The name of the Makefile to generate (default: run_ficture.mk)')

    cmd_params = parser.add_argument_group("Commands", "FICTURE commands to run together")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run the main functions (sorttsv, minibatch, segment, lda, decode, summary)')
    cmd_params.add_argument('--main-ext', action='store_true', default=False, help='Run the main functions with external projections (sorttsv, minibatch, segment, ext-lda, decode, summary)')
    cmd_params.add_argument('--sorttsv', action='store_true', default=False, help='(Main function) Sort the input tsv file')
    cmd_params.add_argument('--minibatch', action='store_true', default=False, help='(Main function) Perform minibatch step')
    cmd_params.add_argument('--segment', action='store_true', default=False, help='(Main function) Perform hexagon segmentation into FICTURE-compatible format')
    cmd_params.add_argument('--init-lda', action='store_true', default=False, help='(Main function) Initialize model with LDA model training')
    cmd_params.add_argument('--init-ext', action='store_true', default=False, help='(Main-ext function) Initialize model from external reference')
    cmd_params.add_argument('--copy-ext-model', action='store_true', default=False, help='(Main-ext function) When running with --init-ext, copy the external model instead of running --init-from-pseudobulk. This requires that the gene lists match exactly')
    cmd_params.add_argument('--projection', action='store_true', default=False, help='(Main function) Perform projection/Transform')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='(Main function) Perform pixel-level decoding')
    cmd_params.add_argument('--summary', action='store_true', default=False, help='(Main function) Generate a JSON file summarizing all fixture parameters for which outputs are available in the <out-dir>.')
    cmd_params.add_argument('--segment-10x', action='store_true', default=False, help='(Additional function) Perform hexagon segmentation into 10x Genomics format')
    cmd_params.add_argument('--viz-per-factor', action='store_true', default=False, help='(Additional function) Generate pixel-level visualization for each factor')
    cmd_params.add_argument('--viz-dotplot', action='store_true', default=False, help='(Additional function) Generate dotplot visualization')
    cmd_params.add_argument('--skip-coarse-report', action='store_true', default=False, help='(Optional) Skip visualization and report generation for init-lda/init-ext and projection steps')

    # three types of external model
    #   - existing lda model: has tw & nf, and all supplement files (such as de) for the model available
    #   - models trained from hexagon-indexed SGE using other tools such as Seurat: has tw & nf without supplement files
    #   - custom models, such as those single-cell RNA-seq data: has nf without tw, and none of the supplement files

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output parameters for FICTURE")
    inout_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    inout_params.add_argument('--out-json', type=str, default=None, help="Output JSON file for summarizing the ficture parameters (default: <out-dir>/ficture.params.json)")
    inout_params.add_argument('--in-transcript', type=str, default=None, help='Path to the input unsorted transcript-indexed SGE file in TSV format (default: <out-dir>/transcripts.unsorted.tsv.gz)')
    inout_params.add_argument('--in-minmax', type=str, default=None, help='Path to the input coordinate minmax TSV file. (default: <out-dir>/coordinate_minmax.tsv)')  
    inout_params.add_argument('--in-feature', type=str, default=None,  help='Path to the input UMI count per gene TSV file.(default: feature.clean.tsv.gz).')
    inout_params.add_argument('--in-cstranscript', type=str, default=None, help='(Optional) Use --in-cstranscript if a coordinate-sorted transcript-indexed SGE file already exists to skip sorting, or to define a custom sorted file name (default: <out-dir>/transcripts.sorted.tsv.gz)')
    inout_params.add_argument('--in-feature-ficture', type=str, default=None, help='(Optional) Use --in-feature-ficture if a feature file for FICTURE analysis (default: <out-dir>/features.ficture.tsv.gz)')

    external_params = parser.add_argument_group("External Model Parameters", """When using --init-ext, provide the pre-trained model, its type, and ID by --ext-path, --ext-type and --ext-id. """)
    external_params.add_argument('--ext-path', type=str, default=None, help='Path to the external model file')
    external_params.add_argument('--ext-type', type=str, default="pseudobulk", help='The type of external model (options: "pseudobulk", "cluster")')
    external_params.add_argument('--ext-id', type=str, default=None, help='ID of the external model') 
    
    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--major-axis', type=str, default="X", choices=['X', 'Y', 'auto'], help='Axis where --in-transcript is sorted. Options: X, Y, auto. If "auto", use the longer axis (default: X)')
    key_params.add_argument('--mu-scale', type=float, default=1.0, help='Scale factor for mu, i.e., pixels per um (default: 1.0)')
    key_params.add_argument('--hexagon-width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in µm) used for creating hexagon-indexed SGE (default: same to train-width)')
    key_params.add_argument('--hexagon-width-10x', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in µm) used for creating hexagon-indexed SGE in 10x Genomics format (default: same to train-width)')
    key_params.add_argument('--train-width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in um) for LDA training (default: None)')
    key_params.add_argument('--n-factor', type=str, default=None, help='Comma-separated list of factor counts for LDA training. Incompatible with --init-ext. Defaults to 12,24 for --init-lda.')
    key_params.add_argument('--anchor-res', type=int, default=6, help='Anchor resolution for decoding (default: 6)')
    key_params.add_argument('--radius-buffer', type=int, default=1, help='Buffer to radius(=anchor_res + radius_buffer) for pixel-level decoding (default: 1)')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--bgzip', type=str, default="bgzip", help='Path to bgzip binary. For faster processing, use "bgzip -@ 4')
    env_params.add_argument('--tabix', type=str, default="tabix", help='Path to tabix binary')
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default="5G", help='Memory size for each process')

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Customizing Auxiliary Parameters", 
                                                    "Auxiliary parameters for customizing features in FICTURE analysis without modifying the original input feature TSV file. This ensures the original feature TSV file is retained in the output JSON file for downstream processing .")
    # given the input sge should be standardized, the csv-delim, csv-colname-feature-name, ftr-delim, ftr-colname-feature-name are not necessary
    aux_ftrfilter_params.add_argument('--out-ficture-feature', type=str, default="features.ficture.tsv.gz", help='File name for the output TSV file of feature used in FICTURE analysis (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-substr', type=str, default=None, help='A substring of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-substr', type=str, default=None, help='A substring of feature/gene names to be excluded (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')
    # type regex
    aux_ftrfilter_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). When --include-feature-type-regex, use --colname-feature-type or --feature-type-ref to provide gene type information.') # (e.g. protein_coding|lncRNA)
    aux_ftrfilter_params.add_argument('--colname-feature-type', type=str, default=None, help='Column name in the --in-transcript that has gene type information (default: None). ')
    aux_ftrfilter_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated reference file to provide gene type information for each each per row (default: None)')
    aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='Column index for gene name in the reference file (default: None).')
    aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='Column index for gene type in the reference file (default: None).')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--colname-count', type=str, default="count", help='Columns from the input transcript file to be used as key')
    # segmentation - ficture
    aux_params.add_argument('--hexagon-n-move', type=int, default=1, help='Level of hexagonal sliding when creating hexagon-indexed SGE in FICTURE compatible format (default: 1)')
    aux_params.add_argument('--hexagon-precision', type=float, default=2, help='Output precision of hexagon coordinates for FICTURE compatible format (default: 2)')
    aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=50, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format (default: 50)')
    # segmentation - 10x
    aux_params.add_argument('--hexagon-n-move-10x', type=int, default=1, help='Level of hexagonal sliding when creating hexagon-indexed SGE in 10x Genomics format (default: 1)')
    aux_params.add_argument('--hexagon-precision-10x', type=float, default=2, help='Output precision of hexagon coordinates for 10x Genomics format (default: 2)')
    aux_params.add_argument('--min-ct-per-unit-hexagon-10x', type=int, default=1, help='Minimum count per hexagon in hexagon segmentation in 10x Genomics format (default: 1)')
    # minibatch
    aux_params.add_argument('--minibatch-size', type=int, default=500, help='Batch size used in minibatch processing (default: 500)')
    aux_params.add_argument('--minibatch-buffer', type=int, default=30, help='Batch buffer used in minibatch processing (default: 30)')
    # train 
    aux_params.add_argument('--train-epoch', type=int, default=3, help='Training epoch for LDA model (default: 3)')
    aux_params.add_argument('--train-epoch-id-len', type=int, default=2, help='Training epoch ID length (default: 2)')
    aux_params.add_argument('--lda-rand-init', type=int, default=10, help='Number of random initialization during model training (default: 10)')
    aux_params.add_argument('--lda-plot-um-per-pixel', type=float, default=1, help='Image resolution for LDA plot (default: 1)')
    # fit 
    aux_params.add_argument('--fit-width',  type=str, default=None, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    aux_params.add_argument('--fit-precision', type=float, default=2, help='Output precision of model fitting (default: 2)')
    aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=20, help='Minimum count per hexagon unit during model fitting (default: 20)')
    aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot (default: 1)')   # in Scopeflow, this is set to 2
    # decode
    aux_params.add_argument('--decode-top-k', type=int, default=3, help='Top K columns to output in pixel-level decoding results (default: 3)')
    aux_params.add_argument('--decode-block-size', type=int, default=100, help='Block size for pixel decoding output (default: 100)')
    aux_params.add_argument('--decode-scale', type=int, default=100, help='Scale parameters for pixel decoding output (default: 100)')
    aux_params.add_argument('--decode-precision', type=float, default=0.01, help='Precision of pixel level decoding (default: 0.01)')
    aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=0.5, help='Image resolution for pixel decoding plot (default: 0.5)')
    # merge_by_pixel
    aux_params.add_argument('--merge-max-dist-um', type=float, default=0.1, help='Maximum distance in um for merging pixel-level decoding results (default: 0.1)') 
    aux_params.add_argument('--merge-max-k', type=int, default=1, help='Maximum number of K columns to output in merged pixel-level decoding results (default: 1)')
    aux_params.add_argument('--merge-max-p', type=int, default=1, help='Maximum number of P columns to output in merged pixel-level decoding results (default: 1)')
    # color map
    aux_params.add_argument('--cmap-name', type=str, default="turbo", help='Name of color map to generate a color map (default: turbo)')
    aux_params.add_argument('--cmap-static', action='store_true', default=False, help='Use a fixed color map for factor visualization')
    aux_params.add_argument('--static-cmap-file', type=str, help='If --cmap-static, define the path to the fixed color map (default: <cartloader_dir>/assets/fixed_color_map_60.tsv)')
    # others parameters shared across steps
    aux_params.add_argument('--min-ct-per-feature', type=int, default=20, help='Minimum count per feature during LDA training, transform and decoding (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='p-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def define_major_axis(args):
    if args.major_axis == "auto":
        args.major_axis = find_major_axis(args.in_minmax, "row")
    return args.major_axis

def run_choose_color(args, fit_tsv, n_factor, prefix):
    if args.cmap_static:
        cmds=[]
    else:
        cmds = cmd_separator([], f" Generate a color map for {n_factor} factors...")
        cmds.append(f"ficture choose_color --input {fit_tsv} --output {prefix} --cmap_name {args.cmap_name}")
    return cmds

def define_cmap(args, model_id):
    model_cmap = os.path.join(args.out_dir, f"{model_id}.rgb.tsv")
    if args.cmap_static:
        cmap_path = args.static_cmap_file 
    else:
        cmap_path = model_cmap
    return cmap_path

def define_lda_runs(args):
    assert args.init_lda, "--init-lda must be ON when running define_lda_runs()"
    assert args.train_width is not None, "When --init-lda is ON, provide at least one train width for LDA training using --train-width"

    if args.n_factor is None:
        print("Warning: --n-factor is not provided for LDA. Using default values: 12,24")
        args.n_factor = "12,24"
    #assert args.n_factor is not None, "When --init-lda is ON, provide at least one n.factor for LDA training using --n-factor"

    train_widths = [int(x) for x in args.train_width.split(",")] if args.train_width else [] #and not (args.use_external_model and args.external_model_type == "custom") else []
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

def define_ext_runs(args):
    assert args.init_ext, "--init-ext must be ON when running define_ext_runs()"
    assert args.train_width is not None, "When --init-ext is ON, provide at least one train width for LDA training using --train-width"
    assert args.n_factor is None, "When --init-ext is ON, --n-factors must NOT be used"
    assert args.ext_path is not None, "When --init-ext is ON, provide external model path using --ext-path"
    assert args.ext_id is not None, "When --init-ext is ON, provide external model type using --ext-id"
    assert "_" not in args.ext_id, "External model ID must not contain '_'"
    #assert args.ext_type == "pseudobulk", "In current implementation, --ext-type must be 'pseudobulk'"

    train_widths = [int(x) for x in args.train_width.split(",")] if args.train_width else [] #and not (args.use_external_model and args.external_model_type == "custom") else []

    # get the number of factors from the model path
    with gzip.open(args.ext_path, "rt") as f:
        header = f.readline()
    n_factor = len(header.strip().split("\t")) - 1
    #print(f"Number of factors in the external model: {n_factor}")

    train_params= [
        {
            "model_type": args.ext_type,
            "model_id": f"t{train_width}_x{args.ext_id}",
            "model_path": f"{args.out_dir}/ref.{args.ext_id}.model_matrix.tsv.gz",
            "train_width": train_width, 
            "n_factor": n_factor,
            **({"factor_map": f"{args.out_dir}/t{train_width}_x{args.ext_id}.factormap.tsv"} if not args.copy_ext_model else {})
        }
        for train_width in train_widths
    ]
    return train_params

def define_training_runs(args):
    if args.init_lda:
        train_params = define_lda_runs(args)
    elif args.init_ext:
        train_params = define_ext_runs(args)
    return train_params

def define_proj_runs(args):
    proj_runs = []
    train_params = define_training_runs(args)
    #print(f"Train params: {train_params}")
    for train_param in train_params:
        # lda/ext
        model_type = train_param["model_type"]
        train_width = train_param["train_width"]
        n_factor = train_param["n_factor"]
        model_id = train_param["model_id"]
        model_prefix = os.path.join(args.out_dir, model_id)
        model_path= f"{model_prefix}.model_matrix.tsv.gz"
        # proj
        fit_widths = [train_width] if args.fit_width is None else [int(x) for x in args.fit_width.split(",")]
        #print(f"fit_width: {fit_widths}")
        for fit_width in fit_widths:
            proj_id = f"{model_id}_p{fit_width}_a{args.anchor_res}"
            cmap_path = define_cmap(args, model_id)
            proj_runs.append({
                "model_type": model_type,
                "model_id": model_id,
                "model_path": model_path,
                "proj_id": proj_id,
                "n_factor": n_factor,
                "fit_width": fit_width,
                "cmap_path": cmap_path,
                "prerequisite_path": f"{model_prefix}.done"
            })
            #print(f"proj_id: {proj_id}")
            #print(f"Projection runs: {proj_runs}")
    #print(f"Projection runs: {proj_runs}")
    return proj_runs

def run_ficture(_args):
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

    # functions
    assert not (args.main and args.main_ext), "Cannot choose both --main and --main-ext"

    if args.main:
        args.sorttsv = True if args.in_transcript is not None else False
        args.minibatch = True
        args.segment=True
        args.init_lda = True
        args.projection = True
        args.decode = True
        args.summary = True
    elif args.main_ext:
        args.sorttsv = True if args.in_transcript is not None else False
        args.minibatch = True
        args.segment=True
        args.init_ext = True
        args.projection = True
        args.decode = True
        args.summary = True

    assert not (args.init_lda and args.init_ext), "Cannot choose both --init-lda and --init-ext"

    # input/output/other files
    # dirs
    os.makedirs(args.out_dir, exist_ok=True)

    # start mm
    mm = minimake()

    # in files
    if args.in_transcript is None:
        args.in_transcript = os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")
    
    if args.in_cstranscript is None:
        args.in_cstranscript = os.path.join(args.out_dir, "transcripts.sorted.tsv.gz")
    
    if args.in_minmax is None:
        args.in_minmax = os.path.join(args.out_dir, "coordinate_minmax.tsv")
    
    if args.in_feature is None:
        args.in_feature = os.path.join(args.out_dir, "feature.clean.tsv.gz")
    
    assert os.path.exists(args.in_transcript) or os.path.exists(args.in_cstranscript), "Provide at least one valid input transcript-indexed SGE file by --in-transcript or --in-cstranscript"
    assert os.path.exists(args.in_minmax), "Provide a valid input coordinate minmax file by --in-minmax"
    assert os.path.exists(args.in_feature), "Provide a valid input feature file by --in-feature"

    # feature customize when enabled 
    if any([args.include_feature_list, args.exclude_feature_list, args.include_feature_substr, args.exclude_feature_substr, args.include_feature_regex, args.exclude_feature_regex, args.include_feature_type_regex]):
        in_feature_ficture = os.path.join(args.out_dir, args.out_ficture_feature)
        in_feature_ficture_record = os.path.join(args.out_dir, args.out_ficture_feature.replace(".tsv.gz", ".record.tsv"))
        cmds = cmd_separator([], f"Customizing features for FICTURE analysis...")
        cmd = " ".join(["cartloader feature_filtering ",
                                    f"--in-csv {args.in_feature}", 
                                    f"--out-csv {in_feature_ficture}", 
                                    f"--out-record  {in_feature_ficture_record}",
                                    f"--include-feature-list {args.include_feature_list}" if args.include_feature_list is not None else "",
                                    f"--exclude-feature-list {args.exclude_feature_list}" if args.exclude_feature_list is not None else "",
                                    f"--include-feature-substr '{args.include_feature_substr}'" if args.include_feature_substr is not None else "",
                                    f"--exclude-feature-substr '{args.exclude_feature_substr}'" if args.exclude_feature_substr is not None else "",
                                    f"--include-feature-regex '{args.include_feature_regex}'" if args.include_feature_regex is not None else "",
                                    f"--exclude-feature-regex '{args.exclude_feature_regex}'" if args.exclude_feature_regex is not None else "",
                                    f"--include-feature-type-regex {args.include_feature_type_regex} --feature-type-ref {args.in_transcript} --feature-type-ref-colname-name gene --feature-type-ref-colname-type {args.colname_feature_type}" if args.include_feature_type_regex is not None and args.colname_feature_type is not None else "",
                                    f"--include-feature-type-regex {args.include_feature_type_regex} --feature-type-ref {args.feature_type_ref} --feature-type-ref-colidx-name {args.feature_type_ref_colidx_name}  --feature-type-ref-colidx-type {args.feature_type_ref_colidx_type}" if args.include_feature_type_regex is not None and args.feature_type_ref is not None else "",
                                    f"--log"
                                    ]) 
        cmds.append(cmd)
        mm.add_target(in_feature_ficture, [args.in_cstranscript, args.in_feature], cmds)
    else:
        in_feature_ficture = args.in_feature


    # out files
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.params.json") 
    major_axis = define_major_axis(args)

    # check static cmap file
    if args.cmap_static:
        if args.static_cmap_file is None:
            scriptdir = os.path.dirname(os.path.realpath(__file__))
            progdir = os.path.dirname(os.path.dirname(scriptdir))
            args.static_cmap_file = os.path.join(progdir, "assets", "fixed_color_map_60.tsv")
        assert os.path.exists(args.static_cmap_file), f"Static color map file {args.static_cmap_file} does not exist"

    # 1. sort 
    if args.sorttsv:
        scheck_app(args.gzip)
        scheck_app(args.sort)
        major_axis = define_major_axis(args)
        major_axis_1stcol = args.colidx_x if major_axis == 'X' else args.colidx_y
        major_axis_2ndcol = args.colidx_y if major_axis == 'X' else args.colidx_x
        sort_keys=f'-k{major_axis_1stcol},{major_axis_1stcol}g  -k{major_axis_2ndcol},{major_axis_2ndcol}g'

        cmds = cmd_separator([], f"Sorting the input transcript-indexed SGE file in tgz: {os.path.basename(args.in_transcript)}")
        cmds.append(f"{args.gzip} -dc {args.in_transcript} | {args.sort} -S {args.sort_mem} {sort_keys} | {args.gzip} -c > {args.in_cstranscript}")
        mm.add_target(args.in_cstranscript, [args.in_transcript], cmds)
    
    # 2. minibatch(create minibatch):
    if args.minibatch:
        scheck_app(args.gzip)
        scheck_app(args.sort)

        batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
        batch_mat_tsv = f"{args.out_dir}/batched.matrix.tsv"

        major_axis = define_major_axis(args)
        major_axis_col = args.colidx_x if major_axis == 'X' else args.colidx_y
        major_axis_col += 1 if major_axis_col != 1 else 0     # make_spatial_minibatch will insert a column on the 2nd position for a random_index.
        sort_keys= f"-k2,2n -k{major_axis_col},{major_axis_col}g"
        
        cmds = cmd_separator([], f"Creating minibatch from {os.path.basename(args.in_cstranscript)}...")
        cmds.append(f"ficture make_spatial_minibatch --input {args.in_cstranscript} --output {batch_mat_tsv} --mu_scale {args.mu_scale} --batch_size {args.minibatch_size} --batch_buff {args.minibatch_buffer} --major_axis {args.major_axis}")
        cmds.append(f"{args.sort} -S {args.sort_mem} {sort_keys} {batch_mat_tsv} | {args.gzip} -c > {batch_mat}")
        cmds.append(f"rm {batch_mat_tsv}")
        mm.add_target(batch_mat, [args.in_cstranscript], cmds)

    # 3. segment
    # segmentation - ficture
    if args.segment:
        scheck_app(args.gzip)
        scheck_app(args.sort)
        if args.hexagon_width is None:
            args.hexagon_width = args.train_width
        assert args.hexagon_width is not None, "When --segment, provide at least one hexagon width for segmentation in FICTURE-compatible format using --hexagon-width"
        hexagon_widths=[int(x) for x in args.hexagon_width.split(",")]

        for hexagon_width in hexagon_widths:
            hexagon_tsv=f"{args.out_dir}/hexagon.d_{hexagon_width}.tsv"
            hexagon=f"{args.out_dir}/hexagon.d_{hexagon_width}.tsv.gz"
            cmds = cmd_separator([], f"Creating hexagon-indexed SGE in FICTURE-compatible format for {hexagon_width}um...")
            cmd = " ".join([
                "ficture", "make_dge",
                f"--input {args.in_cstranscript}",
                f"--major_axis {major_axis}",
                f"--key {args.colname_count}",
                f"--output {hexagon_tsv}",
                f"--hex_width {hexagon_width}",
                f"--n_move {args.hexagon_n_move}",
                f"--mu_scale {args.mu_scale}",
                f"--precision {args.hexagon_precision}",
                f"--min_ct_per_unit {args.min_ct_per_unit_hexagon}"
                ])
            cmds.append(cmd)                
            cmds.append(f"{args.sort} -S {args.sort_mem} -k 1,1n {hexagon_tsv} | {args.gzip} -c > {hexagon}")
            cmds.append(f"rm {hexagon_tsv}")
            mm.add_target(f"{hexagon}", [args.in_cstranscript], cmds)
    
    # segmentation - 10x
    if args.segment_10x:       
        if args.hexagon_width_10x is None:
            args.hexagon_width_10x = args.train_width 
        assert args.hexagon_width_10x is not None, "When --segment-10x, Provide at least one hexagon width for segmentation in 10x Genomics format using --hexagon-width-10x"
        hexagon_widths_10x = [int(x) for x in args.hexagon_width_10x.split(",")]
        
        for hexagon_width in hexagon_widths_10x:
            hexagon_dir=f"{args.out_dir}/hexagon.d_{hexagon_width}.10x"
            cmds=cmd_separator([], f"Creating hexagon-indexed SGE in 10x Genomics format for {hexagon_width}um...")
            cmds.append(f"mkdir -p {hexagon_dir}")
            cmd = " ".join([
                "ficture", "make_sge_by_hexagon",
                f"--input {args.in_cstranscript}",
                f"--feature {in_feature_ficture}",
                f"--major_axis {major_axis}",
                f"--key {args.colname_count}",
                f"--output_path {hexagon_dir}",
                f"--hex_width {hexagon_width}",
                f"--n_move {args.hexagon_n_move_10x}",
                f"--mu_scale {args.mu_scale}",
                f"--precision {args.hexagon_precision_10x}",
                f"--min_ct_per_unit {args.min_ct_per_unit_hexagon_10x}",
                f"--transfer_gene_prefix"
                ])
            cmds.append(cmd)
            cmds.append(f"[ -f {hexagon_dir}/barcodes.tsv.gz ] && [ -f {hexagon_dir}/features.tsv.gz ] && [ -f {hexagon_dir}/matrix.mtx.gz ] && touch {args.out_dir}/hexagon.d_{hexagon_width}.10x.done")
            mm.add_target(f"{args.out_dir}/hexagon.d_{hexagon_width}.10x.done", [args.in_cstranscript, in_feature_ficture], cmds)

    # 4. lda
    if args.init_lda:
        lda_runs= define_lda_runs(args)
        for lda_params in lda_runs:
            # params & prefix
            train_width = lda_params["train_width"]
            n_factor = lda_params["n_factor"]

            model_id = lda_params["model_id"]
            model_prefix = os.path.join(args.out_dir, model_id)

            lda_fillr = int(train_width // 2 + 1)
            # files
            hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
            lda_model_matrix = f"{model_prefix}.model_matrix.tsv.gz"
            lda_fit_tsv = f"{model_prefix}.fit_result.tsv.gz"
            lda_postcount_tsv = f"{model_prefix}.posterior.count.tsv.gz"
            lda_de = f"{model_prefix}.bulk_chisq.tsv"
            # 1) fit model
            cmds = cmd_separator([], f"LDA training for {train_width}um and {n_factor} factors...")
            cmd = " ".join([
                "ficture", "fit_model",
                f"--input {hexagon}",
                f"--feature {in_feature_ficture}",
                f"--output {model_prefix}",
                f"--nFactor {n_factor}",
                f"--epoch {args.train_epoch}",
                f"--epoch_id_length {args.train_epoch_id_len}",
                f"--unit_attr X Y",
                f"--key {args.colname_count}",
                f"--min_ct_per_feature {args.min_ct_per_feature}", 
                f"--test_split 0.5",
                f"--R {args.lda_rand_init}",
                f"--thread {args.threads}",
                ])
            cmds.append(cmd)
            cmds.append(f"[ -f {lda_fit_tsv} ] && [ -f {lda_model_matrix} ] && [ -f {lda_postcount_tsv} ] && touch {model_prefix}.done" )
            mm.add_target(f"{model_prefix}.done", [args.in_cstranscript, hexagon, in_feature_ficture], cmds)

            # 2) choose color (generate a color map no matter skip coarse plot or not, given this color map may be needed for the projection or decode visual)
            cmap = args.static_cmap_file if args.cmap_static else f"{model_prefix}.rgb.tsv"
            cmds = run_choose_color(args, lda_fit_tsv, n_factor, model_prefix)
            if len(cmds) > 0:
                mm.add_target(cmap, [f"{model_prefix}.done"], cmds)

            # 3) coarse plot/DE/report
            if not args.skip_coarse_report:
                cmds = cmd_separator([], f" LDA visualization for {train_width}um and {n_factor} factors...")
                cmds.append(f"ficture plot_base --input {lda_fit_tsv} --output {model_prefix}.coarse --fill_range {lda_fillr} --color_table {cmap} --plot_um_per_pixel {args.lda_plot_um_per_pixel} --plot_discretized")
                mm.add_target(f"{model_prefix}.coarse.png", [f"{model_prefix}.done"], cmds)
            
            # 4) DE
            cmds = cmd_separator([], f" LDA DE/report for {train_width}um and {n_factor} factors...")
            cmds.append(f"ficture de_bulk --input {lda_postcount_tsv} --output {lda_de} --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold}")
            cmds.append(f"ficture factor_report --path {args.out_dir} --pref {model_id} --color_table {cmap}")
            cmds.append(f"[ -f {lda_de} ] && [ -f {model_prefix}.factor.info.html ] && touch {model_prefix}_summary.done")
            mm.add_target(f"{model_prefix}_summary.done", [f"{model_prefix}.done"], cmds)

    # 5. init_ext
    if args.init_ext:
        ext_runs= define_ext_runs(args)

        # get the number of factors from the model path
        with gzip.open(args.ext_path, "rt") as f:
            header = f.readline()
        n_factor = len(header.strip().split("\t")) - 1
        print(f"Number of factors in the external model: {n_factor}")

        ## copy the reference model matrix to the output directory for reproducbility
        print(f"Copying the external model to output directory")
        shutil.copyfile(args.ext_path, f"{args.out_dir}/ref.{args.ext_id}.model_matrix.tsv.gz")

        for ext_params in ext_runs:
            # params & prefix
            train_width = ext_params["train_width"]
            model_id = ext_params["model_id"]
            model_path = ext_params["model_path"]
            model_type = ext_params["model_type"]
            model_prefix = os.path.join(args.out_dir, model_id)

            assert model_type == "pseudobulk", "In current implementation, --ext-type must be 'pseudobulk'"

            ext_fillr = int(train_width // 2 + 1)
            # files
            hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
            init_postcount_tsv = f"{model_prefix}_ref.posterior.count.tsv.gz"
            init_fit_tsv = f"{model_prefix}_ref.fit_result.tsv.gz"

            ext_model_matrix = f"{model_prefix}.model_matrix.tsv.gz"
            ext_postcount_tsv = f"{model_prefix}.posterior.count.tsv.gz"
            ext_fit_tsv = f"{model_prefix}.fit_result.tsv.gz"
            ext_de = f"{model_prefix}.bulk_chisq.tsv"
            ext_factormap_tsv = f"{model_prefix}.factormap.tsv"

            # 1) fit model
            cmds = cmd_separator([], f"init_pseudo: Projecting {model_id} for {train_width}um...")
            if ( args.copy_ext_model ):
                #cmd = f"cp -f {model_path} {model_prefix}.model_matrix.tsv.gz"
                cmd = " ".join([
                "cartloader", "reheader_factors_model",
                f"--input-model {model_path}",
                f"--out {model_prefix}",
                f"--gzip '{args.gzip}' ",
                f"--log"
                ])
                cmds.append(cmd)
                cmd = " ".join([
                    "ficture", "transform",
                    f"--input {args.in_cstranscript}",
                    f"--feature {in_feature_ficture}",
                    f"--output_pref {model_prefix}",
                    f"--model {model_prefix}.model_matrix.tsv.gz",
                    f"--key {args.colname_count}",
                    f"--major_axis {major_axis}",
                    f"--hex_width {train_width}",
                    f"--n_move 1",
                    f"--min_ct_per_unit {args.min_ct_per_unit_fit}",
                    f"--mu_scale {args.mu_scale}",
                    f"--precision {args.fit_precision}",
                    f"--thread {args.threads}",
                    ])
                cmds.append(cmd)
                cmds.append(f"[ -f {ext_fit_tsv} ] && [ -f {ext_model_matrix} ] && [ -f {ext_postcount_tsv} ] && touch {model_prefix}.done" )
            else:
                cmd = " ".join([
                    "ficture", "init_model_from_pseudobulk",
                    f"--input {hexagon}",
                    f"--feature {in_feature_ficture}",
                    f"--output {model_prefix}_ref",
                    f"--model {model_path}",
                    f"--epoch 0",
                    f"--scale_model_rel -1",
                    f"--reorder_factors",
                    f"--key {args.colname_count}",
                    f"--min_ct_per_feature {args.min_ct_per_feature}", 
                    f"--thread 1",  ## use thread 1 because multithreading somehow does not work
                    ])
                cmds.append(cmd)
                
                cmd = " ".join([
                "cartloader", "reheader_factors_tsv",
                f"--fit-tsv {init_fit_tsv}",
                f"--postcount-tsv {init_postcount_tsv}",
                f"--out {model_prefix}",
                f"--gzip '{args.gzip}' ",
                f"--log"
                ])
                cmds.append(cmd)
                cmds.append(f"rm -f {init_fit_tsv} {init_postcount_tsv}")
                cmds.append(f"[ -f {ext_fit_tsv} ] && [ -f {ext_factormap_tsv} ] && [ -f {ext_model_matrix} ] && [ -f {ext_postcount_tsv} ] && touch {model_prefix}.done" )

            mm.add_target(f"{model_prefix}.done", [args.in_cstranscript, hexagon, in_feature_ficture], cmds)

            # 2) choose color 
            # if not args.skip_coarse_report:
            #     cmap = args.static_cmap_file if args.cmap_static else f"{model_prefix}.rgb.tsv"
            #     cmds = run_choose_color(args, ext_fit_tsv, n_factor, model_prefix)
            #     if len(cmds) > 0:
            #         mm.add_target(cmap, [f"{model_prefix}.done"], cmds)
            cmap = args.static_cmap_file if args.cmap_static else f"{model_prefix}.rgb.tsv"
            cmds = run_choose_color(args, ext_fit_tsv, n_factor, model_prefix)
            if len(cmds) > 0:
                mm.add_target(cmap, [f"{model_prefix}.done"], cmds)
            
            # 3) coarse plot
            if not args.skip_coarse_report:
                cmds = cmd_separator([], f"init_pseudo visualization and report for {train_width}um and {n_factor} factors...")
                cmds.append(f"ficture plot_base --input {ext_fit_tsv} --output {model_prefix}.coarse --fill_range {ext_fillr} --color_table {cmap} --plot_um_per_pixel {args.lda_plot_um_per_pixel} --plot_discretized")
                mm.add_target(f"{model_prefix}.coarse.png", [f"{model_prefix}.done", cmap], cmds)

            # 4) DE/report
            cmds = cmd_separator([], f"init_pseudo DE/report for {train_width}um and {n_factor} factors...")
            cmds.append(f"ficture de_bulk --input {ext_postcount_tsv} --output {ext_de} --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold}")
            cmds.append(f"ficture factor_report --path {args.out_dir} --pref {model_id} --color_table {cmap}")
            cmds.append(f"[ -f {ext_de} ] && [ -f {model_prefix}.factor.info.html ] && touch {model_prefix}_summary.done")
            mm.add_target(f"{model_prefix}_summary.done", [f"{model_prefix}.done", cmap], cmds)

    if args.projection:
        scheck_app(args.bgzip)
        scheck_app(args.tabix)
        scheck_app(args.sort)
        batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
        proj_runs = define_proj_runs(args)
        for proj_params in proj_runs:
            #print(proj_params)
            # input
            model_path = proj_params["model_path"]
            cmap_path = proj_params["cmap_path"]
            # prerequisities
            fit_prereq = proj_params["prerequisite_path"]
            # params & prefix
            fit_width = proj_params["fit_width"]
            proj_id = proj_params["proj_id"]
            proj_prefix = os.path.join(args.out_dir, proj_id)
            
            # if ( fit_width % args.anchor_res ) != 0:
            #     raise ValueError(f"fit_width {fit_width} must be divisible by anchor_res {args.anchor_res}")

            fit_n_move = int(fit_width / args.anchor_res)
            fit_fillr = int(args.anchor_res // 2 + 1)
            # files
            proj_fit_tsv = f"{proj_prefix}.fit_result.tsv.gz"
            proj_postcount = f"{proj_prefix}.posterior.count.tsv.gz"
            proj_de = f"{proj_prefix}.bulk_chisq.tsv"
            #1) transform/fit
            cmds=cmd_separator([], f"Creating projection, ID: {proj_id}")
            cmd = " ".join([
                "ficture", "transform",
                f"--input {args.in_cstranscript}",
                f"--feature {in_feature_ficture}",
                f"--output_pref {proj_prefix}",
                f"--model {model_path}",
                f"--key {args.colname_count}",
                f"--major_axis {major_axis}",
                f"--hex_width {fit_width}",
                f"--n_move {fit_n_move}",
                f"--min_ct_per_unit {args.min_ct_per_unit_fit}",
                f"--mu_scale {args.mu_scale}",
                f"--precision {args.fit_precision}",
                f"--thread {args.threads}",
                ])
            cmds.append(cmd)
            cmds.append(f"[ -f {proj_fit_tsv} ] && [ -f {proj_postcount} ] && touch {proj_prefix}.done" )
            mm.add_target(f"{proj_prefix}.done", [fit_prereq, in_feature_ficture], cmds)
            # 3) visualization
            cmds=cmd_separator([], f"Projection visualization, ID: {proj_id}")
            if not args.skip_coarse_report:
                cmd = " ".join([
                    "ficture", "plot_base",
                    f"--input {proj_fit_tsv}",
                    f"--output {proj_prefix}.coarse",
                    f"--fill_range {fit_fillr}",
                    f"--color_table {cmap_path}",
                    f"--plot_um_per_pixel {args.fit_plot_um_per_pixel}",
                    f"--plot_discretized"
                    ])
                cmds.append(cmd)
                mm.add_target(f"{proj_prefix}.coarse.png", [f"{proj_prefix}.done", cmap_path], cmds)
            # 4) DE/report
            cmds=cmd_separator([], f"Projection DE and report, ID: {proj_id}")
            # - transform-DE
            cmd = " ".join([
                "ficture", "de_bulk",
                f"--input {proj_postcount}",
                f"--output {proj_de}",
                f"--min_ct_per_feature {args.min_ct_per_feature}",
                f"--max_pval_output {args.de_max_pval}",
                f"--min_fold_output {args.de_min_fold}"
                ])
            cmds.append(cmd)
            # - transform-report
            cmds.append(f"ficture factor_report --path {args.out_dir} --pref {proj_id} --color_table {cmap_path}")
            # - done & target
            cmds.append(f"[ -f {proj_de} ] && [ -f {proj_prefix}.factor.info.html ] && touch {proj_prefix}_summary.done")
            mm.add_target(f"{proj_prefix}_summary.done", [f"{proj_prefix}.done", cmap_path], cmds)  

    if args.decode:
        scheck_app(args.bgzip)
        scheck_app(args.tabix)
        scheck_app(args.sort)
        major_axis = define_major_axis(args)
        # mini batch file
        batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
        # sort script
        script_path = f"{args.out_dir}/sort_decode.sh"
        
        with open(script_path, "w") as f:
            f.write(r"""#!/bin/bash
input=$1
output=$2
coor=$3
                    
n_factor=$4
bsize=$5
scale=$6
topk=$7
major_axis=$8
                    
bgzip=$9
tabix=${10}
sort=${11}
sort_mem=${12}

# 1) x y limits
while IFS=$'\t' read -r r_key r_val; do
    export "${r_key}"="${r_val}"
done < ${coor}
echo -e "x: ${xmin}, ${xmax}\ny: ${ymin}, ${ymax}"

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
                    

# 2) define the block and sort axis
axis2col_X=2
axis2col_Y=3

if [[ ${major_axis} == "Y" ]]; then
    block_axis="X"
    offblock="${offsetx}"
else
    block_axis="Y"
    offblock="${offsety}"
fi

blockidx0=$( echo "$(eval echo \${axis2col_${block_axis}}) - 1" | bc )   # perl is 0-based
sortidx=$(eval echo \${axis2col_${major_axis}})
                    
# echo
echo -e "block_axis: ${block_axis}\noffblock: ${offblock}\nblockidx in perl: ${blockidx0}\nsortidx: ${sortidx}"
        
header="##K=${n_factor};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=${block_axis};INDEX_AXIS=${major_axis}\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[$bidx]-$offb)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" -bidx="${blockidx0}" -offb="${offblock}"|  ${sort} -S ${sort_mem} -k1,1g -k"${sortidx},${sortidx}g") | ${bgzip} -c > ${output}

${tabix} -f -s1 -b"${sortidx}" -e"${sortidx}" ${output}
                                    
""")
        proj_runs = define_proj_runs(args)
        #print(f"Projection runs: {proj_runs}")
        #print(proj_runs)
        for proj_params in proj_runs:
            model_path = proj_params["model_path"]
            cmap_path  = proj_params["cmap_path"] # None or cmap path
            # params
            radius = int(args.anchor_res + args.radius_buffer)
            proj_id = proj_params["proj_id"]
            proj_prefix = os.path.join(args.out_dir, proj_id)
            decode_id=f"{proj_id}_r{radius}"
            decode_prefix = os.path.join(args.out_dir, decode_id)
            # files
            proj_fit_tsv=f"{proj_prefix}.fit_result.tsv.gz"
            decode_spixel=f"{decode_prefix}.pixel.sorted.tsv.gz"
            decode_postcount=f"{decode_prefix}.posterior.count.tsv.gz"
            decode_anchor=f"{decode_prefix}.anchor.tsv.gz"
            decode_de=f"{decode_prefix}.bulk_chisq.tsv"
            # 1) decode
            cmds=cmd_separator([], f"Performing pixel-level decoding, ID: {decode_id}")
            cmd = " ".join([
                "ficture", "slda_decode",
                f"--input {batch_mat}",
                f"--feature {in_feature_ficture}",
                f"--output {decode_prefix}",
                f"--model {model_path}",
                f"--anchor {proj_fit_tsv}",
                f"--anchor_in_um",
                f"--neighbor_radius {radius}",
                f"--mu_scale {args.mu_scale}",
                f"--key {args.colname_count}",
                f"--precision {args.decode_precision}",
                f"--lite_topk_output_pixel {args.decode_top_k}",
                f"--lite_topk_output_anchor {args.decode_top_k}",
                f"--thread {args.threads}"
                ])
            cmds.append(cmd)
            # - decode-sort
            cmds.append(f"bash {script_path} {decode_prefix}.pixel.tsv.gz {decode_spixel} {args.in_minmax} {proj_params['n_factor']} {args.decode_block_size} {args.decode_scale} {args.decode_top_k} {major_axis} {args.bgzip} {args.tabix} {args.sort} {args.sort_mem}")
            cmds.append(f"[ -f {decode_spixel} ] && [ -f {decode_postcount} ] && [ -f {decode_anchor} ] && touch {decode_prefix}.done")
            #cmds.append(f"rm {decode_prefix}.pixel.tsv.gz")
            mm.add_target(f"{decode_prefix}.done", [batch_mat, in_feature_ficture, f"{proj_prefix}.done"], cmds)
            # 2) decode-visualization
            cmds = cmd_separator([], f"Pixel-level decoding visualization and report, ID: {decode_id}")
            # - decode-pixel 
            cmd = " ".join([
                "ficture", "plot_pixel_full",
                f"--input {decode_spixel}",
                f"--color_table {cmap_path}",
                f"--output {decode_prefix}.pixel.png",
                f"--plot_um_per_pixel {args.decode_plot_um_per_pixel}",
                f"--full"
                ])
            cmds.append(cmd)
            mm.add_target(f"{decode_prefix}.pixel.png", [f"{decode_prefix}.done", cmap_path], cmds)
            # 3) decode-DE/report
            cmds = cmd_separator([], f"Pixel-level decoding DE and report, ID: {decode_id}")
            # - decode-DE
            cmd = " ".join([
                "ficture", "de_bulk",
                f"--input {decode_postcount}",
                f"--output {decode_de}",
                f"--min_ct_per_feature {args.min_ct_per_feature}",
                f"--max_pval_output {args.de_max_pval}",
                f"--min_fold_output {args.de_min_fold}"
                ])
            cmds.append(cmd)
            cmds.append(f"ficture factor_report --path {args.out_dir} --pref {decode_id} --color_table {cmap_path}")
            # done & target
            cmds.append(f"[ -f {decode_prefix}.bulk_chisq.tsv ] && [ -f {decode_prefix}.factor.info.html ] && touch {decode_prefix}_summary.done")
            mm.add_target(f"{decode_prefix}_summary.done", [batch_mat, f"{decode_prefix}.done", cmap_path], cmds)

    if args.summary:
        prerequisities=[]
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
                summary_cmap = args.static_cmap_file if args.cmap_static else f"{model_prefix}.rgb.tsv"
                summary_aux_args_models.append(f"lda,{model_prefix}.model_matrix.tsv.gz,{train_param['model_id']},{train_width},{n_factor},{summary_cmap}")
            summary_aux_args.append(" ".join(summary_aux_args_models))
        elif args.init_ext:
            summary_aux_args_models = ["--ext-model"]
            train_params = define_ext_runs(args)
            for train_param in train_params:
                train_width = train_param["train_width"]
                factor_map = train_param.get("factor_map", None)
                model_type = train_param["model_type"]
                model_path = train_param["model_path"]
                model_prefix = os.path.join(args.out_dir, train_param["model_id"])
                # prerequisities
                prerequisities.append(f"{model_prefix}.done")
                # args
                summary_cmap = args.static_cmap_file if args.cmap_static else f"{model_prefix}.rgb.tsv"
                summary_aux_args_models.append(f"{model_type},{model_path},{train_param['model_id']},{train_width},{factor_map},{summary_cmap}")
            summary_aux_args.append(" ".join(summary_aux_args_models))
        # projection & decode
        if args.projection or args.decode:
            summary_aux_args_projection = ["--projection"]
            summary_aux_args_decode = ["--decode"]
            proj_runs = define_proj_runs(args)
            for proj_params in proj_runs:
                model_type = proj_params["model_type"]
                model_id = proj_params["model_id"]                
                fit_width = proj_params["fit_width"]
                proj_id = proj_params["proj_id"]
                radius = int(args.anchor_res + args.radius_buffer) if args.decode else None
                decode_id=f"{proj_id}_r{radius}" if args.decode else None
                # prerequisities
                if args.projection:
                    prerequisities.append(f"{args.out_dir}/{proj_id}.done")
                if args.decode:
                    prerequisities.append(f"{args.out_dir}/{decode_id}.done")
                # args
                if args.projection:
                    summary_aux_args_projection.append(f"{model_type},{model_id},{proj_id},{fit_width},{args.anchor_res}")
                if args.decode:
                    summary_aux_args_decode.append(f"{model_type},{model_id},{proj_id},{decode_id},{radius}")
            if args.projection and len(summary_aux_args_projection) > 1:
                summary_aux_args.append(" ".join(summary_aux_args_projection))
            if args.decode and len(summary_aux_args_decode) > 1:
                summary_aux_args.append(" ".join(summary_aux_args_decode))
        # summary
        cmds = cmd_separator([], f"Summarizing output into to the <out_json> files...")
        cmd = " ".join([
            "cartloader", "write_json_for_ficture",
                "--merge",
                f"--in-cstranscript {args.in_cstranscript}",
                f"--in-feature {args.in_feature}", # use the original feature file for SGE
                f"--in-feature-ficture {in_feature_ficture}" if in_feature_ficture != args.in_feature else "",
                f"--in-minmax {args.in_minmax}",
                f"--out-dir {args.out_dir}",
                f"--out-json {args.out_json}",
                " ".join(summary_aux_args)
            ])
        cmds.append(cmd)
        mm.add_target(args.out_json, prerequisities, cmds)

    if args.viz_per_factor:
        proj_runs = define_proj_runs(args)
        for proj_params in proj_runs:
            fit_width = proj_params["fit_width"]
            proj_id = proj_params["proj_id"]
            radius = int(args.anchor_res + args.radius_buffer)
            decode_id=f"{proj_id}_r{radius}"
            decode_prefix = os.path.join(args.out_dir, decode_id)
            # files
            decode_spixel=f"{decode_prefix}.pixel.sorted.tsv.gz"
            # cmds
            cmds=cmd_separator([], f"Visualizing single factor at pixel level, ID: {decode_id}")
            cmds.append(f"ficture plot_pixel_single --input {decode_spixel}  --output {decode_prefix}.pixel --plot_um_per_pixel {args.decode_plot_um_per_pixel} --full --all")
            # done & target
            viz_list=[]
            for i in range(n_factor):
                viz_list.append(f"{decode_prefix}.pixel.F_{i}.png")
            check_vizperfactor_cmd = " && ".join([f"[ -f {file} ]" for file in viz_list])
            cmds.append(f"{check_vizperfactor_cmd} && touch {decode_prefix}.viz_per_factor.done")
            mm.add_target(f"{decode_prefix}.viz_per_factor.done", [f"{decode_prefix}.done"], cmds)

    if args.viz_dotplot:
        if args.lda:
            train_params = define_lda_runs(args)
            for train_param in train_params:
                train_width = train_param["train_width"]
                n_factor = train_param["n_factor"]
                model_id = train_param["model_id"]
                model_prefix = os.path.join(args.out_dir, train_param["model_id"])
                cmds = cmd_separator([], f" Dotplot visualization and report for LDA, ID: {model_id}.")
                cmds.append(f"factor_viz dotplot --post_count_file {model_prefix}.posterior.count.tsv.gz --meta_data_file {model_prefix}.factor.info.tsv --output {model_prefix}.dotplot.pdf")
                mm.add_target(f"{model_prefix}.dotplot.pdf", [f"{model_prefix}.done", f"{model_prefix}_summary.done"], cmds)
        if args.projection or args.decode:
            proj_runs = define_proj_runs(args)
            for proj_params in proj_runs:
                fit_width = proj_params["fit_width"]
                proj_id = proj_params["proj_id"]
                radius = int(args.anchor_res + args.radius_buffer)
                decode_id=f"{proj_id}_r{radius}"
                decode_prefix = os.path.join(args.out_dir, decode_id)
                # proj
                cmds = cmd_separator([], f" Dotplot visualization and report for projection, ID: {proj_id}")
                cmds.append(f"factor_viz dotplot --post_count_file {proj_prefix}.posterior.count.tsv.gz --meta_data_file {proj_prefix}.factor.info.tsv --output {proj_prefix}.dotplot.pdf")
                mm.add_target(f"{proj_prefix}.dotplot.pdf", [f"{proj_prefix}.done", f"{proj_prefix}_summary.done"], cmds)
                # decode
                if args.decode:
                    cmds = cmd_separator([], f" Dotplot visualization and report for decoding, ID: {decode_id}")
                    cmds.append(f"factor_viz dotplot --post_count_file {decode_prefix}.posterior.count.tsv.gz --meta_data_file {decode_prefix}.factor.info.tsv --output {decode_prefix}.dotplot.pdf")
                    mm.add_target(f"{decode_prefix}.dotplot.pdf", [f"{decode_prefix}.done", f"{decode_prefix}_summary.done"], cmds)

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
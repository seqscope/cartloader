import sys, os, gzip, argparse, logging, warnings, shutil

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
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run the main functions (sorttsv, minibatch, segment, lda, decode, summary). Note this does NOT include segment-10x and viz-per-factor')
    cmd_params.add_argument('--plus', action='store_true', default=False, help='Run the main functions (sorttsv, minibatch, segment, lda, decode, summary) and additional functions (segment-10x and viz-per-factor)')
    cmd_params.add_argument('--viz-only', action='store_true', default=False, help='For lda, decode functions, only regenerate the visualization and report without running LDA analysis, decode')
    cmd_params.add_argument('--sorttsv', action='store_true', default=False, help='(Main function) Sort the input tsv file')
    cmd_params.add_argument('--minibatch', action='store_true', default=False, help='(Main function) Perform minibatch step')
    cmd_params.add_argument('--segment', action='store_true', default=False, help='(Main function) Perform hexagon segmentation into FICTURE-compatible format')
    cmd_params.add_argument('--lda', action='store_true', default=False, help='(Main function) Perform LDA model training')
    cmd_params.add_argument('--projection', action='store_true', default=False, help='(Main function) Perform projection/Transform')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='(Main function) Perform pixel-level decoding')
    cmd_params.add_argument('--merge-by-pixel', action='store_true', default=False, help='(Main function) Merge pixel-level decoding results with the input TSV SGE into a single file')
    cmd_params.add_argument('--summary', action='store_true', default=False, help='(Main function) Generate a JSON file summarizing all fixture parameters for which outputs are available in the <out-dir>.')
    cmd_params.add_argument('--viz-per-factor', action='store_true', default=False, help='(Additional function) Generate pixel-level visualization for each factor')
    cmd_params.add_argument('--segment-10x', action='store_true', default=False, help='(Additional function) Perform hexagon segmentation into 10x Genomics format')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    key_params.add_argument('--out-json', type=str, default=None, help="Output JSON file for summarizing the ficture parameters. Default: <out-dir>/ficture.params.json ")
    key_params.add_argument('--in-transcript', type=str, default=None, help='Input unsorted transcript-indexed SGE file in TSV format. Default to transcripts.unsorted.tsv.gz in the output directory')
    key_params.add_argument('--in-cstranscript', type=str, default=None, help='(Optional) If a coordinate-sorted transcript-indexed SGE file (sorted by x and y coordinates) exists, specify it by --in-cstranscript to skip the sorting step. Default to transcripts.sorted.tsv.gz in the output directory')
    key_params.add_argument('--in-minmax', type=str, default=None, help='Input coordinate minmax TSV file. Default to coordinate_minmax.tsv in the output directory')
    key_params.add_argument('--in-feature', type=str, default=None,  help='(Optional) Input TSV file that specify which genes to use as input, e.g., feature.clean.tsv.gz. If absent, all genes will be used')
    key_params.add_argument('--major-axis', type=str, default=None, help='Axis where transcripts.tsv.gz are sorted. If not provided, it will be automatically defined by the longer axis. Options: X, Y')
    key_params.add_argument('--mu-scale', type=float, default=1.0, help='Scale factor for mu (pixels per um)')
    key_params.add_argument('--train-width', type=str, default="12", help='Hexagon flat-to-flat width (in um) during training. This width will be used to create hexagon-indexed SGE in FICTURE compatible format for LDA training. Use comma to specify multiple values.')
    key_params.add_argument('--hexagon-width-10x', type=str, default="12", help='Hexagon flat-to-flat width (in µm) used for creating hexagon-indexed SGE in 10x Genomics format. Separate multiple values with commas.')
    key_params.add_argument('--n-factor', type=str, default="12", help='Number of factors to train. Use comma to specify multiple values')
    key_params.add_argument('--anchor-res', type=float, default=4, help='Anchor resolution for decoding')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--csv-colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--csv-colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--key-col', type=str, default="Count", help='Columns from the input file to be used as key')
    # segmentation - ficture
    aux_params.add_argument('--hexagon-n-move', type=int, default=1, help='Level of hexagonal sliding when creating hexagon-indexed SGE in FICTURE compatible format')
    aux_params.add_argument('--hexagon-precision', type=float, default=2, help='Output precision of hexagon coordinates for FICTURE compatible format')
    aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=50, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format')
    # segmentation - 10x
    aux_params.add_argument('--hexagon-n-move-10x', type=int, default=1, help='Level of hexagonal sliding when creating hexagon-indexed SGE in 10x Genomics format')
    aux_params.add_argument('--hexagon-precision-10x', type=float, default=2, help='Output precision of hexagon coordinates for 10x Genomics format')
    aux_params.add_argument('--min-ct-per-unit-hexagon-10x', type=int, default=1, help='Minimum count per hexagon in hexagon segmentation in 10x Genomics format')
    # minibatch
    aux_params.add_argument('--minibatch-size', type=int, default=500, help='Batch size used in minibatch processing')
    aux_params.add_argument('--minibatch-buffer', type=int, default=30, help='Batch buffer used in minibatch processing')
    # train 
    aux_params.add_argument('--train-epoch', type=int, default=3, help='Training epoch for LDA model')
    aux_params.add_argument('--train-epoch-id-len', type=int, default=2, help='Training epoch ID length')
    aux_params.add_argument('--lda-rand-init', type=int, default=10, help='Number of random initialization during model training')
    aux_params.add_argument('--lda-plot-um-per-pixel', type=float, default=1, help='Image resolution for LDA plot')
    # fit 
    aux_params.add_argument('--fit-width',  type=str, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    aux_params.add_argument('--fit-precision', type=float, default=2, help='Output precision of model fitting')
    aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=20, help='Minimum count per hexagon unit during model fitting')
    aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot')   # in Scopeflow, this is set to 2
    # decode
    aux_params.add_argument('--decode-top-k', type=int, default=3, help='Top K columns to output in pixel-level decoding results')
    aux_params.add_argument('--decode-block-size', type=int, default=100, help='Block size for pixel decoding output')
    aux_params.add_argument('--decode-scale', type=int, default=100, help='Scale parameters for pixel decoding output')
    aux_params.add_argument('--decode-precision', type=float, default=0.01, help='Precision of pixel level decoding')
    aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=0.5, help='Image resolution for pixel decoding plot')
    # merge_by_pixel
    aux_params.add_argument('--merge-max-dist-um', type=float, default=0.1, help='Maximum distance in um for merging pixel-level decoding results') 
    aux_params.add_argument('--merge-max-k', type=int, default=1, help='Maximum number of K columns to output in merged pixel-level decoding results')
    aux_params.add_argument('--merge-max-p', type=int, default=1, help='Maximum number of P columns to output in merged pixel-level decoding results')
    # others parameters shared across steps
    aux_params.add_argument('--min-ct-per-feature', type=int, default=20, help='Minimum count per feature during LDA training, transform and decoding')
    aux_params.add_argument('--cmap-name', type=str, default="turbo", help='Name of color map')
    aux_params.add_argument('--extra-cmap', type=str, default=None, help='Extra color map for LDA visualization')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='p-value cutoff for differential expression')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression')
    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools.")
    env_params.add_argument('--spatula', type=str, default=None, help='Path to spatula binary. When not provided, it will use the spatula from the submodules.')    
    env_params.add_argument('--bgzip', type=str, default="bgzip", help='Path to bgzip binary. For faster processing, use "bgzip -@ 4')
    env_params.add_argument('--tabix', type=str, default="tabix", help='Path to tabix binary')
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default="5G", help='Memory size for each process')
    env_params.add_argument('--ficture', type=str, default="ficture", help='Path to ficture repository')
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def define_major_axis(args):
    if args.major_axis is None:
        args.major_axis = find_major_axis(args.in_minmax, "row")
    return args.major_axis

# ==update: replaced by using --csv-colidx-x and --csv-colidx-y==
# def get_col_idx(input_file, column, sep):
#     if input_file.endswith(".gz"):
#         with gzip.open(input_file, "rt") as f:
#             header = f.readline().strip()
#     else:
#         with open(input_file, "r") as f:
#             header = f.readline().strip()
#     headers = header.split(sep)
#     assert column in headers, f"Column {column} is not found in the input file"
#     col_idx = headers.index(column) + 1
#     return col_idx

# def define_sort_keys(input_file, columns, sep):
#     sort_keys = []
#     for col_flag in columns:
#         #print(col_flag)
#         column, flag = col_flag.split(",")
#         col_idx = get_col_idx(input_file, column, sep)
#         sort_key = f"-k{col_idx},{col_idx}{flag}"
#         sort_keys.append(sort_key)
#     return " ".join(sort_keys)
# ==END==

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

    if args.plus:
        args.main = True
        args.segment_10x = True
        args.viz_per_factor = True
        
    if args.main:
        if args.in_cstranscript is not None: ## coordinate-sorted transcript is already available
            args.sorttsv = False
        else:
            args.sorttsv = True
        args.minibatch = True
        args.segment=True
        args.lda = True
        args.projection = True
        args.decode = True
        args.summary = True
        #args.merge_by_pixel = True


    # parse parameters
    train_widths = [int(x) for x in args.train_width.split(",")]
    n_factors = [int(x) for x in args.n_factor.split(",")]
    hexagon_widths_10x = [int(x) for x in args.hexagon_width_10x.split(",")]
    
    # input/output
    # dirs
    os.makedirs(args.out_dir, exist_ok=True)
    # in files
    if args.in_transcript is None:
        args.in_transcript = os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")
    if args.in_cstranscript is None:
        args.in_cstranscript = os.path.join(args.out_dir, "transcripts.sorted.tsv.gz")
    if args.in_minmax is None:
        args.in_minmax = os.path.join(args.out_dir, "coordinate_minmax.tsv")
    # out files 
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.params.json") 
    # no default value for in_feature given that it is optional
    assert os.path.exists(args.in_transcript) or os.path.exists(args.in_cstranscript), "Provide a valid input transcript-indexed SGE file by --in-transcript or --in-cstranscript"
    assert os.path.exists(args.in_minmax), "Provide a valid input coordinate minmax file by --in-minmax"

    # start mm
    mm = minimake()

    # 1. sort 
    if args.sorttsv:
        scheck_app(args.gzip)
        scheck_app(args.sort)
        major_axis=define_major_axis(args)
        #sort_cols=[f"{major_axis},g"]+ [f"{axis},g" for axis in ["X", "Y"] if axis != major_axis]
        #sort_keys=define_sort_keys(args.in_transcript, sort_cols, "\t")
        if major_axis == 'X':
            sort_keys=f'-k{args.csv_colidx_x},{args.csv_colidx_x}g  -k{args.csv_colidx_y},{args.csv_colidx_y}g'
        elif major_axis == 'Y':
            sort_keys=f'-k{args.csv_colidx_y},{args.csv_colidx_y}g  -k{args.csv_colidx_x},{args.csv_colidx_x}g'
        cmds = cmd_separator([], f"Sorting the input transcript-indexed SGE file in tgz...")
        #cmds.append(f"skeys_intsv=$$(cartloader define_keys --sort --input {args.in_transcript} --columns {sort_cols})")
        #cmds.append(f"{args.gzip} -dc {args.in_transcript} | {args.sort} -S {args.sort_mem} $$skeys_intsv | {args.gzip} -c > {args.in_cstranscript}")
        cmds.append(f"{args.gzip} -dc {args.in_transcript} | {args.sort} -S {args.sort_mem} {sort_keys} | {args.gzip} -c > {args.in_cstranscript}")
        mm.add_target(args.in_cstranscript, [args.in_transcript], cmds)
    
    # 2. minibatch(create minibatch):
    if args.minibatch:
        scheck_app(args.gzip)
        scheck_app(args.sort)
        major_axis=define_major_axis(args)
        batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
        batch_mat_tsv = f"{args.out_dir}/batched.matrix.tsv"
        cmds = cmd_separator([], f"Creating minibatch from {args.in_cstranscript}...")
        cmds.append(f"ficture make_spatial_minibatch --input {args.in_cstranscript} --output {batch_mat_tsv} --mu_scale {args.mu_scale} --batch_size {args.minibatch_size} --batch_buff {args.minibatch_buffer} --major_axis {args.major_axis}")
        #major_axis_col = get_col_idx(args.in_cstranscript, major_axis, "\t")
        if major_axis == 'X':
            major_axis_col = args.csv_colidx_x
        elif major_axis == 'Y':
            major_axis_col = args.csv_colidx_y
        # make_spatial_minibatch will insert a column on the 2nd position for a random_index.
        if major_axis_col != 1:
            major_axis_col = major_axis_col + 1
        sort_keys= f"-k2,2n -k{major_axis_col},{major_axis_col}g"
        #cmds.append(f"skeys_minibatch=$$(cartloader define_keys --sort --input {batch_mat_tsv} --columns {sort_cols})")
        #cmds.append(f"{args.sort} -S {args.sort_mem} $$skeys_minibatch {batch_mat_tsv} | {args.gzip} -c > {batch_mat}")
        cmds.append(f"{args.sort} -S {args.sort_mem} {sort_keys} {batch_mat_tsv} | {args.gzip} -c > {batch_mat}")
        cmds.append(f"rm {batch_mat_tsv}")
        mm.add_target(batch_mat, [args.in_cstranscript], cmds)

    # 3. segment
    if args.segment:
        scheck_app(args.gzip)
        scheck_app(args.sort)
        major_axis=define_major_axis(args)
        for train_width in train_widths:
            hexagon_tsv=f"{args.out_dir}/hexagon.d_{train_width}.tsv"
            hexagon=f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
            cmds = cmd_separator([], f"Creating hexagon-indexed SGE for {train_width}um...")
            cmds.append(f"ficture make_dge --key {args.key_col} --input {args.in_cstranscript}  --major_axis {major_axis} --key {args.key_col} --output {hexagon_tsv} --hex_width {train_width} --n_move {args.hexagon_n_move} --mu_scale {args.mu_scale} --precision {args.hexagon_precision} --min_ct_per_unit {args.min_ct_per_unit_hexagon} ")
            cmds.append(f"{args.sort} -S {args.sort_mem} -k 1,1n {hexagon_tsv} | {args.gzip} -c > {hexagon}")
            cmds.append(f"rm {hexagon_tsv}")
            mm.add_target(f"{hexagon}", [args.in_cstranscript], cmds)

    if args.segment_10x:
        major_axis=define_major_axis(args)
        feature_arg = f"--feature {args.in_feature}" if args.in_feature is not None else ""
        for hexagon_width in hexagon_widths_10x:
            hexagon_dir=f"{args.out_dir}/hexagon.d_{hexagon_width}.10x"
            cmds=cmd_separator([], f"Creating hexagon-indexed SGE in 10x Genomics format for {hexagon_width}um...")
            cmds.append(f"mkdir -p {hexagon_dir}")
            cmds.append(f"ficture make_sge_by_hexagon --input {args.in_cstranscript} feature_arg --major_axis {major_axis} --key {args.key_col} --output_path {hexagon_dir} --hex_width {hexagon_width} --n_move {args.hexagon_n_move_10x} --mu_scale {args.mu_scale} --precision {args.hexagon_precision_10x} --min_ct_per_unit {args.min_ct_per_unit_hexagon_10x} --transfer_gene_prefix")
            # done & target
            cmds.append(f"[ -f {hexagon_dir}/barcodes.tsv.gz ] && [ -f {hexagon_dir}/features.tsv.gz ] && [ -f {hexagon_dir}/matrix.mtx.gz ] && touch {args.out_dir}/hexagon.d_{hexagon_width}.10x.done")
            mm.add_target(f"{args.out_dir}/hexagon.d_{hexagon_width}.10x.done", [args.in_cstranscript], cmds)

    # 4. lda
    if args.lda:
        for train_width in train_widths:
            # input
            hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
            feature_arg = f"--feature {args.in_feature}" if args.in_feature is not None else ""
            for n_factor in n_factors:
                # prefix
                lda_basename=f"t{train_width}_f{n_factor}"
                lda_prefix=os.path.join(args.out_dir, lda_basename)
                # 1) fit model
                lda_fit_tsv=f"{lda_prefix}.fit_result.tsv.gz"
                lda_model_tsv=f"{lda_prefix}.model_matrix.tsv.gz"
                lda_postcount_tsv=f"{lda_prefix}.posterior.count.tsv.gz"
                cmds = cmd_separator([], f" LDA training for {train_width}um and {n_factor} factors...")
                cmds.append(f"ficture fit_model --input {hexagon} --output {lda_prefix} {feature_arg} --nFactor {n_factor} --epoch {args.train_epoch} --epoch_id_length {args.train_epoch_id_len} --unit_attr X Y --key {args.key_col} --min_ct_per_feature {args.min_ct_per_feature} --test_split 0.5 --R {args.lda_rand_init} --thread {args.threads}")
                cmds.append(f"[ -f {lda_fit_tsv} ] && [ -f {lda_model_tsv} ] && [ -f {lda_postcount_tsv} ] && touch {lda_prefix}.done" )
                mm.add_target(f"{lda_prefix}.done", [args.in_cstranscript, hexagon], cmds)
                # 2) visualization and report
                lda_fillr = (train_width / 2 + 1)
                cmap=f"{lda_prefix}.rgb.tsv"
                cmds = cmd_separator([], f" LDA visualization and report for {train_width}um and {n_factor} factors...")
                # - choose_color
                if not os.path.exists(cmap):
                    cmds.append(f"ficture choose_color --input {lda_fit_tsv} --output {lda_prefix} --cmap_name {args.cmap_name}")
                # - coarse plot
                cmds.append(f"ficture plot_base --input {lda_fit_tsv} --output {lda_prefix}.coarse --fill_range {lda_fillr} --color_table {cmap} --plot_um_per_pixel {args.lda_plot_um_per_pixel} --plot_discretized")
                # - DE
                cmds.append(f"ficture de_bulk --input {lda_prefix}.posterior.count.tsv.gz --output {lda_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                # - report
                cmds.append(f"ficture factor_report --path {args.out_dir} --pref {lda_basename} --color_table {cmap}")
                # done & target
                cmds.append(f"[ -f {lda_prefix}.coarse.png ] && [ -f {lda_prefix}.bulk_chisq.tsv ] && [ -f {lda_prefix}.factor.info.html ] && touch {lda_prefix}_summary.done")
                mm.add_target(f"{lda_prefix}_summary.done", [f"{lda_prefix}.done"], cmds)
    
    if args.projection:
        scheck_app(args.bgzip)
        scheck_app(args.tabix)
        scheck_app(args.sort)
        major_axis=define_major_axis(args)
        for train_width in train_widths:
            for n_factor in n_factors:
                # prefix
                lda_basename=f"t{train_width}_f{n_factor}"
                lda_prefix=f"{args.out_dir}/{lda_basename}"
                # input
                #hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
                cmap=f"{lda_prefix}.rgb.tsv"
                model_mat=f"{lda_prefix}.model_matrix.tsv.gz"
                batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
                # fit_width
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [int(x) for x in args.fit_width.split(",")]
                for fit_width in fit_widths:
                    # params
                    fit_n_move = int(fit_width / args.anchor_res)
                    fit_fillr = int(args.anchor_res//2+1)
                    # prefix 
                    tsf_basename=f"{lda_basename}_p{fit_width}_a{args.anchor_res}"
                    tsf_prefix = os.path.join(args.out_dir, tsf_basename)
                    # files
                    tsf_fitres=f"{tsf_prefix}.fit_result.tsv.gz"
                    tsf_postcount=f"{tsf_prefix}.posterior.count.tsv.gz"
                    # 1) transform/fit
                    if not args.viz_only:
                        cmds=cmd_separator([], f"Creating projection for {train_width}um and {n_factor} factors, at {fit_width}um")
                        cmds.append(f"ficture transform --input {args.in_cstranscript} --output_pref {tsf_prefix} --model {model_mat} --key {args.key_col} --major_axis {major_axis} --hex_width {fit_width} --n_move {fit_n_move} --min_ct_per_unit {args.min_ct_per_unit_fit} --mu_scale {args.mu_scale} --thread {args.threads} --precision {args.fit_precision}")
                        cmds.append(f"[ -f {tsf_fitres} ] && [ -f {tsf_postcount} ] && touch {tsf_prefix}.done" )
                        mm.add_target(f"{tsf_prefix}.done", [f"{lda_prefix}.done"], cmds)
                    # 2) Transform visualization and report
                    cmds=cmd_separator([], f"Projection visualization and report for {train_width}um and {n_factor} factors, at {fit_width}um")
                    # - transform-cmap if cmap from lda does not exist
                    if not os.path.exists(cmap):
                        cmds.append(f"ficture choose_color --input {tsf_fitres} --output {tsf_prefix} --cmap_name {args.cmap_name}")
                        cmap=f"{tsf_prefix}.rgb.tsv"
                    # - transform-DE
                    cmds.append(f"ficture de_bulk --input {tsf_prefix}.posterior.count.tsv.gz --output {tsf_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                    # - transform-report
                    cmds.append(f"ficture factor_report --path {args.out_dir} --pref {tsf_basename} --color_table {cmap}")
                    # - transform-coarse-plot (add this step to be consistent with Scopeflow and NEDA)
                    cmds.append(f"ficture plot_base --input {tsf_prefix}.fit_result.tsv.gz --output {tsf_prefix}.coarse --fill_range {fit_fillr} --color_table {cmap} --plot_um_per_pixel {args.fit_plot_um_per_pixel} --plot_discretized")
                    # - done & target
                    cmds.append(f"[ -f {tsf_prefix}.bulk_chisq.tsv ] && [ -f {tsf_prefix}.factor.info.html ] && [ -f {tsf_prefix}.coarse.png ] && [ -f {tsf_prefix}.coarse.top.png ] && touch {tsf_prefix}_summary.done")
                    mm.add_target(f"{tsf_prefix}_summary.done", [f"{tsf_prefix}.done"], cmds)  

    if args.decode:
        scheck_app(args.bgzip)
        scheck_app(args.tabix)
        scheck_app(args.sort)
        major_axis=define_major_axis(args)
        script_path = f"{args.out_dir}/sort_decode.sh"
        
        if not args.viz_only:
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
declare -A axis2col
axis2col["X"]=2
axis2col["Y"]=3

if [[ ${major_axis} == "Y" ]]; then
    block_axis="X"
    offblock="${offsetx}"
else
    block_axis="Y"
    offblock="${offsety}"
fi

blockidx0=$( echo "${axis2col[${block_axis}]} - 1" | bc )   # perl is 0-based
sortidx=${axis2col[${major_axis}]}
                    
# echo
echo -e "block_axis: ${block_axis}\noffblock: ${offblock}\nblockidx in perl: ${blockidx0}\nsortidx: ${sortidx}"
        
header="##K=${n_factor};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=${block_axis};INDEX_AXIS=${major_axis}\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[$bidx]-$offb)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" -bidx="${blockidx0}" -offb="${offblock}"|  ${sort} -S ${sort_mem} -k1,1g -k"${sortidx},${sortidx}g") | ${bgzip} -c > ${output}

${tabix} -f -s1 -b"${sortidx}" -e"${sortidx}" ${output}
                                    
""")

        for train_width in train_widths:
            for n_factor in n_factors:
                # prefix
                lda_basename=f"t{train_width}_f{n_factor}"
                lda_prefix=f"{args.out_dir}/{lda_basename}"
                # input
                #hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
                cmap=f"{lda_prefix}.rgb.tsv"
                model_mat=f"{lda_prefix}.model_matrix.tsv.gz"
                batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
                # fit_width
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [int(x) for x in args.fit_width.split(",")]
                for fit_width in fit_widths:
                    # params
                    radius = args.anchor_res + 1
                    # prefix 
                    tsf_basename=f"{lda_basename}_p{fit_width}_a{args.anchor_res}"
                    tsf_prefix = os.path.join(args.out_dir, tsf_basename)
                    decode_basename=f"{tsf_basename}_r{radius}"
                    decode_prefix = os.path.join(args.out_dir, decode_basename)
                    # files
                    tsf_fitres=f"{tsf_prefix}.fit_result.tsv.gz"
                    decode_spixel=f"{decode_prefix}.pixel.sorted.tsv.gz"
                    decode_postcount=f"{decode_prefix}.posterior.count.tsv.gz"
                    decode_anchor=f"{decode_prefix}.anchor.tsv.gz"
                    # 3) decode
                    if not args.viz_only:
                        cmds=cmd_separator([], f"Performing pixel-level decoding for {train_width}um and {n_factor} factors, at {fit_width}um")
                        cmds.append(f"ficture slda_decode --input {batch_mat} --output {decode_prefix} --model {model_mat} --anchor {tsf_fitres} --anchor_in_um --neighbor_radius {radius} --mu_scale {args.mu_scale} --key {args.key_col} --precision {args.decode_precision} --lite_topk_output_pixel {args.decode_top_k} --lite_topk_output_anchor {args.decode_top_k} --thread {args.threads}")
                        # - decode-sort
                        cmds.append(f"bash {script_path} {decode_prefix}.pixel.tsv.gz {decode_spixel} {args.in_minmax} {n_factor} {args.decode_block_size} {args.decode_scale} {args.decode_top_k} {major_axis} {args.bgzip} {args.tabix} {args.sort} {args.sort_mem}")
                        cmds.append(f"[ -f {decode_spixel} ] && [ -f {decode_postcount} ] && [ -f {decode_anchor} ] && touch {decode_prefix}.done")
                        #cmds.append(f"rm {decode_prefix}.pixel.tsv.gz")
                        mm.add_target(f"{decode_prefix}.done", [batch_mat, f"{tsf_prefix}.done"], cmds)

                    # 4) decode-visualization and report
                    cmds=cmd_separator([], f"Pixel-level decoding visualization and report for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture de_bulk --input {decode_prefix}.posterior.count.tsv.gz --output {decode_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                    cmds.append(f"ficture factor_report --path {args.out_dir} --pref {decode_basename} --color_table {cmap}")
                    # - decode-pixel-plot
                    #cmds=cmd_separator(cmds, f"Drawing pixel-level output image for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture plot_pixel_full --input {decode_spixel} --color_table {cmap} --output {decode_prefix}.pixel.png --plot_um_per_pixel {args.decode_plot_um_per_pixel} --full")
                    # done & target
                    cmds.append(f"[ -f {decode_prefix}.bulk_chisq.tsv ] && [ -f {decode_prefix}.factor.info.html ] && [ -f {decode_prefix}.pixel.png ] && touch {decode_prefix}_summary.done")
                    mm.add_target(f"{decode_prefix}_summary.done", [batch_mat, f"{decode_prefix}.done"], cmds)
    
    if args.merge_by_pixel:
        # tools:
        # - spatula
        if args.spatula is None:
            args.spatula = os.path.join(cartloader_repo, "submodules", "spatula", "bin", "spatula")
            print(f"Given spatula is not provided, using the spatula from submodules: {args.spatula}.")
        scheck_app(args.spatula)
        # major axis
        major_axis=define_major_axis(args)
        # collect params
        prerequisities = []
        decode_basename_list=[]
        prerequisities.append(args.in_cstranscript)
        for train_width in train_widths:
            for n_factor in n_factors:
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [float(x) for x in args.fit_width.split(",")]
                for fit_width in fit_widths:
                    radius = args.anchor_res + 1
                    decode_basename=f"t{train_width}_f{n_factor}_p{fit_width}_a{args.anchor_res}_r{radius}"
                    decode_basename_list.append(decode_basename)
                    prerequisities.append(f"{args.out_dir}/{decode_basename}.done")     
        # cmds
        cmds = cmd_separator([], f"Merging pixel-level decoding results with the sorted SGE in TSV format...")
        pix_prefix_tsv_list=[f"--pix-prefix-tsv {decode_basename}__,{args.out_dir}/{decode_basename}.pixel.sorted.tsv.gz" for decode_basename in decode_basename_list]
        out_merge_prefix=f"{args.out_dir}/transcripts_pixel_join_K{args.merge_max_k}P{args.merge_max_p}.sorted"
        cmds.append(f"{args.spatula} join-pixel-tsv --mol-tsv {args.in_cstranscript} {' '.join(pix_prefix_tsv_list)} --out-prefix {out_merge_prefix} --max-dist-um {args.merge_max_dist_um} --out-max-k {args.merge_max_k} --out-max-p {args.merge_max_p}")
        mm.add_target(f"{out_merge_prefix}.tsv.gz", prerequisities, cmds)

    if args.summary:
        # collect prerequisities
        prerequisities = []
        for train_width in train_widths:
            for n_factor in n_factors:
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [int(x) for x in args.fit_width.split(",")]
                for fit_width in fit_widths:
                    # params
                    radius = args.anchor_res + 1
                    # basenames
                    lda_basename=f"t{train_width}_f{n_factor}"
                    tsf_basename=f"{lda_basename}_p{fit_width}_a{args.anchor_res}"
                    decode_basename=f"{tsf_basename}_r{radius}"
                    # prerequisities
                    if args.segment:
                        prerequisities.append(f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz")
                    if args.lda:
                        prerequisities.append(f"{args.out_dir}/{lda_basename}.done")
                    if args.decode:
                        prerequisities.append(f"{args.out_dir}/{tsf_basename}.done")
                        prerequisities.append(f"{args.out_dir}/{decode_basename}.done")
        # cmds
        cmds = cmd_separator([], f"Summarizing output into {args.out_json} files...")
        cmds.append(f"cartloader write_json_for_ficture --out-dir {args.out_dir} --out-json {args.out_json}")
        mm.add_target(args.out_json, prerequisities, cmds)

    if args.viz_per_factor:
        for train_width in train_widths:
            for n_factor in n_factors:
                lda_basename=f"t{train_width}_f{n_factor}"
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [int(x) for x in args.fit_width.split(",")]
                for fit_width in fit_widths:
                    # params
                    radius = args.anchor_res + 1
                    # prefix 
                    # basenames
                    tsf_basename=f"{lda_basename}_p{fit_width}_a{args.anchor_res}"
                    decode_basename=f"{tsf_basename}_r{radius}"
                    decode_prefix = os.path.join(args.out_dir, decode_basename)
                    # files
                    decode_spixel=f"{decode_prefix}.pixel.sorted.tsv.gz"
                    # 1) pixel-plot
                    cmds=cmd_separator([], f"Visualizing single factor at pixel level ({train_width}um and {n_factor} factors, at {fit_width}um)")
                    cmds.append(f"ficture plot_pixel_single --input {decode_spixel}  --output {decode_prefix}.pixel --plot_um_per_pixel {args.decode_plot_um_per_pixel} --full --all")
                    # done & target
                    viz_list=[]
                    for i in range(n_factor):
                        viz_list.append(f"{decode_prefix}.pixel.F_{i}.png")
                    check_vizperfactor_cmd = " && ".join([f"[ -f {file} ]" for file in viz_list])
                    cmds.append(f"{check_vizperfactor_cmd} && touch {decode_prefix}.viz_per_factor.done")
                    mm.add_target(f"{decode_prefix}.viz_per_factor.done", [f"{decode_prefix}.done"], cmds)

    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that ast least one run option was turned on")
        sys.exit(1)

    ## write makefile
    mm.write_makefile(f"{args.out_dir}/{args.makefn}")

    ## run makefile
    if args.dry_run:
        ## run makefile
        os.system(f"make -f {args.out_dir}/{args.makefn} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")
    else:
        os.system(f"make -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")

if __name__ == "__main__":
    # Get the path to the cartloader repository
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
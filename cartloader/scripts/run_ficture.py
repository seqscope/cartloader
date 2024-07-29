import sys, os, gzip, argparse, logging, warnings, shutil

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, find_major_axis, add_param_to_cmd

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader run_ficture", description="Run FICTURE")

    cmd_params = parser.add_argument_group("Commands", "FICTURE commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Run all FICTURE commands (minibatch, segment, lda, decode)')
    cmd_params.add_argument('--sorttsv', action='store_true', default=False, help='Sort the input tsv file')
    cmd_params.add_argument('--minibatch', action='store_true', default=False, help='Perform minibatch step')
    cmd_params.add_argument('--segment', action='store_true', default=False, help='Perform hexagon segmentation into FICTURE-compatible format')
    cmd_params.add_argument('--lda', action='store_true', default=False, help='Perform LDA model training')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='Perform pixel-level decoding')
    cmd_params.add_argument('--summary', action='store_true', default=False, help='Generate a summary yaml file')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default="Makefile", help='The name of the Makefile to generate')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    key_params.add_argument('--in-transcript', type=str, default=None, help='Input unsorted transcript-indexed SGE file in TSV format. Default to transcripts.unsorted.tsv.gz in the output directory')
    key_params.add_argument('--in-cstranscript', type=str, default=None, help='(Optional) If a coordinate-sorted transcript-indexed SGE file (sorted by x and y coordinates) exists, specify it by --in-cstranscript to skip the sorting step. Default to transcripts.sorted.tsv.gz in the output directory')
    key_params.add_argument('--in-minmax', type=str, default=None, help='Input coordinate minmax TSV file. Default to coordinate_minmax.tsv in the output directory')
    key_params.add_argument('--in-feature', type=str, default=None,  help='(Optional) Input TSV file that specify which genes to use as input, e.g., feature.clean.tsv.gz. If absent, all genes will be used')
    key_params.add_argument('--major-axis', type=str, default=None, help='Axis where transcripts.tsv.gz are sorted. If not provided, it will be automatically defined by the longer axis. Options: X, Y')
    key_params.add_argument('--mu-scale', type=float, default=1.0, help='Scale factor for mu (pixels per um)')
    key_params.add_argument('--train-width', type=str, default="12", help='Hexagon flat-to-flat width (in um) during training. Use comma to specify multiple values')
    key_params.add_argument('--n-factor', type=str, default="12", help='Number of factors to train. Use comma to specify multiple values')
    key_params.add_argument('--anchor-res', type=float, default=4, help='Anchor resolution for decoding')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # ficture params
    aux_params.add_argument('--csv-colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--csv-colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--train-epoch', type=int, default=3, help='Training epoch for LDA model')
    aux_params.add_argument('--train-epoch-id-len', type=int, default=2, help='Training epoch ID length')
    aux_params.add_argument('--train-n-move', type=int, default=1, help='Level of hexagonal sliding during training')
    aux_params.add_argument('--sge-n-move', type=int, default=1, help='Level of hexagonal sliding during SGE generation')
    aux_params.add_argument('--fit-width', type=float, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    aux_params.add_argument('--key-col', type=str, default="Count", help='Columns from the input file to be used as key')
    aux_params.add_argument('--minibatch-size', type=int, default=500, help='Batch size used in minibatch processing')
    aux_params.add_argument('--minibatch-buffer', type=int, default=30, help='Batch buffer used in minibatch processing')
    aux_params.add_argument('--min-ct-unit-dge', type=int, default=50, help='Minimum count per hexagon in DGE generation')
    aux_params.add_argument('--min-ct-unit-sge', type=int, default=1, help='Minimum count per hexagon in SGE generation')
    aux_params.add_argument('--min-ct-feature', type=int, default=20, help='Minimum count per feature during LDA training')
    aux_params.add_argument('--min-ct-unit-fit', type=int, default=20, help='Minimum count per hexagon unit during model fitting')
    aux_params.add_argument('--lda-rand-init', type=int, default=10, help='Number of random initialization during model training')
    aux_params.add_argument('--decode-top-k', type=int, default=3, help='Top K columns to output in pixel-level decoding results')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='p-value cutoff for differential expression')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression')
    aux_params.add_argument('--decode-block-size', type=int, default=100, help='Block size for pixel decoding output')
    aux_params.add_argument('--decode-scale', type=int, default=100, help='Scale parameters for pixel decoding output')
    aux_params.add_argument('--cmap-name', type=str, default="turbo", help='Name of color map')
    aux_params.add_argument('--dge-precision', type=float, default=2, help='Output precision of hexagon coordinates')
    aux_params.add_argument('--fit-precision', type=float, default=2, help='Output precision of model fitting')
    aux_params.add_argument('--decode-precision', type=float, default=0.1, help='Precision of pixel level decoding')
    aux_params.add_argument('--lda-plot-um-per-pixel', type=float, default=1, help='Image resolution for LDA plot')
    aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot')   # in Scopeflow, this is set to 2
    aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=0.5, help='Image resolution for pixel decoding plot')
    # summary yaml params
    aux_params.add_argument('--platform', type=str, default=None, help='Provide the information of platform to be added in the summary yaml file')
    aux_params.add_argument('--data-id', type=str, default=None, help='Provide a data id or description to be added in the summary yaml file')
    # applications
    aux_params.add_argument('--bgzip', type=str, default="bgzip", help='Path to bgzip binary. For faster processing, use "bgzip -@ 4')
    aux_params.add_argument('--tabix', type=str, default="tabix", help='Path to tabix binary')
    aux_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    aux_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    aux_params.add_argument('--ficture', type=str, default="ficture", help='Path to ficture repository')
    # mem
    aux_params.add_argument('--sort-mem', type=str, default="5G", help='Memory size for each process')

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

    if args.all:
        args.sorttsv = True
        args.minibatch = True
        args.segment=True
        args.lda = True
        args.decode = True
        args.summary = True

    # parse input parameters
    train_widths = [int(x) for x in args.train_width.split(",")]
    n_factors = [int(x) for x in args.n_factor.split(",")]
    
    # output
    os.makedirs(args.out_dir, exist_ok=True)
    if args.in_transcript is None:
        args.in_transcript = os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")
    if args.in_cstranscript is None:
        args.in_cstranscript = os.path.join(args.out_dir, "transcripts.sorted.tsv.gz")
    if args.in_minmax is None:
        args.in_minmax = os.path.join(args.out_dir, "coordinate_minmax.tsv")
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
            cmds = cmd_separator([], f"Creating DGE for {train_width}um...")
            cmds.append(f"ficture make_dge --key {args.key_col} --input {args.in_cstranscript} --output {hexagon_tsv} --hex_width {train_width} --n_move {args.train_n_move} --min_ct_per_unit {args.min_ct_unit_dge} --mu_scale {args.mu_scale} --precision {args.dge_precision} --major_axis {major_axis}")
            cmds.append(f"{args.sort} -S {args.sort_mem} -k 1,1n {hexagon_tsv} | {args.gzip} -c > {hexagon}")
            cmds.append(f"rm {hexagon_tsv}")
            mm.add_target(f"{hexagon}", [args.in_cstranscript], cmds)

    # 4. lda
    if args.lda:
        for train_width in train_widths:
            # input
            hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
            feature_arg = f"--feature {args.in_feature}" if args.in_feature is not None else ""
            for n_factor in n_factors:
                # prefix
                lda_prefix=f"{args.out_dir}/nF{n_factor}.d_{train_width}"
                # 1) fit model
                cmds = cmd_separator([], f"Creating LDA for {train_width}um and {n_factor} factors...")
                cmds.append(f"ficture fit_model --input {hexagon} --output {lda_prefix} {feature_arg} --nFactor {n_factor} --epoch {args.train_epoch} --epoch_id_length {args.train_epoch_id_len} --unit_attr X Y --key {args.key_col} --min_ct_per_feature {args.min_ct_feature} --test_split 0.5 --R {args.lda_rand_init} --thread {args.threads}")
                # 2) cmap
                lda_fit_tsv=f"{lda_prefix}.fit_result.tsv.gz"
                cmds.append(f"ficture choose_color --input {lda_fit_tsv} --output {lda_prefix} --cmap_name {args.cmap_name}")
                # 3) coarse plot
                lda_fillr = (train_width / 2 + 1)
                cmap=f"{lda_prefix}.rgb.tsv"
                cmds.append(f"ficture plot_base --input {lda_fit_tsv} --output {lda_prefix}.coarse --fill_range {lda_fillr} --color_table {cmap} --plot_um_per_pixel {args.lda_plot_um_per_pixel} --plot_discretized")
                # 4) DE
                cmds.append(f"ficture de_bulk --input {lda_prefix}.posterior.count.tsv.gz --output {lda_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                # 5) report
                lda_basename=os.path.basename(lda_prefix)
                cmds.append(f"ficture factor_report --path {args.out_dir} --pref {lda_basename} --color_table {cmap}")
                # done & target
                cmds.append(f"[ -f {lda_fit_tsv} ] && [ -f {cmap} ] && [ -f {lda_prefix}.coarse.png ] && [ -f {lda_prefix}.model_matrix.tsv.gz ] && [ -f {lda_prefix}.bulk_chisq.tsv ] && [ -f {lda_prefix}.factor.info.html ] && touch {lda_prefix}.done")
                mm.add_target(f"{lda_prefix}.done", [args.in_cstranscript, hexagon], cmds)

    if args.decode:
        scheck_app(args.bgzip)
        scheck_app(args.tabix)
        scheck_app(args.sort)
        major_axis=define_major_axis(args)
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

while IFS=$'\t' read -r r_key r_val; do
    export "${r_key}"="${r_val}"
done < ${coor}
echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
header="##K=${n_factor};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

if [ "${major_axis}" == "X" ]; then
    (echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" | ${sort} -S ${sort_mem} -k1,1g -k2,2g ) | ${bgzip} -c > ${output}
    ${tabix} -f -s1 -b2 -e2 ${output}
else
    (echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" | ${sort} -S ${sort_mem} -k1,1g -k3,3g ) | ${bgzip} -c > ${output}
    ${tabix} -f -s1 -b3 -e3 ${output}
fi
                    
#rm ${input}
""")

        for train_width in train_widths:
            for n_factor in n_factors:
                # prefix
                model_id=f"nF{n_factor}.d_{train_width}"
                lda_prefix=f"{args.out_dir}/{model_id}"
                # input
                #hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
                cmap=f"{lda_prefix}.rgb.tsv"
                model_mat=f"{lda_prefix}.model_matrix.tsv.gz"
                batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
                # fit_width
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [float(x) for x in args.fit_width.split(",")]

                for fit_width in fit_widths:
                    # params
                    fit_nmove = int(fit_width / args.anchor_res)
                    anchor_info=f"prj_{fit_width}.r_{args.anchor_res}"
                    fit_fillr = int(args.anchor_res//2+1)
                    radius = args.anchor_res + 1
                    # prefix 
                    tsf_prefix = f"{args.out_dir}/{model_id}.{anchor_info}"
                    decode_basename=f"{model_id}.decode.{anchor_info}_{radius}"
                    decode_prefix=f"{args.out_dir}/{decode_basename}"
                    # files
                    tsf_fitres=f"{tsf_prefix}.fit_result.tsv.gz"
                    decode_spixel=f"{decode_prefix}.pixel.sorted.tsv.gz"
                    # 1) transform/fit
                    cmds=cmd_separator([], f"Creating projection for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture transform --input {args.in_cstranscript} --output_pref {tsf_prefix} --model {model_mat} --key {args.key_col} --major_axis {major_axis} --hex_width {fit_width} --n_move {fit_nmove} --min_ct_per_unit {args.min_ct_unit_fit} --mu_scale {args.mu_scale} --thread {args.threads} --precision {args.fit_precision}")
                    # - transform-DE
                    cmds.append(f"ficture de_bulk --input {tsf_prefix}.posterior.count.tsv.gz --output {tsf_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                    # - transform-report
                    tsf_basename=os.path.basename(tsf_prefix)
                    cmds.append(f"ficture factor_report --path {args.out_dir} --pref {tsf_basename} --color_table {cmap}")
                    # - transform-coarse-plot (add this step to be consistent with Scopeflow and NEDA)
                    cmds.append(f"ficture plot_base --input {tsf_prefix}.fit_result.tsv.gz --output {tsf_prefix}.coarse --fill_range {fit_fillr} --color_table {cmap} --plot_um_per_pixel {args.fit_plot_um_per_pixel} --plot_discretized")
                    # done & target
                    cmds.append(f"[ -f {tsf_fitres} ] && [ -f {tsf_prefix}.bulk_chisq.tsv ] && [ -f {tsf_prefix}.factor.info.html ] && [ -f {tsf_prefix}.coarse.png ] && [ -f {tsf_prefix}.coarse.top.png ] && touch {tsf_prefix}.done")
                    mm.add_target(f"{tsf_prefix}.done", [f"{lda_prefix}.done"], cmds)
                    #
                    # 2) decode
                    cmds=cmd_separator([], f"Performing pixel-level decoding for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture slda_decode --input {batch_mat} --output {decode_prefix} --model {model_mat} --anchor {tsf_fitres} --anchor_in_um --neighbor_radius {radius} --mu_scale {args.mu_scale} --key {args.key_col} --precision {args.decode_precision} --lite_topk_output_pixel {args.decode_top_k} --lite_topk_output_anchor {args.decode_top_k} --thread {args.threads}")
                    # - decode-sort
                    #cmds=cmd_separator(cmds, f"Creating pixel-level output image for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"bash {script_path} {decode_prefix}.pixel.tsv.gz {decode_spixel} {args.in_minmax} {n_factor} {args.decode_block_size} {args.decode_scale} {args.decode_top_k} {major_axis} {args.bgzip} {args.tabix} {args.sort} {args.sort_mem}")
                    cmds.append(f"rm {decode_prefix}.pixel.tsv.gz")
                    # - decode-de & report
                    #cmds=cmd_separator(cmds, f"Performing pseudo-bulk differential expression analysis for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture de_bulk --input {decode_prefix}.posterior.count.tsv.gz --output {decode_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                    cmds.append(f"ficture factor_report --path {args.out_dir} --pref {decode_basename} --color_table {cmap}")
                    # - decode-pixel-plot
                    #cmds=cmd_separator(cmds, f"Drawing pixel-level output image for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture plot_pixel_full --input {decode_spixel} --color_table {cmap} --output {decode_prefix}.pixel.png --plot_um_per_pixel {args.decode_plot_um_per_pixel} --full")
                    # done & target
                    cmds.append(f"[ -f {decode_spixel} ] && [ -f {decode_prefix}.bulk_chisq.tsv ] && [ -f {decode_prefix}.factor.info.html ] && [ -f {decode_prefix}.pixel.png ] && touch {decode_prefix}.done")
                    mm.add_target(f"{decode_prefix}.done", [batch_mat, f"{lda_prefix}.done"], cmds)

    if args.summary:
        # arg names (aux)
        ans_description=["platform", "data_id"]
        aux_argset = set(item for lst in [ans_description] for item in lst)
        for train_width in train_widths:
            for n_factor in n_factors:
                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [float(x) for x in args.fit_width.split(",")]
                for fit_width in fit_widths:
                    radius = args.anchor_res + 1
                    out_yaml = os.path.join(args.out_dir, f"ficture.nF{args.n_factor}.d_{args.train_width}.decode.prj_{fit_width}.r_{args.anchor_res}_{radius}.yaml") 
                    cmds = cmd_separator([], f"Summarizing output into {out_yaml} files...")
                    yaml_cmds=f"cartloader write_yaml_for_ficture --out-dir {args.out_dir} --out-yaml {out_yaml} --in-transcript {args.in_transcript} --in-cstranscript {args.in_cstranscript} --in-minmax {args.in_minmax} --in-feature {args.in_feature} --train-width {train_width} --n-factor {n_factor} --fit-width {fit_width} --anchor-res {args.anchor_res} "
                    yaml_cmds = add_param_to_cmd(yaml_cmds, args, aux_argset)
                    cmds.append(yaml_cmds)
                    prerequisities = []
                    if args.segment:
                        prerequisities.append(f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz")
                    if args.lda:
                        prerequisities.append(f"{args.out_dir}/nF{n_factor}.d_{train_width}.done")
                    if args.decode:
                        prerequisities.append(f"{args.out_dir}/{model_id}.decode.prj_{fit_width}.r_{args.anchor_res}_{radius}.done")
                    mm.add_target(out_yaml, prerequisities, cmds)

    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
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
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
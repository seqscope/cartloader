import sys, os, gzip, argparse, logging, warnings, shutil

from ficture.utils.minimake import minimake

def parse_arguments(_args):
    if len(_args) == 0:
        parser.print_help()
        return

    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()

    cmd_params = parser.add_argument_group("Commands", "FICTURE commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Run all FICTURE commands (preprocess, segment, lda, decode)')
    cmd_params.add_argument('--preprocess', action='store_true', default=False, help='Perform preprocess step')
    cmd_params.add_argument('--lda', action='store_true', default=False, help='Perform LDA model training')
    cmd_params.add_argument('--decode', action='store_true', default=False, help='Perform pixel-level decoding')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    key_params.add_argument('--in-transcript', type=str, default=None, help='Input unsorted transcript-indexed SGE file in TSV format (e.g. transcript.unsorted.tsv.gz).')
    key_params.add_argument('--in-cstranscript', type=str, default=None, help='(Optional) Input transcript-indexed SGE file that sorted by x and y coordinates in TSV format (e.g. transcript.sorted.tsv.gz)')
    key_params.add_argument('--in-minmax', type=str, default=None, help='Input coordinate minmax TSV file (e.g. coordinate_minmax.tsv). If absent, it will be generated')
    key_params.add_argument('--in-feature', type=str, default=None,  help='Input TSV file (e.g. feature.clean.tsv.gz) that specify which genes to use as input. If absent, it will be use all genes')
    key_params.add_argument('--major-axis', type=str, default='Y', help='Axis where transcripts.tsv.gz are sorted')
    key_params.add_argument('--mu-scale', type=float, default=1.0, help='Scale factor for mu (pixels per um)')
    key_params.add_argument('--train-width', type=str, default="12", help='Hexagon flat-to-flat width (in um) during training. Use comma to specify multiple values')
    key_params.add_argument('--n-factor', type=str, default="12", help='Number of factors to train. Use comma to specify multiple values')
    key_params.add_argument('--anchor-res', type=float, default=4, help='Anchor resolution for decoding')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # ficture params
    aux_params.add_argument('--train-epoch', type=int, default=3, help='Training epoch for LDA model')
    aux_params.add_argument('--train-epoch-id-len', type=int, default=2, help='Training epoch ID length')
    aux_params.add_argument('--train-n-move', type=int, default=2, help='Level of hexagonal sliding during training')
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
    aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=0.5, help='Image resolution for pixel decoding plot')
    # applications
    aux_params.add_argument('--bgzip', type=str, default="bgzip", help='Path to bgzip binary. For faster processing, use "bgzip -@ 4')
    aux_params.add_argument('--tabix', type=str, default="tabix", help='Path to tabix binary')
    aux_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    aux_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    aux_params.add_argument('--ficture', type=str, default="ficture", help='Path to ficture repository')
    # mem
    aux_params.add_argument('--sort-mem', type=str, default="5G", help='Memory size for each process')
    return parser.parse_args()

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
        args.convert = True
        args.preprocess = True
        args.lda = True
        args.decode = True

    # start mm
    mm = minimake()

    # output
    os.makedirs(args.out_dir, exist_ok=True)
    
    # parse input parameters
    train_widths = [int(x) for x in args.train_width.split(",")]
    n_factors = [int(x) for x in args.n_factor.split(",")]
    
    # 1. sort 
    if args.sort:
        scheck_app(args.gzip)
        scheck_app(args.sort)
        cmds = cmd_separator([], f"Sorting the input transcript-indexed SGE file in tgz...")
        cmds.append(f"{args.gzip} -dc {args.in_transcript} | {args.sort} -S {args.sort_mem} -k1,1g -k2,2n | {args.gzip} -c > {args.in_cstranscript}")
        mm.add_target(args.in_cstranscript, [args.in_transcript], cmds)

    # 2. preprocess:
    if args.preprocess:
        batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
        batch_mat_tsv = f"{args.out_dir}/batched.matrix.tsv"
        cmds = cmd_separator([], f"$(info Creating minibatch from {args.in_transcript}...")
        ## create minibatch
        cmds.append(f"ficture make_spatial_minibatch --input {args.in_transcript} --output {batch_mat_tsv} --mu_scale {args.mu_scale} --batch_size {args.minibatch_size} --batch_buff {args.minibatch_buffer} --major_axis {args.major_axis}")
        cmds.append(f"{args.sort} -S {args.sort_mem} -k2,2n -k1,1g {batch_mat_tsv} | {args.gzip} -c > {batch_mat}")
        cmds.append(f"rm {batch_mat_tsv}")
        mm.add_target(batch_mat, [args.in_transcript], cmds)

    # 3. segment
    if args.segment:
        for train_width in train_widths:
            hexagon_tsv=f"{args.out_dir}/hexagon.d_{train_width}.tsv"
            hexagon=f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
            cmds = cmd_separator([], f"Creating DGE for {train_width}um...)")
            cmds.append(f"ficture make_dge --key {args.key_col} --input {args.in_transcript} --output {hexagon_tsv} --hex_width {train_width} --n_move {args.train_n_move} --min_ct_per_unit {args.min_ct_unit_dge} --mu_scale {args.mu_scale} --precision {args.dge_precision} --major_axis {args.major_axis}")
            cmds.append(f"{args.sort} -k 1,1n {hexagon_tsv} | {args.gzip} -c > {hexagon}")
            cmds.append(f"rm {hexagon_tsv}")
            mm.add_target(f"{hexagon}", [args.in_transcript], cmds)

    if args.lda:
        for train_width in train_widths:
            for n_factor in n_factors:
                hexagon = f"{args.out_dir}/hexagon.d_{train_width}.tsv.gz"
                model_id=f"nF{n_factor}.d_{train_width}"
                lda_prefix=f"{args.out_dir}/{model_id}"
                feature_arg = f"--feature {args.in_feature}" if args.in_feature is not None else ""
                cmds = cmd_separator([], f"Creating LDA for {train_width}um and {n_factor} factors...")
                cmds.append(f"ficture fit_model --input {hexagon} --output {lda_prefix} {feature_arg} --nFactor {n_factor} --epoch {args.train_epoch} --epoch_id_length {args.train_epoch_id_len} --unit_attr X Y --key {args.key_col} --min_ct_per_feature {args.min_ct_feature} --test_split 0.5 --R {args.lda_rand_init} --thread {args.threads}")

                fit_tsv=f"{lda_prefix}.fit_result.tsv.gz"
                cmds.append(f"ficture choose_color --input {fit_tsv} --output {lda_prefix} --cmap_name {args.cmap_name}")
                
                fillr = (train_width / 2 + 1)
                cmap=f"{lda_prefix}.rgb.tsv"
                cmds.append(f"ficture plot_base --input {fit_tsv} --output {lda_prefix}.coarse --fill_range {fillr} --color_table {cmap} --plot_um_per_pixel {args.lda_plot_um_per_pixel} --plot_discretized")
                cmds.append(f"touch {lda_prefix}.done")

                mm.add_target(f"{lda_prefix}.done", [args.in_transcript, hexagon], cmds)

    if args.decode:
        scheck_app(args.bgzip)
        scheck_app(args.tabix)

        script_path = f"{args.out_dir}/sort_decode.sh"
        with open(script_path, "w") as f:
            f.write(r"""#!/bin/bash
input=$1
output=$2
coor=$3
model_id=$4
bsize=$5
scale=$6
topk=$7
bgzip=$8
tabix=$9

K=$( echo $model_id | sed 's/nF\([0-9]\{1,\}\)\..*/\1/' )
while IFS=$'\t' read -r r_key r_val; do
    export "${r_key}"="${r_val}"
done < ${coor}
echo -e "${xmin}, ${xmax}; ${ymin}, ${ymax}"

offsetx=${xmin}
offsety=${ymin}
rangex=$( echo "(${xmax} - ${xmin} + 0.5)/1+1" | bc )
rangey=$( echo "(${ymax} - ${ymin} + 0.5)/1+1" | bc )
bsize=2000
scale=100
header="##K=${K};TOPK=${topk}\n##BLOCK_SIZE=${bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n##OFFSET_X=${offsetx};OFFSET_Y=${offsety};SIZE_X=${rangex};SIZE_Y=${rangey};SCALE=${scale}\n#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3"

(echo -e "${header}" && gzip -cd "${input}" | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\t", @F);' -- -bsize="${bsize}" -scale="${scale}" -offx="${offsetx}" -offy="${offsety}" | sort -S 1G -k1,1g -k3,3g ) | ${bgzip} -c > ${output}

${tabix} -f -s1 -b3 -e3 ${output}
rm ${input}
""")
        
        for train_width in train_widths:
            for n_factor in n_factors:
                batch_mat = f"{args.out_dir}/batched.matrix.tsv.gz"
                model_id=f"nF{n_factor}.d_{train_width}"
                lda_prefix=f"{args.out_dir}/{model_id}"
                cmap=f"{lda_prefix}.rgb.tsv"
                model=f"{lda_prefix}.model.p"

                if args.fit_width is None:
                    fit_widths = [train_width]
                else:
                    fit_widths = [float(x) for x in args.fit_width.split(",")]
                
                for fit_width in fit_widths:
                    fit_nmove = int(fit_width / args.anchor_res)
                    anchor_info=f"prj_{fit_width}.r_{args.anchor_res}"
                    radius = args.anchor_res + 1

                    prj_prefix = f"{args.out_dir}/{model_id}.{anchor_info}"
                    cmds=cmd_separator([], f"Creating projection for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture transform --input {args.in_transcript} --output_pref {prj_prefix} --model {model} --key {args.key_col} --major_axis {args.major_axis} --hex_width {fit_width} --n_move {fit_nmove} --min_ct_per_unit {args.min_ct_unit_fit} --mu_scale {args.mu_scale} --thread {args.threads} --precision {args.fit_precision}")

                    decode_basename=f"{model_id}.decode.{anchor_info}_{radius}"
                    decode_prefix=f"{args.out_dir}/{decode_basename}"

                    cmds=cmd_separator(cmds, f"Performing pixel-level decoding for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture slda_decode --input {batch_mat} --output {decode_prefix} --model {model} --anchor {prj_prefix}.fit_result.tsv.gz --anchor_in_um --neighbor_radius {radius} --mu_scale {args.mu_scale} --key {args.key_col} --precision {args.decode_precision} --lite_topk_output_pixel {args.decode_top_k} --lite_topk_output_anchor {args.decode_top_k} --thread {args.threads}")

                    cmds=cmd_separator(cmds, f"Creating pixel-level output image for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"bash {script_path} {decode_prefix}.pixel.tsv.gz {decode_prefix}.pixel.sorted.tsv.gz {args.in_minmax} {model_id} {args.decode_block_size} {args.decode_scale} {args.decode_top_k} {args.bgzip} {args.tabix}")

                    cmds=cmd_separator(cmds, f"Performing pseudo-bulk differential expression analysis for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture de_bulk --input {decode_prefix}.posterior.count.tsv.gz --output {decode_prefix}.bulk_chisq.tsv --min_ct_per_feature {args.min_ct_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold} --thread {args.threads}")
                    cmds.append(f"ficture factor_report --path {args.out_dir} --pref {decode_basename} --color_table {cmap}")

                    cmds=cmd_separator(cmds, f"Drawing pixel-level output image for {train_width}um and {n_factor} factors, at {fit_width}um")
                    cmds.append(f"ficture plot_pixel_full --input {decode_prefix}.pixel.sorted.tsv.gz --color_table {cmap} --output {decode_prefix}.pixel.png --plot_um_per_pixel {args.decode_plot_um_per_pixel} --full")

                    cmds.append(f"touch {decode_prefix}.done")
                    mm.add_target(f"{decode_prefix}.done", [batch_mat, hexagon,f"{lda_prefix}.done"], cmds)

    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)

    ## write makefile
    mm.write_makefile(f"{args.out_dir}/Makefile")

    ## run makefile
    if args.dry_run:
        ## run makefile
        os.system(f"make -f {args.out_dir}/Makefile -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {args.out_dir}/Makefile -j {args.n_jobs}")
    else:
        os.system(f"make -f {args.out_dir}/Makefile -j {args.n_jobs}")

if __name__ == "__main__":
    global cartloader
    cartloader=os.path.dirname(os.path.abspath(__file__))
    sys.path.extend(dir_i for dir_i in [cartloader] if dir_i not in sys.path)

    from utils import cmd_separator, scheck_app
    run_ficture(sys.argv[1:])
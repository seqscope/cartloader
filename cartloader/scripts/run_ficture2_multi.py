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
    run_params.add_argument('--makefn', type=str, default="run_ficture.mk", help='The name of the Makefile to generate (default: run_ficture.mk)')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output parameters for FICTURE")
    inout_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    inout_params.add_argument('--out-json', type=str, default=None, help="Output JSON file for summarizing the ficture parameters (default: <out-dir>/ficture.params.json)")
    inout_params.add_argument('--in-list', type=str, default=None, help='Path to the input list file containing the sample name and input transcript file')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in um) for LDA training (default: None)')
    key_params.add_argument('--n-factor', type=str, default=None, help='Comma-separated list of factor counts for LDA training.')
    key_params.add_argument('--anchor-res', type=int, default=6, help='Anchor resolution for decoding (default: 6)')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default="1G", help='Memory size for each process')
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",
    env_params.add_argument('--ficture2', type=str, required=True,  help='Path to punkst(ficture2) repository')
    env_params.add_argument('--python', type=str, default="python3",  help='Python3 binary')

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Customizing Auxiliary Parameters",
                                                      "Auxiliary parameters for customizing features in FICTURE analysis without modifying the original input feature TSV file. This ensures the original feature TSV file is retained in the output JSON file for downstream processing .")
    aux_ftrfilter_params.add_argument('--min-ct-per-ftr-tile', type=int, default=0, help='Apply a minimum count to filter overlapping feature. Filtering process will be applied if --min-ct-per-overlapftr > 0. (default: 0)')
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--colidx-feature',  type=int, default=3, help='Column index for Y-axis in the --in-transcript (default: 3)')
    aux_params.add_argument('--colidx-count',  type=int, default=4, help='Column index for intensity in the --in-transcript (default: 4)')
    aux_params.add_argument('--tile-size', type=int, default=500, help='Tile size for tiling (default: 500)')
    aux_params.add_argument('--tile-buffer', type=int, default=1000, help='Tile buffer for tiling (default: 1000)')
    aux_params.add_argument('--seed', type=int, default=1, help='Random seed for random number generation (default: 0)')

    # segmentation - ficture
    aux_params.add_argument('--min-count-per-sample', type=int, default=50, help='Minimum count per sample in the tiled SGE (default: 50)')
    aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=50, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format (default: 50)')
    # minibatch
    aux_params.add_argument('--minibatch-size', type=int, default=500, help='Batch size used in minibatch processing (default: 500)')
    # train
    aux_params.add_argument('--train-epoch', type=int, default=2, help='Training epoch for LDA model (default: 2)')
    aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=50, help='Minimum count per hexagon unit during model fitting (default: 20)')
    aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot (default: 1)')   # in Scopeflow, this is set to 2
    # color map
    aux_params.add_argument('--cmap-file', type=str, required=True, help='Define the path to the fixed color map (default: <cartloader_dir>/assets/fixed_color_map_60.tsv)')
    # others parameters shared across steps
    aux_params.add_argument('--min-count-train', type=int, default=50, help='Minimum count for training (default: 50)')
    aux_params.add_argument('--de-min-ct-per-feature', type=int, default=20, help='Minimum count per feature for differential expression test (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='p-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def define_lda_runs(args):
    #assert args.init_lda, "--init-lda must be ON when running define_lda_runs()"
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
        fit_widths = [train_width] 
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

def run_ficture2_multi(_args):
    """Run all functions in FICTURE2 with multi-sample pipeline
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
    
    # start mm
    mm = minimake()

    ficture2bin = os.path.join(args.ficture2, "bin/punkst")
    ficture2report = args.python + " " + os.path.join(args.ficture2, "ext/py/factor_report.py")
    
    # out files
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.params.json")

    if "," in args.width:
        raise ValueError("In multi-sample pipeline, when --train-width is provided, it should be a single value, not a comma-separated list. Use --width instead.")

    ## parse the input list file
    in_samples = []
    in_tsvs = []
    with flexopen(args.in_list, "rt") as f:
        for line in f:
            toks = line.strip().split("\t")
            in_samples.append(toks[0])
            in_tsvs.append(toks[1])

    # step 1. multi-sample tiling and hexagon:
    scheck_app(args.gzip)
    cmds = cmd_separator([], f"Creating tiled tsv from {os.path.basename(args.in_list)}...")
    cmd = " ".join([
        ficture2bin, "multisample-prepare",
        f"--in-tsv-list {args.in_list}",
        f"--out-dir {args.out_dir}",
        f"--out-joint-pref multi",
        f"--icol-x {args.colidx_x-1}",
        f"--icol-y {args.colidx_y-1}",
        f"--icol-feature {args.colidx_feature-1}",
        f"--icol-int {args.colidx_count-1}",
        f"--skip 1",
        f"--temp-dir {args.out_dir}/tmp",
        f"--tile-size {args.tile_size}",
        f"--tile-buffer {args.tile_buffer}",
        f"--threads {args.threads}",
        f"--hex-grid-dist {args.width}",
        f"--min-total-count-per-sample {args.min_count_per_sample}",
        f"--min-count {args.min_ct_per_unit_hexagon}",
        f"--include-feature-regex '{args.include_feature_regex}'" if args.include_feature_regex is not None else "",
        f"--exclude-feature-regex '{args.exclude_feature_regex}'" if args.exclude_feature_regex is not None else "",
    ])
    cmds.append(cmd)
    
    cmds.append(f"[ -f {args.out_dir}/multi.features.tsv ] && [ -f {args.out_dir}/multi.hex_{args.width}.txt ] && touch {args.out_dir}/multi.done" )
    mm.add_target(f"{args.out_dir}/multi.done", [args.in_list], cmds)

    # step 2. multi-sample LDA training
    scheck_app(args.spatula)
    lda_runs = define_lda_runs(args)
    for lda_params in lda_runs:
        # params & prefix
        train_width = lda_params["train_width"]
        n_factor = lda_params["n_factor"]

        model_id = lda_params["model_id"]
        model_prefix = os.path.join(args.out_dir, model_id)

        lda_fillr = int(train_width // 2 + 1)
        # files
        hexagon = os.path.join(args.out_dir, f"multi.hex_{train_width}.txt")
        meta = os.path.join(args.out_dir, f"multi.hex_{train_width}.json")
        lda_model_matrix = f"{model_prefix}.model.tsv"
        lda_fit_tsv = f"{model_prefix}.results.tsv"
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
            f"--minibatch-size {args.minibatch_size}",
            f"--seed {args.seed}",
            f"--n-epochs {args.train_epoch}",
            f"--threads {args.threads}",
            ])
        cmds.append(cmd)
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
        cmds.append(f"rm -f {model_prefix}.unsorted.model.tsv {model_prefix}.unsorted.results.tsv")
        cmds.append(f"[ -f {lda_fit_tsv}.gz ] && [ -f {lda_model_matrix} ] && touch {model_prefix}.done" )
        mm.add_target(f"{model_prefix}.done", [f"{args.out_dir}/multi.done"], cmds)

        # create color table
        out_cmap = f"{model_prefix}.cmap.tsv"
        with open(args.cmap_file, "r") as f:
            with open(out_cmap, "w") as f2:
                for i in range(n_factor+1):
                    f2.write(f.readline())

        # 2) DE
        cmds = cmd_separator([], f" LDA DE/report for {train_width}um and {n_factor} factors...")
        cmds.append(f"{args.spatula} diffexp-model-matrix --tsv1 {lda_model_matrix} --out {lda_de} --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}")
        cmds.append(f"({args.gzip} -cd {lda_de}.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {lda_de}.de.marginal.tsv.gz | tail -n +2 | sort -k 2,2n -k 3,3gr;) > {lda_de}")
        cmds.append(f"rm -f {lda_de}.de.marginal.tsv.gz")
#            cmds.append(f"{ficture2de} --input {lda_model_matrix} --output {lda_de} --feature_label Feature --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold}")
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
        mm.add_target(f"{model_prefix}_summary.done", [f"{model_prefix}.done"], cmds)

        # Perform LDA projection for each sample
        lda_each_targets = []
        for i in range(len(in_samples)):
            sample = in_samples[i]
            cmds = cmd_separator([], f"Creating projection for sample {sample}...")

            hexagon = os.path.join(args.out_dir, "samples", sample, f"{sample}.hex_{train_width}.txt")
            meta = os.path.join(args.out_dir, "samples", sample, f"{sample}.hex_{train_width}.json")
            lda_out_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{model_id}")
            lda_model_matrix = f"{model_prefix}.model.tsv"            
            lda_fit_tsv = f"{lda_out_prefix}.results.tsv"
            cmd = " ".join([
                ficture2bin, "lda4hex",
                f"--in-data {hexagon}",
                f"--in-meta {meta}",
                f"--projection-only", # skip LDA training, only project the data
                f"--model-prior {lda_model_matrix}",
                f"--out-prefix {lda_out_prefix}",
                f"--transform",
                f"--minibatch-size {args.minibatch_size}",
                f"--seed {args.seed}",
                f"--n-epochs {args.train_epoch}",
                f"--threads {args.threads}",
                ])
            cmds.append(cmd)
            cmds.append(f"{args.gzip} -f {lda_fit_tsv}")
            cmds.append(f"[ -f {lda_fit_tsv}.gz ] && touch {lda_out_prefix}.done")
            mm.add_target(f"{lda_out_prefix}.done", [f"{model_prefix}.done", f"{args.out_dir}/multi.done"], cmds)
            lda_each_targets.append(f"{lda_out_prefix}.done")
        cmds=cmd_separator([], f"Finishing LDA projection for each sample for model {model_id}...")
        cmds.append(f"touch {model_prefix}_each.done")
        mm.add_target(f"{model_prefix}_each.done", lda_each_targets, cmds)

    ## step 3. multi-sample pixel-decode
    scheck_app(args.sort)

    ## perform pixel-decode for each sample
    decode_runs = define_decode_runs(args)
    decode_each_targets = []
    for decode_params in decode_runs:            
        # input
        model_prefix = os.path.join(args.out_dir, decode_params["model_id"])
        model_path = decode_params["model_path"]
        cmap_path = decode_params["cmap_path"]
        # params & prefix
        fit_width = decode_params["fit_width"]
        decode_id = decode_params["decode_id"]
        fit_n_move = int(fit_width / args.anchor_res)

        for i in range(len(in_samples)):
            sample = in_samples[i]
            decode_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{decode_id}")
            decode_postcount = f"{decode_prefix}.pseudobulk.tsv"
            decode_fit_tsv = f"{decode_prefix}.tsv"
            decode_de = f"{decode_prefix}.bulk_chisq.tsv"

            #1) transform/fit
            cmds=cmd_separator([], f"Performing pixel-decode, ID {decode_id} for sample {sample}...")
            cmd = " ".join([
                ficture2bin, "pixel-decode",
                f"--model {model_path}",
                f"--in-tsv {args.out_dir}/samples/{sample}/{sample}.tiled.tsv",
                f"--in-index {args.out_dir}/samples/{sample}/{sample}.tiled.index",
                f"--temp-dir {args.out_dir}/tmp/{sample}_{decode_id}",
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
            cmds.append(f"[ -f {decode_fit_tsv} ] && [ -f {decode_postcount} ] && touch {decode_prefix}.tsv.done")
            mm.add_target(f"{decode_prefix}.tsv.done", [cmap_path, f"{args.out_dir}/multi.done", f"{model_prefix}.done"], cmds)  

            cmds=cmd_separator([], f"Performing post-decode tasks, ID {decode_id} for sample {sample}...")
            # - transform-DE
            cmds.append(f"{args.spatula} diffexp-model-matrix --tsv1 {decode_postcount} --out {decode_de} --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}")
            cmds.append(f"({args.gzip} -cd {decode_de}.de.marginal.tsv.gz | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd {decode_de}.de.marginal.tsv.gz | tail -n +2 | sort -k 2,2n -k 3,3gr;) > {decode_de}")
            cmds.append(f"rm -f {decode_de}.de.marginal.tsv.gz")

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

            # 4) visualization
            #cmds=cmd_separator([], f"Decode visualization, ID: {decode_id}")
            cmd = " ".join([
                f"{args.gzip} -dc {decode_fit_tsv}.gz |",
                ficture2bin, "draw-pixel-factors",
                f"--in-tsv /dev/stdin",
                f"--header-json {decode_prefix}.json",
                f"--in-color {cmap_path}",
                f"--out {decode_prefix}.png",
                f"--scale 1",
                f"--range {args.out_dir}/samples/{sample}/{sample}.tiled.coord_range.tsv"
                ])
            cmds.append(cmd)
            #mm.add_target(f"{decode_prefix}.png", [f"{decode_prefix}_summary.done", cmap_path], cmds)

            cmds.append(f"[ -f {decode_fit_tsv}.gz ] && [ -f {decode_postcount} ] && [ -f {decode_de} ] && [ -f {decode_prefix}.factor.info.html ] && [ -f {decode_fit_tsv}.gz ] && [ -f {decode_prefix}.png ] && touch {decode_prefix}.done")
            mm.add_target(f"{decode_prefix}.done", [cmap_path, f"{decode_prefix}.tsv.done", f"{args.out_dir}/multi.done", f"{model_prefix}.done"], cmds)
            decode_each_targets.append(f"{decode_prefix}.done")
        cmds=cmd_separator([], f"Finishing decode, ID: {decode_id}")
        cmds.append(f"touch {args.out_dir}/{decode_id}_each.done")
        mm.add_target(f"{args.out_dir}/{decode_id}_each.done", decode_each_targets, cmds)

    ## step 4. write the output JSON file for each sample
    json_each_targets = []
    for i in range(len(in_samples)):
        sample = in_samples[i]
        cmds = cmd_separator([], f"Writing output JSON file for sample {sample}...")
        sample_out_json = os.path.join(args.out_dir, "samples", sample, "ficture.params.json")
        sample_transcript = in_tsvs[i]
        sample_feature = os.path.join(args.out_dir, "samples", sample, f"{sample}.tiled.features.tsv")
        sample_minmax = os.path.join(args.out_dir, "samples", sample, f"{sample}.tiled.coord_range.tsv")
        sample_out_dir = os.path.join(args.out_dir, "samples", sample)
        
        summary_aux_args = []
        prerequisities = [f"{args.out_dir}/multi.done"]

        summary_aux_args_models = ["--lda-model"]
        lda_runs = define_lda_runs(args)
        for lda_params in lda_runs:
            # params & prefix
            train_width = lda_params["train_width"]
            n_factor = lda_params["n_factor"]
            model_prefix = os.path.join(args.out_dir, decode_params["model_id"])
            summary_cmap = f"{model_prefix}.cmap.tsv"
            prerequisities.append(f"{model_prefix}.done")
            summary_aux_args_models.append(f"lda,{model_prefix}.model.tsv,{lda_params['model_id']},{train_width},{n_factor},{summary_cmap}")
        summary_aux_args.append(" ".join(summary_aux_args_models))

        summary_aux_args_decodes = ["--decode"]
        decode_runs = define_decode_runs(args)
        for decode_params in decode_runs:
            model_type = decode_params["model_type"]
            model_id = decode_params["model_id"]
            fit_width = decode_params["fit_width"]
            decode_id = decode_params["decode_id"]
            decode_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{decode_id}")  ## decode_id contains the sample ID
            prerequisities.append(f"{decode_prefix}.done")
            summary_aux_args_decodes.append(f"{model_type},{model_id},{sample}.{decode_id},{fit_width},{args.anchor_res}")
        summary_aux_args.append(" ".join(summary_aux_args_decodes))

        cmd = " ".join([
            "cartloader", "write_json_for_ficture2",
                "--merge",
                f"--in-transcript {sample_transcript}",
                f"--in-feature {sample_feature}", # use the original feature file for SGE
                f"--in-feature-ficture {sample_feature}",
                f"--in-minmax {sample_minmax}",
                f"--out-dir {sample_out_dir}",
                f"--out-json {sample_out_json}",
                " ".join(summary_aux_args)
            ])
        cmds.append(cmd)
        mm.add_target(sample_out_json, prerequisities, cmds)
        json_each_targets.append(sample_out_json)

    cmds=cmd_separator([], f"Finishing writing the JSON file for each sample...")
    cmds.append(f"touch {args.out_dir}/multi_json_each.done")
    mm.add_target(f"{args.out_dir}/multi_json_each.done", json_each_targets, cmds)

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

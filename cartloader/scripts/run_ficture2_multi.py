import sys, os, gzip, argparse, logging, shutil, subprocess, inspect
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, flexopen, execute_makefile
from cartloader.utils.ficture2_helper import (
    define_lda_runs,
    define_decode_runs,
    add_umap_targets,
)

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

LDA_CONFIG = {
    "default_n_factor": "12",
    "allow_pretrained": True,
    "width_error_msg": "Provide at least one train width for LDA training using --width",
}

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Run FICTURE2")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile, and print commands without executing them')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and start from the beginning')
    run_params.add_argument('--threads', type=int, default=8, help='Maximum number of threads per job (default: 8)')
    run_params.add_argument('--n-jobs', type=int, default=2, help='Number of parallel jobs to run (default: 2)')
    run_params.add_argument('--makefn', type=str, default="run_ficture2_multi.mk", help='File name of Makefile to write (default: run_ficture2_multi.mk)')
    run_params.add_argument('--umap', action='store_true', default=False, help='Create UMAP plots for each LDA model')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output parameters for FICTURE")
    inout_params.add_argument('--out-dir', required=True, type=str, help='Output directory')
    inout_params.add_argument('--out-json', type=str, default=None, help="Path to output JSON file to store analysis parameters (default: <out-dir>/ficture.params.json)")
    inout_params.add_argument('--in-list', type=str, default=None, help='Path to input TSV with one row per sample: sample id and transcript file path')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--width', type=str, required=True, help='Comma-separated hexagon flat-to-flat widths (in um) for LDA training')
    key_params.add_argument('--n-factor', type=str, help='Comma-separated list of factor counts for LDA training.')
    key_params.add_argument('--anchor-res', type=int, default=6, help='Anchor resolution for decoding (default: 6)')
    key_params.add_argument('--cmap-file', type=str, default=os.path.join(repo_dir, "assets", "fixed_color_map_256.tsv"), help='Path to fixed color map TSV (default: <cartloader_dir>/assets/fixed_color_map_256.tsv)')

    # aux params
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    # input column indexes
    aux_params.add_argument('--colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    aux_params.add_argument('--colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    aux_params.add_argument('--colidx-feature',  type=int, default=3, help='Column index for Y-axis in the --in-transcript (default: 3)')
    aux_params.add_argument('--colidx-count',  type=int, default=4, help='Column index for intensity in the --in-transcript (default: 4)')
    # tile
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
    #aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=50, help='Minimum count per hexagon unit during model fitting (default: 20)')
    #aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=1, help='Image resolution for fit coarse plot (default: 1)')  # in Scopeflow, this is set to 2
    aux_params.add_argument('--skip-umap', action='store_true', default=False, help='Skip creating umap')
    aux_params.add_argument('--decode-scale', type=int, default=1, help='Decode scale parameter for plotting')

    # others parameters shared across steps
    aux_params.add_argument('--min-count-train', type=int, default=50, help='Minimum count for training (default: 50)')
    aux_params.add_argument('--de-min-ct-per-feature', type=int, default=20, help='Minimum count per feature for differential expression (default: 20)')
    aux_params.add_argument('--de-max-pval', type=float, default=1e-3, help='P-value cutoff for differential expression (default: 1e-3)')
    aux_params.add_argument('--de-min-fold', type=float, default=1.5, help='Fold-change cutoff for differential expression (default: 1.5)')
    aux_params.add_argument('--redo-pseudobulk-decode', action='store_true', default=False, help='Recompute pseudobulk decode with spatula. If set, the existing pseudobulk decode will be overwritten.')
    aux_params.add_argument('--redo-merge-units', action='store_true', default=False, help='Recompute merge units. If set, the existing mergeed hexagons and LDA results will be overwritten.')
    # project from external model
    aux_params.add_argument('--pretrained-model', type=str, help='Path to a pre-trained model to use for projection. If provided, LDA training will be skipped, and the provided model will be used for projection.')
    aux_params.add_argument('--retrain', action='store_true', default=False, help='If set, retain the pre-trained model. Only applicable when --pretrained-model is set.')

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
    env_params.add_argument('--R', type=str, default="R", help='Path to R binary for UMAP generation (default: R)')


    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def add_multisample_prepare_targets(mm, args, ficture2bin, in_samples):
    """Add Makefile target for multi-sample tiling and hexagon generation."""
    widths = args.width.split(",")

    cmds = cmd_separator([], f"Creating tiled tsv from {os.path.basename(args.in_list)}...")
    cmd = " ".join([
        f"'{ficture2bin}'", "multisample-prepare",
        f"--in-tsv-list '{args.in_list}'",
        f"--out-dir '{args.out_dir}'",
        f"--out-joint-pref multi",
        f"--icol-x {args.colidx_x-1}",
        f"--icol-y {args.colidx_y-1}",
        f"--icol-feature {args.colidx_feature-1}",
        f"--icol-int {args.colidx_count-1}",
        f"--skip 1",
        f"--temp-dir '{args.out_dir}/tmp'",
        f"--tile-size {args.tile_size}",
        f"--tile-buffer {args.tile_buffer}",
        f"--threads {args.threads}",
        f"--hex-grid-dist {args.width.replace(',', ' ')}",
        f"--min-total-count-per-sample {args.min_count_per_sample}",
        f"--min-count {args.min_ct_per_unit_hexagon}",
        f"--include-feature-regex '{args.include_feature_regex}'" if args.include_feature_regex is not None else "",
        f"--exclude-feature-regex '{args.exclude_feature_regex}'" if args.exclude_feature_regex is not None else "",
    ])
    cmds.append(cmd)

    if args.redo_merge_units:
        for width in widths:
            with flexopen(f"{args.out_dir}/multi.hex_{width}.list.tsv", "wt") as wf:
                for sample in in_samples:
                    wf.write(f"{sample}\t{args.out_dir}/samples/{sample}/{sample}.features.tsv\t{args.out_dir}/samples/{sample}/{sample}.hex_{width}.txt\t{args.out_dir}/samples/{sample}/{sample}.hex_{width}.json\t-2\n")
            cmd = " ".join([
                f"'{ficture2bin}'", "merge-units",
                f"--in-list '{args.out_dir}/multi.hex_{width}.list.tsv'",
                f"--min-total-count-per-sample {args.min_count_per_sample}",
                f"--min-count-per-unit {args.min_ct_per_unit_hexagon}",
                f"--out-pref '{args.out_dir}/multi.hex_{width}'",
                f"--threads {args.threads}",
                f"--temp-dir '{args.out_dir}/tmp/multi_hex_{width}'"])
            cmds.append(cmd)

    cmd = f"[ -f '{args.out_dir}/multi.features.tsv' ]" + "".join([f" && [ -f '{args.out_dir}/multi.hex_{width}.txt' ]" for width in widths]) + f" && touch '{args.out_dir}/multi.done'"
    cmds.append(cmd)
    mm.add_target(f"{args.out_dir}/multi.done", [args.in_list], cmds)

def add_lda_training_target(mm, args, ficture2bin, n_factor, train_width, model_prefix, hex_prefix, color_map, ficture2report):
    """Add Makefile target for training (or projecting) an LDA model."""
    cmds = cmd_separator([], f"LDA training for {train_width}um and {n_factor} factors...")
    
    unsorted_prefix = f"{model_prefix}.unsorted"
    lda_model_matrix = f"{model_prefix}.model.tsv"
    lda_fit_tsv = f"{model_prefix}.results.tsv.gz"
    lda_de = f"{model_prefix}.bulk_chisq.tsv"
    
    # copy model
    if args.pretrained_model is not None and not args.retrain:
        cmds.append(f"cp {args.pretrained_model} {unsorted_prefix}.model.tsv")
    
    # configure parameters
    if args.pretrained_model is not None and args.retrain:
        model_prior_arg = f"--model-prior {args.pretrained_model}"
        projection_only_arg = ""  
        n_topics_arg = ""
        sort_topics_arg = ""
    elif args.pretrained_model is not None and not args.retrain:
        model_prior_arg = f"--model-prior {unsorted_prefix}.model.tsv"
        projection_only_arg = "--projection-only"  
        n_topics_arg = ""
        sort_topics_arg = ""
    else:
        model_prior_arg = ""
        projection_only_arg = ""
        n_topics_arg = f"--n-topics {n_factor}"
        sort_topics_arg = "--sort-topics"
    
    # 1) train LDA
    train_cmd = " ".join([
        f"'{ficture2bin}'", "lda4hex",
        f"--in-data '{hex_prefix}.txt'",
        f"--in-meta '{hex_prefix}.json'",
        f"--out-prefix '{unsorted_prefix}'",
        model_prior_arg,
        projection_only_arg,
        n_topics_arg,
        sort_topics_arg,
        "--transform",
        f"--minibatch-size {args.minibatch_size}",
        f"--seed {args.seed}",
        f"--n-epochs {args.train_epoch}",
        f"--threads {args.threads}",
    ])
    cmds.append(train_cmd)

    # 2) append topk
    append_cmd = " ".join([
        f"'{args.spatula}'", "append-topk-tsv",
        f"--in-model '{unsorted_prefix}.model.tsv'",
        f"--in-json '{hex_prefix}.json'",
        f"--out-model '{lda_model_matrix}'",
#        "--reorder",
        f"--in-tsv '{unsorted_prefix}.results.tsv'",
        f"--out-tsv '{lda_fit_tsv}'",
        "--offset-model 1"
    ])
    cmds.append(append_cmd)
    cmds.append(f"cp '{unsorted_prefix}.model.tsv' '{lda_model_matrix}'")
    cmds.append(f"rm -f '{unsorted_prefix}.model.tsv' '{unsorted_prefix}.results.tsv'")
    cmds.append(f"[ -f '{lda_fit_tsv}' ] && [ -f '{lda_model_matrix}' ] && touch '{model_prefix}.done'")
    mm.add_target(f"{model_prefix}.done", [f"{args.out_dir}/multi.done"], cmds)

    # 3) create color table
    cmds = cmd_separator([], f"Generate the color map ")
    cmds.append(f'head -n $(({n_factor} + 1)) "{args.cmap_file}" > "{color_map}"')
    mm.add_target(color_map, [args.cmap_file], cmds)

    # 4) DE
    cmds = cmd_separator([], f" LDA DE/report for {train_width}um and {n_factor} factors...")
    cmds.append(f"'{args.spatula}' diffexp-model-matrix --tsv1 '{lda_model_matrix}' --out '{lda_de}' --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}")
    cmds.append(f"({args.gzip} -cd '{lda_de}.de.marginal.tsv.gz' | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd '{lda_de}.de.marginal.tsv.gz' | tail -n +2 | sort -k 2,2n -k 3,3gr;) > '{lda_de}'")
    cmds.append(f"rm -f '{lda_de}.de.marginal.tsv.gz'")
    #cmds.append(f"{ficture2de} --input {lda_model_matrix} --output {lda_de} --feature_label Feature --min_ct_per_feature {args.min_ct_per_feature} --max_pval_output {args.de_max_pval} --min_fold_output {args.de_min_fold}")
    cmd = " ".join([
        f"{ficture2report}",
        f"--factor_label factor",
        f"--de '{lda_de}'",
        f"--pseudobulk '{lda_model_matrix}'",
        f"--feature_label Feature",
        f"--color_table '{color_map}'",
        f"--output_pref '{model_prefix}.factor.info'"
        ])
    cmds.append(cmd)
    cmds.append(f"[ -f '{lda_de}' ] && [ -f '{model_prefix}.factor.info.html' ] && touch '{model_prefix}_summary.done'")
    mm.add_target(f"{model_prefix}_summary.done", [f"{model_prefix}.done", color_map], cmds)

def add_projection_target_per_sample(mm, args, ficture2bin, model_prefix, model_id, sample, train_width):
    """Add Makefile target that projects a trained LDA model onto a single sample."""
    cmds = cmd_separator([], f"Creating projection for sample {sample}...")

    lda_model_matrix = f"{model_prefix}.model.tsv"
    sample_hex_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.hex_{train_width}")
    sample_lda_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{model_id}")
    sample_lda_fit_tsv = f"{sample_lda_prefix}.results.tsv.gz"

    cmd = " ".join([
        f"{ficture2bin}", "lda4hex",
        f"--in-data '{sample_hex_prefix}.txt'",
        f"--in-meta '{sample_hex_prefix}.json'",
        "--projection-only",
        f"--model-prior '{lda_model_matrix}'",
        f"--out-prefix '{sample_lda_prefix}.unsorted'",
        "--transform",
        f"--minibatch-size {args.minibatch_size}",
        f"--seed {args.seed}",
        f"--n-epochs {args.train_epoch}",
        f"--threads {args.threads}",
    ])
    cmds.append(cmd)

    cmd = " ".join([
        f"'{args.spatula}'", "append-topk-tsv",
        f"--in-model '{lda_model_matrix}'",
        f"--out-model '{sample_lda_prefix}.model.tsv'",
        f"--in-tsv '{sample_lda_prefix}.unsorted.results.tsv'",
        f"--out-tsv '{sample_lda_fit_tsv}'",
        "--offset-model 1",
        "--offset-data 3",
        "--icol-random-key 0"
    ])
    cmds.append(cmd)
    cmds.append(f"rm -f '{sample_lda_prefix}.unsorted.results.tsv'")
    cmds.append(f"[ -f '{sample_lda_fit_tsv}' ] && touch '{sample_lda_prefix}.done'")
    mm.add_target(f"{sample_lda_prefix}.done", [f"{model_prefix}.done", f"{args.out_dir}/multi.done"], cmds)
    return f"{sample_lda_prefix}.done"

def add_pixel_decode_target_per_sample(mm, args, ficture2bin, ficture2report, model_prefix, model_path, cmap_path, decode_id, fit_width, n_factor, fit_n_move, sample):
    """Add Makefile targets to pixel-decode and post-process results for a single sample."""
    decode_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{decode_id}")
    decode_postcount = f"{decode_prefix}.pseudobulk.tsv"
    decode_fit_tsv = f"{decode_prefix}.tsv"
    decode_de = f"{decode_prefix}.bulk_chisq.tsv"

    cmds = cmd_separator([], f"Performing pixel-decode, ID {decode_id} for sample {sample}...")
    cmd = " ".join([
        f"'{ficture2bin}'", "pixel-decode",
        f"--model '{model_path}'",
        f"--in-tsv '{args.out_dir}/samples/{sample}/{sample}.tiled.tsv'",
        f"--in-index '{args.out_dir}/samples/{sample}/{sample}.tiled.index'",
        f"--temp-dir '{args.out_dir}/tmp/{sample}_{decode_id}'",
        f"--out-pref '{decode_prefix}'",
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

    if args.redo_pseudobulk_decode:
        cmds.append(f"rm -f {decode_postcount}")
        cmd = " ".join([
            f"'{args.spatula}'", "pseudobulk-from-decode",
            f"--tsv '{decode_fit_tsv}'",
            f"--out '{decode_postcount}'",
            f"--n-factors {n_factor}"
        ])
        cmds.append(cmd)

    cmds.append(f"{args.gzip} -f '{decode_fit_tsv}'")
    cmds.append(f"{args.gzip} -f '{decode_postcount}'")
    cmds.append(f"[ -f '{decode_fit_tsv}.gz' ] && [ -f '{decode_postcount}.gz' ] && touch '{decode_prefix}.tsv.done'")
    mm.add_target(f"{decode_prefix}.tsv.done", [cmap_path, f"{args.out_dir}/multi.done", f"{model_prefix}.done"], cmds)

    cmds = cmd_separator([], f"Performing post-decode tasks, ID {decode_id} for sample {sample}...")
    cmds.append(f"'{args.spatula}' diffexp-model-matrix --tsv1 '{decode_postcount}.gz' --out '{decode_de}' --min-count {args.de_min_ct_per_feature} --max-pval {args.de_max_pval} --min-fc {args.de_min_fold}")
    cmds.append(f"({args.gzip} -cd '{decode_de}.de.marginal.tsv.gz' | head -1 | sed 's/^Feature/gene/'; {args.gzip} -cd '{decode_de}.de.marginal.tsv.gz' | tail -n +2 | {args.sort} -k 2,2n -k 3,3gr;) > '{decode_de}'")
    cmds.append(f"rm -f '{decode_de}.de.marginal.tsv.gz'")

    cmd = " ".join([
        f"{ficture2report}",
        f"--factor_label factor",
        f"--de '{decode_de}'",
        f"--pseudobulk '{decode_postcount}.gz'",
        f"--feature_label Feature",
        f"--color_table '{cmap_path}'",
        f"--output_pref '{decode_prefix}.factor.info'"
    ])
    cmds.append(cmd)

    cmds.append(f"{args.gzip} -dc '{decode_fit_tsv}.gz' > '{decode_fit_tsv}'")
    cmd = " ".join([
        f"'{ficture2bin}'", "draw-pixel-factors",
        f"--in-tsv '{decode_fit_tsv}'",
        f"--header-json '{decode_prefix}.json'",
        f"--in-color '{cmap_path}'",
        f"--out '{decode_prefix}.png'",
        f"--scale {args.decode_scale}",
        f"--range '{args.out_dir}/samples/{sample}/{sample}.tiled.coord_range.tsv'"
    ])
    cmds.append(cmd)
    cmds.append(f"rm -f '{decode_fit_tsv}'")

    cmds.append(f"[ -f '{decode_de}' ] && [ -f '{decode_prefix}.factor.info.html' ] && [ -f '{decode_prefix}.png' ] && touch '{decode_prefix}.done'")
    mm.add_target(f"{decode_prefix}.done", [cmap_path, f"{decode_prefix}.tsv.done", f"{args.out_dir}/multi.done", f"{model_prefix}.done"], cmds)

    return f"{decode_prefix}.done"

def add_sample_json_target(mm, args, sample, sample_transcript, n_samples):
    """Add Makefile target to write the output JSON for a single sample."""
    cmds = cmd_separator([], f"Writing output JSON file for sample {sample}...")
    sample_out_dir = os.path.join(args.out_dir, "samples", sample)
    sample_out_json = os.path.join(sample_out_dir, "ficture.params.json")
    sample_feature_nohdr = os.path.join(sample_out_dir, f"{sample}.tiled.features.tsv")
    sample_feature_hdr = os.path.join(sample_out_dir, f"{sample}.tiled.features.hdr.tsv")
    sample_minmax = os.path.join(sample_out_dir, f"{sample}.tiled.coord_range.tsv")

    with flexopen(sample_transcript, "rt") as f:
        hdrs = f.readline().strip().split("\t")
        colname_feature = hdrs[args.colidx_feature-1]
        colname_count = hdrs[args.colidx_count-1]

    cmd = f"(printf '{colname_feature}\\t{colname_count}\\n'; cat {sample_feature_nohdr}) > {sample_feature_hdr}"
    cmds.append(cmd)

    summary_aux_args = []
    prerequisities = [f"{args.out_dir}/multi.done"]

    summary_aux_args_models = ["--lda-model"]
    summary_aux_args_umap = ["--umap"] if not args.skip_umap else []
    lda_runs = define_lda_runs(args, **LDA_CONFIG)
    for lda_params in lda_runs:
        train_width = lda_params["train_width"]
        n_factor = lda_params["n_factor"]
        model_id = lda_params["model_id"]
        model_prefix = os.path.join(args.out_dir, model_id)
        sample_prefix = os.path.join(sample_out_dir, f"{sample}.{model_id}")

        model_arg =",".join([
            "lda",
            model_id,
            str(train_width),
            str(n_factor),
            f"{model_prefix}.cmap.tsv",
            f"{model_prefix}.model.tsv",
            f"{model_prefix}.results.tsv.gz",  # shared fit
            f"{sample_prefix}.results.tsv.gz", # sample specific fit 
            f"{model_prefix}.bulk_chisq.tsv",
            f"{model_prefix}.factor.info.tsv"
        ])
        summary_aux_args_models.append(model_arg)
        prerequisities.append(f"{model_prefix}.done")

        if not args.skip_umap:
            umap_tsv = f"{model_prefix}.umap.tsv.gz"
            umap_png = f"{model_prefix}.umap.png"
            umap_single_prob_png = f"{model_prefix}.umap.single.prob.png"
            
            sample_umap_tsv = f"{sample_prefix}.umap.tsv.gz"
            sample_umap_png = f"{sample_prefix}.umap.png"
            sample_umap_single_prob_png = f"{sample_prefix}.umap.single.prob.png"

            prerequisities.extend([umap_tsv, umap_png, umap_single_prob_png, sample_umap_tsv, sample_umap_png, sample_umap_single_prob_png])
            summary_aux_args_umap.append(f"lda,{model_id},{umap_tsv},{umap_png},{umap_single_prob_png},{sample_umap_tsv},{sample_umap_png},{sample_umap_single_prob_png}")
    
    if len(summary_aux_args_models) > 1:
        summary_aux_args.append(" ".join(summary_aux_args_models))
    if not args.skip_umap and len(summary_aux_args_umap) > 1:
        summary_aux_args.append(" ".join(summary_aux_args_umap))

    summary_aux_args_decodes = ["--decode"]
    decode_runs = define_decode_runs(args, **LDA_CONFIG)
    for decode_params in decode_runs:
        decode_id = decode_params["decode_id"]
        decode_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{decode_id}")

        decode_arg = ",".join([
            decode_params["model_type"],
            decode_params["model_id"],
            decode_id,
            str(decode_params["fit_width"]),
            str(args.anchor_res),
            f"{decode_prefix}.tsv.gz",
            f"{decode_prefix}.png",
            f"{decode_prefix}.pseudobulk.tsv.gz",
            f"{decode_prefix}.bulk_chisq.tsv",
            f"{decode_prefix}.factor.info.tsv"
        ])
        summary_aux_args_decodes.append(decode_arg)
        prerequisities.append(f"{decode_prefix}.done")

    if len(summary_aux_args_decodes) > 1:
        summary_aux_args.append(" ".join(summary_aux_args_decodes))

    summary_cmd_parts = [
        "cartloader", "write_json_for_ficture2_multi",
        "--mode append",
        f"--in-transcript '{sample_transcript}'",
        f"--in-feature '{sample_feature_hdr}'",
        f"--in-minmax '{sample_minmax}'",
        f"--out-dir '{sample_out_dir}'",
        f"--out-json '{sample_out_json}'",
        f"--n-samples {n_samples}"
    ]
    if sample_feature_hdr:
        summary_cmd_parts.append(f"--in-feature-ficture {sample_feature_hdr}")
    summary_cmd_parts.extend(arg for arg in summary_aux_args if arg)
    cmd = " ".join(summary_cmd_parts)
    cmds.append(cmd)
    mm.add_target(sample_out_json, prerequisities, cmds)
    return sample_out_json

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

    # validate args
    if args.n_factor is None and args.pretrained_model is None:
        raise ValueError("When --pretrained-model is not provided, --n-factor is required.")
    #if "," in args.width:
    #    raise ValueError("In multi-sample pipeline, when --train-width is provided, it should be a single value, not a comma-separated list. Use --width instead.")

    # input/output/other files
    os.makedirs(args.out_dir, exist_ok=True)

    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, f"ficture.params.json")

    in_samples = []
    in_tsvs = []
    with flexopen(args.in_list, "rt") as f:
        for line in f:
            toks = line.strip().split("\t")
            in_samples.append(toks[0])
            in_tsvs.append(toks[1])

    n_samples = len(in_samples)
    
    assert os.path.exists(args.cmap_file), f"File not found: {args.cmap_file} (--cmap-file)"
    
    # apps/ficture2
    scheck_app(args.spatula)
    if not args.skip_umap:
        scheck_app(args.R)
    scheck_app(args.sort)
    scheck_app(args.gzip)
    
    ficture2bin = os.path.join(args.ficture2, "bin/punkst")
    assert os.path.exists(ficture2bin), f"File not found: {ficture2bin}. FICTURE2 Directory should include bin/punkst (--ficture2)"
    ficture2report = args.python + " '" + os.path.join(args.ficture2, "ext/py/factor_report.py") + "'"
    
    # start mm
    mm = minimake()

    # step 1. multi-sample tiling and hexagon:
    add_multisample_prepare_targets(mm, args, ficture2bin, in_samples)

    # step 2. multi-sample LDA training
    lda_runs = define_lda_runs(args, **LDA_CONFIG)
    for lda_params in lda_runs:
        # params & prefix
        train_width = lda_params["train_width"]
        n_factor = lda_params["n_factor"]

        model_id = lda_params["model_id"]
        model_prefix = os.path.join(args.out_dir, model_id)
        hex_prefix = os.path.join(args.out_dir, f"multi.hex_{train_width}")

        # files
        color_map=f"{model_prefix}.cmap.tsv"

        # 1) shared lda model
        add_lda_training_target(
            mm=mm,
            args=args,
            ficture2bin=ficture2bin,
            n_factor=n_factor,
            train_width=train_width,
            model_prefix=model_prefix,
            hex_prefix=hex_prefix,
            color_map=color_map,
            ficture2report=ficture2report
        )

        # 2) shared UMAP
        if not args.skip_umap:
            add_umap_targets(
                mm=mm,
                input_tsv=f"{model_prefix}.results.tsv.gz",
                color_map=color_map,
                out_prefix=model_prefix,
                subtitle=f"{model_id} - shared"
            )
        
        # 3) Sample-specific LDA projection 
        lda_each_targets = []
        for sample in in_samples:
            sample_lda_prefix = os.path.join(args.out_dir, "samples", sample, f"{sample}.{model_id}")
            target = add_projection_target_per_sample(
                mm=mm,
                args=args,
                ficture2bin=ficture2bin,
                model_prefix=model_prefix,
                model_id=model_id,
                sample=sample,
                train_width=train_width,
            )
            lda_each_targets.append(target)

            # 4) Sample-specific UMAP
            add_umap_targets(mm=mm, 
                            input_tsv=f"{sample_lda_prefix}.results.tsv.gz", 
                            color_map=color_map, 
                            out_prefix=sample_lda_prefix, 
                            subtitle=f"{model_id} - sample specific ({sample})")
        cmds = cmd_separator([], f"Finishing LDA projection for each sample for model {model_id}...")
        cmds.append(f"touch '{model_prefix}_each.done'")
        mm.add_target(f"{model_prefix}_each.done", lda_each_targets, cmds)

    ## step 3. multi-sample pixel-decode (perform pixel-decode for each sample)

    decode_runs = define_decode_runs(args, **LDA_CONFIG)
    for decode_params in decode_runs:            
        model_prefix = os.path.join(args.out_dir, decode_params["model_id"])
        fit_width = decode_params["fit_width"]
        decode_id = decode_params["decode_id"]

        decode_each_targets = []
        for sample in in_samples:
            target = add_pixel_decode_target_per_sample(
                mm=mm,
                args=args,
                ficture2bin=ficture2bin,
                ficture2report=ficture2report,
                model_prefix=model_prefix,
                model_path=decode_params["model_path"],
                cmap_path=decode_params["cmap_path"],
                decode_id=decode_id,
                fit_width=fit_width,
                n_factor=decode_params["n_factor"],
                fit_n_move=int(fit_width / args.anchor_res),
                sample=sample,
            )
            decode_each_targets.append(target)
        cmds = cmd_separator([], f"Finishing decode, ID: {decode_id}")
        cmds.append(f"touch {args.out_dir}/{decode_id}_each.done")
        mm.add_target(f"{args.out_dir}/{decode_id}_each.done", decode_each_targets, cmds)

    ## step 4. write the output JSON file for each sample
    json_each_targets = []
    for sample, sample_transcript in zip(in_samples, in_tsvs):
        sample_out_json = add_sample_json_target(mm, args, sample, sample_transcript, n_samples)
        json_each_targets.append(sample_out_json)

    cmds=cmd_separator([], f"Finishing writing the JSON file for each sample...")
    cmds.append(f"touch '{args.out_dir}/multi_json_each.done'")
    mm.add_target(f"{args.out_dir}/multi_json_each.done", json_each_targets, cmds)

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

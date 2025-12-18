import os
from cartloader.utils.utils import flexopen, cmd_separator

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

create_umap_rscript = f"{repo_dir}/cartloader/r/create_umap.r"
draw_umap_rscript = f"{repo_dir}/cartloader/r/draw_umap.r"
draw_umap_single_rscript = f"{repo_dir}/cartloader/r/draw_umap_single.r"

def _parse_int_csv(value):
    return [int(x) for x in value.split(",")] if value else []


def _fit_widths(args, train_width):
    if hasattr(args, "fit_width") and getattr(args, "fit_width") is not None:
        return _parse_int_csv(getattr(args, "fit_width"))
    return [train_width]


def define_lda_runs(
    args,
    *,
    default_n_factor="12",
    require_init_flag=False,
    allow_pretrained=False,
    width_error_msg=None,
):
    if require_init_flag:
        assert getattr(args, "init_lda"), "--init-lda must be ON when running define_lda_runs()"
    if width_error_msg is None:
        width_error_msg = "When --init-lda is ON, provide at least one train width for LDA training using --width"
    assert args.width is not None, width_error_msg

    train_widths = _parse_int_csv(args.width)

    if allow_pretrained and getattr(args, "pretrained_model", None) is not None:
        if args.n_factor is not None:
            raise ValueError("When --pretrained-model is provided, --n-factor should not be provided.")
        with flexopen(args.pretrained_model) as rf:
            hdrs = rf.readline().rstrip().split("\t")
            n_factor = len(hdrs) - 1
        n_factors = [n_factor]
    else:
        if args.n_factor is None:
            print(f"Warning: --n-factor is not provided for LDA. Using default values: {default_n_factor}")
            args.n_factor = default_n_factor
        n_factors = _parse_int_csv(args.n_factor)

    train_params = [
        {
            "model_type": "lda",
            "train_width": train_width,
            "n_factor": n_factor,
            "model_id": f"t{train_width}_f{n_factor}",
        }
        for train_width in train_widths
        for n_factor in n_factors
    ]
    return train_params


def define_training_runs(args, **kwargs):
    return define_lda_runs(args, **kwargs)

def define_decode_runs(args, **kwargs):
    decode_runs = []
    train_params = define_training_runs(args, **kwargs)
    for train_param in train_params:
        model_type = train_param["model_type"]
        train_width = train_param["train_width"]
        n_factor = train_param["n_factor"]
        model_id = train_param["model_id"]
        model_prefix = os.path.join(args.out_dir, model_id)
        model_path = f"{model_prefix}.model.tsv"
        fit_widths = _fit_widths(args, train_width)
        for fit_width in fit_widths:
            decode_id = f"{model_id}_p{fit_width}_a{args.anchor_res}"
            cmap_path = f"{model_prefix}.cmap.tsv"
            decode_runs.append({
                "model_type": model_type,
                "model_id": model_id,
                "model_path": model_path,
                "decode_id": decode_id,
                "n_factor": n_factor,
                "fit_width": fit_width,
                "cmap_path": cmap_path,
                "prerequisite_path": f"{model_prefix}.done",
            })
    return decode_runs

def add_umap_targets(
    mm,
    input_tsv,
    color_map,
    out_prefix,
    subtitle
):
    """Add Makefile targets for UMAP generation and visualization for a model."""
    umap_tsv = f"{out_prefix}.umap.tsv.gz"
    umap_png = f"{out_prefix}.umap.png"
    umap_single_prob_png = f"{out_prefix}.umap.single.prob.png"

    cmds = cmd_separator([], f"UMAP for ID: {subtitle}...")
    cmd = " ".join([
        f"Rscript {create_umap_rscript}",
        f"--input {input_tsv}",
        f"--out-prefix {out_prefix}"
    ])
    cmds.append(cmd)
    mm.add_target(umap_tsv, [f"{out_prefix}.done", color_map], cmds)

    cmds = cmd_separator([], f"UMAP Visualization for ID: {subtitle}...")
    cmd = " ".join([
        f"Rscript {draw_umap_rscript}",
        f"--input {umap_tsv}",
        f"--out-prefix {out_prefix}",
        f"--cmap {color_map}",
        f'--subtitle "{subtitle}"'
    ])
    cmds.append(cmd)
    mm.add_target(umap_png, [f"{out_prefix}.done", color_map, umap_tsv], cmds)

    cmds = cmd_separator([], f"UMAP Visualization (plot for individual factors; colorized by probability) for ID: {subtitle}...")
    cmd = " ".join([
        f"Rscript {draw_umap_single_rscript}",
        f"--input {umap_tsv}",
        f"--out-prefix {out_prefix}",
        f"--cmap {color_map}",
        f'--subtitle "{subtitle}"',
        f"--mode prob"
    ])
    cmds.append(cmd)
    mm.add_target(umap_single_prob_png, [f"{out_prefix}.done", color_map, umap_tsv], cmds)

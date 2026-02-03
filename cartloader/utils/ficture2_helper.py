import os
from cartloader.utils.utils import flexopen, cmd_separator, factor_id_to_name

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
        f"Rscript '{create_umap_rscript}'",
        f"--input '{input_tsv}'",
        f"--out-prefix '{out_prefix}'"
    ])
    cmds.append(cmd)
    mm.add_target(umap_tsv, [f"{out_prefix}.done", color_map], cmds)

    cmds = cmd_separator([], f"UMAP Visualization for ID: {subtitle}...")
    cmd = " ".join([
        f"Rscript '{draw_umap_rscript}'",
        f"--input '{umap_tsv}'",
        f"--out-prefix '{out_prefix}'",
        f"--cmap '{color_map}'",
        f'--subtitle "{subtitle}"'
    ])
    cmds.append(cmd)
    mm.add_target(umap_png, [f"{out_prefix}.done", color_map, umap_tsv], cmds)

    cmds = cmd_separator([], f"UMAP Visualization (plot for individual factors; colorized by probability) for ID: {subtitle}...")
    cmd = " ".join([
        f"Rscript '{draw_umap_single_rscript}'",
        f"--input '{umap_tsv}'",
        f"--out-prefix '{out_prefix}'",
        f"--cmap '{color_map}'",
        f'--subtitle "{subtitle}"',
        f"--mode prob"
    ])
    cmds.append(cmd)
    mm.add_target(umap_single_prob_png, [f"{out_prefix}.done", color_map, umap_tsv], cmds)

## transform FICTURE parameters to FACTOR assets (new standard)
def ficture2_params_to_factor_assets(params, skip_raster=False, cell_params = None):
    ## model_id
    ## proj_params -> proj_id
    ## proj_params -> decode_params -> decode_id
    suffix_factormap = "-factor-map.tsv"
    suffix_de = "-bulk-de.tsv"
    suffix_info = "-info.tsv"
    suffix_model = "-model.tsv"
    suffix_post = "-pseudobulk.tsv.gz"
    suffix_cells_pmtiles = "-cells.pmtiles"
    suffix_boundaries_pmtiles = "-boundaries.pmtiles"
    suffix_rgb = "-rgb.tsv"
    suffix_hex_coarse = ".pmtiles"
    suffix_raster = "-pixel-raster.pmtiles"
    suffix_umap_tsv = "-umap.tsv.gz"
    suffix_umap_pmtiles = "-umap.pmtiles"
    suffix_umap_png = ".umap.png"
    suffix_umap_ind_png = ".umap.single.prob.png"
    suffix_heatmap_pdf = "-heatmap.pdf"
    suffix_heatmap_tsv = "-heatmap.tsv"

    def _create_decode_asset(decode_param):
        if not "decode_id" in decode_param:
            raise ValueError(f"decode_id is missing from FICTURE parameter for {model_id}")
        decode_id = decode_param["decode_id"].replace("_", "-")

        asset = {
            "id": model_id,
            "name": factor_id_to_name(model_id),
            "model_id": model_id,
            "decode_id": decode_id,
            "de": (decode_id if is_multi_sample else model_id) + suffix_de,
            "info": (decode_id if is_multi_sample else model_id) + suffix_info,
            "model": model_id + suffix_model,
            "post": decode_id + suffix_post,
            "rgb": model_id + suffix_rgb,
            "pmtiles": {
                #"hex_coarse": model_id + suffix_hex_coarse,
                "hex": model_id + suffix_hex_coarse,
                **({"raster": decode_id + suffix_raster} if not skip_raster else {})
            }
        }
        if is_multi_sample:
            asset["shared_de"] = model_id + suffix_de
            asset["shared_post"] = model_id + suffix_model
            asset["shared_info"] = model_id + suffix_info

        if "factor_map" in param:
            asset["factor_map"] = model_id + suffix_factormap
        
        return asset

    out_assets = []
    for param in params: ## train_params is a list of dictionaries
        # print(param)
        # print(param.get("umap", "NO UMAP"))
        if not "model_id" in param:
            raise ValueError(f"model_id is missing from FICTURE parameters")
        model_id = param["model_id"].replace("_", "-")

        len_decode_params = len(param["decode_params"])
        umap_params = param.get("umap",{})
        is_multi_sample = param.get("analysis_type") == "multi-sample"

        # construct out_assets 
        out_asset={}
        if len_decode_params == 0: ## train_param only
            out_asset = {
                "id": model_id,
                "name": factor_id_to_name(model_id),
                "model_id": model_id,
                "de": model_id + suffix_de,
                "info": model_id + suffix_info,
                "model": model_id + suffix_model,
                "rgb": model_id + suffix_rgb,
                "pmtiles": {
                    #"hex_coarse": model_id + suffix_hex_coarse
                    "hex": model_id + suffix_hex_coarse
                }
            }
            if "factor_map" in param:
                out_asset["factor_map"] = model_id + suffix_factormap
        elif len_decode_params == 1:
            decode_param = param["decode_params"][0]
            out_asset = _create_decode_asset(decode_param)
        else: ## multiple decode_params
            for decode_param in param["decode_params"]:
                out_asset = _create_decode_asset(decode_param)
        
        # add analysis info
        if param.get("analysis_type") == "multi-sample":
            out_asset["analysis_type"] = "multi-sample"

        # add umap
        if umap_params:
            if is_multi_sample:
                # Sample-specific UMAP
                out_asset["umap"] = {
                    "tsv": model_id + suffix_umap_tsv,
                    "pmtiles": model_id + suffix_umap_pmtiles,
                    "png": model_id + suffix_umap_png,
                    "ind_png": model_id + suffix_umap_ind_png,
                }
                # Shared UMAP
                out_asset["shared_umap"] = {
                    "tsv": model_id + "-shared" + suffix_umap_tsv,
                    "pmtiles": model_id + "-shared" + suffix_umap_pmtiles,
                    "png": model_id + "-shared" + suffix_umap_png,
                    "ind_png": model_id + "-shared" + suffix_umap_ind_png,
                }
            else:
                out_asset["umap"] = {
                    "tsv": model_id + suffix_umap_tsv,
                    "pmtiles": model_id + suffix_umap_pmtiles,
                    "png": model_id + suffix_umap_png,
                    "ind_png": model_id + suffix_umap_ind_png,
                }
        # append
        out_assets.append(out_asset)

    # process cell_params if provided
    if cell_params is not None:
        for cell_param in cell_params:
            model_id = cell_param["model_id"]
            model_rgb = cell_param["cmap"]
            cell_xy_f = cell_param["cell_xy_path"]
            cell_boundaries_f = cell_param.get("cell_boundaries_path", None)
            cell_clust_f = cell_param.get("cluster_path", None)
            cell_info_f = cell_param.get("cluster_info", None)
            model_manifolds = cell_param.get("manifolds", [])
            cell_de_tsvf = cell_param["cluster_de"]
            cell_post_tsvf = cell_param["cluster_pseudobulk"]
            cell_pixel_tsvf = cell_param["pixel_tsv_path"]
            cell_pixel_pngf = cell_param["pixel_png_path"]
            cell_heatmap_pdf = cell_param["cluster_model_heatmap_pdf"]
            cell_heatmap_tsv = cell_param["cluster_model_heatmap_tsv"]
            out_asset = {
                "id": model_id,
                "name": factor_id_to_name(model_id),
                "model_id": model_id,
                "decode_id": model_id,
                "cells_id": model_id,
                "de": model_id + suffix_de,
                "info": model_id + suffix_info,
                "post": model_id + suffix_post,
                "rgb": model_id + suffix_rgb,
                "heatmap_pdf": model_id + suffix_heatmap_pdf,
                "heatmap_tsv": model_id + suffix_heatmap_tsv,
                "pmtiles": {
                    "cells": model_id + suffix_cells_pmtiles,
                    **({"boundaries": model_id + suffix_boundaries_pmtiles} if cell_boundaries_f is not None else {}),
                    **({"raster": model_id + suffix_raster} if not skip_raster else {})
                }
            }
            if "sample" in model_manifolds and "umap" in model_manifolds["sample"]:
                out_asset["umap"] = {
                    "tsv": model_id + suffix_umap_tsv,
                    "pmtiles": model_id + suffix_umap_pmtiles,
                    "png": model_id + suffix_umap_png,
                }
            if cell_param.get("analysis_type") == "multi-sample":
                out_asset["analysis_type"] = "multi-sample"
                out_asset["shared_de"] = model_id + "-shared" + suffix_de
                out_asset["shared_info"] = model_id + "-shared" + suffix_info
                out_asset["shared_post"] = model_id + "-shared" + suffix_post
                out_asset["shared_heatmap_pdf"] = model_id + "-shared" + suffix_heatmap_pdf
                out_asset["shared_heatmap_tsv"] = model_id + "-shared" + suffix_heatmap_tsv
                if "shared" in model_manifolds and "umap" in model_manifolds["shared"]:
                    out_asset["shared_umap"] = {
                        "tsv": model_id + "-shared" + suffix_umap_tsv,
                        "pmtiles": model_id + "-shared" + suffix_umap_pmtiles,
                        "png": model_id + "-shared" + suffix_umap_png,
                    }
            out_assets.append(out_asset)
    return out_assets

# ## transform FICTURE parameters to FACTOR assets (new standard)
# def ficture2_params_to_factor_assets(params, skip_raster=False):
#     ## model_id
#     ## proj_params -> proj_id
#     ## proj_params -> decode_params -> decode_id
#     suffix_factormap = "-factor-map.tsv"
#     suffix_de = "-bulk-de.tsv"
#     suffix_info = "-info.tsv"
#     suffix_model = "-model.tsv"
#     suffix_post = "-pseudobulk.tsv.gz"
#     suffix_rgb = "-rgb.tsv"
#     suffix_hex_coarse = ".pmtiles"
#     suffix_raster = "-pixel-raster.pmtiles"
#     suffix_umap_tsv = "-umap.tsv.gz"
#     suffix_umap_pmtiles = "-umap.pmtiles"
#     suffix_umap_png = ".umap.png"
#     suffix_umap_ind_png = ".umap.single.prob.png"

#     def _create_decode_asset(decode_param):
#         if not "decode_id" in decode_param:
#             raise ValueError(f"decode_id is missing from FICTURE parameter for {model_id}")
#         decode_id = decode_param["decode_id"].replace("_", "-")

#         asset = {
#             "id": model_id,
#             "name": factor_id_to_name(model_id),
#             "model_id": model_id,
#             "decode_id": decode_id,
#             "de": (decode_id if is_multi_sample else model_id) + suffix_de,
#             "info": (decode_id if is_multi_sample else model_id) + suffix_info,
#             "model": model_id + suffix_model,
#             "post": decode_id + suffix_post,
#             "rgb": model_id + suffix_rgb,
#             "pmtiles": {
#                 "hex_coarse": model_id + suffix_hex_coarse,
#                 **({"raster": decode_id + suffix_raster} if not skip_raster else {})
#             }
#         }
#         if is_multi_sample:
#             asset["shared_de"] = model_id + suffix_de
#             asset["shared_post"] = model_id + suffix_model
#             asset["shared_info"] = model_id + suffix_info

#         if "factor_map" in param:
#             asset["factor_map"] = model_id + suffix_factormap
        
#         return asset

#     out_assets = []
#     for param in params: ## train_params is a list of dictionaries
#         # print(param)
#         # print(param.get("umap", "NO UMAP"))
#         if not "model_id" in param:
#             raise ValueError(f"model_id is missing from FICTURE parameters")
#         model_id = param["model_id"].replace("_", "-")

#         len_decode_params = len(param["decode_params"])
#         umap_params = param.get("umap",{})
#         is_multi_sample = param.get("analysis_type") == "multi-sample"

#         # construct out_assets 
#         out_asset={}
#         if len_decode_params == 0: ## train_param only
#             out_asset = {
#                 "id": model_id,
#                 "name": factor_id_to_name(model_id),
#                 "model_id": model_id,
#                 "de": model_id + suffix_de,
#                 "info": model_id + suffix_info,
#                 "model": model_id + suffix_model,
#                 "rgb": model_id + suffix_rgb,
#                 "pmtiles": {
#                     "hex_coarse": model_id + suffix_hex_coarse
#                 }
#             }
#             if "factor_map" in param:
#                 out_asset["factor_map"] = model_id + suffix_factormap
#         elif len_decode_params == 1:
#             decode_param = param["decode_params"][0]
#             out_asset = _create_decode_asset(decode_param)
#         else: ## multiple decode_params
#             for decode_param in param["decode_params"]:
#                 out_asset = _create_decode_asset(decode_param)
        
#         # add analysis info
#         if param.get("analysis_type") == "multi-sample":
#             out_asset["analysis_type"] = "multi-sample"

#         # add umap
#         if umap_params:
#             if is_multi_sample:
#                 # Sample-specific UMAP
#                 out_asset["umap"] = {
#                     "tsv": model_id + suffix_umap_tsv,
#                     "pmtiles": model_id + suffix_umap_pmtiles,
#                     "png": model_id + suffix_umap_png,
#                     "ind_png": model_id + suffix_umap_ind_png,
#                 }
#                 # Shared UMAP
#                 out_asset["shared_umap"] = {
#                     "tsv": model_id + "-shared" + suffix_umap_tsv,
#                     "pmtiles": model_id + "-shared" + suffix_umap_pmtiles,
#                     "png": model_id + "-shared" + suffix_umap_png,
#                     "ind_png": model_id + "-shared" + suffix_umap_ind_png,
#                 }
#             else:
#                 out_asset["umap"] = {
#                     "tsv": model_id + suffix_umap_tsv,
#                     "pmtiles": model_id + suffix_umap_pmtiles,
#                     "png": model_id + suffix_umap_png,
#                     "ind_png": model_id + suffix_umap_ind_png,
#                 }
#         # append
#         out_assets.append(out_asset)

#     return out_assets


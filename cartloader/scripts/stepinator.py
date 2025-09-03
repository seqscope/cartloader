# TODO: Currently stepinator can only allow dry-run for one step. For example, if requested --run-ficture and --run-cartload-join with --dry-run, it can only process --run-ficture while 
# the dry-run for --run-cartload-join will failed given its input is from run-ficture and run-ficture is not executed.

## 250731: Need to update the cmds related to images.
import os
import sys
import subprocess
import argparse
import yaml
import datetime
from types import SimpleNamespace

from cartloader.utils.utils import add_param_to_cmd, scheck_file
from cartloader.utils.image_helper import update_orient_in_histology
from cartloader.utils.sge_helper import aux_sge_args
from cartloader.utils.execution import write_jobfile, submit_job

timestamp = datetime.datetime.now().strftime("%y%m%d")

img_suffix=[".tif", ".tiff", ".png"]
img_suffix.extend([suffix.upper() for suffix in img_suffix])

aux_env_args = {
        "sge_stitch": ["spatula"],
        "sge_convert": ["spatula", "gzip", "parquet_tools"],
        "north_up": ["gdal_translate", "gdalwarp"],
        "image_stitch": ["gdal_translate", "gdalbuildvrt"],
        "run_ficture": ['bgzip', "tabix", "gzip", "sort", "sort_mem"],
        "run_cartload": ['gzip', 'pmtiles', 'gdal_translate', 'gdaladdo', 'tippecanoe', 'spatula'],
        "image_png2pmtiles": ['pmtiles', 'gdal_translate', 'gdaladdo', "gdalinfo"],
        "upload_aws": ['aws']
}

aux_env_args["run_ficture1"]= aux_env_args["run_ficture"]
aux_env_args["run_ficture2"]= aux_env_args["run_ficture"] + ["ficture2", "python", "spatula"]

aux_params_args = {
    "sge_stitch": ["colname_feature_name", "colname_feature_id", "colname_x", "colname_y"] + ["radius", "quartile", "hex_n_move", "polygon_min_size"],
    "sge_convert": [item for sublist in aux_sge_args.values() for item in sublist] + ["radius", "quartile", "hex_n_move", "polygon_min_size"],
    "north_up": ["srs", "resample"],
    "run_ficture": [ 'anchor_res', 'radius_buffer', 
                     'min_ct_per_unit_hexagon', "min_ct_per_unit_train",
                     'minibatch_size', 
                     'train_epoch',
                     'fit_width', 
                     'min_ct_per_feature', 'de_max_pval', 'de_min_fold',
                     "include_feature_regex", "exclude_feature_regex"]
    }

aux_params_args["run_ficture1"] = aux_params_args["run_ficture"] + ['hexagon_n_move', 'hexagon_precision', 
                                                                    'minibatch_buffer',
                                                                    "hexagon_width_10x",
                                                                    'train_epoch_id_len', 'lda_rand_init', 'lda_plot_um_per_pixel',
                                                                    'fit_precision',  'min_ct_per_unit_fit', 'fit_plot_um_per_pixel',
                                                                    'decode_top_k', 'decode_block_size', 'decode_scale', 'decode_precision', 'decode_plot_um_per_pixel',
                                                                    'merge_max_dist_um', 'merge_max_k', 'merge_max_p',
                                                                    "include_feature_list", "exclude_feature_list", "include_feature_type_regex", "feature_type_ref", "feature_type_ref_colidx_name", "feature_type_ref_colidx_type"
                                                                   ]
aux_params_args["run_ficture2"] = aux_params_args["run_ficture"]

img_keys = [ "path",
            "ome", "lower_thres_quantile", "upper_thres_quantile", "level", "colorize",
            "georeference", "georef_tsv", "georef_bounds", 
            "rotate", "flip"
    ]
img_keys_tiles =  ["row", "col"] + [x for x in img_keys if x not in ["ome", "lower_thres_quantile", "upper_thres_quantile", "level", "colorize"]]

def merge_config(base_config, args, keys, prefix=None):
    """
    Merges parameters from a base configuration dictionary and command-line arguments.
    Args:
        base_config (dict): Dictionary containing default values.
        args (argparse.Namespace): Parsed command-line arguments.
        keys (list): Keys to be merged from args and base_config.
        prefix (str, optional): Prefix to be added to keys in args.
    Returns:
        SimpleNamespace: Merged configuration.
    """
    config = base_config.get(prefix, {}).copy() if prefix else base_config.copy()
    for key in keys:
        val = getattr(args, f"{prefix}_{key}" if prefix else key, None)
        #print(f"{prefix}_{key}" if prefix else key)
        if isinstance(val, str) and val is not None:
            config[key] = val
        elif isinstance(val, list) and len(val) > 0:
            config[key] = val
    return SimpleNamespace(**config)

def cmd_sge_stitch(sgeinfo, args, env, generate_tile_minmax_only=False):
    mkbn = "sge_stitch" if args.mk_id is None else f"sge_stitch_{args.mk_id}"
    # collect tiles 
    stitch_cmd = " ".join([
        "cartloader", "sge_stitch",
        f"--makefn {mkbn}.mk",
        f"--in-tiles {sgeinfo['in_tiles_str']}",
        f"--out-dir {sgeinfo['sge_dir']}",
        f"--units-per-um {sgeinfo.get('units_per_um', None)}" if sgeinfo.get("units_per_um", None) else "",
        f"--colname-count {sgeinfo.get('colname_count', None)}" if sgeinfo.get('colname_count', None) else "",
        f"--sge-visual" if args.sge_visual else "",
        f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        f"--restart" if args.restart else "",
        "--generate-tile-minmax-only" if generate_tile_minmax_only else "",
        f"--filter-by-density --out-filtered-prefix {sgeinfo['filtered_prefix']} --genomic-feature {sgeinfo['colname_count']}" if sgeinfo['filter_by_density'] else "",
        f"--north-up --out-northup-tif {sgeinfo['xy_ntif']}" if args.north_up else "",
    ])

    # add aux stitch parameters
    stitch_aug = merge_config(sgeinfo, args, aux_params_args["sge_stitch"], prefix=None)
    stitch_cmd = add_param_to_cmd(stitch_cmd, stitch_aug, aux_params_args["sge_stitch"])

    # add aux north-up parameters
    if args.north_up:
        northup_aug = merge_config(sgeinfo, args, aux_params_args["north_up"], prefix=None)
        stitch_cmd = add_param_to_cmd(stitch_cmd, northup_aug, aux_params_args["north_up"])
    
    # add aux env
    stitch_cmd = add_param_to_cmd(stitch_cmd, env, aux_env_args["sge_stitch"]+aux_env_args["north_up"]) if args.north_up else add_param_to_cmd(stitch_cmd, env, aux_env_args["sge_stitch"])
    return stitch_cmd

def cmd_sge_convert(sgeinfo, args, env):
    mkbn = "sge_convert" if args.mk_id is None else f"sge_convert_{args.mk_id}"

    platform= sgeinfo.get("platform", None)
    assert platform is not None, "Please provide platform information for sge_convert"

    if platform == "10x_visium_hd":
        assert sgeinfo.get("in_mex", None) is not None, "Please provide --in-mex for 10x_visium_hd"
        assert sgeinfo.get("in_parquet", None) is not None, "Please provide --in-parquet for 10x_visium_hd"
        in_arg= " ".join([
            f"--in-mex {sgeinfo['in_mex']}", 
            f"--in-parquet {sgeinfo['in_parquet']}",
            f"--scale-json {sgeinfo['scale_json']}" if sgeinfo.get("scale_json", None) is not None else "",
        ])
    elif platform == "seqscope":
        assert sgeinfo.get("in_mex", None) is not None, "Please provide --in-mex for seqscope"
        in_arg= f"--in-mex {sgeinfo['in_mex']}"
    elif platform in ["10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st"]:
        assert sgeinfo.get("in_csv", None) is not None, f"Please provide --in-csv for {platform}"
        in_arg= f"--in-csv {sgeinfo['in_csv']}"
    elif platform in ["generic"]:
        assert sgeinfo.get("in_csv", None) is not None, f"Please provide --in-csv for {platform}"
        in_arg= f"--in-csv {sgeinfo['in_csv']}"

    sgeinfo["csv_colnames_others"]= sgeinfo["csv_colnames_others"].split(",") if sgeinfo.get("csv_colnames_others", None) is not None else []

    format_cmd = " ".join([
        "cartloader", "sge_convert",
        f"--makefn {mkbn}.mk",
        f"--platform {sgeinfo['platform']}",
        in_arg,
        f"--out-dir {sgeinfo['sge_dir']}",
        f"--colname-count {sgeinfo['colname_count']}",
        f"--filter-by-density --out-filtered-prefix {sgeinfo['filtered_prefix']} --genomic-feature {sgeinfo['colname_count']}" if sgeinfo['filter_by_density'] else "",
        f"--sge-visual" if args.sge_visual else "",
        f"--n-jobs {args.n_jobs}" if args.n_jobs else ""
        f"--north-up --out-northup-tif {sgeinfo['xy_ntif']}" if args.north_up else "",
    ])
    # add aux parameters
    # remove the "units_per_um" in the aux_params_args["sge_convert"] 
    format_aug = merge_config(sgeinfo, args, aux_params_args["sge_convert"], prefix=None)
    format_cmd = add_param_to_cmd(format_cmd, format_aug, aux_params_args["sge_convert"])

    # add aux north-up parameters
    if args.north_up:
        northup_aug = merge_config(sgeinfo, args, aux_params_args["north_up"], prefix=None)
        format_cmd = add_param_to_cmd(format_cmd, northup_aug, aux_params_args["north_up"])

    # add aux tools
    format_cmd = add_param_to_cmd(format_cmd, env, aux_env_args["sge_convert"]+aux_env_args["north_up"]) if args.north_up else add_param_to_cmd(format_cmd, env, aux_env_args["sge_convert"])

    return format_cmd

def define_sge2fn(filter_by_density, filtered_prefix):
    if filter_by_density:
        assert filtered_prefix is not None, "Error: --filtered-prefix is Required when --filter-by-density is applied"
    tsvfn    = "transcripts.unsorted.tsv.gz" if not filter_by_density else f"{filtered_prefix}.transcripts.unsorted.tsv.gz"
    cstsvfn  = "transcripts.sorted.tsv.gz" if not filter_by_density else f"{filtered_prefix}.transcripts.sorted.tsv.gz"
    ftrfn    = "feature.clean.tsv.gz" if not filter_by_density else f"{filtered_prefix}.feature.lenient.tsv.gz"
    minmaxfn = "coordinate_minmax.tsv" if not filter_by_density else f"{filtered_prefix}.coordinate_minmax.tsv"

    sge2fn={
        "tsv": tsvfn,
        "ftr": ftrfn,
        "minmax": minmaxfn,
        "cstsv": cstsvfn,
    }
    return sge2fn

def define_sge_arg(sgefn, sge_dir, fic_dir, fic_v):
    # if sge_dir is None:
    #     sge_dir = run_dir
    tsv = os.path.join(sge_dir, sgefn["tsv"])
    ftr = os.path.join(sge_dir, sgefn["ftr"])
    minmax = os.path.join(sge_dir, sgefn["minmax"])
    
    cstsv = os.path.join(fic_dir, sgefn["cstsv"])

    if fic_v == "1":
        sge_arg=f"--in-transcript {tsv} --in-feature {ftr} --in-minmax {minmax} --in-cstranscript {cstsv}"
    else:
        sge_arg=f"--in-transcript {tsv} --in-feature {ftr} --in-minmax {minmax}"
    return sge_arg

def link_sge_to_fict(sge_fn, sge_dir, fic_dir):
    # create softlink
    dst2src={}
    for infn in [sge_fn["tsv"], sge_fn["ftr"], sge_fn["minmax"]]:
        src = os.path.join(sge_dir, infn) # source
        dst = os.path.join(fic_dir, infn) # destination
        # if infn == cstsvfn and not os.path.exists(os.path.abspath(src)):
        #     continue
        if not os.path.isfile(dst) and not os.path.exists(dst):
            if os.path.exists(src):
                os.symlink(src, dst)
                print(f"Creating symlink for {infn}")
            else:
                dst2src[dst]=src
    return dst2src

def cmd_run_ficture(run_i, fic_v, args, env):
    ficture_cmds = []
    ext_path = run_i["ext_path"]
    ext_id   = run_i["ext_id"]
    
    # makefile 
    mkbn = f"run_ficture{fic_v}" if ext_path is None else f"run_ficture{fic_v}_{ext_id}"
    mkbn = f"{mkbn}_{args.mk_id}" if args.mk_id is not None else mkbn

    fic_dir=os.path.join(run_i["run_dir"], "ficture") if fic_v == "1" else os.path.join(run_i["run_dir"], "ficture2")

    # sge
    sge2fn = define_sge2fn(run_i["filter_by_density"], run_i["filtered_prefix"])
    sge_arg = define_sge_arg(sge2fn, run_i["sge_dir"], fic_dir, fic_v)
    in_dist = os.path.join(run_i["sge_dir"], "feature.distribution.tsv.gz")

    # * check files
    if not args.sge_convert and not args.sge_stitch:
        if not args.lenient:
            for file in [sge2fn["tsv"], sge2fn["ftr"], sge2fn["minmax"]]:
                scheck_file(os.path.join(run_i["sge_dir"], file))
    # * link files from sge to ficture
    if run_i["sge_dir"] is not None and run_i["sge_dir"] != fic_dir:
        sge_dst2src = link_sge_to_fict(sge2fn, run_i["sge_dir"], fic_dir)
    # * if sge_dst2src is not an empty dict, add the files to the command
    if sge_dst2src:
        for dst, src in sge_dst2src.items():
            ficture_cmds.append(f"ln -s {src} {dst}")
    
    # key params
    assert run_i["train_width"] is not None, "Error: When --run-ficture1 or --run-ficture2, provide width"
    assert not (fic_v == "2" and ext_path), "Error: --run-ficture2 does not support --ext-path"
    if ext_path:
        assert ext_id is not None, "Error: --ext-id is required when running ficture with an external model"
        assert os.path.exists(ext_path), f"Error: --ext-path is defined with a missing file: {ext_path}"
    else:
        assert run_i["n_factor"] is not None, "Error: --n-factor is required when running --run-ficture1/--run-ficture2 with LDA"

    if fic_v == "1":
        # cmap
        if run_i["cmap"] is None:
            cmap_arg = ""
        else:
            if not args.lenient:
                scheck_file(run_i["cmap"])
            cmap_arg = f"--cmap-static --static-cmap-file {run_i['cmap']}"
        # cmd
        ficture_cmd = " ".join([
            "cartloader", "run_ficture",
            f"--makefn {mkbn}.mk",
            f"--out-dir {fic_dir}",
            sge_arg,
            f"--major-axis {run_i['major_axis']}",
            f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
            f"{'--init-ext' if ext_path and args.init_ext else '--main-ext' if ext_path else '--main'}",
            f"--ext-path {ext_path}" if ext_path else "",
            f"--ext-id {ext_id}" if ext_path else "",
            "--copy-ext-model" if ext_path and run_i["copy_ext_model"] else "",
            f"--n-factor {run_i['n_factor']}" if ext_path is None else "",
            f"--train-width {run_i['train_width']}",
            "--skip-coarse-report" if args.skip_coarse_report else "",
            "--segment-10x" if args.segment_10x else "",
            f"--filter-by-overlapping-features --in-feature-dist {in_dist}" if run_i["filter_by_overlapping_features"] else "",
            f"--min-ct-per-ftr-tile {run_i['min_ct_per_ftr_tile']}" if run_i["filter_by_overlapping_features"] and run_i["min_ct_per_ftr_tile"] > 0 else "",
            cmap_arg,
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
            f"--threads {args.threads}"
        ])
    else:
        if run_i["cmap"] is None:
            cmap_arg = ""
        else:
            if not args.lenient:
                scheck_file(run_i["cmap"])
            cmap_arg = f"--cmap-file {run_i['cmap']}" 
        ficture_cmd = " ".join([
            "cartloader", "run_ficture2",
            f"--makefn {mkbn}.mk",
            f"--out-dir {fic_dir}",
            sge_arg,
            f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
            "--main",
            f"--width {run_i['train_width']}",
            f"--n-factor {run_i['n_factor']}",
            f"--filter-by-overlapping-features --in-feature-dist {in_dist}" if run_i["filter_by_overlapping_features"] else "",
            f"--min-ct-per-ftr-tile {run_i['min_ct_per_ftr_tile']}" if run_i["filter_by_overlapping_features"] and run_i["min_ct_per_ftr_tile"] > 0 else "",
            cmap_arg,
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
            f"--threads {args.threads}"
        ])
    
    # add aux env/tools
    ficture_cmd = add_param_to_cmd(ficture_cmd, env, aux_env_args[f"run_ficture{fic_v}"])
    # add aux parameters
    ficture_aug = merge_config(run_i, args, aux_params_args[f"run_ficture{fic_v}"], prefix=None)  # merge auxiliary parameters
    ficture_cmd = add_param_to_cmd(ficture_cmd, ficture_aug, aux_params_args[f"run_ficture{fic_v}"])
    print(ficture_cmd)
    ficture_cmds.append(ficture_cmd)
    return ficture_cmds

def cmd_run_cartload(run_i, cartl_v, args, env):

    mkbn=f"run_cartload{cartl_v}" if args.mk_id is None else f"run_cartload{cartl_v}_{args.mk_id}"

    fic_dir = os.path.join(run_i["run_dir"], ("ficture" if cartl_v == "1" else "ficture2"))
    cartload_dir = os.path.join(run_i["run_dir"], ("cartload" if cartl_v == "1" else "cartload2"))
    os.makedirs(cartload_dir, exist_ok=True)

    if cartl_v == "1":
        cartload_cmd=" ".join([
            "cartloader", "run_cartload_join", 
            f"--makefn {mkbn}.mk",
            f"--fic-dir {fic_dir}", 
            f"--out-dir {cartload_dir}",
            f"--id {run_i['run_id']}",
            f"--major-axis {run_i['major_axis']}",
            f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
            f"--threads {args.threads}"
        ])
    else:
        cartload_cmd=" ".join([
            "cartloader", "run_cartload2", 
            f"--makefn {mkbn}.mk",
            f"--fic-dir {fic_dir}", 
            f"--out-dir {cartload_dir}",
            f"--id {run_i['run_id']}",
            f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
            f"--threads {args.threads}"
        ])
    # add aux env/tools
    cartload_cmd = add_param_to_cmd(cartload_cmd, env, aux_env_args[f"run_cartload"], underscore2dash=False)
    return cartload_cmd

# img_keys and img_keys_tiles are global variables
def parse_image_args(histology_args, by_tiles=False):
    # <image_id>;<path>;<transform>;....;<georeference>;<georef_tsv>;<georef_bounds>;<rotate_degree>;<flip_direction>
    # <image_id>;<tile_info1>;<tile_info2>... where each tile_info is <row>,<col>,<path>,<georeference>,<georef_tsv>,<georef_bounds>,<rotate_degree>,<flip_direction>. 
    img_input = []
    for histology in histology_args:
        img_info = histology.strip().split(";")
        assert len(img_info) >= 2, ("Error: --histology should at least have <image_id> and <histology_path>." if not by_tiles else "Error: --histology-tiles should have <image_id>, and at least one <tile_info>")
        img_dict = {"image_id": img_info[0]}
        if not by_tiles:
            values = img_info[1:]
            for i, key in enumerate(img_keys):
                val = values[i] if i < len(values) else None
                if val in ["", "null", "Null"]:
                    val = None
                img_dict[key] = val          
        else:
            img_dict["in_tiles"] = []
            for tile in img_info[1:]:
                tile_fields = tile.strip().split(",")
                assert len(tile_fields) >= 3, "Error: each <tile_info> must have at least <row>,<col>,<path>"
                tile_dict = {
                    "row": tile_fields[0],
                    "col": tile_fields[1]
                }
                tile_values = tile_fields[2:]
                for i, key in enumerate(img_keys_tiles):
                    val = tile_values[i] if i < len(tile_values) else None
                    if val in ["", "null", "Null"]:
                        val = None
                    tile_dict[key] = val
                img_dict["in_tiles"].append(tile_dict)
        img_input.append(img_dict)
    return img_input

def update_flip_to_vhflip(val):
    vertical = "True" if str(val).lower() == "vertical" else "False"
    horizontal = "True" if str(val).lower() == "horizontal" else "False"
    return vertical, horizontal

def generate_in_tiles_str(in_tiles, param_keys):
    assert len(in_tiles) > 0, "Error: --in-tiles should have at least one tile"
    tile_args = []
    for tile in in_tiles:
        # row = tile.get("row")
        # col = tile.get("col")
        # values =[str(row), str(col)]
        values = []
        for key in param_keys:
            val = tile.get(key, "")
            if key == "georef_bounds":
                val=val.replace("(","")
                val=val.replace(")","")
                val=val.replace(",","_")
            #print(key+":"+str(val))
            if val in [None, "null", "Null"]:
                val = ""
            if key == "flip":
                vertical, horizontal = update_flip_to_vhflip(val)
                values.append(vertical)
                values.append(horizontal)
            elif key == "row" or key == "col":
                values.append(str(val))
            else:
                values.append(str(val))
        tile_args.append(",".join(values))
        #print(tile_args)
    in_tiles_str = " ".join(tile_args)
    #print(in_tiles_str)
    return in_tiles_str

def cmd_image_stitch(imginfo_by_tiles, args, env):
    image_stitch_cmds = []        
    for imginfo_i in imginfo_by_tiles:
        image_id= imginfo_i["image_id"]
        mkbn = f"image_stitch_{image_id}" if args.mk_id is None else f"image_stitch_{image_id}_{args.mk_id}"
        in_tiles_str = generate_in_tiles_str(imginfo_i.get("in_tiles",[]), img_keys_tiles)
        image_stitch_cmd=" ".join([
                "cartloader", "image_stitch",
                f"--output {imginfo_i['path']}", 
                f"--makefn {mkbn}.mk",
                f"--in-offsets {imginfo_i['in_offsets']}",
                f"--in-tiles {in_tiles_str}",
                f"--crop-tile-by-minmax",
                f"--in-minmax {imginfo_i['in_offsets']}",
                f"--restart" if args.restart else "",
                f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        ])
        # add aux tools
        image_stitch_cmd = add_param_to_cmd(image_stitch_cmd, env, aux_env_args["image_stitch"], underscore2dash=False)
        image_stitch_cmds.append(image_stitch_cmd)
    return image_stitch_cmds
  
def cmd_image_png2pmtiles(run_i, cartl_v, args, env):    
    assert len(run_i.get("histology", [])) > 0, "Error: --histology is Required when running --image-png2pmtiles"
    img_cmds=[]
    #cartload_dir=os.path.join(run_i["run_dir"], "cartload")
    cartload_dir = os.path.join(run_i["run_dir"], ("cartload" if cartl_v == "1" else "cartload2"))
    catalog_yaml=os.path.join(cartload_dir, "catalog.yaml")
    for histology in run_i.get("histology", []):
        # 1. histology tif to pmtiles
        img_path = histology["path"]
        # prefix
        img_inname = os.path.basename(img_path).replace("_", "-")
        for suffix in img_suffix:
            if img_inname.endswith(suffix):
                img_inname = img_inname[:-len(suffix)]
                break
        img_prefix = os.path.join(cartload_dir, img_inname)
        
        # update the orientation
        histology = update_orient_in_histology(histology)

        img_cmd = " ".join([
            "cartloader", "import_image",
            # actions
            "--ome2png" if histology.get("ome", False) else "",
            "--png2pmtiles",
            "--georeference" if histology.get("georeference", False) else "",
            "--flip-vertical" if histology["flip"] in ["vertical", "both"] else "",
            "--flip-horizontal" if histology["flip"] in ["horizontal", "both"] else "",
            f"--rotate {histology['rotate']}" if histology.get("rotate", None) is not None else "",
            f"--update-catalog",
            # in/out
            f"--in-img {img_path}",
            f"--out-dir {cartload_dir}",
            # id and color
            f"--img-id {histology['image_id']}",
            f"--colorize {histology.get('colorize', None)}" if histology.get("ome", False) and histology.get("colorize", None) is not None else "",
            # aux params
            f"--upper-thres-quantile {histology.get('upper_thres_quantile', None)}" if histology.get("ome", False) and histology.get("upper_thres_quantile", None) is not None else "",
            f"--lower-thres-quantile {histology.get('lower_thres_quantile', None)}" if histology.get("ome", False) and histology.get("lower_thres_quantile", None) is not None else "",
            f"--upper-thres-intensity {histology.get('upper_thres_intensity', None)}" if histology.get("ome", False) and histology.get("upper_thres_intensity", None) is not None else "",
            f"--lower-thres-intensity {histology.get('lower_thres_intensity', None)}" if histology.get("ome", False) and histology.get("lower_thres_intensity", None) is not None else "",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
            f"--georef-pixel-tsv {histology.get('georef_tsv', None)}" if histology.get("georeference", False) and histology.get("georef_tsv", None) is not None else "",
            f"--georef-bounds {histology.get('georef_bounds', None)}" if histology.get("georeference", False) and histology.get("georef_bounds", None) is not None else "",
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
        ])
        img_cmd = add_param_to_cmd(img_cmd, env, aux_env_args["image_png2pmtiles"], underscore2dash=False)
        img_cmds.append(img_cmd)

    return img_cmds

def cmd_upload_aws(run_i, cartl_v,  args, env, additional_args=None):
    cartload_dir=os.path.join(run_i["run_dir"], "cartload") if cartl_v == "1" else os.path.join(run_i["run_dir"], "cartload2")
    # Option 1: use "." to locate the file and handle by perl -lane
    # aws_cmd="\n".join([
    #     f"out={cartload_dir}",
    #     f"id={run_i['run_id']}",
    #     f"aws_bucket={args.aws_bucket}",
    #     "aws s3 cp ${out}/catalog.yaml s3://${aws_bucket}/${id}/catalog.yaml",
    #     "grep -E '\.' ${out}/catalog.yaml | perl -lane 'print $F[$#F]' | xargs -I {} aws s3 cp ${out}/{} s3://${aws_bucket}/${id}/{}"
    # ])
    # Option 2: use cartloader upload_aws_by_catalog.py
    aws_cmd=" ".join([
        "cartloader", "upload_aws",
        f"--in-dir {cartload_dir}",
        f"--s3-dir \"s3://{args.aws_bucket}/{run_i['run_id']}\"",
        additional_args if additional_args is not None else "",
        f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        f"--restart" if args.restart else ""
        
    ])
    aws_cmd=add_param_to_cmd(aws_cmd, env, aux_env_args["upload_aws"])

    return aws_cmd

def stepinator(_args):
    parser = argparse.ArgumentParser(description="""
    Stepinator: A tool to run steps in cartloader to do SGE format conversion (--sge-convert) or downstream process (--run-ficture*, --run-cartload*, --image-png2pmtiles, --upload-aws) in cartloader.
    The input, actions, parameters, and tools can be provided in two ways: 1) using a YAML file or 2) using individual command-line arguments in IN/OUT Configuration, Environment Configuration, and Auxiliary Parameters for FICTURE.
    It also allows the job execution in two modes: local and slurm.
    """)

    # * mode
    run_params = parser.add_argument_group("Run", "Run mode")
    run_params.add_argument("--dry-run", action="store_true", help="Perform a dry run")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', '-j', type=int, default=None, help='Number of processes to run in parallel (default: 1)')
    run_params.add_argument('--threads', type=int, default=10, help='Maximum number of threads per job (for tippecanoe)')
    run_params.add_argument("--submit", action="store_true", help="Submit the job")
    run_params.add_argument('--submit-mode', type=str, default="local", choices=["slurm", "local"], help='Specify how the job should be executed. Choose "slurm" for SLURM, "local" for local execution')
    run_params.add_argument('--job-id', type=str, default=None, help='Job ID, which will be used as the base name for the job file if --job-file is not provided and the job name for SLURM if --submit-mode is slurm' )
    run_params.add_argument('--mk-id', type=str, default=None, help='Makefile identifier. This helps when multiple sets of parameters will be applied to one function, such as --run_ficture.')
    # run_params.add_argument('--debug', action='store_true', help='For the debug mode')

    # * commands 
    cmd_params = parser.add_argument_group("Commands", "Commands to be performed")
    cmd_params.add_argument("--sge-stitch", action="store_true", help="Run sge-stitch in cartloader")
    cmd_params.add_argument("--list-overlapping-features", action="store_true", help="(Optional if --sge-stitch) List overlapping features in sge-stitch")
    cmd_params.add_argument("--sge-convert", action="store_true", help="Run sge-convert in cartloader")
    cmd_params.add_argument('--sge-visual', action='store_true', help='Plot the SGE from sge-convert or sge-stitch in a PNG file')
    cmd_params.add_argument("--north-up", action="store_true", help="(Optional if --sge-visual) Set the north direction to up in the SGE visual")
    # cmd_params.add_argument("--run-ficture1", action="store_true", help="Run run-ficture in cartloader. Only the main function is executed.")
    # cmd_params.add_argument("--run-cartload1", action="store_true", help="Run run-cartload-join in cartloader")
    cmd_params.add_argument("--run-ficture", action="store_true", help="Run run-ficture in cartloader.")
    cmd_params.add_argument("--run-cartload", action="store_true", help="Run run-cartload in cartloader.")
    cmd_params.add_argument("--image-stitch", action="store_true", help="Run image stitch in cartloader. This is for stitching images.")
    cmd_params.add_argument("--image-png2pmtiles", action="store_true", help="Run image-png2pmtiles in cartloader. This s provide image using --in-yaml or --histology.")
    cmd_params.add_argument("--upload-aws", action="store_true", help="Upload files to AWS S3 bucket. This s provide AWS bucket name using --in-yaml or --aws-bucket.")
    cmd_params.add_argument("--copy-ext-model", action="store_true", help="(Optional if --run-ficture1) Auxiliary action parameters for run-ficture. Copy external model when running FICTURE with an external model")
    cmd_params.add_argument("--init-ext", action="store_true", help="(Optional if --run-ficture1) Auxiliary action parameters for run-ficture. Only Initialize external model without run main-ext")
    cmd_params.add_argument("--skip-coarse-report", action="store_true", help="(Optional if --run-ficture1)Auxiliary action parameters for run-ficture. Skip coarse report")
    cmd_params.add_argument('--segment-10x', action='store_true', help='(Optional if --run-ficture1) Perform hexagon segmentation into 10x Genomics format')
    cmd_params.add_argument("--ficture-version", type=str, default=None, help="Define which ficture version you want to use.", choices=["1","2"])

    # * Key
    key_params = parser.add_argument_group(
    "IN/OUT Configuration",
    "Specify input data and settings using one of the following methods:\n"
    "1. Use --in-yaml (and --run-ids for downstream processing) to provide inputs via a YAML file, enabling batch processing of multiple runs.\n"
    "2. Use individual command-line arguments to specify input files and settings, allowing processing of a single run at a time."
    )
    key_params.add_argument("--in-yaml", '-i', type=str, default=None, help="Input yaml file")
    key_params.add_argument("--run-ids", nargs="*", default=[], help="Run IDs. Required if --run-ficture, --run-cartload-join, --image-png2pmtiles, or --upload-aws is applied. ")
    key_params.add_argument("--image-ids", nargs="*", default=[], help="Image IDs. Specifies the image IDs in the --in-yaml for the --image-png2pmtiles command. If omitted, all images in the input YAML are processed.")
    key_params.add_argument("--out-dir", type=str, default=None, help="Output directory. It will have 3 subdirectories: sge (output from --sge-convert), ficture (output from --run-ficture), and cartload (output from --run-cartload-join and --image-png2pmtiles)")
    # for sge_convert
    key_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st"], help='Required if --sge-convert. Platform of the raw input file to infer the format of the input file. ')
    # - input for 10x_visium_hd, seqscope 
    key_params.add_argument('--in-mex', type=str, default=None, help='Required if --sge-convert on 10x_visium_hd and seqscope datasets. Directory path to input files in Market Exchange (MEX) format.')
    key_params.add_argument('--in-parquet', type=str, default=None, help='Required if --sge-convert on 10x_visium_hd datasets. Path to the input raw parquet file for spatial coordinates (default: None, typical naming convention: tissue_positions.parquet)')
    # - input for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st
    key_params.add_argument('--in-csv', type=str, default=None, help='Required if --sge-convert for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st. Path to the input raw CSV/TSV file (default: None).')
    key_params.add_argument("--in-tiles", nargs="*", default=[], help="Required if --sge-stitch. List of input tiles of each in the format of <transcript_path>,<feature_path>,<minmax_path>,<row>,<col>,<rotate>,<flip>. If feature or minmax is missing, specify its path as null (default: [])")
    key_params.add_argument('--units-per-um', type=float, default=None, help='Applicable if --sge-convert or --sge-stitch. Coordinate unit per um (conversion factor) (default: 1.00)') 
    key_params.add_argument('--scale-json', type=str, default=None, help="Applicable if --sge-convert on 10x_visium_hd datasets. Coordinate unit per um using the scale json file (default: None, typical naming convention: scalefactors_json.json)")
    key_params.add_argument('--precision-um', type=int, default=None, help='Required if --sge-convert. Number of digits of transcript coordinates (default: 2)')
    key_params.add_argument('--filter-by-density', action='store_true', default=False, help='Required if --sge-convert or --run-ficture. If --sge-convert, it enables density-filtering. If --run-ficture, it defines the density-filtered SGE as input. (default: False)')
    key_params.add_argument('--filtered-prefix', type=str, default=None, help='Required if --filtered-by-density. The prefix for density-filtered SGE (default: filtered)')
    key_params.add_argument("--colname-count", type=str, default="count", help="Required if --sge-convert, --sge-stitch, --run-ficture, or --run-cartload-join. Column name that showing the expression count of the genomic feature of interest. (default: gene)")
    # key_params.add_argument("--colnames-other-count", nargs="*", default=[], help="Optional if --sge-convert, --sge-stitch is enabled. It allows to keep other genomic features in the formatted SGE besides the genomic feature of interest. (default: [])")
    # for run_ficture
    key_params.add_argument("--major-axis", type=str, default="X", choices=["X","Y"], help="Major axis (default: X)")
    key_params.add_argument("--ext-path", type=str, default=None, help="Required when --run-ficture1 with an external model. The path for the external model.")
    key_params.add_argument("--ext-id", type=str, default=None, help="Required when --run-ficture1 with an external model. The ID for the external model.")
    key_params.add_argument("--train-width", '-w', type=str, default=None, help="Required if --run-ficture. Train width.")
    key_params.add_argument("--n-factor", '-n', type=str, default=None, help="Required if --run-ficture1 with LDA. Number of factors. ")
    key_params.add_argument("--image", type=str, nargs="?", default=[], help="""
                              (Optional) Provide image info as <image_id>;<path>;<ome>;<lower_thres_quantile>;<upper_thres_quantile>;<level>;<colorize>;<georeference>;<georef_tsv>;<georef_bounds>;<rotate_degree>;<flip_direction>. 
                              Only <image_id> and <path> are required. Supports multiple files; use <image_id> to distinguish them. (Default: [])
                              Define <ome> and <georeference> by True or False. If the image is OME TIFF, define <ome> as True. If georeference=True, specify <georef_tsv> or <georef_bounds>. 
                              Orientation allows rotation (90, 180, 270) and flipping (vertical, horizontal, both). 
                            """)
    key_params.add_argument("--image-tiles", type=str, default=[], nargs="*", help="""
                            (Optional) Provide the image tiles infor for image stitch. The format should be <image_id>;<tile_info1>;<tile_info2>... where each tile_info is <row>,<col>,<path>,<georeference>,<georef_tsv>,<georef_bounds>,<rotate_degree>,<flip_direction>.
                            Note use ";" to separate the image_id and tile_info, and use "," to separate the information in the tile_info. Currently, our image stitch doesn't support transform ome tiles given it's not finalized 
                            """) 
    key_params.add_argument('--aws-bucket', type=str, default=None, help='Required if --upload-aws. AWS bucket name')

    # * tools
    env_params = parser.add_argument_group("Environment Configuration", 
    "Specify paths to tools and environment settings using one of the following methods:\n" 
    "1. Provide all environment settings in the input YAML file (--in-yaml).\n" 
    "2. Provide individual tool paths and environment variables directly as command-line arguments.\n"  
    )
    env_params.add_argument('--slurm-account', type=str, default=None, help='If "--submit-mode slurm", provide a SLURM account')
    env_params.add_argument('--slurm-partition', type=str, default=None, help='If "--submit-mode slurm", provide a SLURM partition')
    env_params.add_argument('--slurm-mail-user', type=str, default=None, help='If "--submit-mode slurm", provide a SLURM mail user')
    env_params.add_argument('--slurm-cpus-per-task', type=int, default=10, help='If "--submit-mode slurm", provide a SLURM cpus per task')
    env_params.add_argument('--slurm-mem-per-cpu', type=str, default="6500mb", help='If "--submit-mode slurm", provide a SLURM memory per cpu')
    # modules & env
    env_params.add_argument('--hpc-modules', type=str, nargs="*", default=[], help='Comma-separated HPC modules to load. When a version is required, use the format: <module>/<version>')
    env_params.add_argument('--conda-base', type=str, default=None, help='Conda base path')
    env_params.add_argument('--conda-env', type=str, default=None, help='Conda environment to activate. If it is installed in the default location, i.e., within the base path, provide the environment name. Otherwise, provide the full path to the environment.')
    # app
    env_params.add_argument('--gzip', type=str, default=None, help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--spatula', type=str, default=None, help='Path to spatula binary.')    
    env_params.add_argument('--parquet-tools', type=str, default=None, help='Path to parquet-tools binary')
    env_params.add_argument('--sort', type=str, default=None, help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default=None, help='Memory size for each process')
    env_params.add_argument('--bgzip', type=str, default=None, help='Path to bgzip binary. For faster processing, use "bgzip -@ 4')
    env_params.add_argument('--tabix', type=str, default=None, help='Path to tabix binary')
    env_params.add_argument('--pmtiles', type=str, default=None, help='Path to pmtiles binary.')
    env_params.add_argument('--tippecanoe', type=str, default=None, help='Path to tippecanoe binary. ')
    env_params.add_argument('--imagemagick', type=str, default=None, help='Path to the imagemagick binary directory')
    env_params.add_argument('--gdal_translate', type=str, default=None, help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=None, help='Path to gdaladdo binary')
    env_params.add_argument('--aws', type=str, default=None, help='Path to aws binary')
    env_params.add_argument('--ficture2', type=str, help='Path to punkst(ficture2) repository') 
    env_params.add_argument('--python', type=str,  help='Python3 binary')

    format_aux_parameter = parser.add_argument_group(
        "Auxiliary Parameters for sge_convert", 
        "Parameters for sge_convert. Required if --sge-convert is used with non-default values."
    )
    # IN-MEX
    format_aux_parameter.add_argument('--icols-mtx', type=str, default=None, help='Input column indices (comma-separated 1-based) in the input matrix file (default: 1,2,3,4,5)')
    format_aux_parameter.add_argument('--icol-bcd-barcode', type=int, default=None, help='(seqscope only) 1-based column index of barcode in the input barcode file (default: 1)')
    format_aux_parameter.add_argument('--icol-bcd-x', type=int, default=None, help='(seqscope only) 1-based column index of x coordinate in the input barcode file (default: 6)')
    format_aux_parameter.add_argument('--icol-bcd-y', type=int, default=None, help='(seqscope only) 1-based column index of y coordinate in the input barcode file (default: 7)')
    format_aux_parameter.add_argument('--icol-ftr-id', type=int, default=None, help='1-based column index of feature ID in the input feature file (default: 1)')
    format_aux_parameter.add_argument('--icol-ftr-name', type=int, default=None, help='1-based column index of feature name in the input feature file (default: 2)')
    format_aux_parameter.add_argument('--pos-colname-barcode', type=str, default=None, help='(10x_visium_hd only) Column name for barcode in the input parquet file (default: barcode)')
    format_aux_parameter.add_argument('--pos-colname-x', type=str, default=None, help='(10x_visium_hd only) Column name for X-axis in the input parquet file (default: pxl_row_in_fullres)')
    format_aux_parameter.add_argument('--pos-colname-y', type=str, default=None, help='(10x_visium_hd only) Column name for Y-axis in the input parquet file (default: pxl_col_in_fullres)')
    format_aux_parameter.add_argument('--pos-delim', type=str, default=None, help='(10x_visium_hd only) Delimiter for the input parquet file (default: ",")')
    # IN-CSV
    format_aux_parameter.add_argument('--csv-delim', type=str, default=None, help='Delimiter for the additional input tsv/csv file (default: "," for 10x_xenium, cosmx_smi, and vizgen_merscope; "\\t" for bgi_stereoseq, pixel_seq, and nova_st) ')
    format_aux_parameter.add_argument('--csv-colname-x', type=str, default=None, help='Column name for X-axis (default: x_location for 10x_xenium; x for bgi_stereoseq; x_local_px for cosmx_smi; global_x for vizgen_merscope; xcoord for pixel_seq; x for nova_st)')
    format_aux_parameter.add_argument('--csv-colname-y', type=str, default=None, help='Column name for Y-axis (default: y_location for 10x_xenium; y for bgi_stereoseq; y_local_px for cosmx_smi; global_y for vizgen_merscope; ycoord for pixel_seq; y for nova_st)')
    format_aux_parameter.add_argument('--csv-colname-feature-name', type=str, default=None, help='Column name for gene name (default: feature_name for 10x_xenium; geneID for bgi_stereoseq; target for cosmx_smi; gene for vizgen_merscope; geneName for pixel_seq; geneID for nova_st)')
    format_aux_parameter.add_argument('--csv-colname-count', type=str, default=None, help='Column name for expression count. If not provided, a count of 1 will be added for a feature in a pixel (default: MIDCounts for bgi_stereoseq; MIDCount for nova_st; None for the rest platforms).')
    format_aux_parameter.add_argument('--csv-colname-feature-id', type=str, default=None, help='Column name for gene id (default: None)')
    format_aux_parameter.add_argument('--csv-colnames-others', nargs='+', default=[], help='Columns names to keep (e.g., cell_id, overlaps_nucleus) (default: None)')
    format_aux_parameter.add_argument('--csv-colname-phredscore', type=str, default=None, help='Column name for Phred-scaled quality value (Q-Score) estimating the probability of incorrect call (default: qv for 10x_xenium and None for the rest platforms).') # qv
    format_aux_parameter.add_argument('--min-phred-score', type=float, default=None, help='Phred-scaled quality score cutoff (default: 20 for 10x_xenium and None for the rest platforms).')
    # feature-filtering
    format_aux_parameter.add_argument('--include-feature-list', type=str, default=None, help='A file provides a list of gene names to include (default: None)')
    format_aux_parameter.add_argument('--exclude-feature-list', type=str, default=None, help='A file provides a list of gene names to exclude (default: None)')
    format_aux_parameter.add_argument('--include-feature-regex', type=str, default=None, help='Regex pattern for gene names to include (default: None)')
    format_aux_parameter.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex pattern for gene names to exclude (default: "^(BLANK|Blank-|NegCon|NegPrb)"). To skip this, use --exclude-feature-regex "".')
    format_aux_parameter.add_argument('--include-feature-type-regex', type=str, default=None, help='Regex pattern for gene type to include (default: None). Requires --csv-colname-feature-type or --feature-type-ref for gene type info') # (e.g. protein_coding|lncRNA)
    format_aux_parameter.add_argument('--csv-colname-feature-type', type=str, default=None, help='If --include-feature-type-regex is used and the input file has gene type, define the column name for gene type info (default: None)')
    format_aux_parameter.add_argument('--feature-type-ref', type=str, default=None, help='Path to a tab-separated gene reference file containing gene type information. The format should be: chrom, start position, end position, gene id, gene name, gene type (default: None)')
    # Polygon-filtering
    #format_aux_parameter.add_argument('--mu-scale', type=int, default=None, help='Scale factor for the polygon area calculation (default: 1.0)')   # Always use --mu-scale 1, since this will be performed after sge conversion, which will apply --units-per-um to generate the output SGE in um
    format_aux_parameter.add_argument('--radius', type=int, default=None, help='Radius for the polygon area calculation (default: 15)')
    format_aux_parameter.add_argument('--quartile', type=int, default=None, help='Quartile for the polygon area calculation (default: 2)')
    format_aux_parameter.add_argument('--hex-n-move', type=int, default=None, help='Sliding step (default: 1)')
    format_aux_parameter.add_argument('--polygon-min-size', type=int, default=None, help='The minimum polygon size (default: 500)')

    stitch_aux_parameter = parser.add_argument_group(
        "Auxiliary Parameters for sge_stitch", 
        "Parameters for sge_stitch. Required if --sge-stitch is used with non-default values."
    )
    stitch_aux_parameter.add_argument('--colname-feature-name', type=str, default=None, help='Feature name column (default: gene)')
    stitch_aux_parameter.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    stitch_aux_parameter.add_argument('--colname-x', type=str, default=None, help='X column name (default: X)')
    stitch_aux_parameter.add_argument('--colname-y', type=str, default=None, help='Y column name (default: Y)')

    ficture_aux_params = parser.add_argument_group(
        "Auxiliary Parameters for FICTURE", 
        "Parameters for run_ficture. Required if --run-ficture is used with non-default values. Default values are recommended."
    )
    # use fic
    ficture_aux_params.add_argument('--filter-by-overlapping-features', action='store_true', default=False, help='Use overlapping features in FICTURE analysis (default: False)')
    ficture_aux_params.add_argument('--min-ct-per-ftr-tile', type=int, default=0, help='Apply a minimum count to filter overlapping feature. Filtering process will be applied if --min-ct-per-overlapftr > 0. (default: 0)')
    ficture_aux_params.add_argument("--cmap", '-c', type=str, default=None, help="Required if the user prefers to use a pre-built color map (Default: None)")
    ficture_aux_params.add_argument('--ficture-feature', type=str, default="features.ficture.tsv.gz", help='File name for the output TSV file of feature used in FICTURE analysis (default: None)')
    ficture_aux_params.add_argument('--fic-include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    ficture_aux_params.add_argument('--fic-exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    ficture_aux_params.add_argument('--fic-include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    ficture_aux_params.add_argument('--fic-exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')
    ficture_aux_params.add_argument('--fic-include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). Requires --csv-colname-feature-type or --feature-type-ref for gene type info') # (e.g. protein_coding|lncRNA)
    # segmentation - ficture
    ficture_aux_params.add_argument('--hexagon-n-move', type=int, default=None, help='Level of hexagonal sliding when creating hexagon-indexed SGE in FICTURE compatible format (default: 1)')
    ficture_aux_params.add_argument('--hexagon-precision', type=float, default=None, help='Output precision of hexagon coordinates for FICTURE compatible format (default: 2)')
    ficture_aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=None, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format (default: 50)')
    # segmentation - 10x
    ficture_aux_params.add_argument('--hexagon-n-move-10x', type=int, default=None, help='Level of hexagonal sliding when creating hexagon-indexed SGE in 10x Genomics format (default: 1)')
    ficture_aux_params.add_argument('--hexagon-precision-10x', type=float, default=None, help='Output precision of hexagon coordinates for 10x Genomics format (default: 2)')
    ficture_aux_params.add_argument('--min-ct-per-unit-hexagon-10x', type=int, default=None, help='Minimum count per hexagon in hexagon segmentation in 10x Genomics format (default: 1)')
    # minibatch
    ficture_aux_params.add_argument('--minibatch-size', type=int, default=None, help='Batch size used in minibatch processing (default: 500)')
    ficture_aux_params.add_argument('--minibatch-buffer', type=int, default=None, help='Batch buffer used in minibatch processing (default: 30)')
    # train 
    ficture_aux_params.add_argument('--train-epoch', type=int, default=None, help='Training epoch for LDA model (default: 3)')
    ficture_aux_params.add_argument('--train-epoch-id-len', type=int, default=None, help='Training epoch ID length (default: 2)')
    ficture_aux_params.add_argument('--lda-rand-init', type=int, default=None, help='Number of random initialization during model training (default: 10)')
    ficture_aux_params.add_argument('--lda-plot-um-per-pixel', type=float, default=None, help='Image resolution for LDA plot (default: 1)')
    # fit 
    ficture_aux_params.add_argument('--fit-width', type=str, default=None, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    ficture_aux_params.add_argument('--fit-precision', type=float, default=None, help='Output precision of model fitting (default: 2)')
    ficture_aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=None, help='Minimum count per hexagon unit during model fitting (default: 20)')
    ficture_aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=None, help='Image resolution for fit coarse plot (default: 1)')   # in Scopeflow, this is set to 2
    # decode
    ficture_aux_params.add_argument('--anchor-res', type=int, default=None, help='Anchor resolution for decoding (default: 4)')
    ficture_aux_params.add_argument('--radius-buffer', type=int, default=None, help='Buffer to radius(=anchor_res + radius_buffer) for pixel-level decoding (default: 1)')
    ficture_aux_params.add_argument('--decode-top-k', type=int, default=None, help='Top K columns to output in pixel-level decoding results (default: 3)')
    ficture_aux_params.add_argument('--decode-block-size', type=int, default=None, help='Block size for pixel decoding output (default: 100)')
    ficture_aux_params.add_argument('--decode-scale', type=int, default=None, help='Scale parameters for pixel decoding output (default: 100)')
    ficture_aux_params.add_argument('--decode-precision', type=float, default=None, help='Precision of pixel level decoding (default: 0.01)')
    ficture_aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=None, help='Image resolution for pixel decoding plot (default: 0.5)')
    # merge_by_pixel
    ficture_aux_params.add_argument('--merge-max-dist-um', type=float, default=None, help='Maximum distance in um for merging pixel-level decoding results (default: 0.1)') 
    ficture_aux_params.add_argument('--merge-max-k', type=int, default=None, help='Maximum number of K columns to output in merged pixel-level decoding results (default: 1)')
    ficture_aux_params.add_argument('--merge-max-p', type=int, default=None, help='Maximum number of P columns to output in merged pixel-level decoding results (default: 1)')
    # others parameters shared across steps
    ficture_aux_params.add_argument('--min-ct-per-feature', type=int, default=None, help='Minimum count per feature during LDA training, projection and decoding (default: 20)')
    ficture_aux_params.add_argument('--de-max-pval', type=float, default=None, help='p-value cutoff for differential expression (default: 1e-3)')
    ficture_aux_params.add_argument('--de-min-fold', type=float, default=None, help='Fold-change cutoff for differential expression (default: 1.5)')

    gis_aux_params = parser.add_argument_group("Auxiliary Parameters for images", "Parameters for images")
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help=' define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--resample', type=str, default='cubic', help='define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')

    dev_params = parser.add_argument_group("Developer", "Developer options")
    dev_params.add_argument('--lenient', action='store_true', help='(Only for development usage) Lenient mode. Skip the required file/app existence check')

    args = parser.parse_args(_args)

    # =========
    #  actions
    # =========
    assert args.sge_stitch or args.sge_convert or args.run_ficture or args.run_cartload or args.image_png2pmtiles or args.upload_aws, "Error: At least one action is required"

    # sge
    assert not (args.sge_convert and args.sge_stitch), "Error: --sge-convert and --sge-stitch cannot be applied together"
    assert not (args.sge_visual and not args.sge_convert and not args.sge_stitch), "Error: --sge-visual can only be applied with --sge-convert or --sge-stitch"
    assert not (args.north_up and not args.sge_visual), "Error: --north-up can only be applied with --sge-visual"   

    # version check
    if args.run_ficture or args.run_cartload or args.image_png2pmtiles or args.upload_aws:
        assert args.ficture_version is not None, "Error: --ficture-version must be specified if any of the following actions is enabled: --run-ficture, --run-cartload, --image-png2pmtiles or --upload-aws"
    
    fic_v = str(args.ficture_version)

    # =========
    #  Read YAML/args
    # =========
    if args.in_yaml is not None:
        with open(args.in_yaml, "r") as f:
            yml = yaml.safe_load(f)

    # output dir
    out_dir = args.out_dir if args.out_dir else yml.get("out_dir", None)
    assert out_dir is not None, "Error: --out-dir is required"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # sge infomation - shared across actions
    if args.in_yaml:
        sgeinfo = yml.get("SGE", {})
        # do density-filtering by default to offer the raw and filtered SGE & xy.png for users to decide.
    else:
        sgeinfo={
            "platform": args.platform,
            "in_mex": args.in_mex,
            "in_parquet": args.in_parquet,
            "in_csv": args.in_csv,
            "units_per_um": args.units_per_um,
            "scale_json": args.scale_json,
            "precision_um": args.precision_um,
            "filter_by_density": args.filter_by_density,
            "filtered_prefix": args.filtered_prefix,
            "colname_count": args.colname_count             # output column name for count
        }
    # add default values for some parameters
    # filter_by_density
    sgeinfo["filter_by_density"] = sgeinfo.get("filter_by_density", False)
    sgeinfo["filtered_prefix"] = sgeinfo.get("filtered_prefix", "filtered")
    # colnames for output
    sgeinfo["colname_count"] = sgeinfo.get("colname_count", "count")
    #sgeinfo["colnames_all_count"] = ",".join([sgeinfo["colname_count"]] + sgeinfo["colnames_other_count"].split(",")) if sgeinfo.get("colnames_other_count", None) is not None else sgeinfo["colname_count"]
    # sge_dir    
    sgeinfo["sge_dir"] = os.path.join(out_dir, "sge")
    os.makedirs(sgeinfo['sge_dir'], exist_ok=True)
    sgeinfo["xy_png"] = "xy.png"
    sgeinfo["xy_ntif"] = "xy_northup.tif"
    # worklog dir
    log_dir = os.path.join(out_dir, "worklog")
    os.makedirs(log_dir, exist_ok=True)

    # =========
    #  env
    # =========
    # * tools (collect from yaml and args) 
    env = merge_config(yml, args, 
                       [item for sublist in aux_env_args.values() for item in sublist] + ["imagemagick"],
                        prefix="env")  

    # * slurm (collect from yaml and args)
    slurm = merge_config(yml, args, ["account", "partition", "mail_user", "cpus_per_task", "mem_per_cpu", "hours"], prefix="slurm")
    
    # * hpc modules
    hpc_modules =  yml.get("env", {}).get("hpc_modules", [])
    if args.hpc_modules:
        hpc_modules = args.hpc_modules.split(",")

    if len(hpc_modules) > 0:
        module_load_cmds="module load "+ " ".join(hpc_modules)
    else:
        module_load_cmds=""

    # * conda env
    conda_base = args.conda_base if args.conda_base else yml.get("env", {}).get("conda", {}).get("base_path", None)
    conda_env  = args.conda_env  if args.conda_env  else yml.get("env", {}).get("conda", {}).get("env_prefix", None)

    conda_cmd=[]
    if conda_base:
        conda_cmd.append(f"source {conda_base}/etc/profile.d/conda.sh")
    if conda_env:
        conda_cmd.append(f"conda activate {conda_env}")

    # * aws 
    if args.aws_bucket is None:
        args.aws_bucket = yml.get("aws_bucket", None)

    # =========
    #  init 
    # =========
    env_cmds=[module_load_cmds] + conda_cmd
    # export path for imagemagick
    if args.sge_stitch or args.run_cartload:
        if getattr(env, "imagemagick", None) is not None:
            env_cmds.append("export PATH=$PATH:"+":".join([env.imagemagick]))
    env_cmds.append("\n")
    cmds = []

    # =========
    #  SGE Convert / Stitch
    # =========
    if args.sge_convert:
        cmds.append(cmd_sge_convert(sgeinfo, args, env))

    if args.sge_stitch:
        # comma-based input per tile
        if args.in_yaml:
            in_tiles=sgeinfo.get("in_tiles", [])
            assert len(in_tiles) > 1, 'When --sge-stitch is enabled, at least two tiles are required.'
            sge_keys_tile = ["in_csv", "in_feature", "in_minmax", "row", "col", "rotate", "flip"]
            sgeinfo["in_tiles_str"] = generate_in_tiles_str(in_tiles, sge_keys_tile)
        else:
            assert len(args.in_tiles) > 1, "When --sge-stitch is enabled, at least two tiles are required"
            #sgeinfo["in_tiles_str"]=args.in_tiles
            in_tiles=[]
            for i, intile in enumerate(args.in_tiles):
                intile_parts = intile.split(",")
                assert len(intile_parts) == 7, f"Error: --in-tiles should be in the format <in_csv>,<in_feature>,<in_minmax>,<row>,<col>,<rotate>,<flip>. Found {len(intile_parts)} parts instead."
                # get the last piece
                intile_flip=intile_parts[-1]
                vertical, horizontal = update_flip_to_vhflip(intile_flip)
                intile_parts[-1] = vertical + "," + horizontal
                intile=",".join(intile_parts)
                in_tiles.append(intile)
            sgeinfo["in_tiles_str"]=" ".join(in_tiles)
        cmds.append(cmd_sge_stitch(sgeinfo, args, env))

    # =========
    #  IMAGE Stitch
    # =========
    # the stitched tif should be shared across run_id in the sge
    if args.image_png2pmtiles or args.image_stitch:
        if args.in_yaml:
            imginfo = yml.get("HISTOLOGY", [])
            if len(args.image_ids) > 0:
                imginfo = [img_i for img_i in imginfo if img_i.get("image_id") in args.image_ids]
        else:
            imginfo_by_tif = parse_image_args(args.image) if len(args.image) > 0 else []
            imginfo_by_tiles = parse_image_args(args.image_tiles, by_tiles=True) if len(args.image_tiles) > 0 else []
            imginfo = imginfo_by_tif + imginfo_by_tiles
        # sanity check: no duplicate image_id in imginfo
        image_ids = [img_i.get("image_id") for img_i in imginfo]
        assert len(image_ids) == len(set(image_ids)), "Error: Duplicate Image IDs found in the input YAML file or command line arguments. Please ensure that each image ID is unique."

    # image stitch
    if args.image_stitch:
        # img_dir to host the output tif from img stitch
        img_dir = os.path.join(out_dir, "hist")
        os.makedirs(img_dir, exist_ok=True)
        # add a cmd to make sure tile_offsets file exists - it was not generated for the early jobs
        tile_offsets_f = os.path.join(sgeinfo["sge_dir"], "coordinate_minmax_per_tile.tsv")
        if not os.path.exists(tile_offsets_f) and not args.sge_stitch:
            cmd.append(cmd_sge_stitch(sgeinfo, args, env, generate_tile_minmax_only=True))
        # extract the img info for img stitch, and add the path
        imginfo_by_tiles = []
        for img_i in imginfo:
            if len(img_i.get("in_tiles", [])) > 0:     # use in_tiles to distinguish the single-tif and tile-tifs
                img_i["path"] = os.path.join(img_dir, img_i["image_id"]+".tif") 
                img_i["in_offsets"] = tile_offsets_f
                imginfo_by_tiles.append(img_i)
        assert len(imginfo_by_tiles) > 0, "Error: --image_stitch is enabled, but no image tile information is provided. Please provide the image tile information in the input YAML file"
        cmds.extend(cmd_image_stitch(imginfo_by_tiles, args, env))

    # =========
    #  Downstream
    # =========
    if args.run_ficture or args.run_cartload or args.image_png2pmtiles or args.upload_aws:
        runinfo=[]
        if args.in_yaml:
            avail_runs=[run_i.get("run_id") for run_i in yml.get("RUNS", [])]
            if len(args.run_ids) == 0 and len(avail_runs) == 1:
                args.run_ids = avail_runs
            assert len(args.run_ids) > 0, "When --run-ficture, --run-cartload-join, --image-png2pmtiles, or --upload-aws is enabled, --run-ids is required"
            
            for run_i in yml.get("RUNS", []):
                if run_i.get("run_id") in args.run_ids:
                    # in/out directory
                    run_i["sge_dir"]=sgeinfo["sge_dir"]
                    run_i["run_dir"]=os.path.join(out_dir, run_i["run_id"])
                    # sge info
                    run_i["filter_by_density"] = run_i.get("filter_by_density", False) # If filter_by_density is required for a run, the user should define it in the run section
                    run_i["filtered_prefix"] = sgeinfo["filtered_prefix"]
                    run_i["major_axis"] = run_i.get("major_axis", "X")
                    run_i["colname_count"] = sgeinfo["colname_count"]
                    # external model when applied
                    run_i["ext_path"]=run_i.get("ext_path", None)
                    run_i["ext_id"]=run_i.get("ext_id", None)
                    run_i["copy_ext_model"]=run_i.get("copy_ext_model", False)
                    run_i["cmap"]=run_i.get("cmap", args.cmap)
                    # analysis parameters
                    run_i["train_width"] = run_i.get("train_width", None)     
                    run_i["n_factor"] = run_i.get("n_factor", None)           # will fill in the default value in lda
                    # customizing feature filtering 
                    if run_i.get("include_feature_type_regex", None) and run_i.get("feature_type_ref", None) is None:
                        if sgeinfo.get("csv_colname_feature_type", None) is None:
                            raise ValueError("Error: --csv-colname-feature-type is required when --include-feature-type-regex is used.")
                        else:
                            run_i["feature_type_ref"] = sgeinfo["csv_colname_feature_type"]
                    runinfo.append(run_i)
        else:
            assert len(args.run_ids) == 1, "When --in-yaml is not provided, only one run ID is allowed"            
            # create a dictionary for the input    
            runinfo=[
                {   
                    # run ID
                    "run_id": args.run_ids,
                    # in/out directory
                    "sge_dir": args.sge_dir,
                    "run_dir": os.path.join(out_dir, args.run_ids), 
                    # sge info
                    "filter_by_density": args.filter_by_density,
                    "filtered_prefix": args.filtered_prefix,
                    "major_axis": args.major_axis,
                    "colname_count": args.colname_count,
                    # external model when applied
                    "ext_id": args.ext_id,
                    "ext_path": args.ext_path,
                    "copy_ext_model": args.copy_ext_model,
                    "cmap": args.cmap,
                    # analysis parameters
                    "train_width": args.train_width,
                    "n_factor": args.n_factor,
                    "histology":[],
                }
            ] 
        #print(runinfo) 
        for run_i in runinfo:
            os.makedirs(run_i["run_dir"], exist_ok=True)

            if args.run_ficture:
                fic_dir = os.path.join(run_i["run_dir"], f"ficture{'' if fic_v == '1' else '2'}")
                os.makedirs(fic_dir, exist_ok=True)
                cmds.extend(cmd_run_ficture(run_i, fic_v, args, env))

            if args.run_cartload:
                cartl_dir = os.path.join(run_i["run_dir"], f"cartload{'' if fic_v == '1' else '2'}")
                os.makedirs(cartl_dir, exist_ok=True)
                cmds.append(cmd_run_cartload(run_i, fic_v, args, env))
                if args.upload_aws:
                    cmds.append(cmd_upload_aws(run_i, fic_v, args, env, "--upload-cartload-only"))

            if args.image_png2pmtiles:
                run_i["histology"] = imginfo
                #print(imginfo)
                cmds.extend(cmd_image_png2pmtiles(run_i, fic_v, args, env))
            
                if args.upload_aws:
                    cmds.append(cmd_upload_aws(run_i, fic_v, args, env, "--upload-basemap-only"))

    #print(cmds)

    if args.dry_run:
        for cmd in cmds:
            print("="*10)
            print(f"Command: {cmd}")
            print("="*10)
            os.system(f"{cmd} --dry-run")

    # TODO: Find a better way to automatically generate the job ID
    def generate_job_id(args, timestamp):
        # Define shorthand encoding
        action_code = "".join([
            "A" if args.sge_convert else "a" if args.sge_stitch else "",
            "B" if args.run_ficture  else "",
            "C" if args.run_cartload else "",
            "D" if args.image_png2pmtiles else "",
            "E" if args.upload_aws else ""
        ])
        # Construct job_id based on run_ids
        if not any([args.run_ficture, args.run_cartload, args.image_png2pmtiles, args.upload_aws]):
            args.job_id = f"sge_convert_{timestamp}" if args.sge_convert else f"sge_stitch_{timestamp}"
        elif len(args.run_ids) == 1:
            args.job_id = f"{args.run_ids[0]}_{action_code}_{timestamp}"
        else:
            shared_text = os.path.commonprefix(args.run_ids)
            indiv_text = [run_id.replace(shared_text, "") for run_id in args.run_ids]
            args.job_id = f"{shared_text}-{'.'.join(indiv_text)}_{action_code}_{timestamp}"
        return args.job_id

    if args.job_id is None:
        args.job_id = generate_job_id(args, timestamp)

    job_f = write_jobfile(args.job_id, log_dir, env_cmds+cmds, slurm, args.submit_mode)

    if args.submit:
        submit_job(job_f, args.submit_mode)



if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    #print(sys.argv)
    func(sys.argv[1:]) 


## 250731: Need to update the cmds related to images.
import os
import sys
import subprocess
import argparse
import yaml
import datetime
from types import SimpleNamespace

from cartloader.utils.utils import add_param_to_cmd, scheck_file
from cartloader.utils.orient_helper import update_orient
from cartloader.utils.sge_helper import aux_sge_args
from cartloader.utils.execution import write_jobfile, submit_job

timestamp = datetime.datetime.now().strftime("%y%m%d")

img_suffix=[".tif", ".tiff", ".png"]
img_suffix.extend([suffix.upper() for suffix in img_suffix])

aux_env_args = {
        "sge_stitch": ["spatula"],
        "sge_convert": ["spatula", "gzip", "pigz", "pigz_threads"],
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

def subset_config(base_config, keys, prefix=None):
    """
    YAML-only mode: build a SimpleNamespace from base_config for the selected keys.
    The args parameter is retained for call-site compatibility but ignored.
    """
    src = base_config.get(prefix, {}) if prefix else base_config
    config = {}
    for key in keys:
        if key in src:
            val = src.get(key)
            if isinstance(val, list) and len(val) == 0:
                continue
            if val is None:
                continue
            config[key] = val
    return SimpleNamespace(**config)

def validate_yaml_config(yml, args):
    """Minimal validation for YAML-only mode.
    - out_dir must be provided (either CLI override or YAML)
    - If SGE-related actions requested, SGE section must exist
    - If run-related actions requested, RUNS must exist and optionally contain run_ids
    - If image actions requested, HISTOLOGY (or in_tiles under SGE for sge_stitch) must exist
    """
    # out_dir
    if args.out_dir is None and yml.get("out_dir", None) is None:
        raise AssertionError("Error: --out-dir or out_dir in YAML is required")

    need_sge = args.sge_convert or args.sge_stitch or args.run_ficture or args.run_cartload
    if need_sge:
        assert "SGE" in yml and isinstance(yml["SGE"], dict), "Error: YAML must contain SGE section for SGE/run actions"

    need_runs = args.run_ficture or args.run_cartload
    if need_runs:
        runs = yml.get("RUNS", [])
        assert isinstance(runs, list) and len(runs) > 0, "Error: YAML must contain RUNS with at least one run when run actions are requested"
        if args.run_ids:
            in_yaml_ids = {str(r.get("run_id")) for r in runs}
            missing = [rid for rid in args.run_ids if str(rid) not in in_yaml_ids]
            assert not missing, f"Error: run_ids not found in YAML RUNS: {', '.join(missing)}"

    need_hist = args.image_png2pmtiles or args.image_stitch
    if need_hist:
        # HISTOLOGY is optional for image_stitch by tiles when using SGE.in_tiles, but required for image_png2pmtiles
        if args.image_png2pmtiles:
            hist = yml.get("HISTOLOGY", [])
            assert isinstance(hist, list) and len(hist) > 0, "Error: YAML must contain HISTOLOGY for image_png2pmtiles"
        # for image_stitch, either HISTOLOGY with in_tiles or SGE.in_tiles will be used downstream
    return True

class StepinatorRunner:
    """Thin orchestrator to hold state and plan commands.
    Keeps handlers in this file; only centralizes plan, env, slurm, and job rendering.
    """
    def __init__(self, args, yml, env_ns: SimpleNamespace, slurm_ns: SimpleNamespace,
                out_dir: str, log_dir: str, env_cmds):
        self.args = args
        self.yml = yml
        self.env = env_ns
        self.slurm = slurm_ns
        self.out_dir = out_dir
        self.log_dir = log_dir
        self.env_cmds = env_cmds
        self.plan: list[str] = []

    def add(self, cmds):
        if cmds is None:
            return
        if isinstance(cmds, str):
            if cmds.strip():
                self.plan.append(cmds)
        elif isinstance(cmds, (list, tuple)):
            for c in cmds:
                self.add(c)
        else:
            # ignore unknown types
            pass

    def print_plan(self):
        print("\n===== Dry Run: Planned Commands =====")
        for i, cmd in enumerate(self.plan, 1):
            print(f"[{i:02d}] {cmd}")
        print("===== End Plan =====\n")

    def _generate_job_id(self, timestamp):
        args = self.args
        action_code = "".join([
            "A" if args.sge_convert else "a" if args.sge_stitch else "",
            "B" if args.run_ficture  else "",
            "C" if args.run_cartload else "",
            "D" if args.image_png2pmtiles else "",
            "E" if args.upload_aws else ""
        ])
        if not any([args.run_ficture, args.run_cartload, args.image_png2pmtiles, args.upload_aws]):
            job_id = f"sge_convert_{timestamp}" if args.sge_convert else f"sge_stitch_{timestamp}"
        elif len(args.run_ids) == 1:
            job_id = f"{args.run_ids[0]}_{action_code}_{timestamp}"
        else:
            shared_text = os.path.commonprefix(args.run_ids)
            indiv_text = [run_id.replace(shared_text, "") for run_id in args.run_ids]
            job_id = f"{shared_text}-{'.'.join(indiv_text)}_{action_code}_{timestamp}"
        return job_id

    def write_and_submit(self):
        args = self.args
        if args.job_id is None:
            args.job_id = self._generate_job_id(timestamp)
        job_f = write_jobfile(args.job_id, self.log_dir, self.env_cmds + self.plan + ["\n"], self.slurm, args.submit_mode)
        if args.submit:
            submit_job(job_f, args.submit_mode)
        return job_f

def cmd_sge_stitch(sgeinfo, args, env, generate_tile_minmax_only=False):
    # collect tiles 
    stitch_cmd = " ".join([
        "cartloader", "sge_stitch",
        f"--makefn sge_stitch.{timestamp}.mk",
        f"--in-tiles {sgeinfo['in_tiles_str']}",
        f"--out-dir {sgeinfo['sge_dir']}",
        f"--units-per-um {sgeinfo.get('units_per_um', None)}" if sgeinfo.get("units_per_um", None) else "",
        f"--colname-count {sgeinfo.get('colname_count', None)}" if sgeinfo.get('colname_count', None) else "",
        f"--sge-visual" if args.sge_visual else "",
        f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        f"--restart" if args.restart else "",
        "--generate-tile-minmax-only" if generate_tile_minmax_only else "",
        f"--filter-by-density --out-filtered-prefix {sgeinfo['filtered_prefix']} --genomic-feature {sgeinfo['colname_count']}" if sgeinfo['filter_by_density'] else "",
        f"--north-up" if args.north_up else "",
    ])

    # add aux stitch parameters
    stitch_aug = subset_config(sgeinfo, aux_params_args["sge_stitch"], prefix=None)
    stitch_cmd = add_param_to_cmd(stitch_cmd, stitch_aug, aux_params_args["sge_stitch"])

    # add aux north-up parameters
    if args.north_up:
        northup_aug = subset_config(sgeinfo, aux_params_args["north_up"], prefix=None)
        stitch_cmd = add_param_to_cmd(stitch_cmd, northup_aug, aux_params_args["north_up"])
    
    # add aux env
    stitch_cmd = add_param_to_cmd(stitch_cmd, env, aux_env_args["sge_stitch"]+aux_env_args["north_up"]) if args.north_up else add_param_to_cmd(stitch_cmd, env, aux_env_args["sge_stitch"])
    return stitch_cmd

def cmd_sge_convert(sgeinfo, args, env):

    platform= sgeinfo.get("platform", None)
    assert platform is not None, "Please provide platform information for sge_convert"

    # define input args
    if platform == "10x_visium_hd":
        assert sgeinfo.get("in_mex", None) is not None, "Please provide in_mex in SGE section"
        assert sgeinfo.get("pos_parquet", None) is not None, "Please provide pos_parquet in SGE section"
        assert sgeinfo.get("scale_json", None) is not None or sgeinfo.get("units_per_um", None) is not None, "Please provide scaling information by units_per_um or scale_json in SGE section."
        in_arg = " ".join([
            f"--in-mex {sgeinfo['in_mex']}",
            f"--pos-parquet {sgeinfo['pos_parquet']}",
            f"--scale-json {sgeinfo['scale_json']}" if sgeinfo.get("scale_json", None) is not None else "",
            f"--units-per-um {sgeinfo['units_per_um']}" if sgeinfo.get("scale_json", None) is None and sgeinfo.get("units_per_um", None) is not None else ""
        ])
    elif platform == "seqscope":
        assert sgeinfo.get("in_mex", None) is not None, "Please provide in_mex for seqscope"
        in_arg = f"--in-mex {sgeinfo['in_mex']}"
    elif platform == "10x_xenium":
        assert sgeinfo.get("in_csv", None) is not None or sgeinfo.get("in_parquet", None) is not None , "Please provide in_csv or in_parquet in SGE section"
        if sgeinfo.get("in_csv", None) is not None:
            in_arg = f"--in-csv {sgeinfo['in_csv']}"
        else:
            in_arg = f"--in-parquet {sgeinfo['in_parquet']}"
    elif platform in ["bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic"]:
        assert sgeinfo.get("in_csv", None) is not None, f"Please provide in_csv in SGE section"
        in_arg = f"--in-csv {sgeinfo['in_csv']}"
    else:
        raise AssertionError(f"Unsupported platform: {platform}")

    # add csv other columns
    if platform in ["bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "10x_xenium", "generic"]:
        csv_colnames_others= sgeinfo["csv_colnames_others"].split(",") if sgeinfo.get("csv_colnames_others", None) is not None else []
        if csv_colnames_others:
            in_arg += " " + " ".join(["--csv-colnames-others"] + csv_colnames_others)

    format_cmd = " ".join([
        "cartloader", "sge_convert",
        f"--makefn sge_convert.{timestamp}.mk",
        f"--platform {sgeinfo['platform']}",
        in_arg,
        f"--out-dir {sgeinfo['sge_dir']}",
        f"--colname-count {sgeinfo['colname_count']}",
        f"--filter-by-density --out-filtered-prefix {sgeinfo['filtered_prefix']} --genomic-feature {sgeinfo['colname_count']}" if sgeinfo['filter_by_density'] else "",
        f"--sge-visual" if args.sge_visual else "",
        f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        f"--north-up" if args.north_up else "",
    ])
    # add aux parameters
    format_aug = subset_config(sgeinfo, aux_params_args["sge_convert"], prefix=None)
    format_cmd = add_param_to_cmd(format_cmd, format_aug, aux_params_args["sge_convert"])

    # add aux north-up parameters
    if args.north_up:
        northup_aug = subset_config(sgeinfo, aux_params_args["north_up"], prefix=None)
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
    tsv = os.path.join(sge_dir, sgefn["tsv"])
    ftr = os.path.join(sge_dir, sgefn["ftr"])
    minmax = os.path.join(sge_dir, sgefn["minmax"])
    
    cstsv = os.path.join(fic_dir, sgefn["cstsv"])

    if fic_v == "1":
        sge_arg=f"--in-transcript {tsv} --in-feature {ftr} --in-minmax {minmax} --in-cstranscript {cstsv}"
    else:
        sge_arg=f"--in-transcript {tsv} --in-feature {ftr} --in-minmax {minmax}"
    return sge_arg

def cmd_run_ficture(run_i, fic_v, args, env):
    ficture_cmds = []
    ext_path = run_i["ext_path"]
    ext_id   = run_i["ext_id"]
    
    # makefile 
    mkbn = f"run_ficture{fic_v}" if ext_path is None else f"run_ficture{fic_v}_{ext_id}"

    # dirs
    fic_dir=run_i["fic_dir"]

    # sge
    sge2fn = define_sge2fn(run_i["filter_by_density"], run_i["filtered_prefix"])
    sge_arg = define_sge_arg(sge2fn, run_i["sge_dir"], fic_dir, fic_v)
    in_dist = os.path.join(run_i["sge_dir"], "feature.distribution.tsv.gz")

    # * check files
    if not args.sge_convert and not args.sge_stitch:
        if not args.lenient:
            for file in [sge2fn["tsv"], sge2fn["ftr"], sge2fn["minmax"]]:
                scheck_file(os.path.join(run_i["sge_dir"], file))
    
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
            f"--makefn {mkbn}.{timestamp}.mk",
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
            # f"--filter-by-overlapping-features --in-feature-dist {in_dist}" if run_i["filter_by_overlapping_features"] else "",
            # f"--min-ct-per-ftr-tile {run_i['min_ct_per_ftr_tile']}" if run_i["filter_by_overlapping_features"] and run_i["min_ct_per_ftr_tile"] > 0 else "",
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
            f"--makefn {mkbn}.{timestamp}.mk",
            f"--out-dir {fic_dir}",
            sge_arg,
            f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
            "--main",
            f"--width {run_i['train_width']}",
            f"--n-factor {run_i['n_factor']}",
            # f"--filter-by-overlapping-features --in-feature-dist {in_dist}" if run_i["filter_by_overlapping_features"] else "",
            # f"--min-ct-per-ftr-tile {run_i['min_ct_per_ftr_tile']}" if run_i["filter_by_overlapping_features"] and run_i["min_ct_per_ftr_tile"] > 0 else "",
            cmap_arg,
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
            f"--threads {args.threads}"
        ])
    
    # add aux env/tools
    ficture_cmd = add_param_to_cmd(ficture_cmd, env, aux_env_args[f"run_ficture{fic_v}"])
    # add aux parameters
    ficture_aug = subset_config(run_i, aux_params_args[f"run_ficture{fic_v}"], prefix=None)  # merge auxiliary parameters
    ficture_cmd = add_param_to_cmd(ficture_cmd, ficture_aug, aux_params_args[f"run_ficture{fic_v}"])
    print(ficture_cmd)
    ficture_cmds.append(ficture_cmd)
    return ficture_cmds

def cmd_run_cartload(run_i, cartl_v, args, env):

    mkbn = f"run_cartload{cartl_v}.{timestamp}"

    fic_dir = run_i["fic_dir"]
    cartload_dir = run_i["cartl_dir"]
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

def update_flip_to_vhflip(val):
    vertical = "True" if str(val).lower() == "vertical" else "False"
    horizontal = "True" if str(val).lower() == "horizontal" else "False"
    return vertical, horizontal

def generate_in_tiles_str(in_tiles, param_keys):
    assert len(in_tiles) > 0, "Error: --in-tiles should have at least one tile"
    tile_args = []
    for tile in in_tiles:
        values = []
        for key in param_keys:
            val = tile.get(key, "")
            if val in [None, "null", "Null"]:
                val = ""
            if key == "georef_bounds" and val:
                val = val.replace("(","").replace(")","").replace(",","_")
            if key == "flip":
                vertical, horizontal = update_flip_to_vhflip(val)
                values.append(vertical)
                values.append(horizontal)
            elif key == "row" or key == "col":
                values.append(str(val))
            else:
                values.append(str(val))
        tile_args.append(",".join(values))
    in_tiles_str = " ".join(tile_args)
    return in_tiles_str

def cmd_image_stitch(imginfo_by_tiles, args, env):
    image_stitch_cmds = []        
    for imginfo_i in imginfo_by_tiles:
        image_id = imginfo_i["image_id"]
        mkbn = f"image_stitch_{image_id}.{timestamp}" 

        # Build in-tiles entries: row,col,path,georef,georef_info,rotate,vertical_flip,horizontal_flip
        tile_strs = []
        for tile in imginfo_i.get("in_tiles", []):
            row = str(tile.get("row", ""))
            col = str(tile.get("col", ""))
            path = str(tile.get("path", ""))
            georef = str(tile.get("georeference", False))
            # georef_info
            georef_info = ""
            if tile.get("georef_bounds", None):
                bounds = str(tile.get("georef_bounds")).replace("(", "").replace(")", "").replace(",", "_")
                georef_info = f"bounds:{bounds}"
            elif tile.get("georef_tsv", None):
                georef_info = f"bounds_pixel_tsv:{tile.get('georef_tsv')}"
            elif tile.get("georef_bounds_tsv", None):
                georef_info = f"bounds_tsv:{tile.get('georef_bounds_tsv')}"

            rotate = str(tile.get("rotate", "")) if tile.get("rotate", None) is not None else ""
            vflip, hflip = update_flip_to_vhflip(tile.get("flip", None))
            tile_strs.append(
                ",".join([row, col, path, georef, georef_info, rotate, vflip, hflip])
            )

        in_tiles_str = " ".join(tile_strs)

        image_stitch_cmd = " ".join([
            "cartloader", "image_stitch",
            f"--output {imginfo_i['path']}",
            f"--makefn {mkbn}.mk",
            f"--in-offsets {imginfo_i['in_offsets']}",
            f"--in-tiles {in_tiles_str}",
            f"--crop-tile-by-minmax",
            f"--in-minmax {imginfo_i['in_minmax']}",
            f"--restart" if args.restart else "",
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        ])
        # add aux tools
        image_stitch_cmd = add_param_to_cmd(image_stitch_cmd, env, aux_env_args["image_stitch"], underscore2dash=False)
        image_stitch_cmds.append(image_stitch_cmd)
    return image_stitch_cmds

# Update the orientation of a histology image
def update_orient_in_histology(histology):
    hist_path = histology["path"]

    rotation = histology.get("rotate", None)
    
    flip = histology.get("flip", None)
    flip_vertical = flip in [ "vertical", "both"]
    flip_horizontal = flip in ["horizontal", "both"]

    new_rot, new_vflip, new_hflip = update_orient(rotation, flip_vertical, flip_horizontal, hist_path)

    # Update the histology dictionary with the best solution
    histology["rotate"] = new_rot
    if new_vflip and new_hflip:
        histology["flip"] = "both"
    elif new_vflip:
        histology["flip"] = "vertical"
    elif new_hflip:
        histology["flip"] = "horizontal"
    else:
        histology["flip"] = None
    return histology

def cmd_image_png2pmtiles(run_i, args, env):    
    assert len(run_i.get("histology", [])) > 0, "Error: --histology is Required when running --image-png2pmtiles"

    cartload_dir = run_i["cartl_dir"]

    img_cmds=[]
    for histology in run_i.get("histology", []):
        # 1. histology tif to pmtiles
        img_path = histology["path"]
        # prefix
        img_inname = os.path.basename(img_path).replace("_", "-")
        for suffix in img_suffix:
            if img_inname.endswith(suffix):
                img_inname = img_inname[:-len(suffix)]
                break
        
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
            f"--georef-pixel-tsv {histology.get('georef_tsv', None)}" if histology.get("georeference", False) and histology.get("georef_tsv", None) is not None else "",
            f"--georef-bounds {histology.get('georef_bounds', None)}" if histology.get("georeference", False) and histology.get("georef_bounds", None) is not None else "",
            f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
            f"--restart" if args.restart else "",
        ])
        img_cmd = add_param_to_cmd(img_cmd, env, aux_env_args["image_png2pmtiles"], underscore2dash=False)
        img_cmds.append(img_cmd)

    return img_cmds

def cmd_upload_aws(run_i, args, env, aws_bucket, additional_args=None):
    # Option 1: use "." to locate the file and handle by perl -lane
    # aws_cmd="\n".join([
    #     f"out={cartload_dir}",
    #     f"id={run_i['run_id']}",
    #     f"aws_bucket={args.aws_bucket}",
    #     "aws s3 cp ${out}/catalog.yaml s3://${aws_bucket}/${id}/catalog.yaml",
    #     "grep -E '\.' ${out}/catalog.yaml | perl -lane 'print $F[$#F]' | xargs -I {} aws s3 cp ${out}/{} s3://${aws_bucket}/${id}/{}"
    # ])
    # Option 2: use cartloader upload_aws.py
    aws_cmd=" ".join([
        "cartloader", "upload_aws",
        f'--in-dir {run_i["cartl_dir"]}',
        f"--s3-dir \"s3://{aws_bucket}/{run_i['run_id']}\"",
        additional_args if additional_args is not None else "",
        f"--n-jobs {args.n_jobs}" if args.n_jobs else "",
        f"--restart" if args.restart else ""
    ])
    aws_cmd=add_param_to_cmd(aws_cmd, env, aux_env_args["upload_aws"])

    return aws_cmd

def stepinator(_args):
    parser = argparse.ArgumentParser(description="""
    Stepinator: A tool to run steps in cartloader to do SGE format conversion (--sge-convert) or downstream process (--run-ficture, --run-cartload, --image-png2pmtiles, --upload-aws) in cartloader.
    Input, actions, parameters, and tools is provided using a YAML file.
    Job execution modes: local and slurm.
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

    # * commands 
    cmd_params = parser.add_argument_group("Commands", "Commands to be performed")
    cmd_params.add_argument("--sge-stitch", action="store_true", help="Run sge-stitch in cartloader")
    cmd_params.add_argument("--list-overlapping-features", action="store_true", help="(Optional if --sge-stitch) List overlapping features in sge-stitch")
    cmd_params.add_argument("--sge-convert", action="store_true", help="Run sge-convert in cartloader")
    cmd_params.add_argument('--sge-visual', action='store_true', help='Plot the SGE from sge-convert or sge-stitch in a PNG file')
    cmd_params.add_argument("--north-up", action="store_true", help="(Optional if --sge-visual) Set the north direction to up in the SGE visual")
    cmd_params.add_argument("--run-ficture", action="store_true", help="Run run-ficture in cartloader.")
    cmd_params.add_argument("--run-cartload", action="store_true", help="Run run-cartload in cartloader.")
    cmd_params.add_argument("--ficture-version", type=str, default=None, help="Define which ficture version you want to use.", choices=["1","2"])
    cmd_params.add_argument("--image-stitch", action="store_true", help="Run image stitch in cartloader. This is for stitching images.")
    cmd_params.add_argument("--image-png2pmtiles", action="store_true", help="Run image-png2pmtiles in cartloader. This s provide image using --in-yaml or --histology.")
    cmd_params.add_argument("--upload-aws", action="store_true", help="Upload files to AWS S3 bucket. This s provide AWS bucket name using --in-yaml.")
    cmd_params.add_argument("--copy-ext-model", action="store_true", help="(Optional if --run-ficture1) Auxiliary action parameters for run-ficture. Copy external model when running FICTURE with an external model")
    cmd_params.add_argument("--init-ext", action="store_true", help="(Optional if --run-ficture1) Auxiliary action parameters for run-ficture. Only Initialize external model without run main-ext")
    cmd_params.add_argument("--skip-coarse-report", action="store_true", help="(Optional if --run-ficture1)Auxiliary action parameters for run-ficture. Skip coarse report")
    cmd_params.add_argument('--segment-10x', action='store_true', help='(Optional if --run-ficture1) Perform hexagon segmentation into 10x Genomics format')

    # * Key
    key_params = parser.add_argument_group(
    "IN/OUT Configuration",
    "Specify input data and settings to provide inputs via a YAML file, enabling batch processing of multiple runs.\n"
    )
    key_params.add_argument("--in-yaml", '-i', type=str, default=None, help="Input yaml file")
    key_params.add_argument("--run-ids", nargs="*", default=[], help="One or more run IDs. Required if --run-ficture, --run-cartload-join, --image-png2pmtiles, or --upload-aws is applied. ")
    key_params.add_argument("--image-ids", nargs="*", default=[], help="One or more image IDs. Required if --image-png2pmtiles is set.")
    # for sge_convert
    dev_params = parser.add_argument_group("Developer", "Developer options")
    dev_params.add_argument('--lenient', action='store_true', help='(Only for development usage) Lenient mode. Skip the required file/app existence check')

    args = parser.parse_args(_args)

    # =========
    #  Sanity check
    # =========
    assert args.sge_stitch or args.sge_convert or args.run_ficture or args.run_cartload or args.image_png2pmtiles or args.upload_aws, "Error: At least one action is required"

    # sge
    assert not (args.sge_convert and args.sge_stitch), "Error: --sge-convert and --sge-stitch cannot be applied together"
    assert not (args.sge_visual and not args.sge_convert and not args.sge_stitch), "Error: --sge-visual can only be applied with --sge-convert or --sge-stitch"
    assert not (args.north_up and not args.sge_visual), "Error: --north-up can only be applied with --sge-visual"   

    # Force YAML-only input mode
    assert args.in_yaml is not None, "Error: YAML-only mode. Provide --in-yaml to specify inputs."

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

    # validate YAML for requested actions
    validate_yaml_config(yml, args)

    # output dir
    out_dir = yml.get("out_dir", None)
    assert out_dir is not None, "Error: --out-dir is required"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    sge_dir = os.path.join(out_dir, "sge")
    os.makedirs(sge_dir, exist_ok=True)

    log_dir = os.path.join(out_dir, "worklog")
    os.makedirs(log_dir, exist_ok=True)

    # =========
    #  env
    # =========
    # * tools (collect from yaml and args) 
    env = subset_config(yml, 
                       [item for sublist in aux_env_args.values() for item in sublist] + ["imagemagick"],
                        prefix="env")  

    # * slurm (collect from yaml and args)
    slurm = subset_config(yml, ["account", "partition", "mail_user", "cpus_per_task", "mem_per_cpu", "hours"], prefix="slurm")
    
    # * hpc modules (YAML-only)
    hpc_modules = yml.get("env", {}).get("hpc_modules", [])

    if len(hpc_modules) > 0:
        module_load_cmds="module load "+ " ".join(hpc_modules)
    else:
        module_load_cmds=""

    # * conda env (YAML-only)
    conda_base = yml.get("env", {}).get("conda", {}).get("base_path", None)
    conda_env  = yml.get("env", {}).get("conda", {}).get("env_prefix", None)

    conda_cmd=[]
    if conda_base:
        conda_cmd.append(f"source {conda_base}/etc/profile.d/conda.sh")
    if conda_env:
        conda_cmd.append(f"conda activate {conda_env}")

    # * aws (YAML-only)
    aws_bucket = yml.get("aws_bucket", None)

    # =========
    #  init 
    # =========
    env_cmds=[module_load_cmds] + conda_cmd
    # export path for imagemagick
    if args.sge_stitch or args.run_cartload:
        if getattr(env, "imagemagick", None) is not None:
            env_cmds.append("export PATH=$PATH:"+":".join([env.imagemagick]))
    env_cmds.append("\n")
    runner = StepinatorRunner(args, yml, env, slurm, out_dir, log_dir, env_cmds)

    # =========
    #  sgeinfo 
    # =========
    # sge infomation - shared across actions
    sgeinfo = yml.get("SGE", {})
    # filter_by_density
    sgeinfo["filter_by_density"] = sgeinfo.get("filter_by_density", False)
    sgeinfo["filtered_prefix"] = sgeinfo.get("filtered_prefix", "filtered")
    # colnames for output
    sgeinfo["colname_count"] = sgeinfo.get("colname_count", "count")
    # sge_dir    
    sgeinfo["sge_dir"] = sge_dir

    # =========
    #  SGE Convert / Stitch
    # =========
    if args.sge_convert:
        runner.add(cmd_sge_convert(sgeinfo, args, env))

    if args.sge_stitch:
        # YAML-only tiles
        in_tiles = sgeinfo.get("in_tiles", [])
        assert len(in_tiles) > 1, 'When --sge-stitch is enabled, at least two tiles are required.'
        sge_keys_tile = ["in_csv", "in_feature", "in_minmax", "row", "col", "rotate", "flip"]
        sgeinfo["in_tiles_str"] = generate_in_tiles_str(in_tiles, sge_keys_tile)
        runner.add(cmd_sge_stitch(sgeinfo, args, env))

    # =========
    #  IMAGE Stitch
    # =========
    # the stitched tif should be shared across run_id in the sge
    if args.image_png2pmtiles or args.image_stitch:
        imginfo = yml.get("HISTOLOGY", [])
        # normalize: accept hist_id alias for image_id
        for img_i in imginfo:
            if "image_id" not in img_i and "hist_id" in img_i:
                img_i["image_id"] = img_i.pop("hist_id")
        if len(args.image_ids) > 0:
            imginfo = [img_i for img_i in imginfo if img_i.get("image_id") in args.image_ids]
        # sanity check: no duplicate image_id in imginfo
        image_ids = [img_i.get("image_id") for img_i in imginfo]
        assert len(image_ids) == len(set(image_ids)), "Error: Found duplicate Image IDs in the input YAML file. Ensure unique image ID(s)"

    # image stitch
    if args.image_stitch:
        # img_dir to host the output tif from img stitch
        img_dir = os.path.join(out_dir, "hist")
        os.makedirs(img_dir, exist_ok=True)
        # add a cmd to make sure tile_offsets file exists - it was not generated for the early jobs
        tile_offsets_f = os.path.join(sgeinfo["sge_dir"], "coordinate_minmax_per_tile.tsv")
        if not os.path.exists(tile_offsets_f) and not args.sge_stitch:
            runner.add(cmd_sge_stitch(sgeinfo, args, env, generate_tile_minmax_only=True))
        # extract the img info for img stitch, and add the path
        imginfo_by_tiles = []
        for img_i in imginfo:
            if len(img_i.get("in_tiles", [])) > 0:
                img_i["path"] = os.path.join(img_dir, img_i["image_id"]+".tif") 
                img_i["in_offsets"] = tile_offsets_f
                img_i["in_minmax"] = tile_offsets_f
                imginfo_by_tiles.append(img_i)
        assert len(imginfo_by_tiles) > 0, "Error: --image_stitch is enabled without image tile information. Provide the image tile in the input YAML file"
        runner.add(cmd_image_stitch(imginfo_by_tiles, args, env))

    # =========
    #  Downstream
    # =========
    if args.run_ficture or args.run_cartload or args.image_png2pmtiles or args.upload_aws:

        # when only 1 run_id exists, no need to specify the run_ids
        avail_runs=[run_i.get("run_id") for run_i in yml.get("RUNS", [])]
        if len(args.run_ids) == 0 and len(avail_runs) == 1:
            args.run_ids = avail_runs
        
        assert len(args.run_ids) > 0, "When --run-ficture, --run-cartload-join, --image-png2pmtiles, or --upload-aws is enabled, --run-ids is required"
        
        for run_i in yml.get("RUNS", []):
            if run_i.get("run_id") in args.run_ids:
                # configure 
                # ----
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
                run_i["cmap"]=run_i.get("cmap", None)
                # analysis parameters
                run_i["train_width"] = run_i.get("train_width", None)     
                run_i["n_factor"] = run_i.get("n_factor", None)           # will fill in the default value in lda
                # customizing feature filtering 
                if run_i.get("include_feature_type_regex", None) and run_i.get("feature_type_ref", None) is None:
                    if sgeinfo.get("csv_colname_feature_type", None) is None:
                        raise ValueError("Error: --csv-colname-feature-type is required when --include-feature-type-regex is used.")
                    else:
                        run_i["feature_type_ref"] = sgeinfo["csv_colname_feature_type"]
                run_i["fic_dir"]    = os.path.join(run_i["run_dir"], f"ficture{'' if fic_v == '1' else '2'}")
                run_i["cartl_dir"]  = os.path.join(run_i["run_dir"], f"cartload{'' if fic_v == '1' else '2'}")

                # execution
                # ----
                os.makedirs(run_i["run_dir"], exist_ok=True)

                if args.run_ficture:
                    os.makedirs(run_i["fic_dir"], exist_ok=True)
                    runner.add(cmd_run_ficture(run_i, fic_v, args, env))

                if args.run_cartload:
                    os.makedirs(run_i["cartl_dir"], exist_ok=True)
                    runner.add(cmd_run_cartload(run_i, fic_v, args, env))
                    if args.upload_aws:
                        runner.add(cmd_upload_aws(run_i, args, env, aws_bucket, "--upload-cartload-only"))

                if args.image_png2pmtiles:
                    run_i["histology"] = imginfo
                    #print(imginfo)
                    runner.add(cmd_image_png2pmtiles(run_i, args, env))
                
                    if args.upload_aws:
                        runner.add(cmd_upload_aws(run_i, args, env, aws_bucket, "--upload-basemap-only"))

    if args.dry_run:
        runner.print_plan()

    runner.write_and_submit()

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    #print(sys.argv)
    func(sys.argv[1:]) 

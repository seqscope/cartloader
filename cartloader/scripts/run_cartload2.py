import sys, os, argparse, logging, subprocess, inspect
import pandas as pd
from pathlib import Path

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file, read_minmax, flexopen, execute_makefile, valid_and_touch_cmd
from cartloader.utils.color_helper import normalize_rgb
from cartloader.utils.ficture2_helper import ficture2_params_to_factor_assets

def parse_arguments(_args):
    """
    Package SGE and optional FICTURE outputs into PMTiles and a catalog for CartoScope (sge → optional ficture2/joins → pmtiles → catalog)
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(
        prog=f"cartloader run_cartload2",
        description="Package SGE and optional FICTURE outputs into PMTiles and a catalog for CartoScope (sge → optional ficture2/joins → pmtiles → catalog)"
    )

    run_params = parser.add_argument_group("Run Options", "Execution controls for generating and running the Makefile")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile but do not execute it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--makefn', type=str, default="run_cartload2.mk", help='Name of the generated Makefile (default: run_cartload2.mk)')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job for tippecanoe (default: 4)')
    run_params.add_argument('--log', action='store_true', default=False, help='Write logs to a file under the output directory')
    run_params.add_argument('--log-suffix', type=str, default=".log", help='Suffix for the log filename; final path is <out_dir>_cartload<suffix> (default: .log)')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Primary input and output locations")
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory (PMTiles, assets JSON, and catalog YAML)')
    inout_params.add_argument('--sge-dir', type=str, help='Path to SGE directory produced by "cartloader sge_convert"; must include SGE files and an assets JSON (see --in-sge-assets)')
    inout_params.add_argument('--fic-dir', type=str, help='Path tp FICTURE results directory produced by "cartloader run_ficture"; must include FICTURE results and a parameter JSON (see --in-fic-params)')
    inout_params.add_argument('--cell-assets', type=str, nargs="+", default=[], help='Optional list of cell asset JSON/YAML files produced by import_xenium_cell or import_visiumhd_cell')
    inout_params.add_argument('--background-assets', type=str, nargs="+", default=[], help='Optional list of background asset specs (JSON/YAML or inline id:path or id1:id2:path)')
    inout_params.add_argument('--square-assets', type=str, nargs="+", default=[], help='Optional list of square asset JSON/YAML files produced by import_visiumhd_square')

    key_params = parser.add_argument_group("Key Parameters", "Metadata")
    key_params.add_argument('--id', type=str, required=True, help='Identifier for the output assets; avoid whitespace (use "-" instead of "_")')
    key_params.add_argument('--title', type=str, help='Human-readable title for the output assets')
    key_params.add_argument('--desc', type=str, help='Short description of the output assets')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Advanced settings; defaults work for most cases")
    aux_params.add_argument('--in-sge-assets', type=str, default="sge_assets.json", help='File name of a SGE assets JSON/YAML under --sge-dir, providing paths to transcript, feature, minmax files (default: sge_assets.json)')
    aux_params.add_argument('--in-fic-params', type=str, default="ficture.params.json", help='File name of FICTURE params JSON/YAML under --fic-dir, providing FICTURE paramaters (default: ficture.params.json)')
    aux_params.add_argument('--out-fic-assets', type=str, default="ficture_assets.json", help='File name of output JSON/YAML for FICTURE asset metadata under --out-dir (default: ficture_assets.json)')
    aux_params.add_argument('--out-catalog', type=str, default="catalog.yaml", help='File name of output catalog YAML under --out-dir (default: catalog.yaml)')
    # sge scale
    # aux_params.add_argument('--sge-scale', type=int, default=1, help='Scale factor from input coordinates to output sge image pixels (default: 1)')
    # tippecanoe/PMTiles 
    aux_params.add_argument('--rename-x', type=str, default='x:lon', help='Column rename mapping for X axis in tippecanoe, format old:new (default: x:lon)')  
    aux_params.add_argument('--rename-y', type=str, default='y:lat', help='Column rename mapping for Y axis in tippecanoe, format old:new (default: y:lat)')  
    aux_params.add_argument('--colname-feature', type=str, default='gene', help='Column name for feature/gene (default: gene)')
    aux_params.add_argument('--colname-count', type=str, default='count', help='Column name for molecule counts (default: count)')
    aux_params.add_argument('--out-molecules-id', type=str, default='genes', help='Base name for output molecules PMTiles files (no directory)')
    aux_params.add_argument('--max-join-dist-um', type=float, default=0.1, help='Max distance (in µm) to associate molecules with decoded pixels (default: 0.1)')
    aux_params.add_argument('--join-tile-size', type=float, default=500, help='Tile size (in µm) when joining molecules with decoded pixels (default: 500)')
    aux_params.add_argument('--bin-count', type=int, default=50, help='Number of bins when splitting input molecules (default: 50)')
    aux_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum tile size in bytes for tippecanoe/PMTiles (default: 5000000)')
    aux_params.add_argument('--max-feature-counts', type=int, default=500000, help='Maximum features per tile for tippecanoe/PMTiles (default: 500000)')
    aux_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Tippecanoe point-density preservation threshold (default: 1024)')
    aux_params.add_argument('--umap-colname-factor', type=str, default='topK', help='Column name encoding the dominant factor assignment in a UMAP TSV (default: topK)')
    aux_params.add_argument('--umap-colname-x', type=str, default='UMAP1', help='Column name for the UMAP X coordinate (default: UMAP1)')
    aux_params.add_argument('--umap-colname-y', type=str, default='UMAP2', help='Column name for the UMAP Y coordinate (default: UMAP2)')
    aux_params.add_argument('--umap-min-zoom', type=int, default=0, help='Minimum zoom for generated UMAP PMTiles (default: 0)')
    aux_params.add_argument('--umap-max-zoom', type=int, default=18, help='Maximum zoom for generated UMAP PMTiles (default: 18)')
    # ?
    aux_params.add_argument('--skip-raster', action='store_true', default=False, help='Skip raster image generation (no GDAL/go-pmtiles required)')
    # tmp
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory (default: <out_dir>/tmp)')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate files instead of cleaning up')
    # img
    aux_params.add_argument('--transparent-below', type=int, help='Set pixels below this value to transparent for dark background (range: 0~255)')
    aux_params.add_argument('--transparent-above', type=int, help='Set pixels above this value to transparent for light background (range: 0~255)')
    aux_params.add_argument('--sge-scale', type=int, default=1, help='scales input coordinates to pixels in the output image (default: 1)')

    env_params = parser.add_argument_group("Env Parameters", "Tool paths (override defaults if needed)")
    # aux_params.add_argument('--magick', type=str, default=f"magick", help='Path to ImageMagick binary') # Disable this function. The user need to add the path to the ImageMagick binary directory to the PATH environment variable
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip-compatible binary (tip: use "pigz -p4" for speed)')
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    env_params.add_argument('--tippecanoe', type=str, default=f"tippecanoe", help='Path to tippecanoe binary') # default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", 
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    
    args=parser.parse_args(_args)

    # dir
    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(args.out_dir, "tmp")
    
    # env
    scheck_app(args.spatula)
    scheck_app(args.tippecanoe)
    if not args.skip_raster:
        scheck_app(args.pmtiles)
        scheck_app(args.gdal_translate)
        scheck_app(args.gdaladdo)

    return args

def pick_sge_inputs(args):
    # Try SGE JSON
    if args.sge_dir:
        sge_json = os.path.join(args.sge_dir, args.in_sge_assets)
        sge = load_file_to_dict(sge_json)
        if all(k in sge for k in ("transcript", "feature", "minmax")):
            return (sge["transcript"], sge["feature"], sge["minmax"], f"{sge_json} (provided by --sge-dir and --in-sge-assets)")
    # fallthrough, try FICTURE
    if args.fic_dir:
        fic_json = os.path.join(args.fic_dir, args.in_fic_params)
        fic = load_file_to_dict(fic_json).get("in_sge", {})
        need = ("in_transcript", "in_feature", "in_minmax")
        if all(k in fic for k in need):
            return (fic["in_transcript"], fic["in_feature"], fic["in_minmax"], f"{fic_json} (provided by --fic-dir and --in-fic-params)")
        # has fic_dir but missing keys, concise error
        missing = ",".join(k for k in need if k not in fic)
        raise KeyError(f"Path not provided for SGE. Missing keys {missing} in FICTURE JSON {fic_json} (provided by --fic-dir and --in-fic-params)")
    # neither source provided, actionable error
    raise KeyError("Path not provided for SGE. Provide using --sge-dir with --in-sge-assets or --fic-dir with --in-fic-params")

def copy_rgb_tsv(in_rgb, out_rgb, restart=False):
    def _get_content(path):
        with open(path, 'r') as f:
            hdrs = f.readline().rstrip().split("\t")
            col2idx = {hdr: i for i, hdr in enumerate(hdrs)}
            lines = ["\t".join(["Name", "Color_index", "R", "G", "B"])]
            for line in f:
                toks = line.rstrip().split("\t")
                if len(toks) != len(hdrs):
                    raise ValueError(f"Input RGB file {path} has inconsistent number of columns")
                rgb_r = float(toks[col2idx["R"]])
                rgb_g = float(toks[col2idx["G"]])
                rgb_b = float(toks[col2idx["B"]])
                rgb_r, rgb_g, rgb_b = normalize_rgb(rgb_r, rgb_g, rgb_b)
                name = toks[col2idx["Name"]]
                lines.append(f"{name}\t{name}\t{rgb_r:.5f}\t{rgb_g:.5f}\t{rgb_b:.5f}")
        return "\n".join(lines) + "\n"

    # Desired behavior:
    # - If restart is True OR output file does not exist: (re)generate the file.
    # - Otherwise: compare to expected output and skip rewrite if identical.
    expected_content = _get_content(in_rgb)

    if restart or not os.path.exists(out_rgb):
        with open(out_rgb, 'w') as f:
            f.write(expected_content)
        return

    try:
        with open(out_rgb, 'r') as f:
            existing_content = f.read()
        if existing_content == expected_content:
            return  # up-to-date; no rewrite needed
    except Exception:
        # On read/compare failure, fall through to regenerate
        pass

    with open(out_rgb, 'w') as f:
        f.write(expected_content)

def process_umap(umap, mm, args, out_prefix, model_id, fic_jsonf):
    """Add Makefile target to convert a UMAP bundle into PMTiles and copies."""
    if not umap:
        return None

    umap_tsv = umap.get("tsv", None)
    umap_png = umap.get("png", None)
    umap_idv_png = umap.get("ind_png", None)
    for f in (umap_tsv, umap_png, umap_idv_png):
        assert f is not None, f"UMAP entry incomplete for model {model_id} in FICTURE params {fic_jsonf} (provided by --fic-dir and --in-fic-params)"
        assert os.path.exists(f), f"UMAP file not found: {f} for model {model_id} in FICTURE params {fic_jsonf} (provided by --fic-dir and --in-fic-params)"

    cmds = cmd_separator([], f"Converting UMAP for {model_id} into PMTiles and copying relevant files..")
    prerequisites = [umap_tsv, umap_png, umap_idv_png]
    outfiles=[]

    umap_tsv_out = f"{out_prefix}-umap.tsv.gz"
    if umap_tsv.endswith(".gz"):
        cmds.append(f"cp {umap_tsv} {umap_tsv_out}")
    else:
        cmds.append(f"{args.gzip} -c '{umap_tsv}' > '{umap_tsv_out}'")
    outfiles.append(umap_tsv_out)

    cmds.append(f"cp {umap_png} {out_prefix}.umap.png")
    outfiles.append(f"{out_prefix}.umap.png")

    cmds.append(f"cp {umap_idv_png} {out_prefix}.umap.single.prob.png")
    outfiles.append(f"{out_prefix}.umap.single.prob.png")

    umap_ndjson = f"{out_prefix}-umap.ndjson"
    umap_pmtiles = f"{out_prefix}-umap.pmtiles"
    convert_cmd = " ".join([
        "cartloader", "render_umap",
        f"--input {umap_tsv_out}",
        f"--out {umap_ndjson}",
        f"--colname-factor {args.umap_colname_factor}",
        f"--colname-x {args.umap_colname_x}",
        f"--colname-y {args.umap_colname_y}"
    ])
    cmds.append(convert_cmd)

    tippecanoe_cmd = " ".join([
        f"TIPPECANOE_MAX_THREADS={args.threads}",
        f"'{args.tippecanoe}'",
        f"-t {args.tmp_dir}",
        f"-o {umap_pmtiles}",
        "-Z", str(args.umap_min_zoom),
        "-z", str(args.umap_max_zoom),
        "-l", "umap",
        "--force",
        "--drop-densest-as-needed",
        "--extend-zooms-if-still-dropping",
        "--no-duplication",
        f"--preserve-point-density-threshold={args.preserve_point_density_thres}",
        umap_ndjson
    ])
    cmds.append(tippecanoe_cmd)
    if not args.keep_intermediate_files:
        cmds.append(f"rm -f {umap_ndjson}")
    outfiles.append(umap_pmtiles)

    touch_flag_cmd=valid_and_touch_cmd(outfiles, f"{out_prefix}-umap.done") # this only touch the flag file when all output files exist
    cmds.append(touch_flag_cmd)
    target = f"{out_prefix}-umap.done"
    mm.add_target(target, prerequisites, cmds)
    return target

def run_cartload2(_args):
    """
    Build resources for CartoScope
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_dir + "_cartload" + args.log_suffix if args.log else None)
    logger.info("Analysis started")

    # start mm
    mm = minimake()

    # create output directory if needed
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.tmp_dir, exist_ok=True)
    
    # output files/prefix
    out_catalog_f = os.path.join(args.out_dir,args.out_catalog)
    out_assets_f = os.path.join(args.out_dir, args.out_fic_assets)
    out_molecules_prefix=os.path.join(args.out_dir,args.out_molecules_id)

    # 1. Load SGE metadata or FICTURE metadata to define and 
    
    in_molecules, in_features, in_minmax, src_hint = pick_sge_inputs(args)
    assert os.path.exists(in_molecules), f"File not found: {in_molecules} (transcript) {src_hint}"
    assert os.path.exists(in_features), f"File not found: {in_features} (feature) {src_hint}"
    assert os.path.exists(in_minmax), f"File not found: {in_minmax} (minmax) {src_hint}"
    
    # 2. deploy SGE
    if not args.skip_raster:
        ## create raster mono pmtiles for SGE
        cmds = cmd_separator([], f"Converting SGE counts into PMTiles")
        cmd = " ".join([
            "cartloader", "run_tsv2mono",
            "--in-tsv", in_molecules,
            "--in-minmax", in_minmax,
            "--out-prefix", f"{args.out_dir}/sge-mono",
            "--colname-count", args.colname_count,
            "--main",
            "--units-per-pixel", str(args.sge_scale),
            f"--pmtiles '{args.pmtiles}'",
            f"--gdal_translate '{args.gdal_translate}'",
            f"--gdaladdo '{args.gdaladdo}'",
            f"--spatula '{args.spatula}'",
            "--keep-intermediate-files" if args.keep_intermediate_files else "",
            # f"--sge-scale {args.sge_scale}" if args.sge_scale else "",
        ])
        cmds.append(cmd)
        # Use a flag to make sure both light and dark pmtiles are done 
        tsv2mono_flag = f"{args.out_dir}/sge-mono.done" 
        cmds.append(f"[ -f {args.out_dir}/sge-mono-dark.pmtiles.done ] && [ -f {args.out_dir}/sge-mono-light.pmtiles.done ] && touch {tsv2mono_flag}")
        mm.add_target(tsv2mono_flag, [in_molecules, in_minmax], cmds)

    ## 3. deploy FICTURE results
    join_pixel_tsvs = []
    join_pixel_ids = []
    if args.fic_dir is not None:
        fic_jsonf = os.path.join(args.fic_dir, args.in_fic_params)
        fic_data = load_file_to_dict(fic_jsonf)
        in_fic_params = fic_data.get("train_params", [])
        if len(in_fic_params) == 0:  # parameters are empty
            logger.error(f"FICTURE 'train_params' is empty after loading {fic_jsonf} (provided by --fic-dir and --in-fic-params)")

        # create the output assets json
        out_fic_assets = ficture2_params_to_factor_assets(in_fic_params, args.skip_raster)

        train_targets = []

        minmax = read_minmax(in_minmax, "row")
        xmin = minmax["xmin"]
        xmax = minmax["xmax"]
        ymin = minmax["ymin"]
        ymax = minmax["ymax"]

        for train_param in in_fic_params:
            model_id = train_param["model_id"]
            model_rgb = train_param["cmap"]
            train_width = train_param["train_width"]
            in_prefix = f"{args.fic_dir}/{model_id}"

            ## ?? what is factormap?
            factormap_path = train_param.get("factor_map", None)
            
            out_id = model_id.replace("_", "-")
            out_prefix = f"{args.out_dir}/{out_id}"

            out_train_id = out_id
            out_train_prefix = out_prefix

            train_inout ={
                "model":{
                    "required": True,
                    "in": train_param.get("model_path", f"{in_prefix}.model.tsv"),
                    "out": f"{out_prefix}-model.tsv"
                },
                "rgb":{
                    "required": True,
                    "in":  model_rgb,
                    "out": f"{out_prefix}-rgb.tsv"
                },
                "de":{
                    "required": False,
                    "in": train_param.get("de_path", f"{in_prefix}.bulk_chisq.tsv"),
                    "out": f"{out_prefix}-bulk-de.tsv"
                },
                "info":{
                    "required": False,
                    "in": train_param.get("info_path", f"{in_prefix}.factor.info.tsv"),
                    "out": f"{out_prefix}-info.tsv"
                }
            }

            if factormap_path is not None:
                train_inout["factor_map"] = {
                    "required": False,
                    "in":  factormap_path,
                    "out": f"{out_prefix}-factor-map.tsv"
                }

            # umap
            umap = train_param.get("umap", {})
            # if umap is a dict,
            if train_param.get("analysis") == "multi-sample":
                process_umap(umap.get("shared"), mm, args, out_prefix+"-shared", model_id, fic_jsonf)
                process_umap(umap.get("sample"), mm, args, out_prefix, model_id, fic_jsonf)
            else:
                process_umap(umap, mm, args, out_prefix, model_id, fic_jsonf)
                
            cmds = cmd_separator([], f"Converting LDA-trained factors {model_id} into PMTiles and copying relevant files..")
            
            prerequisites = []
            outfiles=[]

            # fit_results
            in_fit_tsvf = train_param.get("fit_path", f"{in_prefix}.results.tsv.gz")
            if os.path.exists(in_fit_tsvf):
                cmd = " ".join([
                    "cartloader", "run_hex2pmtiles",
                    "--in-tsv", in_fit_tsvf,
                    "--out-prefix", out_prefix,
                    f"--hex-width {train_width}",
                    "--rename-column", args.rename_x, args.rename_y,
                    "--threads", str(args.threads),
                    "--max-tile-bytes", str(args.max_tile_bytes),
                    "--max-feature-counts", str(args.max_feature_counts),
                    "--preserve-point-density-thres", str(args.preserve_point_density_thres),
                    f"--tippecanoe '{args.tippecanoe}'",
                    f"--log --log-suffix '{args.log_suffix}'" if args.log else "",
                    f"--tmp-dir '{args.tmp_dir}'",
                    "--keep-intermediate-files" if args.keep_intermediate_files else ""
                ])
                cmds.append(cmd)
                prerequisites.append(in_fit_tsvf)
                outfiles.append(f"{out_prefix}.pmtiles")
            
            # mode/rgb/de/posterior/info
            for key, val in train_inout.items():
                if val["required"] or os.path.exists(val["in"]):
                    prerequisites.append(val["in"])
                    if key == "rgb":
                        copy_rgb_tsv(val["in"], val["out"], restart=args.restart)
                    else:
                        cmds.append(f"cp {val['in']} {val['out']}")
                    outfiles.append(val["out"])
            
            touch_flag_cmd=valid_and_touch_cmd(outfiles, f"{out_prefix}.done") # this only touch the flag file when all output files exist
            cmds.append(touch_flag_cmd)

            mm.add_target(f"{out_prefix}.done", prerequisites, cmds)

            ## add to the output asset json
            # out_decode_assets = [] ## NOT USED?

            sources = [f"{out_prefix}.done"]

            for decode_param in train_param["decode_params"]:
                in_id = decode_param["decode_id"]
                in_prefix = f"{args.fic_dir}/{in_id}"
                in_pixel_tsvf = decode_param.get("pixel_tsv_path", f"{in_prefix}.tsv.gz")
                in_pixel_png = decode_param.get("pixel_png_path", f"{in_prefix}.png")
                in_de_tsvf  = decode_param.get("de_tsv_path", f"{in_prefix}.bulk_chisq.tsv")
                in_post_tsvf = decode_param.get("pseudobulk_tsv_path", f"{in_prefix}.pseudobulk.tsv.gz")
                in_info_tsvf = decode_param.get("info_tsv_path", f"{in_prefix}.factor.info.tsv")

                out_id = in_id.replace("_", "-")
                out_prefix = os.path.join(args.out_dir, out_id)

                join_pixel_tsvs.append(in_pixel_tsvf)
                join_pixel_ids.append(out_id)

                cmds = cmd_separator([], f"Converting decoded factors {in_id} into PMTiles and copying relevant files.")
                outfiles=[]
                if not args.skip_raster:
                    cmd = " ".join([
                        "cartloader", "image_png2pmtiles",
                        "--georeference", "--geotif2mbtiles", "--mbtiles2pmtiles",
                        "--in-img", in_pixel_png,
                        f'--georef-bounds="{xmin},{ymin},{xmax},{ymax}"',
                        "--out-prefix", f"{out_prefix}-pixel-raster",
                        f"--pmtiles '{args.pmtiles}'",
                        f"--gdal_translate '{args.gdal_translate}'",
                        f"--gdaladdo '{args.gdaladdo}'",
                        "--keep-intermediate-files" if args.keep_intermediate_files else ""
                    ])
                    cmds.append(cmd)
                    outfiles.append(f"{out_prefix}-pixel-raster.pmtiles")

                cmds.append(f"cp {in_de_tsvf} {out_prefix}-bulk-de.tsv")
                cmds.append(f"cp {in_post_tsvf} {out_prefix}-pseudobulk.tsv.gz")
                cmds.append(f"cp {in_info_tsvf} {out_prefix}-info.tsv")
                
                outfiles.extend([f"{out_prefix}-bulk-de.tsv", f"{out_prefix}-pseudobulk.tsv.gz", f"{out_prefix}-info.tsv"])
    
                touch_flag_cmd=valid_and_touch_cmd(outfiles, f"{out_prefix}.done") # this only touch the flag file when all output files exist
                cmds.append(touch_flag_cmd)

                mm.add_target(f"{out_prefix}.done", [in_pixel_tsvf, in_pixel_png, in_de_tsvf, in_post_tsvf, model_rgb, in_info_tsvf], cmds)
                sources.append(f"{out_prefix}.done")

            cmds = cmd_separator([], f"Finishing up for train parameters {out_train_id}")
            cmds.append(f"touch {out_train_prefix}.alldone")
            mm.add_target(f"{out_train_prefix}.alldone", sources, cmds)
            train_targets.append(f"{out_train_prefix}.alldone")

        ## Touch a flag file for all ficture deployment
        cmds = cmd_separator([], f"Finalizing all conversions")
        cmds.append(f"touch {args.out_dir}/ficture.done")
        mm.add_target(f"{args.out_dir}/ficture.done", train_targets, cmds)

        ## Write the output asset json for FICTURE output
        logger.info("Writing a json file for expected FICTURE output assets")
        write_dict_to_file(out_fic_assets, out_assets_f, check_equal=True)

    ## 4. If FICTURE is provided with decoding results, join pixel-level TSVs
    molecules_f = in_molecules
    if ( len(join_pixel_tsvs) > 0 ):
        cmds = cmd_separator([], f"Joining pixel-level TSVs")
        out_join_pixel_prefix = f"{args.out_dir}/transcripts_pixel_joined"
        # cmd = " ".join([
        #         f"'{args.spatula}'", "paste-pixel-tsv",
        #         f"--out-tsv {out_join_pixel_prefix}.tsv.gz",
        #         f"--colname-x", args.rename_x.split(":")[0],
        #         f"--colname-y", args.rename_y.split(":")[0]
        #     ]
        #     + [ f"--pix-prefix-tsv {join_pixel_ids[i]}_,{join_pixel_tsvs[i]}" for i in range(len(join_pixel_tsvs)) ]
        # )
        cmd = " ".join([
                f"'{args.spatula}'", "join-pixel-decode",
                f"--out-prefix {out_join_pixel_prefix}",
                f"--mol-tsv {in_molecules}",
                f"--threads {args.threads}",
                f"--max-dist {args.max_join_dist_um}",
                f"--tile-size {args.join_tile_size}"
            ]
            + [ f"--decode-prefix-tsv {join_pixel_ids[i]}_,{join_pixel_tsvs[i]}" for i in range(len(join_pixel_tsvs)) ]
        )
        cmds.append(cmd)
        cmds.append(f"{args.gzip} -f {out_join_pixel_prefix}.tsv")
        mm.add_target(f"{out_join_pixel_prefix}.tsv.gz", [in_molecules]+join_pixel_tsvs, cmds)
        molecules_f = f"{out_join_pixel_prefix}.tsv.gz"

    ## 5. run tsv2pmtiles for the convert the joined pixel-level TSV to PMTiles
    cmds = cmd_separator([], f"Converting the joined pixel-level TSV to PMTiles")
    cmd = " ".join([
        "cartloader", "run_tsv2pmtiles",
        "--in-molecules", molecules_f,
        "--in-features", in_features,
        "--out-prefix", f"{out_molecules_prefix}",
        "--threads", str(args.threads),
        "--col-rename", args.rename_x, args.rename_y, f"feature:{args.colname_feature}", f"ct:{args.colname_count}",
        "--colname-feature", args.colname_feature,
        "--colname-count", args.colname_count,
        "--max-tile-bytes", str(args.max_tile_bytes),
        "--max-feature-counts", str(args.max_feature_counts),
        "--preserve-point-density-thres", str(args.preserve_point_density_thres),
        "--bin-count", str(args.bin_count),
        "--all",
        "--n-jobs", str(args.n_jobs),
        f"--log --log-suffix '{args.log_suffix}'" if args.log else "",
        f"--tippecanoe '{args.tippecanoe}'",
        f"--tmp-dir '{args.tmp_dir}'",
        "--keep-intermediate-files" if args.keep_intermediate_files else ""
    ])
    cmds.append(cmd)
    mm.add_target(f"{out_molecules_prefix}_pmtiles_index.tsv", [molecules_f], cmds)

    # 6. Create a yaml for all output assets
    cmds = cmd_separator([], f"Writing catalog YAML for all output assets")
    cmd = " ".join([
        "cartloader", "write_catalog_for_assets",
        "--sge-index", f"{out_molecules_prefix}_pmtiles_index.tsv", 
        "--sge-counts", f"{out_molecules_prefix}_bin_counts.json", 
        f"--fic-assets {out_assets_f}" if args.fic_dir else "",
        "--out-catalog", f"{out_catalog_f}",
        f"--log --log-suffix '{args.log_suffix}'" if args.log else ""
        ] + (["--overview", f"sge-mono-dark.pmtiles", "--basemap", f"sge:dark:sge-mono-dark.pmtiles", f"sge:light:sge-mono-light.pmtiles"] if not args.skip_raster else []
        ) + (["--background-assets"] + args.background_assets if args.background_assets else []
        ) + ([f"--cell-assets"] + args.cell_assets if args.cell_assets else []
        ) + ([f"--square-assets"] + args.square_assets if args.square_assets else [])
    )

    if ( args.id is not None ):
        cmd += f" --id {args.id}"

    if ( args.title is None ): ## if title is empty, use ID as title
        cmd += f" --title {args.id}"
    else:
        cmd += f" --title {args.title}"

    if ( args.desc is not None ):
        cmd += f" --desc {args.desc}"

    prerequisites_yaml=[f"{out_molecules_prefix}_pmtiles_index.tsv"]
    if args.fic_dir:
        prerequisites_yaml.append(out_assets_f)
    if not args.skip_raster:
        prerequisites_yaml.append(tsv2mono_flag)
    if args.background_assets:
        prerequisites_yaml.extend(args.background_assets)
    if args.cell_assets:
        prerequisites_yaml.extend(args.cell_assets)
    if args.square_assets:
        prerequisites_yaml.extend(args.square_assets)

    cmds.append(cmd)
    mm.add_target(f"{out_catalog_f}", prerequisites_yaml, cmds)

    if len(mm.targets) == 0:
        logger.error("No tasks were generated. Check inputs and parameters.")
        sys.exit(1)

    ## write makefile
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


def run_cartload2_generic(_args):
    """Backward compatibility: old command name forwards to run_cartload2."""
    return run_cartload2(_args)

    logger.info("Cartload2-generic finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

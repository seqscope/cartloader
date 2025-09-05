import sys, os, argparse, logging, subprocess, inspect
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file, ficture2_params_to_factor_assets, read_minmax, flexopen, execute_makefile, valid_and_touch_cmd

def parse_arguments(_args):
    """
    Build resources for CartoScope, joining the pixel-level results from FICTURE
    Key Parameters:
    - tsv-dir: Directory containing the TSV files of convert-sge output
    - fic-dir: Directory containing the FICTURE output
    - out-dir: Output directory
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Package SGE and FICTURE outputs into PMTiles and a catalog for CartoScope (ficture2 → joins → pmtiles → catalog)"
    )

    run_params = parser.add_argument_group("Run Options", "Execution controls for generating and running the Makefile")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate Makefile and print commands without executing')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_convert.mk", help='File name of Makefile to write (default: sge_convert.mk)')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job for tippecanoe')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output paths")
    inout_params.add_argument('--fic-dir', type=str, required=True, help='Directory containing FICTURE outputs (models, decodes, metadata)')
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output root directory (PMTiles, assets JSON, and catalog YAML)')
    inout_params.add_argument('--background-assets', type=str, nargs="+", help='One of more background asset in the formal of id:path or id1:id2:path')

    key_params = parser.add_argument_group("Key Parameters", "Metadata and logging")
    key_params.add_argument('--id', type=str, required=True, help='Identifier for the output assets; no whitespace; prefer "-" over "_"')
    key_params.add_argument('--title', type=str, help='Human-readable title for the output assets. Quote if contains spaces, e.g., "Example Title"')
    key_params.add_argument('--desc', type=str, help='Description of the output assets. Quote if contains spaces, e.g., "Example short description"')
    key_params.add_argument('--log', action='store_true', default=False, help='Write logs to a file under the output directory')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='Suffix for the log filename; final path is <out_dir>_cartload<suffix> (default: .log)')

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p4"')
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    env_params.add_argument('--tippecanoe', type=str, default=f"tippecanoe", help='Path to tippecanoe binary') # default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", 
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Advanced settings; defaults work for most cases")
    aux_params.add_argument('--in-fic-params', type=str, default="ficture.params.json", help='File name of YAML/JSON file with SGE paths and FICTURE parameters under --fic-dir (default: ficture.params.json)')
    aux_params.add_argument('--out-fic-assets', type=str, default="ficture_assets.json", help='File name of output JSON/YAML for FICTURE asset metadata under --out-dir (default: ficture_assets.json)')
    aux_params.add_argument('--out-catalog', type=str, default="catalog.yaml", help='File name of output catalog YAML under --out-dir (default: catalog.yaml)')
    aux_params.add_argument('--rename-x', type=str, default='x:lon', help='Column rename mapping for X axis in tippecanoe, format old:new (default: x:lon)')  
    aux_params.add_argument('--rename-y', type=str, default='y:lat', help='Column rename mapping for Y axis in tippecanoe, format old:new (default: y:lat)')  
    aux_params.add_argument('--colname-feature', type=str, default='gene', help='Column name for feature (default: gene)')
    aux_params.add_argument('--colname-count', type=str, default='count', help='Column name for UMI counts (default: count)')
    aux_params.add_argument('--out-molecules-id', type=str, default='genes', help='Base name for output molecules PMTiles files (no directory)')
    aux_params.add_argument('--max-join-dist-um', type=float, default=0.1, help='Max distance (in µm) to associate molecules with decoded pixels (default: 0.1)')
    aux_params.add_argument('--bin-count', type=int, default=50, help='Number of bins when splitting input molecules (default: 50)')
    aux_params.add_argument('--join-tile-size', type=float, default=500, help='Tile size (in µm) when joining molecules with decoded pixels (default: 500)')
    aux_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum tile size in bytes for tippecanoe/PMTiles (default: 5000000)')
    aux_params.add_argument('--max-feature-counts', type=int, default=500000, help='Maximum features per tile for tippecanoe/PMTiles (default: 500000)')
    aux_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Tippecanoe point-density preservation threshold (default: 1024)')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate files instead of cleaning up')
    aux_params.add_argument('--skip-raster', action='store_true', default=False, help='Skip raster image generation (no GDAL/go-pmtiles required)')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory (default: <out_dir>/tmp)')
    aux_params.add_argument('--transparent-below', type=int,  help='Set pixels below this value to transparent for dark background (optional)')
    aux_params.add_argument('--transparent-above', type=int,  help='Set pixels above this value to transparent for light background (optional)')
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def copy_rgb_tsv(in_rgb, out_rgb):
    with open(in_rgb, 'r') as f:
        hdrs = f.readline().rstrip().split("\t")
        col2idx = {hdr: i for i, hdr in enumerate(hdrs)}
        with open(out_rgb, 'w') as outf:
            outf.write("\t".join(["Name", "Color_index", "R", "G", "B"]) + "\n")
            for line in f:
                toks = line.rstrip().split("\t")
                if len(toks) != len(hdrs):
                    raise ValueError(f"Input RGB file {in_rgb} has inconsistent number of columns")
                rgb_r = int(toks[col2idx["R"]])/255
                rgb_g = int(toks[col2idx["G"]])/255
                rgb_b = int(toks[col2idx["B"]])/255
                name = toks[col2idx["Name"]]
                outf.write(f"{name}\t{name}\t{rgb_r:.5f}\t{rgb_g:.5f}\t{rgb_b:5f}\n")

def run_cartload2(_args):
    """
    Build resources for CartoScope
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_dir + "_cartload" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    # start mm
    mm = minimake()

    scheck_app(args.spatula)
    scheck_app(args.tippecanoe)
    if not args.skip_raster:
        scheck_app(args.pmtiles)
        scheck_app(args.gdal_translate)
        scheck_app(args.gdaladdo)
        #scheck_app(args.magick)

    # create output directory if needed
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(args.out_dir, "tmp")
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir, exist_ok=True)
    
    # output files/prefix
    out_assets_f = os.path.join(args.out_dir, args.out_fic_assets)
    out_catalog_f = os.path.join(args.out_dir,args.out_catalog)
    out_molecules_prefix=os.path.join(args.out_dir,args.out_molecules_id)

    # 2. Load FICTURE output metadata
    in_data = {}
    in_jsonf=f"{args.fic_dir}/{args.in_fic_params}"
    logger.info(f"Loading FICTURE metadata: {in_jsonf}")
    in_data = load_file_to_dict(in_jsonf)

    ## sge 
    in_sge = in_data.get("in_sge", {})
    in_molecules = in_sge.get("in_transcript", None)
    assert in_molecules is not None and os.path.exists(in_molecules), "Provide a valid input molecules file in the json file"
    in_features = in_sge.get("in_feature", None)
    assert in_features is not None and os.path.exists(in_features), "Provide a valid input features file in the json file"
    in_minmax = in_sge.get("in_minmax", None)
    assert in_minmax is not None and os.path.exists(in_minmax), "Provide a valid input minmax file in the json file"
    
    if not args.skip_raster:
        ## create raster mono pmtiles for SGE
        cmds = cmd_separator([], f"Converting SGE counts into PMTiles")
        cmd = " ".join([
            "cartloader", "run_tsv2mono",
            "--in-tsv", in_molecules,
            "--in-minmax", in_minmax,
            "--out-prefix", f"{args.out_dir}/sge-mono",
            "--colname-count", args.colname_count,
            f"--transparent-below {args.transparent_below}" if args.transparent_below else "",
            f"--transparent-above {args.transparent_above}" if args.transparent_above else "",
            "--main",
            f"--pmtiles '{args.pmtiles}'",
            f"--gdal_translate '{args.gdal_translate}'",
            f"--gdaladdo '{args.gdaladdo}'",
            f"--spatula '{args.spatula}'",
            "--keep-intermediate-files" if args.keep_intermediate_files else ""
        ])
        cmds.append(cmd)
        
        tsv2mono_flag = f"{args.out_dir}/sge-mono.done" # Update: Use a flag to make sure both light and dark pmtiles are done 
        cmds.append(f"[ -f {args.out_dir}/sge-mono-dark.pmtiles.done ] && [ -f {args.out_dir}/sge-mono-light.pmtiles.done ] && touch {tsv2mono_flag}")
        mm.add_target(tsv2mono_flag, [in_molecules, in_minmax], cmds)

    ## fic
    in_fic_params = in_data.get("train_params", [])
    #print(in_fic_params)
    if len(in_fic_params) == 0: ## parameters are empty
        logger.error(f"No training parameters found in {args.fic_dir}/{args.in_fic_params}")

    # create the output assets json
    out_fic_assets = ficture2_params_to_factor_assets(in_fic_params, args.skip_raster)

    train_targets = []

    join_pixel_tsvs = []
    join_pixel_ids = []

    minmax = read_minmax(in_minmax, "row")
    xmin = minmax["xmin"]
    xmax = minmax["xmax"]
    ymin = minmax["ymin"]
    ymax = minmax["ymax"]

    for train_param in in_fic_params:
        model_id = train_param["model_id"]
        # model_path = train_param.get("model_path", None)
        # fit_path = train_param.get("fit_path", None)
        # de_path = train_param.get("de_path", None)
        # info_path = train_param.get("info_path", None)
        factormap_path = train_param.get("factor_map", None)
        model_rgb = train_param["cmap"]

        in_prefix = f"{args.fic_dir}/{model_id}"
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
                "in": train_param.get("de_path",f"{in_prefix}.bulk_chisq.tsv"),
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

        cmds = cmd_separator([], f"Converting LDA-trained factors {model_id} into PMTiles and copying relevant files.")
        
        prerequisites = []
        outfiles=[]

        in_fit_tsvf = train_param.get("fit_path", f"{in_prefix}.results.tsv.gz")
        if os.path.exists(in_fit_tsvf):
            cmd = " ".join([
                "cartloader", "convert_generic_tsv_to_pmtiles",
                #"--in-tsv", out_fit_tsvf, 
                "--in-tsv", in_fit_tsvf,
                "--out-prefix", out_prefix,
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
                outfiles.append(val["out"])
                if key == "rgb":
                    copy_rgb_tsv(val["in"], val["out"])
                else:
                    cmds.append(f"cp {val['in']} {val['out']}")
        
        touch_flag_cmd=valid_and_touch_cmd(outfiles, f"{out_prefix}.done") # this only touch the flag file when all output files exist
        cmds.append(touch_flag_cmd)
        
        mm.add_target(f"{out_prefix}.done", prerequisites, cmds)

        ## add to the output asset json
        out_decode_assets = []

        sources = [f"{out_prefix}.done"]

        for decode_param in train_param["decode_params"]:
            in_id = decode_param["decode_id"]
            in_prefix = f"{args.fic_dir}/{in_id}"
            in_pixel_tsvf = decode_param.get("pixel_tsv_path",f"{in_prefix}.tsv.gz")
            in_pixel_png = decode_param.get("pixel_png_path",f"{in_prefix}.png")
            in_de_tsvf  = decode_param.get("de_tsv_path",f"{in_prefix}.bulk_chisq.tsv")
            in_post_tsvf = decode_param.get("pseudobulk_tsv_path",f"{in_prefix}.pseudobulk.tsv.gz")
            in_info_tsvf = decode_param.get("info_tsv_path",f"{in_prefix}.factor.info.tsv")

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

    ## Join pixel-level TSVs
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
        mm.add_target(f"{out_join_pixel_prefix}.tsv.gz", [in_molecules] + join_pixel_tsvs, cmds)
        molecules_f = f"{out_join_pixel_prefix}.tsv.gz"

    ## run tsv2pmtiles for the convert the joined pixel-level TSV to PMTiles
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

    cmds = cmd_separator([], f"Finalizing all conversions")
    cmds.append(f"touch {args.out_dir}/ficture.done")
    mm.add_target(f"{args.out_dir}/ficture.done", train_targets, cmds)

    ## Write the output asset json for FICTURE output
    logger.info("Writing FICTURE asset metadata JSON")
    write_dict_to_file(out_fic_assets, out_assets_f, check_equal=True)

    # 3. Create a yaml for all output assets
    cmds = cmd_separator([], f"Writing catalog YAML for all output assets")
    cmd = " ".join([
        "cartloader", "write_catalog_for_assets",
        "--sge-index", f"{out_molecules_prefix}_pmtiles_index.tsv", 
        "--sge-counts", f"{out_molecules_prefix}_bin_counts.json", 
        "--fic-assets", out_assets_f,
        "--out-catalog", f"{out_catalog_f}",
        f"--log --log-suffix '{args.log_suffix}'" if args.log else ""
    ] + (["--overview", f"sge-mono-dark.pmtiles",
          "--basemap", f"sge:dark:sge-mono-dark.pmtiles", 
          f"sge:light:sge-mono-light.pmtiles"] if not args.skip_raster else []
        ) + (args.background_assets if args.background_assets is not None else []))
      
    if ( args.id is not None ):
        cmd += f" --id {args.id}"

    if ( args.title is None ): ## if title is empty, use ID as title
        cmd += f" --title {args.id}"
    else:
        cmd += f" --title {args.title}"

    if ( args.desc is not None ):
        cmd += f" --desc {args.desc}"

    prerequisites_yaml=[f"{out_molecules_prefix}_pmtiles_index.tsv", out_assets_f]
    if not args.skip_raster:
        prerequisites_yaml.append(tsv2mono_flag)
    
    cmds.append(cmd)
    mm.add_target(f"{out_catalog_f}", prerequisites_yaml, cmds)

    if len(mm.targets) == 0:
        logger.error("No tasks were generated. Check inputs and parameters.")
        sys.exit(1)

    ## write makefile
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

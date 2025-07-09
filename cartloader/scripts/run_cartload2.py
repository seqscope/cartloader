import sys, os, argparse, logging, subprocess
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file, ficture2_params_to_factor_assets, read_minmax, flexopen

def parse_arguments(_args):
    """
    Build resources for CartoScope, joining the pixel-level results from FICTURE
    Key Parameters:
    - tsv-dir: Directory containing the TSV files of convert-sge output
    - fic-dir: Directory containing the FICTURE output
    - out-dir: Output directory
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_cartload2", description="Build resources for CartoScope, joining the pixel-level results from FICTURE")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default="run_cartload2.mk", help='The name of the Makefile to generate')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--fic-dir', type=str, required=True, help='Input directory containing FICTURE output')
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--id', type=str, required=True, help='The identifier of the output assets')
    key_params.add_argument('--title', type=str, help='The title of the output assets')
    key_params.add_argument('--desc', type=str, help='The description of output assets')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    # aux_params.add_argument('--magick', type=str, default=f"magick", help='Path to ImageMagick binary') # Disable this function. The user need to add the path to the ImageMagick binary directory to the PATH environment variable
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p4"')
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    env_params.add_argument('--tippecanoe', type=str, default=f"tippecanoe", help='Path to tippecanoe binary') # default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", 
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--in-fic-params', type=str, default="ficture.params.json", help='The YAML/JSON file containing both the SGE files and FICTURE parameters')
    aux_params.add_argument('--out-fic-assets', type=str, default="ficture_assets.json", help='The YAML/JSON file containing FICTURE output assets')
    aux_params.add_argument('--out-catalog', type=str, default="catalog.yaml", help='The YAML file containing the output catalog')
    aux_params.add_argument('--background-assets', type=str, nargs="+", help='The JSON/YAML file containing background assets in the form of [id:file] or [id1:id2:file]')
    aux_params.add_argument('--rename-x', type=str, default='x:lon', help='tippecanoe parameters to rename X axis')  
    aux_params.add_argument('--rename-y', type=str, default='y:lat', help='tippecanoe parameters to rename Y axis')  
    aux_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    aux_params.add_argument('--colname-count', type=str, default='count', help='Column name for feature counts')
    aux_params.add_argument('--out-molecules-id', type=str, default='genes', help='Prefix of output molecules PMTiles files. No directory path should be included')
    aux_params.add_argument('--max-join-dist-um', type=float, default=0.1, help='Maximum distance allowed to join molecules and pixel in micrometers')
    aux_params.add_argument('--join-tile-size', type=float, default=500, help='Tile size for joining molecules and pixel in micrometers')
    aux_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles')
    aux_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles')
    aux_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--skip-raster', action='store_true', default=False, help='Skip processing raster files, removing dependency to gdal, go-pmtiles')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory to be used (default: {out-dir}/tmp')
    aux_params.add_argument('--bin-count', type=int, default=50, help='Number of bins for splitting the input molecules')
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
    print(in_jsonf)
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
            "--main",
            f"--pmtiles '{args.pmtiles}'",
            f"--gdal_translate '{args.gdal_translate}'",
            f"--gdaladdo '{args.gdaladdo}'",
            f"--spatula '{args.spatula}'",
            "--keep-intermediate-files" if args.keep_intermediate_files else ""
        ])
        cmds.append(cmd)
        mm.add_target(f"{args.out_dir}/sge-mono-dark.pmtiles.done", [in_molecules, in_minmax], cmds)

    ## fic
    in_fic_params = in_data.get("train_params", [])
    #print(in_fic_params)
    if len(in_fic_params) == 0: ## parameters are empty
        logging.error(f"The parameters are empty after loading {args.fic_dir}/{args.in_fic_params}")

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
                "in":  f"{in_prefix}.model.tsv",
                "out": f"{out_prefix}-model.tsv"
            },
            "rgb":{
                "required": True,
                "in":  model_rgb,
                "out": f"{out_prefix}-rgb.tsv"
            },
            "de":{
                "required": False,
                "in":  f"{in_prefix}.bulk_chisq.tsv",
                "out": f"{out_prefix}-bulk-de.tsv"
            },
            "info":{
                "required": False,
                "in":  f"{in_prefix}.factor.info.tsv",
                "out": f"{out_prefix}-info.tsv"
            }
        }
        if factormap_path is not None:
            train_inout["factor_map"] = {
                "required": False,
                "in":  factormap_path,
                "out": f"{out_prefix}-factor-map.tsv"
            }

        prerequisites = []
        cmds = cmd_separator([], f"Converting LDA-trained factors {model_id} into PMTiles and copying relevant files..")

        # fit_results
        #print(f"--tippecanoe '{args.tippecanoe}'")
        in_fit_tsvf = f"{in_prefix}.results.tsv.gz"
        #out_fit_tsvf = f"{out_prefix}.results.tsv.gz"
        if os.path.exists(in_fit_tsvf):
            # cmd = " ".join([
            #     args.spatula, "append-topk-tsv",
            #     "--in-tsv", in_fit_tsvf,
            #     "--out-tsv", out_fit_tsvf,
            #     "--icol-beg", "2"
            # ])
            # cmds.append(cmd)
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
            #cmds.append(f"rm -f {out_fit_tsvf}")
            prerequisites.append(in_fit_tsvf)

        # mode/rgb/de/posterior/info
        for key, val in train_inout.items():
            if val["required"] or os.path.exists(val["in"]):
                prerequisites.append(val["in"])
                if key == "rgb":
                    copy_rgb_tsv(val["in"], val["out"])
                else:
                    cmds.append(f"cp {val['in']} {val['out']}")
        
        cmds.append(f"touch {out_prefix}.done")
        mm.add_target(f"{out_prefix}.done", prerequisites, cmds)

        ## add to the output asset json
        out_decode_assets = []

        sources = [f"{out_prefix}.done"]

        for decode_param in train_param["decode_params"]:
            in_id = decode_param["decode_id"]
            in_prefix = f"{args.fic_dir}/{in_id}"
            in_pixel_tsvf = f"{in_prefix}.tsv.gz"
            in_pixel_png = f"{in_prefix}.png"
            in_de_tsvf  = f"{in_prefix}.bulk_chisq.tsv"
            in_post_tsvf = f"{in_prefix}.pseudobulk.tsv"
            in_info_tsvf = f"{in_prefix}.factor.info.tsv"

            out_id = in_id.replace("_", "-")
            out_prefix = os.path.join(args.out_dir, out_id)

            join_pixel_tsvs.append(in_pixel_tsvf)
            join_pixel_ids.append(out_id)

            cmds = cmd_separator([], f"Converting decoded factors {in_id} into PMTiles and copying relevant files..")
            if not args.skip_raster:
                cmd = " ".join([
                    "cartloader", "run_fig2pmtiles",
                    "--georeference", "--geotif2mbtiles", "--mbtiles2pmtiles",
                    "--in-fig", in_pixel_png,
                    f"--in-bounds={xmin},{ymin},{xmax},{ymax}",
                    "--out-prefix", f"{out_prefix}-pixel-raster",
                    f"--pmtiles '{args.pmtiles}'",
                    f"--gdal_translate '{args.gdal_translate}'",
                    f"--gdaladdo '{args.gdaladdo}'",
                    "--keep-intermediate-files" if args.keep_intermediate_files else ""
                ])
                cmds.append(cmd)

            cmds.append(f"cp {in_de_tsvf} {out_prefix}-bulk-de.tsv")
            cmds.append(f"cp {in_post_tsvf} {out_prefix}-pseudobulk.tsv")
            cmds.append(f"cp {in_info_tsvf} {out_prefix}-info.tsv")
            cmds.append(f"touch {out_prefix}.done")
            mm.add_target(f"{out_prefix}.done", [in_pixel_tsvf, in_pixel_png, in_de_tsvf, in_post_tsvf, model_rgb, in_info_tsvf], cmds)
            sources.append(f"{out_prefix}.done")

        cmds = cmd_separator([], f"Finishing up for train parameters {out_train_id}")
        cmds.append(f"touch {out_train_prefix}.alldone")
        mm.add_target(f"{out_train_prefix}.alldone", sources, cmds)
        train_targets.append(f"{out_train_prefix}.alldone")

    ## Join pixel-level TSVs
    if ( len(join_pixel_tsvs) > 0 ):
        cmds = cmd_separator([], f"Pasting pixel-level TSVs")
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
        mm.add_target(f"{out_join_pixel_prefix}.tsv.gz", [in_molecules], cmds)

    ## run tsv2pmtiles for the convert the joined pixel-level TSV to PMTiles
    cmds = cmd_separator([], f"Converting the joined pixel-level TSV to PMTiles")
    cmd = " ".join([
        "cartloader", "run_tsv2pmtiles",
        "--in-molecules", f"{out_join_pixel_prefix}.tsv.gz",
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
    mm.add_target(f"{out_molecules_prefix}_pmtiles_index.tsv", [f"{out_join_pixel_prefix}.tsv.gz"], cmds)

    cmds = cmd_separator([], f"Finishing up for all conversions")
    cmds.append(f"touch {args.out_dir}/ficture.done")
    mm.add_target(f"{args.out_dir}/ficture.done", train_targets, cmds)

    ## Write the output asset json for FICTURE output
    logger.info("Writing a json file for expected FICTURE output assets")
    write_dict_to_file(out_fic_assets, out_assets_f, check_equal=True)

    # 3. Create a yaml for all output assets
    cmds = cmd_separator([], f"Writing a YAML file for all output assets")
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

    cmds.append(cmd)
    mm.add_target(f"{out_catalog_f}", [f"{out_molecules_prefix}_pmtiles_index.tsv", out_assets_f], cmds)

    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)

    ## write makefile
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    if args.dry_run:
        os.system(f"make -f {make_f} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs}")
    else:
        exe_cmd=f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}"
        result = subprocess.run(exe_cmd, shell=True)
        if result.returncode != 0:
            print(f"Error in executing: {exe_cmd}")
            sys.exit(1)

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

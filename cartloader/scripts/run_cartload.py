import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, json
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file

def parse_arguments(_args):
    """
    Build resources for CartoScope
    Key Parameters:
    - tsv-dir: Directory containing the TSV files of convert-sge output
    - fic-dir: Directory containing the FICTURE output
    - out-dir: Output directory
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_cartload", description="Build resources for CartoScope")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-molecules', type=str, help='Input Long Format TSV/CSV (possibly gzipped) file containing the X/Y coordinates and gene expression counts per spot')
    inout_params.add_argument('--in-features', type=str, help='Input TSV/CSV (possibly gzipped) file containing the gene name and total count for each gene')
    inout_params.add_argument('--fic-dir', type=str, required=True, help='Input directory containing FICTURE output')
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--id', type=str, required=True, help='The identifier of the output assets')
    key_params.add_argument('--title', type=str, help='The title of the output assets')
    key_params.add_argument('--desc', type=str, help='The description of output assets')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--rename-x', type=str, default='X:lon', help='tippecanoe parameters to rename X axis')  
    aux_params.add_argument('--rename-y', type=str, default='Y:lat', help='tippecanoe parameters to rename Y axis')  
    aux_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    aux_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    aux_params.add_argument('--out-molecules-prefix', type=str, default='genes', help='Prefix of output molecules PMTiles files')
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--in-fic-params', type=str, default="ficture.params.json", help='The YAML/JSON file containing FICTURE parameters')
    aux_params.add_argument('--out-fic-assets', type=str, default="ficture_assets.json", help='The YAML/JSON file containing FICTURE output assets')
    aux_params.add_argument('--out-catalog', type=str, default="catalog.yaml", help='The YAML file containing the output catalog')
    aux_params.add_argument('--background-assets', type=str, help='The JSON/YAML file containing background assets')
    aux_params.add_argument('--makefn', type=str, default="Makefile", help='The name of the Makefile to generate')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_cartload(_args):
    """
    Build resources for CartoScope
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_dir + "_cartload" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    # start mm
    mm = minimake()

    # create output directory if needed
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)
    
    # 1. Run tsv2pmtiles
    cmds = cmd_separator([], f"Running tsv2pmtiles TSV file")
    cmd = " ".join([
        "cartloader", "run_tsv2pmtiles",
        "--in-molecules", args.in_molecules,
        "--in-features", args.in_features,
        "--out-prefix", f"{args.out_dir}/{args.out_molecules_prefix}",
        "--all",
        "--n-jobs", str(args.n_jobs),
    ])
    cmds.append(cmd)
    mm.add_target(f"{args.out_dir}/{args.out_molecules_prefix}_bin_counts.json", [args.in_molecules, args.in_features], cmds)
    logger.info("Wrote commands for tsv2pmtiles")

    # 2. Load FICTURE output metadata
    # Assume that the parameters have 
    in_fic_params = {}
    out_fic_assets = {}
    in_fic_params = load_file_to_dict(f"{args.fic_dir}/{args.in_fic_params}")
    if len(in_fic_params) == 0: ## parameters are empty
        logging.error(f"The parameters are empty after loading {args.fic_dir}/{args.in_fic_params}")

    train_targets = []
    out_train_assets = []
    out_fic_assets["train_assets"] = out_train_assets
    for train_param in in_fic_params["train_params"]:
        train_width = train_param["train_width"]
        n_factor = train_param["n_factor"]

        in_prefix = f"{args.fic_dir}/nF{n_factor}.d_{train_width}"
        in_fit_tsvf = f"{in_prefix}.fit_result.tsv.gz"
        in_de_tsvf  = f"{in_prefix}.bulk_chisq.tsv"
        in_post_tsvf = f"{in_prefix}.posterior.count.tsv.gz"
        in_model_tsvf = f"{in_prefix}.model_matrix.tsv.gz"
        in_rgb_tsvf = f"{in_prefix}.rgb.tsv"
        in_info_tsvf = f"{in_prefix}.factor.info.tsv"
        out_prefix = f"factor-hexagon-nF{n_factor}-d{train_width}"
        
        cmds = cmd_separator([], f"Converting LDA-trained factors {in_prefix} into PMTiles and copying relevant files..")
        cmd = " ".join([
            "cartloader", "convert_generic_tsv_to_pmtiles",
            "--in-tsv", in_fit_tsvf, 
            "--out-prefix", f"{args.out_dir}/{out_prefix}",
            "--rename-column", "x:lon", "y:lat",
            "--colname-feature", args.colname_feature,
            "--colname-count", args.colname_count,
            "--log"
        ])
        cmds.append(cmd)
        cmds.append(f"cp {in_de_tsvf} {args.out_dir}/{out_prefix}-bulk-de.tsv")
        cmds.append(f"cp {in_post_tsvf} {args.out_dir}/{out_prefix}-posterior-counts.tsv.gz")
        cmds.append(f"cp {in_model_tsvf} {args.out_dir}/{out_prefix}-model-matrix.tsv.gz")
        cmds.append(f"cp {in_rgb_tsvf} {args.out_dir}/{out_prefix}-rgb.tsv")
        cmds.append(f"cp {in_info_tsvf} {args.out_dir}/{out_prefix}-info.tsv")
        mm.add_target(f"{args.out_dir}/{out_prefix}.pmtiles", [in_fit_tsvf, in_de_tsvf, in_post_tsvf, in_rgb_tsvf, in_info_tsvf], cmds)

        ## add to the output asset json
        out_proj_assets = []
        out_train_assets.append({
            "prefix": out_prefix,
            "pmtiles" : f"{out_prefix}.pmtiles",
            "post" : f"{out_prefix}-posterior-counts.tsv.gz",
            "model" : f"{out_prefix}-model-matrix.tsv.gz",
            "de" : f"{out_prefix}-bulk-de.tsv",
            "rgb" : f"{out_prefix}-rgb.tsv",
            "info" : f"{out_prefix}-info.tsv",
            "proj_assets": out_proj_assets
        })

        sources = [f"{args.out_dir}/{out_prefix}.pmtiles"]

        for proj_param in train_param["proj_params"]:
            fit_width = proj_param["fit_width"]
            anchor_res = proj_param["anchor_res"]
            in_prefix = f"{args.fic_dir}/nF{n_factor}.d_{train_width}.prj_{fit_width}.r_{anchor_res}"
            in_fit_tsvf = f"{in_prefix}.fit_result.tsv.gz"
            in_de_tsvf  = f"{in_prefix}.bulk_chisq.tsv"
            in_post_tsvf = f"{in_prefix}.posterior.count.tsv.gz"
            in_info_tsvf = f"{in_prefix}.factor.info.tsv"
            out_prefix = f"factor-hexagon-nF{n_factor}-d{train_width}-prj{fit_width}-r{anchor_res}"

            cmds = cmd_separator([], f"Converting projected factors {in_prefix} into PMTiles and copying relevant files..")
            cmd = " ".join([
                "cartloader", "convert_generic_tsv_to_pmtiles",
                "--in-tsv", in_fit_tsvf, 
                "--out-prefix", f"{args.out_dir}/{out_prefix}",
                "--rename-column", "x:lon", "y:lat",
                "--log"
            ])
            cmds.append(cmd)
            cmds.append(f"cp {in_de_tsvf} {args.out_dir}/{out_prefix}-bulk-de.tsv")
            cmds.append(f"cp {in_post_tsvf} {args.out_dir}/{out_prefix}-posterior-counts.tsv.gz")
            cmds.append(f"cp {in_info_tsvf} {args.out_dir}/{out_prefix}-info.tsv")
            mm.add_target(f"{args.out_dir}/{out_prefix}.pmtiles", [in_fit_tsvf, in_de_tsvf, in_post_tsvf, in_rgb_tsvf, in_info_tsvf], cmds)
            sources.append(f"{args.out_dir}/{out_prefix}.pmtiles")

            out_decode_assets = []
            out_proj_assets.append({
                "prefix": out_prefix,
                "pmtiles" : f"{out_prefix}.pmtiles",
                "post" : f"{out_prefix}-posterior-counts.tsv.gz",
                "de" : f"{out_prefix}-bulk-de.tsv",
                "info" : f"{out_prefix}-info.tsv",
                "decode_assets": out_decode_assets
            })


            for decode_param in proj_param["decode_params"]:
                radius = decode_param["radius"]
                in_prefix = f"{args.fic_dir}/nF{n_factor}.d_{train_width}.decode.prj_{fit_width}.r_{anchor_res}_{radius}"
                in_pixel_tsvf = f"{in_prefix}.pixel.sorted.tsv.gz"
                in_de_tsvf = f"{in_prefix}.bulk_chisq.tsv"
                in_post_tsvf = f"{in_prefix}.posterior.count.tsv.gz"
                in_info_tsvf = f"{in_prefix}.factor.info.tsv"
                out_prefix = f"factor-pixel-nF{n_factor}-d{train_width}-prj{fit_width}-r{anchor_res}-{radius}"

                cmds = cmd_separator([], f"Converting pixel-level factors {in_prefix} into PMTiles and copying relevant files..")
                cmd = " ".join([
                    "cartloader", "convert_ficture_pixel_to_pmtiles",
                    "--in-pixel", in_pixel_tsvf, 
                    "--out-prefix", f"{args.out_dir}/{out_prefix}",
                    "--log"
                ])
                cmds.append(cmd)
                cmds.append(f"cp {in_de_tsvf} {args.out_dir}/{out_prefix}-bulk-de.tsv")
                cmds.append(f"cp {in_post_tsvf} {args.out_dir}/{out_prefix}-posterior-counts.tsv.gz")
                cmds.append(f"cp {in_info_tsvf} {args.out_dir}/{out_prefix}-info.tsv")
                mm.add_target(f"{args.out_dir}/{out_prefix}.pmtiles", [in_fit_tsvf, in_de_tsvf, in_post_tsvf, in_rgb_tsvf, in_info_tsvf], cmds)
                sources.append(f"{args.out_dir}/{out_prefix}.pmtiles")

                out_decode_assets.append({
                    "prefix": out_prefix,
                    "pmtiles" : f"{out_prefix}.pmtiles",
                    "post" : f"{out_prefix}-posterior-counts.tsv.gz",
                    "de" : f"{out_prefix}-bulk-de.tsv",
                    "info" : f"{out_prefix}-info.tsv"
                })


        out_prefix = f"factor-hexagon-nF{n_factor}-d{train_width}"
        cmds = cmd_separator([], f"Finishing up for train parameters {out_prefix}")
        cmds.append(f"touch {args.out_dir}/{out_prefix}.done")
        mm.add_target(f"{args.out_dir}/{out_prefix}.done", sources, cmds)
        train_targets.append(f"{args.out_dir}/{out_prefix}.done")

    cmds = cmd_separator([], f"Finishing up for all train parameters")
    cmds.append(f"touch {args.out_dir}/ficture.done")
    mm.add_target(f"{args.out_dir}/ficture.done", train_targets, cmds)

    ## Write the output asset json for FICTURE output
    logger.info("Writing a json file for expected FICTURE output assets")

    write_dict_to_file(out_fic_assets, f"{args.out_dir}/{args.out_fic_assets}")

    # 3. Create output asset json for all output assets
    cmds = cmd_separator([], f"Writing a YAML file for all output assets")
    cmd = " ".join([
        "cartloader", "write_yaml_for_assets",
        "--sge-index", f"{args.out_dir}/{args.out_molecules_prefix}_pmtiles_index.tsv", 
        "--sge-counts", f"{args.out_dir}/{args.out_molecules_prefix}_bin_counts.json", 
        "--fic-assets", f"{args.out_dir}/{args.out_fic_assets}",
        "--out-catalog", f"{args.out_dir}/{args.out_catalog}",
        "--log"
    ])
    if ( args.background_assets is not None ):
        cmd += f" --background-assets {args.background_assets}"
    if ( args.id is not None ):
        cmd += f" --id {args.id}"

    if ( args.title is None ): ## if title is empty, use ID as title
        cmd += f" --title {args.id}"
    else:
        cmd += f" --title {args.title}"

    if ( args.desc is not None ):
        cmd += f" --desc {args.desc}"

    cmds.append(cmd)
    mm.add_target(f"{args.out_dir}/{args.out_catalog}", [f"{args.out_dir}/{args.out_molecules_prefix}_bin_counts.json", f"{args.out_dir}/{args.out_fic_assets}"], cmds)

    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)

    ## write makefile
    mm.write_makefile(f"{args.out_dir}/{args.makefn}")

    ## run makefile
    if args.dry_run:
        ## run makefile
        os.system(f"make -f {args.out_dir}/{args.makefn} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")
    else:
        os.system(f"make -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, json, inspect
import pandas as pd
import yaml

from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file, update_and_copy_paths

def parse_arguments(_args):
    """
    Write a new YAML file or update an existing YAML file for all output assets
    """
    
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Write YAML for all output assets")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directories and files")
    inout_params.add_argument('--out-catalog', type=str, required=True, help='JSON/YAML file containing the output assets')
    inout_params.add_argument('--write-mode', type=str, default="append", choices=["append", "write"], help='Write mode for the output catalog. Default: append. If write, a new file will be created based on the arguments provided. If append, the new assets will be added or updated in the existing catalog file.')
    inout_params.add_argument('--sge-index', type=str, help='Path to an index TSV file containing the SGE output converted to PMTiles.')
    inout_params.add_argument('--sge-counts', type=str, help='Path to a JSON file containing SGE counts per gene')
    inout_params.add_argument('--fic-assets', type=str,  help='Path to a JSON/YAML file containing FICTURE output assets')
    inout_params.add_argument('--background-assets', type=str, nargs="+", default=[], help='Path(s) of one or more JSON/YAML file(s) containing background assets, if present')
    inout_params.add_argument('--cell-assets', type=str, nargs="+", default=[], help='Path(s) to one or more JSON/YAML file(s) containing cell segmentation assets, if present')
    inout_params.add_argument('--basemap', type=str, nargs="+", default=[], help='One or more basemap assets. Each must be in the format "id:filename" or "id1:id2:filename", where each ID represents a basemap identifier.')    
    inout_params.add_argument('--basemap-dir', type=str, default=None, help='Directory containing the basemap files. By default, the directory of the out-catalog file is used as the basemap directory.')
    inout_params.add_argument('--overview', type=str, help='Specify one of these basemaps as the overview asset, using its filename.')
    inout_params.add_argument('--check-equal', action="store_true", default=False, help="If enabled, the script checks for an existing file with matching content and skips writing the JSON/YAML unless differences are found.")

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--id', type=str, help='The identifier of the output assets')
    key_params.add_argument('--title', type=str, help='The title of the output assets')
    key_params.add_argument('--desc', type=str, help='The description of output assets')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def write_catalog_for_assets(_args):
    """
    Build resources for CartoScope
    """

    # parse argument
    args=parse_arguments(_args)
    out_dir = os.path.dirname(args.out_catalog)
    os.makedirs(out_dir, exist_ok=True)

    logger = create_custom_logger(__name__, args.out_catalog + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    flags = []

    ## if the output catalog file exists, check the write mode
    if os.path.exists(args.out_catalog) and args.write_mode == "append":
        logger.info(f"Found an existing catalog file {args.out_catalog}. This file will be appended.")
        with open(args.out_catalog, 'r') as stream:
            catalog_dict = yaml.load(stream, Loader=yaml.FullLoader)
    else:
        catalog_dict = {}

    ## 1st - level keys
    if ( args.id is not None ):
        catalog_dict["id"] = args.id
    if ( args.title is not None ):
        catalog_dict["title"] = args.title
    if ( args.desc is not None ):
        catalog_dict["description"] = args.desc
    if "assets" not in catalog_dict:
        catalog_dict["assets"] = {}
    
    # 2nd - level keys
    # * assets
    if args.sge_index is not None or args.sge_counts is not None:
        logger.info(f"Reading the SGE information (sge index file: {args.sge_index}; sge counts file: {args.sge_counts})")

        if "sge" not in catalog_dict["assets"]:
            sge_dict={}
        else:
            sge_dict = catalog_dict["assets"]["sge"]
        ## read the SGE index file
        if args.sge_index is not None:
            sge_df = pd.read_csv(args.sge_index, sep="\t")
            all_pmtiles = sge_df.loc[sge_df['bin_id'] == 'all']['pmtiles_path'].to_list()[0]
            bin_pmtiles = sge_df.loc[sge_df['bin_id'] != 'all']['pmtiles_path'].to_list()
            sge_dict["all"] = all_pmtiles
            sge_dict["bins"] = bin_pmtiles
        if args.sge_counts is not None:
            sge_dict["counts"] = os.path.basename(args.sge_counts)
        catalog_dict["assets"]["sge"] = sge_dict

    # - overview
    if args.overview is not None:
        logger.info(f"Updating the overview parameters {args.overview}")
        catalog_dict["assets"]["overview"]=args.overview

    # - factors
    if args.fic_assets is not None or len(args.cell_assets)>0:
        logger.info(f"Reading the FICTURE assets params {args.fic_assets}")
        if "factors" not in catalog_dict["assets"]:
            factors_list=[]
        else:
            factors_list = catalog_dict["assets"]["factors"]
        
        if args.fic_assets is not None: 
            fic_assets = load_file_to_dict(args.fic_assets)
            factors_list.extend(fic_assets)

        if len(args.cell_assets)>0:
            for cell_assets_f in args.cell_assets:
                cell_assets=load_file_to_dict(cell_assets_f)
                # update_and_copy_paths will update the path to be only filename, which fits the needs of catalog.yaml
                cell_assets=update_and_copy_paths(cell_assets, out_dir, skip_keys=["id", "name", "cells_id"], exe_copy=True)
                factors_list.append(cell_assets)
                flags.append(f"{cell_assets_f}.done")

        ## add files to the catalog
        unique_factors = {item['id']: item for item in factors_list}         # Use a dictionary to keep the latest entry by 'id'
        factors_list = list(unique_factors.values())                         # Get the deduplicated list
        catalog_dict["assets"]["factors"]=factors_list

    # - basemap
    #   * load basemap dict
    if len(args.background_assets)>0 or len(args.basemap)>0 :
        # define basemap_dict
        if "basemap" not in catalog_dict["assets"]:
            basemap_dict = {}
        else:
            basemap_dict = catalog_dict["assets"]["basemap"]
        
        if ( len(args.basemap)>0 ):
            logger.info(f"Updating the basemap parameters {args.basemap}")

            if args.basemap_dir is None:
                args.basemap_dir = os.path.dirname(args.out_catalog)

            for basemap in args.basemap:
                logger.info(f"  - Adding the basemap {basemap}")
                toks = basemap.split(":")
                if ( len(toks) == 2 ):
                    (basemap_type, basemap_fn) = toks
                    basemap_dict[basemap_type] = basemap_fn
                elif ( len(toks) == 3 ):
                    (basemap_type, basemap_id, basemap_fn) = toks
                    if ( basemap_type not in basemap_dict):
                        basemap_dict[basemap_type] = {}
                        basemap_dict[basemap_type]["default"] = basemap_id
                    basemap_dict[basemap_type][basemap_id] = basemap_fn
                else:
                    logger.error(f"Invalid basemap format {basemap}")
                    sys.exit(1)
                if basemap_type != "sge":
                    flags.append(os.path.join(args.basemap_dir, f"{basemap_fn}.yaml.done"))
    
        if len(args.background_assets)>0:
            for background_assets_f in args.background_assets:
                ## each backgroudn_assets file should contains only one key:value pair 
                background_assets = load_file_to_dict(background_assets_f)
                background_assets = update_and_copy_paths(background_assets, out_dir, skip_keys=[], exe_copy=True)
                basemap_dict.update(background_assets)
                flags.append(f"{background_assets_f}.done")
            
        catalog_dict["assets"]["basemap"] = basemap_dict

    ## write the output catalog
    logger.info(f"Writing the output catalog file {args.out_catalog}")
    write_dict_to_file(catalog_dict, args.out_catalog, check_equal=args.check_equal)

    if ( len(flags) > 0 ):
        for flag_file in flags:
            subprocess.run(["touch", flag_file])

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

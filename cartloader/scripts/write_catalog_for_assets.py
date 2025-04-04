import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, json
import pandas as pd
import yaml

from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file

def parse_arguments(_args):
    """
    Write a new YAML file or update an existing YAML file for all output assets
    Key Parameters:
    - sge-index: TSV file containing the index of the SGE output converted to PMTiles
    - fic-assets: JSON/YAML file containing the FICTURE output assets
    - background-assets: JSON/YAML file containing the background assets
    - out: Output YAML file
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader write_yaml_for_assets", description="Write YAML for all output assets")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--out-catalog', type=str, required=True, help='JSON/YAML file containing the output assets')
    inout_params.add_argument('--write-mode', type=str, default="append", choices=["append", "write"], help='Write mode for the output catalog. Default: append. If write, a new file will be created based on the arguments provided. If append, the new assets will be added or updated in the existing catalog file.')
    inout_params.add_argument('--sge-index', type=str, help='Index TSV file containing the SGE output converted to PMTiles.')
    inout_params.add_argument('--sge-counts', type=str, help='JSON file containing SGE counts per gene')
    inout_params.add_argument('--fic-assets', type=str,  help='JSON/YAML file containing FICTURE output assets')
#    inout_params.add_argument('--background-assets', type=str, help='JSON/YAML file containing background assets, if exists')
    inout_params.add_argument('--overview', type=str, help='File containing the overview assets')
    inout_params.add_argument('--basemap', type=str, nargs="+", default=[], help='[type:filename] or [type:id:filename] containing the basemap assets, where the id is the identifier of the basemap, for example hist_id in the yaml file.')    
    inout_params.add_argument('--basemap-dir', type=str, default=None, help='Directory containing the basemap files. By default, use the directory of the out-catalog file')

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
    if args.fic_assets is not None:
        logger.info(f"Reading the FICTURE assets params {args.fic_assets}")
        if "factors" not in catalog_dict["assets"]:
            factors_list=[]
        else:
            factors_list = catalog_dict["assets"]["factors"]
        ## load json of FICTURE assets
        fic_assets = load_file_to_dict(args.fic_assets)
        # updat factors_list by fic_assets 
        factors_combined = factors_list + fic_assets
        unique_factors = {item['id']: item for item in factors_combined}         # Use a dictionary to keep the latest entry by 'id'
        factors_combined = list(unique_factors.values())                                # Get the deduplicated list
        catalog_dict["assets"]["factors"]=factors_combined
        #print(f"Factors list: {factors_combined}")

    # - basemap
    if ( args.basemap is not [] ):
        logger.info(f"Updating the basemap parameters {args.basemap}")
        # update directory
        if args.basemap_dir is None:
            args.basemap_dir = os.path.dirname(args.out_catalog)

        # update basemap
        if "basemap" not in catalog_dict["assets"]:
            basemap_dict = {}
        else:
            basemap_dict = catalog_dict["assets"]["basemap"]
        
        basemap_flags =[]
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
                basemap_flags.append(os.path.join(args.basemap_dir, f"{basemap_fn}.yaml.done"))
        catalog_dict["assets"]["basemap"] = basemap_dict

    ## write the output catalog
    logger.info(f"Writing the output catalog file {args.out_catalog}")
    write_dict_to_file(catalog_dict, args.out_catalog)

    if ( len(basemap_flags) > 0 ):
        for flag_file in basemap_flags:
            subprocess.run(["touch", flag_file])

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
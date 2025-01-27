import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, json
import pandas as pd

from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file

def parse_arguments(_args):
    """
    Write a YAML file for all output assets
    Key Parameters:
    - sge-index: TSV file containing the index of the SGE output converted to PMTiles
    - fic-assets: JSON/YAML file containing the FICTURE output assets
    - background-assets: JSON/YAML file containing the background assets
    - out: Output YAML file
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader write_yaml_for_assets", description="Write YAML for all output assets")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--sge-index', type=str, required=True, help='Index TSV file containing the SGE output converted to PMTiles')
    inout_params.add_argument('--sge-counts', type=str, required=True, help='JSON file containing SGE counts per gene')
    inout_params.add_argument('--fic-assets', type=str, required=True, help='JSON/YAML file containing FICTURE output assets')
#    inout_params.add_argument('--background-assets', type=str, help='JSON/YAML file containing background assets, if exists')
    inout_params.add_argument('--overview', type=str, help='File containing the overview assets')
    inout_params.add_argument('--basemap', type=str, nargs="+", help='[type:filename] or [type:id:filename] containing the basemap assets, where the id is the identifier of the basemap, for example hist_id in the yaml file.')    
    inout_params.add_argument('--basemap-dir', type=str, default=None, help='Directory containing the basemap files. By default, use the directory of the out-catalog file')
    inout_params.add_argument('--out-catalog', type=str, required=True, help='JSON/YAML file containing the output assets')

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

    logger = create_custom_logger(__name__, args.out_catalog + args.log_suffix if args.log else None)
    logger.info("Analysis Started")


    ## read the SGE index file
    logger.info(f"Reading the SGE index file {args.sge_index}")

    sge_df = pd.read_csv(args.sge_index, sep="\t")
    all_pmtiles = sge_df.loc[sge_df['bin_id'] == 'all']['pmtiles_path'].to_list()[0]
    bin_pmtiles = sge_df.loc[sge_df['bin_id'] != 'all']['pmtiles_path'].to_list()

    sge_dict = {
        "all": all_pmtiles,
        "bins": bin_pmtiles,
        "counts": os.path.basename(args.sge_counts)
    }

    ## load json of FICTURE assets
    logger.info(f"Reading the FICTURE assets file {args.fic_assets}")
    fic_assets = load_file_to_dict(args.fic_assets)

    ## load json of background assets
    # background_assets = None

    ## create output directory
    out_dict = {}
    if ( args.id is not None ):
        out_dict["id"] = args.id
    if ( args.title is not None ):
        out_dict["title"] = args.title
    if ( args.desc is not None ):
        out_dict["description"] = args.desc

    out_dict["assets"] = {
        "sge": sge_dict,
        "factors": fic_assets
    }

    if ( args.overview is not None ):
        out_dict["assets"]["overview"] = args.overview
    
    if args.basemap_dir is None:
        args.basemap_dir = os.path.dirname(args.out_catalog)

    basemap_flags =[]
    if ( args.basemap is not None ):
        basemap_dict = {}
        for basemap in args.basemap:
            toks = basemap.split(":")
            if ( len(toks) == 2 ):
                (basemap_type, basemap_fn) = toks
                if ( basemap_type in basemap_dict ):
                    logger.error(f"Duplicate basemap id {basemap_type}")
                    sys.exit(1)
                basemap_dict[basemap_type] = basemap_fn
            elif ( len(toks) == 3 ):
                (basemap_type, basemap_id, basemap_fn) = toks
                if ( basemap_type not in basemap_dict ):
                    basemap_dict[basemap_type] = {}
                    basemap_dict[basemap_type]["default"] = basemap_id
                if ( basemap_id in basemap_dict[basemap_type] ):
                    logger.error(f"Duplicate basemap 'type:id' : {basemap_type}:{basemap_id}")
                    sys.exit(1)
                basemap_dict[basemap_type][basemap_id] = basemap_fn
            else:
                logger.error(f"Invalid basemap format {basemap}")
                sys.exit(1)
            if basemap_type != "sge":
                basemap_flags.append(os.path.join(args.basemap_dir, f"{basemap_fn}.yaml.done"))
        out_dict["assets"]["basemap"] = basemap_dict


    # if ( args.background_assets is not None ):
    #     out_dict["assets"]["background"] = args.background_assets

    ## write the output catalog
    logger.info(f"Writing the output catalog file {args.out_catalog}")
    write_dict_to_file(out_dict, args.out_catalog)

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
# !!! This script is not yet finished & implemented !!

import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, yaml
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Convert SGE into a tsv format.")

    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    key_params.add_argument('--yaml-file', default=None, type=str, help='If missing, the yaml file will be named as catalog.yaml and saved in the output directory')
    key_params.add_argument('--lda-ids', nargs='*', type=str, help='Training IDs for LDA model')
    key_params.add_argument('--decode-ids', nargs='*', type=str, help='Projection IDs for pixel-level decoding')
    key_params.add_argument('--all-ids', type=str, help='Write down all results into the yaml file')

    desc_params = parser.add_argument_group("Description Parameters", "Description parameters for the input data")
    desc_params.add_argument('--stac-version', type=str, default="1.0.0", help='STAC version')
    desc_params.add_argument('--id', type=str, help='ID of the data')
    desc_params.add_argument('--description', type=str, help='Description of the data')
    desc_params.add_argument('--title', type=str, help='Title of the data')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def ini_yaml(_args):
    args=parse_arguments(_args)
    yaml_content = {
        "stac_version": args.stac_version,
        "id": args.id,
        "description": args.description,
        "title": args.title,
        "assets":{
            "factor": {
                "DE": [],
                "hexagon": [],
                "pixel": [],
            },
            "segment": {
                "fullDepth": None,
                "perBin": None,
                "statBin": None,

            },
            "basemap": []
        }
    }
    return yaml_content

def write_yaml_for_zenodo(_args):
    args = parse_arguments(_args)
    # in/out
    assert os.path.exists(args.out_dir), "Provide a valid output directory"
    if args.yaml_file is None:
        yaml_file = os.path.join(args.out_dir, "catalog.yaml")
    else:
        os.makedirs(os.path.dirname(args.yaml_file), exist_ok=True)
        yaml_file = args.yaml_file
    # yaml content
    yaml_content = ini_yaml(args)
    de_list = []
    hex_df = []
    pixel_df = []
    #for lda_id in args.lda_ids:



if __name__ == "__main__":
    # get the cartloader path
    global cartloader_repo
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
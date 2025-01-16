import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, json, yaml
import pandas as pd

from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, load_file_to_dict, write_dict_to_file

def parse_arguments(_args):
    """
    Update the YAML file to add basemap in the assets
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader update_yaml_for_basemap", description="Update the yaml file to include the pmtiles to basemap in the assets")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-yaml', type=str, required=True, help='Input JSON/YAML file containing the assets')
    inout_params.add_argument('--out-yaml', type=str,  default=None, help='Output JSON/YAML file containing the assets. If not provided, the input file will be overwritten')
    # add a list of basemaps pmtiles
    inout_params.add_argument('--pmtiles', type=str, nargs='+', required=True, help='List of pmtiles to add to the basemap')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)


def update_yaml_for_basemap(_args):

    # parse argument
    args = parse_arguments(_args)

    # Load the catalog YAML file while preserving order
    with open(args.in_yaml, 'r') as stream:
        catalog = yaml.load(stream, Loader=yaml.FullLoader)  # Preserves order

    # Add the pmtiles to the basemap
    basemap_dict = catalog["assets"].get("basemap", {})
    for pmtile in args.pmtiles:
        toks = pmtile.split(":")
        if len(toks) == 2:
            (basemap_id, basemap_file) = toks
            if ( basemap_id in basemap_dict ):
                raise ValueError(f"Duplicate basemap id {basemap_id}")
            basemap_dict[basemap_id] = basemap_file
        elif ( len(toks) == 3 ):
            (basemap_id1, basemap_id2, basemap_file) = toks
            if ( basemap_id1 not in basemap_dict ):
                basemap_dict[basemap_id1] = {}
                basemap_dict[basemap_id1]["default"] = basemap_id2
            if ( basemap_id2 in basemap_dict[basemap_id1] ):
                raise ValueError (f"Duplicate basemap id {basemap_id1}:{basemap_id2}")
                sys.exit(1)
            basemap_dict[basemap_id1][basemap_id2] = basemap_file
        else:
            raise ValueError(f"Invalid basemap format {pmtile}")

    #pmtile_bn = [os.path.basename(pmtile) for pmtile in args.pmtiles]
    #basemap_list.extend(pmtile_bn)
    catalog["assets"]["basemap"] = basemap_dict

    # Write the updated catalog YAML file, preserving order
    if args.out_yaml is None:
        args.out_yaml = args.in_yaml
    with open(args.out_yaml, 'w') as stream:
        yaml.dump(catalog, stream, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
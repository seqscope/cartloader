# UPDATED: 2025-04-03
# USE write_catalog_for_assets instead.


import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, json, yaml
import pandas as pd

def parse_arguments(_args):
    """
    Update the YAML file to add basemap in the assets
    """
    #repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader update_yaml_for_basemap", description="Update the yaml file to include the basemaps to basemap in the assets")

    parser.add_argument('--in-yaml', type=str, required=True, help='Input JSON/YAML file containing the assets')
    parser.add_argument('--out-yaml', type=str,  default=None, help='Output JSON/YAML file containing the assets. If not provided, the input file will be overwritten')
    # add a list of basemaps basemaps
    parser.add_argument('--basemap', type=str, nargs="+", help='[type:filename] or [type:id:filename] containing the basemap assets, where the id is the identifier of the basemap, for example hist_id in the yaml file.')    
    parser.add_argument('--basemap-dir', type=str, default=None, help='Directory containing the basemap files. By default, use the directory of the in-yaml file')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite the existing basemap if it exists')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)


def update_catalog_for_basemap(_args):

    # parse argument
    args = parse_arguments(_args)

    # Load the catalog YAML file while preserving order
    with open(args.in_yaml, 'r') as stream:
        catalog = yaml.load(stream, Loader=yaml.FullLoader)  # Preserves order

    if args.basemap_dir is None:
        args.basemap_dir = os.path.dirname(args.in_yaml)

    # Add the basemaps to the basemap
    # Recommend to use toks==3, given it may have more than one histology file
    basemap_dict = catalog["assets"].get("basemap", {})
    flag_files=[]
    for basemap in args.basemap:
        toks = basemap.split(":")
        if len(toks) == 2:
            (basemap_type, basemap_fn) = toks
            basemap_dict[basemap_type] = basemap_fn
        elif ( len(toks) == 3 ):
            (basemap_id1, basemap_id2, basemap_fn) = toks
            if ( basemap_id1 not in basemap_dict ):
                basemap_dict[basemap_id1] = {}
                basemap_dict[basemap_id1]["default"] = basemap_id2
            basemap_dict[basemap_id1][basemap_id2] = basemap_fn
        else:
            raise ValueError(f"Invalid basemap format {basemap}")
        
        flag_files.append(os.path.join(args.basemap_dir, f"{basemap_fn}.yaml.done"))

    catalog["assets"]["basemap"] = basemap_dict

    # Write the updated catalog YAML file, preserving order
    if args.out_yaml is None:
        args.out_yaml = args.in_yaml
    with open(args.out_yaml, 'w') as stream:
        yaml.dump(catalog, stream, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)

    # touch empty done files
    for flag_file in flag_files:
        subprocess.run(["touch", flag_file])

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
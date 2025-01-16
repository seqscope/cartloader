import os
import subprocess
import argparse
import sys
import yaml


from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger

def run_command(command, dry_run):
    if dry_run:
        print("Dry-run command:", " ".join(command))
    else:
        subprocess.run(command, check=True)
def traverse_dict(d, parent_key=''):
    pairs = []
    for key, value in d.items():
        new_key = f"{parent_key}.{key}" if parent_key else key
        if isinstance(value, dict):
            pairs.extend(traverse_dict(value, new_key))  # Recurse into dictionaries and extend pairs list
        elif isinstance(value, list):
            for i, item in enumerate(value):
                pairs.extend(traverse_dict({f"item_{i}": item}, new_key))  # Recurse into list items
        else:
            pairs.append((new_key, value))  # Append key-value pair
    return pairs

def upload_aws_by_catalog(_args):
    parser = argparse.ArgumentParser(description="Upload files to S3 as specified in catalog.yaml. All files must be in the input directory.")
    parser.add_argument("--in-dir", help="Path to the input directory (e.g., /<in_dir>/cartload/).")
    parser.add_argument("--s3-dir", help="S3 directory path (e.g., s3://<bucket_name>/<ID>).")
    parser.add_argument("--catalog-yaml", default=None, help="Path to the catalog.yaml file (default: /<in_dir>/catalog.yaml).")
    parser.add_argument("--histology", nargs="*", default=None, help="List of histology pmtiles to upload (default: None).")
    parser.add_argument("--dry-run", action="store_true", help="Show commands without executing them (default: False).")

    args = parser.parse_args(_args)

    if args.catalog_yaml is None:
        catalog_f= os.path.join(args.in_dir, "catalog.yaml")
    else:
        catalog_f = args.catalog_yaml

    commands = []
    s3_catalog_f = f"{args.s3_dir}/catalog.yaml"

    # start mm
    mm = minimake()

    # Step 1: By default, the histology should already registered in the catalog.yaml. 
    # Added this --histology option to faciliate additional histology files.
    if args.histology is not None:
        for hist_path in args.histology:
            hist_bn = os.path.basename(hist_path)
            hist_copy = os.path.join(args.in_dir, hist_bn)

            # Add commands for uploading pmtiles and updating catalog.yaml
            cmds = cmd_separator([], f"Registering histology files to catalogy yaml...'")
            commands.append(f"cp {hist_path} {hist_copy}")
            commands.append([
                "cartloader", "update_yaml_for_basemap",
                "--catalog-yaml", catalog_f,
                "--pmtiles", hist_path
            ])

    commands.append(f"echo '============================================================================'")
    commands.append(f"echo 'Uploading cartload files to AWS...'")
    commands.append(f"echo '============================================================================'")

    # Option 1: use "." to locate the files
    # with open(catalog_f, "r") as catalog_file:
    #     for line in catalog_file:
    #         if "." in line:
    #             # Extract the file name
    #             file_name = line.strip().split()[-1]
    #             file_path = os.path.join(args.in_dir, file_name)
    #             s3_file_path = f"{args.s3_dir}/{file_name}"

    #             commands.append(f"echo 'Uploading {file_name} to S3...'")
    #             commands.append(["aws", "s3", "cp", file_path, s3_file_path])

    # Option 2: use traverse_dict to locate the files.
    # Note, if catalog.yaml file structure changes, this keys_id may need to be updated.
    keys_id=["id", "title", "name", "model_id", "proj_id", "decode_id"]
    with open(catalog_f, "r") as catalog_file:
        catalog = yaml.safe_load(catalog_file)
        catalog_pair=traverse_dict(catalog)
        for key, value in catalog_pair:
            subkey=key.split(".")[-1] if "." in key else key
            if subkey not in keys_id and key != "assets.basemap.sge.default":
                # now check if it exists in the input directory
                file_path=os.path.join(args.in_dir, value)
                if os.path.exists(file_path):
                    s3_file_path = f"{args.s3_dir}/{value}"
                    commands.append(["aws", "s3", "cp", file_path, s3_file_path])
                else:
                    raise ValueError (f"The file {value} does not exist in the input directory. Full path: {file_path}")

    # Step 3: Upload catalog.yaml to S3
    commands.append(f"echo '============================================================================'")
    commands.append(f"echo 'Uploading catalog to AWS...'")
    commands.append(f"echo '============================================================================'")
    commands.append(["aws", "s3", "cp", catalog_f, s3_catalog_f])

    commands.append(f"echo 'S3 catalog.yaml: {s3_catalog_f}'")
    # Execute or print all commands
    for cmd in commands:
        if isinstance(cmd, str):
            print(cmd)  # Always print echo statements
        else:
            run_command(cmd, args.dry_run)

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
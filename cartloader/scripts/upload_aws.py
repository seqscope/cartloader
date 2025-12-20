import os
import subprocess
import argparse
import sys
import yaml
from typing import Dict, List, Tuple
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, execute_makefile

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

def collect_files_from_yaml(catalog_f):
    print(f"Extracting input files from catalog file: {catalog_f}")
    with open(catalog_f, "r") as catalog_file:
        catalog = yaml.safe_load(catalog_file)
    catalog_pair=traverse_dict(catalog)

    keys_id = {"id", "title", "name", "model_id", "proj_id", "decode_id", "cells_id", "cell_id", "square_id", "raw_pixel_col", "analysis"}
    cartload_required_files = [] 
    cartload_optional_files = [] 
    basemap_files = []

    for key, value in catalog_pair:
        if key.startswith("assets.basemap"):
            if key.startswith("assets.basemap.sge") and not key.endswith("default"):
                cartload_required_files.append(value)
            elif not key.startswith("assets.basemap.sge") and not key.endswith("default"):
                basemap_files.append(value)
        else:
            subkey = key.split(".")[-1] if "." in key else key
            if subkey not in keys_id:
                if key.endswith("umap.tsv") or key.endswith("umap.png") or key.endswith("umap.ind_png") or key.endswith("umap.pmtiles") or subkey == "alias":
                    cartload_optional_files.append(value)
                    # print(repr(key)+"\t"+"cartload_optional_files")
                else:
                    cartload_required_files.append(value)
                    # print(repr(key)+"\t"+"cartload_required_files")


    # Deduplicate cartload_files
    cartload_required_files = list(set(cartload_required_files))
    cartload_optional_files = list(set(cartload_optional_files))

    return cartload_required_files, cartload_optional_files, basemap_files

def upload_aws(_args):
    parser = argparse.ArgumentParser(
        description=(
            "Upload files to S3 as specified in catalog.yaml. "
            "Supports single dataset upload (default) or collection upload via --in-list. "
            "For single uploads, use --in-dir/--s3-dir (and optional --catalog-yaml). "
            "For collection uploads, provide --in-list and point --in-dir/--s3-dir to parent directories."
        )
    )
    run_params = parser.add_argument_group("Run Options", "")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate Makefile and print commands; do not execute')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and rebuild from scratch')
    run_params.add_argument('--makefn', type=str, default="upload_aws.mk", help='Makefile name to write (default: upload_aws.mk)')
    run_params.add_argument('--n-jobs', type=int, default=2, help='Number of parallel jobs (default: 2)')

    key_params = parser.add_argument_group("Key Parameters")
    key_params.add_argument('--in-dir', required=True, help='Input dir (single) or parent dir with per-sample subdirs (collection)')
    key_params.add_argument('--s3-dir', required=True, help='S3 destination (single) or parent prefix for per-sample outputs (collection)')
    key_params.add_argument('--catalog-yaml', default=None, help='Path to catalog.yaml (single mode only; default: <in_dir>/catalog.yaml)')
    key_params.add_argument('--upload-basics-only', action="store_true", default=False, help='Upload only the cartload mandatory files from run_cartload*')
    key_params.add_argument('--upload-optional-only', action="store_true", default=False, help='Upload only the cartload optional files, i.e., umap, alias ')
    key_params.add_argument('--upload-basemap-only', action="store_true", default=False, help='Upload only additional basemap pmtiles, i.e., histology ')
    key_params.add_argument('--aws', type=str, default="aws", help='aws CLI executable (default: aws)')
    key_params.add_argument('--profile', type=str, default=None, help='AWS CLI profile to use (optional)')

    multi_params = parser.add_argument_group(
        "Optional Collection Upload Parameters",
        "Providing --in-list enables collection mode (multiple datasets). "
        "This is designed as a shortcut to upload files generated by run_cartload2_multi.",
    )
    multi_params.add_argument('--in-list', type=str, default=None, help='TSV of sample IDs (no header); enables collection mode')
    args = parser.parse_args(_args)

    def upload_single_to_aws(in_dir, s3_dir, catalog_f, args, mm):
        s3_catalog_f = f"{s3_dir}/catalog.yaml"
        aws_base = f"{args.aws}" + (f" --profile {args.profile}" if args.profile else "")

        # get the cartload output and basemaps files from catalog
        basics_files, optional_files, basemap_files = collect_files_from_yaml(catalog_f)
            
        # step 1. Upload cartload files to AWS (all files, except the non-SGE in basemaps)
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
        # Option 2: use traverse_dict to locate the files. (function: collect_files_from_yaml)
        # Note, if catalog.yaml file structure changes, this keys_id may need to be updated.
        
        upload_basics = True if not args.upload_basemap_only and not args.upload_optional_only else False
        upload_opt = True if not args.upload_basics_only and not args.upload_basemap_only else False
        upload_bm = True if not args.upload_basics_only and not args.upload_optional_only else False

        # print(";".join(basics_files)+"\n")
        # upload cartload basics
        if upload_basics:
            cmds=cmd_separator([], f"Uploading cartload files to AWS...")
            basics_prerequisites=[os.path.join(in_dir, filename) for filename in basics_files]
            basics_prerequisites.append(catalog_f)
            
            for filename in basics_files:
                file_path=os.path.join(in_dir, filename)
                s3_file_path = os.path.join(s3_dir, filename) 
                cmds.append(f"{aws_base} s3 cp {file_path} {s3_file_path}")
            
            cmds.append(f"{aws_base} s3 cp {catalog_f} {s3_catalog_f}")
            cartload_flag=os.path.join(in_dir, "cartload.aws.done")
            cmds.append(f"touch {cartload_flag}")
            mm.add_target(cartload_flag, basics_prerequisites, cmds)

        # print(";".join(optional_files)+"\n")
        # upload cartload optional
        if upload_opt:
            for filename in optional_files:
                cmds=cmd_separator([], f"Uploading optional cartload files to AWS...")
                file_path=os.path.join(in_dir, filename)
                s3_file_path = os.path.join(s3_dir, filename) 
                cmds.append(f'{aws_base} s3 cp "{file_path}" "{s3_file_path}" && {aws_base} s3 cp "{catalog_f}" "{s3_catalog_f}" && touch {file_path}.aws.done')
                mm.add_target(f"{file_path}.aws.done", [file_path], cmds)
        
        # print(";".join(basemap_files)+"\n")
        # Upload basemap files to AWS besides cartload
        if upload_bm:
            for filename in basemap_files:
                cmds=cmd_separator([], f"Uploading additional basemaps to AWS: {filename}...")
                file_path=os.path.join(in_dir, filename)
                s3_file_path = os.path.join(s3_dir, filename) 
                cmds.append(f'{aws_base} s3 cp "{file_path}" "{s3_file_path}" && {aws_base} s3 cp "{catalog_f}" "{s3_catalog_f}" && touch {file_path}.aws.done')
                mm.add_target(f"{file_path}.aws.done", [file_path], cmds)

    # check aws cli
    scheck_app(args.aws)

    # start mm
    mm = minimake()

    input_mode="single" if args.in_list is None else "multi"

    if input_mode == "single":
        in_dir = args.in_dir
        s3_dir = args.s3_dir
        assert os.path.exists(in_dir), f"Input directory not found: {in_dir} (--in-dir)"

        catalog_f = args.catalog_yaml if args.catalog_yaml is not None else os.path.join(args.in_dir, "catalog.yaml")
        assert os.path.exists(catalog_f), f"File not found: {catalog_f}" + ("(--catalog-yaml)" if args.catalog_yaml is not None else "(defined using --in-dir)")

        upload_single_to_aws(in_dir, s3_dir, catalog_f, args, mm)

    if input_mode == "multi":
        assert os.path.exists(args.in_list), f"File not found: {args.in_list} (--in-list)"
        # read the in_list
        df=pd.read_csv(args.in_list, sep="\t", dtype=str, header=None)
        
        if df.shape[1] < 1:
            raise ValueError(f"Input list {args.in_list} must have at least 1 columns: sample id")
        # extract the first column as input_ids
        input_ids=df[0].tolist()

        for input_id in input_ids:
            print(f"Processing input ID: {input_id}")
            
            in_dir = os.path.join(args.in_dir, input_id) 
            assert os.path.exists(in_dir), f"Input directory not found: {in_dir} (--in-dir with --in-list)"

            assert args.catalog_yaml is None, "When --in-list is provided, --catalog-yaml should not be specified. Each input directory must contain its own catalog.yaml file."
            catalog_f = os.path.join(in_dir, "catalog.yaml")
            assert os.path.exists(catalog_f), f"File not found: {catalog_f}" + ("(--catalog-yaml)" if args.catalog_yaml is not None else "(defined using --in-dir with --in-list)")
            
            # s3 dir
            with open(catalog_f, 'r') as stream:
                catalog_data = yaml.load(stream, Loader=yaml.FullLoader)
            if "id" in catalog_data:
                print(f"Extract s3 ID from catalog.yaml: {catalog_data['id']}; using this ID for s3 subdirectory name.")
                s3_id = catalog_data["id"]
            else:
                print(f"No 'id' found in catalog.yaml. Generating s3 ID from input ID: {input_id}")
                s3_id = input_id.replace("_", "-").lower()
            s3_dir = os.path.join(args.s3_dir, s3_id)
            
            # upload single data
            upload_single_to_aws(in_dir, s3_dir, catalog_f, args, mm)

    ## write makefile
    make_f = os.path.join(args.in_dir, args.makefn)
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

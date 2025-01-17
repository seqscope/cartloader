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
    # parser.add_argument("--histology", nargs="*", default=None, help="List of histology pmtiles to upload (default: None).")
    parser.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    parser.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    parser.add_argument('--makefn', type=str, default="upload_aws_by_catalog.mk", help='The file name of the Makefile to generate (default: {out-prefix}.mk)')
    parser.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')
    args = parser.parse_args(_args)

    # catalog files
    if args.catalog_yaml is None:
        catalog_f= os.path.join(args.in_dir, "catalog.yaml")
    else:
        catalog_f = args.catalog_yaml

    s3_catalog_f = f"{args.s3_dir}/catalog.yaml"

    # catalog info
    with open(catalog_f, "r") as catalog_file:
        catalog = yaml.safe_load(catalog_file)
    catalog_pair=traverse_dict(catalog)

    # start mm
    mm = minimake()

    # step 1. Upload cartload files to AWS （all files, except the nonsge in basemaps 
    cmds=cmd_separator([], f"Uploading cartload files to AWS...")
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
    cartload_files=[]
    basemap_files=[]

    for key, value in catalog_pair:
        # skip all basemaps except the sge
        if key.startswith("assets.basemap"):
            if key.startswith("assets.basemap.sge") and not key.endswith("default"):
                cartload_files.append(value)
            elif not key.startswith("assets.basemap.sge") and not key.endswith("default"):
                basemap_files.append(value)
        else:
            subkey=key.split(".")[-1] if "." in key else key
            if subkey not in keys_id:
                cartload_files.append(value)

    cartload_files=list(set(cartload_files))
    cartload_prerequisites=[os.path.join(args.in_dir, filename) for filename in cartload_files]
    cartload_prerequisites.append(catalog_f)
    for filename in cartload_files:
        file_path=os.path.join(args.in_dir, filename)
        s3_file_path = os.path.join(args.s3_dir, filename) 
        cmds.append(f"aws s3 cp {file_path} {s3_file_path}")
    cartload_flag=os.path.join(args.in_dir, "cartload.aws.done")
    # upload catalog.yaml to AWS
    cmds.append(f"aws s3 cp  {catalog_f} {s3_catalog_f}")
    cmds.append(f"touch {cartload_flag}")
    mm.add_target(cartload_flag, cartload_prerequisites, cmds)

    # step 2. Upload basemap files to AWS besides sge
    for filename in basemap_files:
        cmds=cmd_separator([], f"Uploading additional basemaps to AWS: {filename}...")
        file_path=os.path.join(args.in_dir, filename)
        s3_file_path = os.path.join(args.s3_dir, filename) 
        cmds.append(f"aws s3 cp {file_path} {s3_file_path}")
         # upload catalog.yaml to AWS
        cmds.append(f"aws s3 cp  {catalog_f} {s3_catalog_f}")
        cmds.append(f"touch {file_path}.aws.done")
        mm.add_target(f"{file_path}.aws.done", [file_path], cmds)

    # # # Step 2.2: By default, the histology should already registered in the catalog.yaml. 
    # # # Added this --histology option to faciliate additional histology files.
    # # if args.histology is not None:
    # #     for hist_f in args.histology:
    # #         hist_bn = os.path.basename(hist_f)
    # #         hist_copy = os.path.join(args.in_dir, hist_bn)
    # #         s3_hist_f=f"{args.s3_dir}/{hist_bn}"
    # #         hist_yaml_flag = f"{hist_copy}.yaml.done"
    # #         hist_aws_flag = f"{hist_copy}.aws.done"

    # #         if hist_f != hist_copy:
    # #             cmds=cmd_separator([], f"Copying histology files to input directory...")
    # #             cmds.append(f"cp {hist_f} {hist_copy}")
    # #             mm.add_target(f"{hist_copy}", [hist_f], cmds)

    # #         # # updating catalog.yaml
    # #         # cmds = cmd_separator([], f"Registering histology files to catalogy yaml...'")
    # #         # cmds.append(" ".join ([
    # #         #     "cartloader", "update_yaml_for_basemap",
    # #         #     "--catalog-yaml", catalog_f,
    # #         #     "--pmtiles", hist_f
    # #         # ]))
    # #         # cmds.append(f"touch hist_yaml_flag")
    # #         # mm.add_target(hist_yaml_flag, [hist_copy], cmds)
            
    # #         # uploading histology files to AWS
    # #         cmds.append(f"aws s3 cp {hist_copy} {s3_hist_f}")
    # #         cmds.append(f"touch {hist_flag}")
    # #         mm.add_target(hist_flag, [hist_copy], cmds)

    ## write makefile
    make_f=os.path.join(args.in_dir, args.makefn)
    mm.write_makefile(make_f)

    if args.dry_run:
        os.system(f"make -f {make_f} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs}")
    else:
        result = subprocess.run(f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}", shell=True)
        if result.returncode != 0:
            print(f"Error in executing the pipeline. Check the log files for more details.")
            sys.exit(1)

    # # Execute or print all commands
    # for cmd in commands:
    #     if isinstance(cmd, str):
    #         print(cmd)  # Always print echo statements
    #     else:
    #         run_command(cmd, args.dry_run)

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
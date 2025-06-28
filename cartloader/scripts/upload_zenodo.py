

import requests
import json, os, sys, argparse, inspect, glob
import yaml
from cartloader.scripts.upload_aws import collect_files_from_yaml

def get_existing_files(bucket_url, ACCESS_TOKEN):
    r = requests.get(
        bucket_url,
        params={'access_token': ACCESS_TOKEN},
    )
    print(f"Bucket Status: {r.status_code}")
    if r.status_code == 200:
        files = [file['filename'] for file in r.json()['contents']]
        return files
    elif r.status_code == 500:
        return []
    else:
        print(f"Error fetching existing files: {r.status_code} - {r.text}")
        exit(1)
        return []

def upload_file_to_bucket(filepath, bucket_url, ACCESS_TOKEN):
    filename = os.path.basename(filepath)
    with open(filepath, "rb") as fp:
        r = requests.put(
            f"{bucket_url}/{filename}",
            data=fp,
            params={'access_token': ACCESS_TOKEN},
        )
    return r

# use collect_files_from_yaml from upload aws
# # define a function to check if a variable is a list or a dict
# def collect_files_from_yaml(yaml_dat):
#     file_set = set()
#     def recurse_items(items):
#         if isinstance(items, dict):
#             for key, value in items.items():
#                 recurse_items(value)
#         elif isinstance(items, list):
#             for item in items:
#                 recurse_items(item)
#         elif isinstance(items, str):
#             file_set.add(items)
#     recurse_items(yaml_dat)
#     # sort the file list
#     file_list=sorted(list(file_set))
#     return file_list

def upload_zenodo(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Upload results files to Zenodo for cartoscope.")
    parser.add_argument('--zenodo-id', type=str, required=True, help='The Zenodo deposition ID.')
    parser.add_argument('--zenodo-token', type=str, required=True, help='The path of file for your access token for Zenodo.')
    parser.add_argument('--in-dir', type=str, required=True, help='The input directory where the files are located.')
    parser.add_argument('--upload-method', type=str, default="user_list", choices=["all", "catalog", "user_list"], help='Specify which files to upload to Zenodo. Choose "all" to upload all files in the input directory, "catalog" to upload files listed in the catalog YAML file, or "user_list" to upload files specified by the user.')
    parser.add_argument('--in-list', type=str, nargs='+', default=[], help='If --upload-method is "user_list", provide the name of files to be uploaded to Zenodo.')
    parser.add_argument('--catalog-yaml', type=str, default=None, help='If --upload-method is "catalog_yaml", provide the catalog YAML file that lists the files to be uploaded to Zenodo.')
    parser.add_argument('--overwrite', action='store_true', default=False, help='If set, overwrite existing files in the Zenodo bucket. If not set, skip existing files.')
    parser.add_argument('--dry-run', action='store_true', default=False, help='If set, overwrite existing files in the Zenodo bucket. If not set, skip existing files.')
    args = parser.parse_args(_args)

    # 1. Obtain info from zenodo deposition
    print(f"="*10)
    print(f" 1. Accessing Zenodo ")
    print(f"="*10)
    print(f" * Zenodo Token File: {args.zenodo_token}")
    ACCESS_TOKEN = open(args.zenodo_token).read().strip()

    print(f" * Accessing the deposition ID: {args.zenodo_id}")
    r = requests.get(f'https://zenodo.org/api/deposit/depositions/{args.zenodo_id}',
                    params={'access_token': ACCESS_TOKEN})
    
    if r.status_code != 200:
        print(f"Error fetching deposition: {r.status_code} - {r.text}")
        sys.exit(1)

    bucket_url = r.json()["links"]["bucket"]

    # check if any file exists in the bucket
    print(f" * Checking files in the bucket_URL {bucket_url}")
    existing_files = get_existing_files(bucket_url, ACCESS_TOKEN)
    n_exfiles=len(existing_files)
    print(f"    - Found {n_exfiles} existing files")

    # define input by the upload_method
    print(f"="*10)
    print(f" 2. Accessing input ")
    print(f"="*10)
    if args.upload_method == "all":
        in_files_raw = glob.glob(f"{args.in_dir}/*")
    elif args.upload_method == "catalog":
        if args.catalog_yaml is None:
            catalog_f= os.path.join(args.in_dir, "catalog.yaml")
        else:
            catalog_f = args.catalog_yaml
        assert os.path.exists(catalog_f), "Provide an invalid catalog yaml file when using --upload-method catalog"
        # with open(args.catalog_yaml, 'r') as file:
        #     yaml_dat = yaml.safe_load(file)
        # fn_list = collect_files_from_yaml(yaml_dat["assets"])
        cartload_files, basemap_files = collect_files_from_yaml(catalog_f)
        fn_list=list(cartload_files)+basemap_files
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in fn_list]
    elif args.upload_method == "user_list":
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in args.in_list]
    
    if len(in_files_raw) == 0:
        print("Error: No file to be uploaded.")
        sys.exit(1)

    if len(existing_files)>0:
        in_files_existing = [f for f in in_files_raw if os.path.basename(f) in existing_files]
        if len(in_files_existing) > 0:
            print(f" * Found {len(in_files_existing)} files that already exist in the Zenodo bucket:")
            for in_file in in_files_existing:
                print(f"    - {in_file}")
            if args.overwrite:
                print(" * Will overwrite existing files in the Zenodo bucket.")
                in_files = in_files_raw
            else:
                print(" * Will skip existing files in the Zenodo bucket.")
                in_files = [f for f in in_files_raw if f not in in_files_existing]
        else:
            print(f" * No input file exists in the Zenodo bucket.")
            in_files = in_files_raw
    else:
        in_files = in_files_raw
    
    n_infiles=len(in_files)
    print(f" * Found {n_infiles} to be uploaded...")
    # upload the files
    if args.dry_run:
        print(f"Showing the files to be uploaded to Zenodo")
        for in_file in in_files:
            print(f" * {in_file}")
    else:
        print(f"Uploading: Uploading files to Zenodo bucket {bucket_url}...")
        for in_file in in_files:
            try:
                res = upload_file_to_bucket(in_file, bucket_url, ACCESS_TOKEN)
                if res.status_code not in (200, 201):
                    print(f"    - Failed to upload {in_file}: {res.status_code} - {res.text}")
                else:
                    print(f"    - Successfully uploaded {in_file}")
            except Exception as e:
                print(f"    - Error uploading {in_file}: {e}")

        # check the status
        print(f"Check the status of the deposition {args.zenodo_id}")
        print(r.json())

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
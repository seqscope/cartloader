

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
        sys.exit(1)
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

def upload_zenodo(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Upload results files to Zenodo for cartoscope.")
    parser.add_argument('--zenodo-id', type=str, required=True, help='The Zenodo deposition ID.')
    parser.add_argument('--zenodo-token', type=str, required=True, help='The path of file for your access token for Zenodo.')
    parser.add_argument('--in-dir', type=str, required=True, help='The input directory where the files are located.')
    parser.add_argument('--upload-method', type=str, default="user_list", choices=["all", "catalog", "user_list"], help='Specify which files to upload to Zenodo. Choose "all" to upload all files in the input directory, "catalog" to upload files listed in the catalog YAML file, or "user_list" to upload files specified by the user.')
    parser.add_argument('--in-list', type=str, nargs='+', default=[], help='If --upload-method is "user_list", provide the name of files to be uploaded to Zenodo.')
    parser.add_argument('--catalog-yaml', type=str, default=None, help='If --upload-method is "catalog", provide the catalog YAML file that lists the files to be uploaded to Zenodo.')
    parser.add_argument('--overwrite', action='store_true', default=False, help='If set, overwrite existing files in the Zenodo bucket. If not set, skip existing files.')
    parser.add_argument('--dry-run', action='store_true', default=False, help='If set, simulate the upload without actually uploading files.')
    args = parser.parse_args(_args)

    # ==========
    # 1. Obtain info from zenodo deposition
    # ==========
    print(f"="*10)
    print(f" 1. Accessing Zenodo ")
    print(f"="*10)

    # get token file
    print(f" * Zenodo Token File: {args.zenodo_token}")
    ACCESS_TOKEN = open(args.zenodo_token).read().strip()
    if not ACCESS_TOKEN:
        print("ERROR: Zenodo access token file is empty.")
        sys.exit(1)

    # check deposition
    print(f" * Accessing the deposition ID: {args.zenodo_id}")
    r = requests.get(f'https://zenodo.org/api/deposit/depositions/{args.zenodo_id}',
                    params={'access_token': ACCESS_TOKEN})
    
    if r.status_code != 200:
        print(f"Error fetching deposition: {r.status_code} - {r.text}")
        sys.exit(1)
    
    # extract bucket url
    bucket_url = r.json()["links"]["bucket"]


    # check if any file exists in the bucket
    print(f" * Checking files in the bucket_URL {bucket_url}")
    existing_files = get_existing_files(bucket_url, ACCESS_TOKEN)
    n_exfiles=len(existing_files)
    print(f"    - Found {n_exfiles} existing files")

    # ==========
    # 2. Define input by the upload_method
    # ==========
    print("="*10)
    print(" 2. Accessing Input Files")
    print("="*10)
    if args.upload_method == "all":
        in_files_raw = glob.glob(os.path.join(args.in_dir, "*"))
    elif args.upload_method == "catalog":
        catalog_f = args.catalog_yaml or os.path.join(args.in_dir, "catalog.yaml")
        assert os.path.exists(catalog_f), f" * Missing catalog file: {catalog_f}"
        cartload_files, basemap_files = collect_files_from_yaml(catalog_f)
        fn_list = list(cartload_files) + basemap_files
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in fn_list]
    elif args.upload_method == "user_list":
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in args.in_list]
    else:
        raise ValueError(f"Unsupported upload method: {args.upload_method}")

    # raw list of input files
    if not in_files_raw:
        print(f" * No input files found using method: {args.upload_method}")
        sys.exit(1)
    print(f" * Initially, located {len(in_files_raw)} input file(s) for upload.")

    in_files_invalid = [f for f in in_files_raw if not os.path.exists(f)]
    if in_files_invalid:
        print(" * Error: The following input file(s) do not exist:")
        for f in in_files_invalid:
            print(f"    - {f}")
        sys.exit(1)

    # ==========
    # 3. Compare input with existing Zenodo files
    # ==========
    print("\n" + "="*10)
    print(" 3. Zenodo File Comparison")
    print("="*10)

    input_fnames = {os.path.basename(f): f for f in in_files_raw}
    existing_fnames = set(existing_files)

    input_overlap = [input_fnames[f] for f in input_fnames if f in existing_fnames]
    input_new = [input_fnames[f] for f in input_fnames if f not in existing_fnames]
    existing_only = [f for f in existing_files if f not in input_fnames]

    in_files = in_files_raw if args.overwrite else input_new

    if input_new:
        print(f"\n * New input file(s) (not in deposition): {len(input_new)}")
        for f in input_new:
            print(f"    - {f}")
    else:
        print("\n * No new input files.")

    if input_overlap:
        if args.overwrite:
            print(f"\n * Overlapping input file(s) (already in deposition): {len(input_overlap)} (will be overwritten)")
        else:
            print(f"\n * Overlapping input file(s) (already in deposition): {len(input_overlap)} (will be skipped)")
        print(f"    - {f}")
    else:
        print("\n * No overlapping input files.")

    if existing_only:
        print(f"\n * Existing file(s) in Zenodo not listed in input: {len(existing_only)}")
        for f in existing_only:
            print(f"    - {f}")
    else:
        print("\n * No additional existing files in the deposition.")

    files_after_upload = existing_fnames.union({os.path.basename(f) for f in in_files})
    n_expected = len(files_after_upload)
    n_infiles=len(in_files)

    print(f"\n * {n_infiles} file(s) will be uploaded.")
    print(f" * Expected total files in Zenodo bucket after upload: {n_expected}")

    if n_infiles == 0:
        print("\nWARNING: No files to be uploaded. Aborting upload.")
        return

    if n_expected > 100:
        print("\nWARNING: Total number of expected files exceeds 100 (Doesn't support by Zenodo). Aborting upload.")
        sys.exit(1)

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
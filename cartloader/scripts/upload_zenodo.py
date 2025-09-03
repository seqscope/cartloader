
import requests, json, os, sys, argparse, inspect, glob, yaml
from pathlib import Path
from datetime import date
from cartloader.scripts.upload_aws import collect_files_from_yaml

# === FILES ==
##  Zenodo’s bucket URL is actually a direct S3 endpoint, not a RESTful Zenodo API endpoint -- it cannot get the files that actually exist
# def list_existing_files(bucket_url, access_token):
#     response = requests.get(bucket_url, params={'access_token': access_token})
#     print(f"Bucket Status: {response.status_code}")
#     if response.status_code == 200:
#         return [file['filename'] for file in response.json()['contents']]
#     elif response.status_code == 500:
#         return []
#     else:
#         print(f"Error fetching existing files: {response.status_code} - {response.text}")
#         sys.exit(1)

def list_existing_files(deposition_id, access_token):
    url = f"https://zenodo.org/api/deposit/depositions/{deposition_id}/files"
    response = requests.get(url, params={'access_token': access_token})
    print(f"File list status: {response.status_code}")
    if response.status_code == 200:
        return [file['filename'] for file in response.json()]
    else:
        print(f"Error fetching files: {response.status_code} - {response.text}")
        return []

def upload_file(filepath, bucket_url, access_token):
    filename = os.path.basename(filepath)
    with open(filepath, "rb") as fp:
        response = requests.put(
            f"{bucket_url}/{filename}",
            data=fp,
            params={'access_token': access_token},
        )
    return response

# === META DATA ==
def build_metadata(args, existing_metadata=None):
    # existing_metadata should be provided if the deposition exists and existing metadata should be used. 
    if existing_metadata is None:
        if not args.title:
            raise ValueError("--title is required when creating a new deposition")
        if not args.creators:
            raise ValueError("--creators is required and must include at least one person")

        creators = [{"name": name} for name in args.creators]

        metadata = {
            "title": args.title,
            "upload_type": args.upload_type,
            "description": args.description or "No description provided.",
            "creators": creators,
            "publication_date": date.today().isoformat()
        }
    else:
        # Start with a copy of the existing metadata
        metadata = existing_metadata.copy()

        # Update with any provided args
        if args.title:
            metadata["title"] = args.title
        if args.upload_type:
            metadata["upload_type"] = args.upload_type
        if args.description:
            metadata["description"] = args.description
        if args.creators:
            metadata["creators"] = [{"name": name} for name in args.creators]
        if "publication_date" not in metadata:
            metadata["publication_date"] = date.today().isoformat()

    return metadata

def update_deposition_metadata(deposition_id, access_token, metadata):
    url = f"https://zenodo.org/api/deposit/depositions/{deposition_id}"
    response = requests.put(url, params={'access_token': access_token}, json={"metadata": metadata}, headers={"Content-Type": "application/json"})
    if response.status_code == 200:
        print(f" * Successfully updated metadata for deposition {deposition_id}")
    else:
        raise RuntimeError(f"Failed to update deposition metadata. Please manually add that metadata: {response.status_code} - {response.text}")

# === DEPOSITION ==
def create_deposition(access_token, metadata=None):
    headers = {"Content-Type": "application/json"}
    url = "https://zenodo.org/api/deposit/depositions"
    payload = {"metadata": metadata} if metadata else {}

    response = requests.post(url, params={"access_token": access_token}, json=payload, headers=headers)

    if response.status_code == 201:
        deposition = response.json()
        deposition_id = str(deposition["id"])
        bucket_url = deposition["links"]["bucket"]
        print(f" * Created new deposition: ID = {deposition_id}")
        return deposition_id, bucket_url
    else:
        raise RuntimeError(f"Error creating deposition: {response.status_code} - {response.text}")

def publish_deposition(deposition_id, access_token):
    url = f"https://zenodo.org/api/deposit/depositions/{deposition_id}/actions/publish"
    response = requests.post(url, params={'access_token': access_token})
    if response.status_code == 202:
        print("Deposition published successfully.")
    else:
        print(f"Failed to publish. Please publish it manually: {response.status_code} - {response.text}")

# == NEW VERSION ==
def create_new_version(old_id, token):
    url = f"https://zenodo.org/api/deposit/depositions/{old_id}/actions/newversion"
    r = requests.post(url, params={"access_token": token})
    if r.status_code != 201:
        raise RuntimeError(f"Failed to create new version: {r.status_code} - {r.text}")
    print("New version draft created.")
    return r.json()["links"]["latest_draft"]


def upload_zenodo(_args):
    # Categorized argument groups
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Upload files to Zenodo.")

    # I/O arguments
    group_io = parser.add_argument_group("Input/Output")
    group_io.add_argument('--in-dir', type=str, required=True, help='(Required) Path to the input directory containing files to upload.')
    group_io.add_argument('--upload-method', type=str, default="all", choices=["all", "catalog", "user_list"], help='Method to determine which files to upload: "all" uploads every file in --in-dir; "catalog" uses filenames listed in a catalog YAML file; "user_list" uploads only files specified via --in-list.')
    group_io.add_argument('--in-list', type=str, nargs='+', default=[], help='(Required if --upload-method is "user_list") One or more filenames to upload, e.g., "--in-list fileA.tif fileB.tif".')
    group_io.add_argument('--catalog-yaml', type=str, default=None, help='(Required if --upload-method is "catalog") Path to the catalog YAML file listing files to upload (defaults to <in_dir>/catalog.yaml)')

    # Zenodo credentials and config
    group_zenodo = parser.add_argument_group("Zenodo Configuration")
    group_zenodo.add_argument('--zenodo-token', type=str, required=True, help='Path to a file containing your Zenodo access token.')
    group_zenodo.add_argument('--zenodo-deposition-id', type=str, default=None,  help='(Optional) ID of an existing Zenodo deposition to upload files into.'
                                                                                        'If the deposition is published, a new version will be automatically created. '
                                                                                        'If not provided, a new deposition will be created.')
    # group_zenodo.add_argument('--create-new-deposition', action='store_true', help='Create a new Zenodo deposition instead of using an existing one. Must specify exactly one of --zenodo-deposition-id or --create-new-deposition.')
    # group_zenodo.add_argument('--create-new-version', action='store_true', help='If set, create a new version of an existing published Zenodo deposition. Requires --zenodo-deposition-id to specify the DOI or ID of the published deposition.')

    # Metadata fields
    group_meta = parser.add_argument_group("Deposition Metadata")
    group_meta.add_argument('--title', type=str, default=None, help='Title of the deposition. Required if creating a new deposition or if the existing deposition does not have a title.')
    group_meta.add_argument('--upload-type', type=str, default='dataset', choices=['dataset', 'software', 'publication', 'poster', 'presentation', 'image', 'video', 'lesson', 'other'], help='Type of upload for the Zenodo deposition (default: dataset). Required if creating a new deposition.')
    group_meta.add_argument('--creators', type=str, nargs='+', default=[], help='List of creators in "Lastname, Firstname" format. Each name should be quoted. Required if creating a new deposition or if the existing deposition does not have creator information.')
    group_meta.add_argument('--description', type=str, default=None, help='(Optional) Description of the deposition.')

    # Behavior flags
    group_flags = parser.add_argument_group("Options")
    group_flags.add_argument('--publish', action='store_true', default=False, help='If set, publish the deposition automatically after upload. '
                                                                                    'Recommended to leave this DISABLED and publish manually after verifying the deposition via the Zenodo web interface.'
                                                                                    )
    group_flags.add_argument('--overwrite', action='store_true', default=False, help='Overwrite existing files with the same name.')
    group_flags.add_argument('--dry-run', action='store_true', default=False, help='Simulate the upload process without making changes.')

    args = parser.parse_args(_args)

    # ==========
    # 1. Obtain info from zenodo deposition
    # ==========
    print(f"-"*10)
    print(f" 1. Accessing Zenodo ")
    print(f"-"*10)

    # get token file
    print(f" * Zenodo Token File: {args.zenodo_token}")
    
    with open(args.zenodo_token) as f:
        ACCESS_TOKEN = f.read().strip()

    if not ACCESS_TOKEN:
        print("ERROR: Zenodo access token file is empty.")
        sys.exit(1)

    if args.zenodo_deposition_id:
        # print(f" * Accessing the deposition ID: {args.zenodo_deposition_id}")
        # if args.create_new_version:
        #     print(f"    - Create a new version: {args.zenodo_deposition_id}")
        #     draft_url = create_new_version(args.zenodo_deposition_id, ACCESS_TOKEN)
        #     response = requests.get(draft_url, params={"access_token": ACCESS_TOKEN})
        #     if response.status_code != 200:
        #         raise RuntimeError(f"Failed to retrieve new draft: {response.status_code} - {response.text}")
        #     new_deposition_id = response.json()["id"]
        #     print(f"    - New version is created. New deposition ID: {new_deposition_id}")
        #     args.zenodo_deposition_id = new_deposition_id
        # else:
        #     response = requests.get(f'https://zenodo.org/api/deposit/depositions/{args.zenodo_deposition_id}',
        #                     params={'access_token': ACCESS_TOKEN})
        #     if response.status_code != 200:
        #         print(f"Error fetching deposition: {response.status_code} - {response.text}")
        #         sys.exit(1)
        
        # Update: automatically detect whether it should create a new version -- if the status is submitted, use a new version
        print(f" * Using an existing Zenodo deposition ({args.zenodo_deposition_id})...")

        response = requests.get(
            f'https://zenodo.org/api/deposit/depositions/{args.zenodo_deposition_id}',
            params={'access_token': ACCESS_TOKEN}
        )

        if response.status_code != 200:
            print(f"Error fetching deposition: {response.status_code} - {response.text}")
            sys.exit(1)

        deposition_metadata = response.json()
        is_published = deposition_metadata.get("submitted", False)  # Or check "state" == "done"

        if is_published:
            print(f"    - Detected published deposition. Creating a new version of {args.zenodo_deposition_id}")
            draft_url = create_new_version(args.zenodo_deposition_id, ACCESS_TOKEN)
            response = requests.get(draft_url, params={"access_token": ACCESS_TOKEN})
            if response.status_code != 200:
                raise RuntimeError(f"Failed to retrieve new draft: {response.status_code} - {response.text}")
            args.zenodo_deposition_id = response.json()["id"]
            print(f"    - New version created. New deposition ID: {args.zenodo_deposition_id }")
        else:
            print(f"    - Deposition is still a draft. Using existing deposition ID: {args.zenodo_deposition_id}")

        # update metadata to deposition
        print(f" * Updating metadata")
        raw_metadata=response.json().get("metadata", {})
        metadata = build_metadata(args, existing_metadata=raw_metadata)
        update_deposition_metadata(args.zenodo_deposition_id, ACCESS_TOKEN, metadata)

        # check if any file exists in the bucket
        bucket_url = response.json()["links"]["bucket"]
        print(f" *  Checking files in the bucket_URL {bucket_url}")
        existing_files = list_existing_files(args.zenodo_deposition_id, ACCESS_TOKEN)
        print(f"      Found {len(existing_files)} existing files")
    else:
        print(" * Creating a new Zenodo deposition...")
        print(f"    - Building metadata")
        metadata = build_metadata(args)
        args.zenodo_deposition_id, bucket_url = create_deposition(ACCESS_TOKEN, metadata)
        print(f"    - New deposition ID: {args.zenodo_deposition_id}")

        existing_files=[]
    
    # ==========
    # 2. Define input by the upload_method
    # ==========
    print("-"*10)
    print(" 2. Accessing Input Files")
    print("-"*10)
    
    print(f" * Upload method: {args.upload_method}")
    if args.upload_method == "all":
        in_files_raw = glob.glob(os.path.join(args.in_dir, "*"))
    elif args.upload_method == "catalog":
        catalog_f = args.catalog_yaml or os.path.join(args.in_dir, "catalog.yaml")
        assert os.path.exists(catalog_f), f"File not found: {catalog_f} (--catalog-yaml)"
        cartload_files, basemap_files = collect_files_from_yaml(catalog_f)
        fn_list = list(cartload_files) + basemap_files
        cartload_files_raw = [os.path.join(args.in_dir, fn) for fn in list(cartload_files)]
        basemap_files_raw = [os.path.join(args.in_dir, fn) for fn in list(basemap_files)]
        in_files_raw = list(set(cartload_files_raw + basemap_files_raw))
    elif args.upload_method == "user_list":
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in args.in_list]
    else:
        raise ValueError(f"Unsupported upload method: {args.upload_method}")

    print(f" *  Checking invalid input files")
    print(f"    -  Initially, located {len(in_files_raw)} input file(s) for upload.")
    in_files_invalid = [f for f in in_files_raw if not os.path.exists(f)]
    if in_files_invalid:
        print("    - Error: The following input file(s) do not exist:")
        for f in in_files_invalid:
            print(f"    - {f}")
        sys.exit(1)
    
    print(f" * Checking Zenodo capacity")
    existing_fnames = set(existing_files)
    files_after_upload = existing_fnames.union({os.path.basename(f) for f in in_files_raw})
    n_expected = len(files_after_upload)
    print(f"    - Expected total files in Zenodo bucket after upload: {n_expected}")
    if n_expected > 100:
        print("\n    - WARNING: Total number of expected files exceeds 100 (Doesn't support by Zenodo). Aborting upload.")
        sys.exit(1)

    print("-"*10)
    print(" 3. Uploading files")
    print("-"*10)
    
    def uploading(in_files, dry_run, touch_flag=False, flag_suffix="zenodo.done"):
        failed_list=[]
        if dry_run:
            print(f"\nShowing the files to be uploaded to Zenodo")
            for in_file in in_files:
                print(f" * {in_file}")
        else:
            print(f"\nUploading files to Zenodo bucket {bucket_url}...")
            for in_file in in_files:
                try:
                    res = upload_file(in_file, bucket_url, ACCESS_TOKEN)
                    if res.status_code not in (200, 201):
                        print(f"    - Failed to upload {in_file}: {res.status_code} - {res.text}")
                        failed_list.append(in_file)
                    else:
                        print(f"    - Successfully uploaded {in_file}")
                        if touch_flag:
                            flag_f = f"{in_file}.{flag_suffix}"
                            Path(flag_f).touch(exist_ok=True)
                except Exception as e:
                    print(f"    - Error uploading {in_file}: {e}")
        return failed_list

    def process_uploading_by_list(in_files_raw, existing_files, force_upload_files=[], touch_flag=False, flag_suffix="zenodo.done", overwrite=False, dry_run=False):
        if not in_files_raw:
            raise ValueError(f" * No input files found")

        failed_list = []
        print(f" * Checking overlapping with existing files in the Zenodo deposition (overwrite mode: {'On' if args.overwrite else 'Off'})")
        input_fnames = {os.path.basename(f): f for f in in_files_raw}
        existing_only = [f for f in existing_files if f not in input_fnames]
        input_overlap = [input_fnames[f] for f in input_fnames if f in existing_fnames]
        input_only = [input_fnames[f] for f in input_fnames if f not in existing_fnames]
        in_files = in_files_raw if overwrite else input_only
        
        if len(force_upload_files) > 0 :
            print(f"    - Files to upload regardless of existence. {";".join(force_upload_files)}")
            in_files.extend(force_upload_files)

        if len(in_files) == 0:
            print("WARNING: No files to be uploaded. Aborting upload.")
        else:
            print(f" * Uploading (dry-run mode)") if dry_run else print(f"2) Uploading") 
            if input_only:
                print(f"\n    - New input file(s) (not in deposition): N={len(input_only)}")
                failed_list1=uploading(input_only, dry_run, touch_flag=False, flag_suffix="zenodo.done")
                failed_list.extend(failed_list1)
            else:
                print("\n    - New input file(s) (not in deposition): N=0")

            if input_overlap:
                if args.overwrite:
                    print(f"\n    - Overlapping input file(s) (already in deposition): N={len(input_overlap)} (will be overwritten)")
                    failed_list2=uploading(input_overlap, touch_flag=False, flag_suffix="zenodo.done")
                    failed_list.extend(failed_list2)
                else:
                    print(f"\n    - Overlapping input file(s) (already in deposition): N={len(input_overlap)} (will be skipped)")
                    for f in input_overlap:
                        print(f"        * {f}")
            else:
                print("\n * Overlapping file(s): N=0")

            print(f" * Showing existing files in the deposition") 
            if existing_only:
                print(f"\n    - Existing file(s) in Zenodo not listed in input: N={len(existing_only)}")
                for f in existing_only:
                    print(f"        * {f}")
            else:
                print(f"\n    - Existing file(s) in Zenodo not listed in input: N=0")
        return failed_list

    failed_list=[]
    if args.upload_method == "all" or args.upload_method == "user_list":
        failed_sublist=process_uploading_by_list(in_files_raw, existing_files, force_upload_files=[], touch_flag=False, flag_suffix="zenodo.done", overwrite=args.overwrite, dry_run=args.dry_run)
        failed_list.extend(failed_sublist)
    elif args.upload_method == "catalog":
        print(f"\n1) Upload: tiled map data for SGE (with or without FICTURE)\n")
        failed_sublist=process_uploading_by_list(cartload_files_raw, existing_files, force_upload_files=[], touch_flag=False, flag_suffix="zenodo.done", overwrite=args.overwrite, dry_run=args.dry_run)
        failed_list.extend(failed_sublist)
        cartload_flag=os.path.join(args.in_dir, "cartload.zenodo.done")
        Path(cartload_flag).touch(exist_ok=True)
        if len(basemap_files_raw) > 0:
            print(f"\n2) Upload: tiled map data for background images, such as histology images\n")
            failed_sublist=process_uploading_by_list(basemap_files_raw, existing_files, force_upload_files=[], touch_flag=True, flag_suffix="zenodo.done", overwrite=args.overwrite, dry_run=args.dry_run)
            failed_list.extend(failed_sublist)
    if args.publish and not args.dry_run:
        print(f"\n Publishing the deposition {args.zenodo_deposition_id} ...")
        failed_sublist=publish_deposition(args.zenodo_deposition_id, ACCESS_TOKEN)
        failed_list.extend(failed_sublist)

    if failed_list:
        print("-"*10)
        print(" 4. Files failed to be uploaded")
        print("-"*10)
        print(f"N={len(failed_list)}")
        for failed_f in failed_list:
            print(f" * {failed_f}")

        print("\nTry the following commands to re-upload those failed files:")
        zenodo_cmd = " ".join([
            "cartloader", "upload_zenodo",
            f"--in-dir {args.in_dir}",
            f"--upload-method list",
            f"--zenodo-token {args.zenodo_token}",
            f"--zenodo-deposition-id {args.zenodo_deposition_id}",
            f"--in-list {" ".join(failed_list)}"
        ])
        print(zenodo_cmd)
    

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

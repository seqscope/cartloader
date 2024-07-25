

import requests
import requests
import json, os, sys, argparse, inspect, glob
import yaml

def get_existing_files(bucket_url, ACCESS_TOKEN):
    """
    Gets a list of files that already exist in the specified bucket.
    :param bucket_url: The URL of the bucket to check.
    :param ACCESS_TOKEN: The access token used for authentication.
    :return: A list of filenames that already exist in the bucket.
    """
    r = requests.get(
        bucket_url,
        params={'access_token': ACCESS_TOKEN},
    )
    if r.status_code == 200:
        files = [file['filename'] for file in r.json()['contents']]
        return files
    else:
        print(f"Error fetching existing files: {r.status_code} - {r.text}")
        exit(1)
        return []

def upload_file_to_bucket(r, filepath, bucket_url, ACCESS_TOKEN):
    """
    Uploads a given file to the specified bucket URL using an access token.
    :param filepath: The path to the file to be uploaded.
    :param bucket_url: The URL of the bucket to which the file will be uploaded.
    :param ACCESS_TOKEN: The access token used for authentication.
    """
    filename = os.path.basename(filepath)
    with open(filepath, "rb") as fp:
        r = requests.put(
            f"{bucket_url}/{filename}",
            data=fp,
            params={'access_token': ACCESS_TOKEN},
        )
    return r

# define a function to check if a variable is a list or a dict
def collect_files_from_yaml(yaml_dat):
    """
    Collects and returns a list of unique file paths from the given YAML file.
    :param yaml_file: Path to the YAML file.
    :return: A list of unique file paths.
    """
    file_set = set()
    def recurse_items(items):
        if isinstance(items, dict):
            for key, value in items.items():
                recurse_items(value)
        elif isinstance(items, list):
            for item in items:
                recurse_items(item)
        elif isinstance(items, str):
            file_set.add(items)
    recurse_items(yaml_dat)
    # sort the file list
    file_list=sorted(list(file_set))
    return file_list

def upload_zenodo(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Upload results files to Zenodo for cartoscope.")
    parser.add_argument('--zenodo-id', type=str, required=True, help='The Zenodo deposition ID.')
    parser.add_argument('--zenodo-token', type=str, required=True, help='The path of file for your access token for Zenodo.')
    parser.add_argument('--in-dir', type=str, required=True, help='The input directory where the files are located.')
    parser.add_argument('--token-path', type=str, required=True, help='The path to your Zenodo access token file.')
    parser.add_argument('--upload-method', type=str, default="user_list", choices=["all", "yaml_list", "user_list"], help='Specify which files to upload to Zenodo. Choose "all_files" to upload all files in the input directory, "yaml_list" to upload files listed in a YAML file, or "user_list" to upload files specified by the user.')
    parser.add_argument('--in-list', type=str, nargs='+', default=[], help='If --upload-method is "user_list", provide the name of files to be uploaded to Zenodo.')
    parser.add_argument('--in-yaml', type=str, default=None, help='If --upload-method is "yaml_list", provide the YAML file that lists the files to be uploaded to Zenodo.')
    parser.add_argument('--overwrite', action='store_true', default=False, help='If set, overwrite existing files in the Zenodo bucket. If not set, skip existing files.')
    args = parser.parse_args(_args)

    # define input by the upload_method
    if args.upload_method == "all_files":
        in_files_raw = glob.glob(f"{args.input_dir}/*")
    elif args.upload_method == "yaml_list":
        assert args.in_yaml is not None, "Please provide a YAML file that lists the files to be uploaded to Zenodo."
        with open(args.in_yaml, 'r') as file:
            yaml_dat = yaml.safe_load(file)
        fn_list = collect_files_from_yaml(yaml_dat["assets"])
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in fn_list]
    elif args.upload_method == "user_list":
        in_files_raw = [os.path.join(args.in_dir, fn) for fn in args.in_list]
    
    if len(in_files_raw) == 0:
        print("Error: No file to be uploaded.")
        sys.exit(1)
    
    print(f"Input: Find {len(in_files_raw)} files to upload to Zenodo:")
    for in_file in in_files_raw:
        print(f"    - {in_file}")

    # read the access token from a txt file /home/weiqiuc/zenodo_token.txt
    ACCESS_TOKEN = open(args.zenodo_token).read().strip()

    # get the deposition
    print(f"Get the deposition {args.id}")
    r = requests.get(f'https://zenodo.org/api/deposit/depositions/{args.id}',
                    params={'access_token': ACCESS_TOKEN})
    
    if r.status_code != 200:
        print(f"Error fetching deposition: {r.status_code} - {r.text}")
        sys.exit(1)

    bucket_url = r.json()["links"]["bucket"]

    # locate the existing files in the bucket
    existing_files = get_existing_files(bucket_url, ACCESS_TOKEN)
    in_files_existing = [f for f in in_files_raw if os.path.basename(f) in existing_files]
    print(f"Check Overlapping: ")
    if len(in_files_existing) > 0:
        print(f" - Found {len(in_files_existing)} files that already exist in the Zenodo bucket:")
        for in_file in in_files_existing:
            print(f"    - {in_file}")
        if args.overwrite:
            print(" - Will overwrite existing files in the Zenodo bucket.")
            in_files = in_files_raw
        else:
            print(" - Will skip existing files in the Zenodo bucket.")
            in_files = [f for f in in_files_raw if f not in in_files_existing]
    else:
        print(f" - None of the files exists in the Zenodo bucket.")
        in_files = in_files_raw

    # upload the files
    print(f"Uploading: Uploading files to Zenodo bucket {bucket_url}...")
    for in_file in in_files:
        try:
            r = upload_file_to_bucket(in_file, bucket_url, ACCESS_TOKEN)
            if r.status_code == 200:
                print(f"    - Successfully uploaded: {in_file}")
            else:
                print(f"    - Failed to upload {in_file}: {r.status_code} - {r.text}")
        except Exception as e:
            print(f"    - Error uploading {in_file}: {e}")

    # check the status
    print(f"Check the status of the deposition {args.id}")
    print(r.json())

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
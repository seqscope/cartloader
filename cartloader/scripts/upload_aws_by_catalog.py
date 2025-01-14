import os
import subprocess
import argparse

def run_command(command, dry_run):
    if dry_run:
        print("Dry-run command:", " ".join(command))
    else:
        subprocess.run(command, check=True)

def upload_catalog(in_dir, catalog, s3_dir, histology, dry_run):
    commands = []
    catalog_path = os.path.join(in_dir, catalog)
    s3_catalog_path = f"s3://{s3_dir}/{catalog}"

    commands.append(f"echo '============================================================================'")
    commands.append(f"echo 'Uploading FICTURE results to AWS...'")
    commands.append(f"echo '============================================================================'")
    with open(catalog_path, "r") as catalog_file:
        for line in catalog_file:
            if "." in line:
                # Extract the file name
                file_name = line.strip().split()[-1]
                file_path = os.path.join(in_dir, file_name)
                s3_file_path = f"s3://{s3_dir}/{file_name}"

                commands.append(f"echo 'Uploading {file_name} to S3...'")
                commands.append(["aws", "s3", "cp", file_path, s3_file_path])

    # Step 2: If histology is provided, add additional commands
    if histology is not None:
        for hist_path in histology:
            hist_bn = os.path.basename(hist_path)
            s3_hist_path = f"s3://{s3_dir}/{hist_bn}"

            # Add commands for uploading pmtiles and updating catalog.yaml
            commands.append(f"echo '============================================================================'")
            commands.append(f"echo 'Updating histology file to AWS...'")
            commands.append(f"echo '============================================================================'")
            commands.append(f"echo 'Uploading {hist_bn}.pmtiles to S3...'")
            commands.append(["aws", "s3", "cp", hist_path, s3_hist_path])
            commands.append([
                "cartloader", "update_yaml_for_basemap",
                "--yaml", catalog_path,
                "--pmtiles", hist_path
            ])
            commands.append(f"echo 'Re-uploading updated {catalog} to S3...'")

    # Step 3: Upload catalog.yaml to S3
    commands.append(f"echo '============================================================================'")
    commands.append(f"echo 'Uploading {catalog} to AWS...'")
    commands.append(f"echo '============================================================================'")
    commands.append(["aws", "s3", "cp", catalog_path, s3_catalog_path])

    # Execute or print all commands
    for cmd in commands:
        if isinstance(cmd, str):
            print(cmd)  # Always print echo statements
        else:
            run_command(cmd, dry_run)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Upload files to S3")
    parser.add_argument("--in_dir", help="Directory containing the files")
    parser.add_argument("--catalog", default="catalog.yaml", help="catalog yaml file")
    parser.add_argument("--s3_dir", help="S3 directory name. Typically, s3://<bucket>/<ID>")
    parser.add_argument("--histology", nargs="*", default=None, help="List of histology pmtiles to upload")
    parser.add_argument("--dry-run", action="store_true", help="Print commands instead of running them. Default: False")

    args = parser.parse_args()
    upload_catalog(args.in_dir, args.catalog, args.s3_dir, args.histology, args.dry_run)


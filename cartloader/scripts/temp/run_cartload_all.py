# Purpose:
# This should help to run sge_convert, run_ficture, run_cartload_join  and upload to aws
# This only provide the minimal argument.

import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd
from cartloader.utils.utils import cmd_separator, scheck_app, find_major_axis, add_param_to_cmd

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""This script will take the raw file from different platform and generate a script to apply sge_convert, run_fixture, run_cartload_join and upload_aws_by_catalog.""")
    # run params
    run_params = parser.add_argument_group("Run Options", "")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--jobfn', type=str, default="from_raw_to_aws.job", help='Job filename (default: from_raw_to_aws.job)')
    
    # input params
    inout_params = parser.add_argument_group("Input Options", "")
    inout_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st"], required=True, help='Platform of the raw input file to infer the format of the input file.')
    inout_params.add_argument('--in-parquet', type=str, default="tissue_positions.parquet", help='Path to the input raw parquet file for spatial coordinates. Required for 10x_visium_hd platform (default: tissue_positions.parquet).') # naming convention: tissue_positions.parquet
    inout_params.add_argument('--in-sge', type=str, default=None, help='Path to the input SGE directory. Required for 10x_visium_hd and seqscope platform. ')
    inout_params.add_argument('--in-csv', type=str, default=None, help='Path to the input raw CSV/TSV file. Required for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st platforms (default: None).')
    inout_params.add_argument('--out-dir', type=str, default=None, help='Path to the output directory.')
    
    # sge_convert params
    sge_params = parser.add_argument_group("sge_convert Options", "")
    sge_params.add_argument('--filter-by-density', action='store_true', default=False, help='Filter the transcript-indexed SGE file by density (default: False)')
    sge_params.add_argument('--units-per-um', type=float, default=1.00, help='Coordinate unit per um (conversion factor) (default: 1.00)') 
    sge_params.add_argument('--scale-json', type=str, default=None, help="For 10x_visium_hd datasets, users could use --scale-json to provide the path to the scale json file for calculating units-per-um (default: None). Typical naming convention: scalefactors_json.json")
    sge_params.add_argument('--csv-colnames-count', type=str, default=None, help='Column name for feature expression count. If not provided, a count of 1 will be added for a feature in a pixel (default: MIDCounts for bgi_stereoseq; MIDCount for nova_st; None for the rest platforms).')
    sge_params.add_argument('--exclude-feature-regex', type=str, default="^(BLANK|Blank-|NegCon|NegPrb)", help='A regex pattern of feature/gene names to be excluded (default: "^(BLANK|Blank-|NegCon|NegPrb)"). To avoid filtering features by name using regex, use --exclude-feature-regex "". ')
    sge_params.add_argument('--print-removed-transcripts', action='store_true', default=False, help='For debugging purposes, print the list of features removed based on the filtering criteria (default: False). Currently it only works for datasets from 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope or pixel_seq.')
    
    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools.")
    env_params.add_argument('--spatula', type=str, help='Path to spatula binary.')    

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_cartload_sge2aws:
    # Parse command-line arguments
    args = parse_arguments(_args)

    # Check the availability of the required tools
    scheck_app(args.spatula, "spatula")

    # Create the output directory if it does not exist
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # start mm
    mm = minimake()




if __name__ == "__main__":
    # Get the path to the cartloader repository
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd
from cartloader.utils.utils import cmd_separator, scheck_app, find_major_axis, add_param_to_cmd

# This should help to run sge_convert, run_ficture, run_cartload_join  and upload to aws
# This only provide the minimal argument.

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
    sge_params = parser.add_argument_group("Input Options", "")
    sge_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st"], required=True, help='Platform of the raw input file to infer the format of the input file.')
    sge_params.add_argument('--units-per-um', type=float, default=1.00, help='Coordinate unit per um (conversion factor) (default: 1.00)') 
    sge_params.add_argument('--scale-json', type=str, default=None, help="For 10x_visium_hd datasets, users could use --scale-json to provide the path to the scale json file for calculating units-per-um (default: None). Typical naming convention: scalefactors_json.json")
    sge_params.add_argument('--precision-um', type=int, default=2, help='Number of digits to store the transcript coordinates (only if --px_to_um is in use). Set it to 0 to round to integer (default: 2)')
    sge_params.add_argument('--filter-by-density', action='store_true', default=False, help='Filter the transcript-indexed SGE file by density (default: False)')
    sge_params.add_argument('--in-parquet', type=str, default="tissue_positions.parquet", help='Path to the input raw parquet file for spatial coordinates. Required for 10x_visium_hd platform (default: tissue_positions.parquet).') # naming convention: tissue_positions.parquet
    sge_params.add_argument('--in-sge', type=str, default=os.getcwd(), help='Path to the input SGE directory. Required for 10x_visium_hd and seqscope platform. Defaults to the current working directory.')
    sge_params.add_argument('--in-csv', type=str, default=None, help='Path to the input raw CSV/TSV file. Required for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st platforms (default: None).')

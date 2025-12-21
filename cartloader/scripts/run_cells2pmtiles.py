import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, inspect
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, execute_makefile, flexopen

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Split and convert transcripts TSV file into pmtiles")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-tsv', type=str, help='Input Format TSV/CSV (possibly gzipped) file containing the X/Y coordinates and gene expression counts per spot')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='The output prefix. New directory will be created if needed')  
    inout_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    conv_params = parser.add_argument_group("Parameters for pmtiles conversion")
    conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level')
    conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level')
    conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles')
    conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles')
    conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory to be used (default: out-dir/tmp; specify /tmp if needed)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_cells2pmtiles(_args):
    """
    Convert cross-platform molecules (TSV) file to pmtiles
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_cell2pmtiles" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    ## rename X/Y to lon/lat for tippecanoe compatibility
    if args.col_rename is None:
        args.col_rename = []
        args.col_rename.append("X:lon")
        args.col_rename.append("Y:lat")
        args.col_rename.append("cluster:topK")

    # start mm
    mm = minimake()

    # create output directory if needed
    out_dir = os.path.dirname(args.out_prefix)
    out_base = os.path.basename(args.out_prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(out_dir, "tmp")
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir, exist_ok=True)

    ## convert TSV into CSV files for pmtiles conversion
    col_rename_dict = {}
    for col_rename in args.col_rename:
        old_name, new_name = col_rename.split(":")
        col_rename_dict[old_name] = new_name

    logger.info("  * Converting input TSV into CSV for pmtiles conversion")
    with flexopen(args.in_tsv, 'rt') as rf, open(f"{args.out_prefix}.csv", 'wt') as wf:
        in_hdrs = rf.readline().strip().split('\t')
        out_hdrs = []
        for h in in_hdrs:
            if h in col_rename_dict:
                out_hdrs.append(col_rename_dict[h])
            else:
                out_hdrs.append(h)
        wf.write(",".join(out_hdrs))
        wf.write("\n")
        for line in rf:
            toks = line.strip().split('\t')
            wf.write(",".join(toks))
            wf.write("\n")

    ## convert CSV into PMtiles
    tippecanoe_cmd = " ".join([
        f"TIPPECANOE_MAX_THREADS={args.threads}",
        f"'{args.tippecanoe}'",
        f"-t {args.tmp_dir}",
        f"-o {args.out_prefix}.pmtiles",
        f"-Z {args.min_zoom} -z {args.max_zoom} --force",
        f"-s EPSG:3857 -M {args.max_tile_bytes} -O {args.max_feature_counts}",
        f"--drop-densest-as-needed",
        f"--extend-zooms-if-still-dropping",
        f"'--preserve-point-density-threshold={args.preserve_point_density_thres}'",
        f"--no-duplication",
        f"--no-clipping",
        f"--buffer 0",
        f"{args.out_prefix}.csv"
    ])

    logger.info(f"  * Running tippecanoe command: {tippecanoe_cmd}")
    result = subprocess.run(tippecanoe_cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        logger.error(f"Command {tippecanoe_cmd}\nfailed with error: {result.stderr.decode()}")
        sys.exit(1)
    else:
        logger.info("   * PMTiles creation command completed successfully")
    
    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

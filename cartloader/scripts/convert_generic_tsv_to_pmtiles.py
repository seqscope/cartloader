import sys, os, re, gzip, logging, argparse, inspect, subprocess, gzip
import pandas as pd
import numpy as np

from collections import Counter
from cartloader.utils.utils import create_custom_logger

def convert_generic_tsv_to_pmtiles(_args):
    """
    Convert a generic TSV file into PMTiles format
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Convert a generic TSV file into PMTiles format that can be converted to pmtiles")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-tsv', type=str, required=True, help='Genetic TSV file to convert to PMTiles. Typically named as ....tsv.gz, and should have spatial coordinates')
    inout_params.add_argument('--out-prefix', type=str, required=True, help='Prefix of output files. New directory will be created if needed')

    ## columns to add, remove, or rename
    iocol_params = parser.add_argument_group("Input/Output Columns Parameters", "Input/output column parameters .")
    iocol_params.add_argument('--remove-column', type=str, nargs='+', help='List of column names to remove')
    iocol_params.add_argument('--rename-column', type=str, nargs='+', help='List of columns to rename in the format of [old_name1:new_name1] [old_name2:new_name2] .... Note that lon/lat must exist in the output columns')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters frequently used by users")
    aux_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    aux_params.add_argument('--skip-pmtiles', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    aux_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')
    aux_params.add_argument('--in-delim', type=str, default='\t', help="Delimiter of the input file")
    aux_params.add_argument('--out-delim', type=str, default=',', help="Delimiter of the output file")
    aux_params.add_argument('--out-csv-suffix', type=str, default='.csv', help="Suffix of output CSV file")
    aux_params.add_argument('--out-pmtiles-suffix', type=str, default='.pmtiles', help="Suffix of output PMTiles file")
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level')
    aux_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level')
    aux_params.add_argument('--max-tile-bytes', type=int, default=1500000, help='Maximum bytes for each tile')
    aux_params.add_argument('--max-feature', type=int, default=200000, help='Max feature limits per tile')
    aux_params.add_argument('--preserve-point-density-thres', type=int, default=16384, help='Threshold for preserving point density for tippercanoe')

    args = parser.parse_args(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_convert_generic_tsv_to_pmtiles" + args.log_suffix if args.log else None)

    logger.info("Reading the metadata")

    ## parse the header information
    hdr_cols_input = []
    hdr_cols_output = []
    hdr_index_output = []
    with gzip.open(args.in_tsv, 'rt', encoding='utf-8') as rf:
        hdr_cols_input = rf.readline().rstrip().split(args.in_delim)
    
    if args.remove_column:
        set_cols_to_remove = set(args.remove_column)
        if len(set_cols_to_remove) != len(args.remove_column):
            ## identify duplicate column names from args.remove_column
            duplicates = list(filter(lambda x: Counter(args.remove_column)[x] > 1, set_cols_to_remove))
            logger.error(f"Duplicate column names {duplicates} in --remove-column")
            sys.exit(1)
    else:
        set_cols_to_remove = set()
    
    if args.rename_column:
        dict_cols_to_rename = dict([rename.split(':') for rename in args.rename_column])
        ## make sure that the keys are unique
        if len(dict_cols_to_rename) != len(args.rename_column):
            ## identify duplicate column names from args.remove_column
            rename_keys = [rename.split(':')[0] for rename in args.rename_column]
            duplicates = list(filter(lambda x: Counter(rename_keys)[x] > 1, dict_cols_to_rename.keys()))
            logger.error(f"Duplicate column names {duplicates} in --rename-column")
            sys.exit(1)
    else:
        dict_cols_to_rename = {}

    ## make sure that set_cols_to_remove and dict_cols_to_rename are disjoint
    if set_cols_to_remove & set(dict_cols_to_rename.keys()):
        ## identify overlapping columns
        overlapping_cols = set_cols_to_remove & set(dict_cols_to_rename.keys())
        logger.error(f"Overlapping column names {overlapping_cols} between --remove-column and --rename-column")
        sys.exit(1)

    n_removed = 0
    n_renamed = 0    
    for i, colname in enumerate(hdr_cols_input):
        if colname not in set_cols_to_remove:
            if colname in dict_cols_to_rename:
                hdr_cols_output.append(dict_cols_to_rename[colname])
                n_renamed += 1
            else:
                hdr_cols_output.append(colname)
            hdr_index_output.append(i)
        else:
            n_removed += 1

    ## make sure that the output columns are unique
    if len(hdr_cols_output) != len(set(hdr_cols_output)):
        ## identify duplicate column names
        duplicates = list(filter(lambda x: Counter(hdr_cols_output)[x] > 1, hdr_cols_output))
        logger.error(f"Duplicate column names {duplicates} in the output columns")
        sys.exit(1)

    if ( args.remove_column and len(args.remove_column) != n_removed ):
        ## identify columns that are not found
        missing_cols = set_cols_to_remove - set(hdr_cols_input)
        logger.error(f"Error in removing columns. The following columns are missing {missing_cols}")
        sys.exit(1)

    if ( args.rename_column and len(args.rename_column) != n_renamed ):
        ## identify columns that are not found
        missing_cols = dict_cols_to_rename - set(dict_cols_to_rename.keys())
        logger.error(f"Error in renaming columns. The following columns are missing {missing_cols}")
        sys.exit(1)
        
    logger.info("Converting the input data to CSV/TSV format")

    n_chunks = 0
    for chunk in pd.read_csv(args.in_tsv, sep=args.in_delim, chunksize=args.chunk_size, usecols=hdr_index_output):
        chunk.columns = hdr_cols_output

        ## write to a CSV file
        if n_chunks == 0:
            chunk.to_csv(f"{args.out_prefix}{args.out_csv_suffix}", sep=args.out_delim, index=False, header=True)
        else:
            chunk.to_csv(f"{args.out_prefix}{args.out_csv_suffix}", sep=args.out_delim, index=False, header=False, mode='a')

        n_chunks += 1
        logger.info(f"Finished processing chunk {n_chunks} of size {args.chunk_size}...")

    ## Converting the output to pmtiles
    if not args.skip_pmtiles:
        logger.info("Converting the CSV file to PMTiles format")
        cmd = f"{args.tippecanoe} -o {args.out_prefix}{args.out_pmtiles_suffix} -Z {args.min_zoom} -z {args.max_zoom} --force -s EPSG:3857 -M {args.max_tile_bytes} --drop-densest-as-needed --extend-zooms-if-still-dropping '--preserve-point-density-threshold={args.preserve_point_density_thres}' --no-duplication --no-clipping --no-tile-size-limit --buffer 0 {args.out_prefix}{args.out_csv_suffix}"
        print("Command to run:" + cmd)
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            logger.error("Error in converting the input TSV file into output CSV/TSV file")
            sys.exit(1)

    ## Remove the intermediate CSV file
    if not args.keep_intermediate_files:
        logger.info(f"Removing the intermediate CSV file {args.out_prefix}{args.out_csv_suffix}...")
        os.remove(f"{args.out_prefix}{args.out_csv_suffix}")

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
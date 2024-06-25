import sys, os, re, gzip, logging, argparse, inspect
import pandas as pd
import numpy as np
from collections import Counter

from cartloader.utils.utils import create_custom_logger

# Function to get log2 bins from the first file
def get_log2_bins(multiplier, in_features, out_prefix, out_features_suffix, delim, colname_feature, colname_count):
    df = pd.read_csv(in_features, sep=delim)
    df['log2_bin'] = np.floor(multiplier * np.log2(df[colname_count])).astype(int)
    log2_bins = df.set_index(colname_feature)['log2_bin'].to_dict()

    for log2_bin, group in df.groupby('log2_bin'):
        output_file = f"{out_prefix}.bin_{log2_bin}.{out_features_suffix}"
        group.drop('log2_bin', axis=1).to_csv(output_file, sep='\t', index=False)
    
    return log2_bins

# Function to process chunks of the second file
def process_chunk(chunk, log2_bins, colname_feature, out_prefix, out_tsv_suffix, out_tsv_delim, rename_dict):
    bin2nmols = {}
    for log2_bin, group in chunk.groupby(chunk[colname_feature].map(log2_bins)):
        bin2nmols[log2_bin] = group.shape[0]
        output_file = f"{out_prefix}.bin_{log2_bin}.{out_tsv_suffix}"
        mode = 'a' if os.path.exists(output_file) else 'w'
        group.rename(columns=rename_dict).to_csv(output_file, sep=out_tsv_delim, index=False, mode=mode, header=(mode == 'w'))
    return bin2nmols

def split_molecule_counts(_args):
    """
    Split the cross-platform TSV/CSV file by groups based on total number of molecular observations
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Split the cross-platform TSV/CSV file by groups based on total number of expression")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-molecules', type=str, help='Input Long Format TSV/CSV (possibly gzipped) file containing the X/Y coordinates and gene expression counts per spot')
    inout_params.add_argument('--in-features', type=str, help='Input TSV/CSV (possibly gzipped) file containing the gene name and total count for each gene')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='The output prefix. New directory will be created if needed')  
    inout_params.add_argument('--in-molecules-delim', type=str, default='\t', help='Delimiter used in the input molecules files. Default is tab.')  
    inout_params.add_argument('--in-features-delim', type=str, default='\t', help='Delimiter used in the input feature files. Default is tab.')  
    inout_params.add_argument('--out-molecules-delim', type=str, default='\t', help='Delimiter used in the output molecule TSV/CSV files. Default is ,')
    inout_params.add_argument('--out-features-delim', type=str, default='\t', help='Delimiter used in the output feature files. Default is tab.')
    inout_params.add_argument('--out-molecules-suffix', type=str, default="molecules.tsv.gz", help='The output file to store individual molecule count matrix. Each file name will be [out-prefix].split.[bin_id].[out-molecules-suffix]. Default: molecules.tsv.gz')
    inout_params.add_argument('--out-features-suffix', type=str, default="features.tsv.gz", help='The output file to store feature count matrix. Each file name will be [out-prefix].split.[bin_id].[out-features-suffix]. Default: features.tsv.gz')

    iocol_params = parser.add_argument_group("Input/Output Columns Parameters", "Input/output column parameters .")
    iocol_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    iocol_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    iocol_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--bin-multiplier', type=float, default=1.0, help='Multiplier used to determine the bin. Bin is determined as [multiplier] * log2(total_count). Default is 1.0')
    key_params.add_argument('--dummy-genes', type=str, default='', help='A single name or a regex describing the names of negative control probes')
    key_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    key_params.add_argument('--skip-original', action='store_true', default=False, help='Skip writing the original file')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    args = parser.parse_args(_args)

    logger = create_custom_logger(__name__, args.out_prefix + args.log_suffix if args.log else None)

    logger.info("Reading the feature counts and splitting into bins")
    # Get log2 bins from the tsv file
    log2_bins = get_log2_bins(args.bin_multiplier, args.in_features, args.out_prefix, args.out_features_suffix, args.in_features_delim, args.colname_feature, args.colname_count)

    # set the column names to be renamed
    rename_dict = {}
    if args.col_rename:
        for rename in args.col_rename:
            old_name, new_name = rename.split(':')
            rename_dict[old_name] = new_name

    logger.info("Splitting the TSV file by chunks and writing to individual bins")  
    # Process the second file in chunks
    nchunks = 0
    bin2mols = {}
    for chunk in pd.read_csv(args.in_molecules, sep=args.in_molecules_delim, chunksize=args.chunk_size):
        bin2mols_chunk = process_chunk(chunk, log2_bins, args.colname_feature, args.out_prefix, args.out_molecules_suffix, args.out_molecules_delim, rename_dict)
        for bin_id, nmols in bin2mols_chunk.items():
            bin2mols[bin_id] = bin2mols.get(bin_id, 0) + nmols
            bin2mols["all"] = bin2mols.get("all", 0) + nmols
        if not args.skip_original:
            chunk.rename(columns=rename_dict).to_csv(f"{args.out_prefix}.bin_all.{args.out_molecules_suffix}", sep=args.out_molecules_delim, index=False, mode='a', header=(nchunks == 0))
        nchunks += 1
        logger.info(f"Finished processing chunk {nchunks} of size {args.chunk_size}...")

    bin2nftrs = Counter(log2_bins.values())
    bin2nftrs["all"] = len(log2_bins)

    ## write the manifest file
    logger.info("Writing the index file...")  
    with open(f"{args.out_prefix}.index.tsv", 'w') as wf:
        wf.write(f"bin_id\tmolecules_count\tfeatures_count\tmolecules_path\tfeatures_path\n")
        out_basename = os.path.basename(args.out_prefix)
        if not args.skip_original:
            nmols = bin2mols.get("all", 0)
            nftrs = bin2nftrs.get("all", 0)
            wf.write(f"all\t{nmols}\t{nftrs}\t{out_basename}.bin_all.{args.out_molecules_suffix}\t{out_basename}.bin_all.{args.out_features_suffix}\n")
        for bin_id in sorted(list(set(log2_bins.values()))):
            nmols = bin2mols.get(bin_id, 0)
            nftrs = bin2nftrs.get(bin_id, 0)
            wf.write(f"{bin_id}\t{nmols}\t{nftrs}\t{out_basename}.bin_{bin_id}.{args.out_molecules_suffix}\t{out_basename}.bin_{bin_id}.{args.out_features_suffix}\n")

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
import sys, os, re, gzip, logging, argparse, inspect
import pandas as pd
import numpy as np
from collections import Counter

from cartloader.utils.utils import create_custom_logger

# Function to get equal bins from the first file
def get_equal_bins(n_bins, in_features, out_prefix, out_features_suffix, delim, colname_feature, colname_count, skip_original):
    ## read the feature data frame
    df = pd.read_csv(in_features, sep=delim)

    ## sort by smallest to largest
    df = df.sort_values(by=colname_count, ascending=False).reset_index(drop=True)

    # Calculate the cumulative sum of 'count'
    df['cumulative_sum'] = df[colname_count].cumsum()

    # Define the total sum and target sum for each bin
    total_sum = df[colname_count].sum()
    target_sum = total_sum / n_bins
    
    df['bin'] = 0
    current_sum = 0
    current_bin = 1
    for index, row in df.iterrows():
#        print("***", total_sum, target_sum, current_sum, current_bin, row[colname_count])
        if current_sum + row[colname_count] < target_sum: # keep adding more features
            current_sum += row[colname_count]
            df.at[index, 'bin'] = current_bin
        else: # stop here
            current_sum += row[colname_count]
            df.at[index, 'bin'] = current_bin
            current_bin += 1
            ## update target_sum
            total_sum -= current_sum
            if current_bin <= n_bins:
              target_sum = total_sum / (n_bins - current_bin + 1)
              current_sum = 0
    

    equal_bins = df.set_index(colname_feature)['bin'].to_dict()

    df = df.drop(columns=['cumulative_sum'])

    if not skip_original:
        df.to_csv(f"{out_prefix}_all_{out_features_suffix}", sep='\t', index=False, na_rep='NA')
    df.to_json(f"{out_prefix}_bin_counts.json", orient='records')

    for equal_bin, group in df.groupby('bin'):
        output_file = f"{out_prefix}_bin{equal_bin}_{out_features_suffix}"
        group.drop('bin', axis=1).to_csv(output_file, sep='\t', index=False, na_rep='NA')
    
    return equal_bins

# Function to get log2 bins from the first file
def get_log2_bins(multiplier, in_features, out_prefix, out_features_suffix, delim, colname_feature, colname_count, skip_original):
    #print(f"delim = {delim} {len(delim)}")
    df = pd.read_csv(in_features, sep=delim)
    df['bin'] = np.floor(multiplier * np.log2(df[colname_count])).astype(int)
    log2_bins = df.set_index(colname_feature)['bin'].to_dict()

    if not skip_original:
        df.to_csv(f"{out_prefix}_all_{out_features_suffix}", sep='\t', index=False)
    df.to_json(f"{out_prefix}_bin_counts.json", orient='records')

    for log2_bin, group in df.groupby('bin'):
        output_file = f"{out_prefix}_bin{log2_bin}_{out_features_suffix}"
        group.drop('bin', axis=1).to_csv(output_file, sep='\t', index=False)
    
    return log2_bins

# Function to process chunks of the second file
def process_chunk(chunk, bins, colname_feature, out_prefix, out_tsv_suffix, out_tsv_delim, rename_dict):
    bin2nmols = {}
    chunk.rename(columns=rename_dict, inplace=True)
    for bin_id, group in chunk.groupby(chunk[colname_feature].map(bins)):
        bin2nmols[bin_id] = group.shape[0]
        # bin_id may be 1.0, 2.0, etc. Convert to int and then to str
        output_file = f"{out_prefix}_bin{str(int(bin_id))}_{out_tsv_suffix}"
        mode = 'a' if os.path.exists(output_file) else 'w'
        #group.rename(columns=rename_dict).to_csv(output_file, sep=out_tsv_delim, index=False, mode=mode, header=(mode == 'w'), na_rep='NA')
        group.to_csv(output_file, sep=out_tsv_delim, index=False, mode=mode, header=(mode == 'w'), na_rep='NA')
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
    inout_params.add_argument('--out-molecules-suffix', type=str, default="molecules.tsv.gz", help='The output file to store individual molecule count matrix. Each file name will be [out-prefix]_[bin_id]_[out-molecules-suffix]. Default: molecules.tsv.gz')
    inout_params.add_argument('--out-features-suffix', type=str, default="features.tsv.gz", help='The output file to store feature count matrix. Each file name will be [out-prefix]_[bin_id]_[out-features-suffix]. Default: features.tsv.gz')

    iocol_params = parser.add_argument_group("Input/Output Columns Parameters", "Input/output column parameters .")
    iocol_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    iocol_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    iocol_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--equal-bins', action='store_true', default=False, help='Use equal bins instead of log2 bins')
    key_params.add_argument('--bin-count', type=int, default=50, help='When --equal-bins is used, determine the number of bins to split the data into the same bin')
    key_params.add_argument('--log2-multiplier', type=float, default=1.0, help='Multiplier used to determine the log2 bin. Bin is determined as [multiplier] * log2(total_count). Default is 1.0')
    key_params.add_argument('--dummy-genes', type=str, default='', help='A single name or a regex describing the names of negative control probes')
    key_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    key_params.add_argument('--skip-original', action='store_true', default=False, help='Skip writing the original file')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    args = parser.parse_args(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_split_molecule_counts" + args.log_suffix if args.log else None)

    logger.info("Reading the feature counts and splitting into bins")
    # Get log2 bins from the tsv file
    if args.equal_bins:
        bins = get_equal_bins(args.bin_count, args.in_features, args.out_prefix, args.out_features_suffix, args.in_features_delim, args.colname_feature, args.colname_count, args.skip_original)
    else:
        bins = get_log2_bins(args.log2_multiplier, args.in_features, args.out_prefix, args.out_features_suffix, args.in_features_delim, args.colname_feature, args.colname_count, args.skip_original)

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
    for chunk in pd.read_csv(args.in_molecules, sep=args.in_molecules_delim, chunksize=args.chunk_size, keep_default_na=False, na_values=[]):
        bin2mols_chunk = process_chunk(chunk, bins, args.colname_feature, args.out_prefix, args.out_molecules_suffix, args.out_molecules_delim, rename_dict)
        for bin_id, nmols in bin2mols_chunk.items():
            bin2mols[bin_id] = bin2mols.get(bin_id, 0) + nmols
            bin2mols["all"] = bin2mols.get("all", 0) + nmols
        if not args.skip_original:
            chunk.rename(columns=rename_dict).to_csv(f"{args.out_prefix}_all_{args.out_molecules_suffix}", sep=args.out_molecules_delim, index=False, mode='a', header=(nchunks == 0),na_rep='NA')
        nchunks += 1
        logger.info(f"Finished processing chunk {nchunks} of size {args.chunk_size}...")

    bin2nftrs = Counter(bins.values())
    bin2nftrs["all"] = len(bins)

    ## write the index file
    logger.info("Writing the index file...")  
    with open(f"{args.out_prefix}_index.tsv", 'w') as wf:
        wf.write(f"bin_id\tmolecules_count\tfeatures_count\tmolecules_path\tfeatures_path\n")
        out_basename = os.path.basename(args.out_prefix)
        if not args.skip_original:
            nmols = bin2mols.get("all", 0)
            nftrs = bin2nftrs.get("all", 0)
            wf.write(f"all\t{nmols}\t{nftrs}\t{out_basename}_all_{args.out_molecules_suffix}\t{out_basename}_all_{args.out_features_suffix}\n")
        for bin_id in sorted(list(set(bins.values()))):
            nmols = bin2mols.get(bin_id, 0)
            nftrs = bin2nftrs.get(bin_id, 0)
            wf.write(f"{bin_id}\t{nmols}\t{nftrs}\t{out_basename}_bin{bin_id}_{args.out_molecules_suffix}\t{out_basename}_bin{bin_id}_{args.out_features_suffix}\n")

    ## write json file of 

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
import sys, os, gzip, logging, argparse, inspect, json
import polars as pl
import numpy as np
from collections import Counter

from cartloader.utils.utils import create_custom_logger


def _write_csv(df, path, separator='\t', null_value=''):
    """Write a polars DataFrame to CSV, supporting gzip compression."""
    if path.endswith('.gz'):
        csv_bytes = df.write_csv(separator=separator, null_value=null_value).encode('utf-8')
        with gzip.open(path, 'wb') as f:
            f.write(csv_bytes)
    else:
        df.write_csv(path, separator=separator, null_value=null_value)


# Function to get equal bins from the first file
def get_equal_bins(n_bins, in_features, out_prefix, out_features_suffix, delim, colname_feature, colname_count, skip_original):
    ## read the feature data frame
    df = pl.read_csv(in_features, separator=delim)

    ## sort by largest to smallest
    df = df.sort(colname_count, descending=True)

    # Fast bin assignment using numpy arrays instead of iterrows
    counts = df[colname_count].to_numpy().astype(np.float64)
    n = len(counts)
    bins_arr = np.zeros(n, dtype=np.int64)

    total_sum = float(counts.sum())
    target_sum = total_sum / n_bins
    current_sum = 0.0
    current_bin = 1

    for i in range(n):
        c = counts[i]
        if current_sum + c < target_sum or current_bin >= n_bins:
            current_sum += c
            bins_arr[i] = current_bin
        else:
            current_sum += c
            bins_arr[i] = current_bin
            current_bin += 1
            total_sum -= current_sum
            if current_bin <= n_bins:
                target_sum = total_sum / (n_bins - current_bin + 1)
                current_sum = 0.0

    df = df.with_columns(pl.Series("bin", bins_arr))

    equal_bins = dict(zip(df[colname_feature].to_list(), df["bin"].to_list()))

    if not skip_original:
        _write_csv(df, f"{out_prefix}_all_{out_features_suffix}", separator='\t', null_value='NA')

    with open(f"{out_prefix}_bin_counts.json", 'w') as f:
        json.dump(df.to_dicts(), f, separators=(',', ':'))

    for bin_val in sorted(df["bin"].unique().to_list()):
        group = df.filter(pl.col("bin") == bin_val).drop("bin")
        output_file = f"{out_prefix}_bin{bin_val}_{out_features_suffix}"
        _write_csv(group, output_file, separator='\t', null_value='NA')

    return equal_bins

# Function to get log2 bins from the first file
def get_log2_bins(multiplier, in_features, out_prefix, out_features_suffix, delim, colname_feature, colname_count, skip_original):
    df = pl.read_csv(in_features, separator=delim)
    df = df.with_columns(
        (pl.col(colname_count).cast(pl.Float64) + 1).log(2.0).mul(multiplier).floor().cast(pl.Int64).alias("bin")
    )
    log2_bins = dict(zip(df[colname_feature].to_list(), df["bin"].to_list()))

    if not skip_original:
        _write_csv(df, f"{out_prefix}_all_{out_features_suffix}", separator='\t')

    with open(f"{out_prefix}_bin_counts.json", 'w') as f:
        json.dump(df.to_dicts(), f, separators=(',', ':'))

    for bin_val in sorted(df["bin"].unique().to_list()):
        group = df.filter(pl.col("bin") == bin_val).drop("bin")
        output_file = f"{out_prefix}_bin{bin_val}_{out_features_suffix}"
        _write_csv(group, output_file, separator='\t')

    return log2_bins

def split_molecule_counts(_args):
    """
    Split the cross-platform TSV/CSV file by groups based on total number of molecular observations
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Split the cross-platform TSV/CSV file by groups based on total number of expression")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-molecules', type=str, help='Input Long Format TSV/CSV (possibly gzipped) file containing the X/Y coordinates and gene expression counts per spot')
    inout_params.add_argument('--in-features', type=str, help='Input TSV/CSV (possibly gzipped) file containing the gene name and total count for each gene')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='Output prefix. New directory will be created if needed')
    inout_params.add_argument('--in-molecules-delim', type=str, default='\t', help='Delimiter used in the input molecules files. Default is tab.')
    inout_params.add_argument('--in-features-delim', type=str, default='\t', help='Delimiter used in the input feature files. Default is tab.')
    inout_params.add_argument('--out-molecules-delim', type=str, default='\t', help='Delimiter used in the output molecule TSV/CSV files. Default is ,')
    inout_params.add_argument('--out-features-delim', type=str, default='\t', help='Delimiter used in the output feature files. Default is tab.')
    inout_params.add_argument('--out-molecules-suffix', type=str, default="molecules.tsv.gz", help='Name prefix of a output file to store individual molecule count matrix. Each file name will be [out-prefix]_[bin_id]_[out-molecules-suffix]. Default: molecules.tsv.gz')
    inout_params.add_argument('--out-features-suffix', type=str, default="features.tsv.gz", help='Name prefix of a output file to store feature count matrix. Each file name will be [out-prefix]_[bin_id]_[out-features-suffix]. Default: features.tsv.gz')

    iocol_params = parser.add_argument_group("Input/Output Columns Parameters", "Input/output column parameters .")
    iocol_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    iocol_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    iocol_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--equal-bins', action='store_true', default=False, help='Use equal bins instead of log2 bins')
    key_params.add_argument('--bin-count', type=int, default=50, help='When --equal-bins is used, determine the number of bins to split the data into the same bin')
    key_params.add_argument('--log2-multiplier', type=float, default=1.0, help='Multiplier used to determine the log2 bin. Bin is determined as [multiplier] * log2(total_count). Default is 1.0')
    key_params.add_argument('--dummy-genes', type=str, default='', help='A single name or a regex describing the names of negative control probes')
    key_params.add_argument('--chunk-size', type=int, default=1000000, help='Kept for backward compatibility. Polars reads the full file at once for efficiency.')
    key_params.add_argument('--skip-original', action='store_true', default=False, help='Skip writing the original file')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    args = parser.parse_args(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_split_molecule_counts" + args.log_suffix if args.log else None)

    logger.info("Reading the feature counts and splitting into bins")
    # Get bins from the feature file
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

    # Create bins lookup DataFrame for efficient join-based bin mapping
    bins_df = pl.DataFrame({
        args.colname_feature: list(bins.keys()),
        "__bin__": [int(v) for v in bins.values()],
    })

    logger.info("Reading the molecules file...")
    df_mol = pl.read_csv(args.in_molecules, separator=args.in_molecules_delim)

    if rename_dict:
        df_mol = df_mol.rename(rename_dict)

    logger.info(f"Read {df_mol.height} molecules. Joining with bin assignments and splitting...")

    # Join with bins lookup for vectorized bin mapping (replaces per-row dict.map)
    df_mol = df_mol.join(bins_df, on=args.colname_feature, how="left")

    # Write per-bin output files
    bin2mols = {}
    for group in df_mol.partition_by("__bin__", maintain_order=True):
        bin_id_val = group["__bin__"][0]
        if bin_id_val is None:
            continue
        bin_id = int(bin_id_val)
        group = group.drop("__bin__")
        bin2mols[bin_id] = group.height
        bin2mols["all"] = bin2mols.get("all", 0) + group.height
        output_file = f"{args.out_prefix}_bin{bin_id}_{args.out_molecules_suffix}"
        _write_csv(group, output_file, separator=args.out_molecules_delim, null_value='NA')

    if not args.skip_original:
        _write_csv(
            df_mol.drop("__bin__"),
            f"{args.out_prefix}_all_{args.out_molecules_suffix}",
            separator=args.out_molecules_delim,
            null_value='NA',
        )

    logger.info(f"Finished splitting molecules into {len(bin2mols) - (1 if 'all' in bin2mols else 0)} bins")

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

    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
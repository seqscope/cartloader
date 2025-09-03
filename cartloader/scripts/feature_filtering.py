import argparse, os, re, sys, logging, gzip, inspect
import pandas as pd
import numpy as np

from cartloader.utils.utils import create_custom_logger, log_info

def filter_feature_by_type(include_feature_type_regex, feature_type_ref, feature_type_ref_delim="\t",
                           feature_type_ref_colname_name=None, feature_type_ref_colname_type=None,
                           feature_type_ref_colidx_name=None, feature_type_ref_colidx_type=None,
                           logger=None, chunksize=50000):
    assert feature_type_ref is not None, "Provide --feature-type-ref when --include-feature-type-regex is enabled."
    assert (feature_type_ref_colidx_name is not None and feature_type_ref_colidx_type is not None) or (feature_type_ref_colname_name is not None and feature_type_ref_colname_type is not None), "Provide either --feature-type-ref-colidx-name and --feature-type-ref-colidx-type or --feature-type-ref-colname-name and --feature-type-ref-colname-type"

    ftrlist_keep_type = []
    n_ftr = 0
    # log
    log_info(f"  - include_feature_type_regex: {include_feature_type_regex} using --feature-type-ref {feature_type_ref} ", logger)
    if feature_type_ref_colname_name is not None and feature_type_ref_colname_type is not None:
        log_info(f"  - --feature-type-ref-colname-name: {feature_type_ref_colname_name} and --feature-type-ref-colname-type: {feature_type_ref_colname_type}", logger)
    elif feature_type_ref_colidx_name is not None and feature_type_ref_colidx_type is not None:
        log_info(f"  - --feature-type-ref-colidx-name: {feature_type_ref_colidx_name} and --feature-type-ref-colidx-type: {feature_type_ref_colidx_type}", logger)
    # processing
    for chunk in pd.read_csv(feature_type_ref, chunksize=chunksize, index_col=None, sep=feature_type_ref_delim, comment="#"):
        if feature_type_ref_colname_name is not None and feature_type_ref_colname_type is not None:
            chunk = chunk[[feature_type_ref_colname_name, feature_type_ref_colname_type]]
        elif feature_type_ref_colidx_name is not None and feature_type_ref_colidx_type is not None:
            chunk = chunk.iloc[:, [feature_type_ref_colidx_name, feature_type_ref_colidx_type]]
        chunk.columns = ["ftr", "ftr_type"]
        n_ftr += chunk.shape[0]
        ftrlist_keep_type.extend(chunk[chunk["ftr_type"].str.contains(include_feature_type_regex, regex=True)]["ftr"].tolist())
    ftrlist_keep_type = list(set(ftrlist_keep_type))
    log_info(f"  - N features in the reference file: {n_ftr}", logger)
    log_info(f"  - N features with the specified feature type: {len(ftrlist_keep_type)}", logger)
    return ftrlist_keep_type

def map_filter_per_feature(df_ftrinfo, 
                           include_feature_list=None, exclude_feature_list=None, 
                           include_feature_substr=None, exclude_feature_substr=None, 
                           include_feature_regex=None, exclude_feature_regex=None,
                           include_feature_type_regex=None, feature_type_ref=None, feature_type_ref_delim='\t', feature_type_ref_colname_name=None, feature_type_ref_colname_type=None, feature_type_ref_colidx_name=None, feature_type_ref_colidx_type=None,
                           logger=None, chunksize=50000):
    """
    df_ftrinfo: pd.DataFrame with a column named "feature" containing feature names.
    """
    df_ftrinfo.drop_duplicates(inplace=True)
    n_raw_ftr = df_ftrinfo.shape[0]
    df_ftrinfo.columns = ["feature"]

    # ====
    # 1. Examining each feature for filtering criteria
    # ====

    # Filtering based on provided lists
    log_info(" Examining each feature for filtering criteria", logger) 
    if include_feature_list is not None:
        log_info(f"  - include_feature_list: {include_feature_list}", logger)
        ftr_keep_list = pd.read_csv(include_feature_list, header=None, names=["feature"])['feature'].tolist()
        df_ftrinfo["rm_by_include_feature_list"] = df_ftrinfo["feature"].apply(lambda x: "include_feature_list" if x not in ftr_keep_list else np.nan)

    if exclude_feature_list is not None:
        log_info(f"  - exclude_feature_list: {exclude_feature_list}", logger)
        ftr_exclude_list = pd.read_csv(exclude_feature_list, header=None, names=["feature"])['feature'].tolist()
        df_ftrinfo["rm_by_exclude_feature_list"] = df_ftrinfo["feature"].apply(lambda x: "exclude_feature_list" if x in ftr_exclude_list else np.nan)
        
    # Filtering based on substrings and regex patterns
    filters = {
            "include_feature_substr": (include_feature_substr, True, False),
            "exclude_feature_substr": (exclude_feature_substr, False, False),
            "include_feature_regex": (include_feature_regex, True, True),
            "exclude_feature_regex": (exclude_feature_regex, False, True),
        }

    for filter_name, (filter_value, include, is_regex) in filters.items():
            if filter_value:  # Skip if no filtering value is provided
                log_info(f"  - Applying {filter_name}: {filter_value}", logger)
                # Ensure regex pattern does not introduce capture groups
                pattern = filter_value if is_regex else re.escape(filter_value)
                # Create a boolean mask based on the filtering condition
                mask = df_ftrinfo["feature"].str.contains(pattern, regex=is_regex, na=False) ^ include
                # Assign filtering reason where the condition is met
                df_ftrinfo.loc[mask, f"rm_by_{filter_name}"] = filter_name

    # Feature type filtering
    if include_feature_type_regex is not None:
        ftrlist_keep_type= filter_feature_by_type(include_feature_type_regex, feature_type_ref, 
                                                  feature_type_ref_delim=feature_type_ref_delim, 
                                                  feature_type_ref_colname_name=feature_type_ref_colname_name, 
                                                  feature_type_ref_colname_type=feature_type_ref_colname_type, 
                                                  feature_type_ref_colidx_name=feature_type_ref_colidx_name, 
                                                  feature_type_ref_colidx_type=feature_type_ref_colidx_type, 
                                                  logger=logger, chunksize=chunksize)
        df_ftrinfo["rm_by_include_feature_type"] = df_ftrinfo["feature"].apply(lambda x: "include_feature_type_regex" if x not in ftrlist_keep_type else np.nan)

    # ====
    # 2. Summarize the filtering criteria & N removal features
    # ====
    log_info(" Summarizing filtering per feature ...", logger)
    log_info(f"  - N raw features: {n_raw_ftr}", logger)
    comment_cols = [x for x in df_ftrinfo.columns if x.startswith("rm_by_")]
    for cmtcol in comment_cols:
        rm_ftr_reason = cmtcol.replace("rm_by_", "Number of features will be removed by --").replace("_", "-")
        df_ftrinfo[cmtcol] = df_ftrinfo[cmtcol].replace("nan", np.nan)
        df_ftrinfo[cmtcol] = df_ftrinfo[cmtcol].replace({pd.NA: np.nan}).dropna()
        rm_ftr_count = df_ftrinfo[cmtcol].dropna().shape[0]
        log_info(f"  - {rm_ftr_reason}: {rm_ftr_count}", logger)

    # df_ftrinfo 
    df_ftrinfo["filtering"] = df_ftrinfo[comment_cols].apply(lambda x: ";".join(x.dropna()) if len(x.dropna()) > 1 else x.dropna().values[0] if len(x.dropna()) == 1 else None, axis=1)
    df_ftrinfo = df_ftrinfo[["feature", "filtering"]]
    
    log_info(f'  - N keep features: {df_ftrinfo[df_ftrinfo["filtering"].isna()].shape[0]}', logger)
    return df_ftrinfo


def feature_filtering(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description=(
            "Filter features (e.g., genes) by lists, substrings, regex, or type. "
            "Writes an annotated record and/or a filtered CSV/TSV."
        ),
    )
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output paths")
    inout_params.add_argument('--in-csv', type=str, required=True, help='Path to input CSV/TSV containing feature names')
    inout_params.add_argument('--out-csv', type=str, required=True, help='Path to compressed output CSV/TSV of rows kept after filtering (writes gzip)')
    inout_params.add_argument('--out-record', type=str, default=None, help='Path to output TSV of feature and filtering status (default: None)')
    inout_params.add_argument('--csv-delim', type=str, default='\t', help='Delimiter for input/output CSV/TSV (default: tab)')
    inout_params.add_argument('--csv-colname-feature-name', type=str, default="gene", help='Column name of feature in input CSV/TSV (default: gene)')
    inout_params.add_argument('--chunksize', type=int, default=50000, help='Chunk size for streaming CSV reads (default: 50000)')
    inout_params.add_argument('--log', action='store_true', default=False, help='Enable logging to a .log file next to outputs (default: False)')

    # Filtering parameters (provide at least one criterion)
    list_params = parser.add_argument_group("List Filters", "Include/exclude by explicit feature lists")
    list_params.add_argument('--include-feature-list', type=str, default=None, help='Path to file with feature names/IDs to include (default: None)')
    list_params.add_argument('--exclude-feature-list', type=str, default=None, help='Path to file with feature names/IDs to exclude (default: None)')

    substr_params = parser.add_argument_group("Substring Filters", "Include/exclude by substring match")
    substr_params.add_argument('--include-feature-substr', type=str, default=None, help='Substring of feature names to include (default: None)')
    substr_params.add_argument('--exclude-feature-substr', type=str, default=None, help='Substring of feature names to exclude (default: None)')

    regex_params = parser.add_argument_group("Regex Filters", "Include/exclude by regex match")
    regex_params.add_argument('--include-feature-regex', type=str, default=None, help='Regex of feature names to include (default: None)')
    regex_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude (default: None)')

    type_params = parser.add_argument_group("Type Filters", "Filter by feature/gene type using a reference file; provide --feature-type-ref and specify name/type columns via --feature-type-ref-colname-* or --feature-type-ref-colidx-*")
    type_params.add_argument('--include-feature-type-regex', type=str, default=None, help='Regex of feature type to include (default: None)')
    type_params.add_argument('--feature-type-ref', type=str, default=None, help='Path to a reference file with feature name and type columns (default: None)')
    type_params.add_argument('--feature-type-ref-delim', type=str, default='\t', help='Delimiter for reference file (default: tab)')
    type_params.add_argument('--feature-type-ref-colname-name', type=str, default=None, help='Column name for feature name in reference file (default: None)')
    type_params.add_argument('--feature-type-ref-colname-type', type=str, default=None, help='Column name for feature type in reference file (default: None)')
    type_params.add_argument('--feature-type-ref-colidx-name', type=int, default=None, help='0-based column index for feature name in reference file (default: None)')
    type_params.add_argument('--feature-type-ref-colidx-type', type=int, default=None, help='0-based column index for feature type in reference file (default: None)')

    args = parser.parse_args(_args)

    # sanity check 
    assert os.path.exists(args.in_csv), f"Input file does not exist: {args.in_csv}"
    assert any([
        args.include_feature_list,
        args.exclude_feature_list,
        args.include_feature_substr,
        args.exclude_feature_substr,
        args.include_feature_regex,
        args.exclude_feature_regex,
        args.include_feature_type_regex,
    ]), "At least one filtering criterion is required."

    # output
    out_dir = os.path.dirname(args.out_csv) or '.'
    os.makedirs(out_dir, exist_ok=True)

    out_log = os.path.join(out_dir, "feature_filtering.log")
    logger = create_custom_logger(__name__, out_log if args.log else None)

    # ===============================
    #  Annotate Feature filtering information
    # ===============================

    # Read all input feature names
    logger.info(" Reading input --in-csv: "+args.in_csv)
    assert os.path.exists(args.in_csv), f"File not found: {args.in_csv} (--in-csv)"
    df_ftrinfo = pd.read_csv(args.in_csv, usecols=[args.csv_colname_feature_name], sep=args.csv_delim, header=0, comment="#")
    df_ftrinfo = map_filter_per_feature(df_ftrinfo,
                                        include_feature_list=args.include_feature_list, exclude_feature_list=args.exclude_feature_list,
                                        include_feature_substr=args.include_feature_substr, exclude_feature_substr=args.exclude_feature_substr,
                                        include_feature_regex=args.include_feature_regex, exclude_feature_regex=args.exclude_feature_regex,
                                        include_feature_type_regex=args.include_feature_type_regex, feature_type_ref=args.feature_type_ref, feature_type_ref_delim=args.feature_type_ref_delim, 
                                        feature_type_ref_colname_name=args.feature_type_ref_colname_name, feature_type_ref_colname_type=args.feature_type_ref_colname_type, feature_type_ref_colidx_name=args.feature_type_ref_colidx_name, feature_type_ref_colidx_type=args.feature_type_ref_colidx_type,
                                        logger=logger, chunksize=args.chunksize)

    # Save filtering records
    if args.out_record is not None:
        df_ftrinfo.to_csv(args.out_record, sep=args.csv_delim, index=False)
        logger.info(f" - Saved feature filtering records to: {args.out_record}")

    # ===============================
    #  filtering
    # ===============================
    # std out_csv

    if not args.out_csv.endswith('.gz'):
        logger.info("--out-csv does not end with .gz; appending .gz to match gzip compression")
        args.out_csv = args.out_csv + '.gz'

    if os.path.exists(args.out_csv):
        logger.info(f"  - Existing --out-csv found; overwriting: {args.out_csv}")
        os.remove(args.out_csv)

    # read-in
    logger.info("Applying filtering...")
    df_ftrinfo = df_ftrinfo[df_ftrinfo["filtering"].isna()]
    
    n_raw_csv = 0
    n_filtered_csv = 0
    for chunk in pd.read_csv(args.in_csv, chunksize=args.chunksize, index_col=None, sep=args.csv_delim, comment="#"):
        n_raw_csv += chunk.shape[0]
        chunk = chunk[chunk[args.csv_colname_feature_name].isin(df_ftrinfo["feature"])]
        n_filtered_csv += chunk.shape[0]
        # write down
        if not os.path.exists(args.out_csv):
            chunk.to_csv(args.out_csv, sep=args.csv_delim, index=False, compression='gzip')
        else:
            chunk.to_csv(args.out_csv, sep=args.csv_delim, index=False, compression='gzip', mode='a', header=False)
    logger.info(f'  - N raw rows in --in-csv: {n_raw_csv}. N row may differ from N features due to the scenario of same feature name with different feature IDs')
    logger.info(f'  - N filtered rows in --out-csv: {n_filtered_csv}')
    logger.info(f"  - Saved filtered CSV to: {args.out_csv}")

if __name__ == "__main__":
    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], os.path.splitext(os.path.basename(__file__))[0])
    func(sys.argv[1:])

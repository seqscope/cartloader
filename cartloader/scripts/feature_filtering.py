import argparse, os, re, sys, logging, gzip, inspect
import pandas as pd
import numpy as np

from cartloader.utils.utils import create_custom_logger, log_info
from cartloader.scripts.feature_filtering_by_type import filter_feature_by_type

def map_filter_per_feature(df_ftrinfo, 
                           include_feature_list=None, exclude_feature_list=None, 
                           include_feature_substr=None, exclude_feature_substr=None, 
                           include_feature_regex=None, exclude_feature_regex=None,
                           include_feature_type_regex=None, feature_type_ref=None, feature_type_ref_delim='\t', feature_type_ref_colname_name=None, feature_type_ref_colname_type=None, feature_type_ref_colidx_name=None, feature_type_ref_colidx_type=None,
                           logger=None):
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
        ftrlist_keep_type= filter_feature_by_type(include_feature_type_regex, feature_type_ref, feature_type_ref_delim=feature_type_ref_delim, 
                                                  feature_type_ref_colname_name=feature_type_ref_colname_name, feature_type_ref_colname_type=feature_type_ref_colname_type, 
                                                  feature_type_ref_colidx_name=feature_type_ref_colidx_name, feature_type_ref_colidx_type=feature_type_ref_colidx_type, 
                                                  logger=logger)
        df_ftrinfo["rm_by_include_feature_type"] = df_ftrinfo["feature"].apply(lambda x: "include_feature_type_regex" if x not in ftrlist_keep_type else np.nan)

    # ====
    # 2. Summarize the filtering criteria & N removal genes
    # ====
    log_info(" Summarizing filtering per gene ...", logger)
    log_info(f"  - N raw features: {n_raw_ftr}", logger)
    comment_cols = [x for x in df_ftrinfo.columns if x.startswith("rm_by_")]
    for cmtcol in comment_cols:
        rm_ftr_reason = cmtcol.replace("rm_by_", "Number of genes will be removed by --").replace("_", "-")
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
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                        description="""
                                        This script filters features based on specified criteria. 
                                        If --out-record is specified, the script saves the features along with their filtering status in an output record file.
                                        If --out-csv is specified, the script applies filtering and saves the remaining features in a compressed output file.
                                        """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input and output files/directories.")
    inout_params.add_argument('--in-csv', type=str, help='Path to the input transcript file containing feature names (e.g., gene names).')
    inout_params.add_argument('--out-record', type=str, default=None, help='Path to an output tsv file with feature name and filtering status. Only generated if provided.')
    inout_params.add_argument('--out-csv', type=str, default=None, help='Path to an compressed output file with remaining features aftering filtering. Only generated if provided.')
    inout_params.add_argument('--csv-delim', type=str, default='\t', help='Delimiter used in the input transcript file (default: tab).')
    inout_params.add_argument('--csv-colname-feature-name', type=str, default="gene", help='Column name containing gene names in the input transcript file (default: gene).')
    inout_params.add_argument('--log', action='store_true', default=False, help='Enable logging.')

    ftr_params = parser.add_argument_group("Feature Filtering Parameters", "Filtering criteria for features. Must provide at least one filtering criteria.")
    ftr_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    ftr_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    ftr_params.add_argument('--include-feature-substr', type=str, default=None, help='A substring of feature/gene names to be included (default: None)')
    ftr_params.add_argument('--exclude-feature-substr', type=str, default=None, help='A substring of feature/gene names to be excluded (default: None)')
    ftr_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    ftr_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')
    ftr_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). When --include-feature-type-regex, provide --feature-type-ref with either colidx or colname for name and type') # (e.g. protein_coding|lncRNA) # (e.g. protein_coding|lncRNA)
    ftr_params.add_argument('--feature-type-ref', type=str, default=None, help='If --include-feature-type-regex, specify the path to a tab-separated reference file to provide gene type information for each each per row. (default: None)')
    ftr_params.add_argument('--feature-type-ref-delim', type=str, default='\t', help='Delimiter used in the reference file (default: tab).')
    ftr_params.add_argument('--feature-type-ref-colname-name', type=str, default=None, help='Column name for gene name in the reference file (default: None).')
    ftr_params.add_argument('--feature-type-ref-colname-type', type=str, default=None, help='Column name for gene type in the reference file (default: None).')
    ftr_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='Column index for gene name in the reference file (default: None).')
    ftr_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='Column index for gene type in the reference file (default: None).')

    args = parser.parse_args(_args)

    # sanity check 
    assert os.path.exists(args.in_csv), f"Input file does not exist: {args.in_csv}"
    assert args.out_record is not None or args.out_csv is not None, "Provide either --out-record or --out-csv to save the filtered features."
    assert any([args.include_feature_list, args.exclude_feature_list, args.include_feature_substr, args.exclude_feature_substr, args.include_feature_regex, args.exclude_feature_regex, args.include_feature_type_regex]), "At least one filtering criteria is required."

    # output
    if args.out_record:
        out_dir = os.path.dirname(args.out_record)
        out_log = args.out_record.replace(".tsv", ".log")
    elif args.out_csv:
        out_dir = os.path.dirname(args.out_csv)
        out_log = args.out_csv.replace(".tsv.gz", ".log")

    os.makedirs(out_dir, exist_ok=True)

    logger = create_custom_logger(__name__, out_log if args.log else None)

    # ===============================
    #  Annotate Feature filtering information
    # ===============================

    # Read all input feature names
    logger.info(" Reading input --in-csv: "+args.in_csv)
    df_ftrinfo = pd.read_csv(args.in_csv, usecols=[args.csv_colname_feature_name], sep=args.csv_delim, header=0, comment="#")
    df_ftrinfo = map_filter_per_feature(df_ftrinfo,
                                        include_feature_list=args.include_feature_list, exclude_feature_list=args.exclude_feature_list,
                                        include_feature_substr=args.include_feature_substr, exclude_feature_substr=args.exclude_feature_substr,
                                        include_feature_regex=args.include_feature_regex, exclude_feature_regex=args.exclude_feature_regex,
                                        include_feature_type_regex=args.include_feature_type_regex, feature_type_ref=args.feature_type_ref, feature_type_ref_delim=args.feature_type_ref_delim, 
                                        feature_type_ref_colname_name=args.feature_type_ref_colname_name, feature_type_ref_colname_type=args.feature_type_ref_colname_type, feature_type_ref_colidx_name=args.feature_type_ref_colidx_name, feature_type_ref_colidx_type=args.feature_type_ref_colidx_type,
                                        logger=logger)

    # Save the filtered features
    if args.out_record is not None:
        df_ftrinfo.to_csv(args.out_record, sep=args.csv_delim, index=False)
        logger.info(f" - Saved feature filtering records to: {args.out_record}")

    # ===============================
    #  filtering
    # ===============================
    if args.out_csv is not None:
        # read-in
        logger.info(" Applying filtering...")
        df_ftrinfo = df_ftrinfo[df_ftrinfo["filtering"].isna()]
        n_raw_csv = 0
        n_filtered_csv = 0
        for chunk in pd.read_csv(args.in_csv, chunksize=50000, index_col=None, sep=args.csv_delim, comment="#"):
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
        logger.info(f"  - Saved filtered csv to: {args.out_csv}")
    
if __name__ == "__main__":#
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
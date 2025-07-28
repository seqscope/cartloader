import argparse, os, re, sys, logging, gzip, inspect
import pandas as pd
import numpy as np

from cartloader.utils.utils import create_custom_logger, log_info

def filter_feature_by_type(include_feature_type_regex, feature_type_ref, feature_type_ref_delim="\t", 
                           feature_type_ref_colname_name=None, feature_type_ref_colname_type=None, 
                           feature_type_ref_colidx_name=None, feature_type_ref_colidx_type=None, 
                           logger=None):    
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
    for chunk in pd.read_csv(feature_type_ref, chunksize=50000, index_col=None, sep=feature_type_ref_delim, comment="#"):
        if feature_type_ref_colname_name is not None and feature_type_ref_colname_type is not None:
            chunk = chunk[[feature_type_ref_colname_name, feature_type_ref_colname_type]]
        elif feature_type_ref_colidx_name is not None and feature_type_ref_colidx_type is not None:
            chunk = chunk.iloc[:, [feature_type_ref_colidx_name, feature_type_ref_colidx_type]]
        chunk.columns = ["ftr", "ftr_type"]
        n_ftr += chunk.shape[0]
        ftrlist_keep_type.extend(chunk[chunk["ftr_type"].str.contains(include_feature_type_regex, regex=True)]["ftr"].tolist())
    ftrlist_keep_type = list(set(ftrlist_keep_type)) 
    log_info(f"  - N genes in the reference file: {n_ftr}", logger)
    log_info(f"  - N genes with the specified feature type: {len(ftrlist_keep_type)}", logger)
    return ftrlist_keep_type

def feature_filtering_by_type(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                        description="""Filter features by gene type using a reference file.
                                        The reference file should contain gene name and gene type information.
                                        """)
    ftr_params = parser.add_argument_group("Feature Filtering Parameters", "Filtering criteria for features.")
    ftr_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). When --include-feature-type-regex, provide --feature-type-ref with either colidx or colname for name and type') # (e.g. protein_coding|lncRNA) # (e.g. protein_coding|lncRNA)
    ftr_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated reference file to provide gene type information for each each per row. (default: None)')
    ftr_params.add_argument('--feature-type-ref-delim', type=str, default='\t', help='Delimiter used in the reference file (default: tab).')
    ftr_params.add_argument('--feature-type-ref-colname-name', type=str, default=None, help='Column name for gene name in the reference file (default: None).')
    ftr_params.add_argument('--feature-type-ref-colname-type', type=str, default=None, help='Column name for gene type in the reference file (default: None).')
    ftr_params.add_argument('--feature-type-ref-colidx-name', type=int, default=None, help='Column index for gene name in the reference file (default: None).')
    ftr_params.add_argument('--feature-type-ref-colidx-type', type=int, default=None, help='Column index for gene type in the reference file (default: None).')
    ftr_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    ftr_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    ftr_params.add_argument('--out-feature', type=str, default=None, help='Output file path to save the filtered features (default: None)')
    ftr_params.add_argument('--log', action='store_true', default=False, help='Enable logging.')
    args = parser.parse_args(_args)

    assert args.out_feature is not None, "Provide --out-feature to save the filtered feature list."
    out_dir=os.path.dirname(args.out_feature)
    os.makedirs(out_dir, exist_ok=True)
    
    logger = create_custom_logger(__name__, os.path.join(out_dir, args.out_feature + ".log") if args.log else None)

    # log the sys.argv
    log_info(f"Running: {' '.join(sys.argv)}", logger)
    ftrlist_keep_type= filter_feature_by_type(args.include_feature_type_regex, 
                                              feature_type_ref = args.feature_type_ref, feature_type_ref_delim=args.feature_type_ref_delim, 
                                              feature_type_ref_colidx_name=args.feature_type_ref_colidx_name, feature_type_ref_colidx_type=args.feature_type_ref_colidx_type, 
                                              feature_type_ref_colname_name=args.feature_type_ref_colname_name, feature_type_ref_colname_type=args.feature_type_ref_colname_type,
                                              logger=logger)
    
    # if --include-feature-list is provided, merge the gene list with the gene type list
    if args.include_feature_list is not None:
        df_ftrlist_raw = pd.read_csv(args.include_feature_list, header=None, names=["gene_name"])
        ftrlist_include = df_ftrlist_raw["gene_name"].tolist()
        log_info(f"  - N genes in --include-feature-list : {len(ftrlist_include)}", logger)
        ftrlist_include = list(set(ftrlist_include) & set(ftrlist_keep_type))
        log_info(f"  - N genes after applying --include-feature-list : {len(ftrlist_include)}", logger)
    if args.exclude_feature_list is not None:
        df_ftrlist_raw = pd.read_csv(args.exclude_feature_list, header=None, names=["gene_name"])
        ftrlist_exclude = df_ftrlist_raw["gene_name"].tolist()
        log_info(f"  - N genes in --exclude-feature-list : {len(ftrlist_exclude)}", logger)
        ftrlist_keep_type = list(set(ftrlist_keep_type) - set(ftrlist_exclude))
        log_info(f"  - N genes after applying --exclude-feature-list : {len(ftrlist_keep_type)}", logger)
    
    # convert a pd.dataframe to a file
    df_ftr = pd.DataFrame(ftrlist_keep_type, columns=["gene_name"])
    df_ftr.to_csv(args.out_feature, index=False, header=False, sep="\t")
    log_info(f"  - Saved the filtered feature list to {args.out_feature}", logger)
    
if __name__ == "__main__":#
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
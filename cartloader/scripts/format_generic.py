import argparse, os, re, sys, logging, gzip, inspect
import pandas as pd
import numpy as np
import subprocess
import hashlib

from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, log_info
from cartloader.scripts.feature_filtering import map_filter_per_feature

def format_generic(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Standardize Spatial Transcriptomics (ST) datasets from various platforms, including 10X Xenium, BGI Stereoseq, Cosmx SMI, Vizgen Merscope, Pixel-Seq. 
                                     It will generate output in micrometer precision, including a transcript-indexed SGE file in TSV format, a feature file couning UMIs per gene, and a minmax file for X Y coordinates.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--input', type=str, help='Specify the input transcript-indexed SGE file in TSV or CSV format')
    inout_params.add_argument('--out-dir', required= True, type=str, help='The output directory')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv", help='The output transcript-indexed SGE file in TSV format. Default: transcripts.unsorted.tsv')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='The output TSV file for min and max X and Y coordinates. Default: coordinate_minmax.tsv')
    inout_params.add_argument('--out-feature', type=str, default="features.clean.tsv.gz", help='The output file collects UMI counts on a per-gene basis. Default: features.clean.tsv.gz')
    
    key_params = parser.add_argument_group("Key Parameters", "Key parameters, such as filtering cutoff.")
    key_params.add_argument('--precision-um', type=int, default=2, help='Number of digits to store the transcript coordinates in micrometer')
    key_params.add_argument('--units-per-um', type=float, default=1, help='Units per micrometer (default: 1)')
    key_params.add_argument('--log', type=str, default=None, help='Specify the log file path (default: None)')

    incol_params = parser.add_argument_group("Input Columns Parameters", "Input column parameters.")
    incol_params.add_argument('--csv-comment',  action='store_true', default=False, help='Lines starts with # is comment lines in the input file (default: False).')
    incol_params.add_argument('--csv-delim', type=str, default='\t', help='Delimiter for the additional barcode position file (default: \\t).')
    incol_params.add_argument('--csv-colname-x',  type=str, default=None, required=True, help='Specify the input column name for X-axis')
    incol_params.add_argument('--csv-colname-y',  type=str, default=None, required=True, help='Specify the input column name for Y-axis')
    incol_params.add_argument('--csv-colname-feature-name', type=str, default=None, required=True, help='Specify the input column name for gene name')
    incol_params.add_argument('--csv-colnames-count', type=str, default=None, help='Specify column name(s) for UMI count. For multiple columns, separate names with commas. If not provided, each feature will be assigned a default count of 1 per pixel.')
    incol_params.add_argument('--csv-colnames-others', nargs='+', default=[], help='Specify the input column names to keep (e.g. cell_id, overlaps_nucleus). Please note this will affect how the count is aggregated (default: [])')

    outcol_params = parser.add_argument_group("Output Columns Parameters", "Output column parameters .")
    outcol_params.add_argument('--colname-x', type=str, default='X', help='Specify the output column name for X (default: X)')
    outcol_params.add_argument('--colname-y', type=str, default='Y', help='Specify the output column name for Y (default: Y)')
    outcol_params.add_argument('--colname-feature-name', type=str, default='gene', help='Specify the output column name for feature(gene) name')
    outcol_params.add_argument('--colnames-count', type=str, default='count', help='Output column name(s) for UMI count. If multiple input columns are specified in `--csv-colnames-count`, provide corresponding output column names (default: count).')

    aux_params = parser.add_argument_group("Auxiliary Parameters", 
                                           """
                                           Auxiliary parameters. 
                                           1) Use --csv-colname-feature-id and --colname-feature-id to add a gene ID column to the output transcript-indexed SGE tsv file.
                                           2) Use --add-molecule-id to add a molecule ID to the output transcript-indexed SGE tsv file.
                                           3) Use --csv-colname-phredscore with --min-phred-score to filter out low-quality reads. 
                                           4) Use --*-feature-list or --*-feature-substr or --*-feature-regex to filter out features by feature name.
                                           5) Use --include-feature-type-regex with --csv-colname-feature-type or --feature-type-ref to filter out features by gene type.
                                           """)
    aux_params.add_argument('--csv-colname-feature-id', type=str, default=None, help='Specify the input column name for gene ID if available (default: None)')
    aux_params.add_argument('--colname-feature-id', type=str, default='gene_id', help='Specify the output column name for gene ID when --csv-colname-feature-id is enabled (default: gene_id)')
    aux_params.add_argument('--add-molecule-id', action='store_true', default=False, help='If needed, use --add-molecule-id to generate a column of "molecule_id" in the output file to track the index of the original input. If the input file already has a column of molecule ID, use --csv-colnames-others instead of --add-molecule-id. (default: False)')
    aux_params.add_argument('--csv-colname-phredscore', type=str, default=None, help='Specify the input column name for Phred-scaled quality score (Q-Score) estimating the probability of incorrect call. If provided,  (platform: 10x_xenium; default: None)')
    aux_params.add_argument('--min-phred-score', type=float, default=None, help='Specify the Phred-scaled quality score cutoff') 
    aux_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    aux_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    aux_params.add_argument('--include-feature-substr', type=str, default=None, help='A substring of feature/gene names to be included (default: None)')
    aux_params.add_argument('--exclude-feature-substr', type=str, default=None, help='A substring of feature/gene names to be excluded (default: None)')
    aux_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')
    aux_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). It must be used with --csv-colname-feature-type or --feature-type-ref') # (e.g. protein_coding|lncRNA)
    aux_params.add_argument('--csv-colname-feature-type', type=str, default=None, help='If your input file has a column of gene type information, specify the column name (default: None)')
    aux_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated reference file to provide gene type information for each each per row (default: None)')
    aux_params.add_argument('--feature-type-ref-delim', type=str, default='\t', help='If --feature-type-ref, define delimiter used in the reference file (default: tab).')
    aux_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='If --feature-type-ref, define the column index for gene name in the reference file (default: None).')
    aux_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='If --feature-type-ref, define the column index for gene type in the reference file (default: None).')
    aux_params.add_argument('--print-removed-transcripts', action='store_true', default=False, help='Print the list of removed transcript with corresponding filtering criteria (default: False)')
    args = parser.parse_args(_args)

    # output
    os.makedirs(args.out_dir, exist_ok=True)
    out_transcript_path=os.path.join(args.out_dir, args.out_transcript)
    out_feature_path=os.path.join(args.out_dir, args.out_feature)
    out_minmax_path=os.path.join(args.out_dir, args.out_minmax)    

    # initial output
    # 1) feature
    feature=pd.DataFrame()

    # 2) minmax
    xmin = np.inf
    ymin = np.inf
    xmax = -np.inf
    ymax = -np.inf
    # xmin=sys.maxsize
    # xmax=0
    # ymin=sys.maxsize
    # ymax=0

    # 3) transcript (header)
    csv_comment="#" if args.csv_comment else None
    print(f"csv_comment: {csv_comment}")
    
    feature_cols=[args.colname_feature_name] if args.csv_colname_feature_id is None else [args.colname_feature_name, args.colname_feature_id]
    other_cols = args.csv_colnames_others + ["molecule_id"] if args.add_molecule_id else args.csv_colnames_others 
    
    iheader = pd.read_csv(args.input, nrows=1, sep=args.csv_delim, comment=csv_comment).columns.tolist()
    icols = [args.csv_colname_x, args.csv_colname_y, args.csv_colname_feature_name] + args.csv_colnames_others
    for j in [args.csv_colname_feature_id, args.csv_colname_feature_type, args.csv_colname_phredscore]:
        if j is not None:
            icols.append(j)
    
    countcols_out= args.colnames_count.split(",")
    if args.csv_colnames_count is not None:
        countcols_in = args.csv_colnames_count.split(",")
        icols += countcols_in

        # check if the number of input and output column names for UMI count are the same
        assert len(countcols_in) == len(countcols_out), "The number of input and output column names for UMI count must be the same."
        countcols_in2out={countcols_in[i]:countcols_out[i] for i in range(len(countcols_in))}
    else:
        countcols_in = []
        countcols_in2out = {}

    for icol in icols:
        assert icol in iheader, f"Column {icol} is not found in the input file."
    
    unit_info = [args.colname_x, args.colname_y] + feature_cols + other_cols
    oheader = unit_info + countcols_out
    with open(out_transcript_path, 'w') as wf:
        _ = wf.write('\t'.join(oheader)+'\n')
    
    # params
    # 1) float_format
    float_format="%.2f"
    if args.precision_um >= 0:
        float_format = f"%.{args.precision_um}f"
    
    # 2) renames: column names from old to new
    col_dict={args.csv_colname_x:args.colname_x,
              args.csv_colname_y:args.colname_y,
              args.csv_colname_feature_name:args.colname_feature_name}
    
    if args.csv_colname_feature_id is not None:
        col_dict[args.csv_colname_feature_id] = args.colname_feature_id

    # 3) feature preprocessing, read all features from the input and apply all ftr-related filters. 
    #   This returns a df with columns [feature, filtering]
    #   For a keep feature, the filtering column is na
    if any([args.include_feature_list, args.exclude_feature_list, args.include_feature_substr, args.exclude_feature_substr, args.include_feature_regex, args.exclude_feature_regex, args.include_feature_type_regex]):
        # ftrinfo_id = hashlib.md5(";".join([
        #     str(i) if i is not None else ""  
        #     for i in [args.include_feature_list, args.exclude_feature_list, args.include_feature_substr, args.exclude_feature_substr, args.include_feature_regex, args.exclude_feature_regex, args.include_feature_type_regex, args.csv_colname_feature_type, args.feature_type_ref]]).encode()).hexdigest()[:10]
        df_ftrinfo = pd.read_csv(args.input, usecols=[args.csv_colname_feature_name], sep=args.csv_delim, index_col=None, header=0, comment=csv_comment)
        feature_filter_args = {
            "include_feature_list": args.include_feature_list,
            "exclude_feature_list": args.exclude_feature_list,
            "include_feature_substr": args.include_feature_substr,
            "exclude_feature_substr": args.exclude_feature_substr,
            "include_feature_regex": args.include_feature_regex,
            "exclude_feature_regex": args.exclude_feature_regex,
            "include_feature_type_regex": args.include_feature_type_regex,
            "logger": None,
        }
        if args.include_feature_type_regex is not None:
            if args.csv_colname_feature_type is not None:
                feature_filter_args.update({
                    "feature_type_ref": args.input,
                    "feature_type_ref_delim": args.csv_delim,
                    "feature_type_ref_colname_name": args.csv_colname_feature_name,
                    "feature_type_ref_colname_type": args.csv_colname_feature_type,
                })
            elif args.feature_type_ref is not None:
                feature_filter_args.update({
                    "feature_type_ref": args.feature_type_ref,
                    "feature_type_ref_delim": args.feature_type_ref_delim,
                    "feature_type_ref_colidx_name": args.feature_type_ref_colidx_name,
                    "feature_type_ref_colidx_type": args.feature_type_ref_colidx_type,
                })
        df_ftrinfo = map_filter_per_feature(df_ftrinfo, **feature_filter_args)
        ftr2filter = df_ftrinfo.set_index("feature")["filtering"].to_dict()
    else:
        ftr2filter = {}

    # 4) phred score filtering
    if args.csv_colname_phredscore is not None and args.min_phred_score is None:
        print(f"Warning: While the --csv-colname-phredscore is enabled, the --min-phred-score is not provided. carloader will SKIP filtering SGE by phred score.")

    # 5) ini a list for collecting rm rows
    filtered_out_rows = []

    # processing
    for chunk in pd.read_csv(args.input, header=0, chunksize=500000, index_col=None, sep=args.csv_delim, comment=csv_comment):
        # drop the lines starts with '#' (novast start with '#')
        #chunk = chunk[~chunk[args.csv_colname_x].str.startswith('#')]

        Rraw=chunk.shape[0]
        if len(countcols_in) > 0:
            Craw="s;".join([col+":"+str(chunk[col].sum()) for col in countcols_in])
        else:
            Craw=Rraw
        
        # filter by feature 
        # if ftr2filter is not empty
        if ftr2filter:
            # add ftr2filter value to removed by key == csv_colname_feature_name
            removed = chunk.copy()
            removed["reason"] = removed[args.csv_colname_feature_name].map(ftr2filter)
            # drop the rows that the reason is na
            removed = removed[removed["reason"].notna()]
            filtered_out_rows.append(removed)
            chunk = chunk[chunk[args.csv_colname_feature_name].map(ftr2filter).isna()]

        # filter by phred scores (low-quality reads)
        if args.csv_colname_phredscore is not None and args.min_phred_score is not None:
            removed = chunk[chunk[args.csv_colname_phredscore] < args.min_phred_score].copy()
            removed['reason'] = 'min_phred_score'
            filtered_out_rows.append(removed)
            chunk = chunk[chunk[args.csv_colname_phredscore] >= args.min_phred_score]
        
        # rename columns
        chunk.rename(columns = col_dict, inplace=True)
        
        # coord conversion (applied float_format before aggregate the count)
        if args.units_per_um != 1:
            chunk[args.colname_x] = chunk[args.colname_x] / args.units_per_um
            chunk[args.colname_y] = chunk[args.colname_y] / args.units_per_um
        chunk[args.colname_x] = chunk[args.colname_x].map(lambda x: float(float_format % x))
        chunk[args.colname_y] = chunk[args.colname_y].map(lambda y: float(float_format % y))
        
        # add molecule id if provided
        if args.add_molecule_id:
            chunk["molecule_id"] = chunk.index.values
        
        # umi count aggregation
        # * Applying args.csv_colnames_others may introduce NaN, and NaN in any column in unit_info may reduce the number of groups.
        #   So, replace NaN with a placeholder before grouping to ensure they are not causing issues in the grouping process.
        # 1) add a count column
        if len(countcols_in) > 0:
            for k,v in countcols_in2out.items():
                chunk.rename(columns = {k:v}, inplace=True)
        else:
            chunk[args.colnames_count] = 1
        # 2) drop unnecessary columns
        chunk=chunk[unit_info + countcols_out]
        # 3) replace NaN with a placeholder
        for col in chunk.columns:
            #chunk[col].fillna('NA', inplace=True) # avoiding chained assignments due to FutureWarning
            chunk[col] = chunk[col].fillna('NA')
        # 4) group by unit_info and aggregate each count
        # for all cols in countcols_out, aggregate the count
        # chunk = chunk.groupby(by = unit_info).agg({args.colnames_count:'sum'}).reset_index()
        chunk = chunk.groupby(by=unit_info).agg({col: 'sum' for col in countcols_out}).reset_index()

        # 5) replace back
        for col in chunk.columns:
            #chunk[col].replace('NA', np.nan, inplace=True) # avoiding chained assignments due to FutureWarning
            chunk[col] = chunk[col].replace('NA', np.nan)
        
        # write down
        chunk[oheader].to_csv(out_transcript_path, sep='\t',mode='a',index=False,header=False, float_format=float_format, na_rep="NA")
        
        # log
        Rqc=chunk.shape[0]
        #Cqc=chunk[args.colnames_count].sum()
        Cqc = " ;".join([col+":"+str(chunk[col].sum()) for col in countcols_out])

        print(f"processing: 1) rows {Rraw} -> {Rqc}\t2) counts {Craw} -> {Cqc}")
        
        # feature
        #feature = pd.concat([feature, chunk.groupby(by=feature_cols).agg({args.colnames_count:"sum"}).reset_index()])
        feature = pd.concat([feature, chunk.groupby(by=feature_cols).agg({col: 'sum' for col in countcols_out}).reset_index()])
        
        # minmax chunk
        x0 = chunk[args.colname_x].min()
        x1 = chunk[args.colname_x].max()
        y0 = chunk[args.colname_y].min()
        y1 = chunk[args.colname_y].max()
        xmin = min(xmin, x0)
        xmax = max(xmax, x1)
        ymin = min(ymin, y0)
        ymax = max(ymax, y1)

    # print out the removed features

    if args.print_removed_transcripts:
        if len(filtered_out_rows) > 0:
            filtered_out_df = pd.concat(filtered_out_rows, ignore_index=True)
            filtered_out_df.to_csv(os.path.join(args.out_dir, "transcripts.removed.tsv"), sep='\t', index=False)
        else:
            # touch an empty file
            with open(os.path.join(args.out_dir, "transcripts.removed.tsv"), 'w') as wf:
                pass

    # feature
    #feature = feature.groupby(by=feature_cols).agg({args.colnames_count:"sum"}).reset_index()
    feature = feature.groupby(by=feature_cols).agg({col: 'sum' for col in countcols_out}).reset_index()
    feature.to_csv(out_feature_path, sep='\t',index=False)

    # minmax
    with open(out_minmax_path, 'w') as wf:
        wf.write(f"xmin\t{xmin:.2f}\n")
        wf.write(f"xmax\t{xmax:.2f}\n")
        wf.write(f"ymin\t{ymin:.2f}\n")
        wf.write(f"ymax\t{ymax:.2f}\n")

if __name__ == "__main__":#
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
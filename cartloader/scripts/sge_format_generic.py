import argparse, os, re, sys, logging, gzip, inspect
import pandas as pd
import numpy as np
import subprocess
import hashlib

from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, log_info
from cartloader.scripts.feature_filtering import map_filter_per_feature

def sge_format_generic(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Standardize Spatial Transcriptomics (ST) datasets in CSV/TSV format.
                                     Output: a transcript-indexed SGE file in TSV format, a feature file couning UMIs per gene, and a minmax file for X Y coordinates.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--input', type=str, help='Path to input transcript-indexed SGE TSV/CSV')
    inout_params.add_argument('--out-dir', required= True, type=str, help='Path to output directory')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv", help='File name of output transcript-indexed SGE TSV (default: transcripts.unsorted.tsv)')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='File name of output coordinate min/max TSV (default: coordinate_minmax.tsv)')
    inout_params.add_argument('--out-feature', type=str, default="features.clean.tsv.gz", help='File name of output per-feature UMI count TSV (default: features.clean.tsv.gz)')
    
    key_params = parser.add_argument_group("Key Parameters")
    key_params.add_argument('--precision-um', type=int, default=2, help='Number of digits for transcript coordinates in micrometers (default: 2)')
    key_params.add_argument('--units-per-um', type=float, default=1, help='Units per micrometer (default: 1)')
    key_params.add_argument('--min-phred-score', type=float, default=None, help='Minimum Phred-scaled quality score cutoff (used with --csv-colname-phredscore)') 
    key_params.add_argument('--print-removed-transcripts', action='store_true', default=False, help='Print list of removed transcripts with reasons')
    key_params.add_argument('--add-molecule-id', action='store_true', default=False, help='Enable to create a "molecule_id" column in the output to track the index of the original input. If the input file already has a column of molecule ID, use --csv-colnames-others instead of --add-molecule-id.')

    incol_params = parser.add_argument_group("Input Columns Parameters")
    incol_params.add_argument('--csv-comment',  action='store_true', default=False, help='Treat lines starting with # as comments in --input')
    incol_params.add_argument('--csv-delim', type=str, default='\t', help='Delimiter in --input (default: "\t")')
    incol_params.add_argument('--csv-colname-x',  type=str, default=None, required=True, help='Column name of X coordinate in --input')
    incol_params.add_argument('--csv-colname-y',  type=str, default=None, required=True, help='Column name of Y coordinate in --input')
    incol_params.add_argument('--csv-colname-feature-name', type=str, default=None, required=True, help='Column name of gene name in --input')
    incol_params.add_argument('--csv-colname-count', type=str, default=None, help='Column name of UMI count in --input. If not provided, each feature is assigned a count of 1 per pixel')
    incol_params.add_argument('--csv-colnames-others', nargs='+', default=[], help='Column names to keep from --input (e.g., cell_id, overlaps_nucleus). Affects aggregation')
    incol_params.add_argument('--csv-colname-feature-id', type=str, default=None, help='Column name of gene ID in --input, if available')
    incol_params.add_argument('--csv-colname-phredscore', type=str, default=None, help='Column name of Phred-scaled quality score (Q-Score) in --input. Use with --min-phred-score (platform: 10x_xenium; default: None)')

    outcol_params = parser.add_argument_group("Output Columns Parameters")
    outcol_params.add_argument('--colname-x', type=str, default='X', help='Column name for X in the output (default: X)')
    outcol_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y in the output (default: Y)')
    outcol_params.add_argument('--colname-count', type=str, default='count', help='Column name for UMI count in the output (default: count)')
    outcol_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for gene name in the output (default: gene)')
    outcol_params.add_argument('--colname-feature-id', type=str, default='gene_id', help='Column name for gene ID in the output (used with --csv-colname-feature-id; default: gene_id)')

    aux_params = parser.add_argument_group("Auxiliary Feature Filtering Parameters")
    aux_params.add_argument('--include-feature-list', type=str, default=None, help='Path to file listing genes to include')
    aux_params.add_argument('--exclude-feature-list', type=str, default=None, help='Path to file listing genes to exclude')
    aux_params.add_argument('--include-feature-regex', type=str, default=None, help='Regex of feature names to include')
    aux_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude')
    aux_params.add_argument('--include-feature-type-regex', type=str, default=None, help='Regex of feature type to include (default: None). Use with --csv-colname-feature-type or --feature-type-ref')
    aux_params.add_argument('--csv-colname-feature-type', type=str, default=None, help='Column name of feature type in --input')
    aux_params.add_argument('--feature-type-ref', type=str, default=None, help='Path to tab-separated reference with gene type per row')
    aux_params.add_argument('--feature-type-ref-delim', type=str, default='\t', help='Delimiter in --feature-type-ref (default: "\t")')
    aux_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='Column index for feature name in --feature-type-ref')
    aux_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='Column index for feature type in --feature-type-ref')
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
    # update: the first four columns should consistently be: x, y, feature_name, countcol

    # * input header and input columns
    csv_comment="#" if args.csv_comment else None
    print(f"csv_comment: {csv_comment}")
    
    icols_count = [args.csv_colname_count] if args.csv_colname_count is not None else []
    icols = [args.csv_colname_x, args.csv_colname_y, args.csv_colname_feature_name, *args.csv_colnames_others, *icols_count] + [c for c in [args.csv_colname_feature_id, args.csv_colname_feature_type, args.csv_colname_phredscore] if c]
    
    iheader = pd.read_csv(args.input, nrows=1, sep=args.csv_delim, comment=csv_comment).columns.tolist()
    
    # * output header and output columns
    ocols_count = [args.colname_count]
    ocols_ftrs = [args.colname_feature_name] if args.csv_colname_feature_id is None else [args.colname_feature_name, args.colname_feature_id]
    ocols_others = args.csv_colnames_others + ([args.csv_colname_feature_id] if args.csv_colname_feature_id is not None else []) + (["molecule_id"] if args.add_molecule_id else [])
    
    unit_info = [args.colname_x, args.colname_y, args.colname_feature_name] + ocols_others

    oheader = [args.colname_x, args.colname_y, args.colname_feature_name] + ocols_count + ocols_others

    with open(out_transcript_path, 'w') as wf:
        _ = wf.write('\t'.join(oheader)+'\n')
    
    # * input output dict:
    countcols_in2out={icols_count[i]:ocols_count[i] for i in range(len(icols_count))} if args.csv_colname_count is not None else {}

    col_dict = {
        args.csv_colname_x: args.colname_x,
        args.csv_colname_y: args.colname_y,
        args.csv_colname_feature_name: args.colname_feature_name,
        **({args.csv_colname_feature_id: args.colname_feature_id} if args.csv_colname_feature_id else {})
    }

    # * sanity check for input and output columns
    #   - all specified input columns should be available in the input tsv
    for icol in icols:
        assert icol in iheader, f"Missing input column(s) ({icol}) in the header of the input file: {iheader}."
    #   - len of column should be consistent between input and output
    if args.csv_colname_count is not None:
        assert len(icols_count) == len(ocols_count), f"Mismatch in the number of UMI count columns in input(N={len(icols_count)}) and output (N={len(ocols_count)})"

    # * float_format
    float_format="%.2f"
    if args.precision_um >= 0:
        float_format = f"%.{args.precision_um}f"

    # 4) feature preprocessing, read all features from the input and apply all ftr-related filters. 
    #   This returns a df with columns [feature, filtering]. For a keep feature, the filtering column is na
    if any([args.include_feature_list, args.exclude_feature_list, args.include_feature_regex, args.exclude_feature_regex, args.include_feature_type_regex]):
        df_ftrinfo = pd.read_csv(args.input, usecols=[args.csv_colname_feature_name], sep=args.csv_delim, index_col=None, header=0, comment=csv_comment)
        feature_filter_args = {
            "include_feature_list": args.include_feature_list,
            "exclude_feature_list": args.exclude_feature_list,
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

    # 5) phred score filtering
    if args.csv_colname_phredscore is not None and args.min_phred_score is None:
        print(f"Warning: While the --csv-colname-phredscore is enabled, the --min-phred-score is not provided. carloader will SKIP filtering SGE by phred score.")

    # 6) ini a list for collecting rm rows
    filtered_out_rows = []

    # processing
    for chunk in pd.read_csv(args.input, header=0, chunksize=500000, index_col=None, sep=args.csv_delim, comment=csv_comment):

        Rraw=chunk.shape[0]
        if len(icols_count) > 0:
            Craw=";".join([col+":"+str(chunk[col].sum()) for col in icols_count])
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
        if len(icols_count) > 0:
            for k,v in countcols_in2out.items():
                chunk.rename(columns = {k:v}, inplace=True)
            # filter by counts: if the sum of icols_count is 0, drop the row
            chunk = chunk[chunk[ocols_count].sum(axis=1) > 0]
        else:
            chunk[args.colname_count] = 1
        # 2) drop unnecessary columns
        chunk=chunk[unit_info + ocols_count]
        # 3) replace NaN with a placeholder
        for col in chunk.columns:
            chunk[col] = chunk[col].fillna('NA')
        # 4) Group by unit_info and Aggregate each count. Allows >=1 count column(s).
        chunk = chunk.groupby(by=unit_info).agg({col: 'sum' for col in ocols_count}).reset_index()
        # 5) replace back
        for col in chunk.columns:
            chunk[col] = chunk[col].replace('NA', np.nan)
        
        # write down
        chunk[oheader].to_csv(out_transcript_path, sep='\t',mode='a',index=False,header=False, float_format=float_format, na_rep="NA")
        
        # log
        Rqc=chunk.shape[0]
        Cqc = " ;".join([col+":"+str(chunk[col].sum()) for col in ocols_count])
        print(f"processing: 1) rows {Rraw} -> {Rqc}\t2) counts {Craw} -> {Cqc}")
        
        # feature
        feature = pd.concat([feature, chunk.groupby(by=ocols_ftrs).agg({col: 'sum' for col in ocols_count}).reset_index()])

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
    # TBC: reorder the feature columns, make sure that first column is feature name and the second column should be count
    feature = feature.groupby(by=ocols_ftrs).agg({col: 'sum' for col in ocols_count}).reset_index()
    oheader_ftr = ocols_ftrs + ocols_count
    feature = feature[ oheader_ftr ]
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

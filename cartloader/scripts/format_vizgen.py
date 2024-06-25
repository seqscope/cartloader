import sys, os, re, copy, gzip, time, logging, pickle, argparse, inspect
import numpy as np
import pandas as pd

def format_vizgen(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Format transcript file from Vizgen MERSCOPE format.")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--input', type=str, help='Input transcript file from Vizgen MERSCOPE, likely named like detected_transcripts.csv.gz')
    inout_params.add_argument('--out-dir', required= True, type=str, help='The output directory.')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv", help='The output transcript-indexed SGE file in TSV format. Default: transcripts.unsorted.tsv')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='The output coordinate minmax TSV file. Default: coordinate_minmax.tsv')
    inout_params.add_argument('--out-feature', type=str, default="features.clean.tsv.gz", help='The output files for gene. Default: features.clean.tsv.gz')

    incol_params = parser.add_argument_group("Input Columns Parameters", "Input column parameters .")
    incol_params.add_argument('--tsv-colname-x',  type=str, default='global_x', help='Column name for X-axis (default: global_x)')
    incol_params.add_argument('--tsv-colname-y',  type=str, default='global_y', help='Column name for Y-axis (default: global_y)')
    incol_params.add_argument('--tsv-colname-feature-name', type=str, default=None, help='Column name for gene name (e.g.: gene)')
    incol_params.add_argument('--tsv-colname-feature-id', type=str, default=None, help='Column name for gene id (e.g.: transcript_id)')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters, such as filtering cutoff.")
    key_params.add_argument('--dummy-genes', type=str, default='', help='A single name or a regex describing the names of negative control probes')
    key_params.add_argument('--precision-um', type=int, default=2, help='Number of digits to store the transcript coordinates in micrometer')
    key_params.add_argument('--debug', action="store_true", help='')
    
    outcol_params = parser.add_argument_group("Output Columns Parameters", "Output column parameters for TSV.")
    outcol_params.add_argument('--colname-x', type=str, default='X', help='Column name for X (default: X)')
    outcol_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y (default: Y)')
    outcol_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for feature/gene name (default: None)')
    outcol_params.add_argument('--colname-feature-id', type=str, default='gene_id', help='Column name for feature/gene name (default: None)')
    outcol_params.add_argument('--colnames-count', type=str, default='gn', help='Comma-separate column names for Count (default: gn)')
    outcol_params.add_argument('--colname-molecule-id', type=str, default='MoleculeID', help='Column name for MoleculeID (default: MoleculeID)')
 
    args = parser.parse_args(_args)

    logging.basicConfig(level=logging.INFO)

    # output
    os.makedirs(args.out_dir, exist_ok=True)
    out_transcript_path=os.path.join(args.out_dir, args.out_transcript)
    out_feature_path=os.path.join(args.out_dir, args.out_feature)
    out_minmax_path=os.path.join(args.out_dir, args.out_minmax)

    # params
    float_format = f"%.{args.precision_um}f"

    # feature
    feature=pd.DataFrame()
    feature_cols = [col for col in [args.tsv_colname_feature_name, args.tsv_colname_feature_id] if col is not None]
    if len(feature_cols) == 0:
        logging.error("Please provide at least one of the feature columns.")
        sys.exit(1)
    
    # minmax
    xmin = np.inf
    ymin = np.inf
    xmax = -np.inf
    ymax = -np.inf

    # transcript
    unit_info=[args.colname_x, args.colname_y] + feature_cols + [args.colname_molecule_id]
    oheader = unit_info + [args.colname_count]
    with open(out_transcript_path, 'w') as wf:
        _ = wf.write('\t'.join(oheader)+'\n')
    
    for chunk in pd.read_csv(args.input,header=0,chunksize=500000,index_col=0):
        # filter
        if args.dummy_genes != '':
            chunk = chunk[~chunk[args.tsv_colname_feature_name].str.contains(args.dummy_genes, flags=re.IGNORECASE, regex=True)]    
        # rename
        chunk.rename(columns = {args.tsv_colname_x:args.colname_x, args.tsv_colname_y:args.colname_y}, inplace=True)
        if args.tsv_colname_feature_name is not None:
            chunk.rename(columns = {args.tsv_colname_feature_name:args.colname_feature_name}, inplace=True)
        if args.tsv_colname_feature_id is not None:
            chunk.rename(columns = {args.tsv_colname_feature_id:args.colname_feature_id}, inplace=True)
        # count
        chunk[args.colname_count] = 1
        chunk[args.colname_molecule_id] = chunk.index.values
        x,y = chunk[[args.colname_x, args.colname_y]].values.min(axis = 0)
        xmin = min(xmin,x)
        ymin = min(ymin,y)
        x,y = chunk[[args.colname_x, args.colname_y]].values.max(axis = 0)
        xmax = max(xmax,x)
        ymax = max(ymax,y)
        chunk[oheader].to_csv(out_transcript_path,sep='\t',mode='a',index=False,header=False,float_format=float_format)
        logging.info(f"{chunk.shape[0]}")
        #feature = pd.concat([feature, chunk.groupby(by=['gene','transcript_id']).agg({'Count':sum}).reset_index()])
        feature = pd.concat([feature, chunk.groupby(by=feature_cols).agg({args.colname_count:"sum"}).reset_index()])
        if args.debug:
            break

    # feature
    #feature = feature.groupby(by=['gene','transcript_id']).agg({'Count':sum}).reset_index()
    #feature.to_csv(out_feature_path,sep='\t',index=False)
    feature = feature.groupby(by=feature_cols).agg({args.colname_count:"sum"}).reset_index()
    feature.to_csv(out_feature_path,sep='\t',index=False)

    # minmax
    line = f"xmin\t{xmin:.{args.precision_um}f}\nxmax\t{xmax:.{args.precision_um}f}\nymin\t{ymin:.{args.precision_um}f}\nymax\t{ymax:.{args.precision_um}f}\n"
    with open(out_minmax_path, 'w') as wf:
        _ = wf.write(line)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
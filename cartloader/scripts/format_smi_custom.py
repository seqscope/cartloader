import sys, os, re, copy, gzip, time, logging, pickle, argparse
import numpy as np
import pandas as pd

def format_smi():

    parser = argparse.ArgumentParser()
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--input', type=str, help='Input transcript file from Vizgen MERSCOPE, likely named like detected_transcripts.csv.gz')
    inout_params.add_argument('--out-dir', required= True, type=str, help='The output directory.')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv", help='The output transcript-indexed SGE file in TSV format. Default: transcripts.unsorted.tsv')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='The output coordinate minmax TSV file. Default: coordinate_minmax.tsv')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='The output files for gene. Default: feature.clean.tsv.gz')

    incol_params = parser.add_argument_group("Input Columns Parameters", "Input column parameters for CSV.")
    incol_params.add_argument('--csv-colname-x',  type=str, default='global_x', help='Column name for X-axis (default: global_x)')
    incol_params.add_argument('--csv-colname-y',  type=str, default='global_y', help='Column name for Y-axis (default: global_y)')
    incol_params.add_argument('--csv-colname-feature-name', type=str, default='target', help='Column name for gene name (e.g.: target)')
    incol_params.add_argument('--csv-colnames-others', nargs='+', default=[], help='Columns names to keep(e.g. cell_ID, CellComp).')

    outcol_params = parser.add_argument_group("Output Columns Parameters", "Output column parameters for CSV.")
    outcol_params.add_argument('--colname-x', type=str, default='X', help='Column name for X (default: X)')
    outcol_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y (default: Y)')
    outcol_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for feature/gene name (default: None)')
    # outcol_params.add_argument('--colname-feature-id', type=str, default='gene_id', help='Column name for feature/gene name (default: None)')
    outcol_params.add_argument('--colnames-count', type=str, default='gn', help='Comma-separate column names for Count (default: gn)')
 
    key_params = parser.add_argument_group("Key Parameters", "Key parameters, such as filtering cutoff, for TSV.")
    #key_params.add_argument('--px_to_um', type=float, help='Convert pixel unit as used in x_local_px and x_global_px to micrometer, this number should be found in the README of your SMI output, it is likely 0.12~0.18', default=1)
    key_params.add_argument('--units-per-um', type=float, default=1.00, help=' Coordinate unit per um (conversion factor) (default: 1.00)') 
    key_params.add_argument('--precision-um', type=int, default=-1, help='Number of digits to store the transcript coordinates (only if --px_to_um is in use). Set it to 0 to round to integer. Default is -1, without rounding.')
    key_params.add_argument('--dummy-genes', type=str, default='NegPrb', help='A single name or a regex describing the names of negative control probes')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    
    # output
    os.makedirs(args.out_dir, exist_ok=True)
    out_transcript_path=os.path.join(args.out_dir, args.out_transcript)
    out_feature_path=os.path.join(args.out_dir, args.out_feature)
    out_minmax_path=os.path.join(args.out_dir, args.out_minmax)

    # params
    float_format="%.2f"
    if args.precision_um >= 0:
        float_format = f"%.{args.precision_um}f"

    # feature
    feature=pd.DataFrame()

    # minmax
    xmin=sys.maxsize
    xmax=0
    ymin=sys.maxsize
    ymax=0

    # transcript header
    unit_info=[args.colname_x, args.colname_y, args.colname_feature_name] + args.csv_colnames_others
    oheader = unit_info + [args.colname_count]
    with open(out_transcript_path, 'w') as wf:
        _ = wf.write('\t'.join(oheader)+'\n')
    
    for chunk in pd.read_csv(args.input,header=0,chunksize=500000):
        # filter
        if args.dummy_genes != '':
            chunk = chunk[~chunk[args.csv_colname_feature_name].str.contains(args.dummy_genes, flags=re.IGNORECASE, regex=True)]
        # rename
        chunk.rename(columns = {args.csv_colname_x:args.colname_x, 
                                args.csv_colname_y:args.colname_y, 
                                args.csv_colname_feature_name:args.colname_feature_name}, inplace=True)
        # count
        chunk[args.colname_count] = 1
        chunk = chunk.groupby(by = unit_info).agg({args.colname_count:'sum'}).reset_index()
        # conversion
        #if args.px_to_um != 1:
            # chunk[args.csv_colname_x] *= args.px_to_um
            # chunk[args.csv_colname_y] *= args.px_to_um
        if args.units_per_um != 1:
            chunk[args.csv_colname_x] = chunk[args.csv_colname_x] / args.units_per_um
            chunk[args.csv_colname_y] = chunk[args.csv_colname_y] / args.units_per_um
        chunk[oheader].to_csv(out_transcript_path,sep='\t',mode='a',index=False,header=False,float_format=float_format)
        
        logging.info(f"{chunk.shape[0]}")
        feature = pd.concat([feature, chunk.groupby(by=args.colname_feature_name).agg({args.colname_count:"sum"}).reset_index()])
        x0 = chunk[args.colname_x].min()
        x1 = chunk[args.colname_x].max()
        y0 = chunk[args.colname_y].min()
        y1 = chunk[args.colname_y].max()
        xmin = min(xmin, x0)
        xmax = max(xmax, x1)
        ymin = min(ymin, y0)
        ymax = max(ymax, y1)

    # feature
    feature = feature.groupby(by=args.colname_feature_name).agg({args.colname_count:"sum"}).reset_index()
    feature.to_csv(out_feature_path,sep='\t',index=False)

    # minmax
    with open(out_minmax_path, 'w') as wf:
        wf.write(f"xmin\t{xmin:.2f}\n")
        wf.write(f"xmax\t{xmax:.2f}\n")
        wf.write(f"ymin\t{ymin:.2f}\n")
        wf.write(f"ymax\t{ymax:.2f}\n")

if __name__ == '__main__':
    format_smi()
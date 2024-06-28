import sys, os, re, copy, gzip, time, logging, pickle, argparse, inspect, warnings
import numpy as np
import pandas as pd

<<<<<<< HEAD
def parse_arguments(_args):
    parser = argparse.ArgumentParser()
    # input/output
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--input', type=str, help='Input transcript file from Xenium output, likely named transcripts.csv.gz')
    inout_params.add_argument('--out-dir', required= True, type=str, help='The output directory.')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.csv", help='The output transcript-indexed SGE file in csv format. Default: transcripts.unsorted.csv')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.csv", help='The output coordinate minmax csv file. Default: coordinate_minmax.csv')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.csv.gz", help='The output files for gene. Default: feature.clean.csv.gz')
    # input columns
=======
def format_xenium(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Format transcript file from 10X Xenium format.")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--input', type=str, help='Input transcript file from Xenium output, likely named transcripts.csv.gz')
    inout_params.add_argument('--out-dir', required= True, type=str, help='The output directory.')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv", help='The output transcript-indexed SGE file in TSV format. Default: transcripts.unsorted.tsv')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='The output coordinate minmax TSV file. Default: coordinate_minmax.tsv')
    inout_params.add_argument('--out-feature', type=str, default="features.clean.tsv.gz", help='The output files for gene. Default: features.clean.tsv.gz')

>>>>>>> origin/hyun
    incol_params = parser.add_argument_group("Input Columns Parameters", "Input column parameters .")
    incol_params.add_argument('--csv-colname-x',  type=str, default='x_location', help='Column name for X-axis (default: x_location)')
    incol_params.add_argument('--csv-colname-y',  type=str, default='y_location', help='Column name for Y-axis (default: y_location)')
    incol_params.add_argument('--csv-colname-feature-name', type=str, default='feature_name', help='Column name for gene name(default: feature_name)')
    incol_params.add_argument('--csv-colname-phredscore', type=str, default='qv', help='Column name for Phred-scaled quality value (Q-Score) estimating the probability of incorrect call(default: qv)')
    incol_params.add_argument('--csv-colnames-others', nargs='+', default=[], help='Columns names to keep (e.g. cell_id, overlaps_nucleus).')
    # key parameters
    key_params = parser.add_argument_group("Key Parameters", "Key parameters, such as filtering cutoff.")
    key_params.add_argument('--min-phred-score', type=float, default=13, help='Quality score cutoff')
    key_params.add_argument('--dummy-genes', type=str, default='', help='A single name or a regex describing the names of negative control probes')
    key_params.add_argument('--precision-um', type=int, default=-1, help='Number of digits to store the transcript coordinates (only if --px_to_um is in use). Set it to 0 to round to integer. Default is -1, without rounding.')
    # output columns
    outcol_params = parser.add_argument_group("Output Columns Parameters", "Output column parameters for csv.")
    outcol_params.add_argument('--colname-x', type=str, default='X', help='Output Options. Column name for X (default: X)')
    outcol_params.add_argument('--colname-y', type=str, default='Y', help='Output Options. Column name for Y (default: Y)')
    outcol_params.add_argument('--colname-feature-name', type=str, default='gene', help='Output Options. Column name for feature/gene name (default: gene)')
<<<<<<< HEAD
    outcol_params.add_argument('--colnames-count', type=str, default='gn', help='Output Options. Comma-separate column names for Count (default: gn)')
    args = parser.parse_args()
    return args
=======
    outcol_params.add_argument('--colname-count', type=str, default='gn', help='Output Options. Column name for Count (default: gn)')

    args = parser.parse_args(_args)
>>>>>>> origin/hyun

def format_xenium(_args):
    args=parse_arguments(_args)
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
    # transcript
    # - header
    unit_info=[args.colname_x, args.colname_y, args.colname_feature_name] + args.csv_colnames_others
    oheader = unit_info + [args.colnames_count]
    with open(out_transcript_path, 'w') as wf:
        _ = wf.write('\t'.join(oheader)+'\n')
    # - chunks
    for chunk in pd.read_csv(args.input,header=0,chunksize=500000):
        # filter
        chunk = chunk.loc[chunk[args.csv_colname_phredscore] > args.min_phred_score]
        if args.dummy_genes != '':
<<<<<<< HEAD
            chunk = chunk[~chunk[args.csv_colname_feature_name].str.contains(args.dummy_genes, flags=re.IGNORECASE, regex=True)]
=======
            with warnings.catch_warnings(): ## to suppress the warning message for regex input
                warnings.simplefilter("ignore", UserWarning)
                chunk = chunk[~chunk[args.tsv_colname_feature_name].str.contains(args.dummy_genes, flags=re.IGNORECASE, regex=True)]
        
>>>>>>> origin/hyun
        # rename
        chunk.rename(columns = {args.csv_colname_x:args.colname_x, 
                                args.csv_colname_y:args.colname_y, 
                                args.csv_colname_feature_name:args.colname_feature_name}, inplace=True)
        # count
        chunk[args.colnames_count] = 1
        chunk = chunk.groupby(by = unit_info).agg({args.colnames_count:'sum'}).reset_index()
        # conversion (skip given that the input x y are in um, no need to convert.)
        # if args.units_per_um != 1:
        #     chunk[args.csv_colname_x] = chunk[args.csv_colname_x] / args.units_per_um
        #     chunk[args.csv_colname_y] = chunk[args.csv_colname_y] / args.units_per_um
        # write down
        chunk[oheader].to_csv(out_transcript_path, sep='\t',mode='a',index=False,header=False,float_format=float_format)
        # feature
        logging.info(f"{chunk.shape[0]}")
        feature = pd.concat([feature, chunk.groupby(by=args.colname_feature_name).agg({args.colnames_count:"sum"}).reset_index()])
        #Update bounds
        x0, x1 = chunk[args.colname_x].min(), chunk[args.colname_x].max()
        y0, y1 = chunk[args.colname_y].min(), chunk[args.colname_y].max()
        xmin, xmax = min(xmin, x0), max(xmax, x1)
        ymin, ymax = min(ymin, y0), max(ymax, y1)
    # write down
    # - feature
    feature = feature.groupby(by=args.colname_feature_name).agg({args.colnames_count:"sum"}).reset_index()
    feature.to_csv(out_feature_path,sep='\t',index=False)
    # - minmax
    with open(out_minmax_path, 'w') as wf:
        wf.write(f"xmin\t{xmin:.2f}\n")
        wf.write(f"xmax\t{xmax:.2f}\n")
        wf.write(f"ymin\t{ymin:.2f}\n")
        wf.write(f"ymax\t{ymax:.2f}\n")

<<<<<<< HEAD
if __name__ == '__main__':
    format_xenium(sys.argv[1:])
=======
if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
>>>>>>> origin/hyun

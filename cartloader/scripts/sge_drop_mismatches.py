import argparse, os, re, sys, logging, gzip, inspect
import pandas as pd
import numpy as np


from cartloader.utils.utils import run_command

def get_xymin(minmax_path):
    with open(minmax_path, 'r') as file:
        minmax_values = {}
        for line in file:
            key, value = line.strip().split("\t")
            if key == "xmin":
                xmin = float(value)
            elif key == "ymin":
                ymin = float(value)
    return xmin, ymin

def sge_drop_mismatches(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                        Drop mismatches from the transcript and feature files for SeqScope datasets.
                                        Mismatch is defined as x=0 and y=0.
                                        If mismatch exists, the script will drop the mismatches from the transcript and feature files and recalculate the min and max coordinates.
                                        The new files will be saved with the same name as the input files.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-dir', required= True, type=str, help='The input directory')
    inout_params.add_argument('--transcript', type=str, default="transcripts.unsorted.tsv.gz", help='File name of input compressed transcript-indexed SGE file in TSV format from spatula. Default: transcripts.unsorted.tsv.gz')
    inout_params.add_argument('--minmax', type=str, default="coordinate_minmax.tsv", help='File name of input TSV file for min and max X and Y coordinates from spatula. Default: coordinate_minmax.tsv')
    inout_params.add_argument('--feature', type=str, default="features.clean.tsv.gz", help='File name of input compressed TSV file collects UMI counts on a per-gene basis from spatula. Default: features.clean.tsv.gz')

    temp_params = parser.add_argument_group("TEMP Parameters", "Temporary files for processing.")
    temp_params.add_argument('--temp-mismatches', type=str, default="feature.mismatches.tsv.gz", help='Temporary file to store mismatches. Default: feature.mismatches.tsv.gz')
    temp_params.add_argument('--temp-transcript', type=str, default="transcripts.unsorted.withmismatches.tsv.gz", help='Temporary file to store transcripts with mismatches. Default: transcripts.unsorted.withmismatches.tsv.gz')
    temp_params.add_argument('--temp-feature', type=str, default="feature.clean.withmismatches.tsv.gz", help='Temporary file to store features with mismatches. Default: feature.clean.withmismatches.tsv.gz')
    temp_params.add_argument('--temp-minmax', type=str, default="coordinate_minmax.withmismatches.tsv", help='Temporary file to store min and max coordinates with mismatches. Default: coordinate_minmax.withmismatches.tsv')

    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4".')
    args = parser.parse_args(_args)

    # out files
    minmax_path     = os.path.join(args.in_dir, args.minmax)
    transcript_path = os.path.join(args.in_dir, args.transcript)
    feature_path    = os.path.join(args.in_dir, args.feature)

    # check if mismatch exists by xmin and ymin
    minmax_path=os.path.join(args.in_dir, args.minmax)
    xmin, ymin=get_xymin(minmax_path)
    print(f"\t - xmin {xmin};")
    print(f"\t - ymin {ymin};")

    if xmin == 0 and ymin == 0:
        print("Mismatches found, dropping mismatches ...")
        # temp files
        temp_mismatch_path  = os.path.join(args.in_dir, args.temp_mismatches)
        temp_transcript_path = os.path.join(args.in_dir, args.temp_transcript)
        temp_feature_path   = os.path.join(args.in_dir, args.temp_feature)
        temp_minmax_path    = os.path.join(args.in_dir, args.temp_minmax)

        # Move files
        run_command(f"mv {transcript_path} {temp_transcript_path}")
        run_command(f"mv {feature_path} {temp_feature_path}")
        run_command(f"mv {minmax_path} {temp_minmax_path}")

        # Process mismatches
        mismatch_cmd = f"{args.gzip} -cd {temp_transcript_path} | awk '$1==0 && $2==0' | awk '{{sum[$3] += $4}} END {{for (gene in sum) print gene\"\\t\"sum[gene]}}' | {args.gzip} -c > {temp_mismatch_path}"
        run_command(mismatch_cmd)

        # Filter and recompress transcript and feature files
        tsv_cmd = f"{args.gzip} -cd {temp_transcript_path} | awk '$1!=0 || $2!=0' | {args.gzip} -c > {transcript_path}"
        run_command(tsv_cmd)

        ftr_cmd = f"{args.gzip} -cd {temp_mismatch_path} | awk 'BEGIN{{OFS=\"\\t\"}} NR==FNR {{a[$1]=$2; next}} {{if($1 in a) $3=$3-a[$1]; print}}' - <(zcat {temp_feature_path}) | {args.gzip} -c > {feature_path}"
        #ftr_cmd = f"{args.gzip} -cd {temp_mismatch_path} | awk 'NR==FNR {{a[$1]=$2; next}} {{if($1 in a) $3=$3-a[$1]; print}}' - <(zcat {temp_feature_path}) | {args.gzip} -c > {feature_path}"
        run_command(ftr_cmd, use_bash=True)

        # Recalculate min and max
        minmax_cmd = f"{args.gzip} -cd {transcript_path} | awk 'NR==2 {{xmin=$1; xmax=$1; ymin=$2; ymax=$2}} {{if($1<xmin) xmin=$1; if($1>xmax) xmax=$1; if($2<ymin) ymin=$2; if($2>ymax) ymax=$2}} END {{print \"xmin\", xmin; print \"xmax\", xmax; print \"ymin\", ymin; print \"ymax\", ymax}}' OFS=\"\\t\" > {minmax_path}"
        run_command(minmax_cmd)
    else:
        print(f"No mismatch found, skipping the mismatch removal step.")


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
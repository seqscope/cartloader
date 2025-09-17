import argparse, os, sys, gzip, inspect
import subprocess, shlex
from cartloader.utils.utils import run_command

# def detect_feature_count_col_by_header(feature_gz_path: str) -> int:
#     """
#     Read the first non-empty line in the (gzipped) features file as header,
#     and return the 1-based index of the LAST column (count column).
#     Expected header forms:
#       2 columns: gene, count
#       3 columns: gene, gene_id, count
#     """
#     with gzip.open(feature_gz_path, "rt") as f:
#         for line in f:
#             header = line.strip()
#             if not header:
#                 continue
#             col_count = len(header.split("\t"))
#             if col_count < 2:
#                 raise ValueError(
#                     f"Unexpected features header with <2 columns in {feature_gz_path}: {header}"
#                 )
#             return col_count  # 1-based index of last column
#     raise ValueError(f"Empty features file (no header found): {feature_gz_path}")

## SLOW, currently use grep.
# def find_mismatch_transcripts(transcript_gz_path: str, gzip_bin: str = "gzip") -> bool:
#     """
#     Fast check: return True if any data row has x==0 && y==0.
#     Assumes the transcript file always has a header in the first line.
#     """
#     # Skip header with NR>1 and exit on first match
#     awk_prog = r'NR>1 && $1==0 && $2==0 { print 1; exit }'
#     cmd = f'{gzip_bin} -cd {shlex.quote(transcript_gz_path)} | awk \'{awk_prog}\''
#     res = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
#     return (res.stdout or "").strip() == "1"

# def find_mismatch_transcripts(transcript_gz_path: str, gzip_bin: str = "gzip") -> bool:
#     """
#     Return True if any data row has x==0 && y==0.
#     Assumes the transcript file always has a header on the first line.
#     Uses awk exit codes for fastest early-exit behavior.
#     """
#     # Skip header (NR==1), then exit 0 on first match; otherwise exit 1 at EOF
#     awk_prog = r'NR==1{next} { if ($1==0 && $2==0) exit 0 } END{ exit 1 }'
#     # Use LC_ALL=C for slightly faster parsing in awk
#     cmd = f'LC_ALL=C {gzip_bin} -cd {shlex.quote(transcript_gz_path)} | awk -F"\t" \'{awk_prog}\''
#     # Don't capture output; just check returncode. shell=True is required for the pipe.
#     res = subprocess.run(cmd, shell=True)
#     return res.returncode == 0

def find_mismatch_transcripts(transcript_gz_path: str, gzip_bin: str = "gzip") -> bool:
    """
    True if any data row has x==0 && y==0. Assumes header on line 1.
    Uses grep for fastest early-exit scan.
    """
    # Skip header (tail -n +2), anchor to row-start, require tab delimiter
    # -m1: stop after 1 match; -q: quiet (no output)
    cmd = (
        f"LC_ALL=C {gzip_bin} -cd {shlex.quote(transcript_gz_path)} "
        f"| tail -n +2 "
        f"| grep -m1 -q -E '^(0(\\.0+)?)[[:space:]]+(0(\\.0+)?)([[:space:]]|$)'"
    )
    res = subprocess.run(cmd, shell=True)
    return res.returncode == 0

def read_xmin_ymin(minmax_tsv_path: str):
    """Return (xmin, ymin) from coordinate_minmax.tsv, raising if missing."""
    xmin = ymin = None
    with open(minmax_tsv_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, value = line.split("\t")
            if key == "xmin":
                xmin = float(value)
            elif key == "ymin":
                ymin = float(value)
    if xmin is None or ymin is None:
        raise ValueError(f"Missing xmin/ymin in {minmax_tsv_path}")
    return xmin, ymin

def sge_drop_mismatches(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Removes entries with coordinates x=0, y=0 from SeqScope transcript and feature files. These mismatches occur when reads fail to map to a valid spatial position, often due to barcode or anchor errors, and are incorrectly assigned to the origin. 
                                     When detected, the script filters them from the transcript and feature tables, subtracts their counts from per-gene totals, and recalculates the coordinate ranges. 
                                     The cleaned files are written back under the original filenames, replacing the mismatched versions.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-dir', required= True, type=str, help='The input directory')
    inout_params.add_argument('--transcript', type=str, default="transcripts.unsorted.tsv.gz", help='File name of input compressed transcript-indexed SGE file in TSV format from spatula. Default: transcripts.unsorted.tsv.gz')
    inout_params.add_argument('--minmax', type=str, default="coordinate_minmax.tsv", help='File name of input TSV file for min and max X and Y coordinates from spatula. Default: coordinate_minmax.tsv')
    inout_params.add_argument('--feature', type=str, default="features.clean.tsv.gz", help='File name of input compressed TSV file collects UMI counts on a per-gene basis from spatula. Default: features.clean.tsv.gz')

    temp_params = parser.add_argument_group("TEMP Parameters", "Temporary file names.")
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

    def ensure_exists(paths):
        for p in paths:
            if not os.path.exists(p):
                raise FileNotFoundError(p)

    ensure_exists([transcript_path, feature_path, minmax_path])

    # temp files
    temp_mismatch_path  = os.path.join(args.in_dir, args.temp_mismatches)
    temp_transcript_path = os.path.join(args.in_dir, args.temp_transcript)
    temp_feature_path   = os.path.join(args.in_dir, args.temp_feature)
    temp_minmax_path    = os.path.join(args.in_dir, args.temp_minmax)

    # Read min/min to find if xmin==0 and ymin==0
    xmin, ymin = read_xmin_ymin(minmax_path)
    print(f"\t- xmin: {xmin}")
    print(f"\t- ymin: {ymin}")

    # update: NovaScope & Scopeflow combined outputs may set xmin/ymin to 0 even without true (0,0) rows.
    # Use a stricter check: actually look for any (x==0 && y==0) records in the transcript file.
    has_mismatch = False
    # update_minmax = False 
    if xmin == 0 and ymin == 0:
        if find_mismatch_transcripts(transcript_path):
            has_mismatch = True
        # else:
        #     update_minmax = True
        
    if has_mismatch == True:
        print("Mismatches found, dropping mismatches ...")

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
    # elif update_minmax == True:
    #     print("No mismatches found, but xmin/ymin are 0. Recomputing min/max …")
    #     run_command(f"mv {minmax_path} {temp_minmax_path}")
    #     # Recalculate min and max
    #     minmax_cmd = f"{args.gzip} -cd {transcript_path} | awk 'NR==2 {{xmin=$1; xmax=$1; ymin=$2; ymax=$2}} {{if($1<xmin) xmin=$1; if($1>xmax) xmax=$1; if($2<ymin) ymin=$2; if($2>ymax) ymax=$2}} END {{print \"xmin\", xmin; print \"xmax\", xmax; print \"ymin\", ymin; print \"ymax\", ymax}}' OFS=\"\\t\" > {minmax_path}"
    #     run_command(minmax_cmd)
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
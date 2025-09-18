import argparse, os, sys, gzip, inspect
import subprocess, shlex
from cartloader.utils.utils import parquet_to_csv_with_polars_pigz

def parquet_to_csv_rapid(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Rapid conversion from parquet to compressed CSV using polars and pigz.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-parquet', required= True, type=str, help='The input parquet file')
    inout_params.add_argument('--out-csv-gz', type=str, required=True, help='Output compressed CSV file name (gzip).')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Other parameters")
    aux_params.add_argument('--threads', type=int, default=24, help='Number of threads for pigz. Default: 24')
    aux_params.add_argument('--pigz', type=str, default="pigz", help='Path to pigz binary. Default: pigz')
    aux_params.add_argument('--batch-size', type=int, default=131072, help='Batch size for processing. Default: 131072')
    aux_params.add_argument('--compress-level', type=int, default=6, help='Compression level for pigz. Default: 6')

    args = parser.parse_args(_args)

    parquet_to_csv_with_polars_pigz(
        parquet_file=args.in_parquet,
        out_path=args.out_csv_gz,
        batch_size=args.batch_size,
        compress_level=args.compression_level,
        pigz_path=args.pigz,
        pigz_threads=args.threads
    )


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_tsv2pmtiles", description="Split and convert transcripts TSV file into pmtiles")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Run all commands (split, convert)')
    cmd_params.add_argument('--split', action='store_true', default=False, help='Split molecules TSV file into group-wise CSVs')
    cmd_params.add_argument('--convert', action='store_true', default=False, help='Convert the CSV files into pmtiles')
    cmd_params.add_argument('--clean', action='store_true', default=False, help='Clean intermedate files')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-molecules', type=str, help='Input Long Format TSV/CSV (possibly gzipped) file containing the X/Y coordinates and gene expression counts per spot')
    inout_params.add_argument('--in-features', type=str, help='Input TSV/CSV (possibly gzipped) file containing the gene name and total count for each gene')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='The output prefix. New directory will be created if needed')  
    inout_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    inout_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    inout_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--bin-count', type=int, default=50, help='Number of bins to equally divide the genes into (default: 50)')
#    key_params.add_argument('--log2-multiplier', type=float, default=1.0, help='Multiplier used to determine the log2 bin. Bin is determined as [multiplier] * log2(total_count). Default is 1.0')
    key_params.add_argument('--dummy-genes', type=str, default='', help='A single name or a regex describing the names of negative control probes')
    key_params.add_argument('--chunk-size', type=int, default=1000000, help='Number of rows to read at a time. Default is 1000000')
    key_params.add_argument('--skip-original', action='store_true', default=False, help='Skip writing the original file')
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    conv_params = parser.add_argument_group("Parameters for pmtiles conversion")
    conv_params.add_argument('--min-zoom', type=int, default=10, help='Minimum zoom level')
    conv_params.add_argument('--max-zoom', type=int, default=18, help='Maximum zoom level')
    conv_params.add_argument('--max-tile-bytes', type=int, default=1500000, help='Maximum bytes for each tile in PMTiles')
    conv_params.add_argument('--max-feature-counts', type=int, default=200000, help='Max feature limits per tile in PMTiles')
    conv_params.add_argument('--preserve-point-density-thres', type=int, default=16384, help='Threshold for preserving point density in PMTiles')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--in-molecules-delim', type=str, default='\t', help='Delimiter used in the input molecules files. Default is tab.')  
    aux_params.add_argument('--in-features-delim', type=str, default='\t', help='Delimiter used in the input feature files. Default is tab.')  
    aux_params.add_argument('--out-molecules-delim', type=str, default=',', help='Delimiter used in the output molecule TSV/CSV files. Default is ,')
    aux_params.add_argument('--out-features-delim', type=str, default='\t', help='Delimiter used in the output feature files. Default is tab.')
    aux_params.add_argument('--out-molecules-suffix', type=str, default="molecules.csv", help='The output file to store individual molecule count matrix. Each file name will be [out-prefix].split.[bin_id].[out-molecules-suffix]. Default: molecules.tsv.gz')
    aux_params.add_argument('--out-features-suffix', type=str, default="features.tsv.gz", help='The output file to store feature count matrix. Each file name will be [out-prefix].split.[bin_id].[out-features-suffix]. Default: features.tsv.gz')
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_tsv2pmtiles(_args):
    """
    Convert cross-platform molecules (TSV) file to pmtiles
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_tsv2pmtiles" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    if args.all:
        args.split = True
        args.convert = True
        args.clean = True

    ## rename X/Y to lon/lat for tippecanoe compatibility
    if args.col_rename is None:
        args.col_rename = []
    args.col_rename.append("X:lon")
    args.col_rename.append("Y:lat")

    # start mm
    mm = minimake()

    # create output directory if needed
    out_dir = os.path.dirname(args.out_prefix)
    out_base = os.path.basename(args.out_prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    
    # 1. Perform split without running makefile
    if args.split:
        logger.info("Splitting the input cross-platform TSV file into CSV files")

        # cmd = [
        #     "cartloader", "split_molecule_counts",
        #     "--in-molecules", args.in_molecules,
        #     "--in-features", args.in_features,
        #     "--out-prefix", args.out_prefix,
        #     "--in-molecule-delim", repr(args.in_molecules_delim),
        #     "--in-features-delim", repr(args.in_features_delim),
        #     "--out-molecule_suffix", args.out_molecules_suffix,
        #     "--out-molecule-delim", repr(args.out_molecules_delim),
        #     "--out-features-suffix", args.out_features_suffix,
        #     "--out-features-delim", repr(args.out_features_delim),
        #     "--colname-feature", args.colname_feature,
        #     "--colname-count", args.colname_count,
        #     "--col_rename", " ".join(args.col_rename),
        #     "--bin_multiplier", str(args.bin_multiplier),
        #     "--chunk-size", str(args.chunk_size),
        #     "--log-suffix", args.log_suffix
        # ]

        # if args.dummy_genes != '':
        #     cmd.append("--dummy-genes", args.dummy_genes)
        # if args.log:
        #     cmd.append("--log")
        # if args.skip_original:
        #     cmd.append("--skip-original")

        cmd = f"""cartloader split_molecule_counts \\
                --in-molecules {args.in_molecules} \\
                --in-features {args.in_features} \\
                --out-prefix {args.out_prefix} \\
                --in-molecules-delim \"{args.in_molecules_delim}\" \\
                --in-features-delim \"{args.in_features_delim}\" \\
                --out-molecules-suffix {args.out_molecules_suffix} \\
                --out-molecules-delim \"{args.out_molecules_delim}\" \\
                --out-features-suffix {args.out_features_suffix} \\
                --out-features-delim \"{args.out_features_delim}\" \\
                --colname-feature {args.colname_feature} \\
                --colname-count {args.colname_count} \\
                --col-rename {" ".join(args.col_rename)} \\
                --bin-count {args.bin_count} \\
                --equal-bins \\
                --dummy-genes \"{args.dummy_genes}\" \\
                --chunk-size {args.chunk_size} \\
                --log-suffix {args.log_suffix} \\
        """ + ("--log" if args.log else "") + ("--skip-original" if args.skip_original else "")

        print(cmd)
        result = subprocess.run(cmd, shell=True)

#        cmd_str = " ".join(cmd)
#        result = subprocess.run(cmd_str, shell=True)
        if result.returncode != 0:
            logger.error("Error in splitting the input TSV file into CSV files")
            sys.exit(1)

    # 2. Perform conversion:
    if args.convert:
        ## open index file
        df = pd.read_csv(f"{args.out_prefix}_index.tsv", sep="\t")

        ## add targets for each bin
        for i, row in df.iterrows():
            bin_id = row["bin_id"]
            csv_path = out_dir + "/" + row["molecules_path"]
            if bin_id == "all":
                pmtiles_path = args.out_prefix + "_all.pmtiles"
            else:
                pmtiles_path = args.out_prefix + "_bin" + bin_id + ".pmtiles"
            cmds = cmd_separator([], f"Converting bin {bin_id} to pmtiles")
            cmds.append(f"{args.tippecanoe} -o {pmtiles_path} -Z {args.min_zoom} -z {args.max_zoom} --force -s EPSG:3857 -M {args.max_tile_bytes} --drop-densest-as-needed --extend-zooms-if-still-dropping '--preserve-point-density-threshold={args.preserve_point_density_thres}' --no-duplication --no-clipping --no-tile-size-limit --buffer 0 {csv_path}")
            mm.add_target(pmtiles_path, [csv_path], cmds)

        if len(mm.targets) == 0:
            logging.error("There is no target to run. Please make sure that at least one run option was turned on")
            sys.exit(1)

        ## write makefile
        mm.write_makefile(f"{args.out_prefix}.Makefile")

        logger.info("Running makefile to convert the CSV files to pmtiles")

        result = subprocess.run(f"make -f {args.out_prefix}.Makefile -j {args.n_jobs} {'-B' if args.restart else ''}", shell=True)
        if result.returncode != 0:
            logger.error("Error in splitting the input TSV file into CSV files")
            sys.exit(1)

    # 3. clean the intermediate files and write new output files
    logger.info("Cleaning intermediate files")

    df = pd.read_csv(f"{args.out_prefix}_index.tsv", sep="\t")
    for i, row in df.iterrows():
        bin_id = row["bin_id"]
        csv_path = out_dir + "/" + row["molecules_path"]
        if args.clean and os.path.exists(csv_path):
            os.remove(csv_path)
        ftr_path = out_dir + "/" + row["features_path"]
        if args.clean and os.path.exists(ftr_path):
            os.remove(ftr_path)

    df_out = df.drop(columns=['molecules_path','features_path'])
    df_out['pmtiles_path'] = [f"{out_base}_all.pmtiles" if x == "all" else f"{out_base}_bin{x}.pmtiles" for x in df_out['bin_id']]
    df_out.to_csv(f"{args.out_prefix}_pmtiles_index.tsv", sep="\t", index=False)
    
    logger.info("Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, inspect
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, execute_makefile

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Split and convert transcripts TSV file into pmtiles")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Run all commands (split, convert, clean)')
    cmd_params.add_argument('--split', action='store_true', default=False, help='Split molecules TSV file into group-wise CSVs')
    cmd_params.add_argument('--convert', action='store_true', default=False, help='Convert the CSV files into pmtiles')
    cmd_params.add_argument('--clean', action='store_true', default=False, help='Clean intermediate files')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-molecules', type=str, help='Input Long Format TSV/CSV (possibly gzipped) file containing the X/Y coordinates and gene expression counts per spot')
    inout_params.add_argument('--in-features', type=str, help='Input TSV/CSV (possibly gzipped) file containing the gene name and total count for each gene')
    inout_params.add_argument('--in-bin-json', type=str, default=None, help='Optional precomputed gene->bin assignment JSON (spatula assign-feature2bin output, i.e. a _bin_counts.json). When provided, the gene-to-bin assignment is reused from this file instead of being derived from --in-features, so assignment is identical across all datasets that share it. A copy is written to <out-prefix>_bin_counts.json.')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='The output prefix. New directory will be created if needed')  
    inout_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    inout_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    inout_params.add_argument('--col-rename', type=str, nargs='+', help='Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...')
    # Column names AS THEY APPEAR IN --in-molecules. split-mol2bin looks its columns up
    # by the input name (--col-rename only rewrites the output header), so an input whose
    # header differs from the defaults (X/Y/gene) must name them here. E.g. a punkst tiled
    # TSV has the header "#X Y Feature count", so its X/Y/count already match but its
    # feature column needs --in-colname-feature Feature. Unset = leave split-mol2bin's
    # own defaults in place.
    inout_params.add_argument('--in-colname-x', type=str, default=None, help='Column name for X in --in-molecules (default: split-mol2bin default, X)')
    inout_params.add_argument('--in-colname-y', type=str, default=None, help='Column name for Y in --in-molecules (default: split-mol2bin default, Y)')
    inout_params.add_argument('--in-colname-feature', type=str, default=None, help='Column name for the feature/gene in --in-molecules (default: Feature with --use-pmpoint, else split-mol2bin default, gene)')

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
    conv_params.add_argument('--max-tile-bytes', type=int, default=5000000, help='Maximum bytes for each tile in PMTiles')
    conv_params.add_argument('--max-feature-counts', type=int, default=500000, help='Max feature limits per tile in PMTiles')
    conv_params.add_argument('--preserve-point-density-thres', type=int, default=1024, help='Threshold for preserving point density in PMTiles')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--in-molecules-delim', type=str, default='\t', help='Delimiter used in the input molecules files. Default is tab.')  
    aux_params.add_argument('--in-features-delim', type=str, default='\t', help='Delimiter used in the input feature files. Default is tab.')  
    aux_params.add_argument('--out-molecules-delim', type=str, default=',', help='Delimiter used in the output molecule TSV/CSV files. Default is ,')
    aux_params.add_argument('--out-features-delim', type=str, default='\t', help='Delimiter used in the output feature files. Default is tab.')
    aux_params.add_argument('--out-molecules-suffix', type=str, default="molecules.csv", help='The output file to store individual molecule count matrix. Each file name will be [out-prefix].split.[bin_id].[out-molecules-suffix]. Default: molecules.tsv.gz')
    aux_params.add_argument('--out-features-suffix', type=str, default="features.tsv.gz", help='The output file to store feature count matrix. Each file name will be [out-prefix].split.[bin_id].[out-features-suffix]. Default: features.tsv.gz')
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--pmpoint', type=str, default=f"{repo_dir}/submodules/pmpoint/bin/pmpoint", help='Path to pmpoint binary')
    aux_params.add_argument('--spatula', type=str, default=f"{repo_dir}/submodules/spatula/bin/spatula", help='Path to spatula binary')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--use-pmpoint', action='store_true', default=False, help='Use pmpoint/MLT instead of tippecanoe for point PMTiles generation (requires --pmpoint)')
    aux_params.add_argument('--tile-format-pmpoint', choices=['MLT', 'MVT'], default='MVT', help='Tile format to use when --use-pmpoint is enabled (default: MVT)')
    aux_params.add_argument('--pmpoint-compression-scale', type=float, default=10.0, help='Additional compression scale for pmpoint when --use-pmpoint is turned on. Default: 10.0')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory to be used (default: out-dir/tmp; specify /tmp if needed)')

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

    if args.use_pmpoint and "Feature:gene" not in args.col_rename:
        args.col_rename.append("Feature:gene")

    # start mm
    mm = minimake()

    # create output directory if needed
    out_dir = os.path.dirname(args.out_prefix)
    out_base = os.path.basename(args.out_prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(out_dir, "tmp")
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir, exist_ok=True)
    
    # 1. Perform split without running makefile
    if args.split:
        logger.info("Splitting the input cross-platform TSV file into per-bin files")

        in_colname_feature = args.in_colname_feature or ("Feature" if args.use_pmpoint else None)
        pmpoint_arg = f"--colname-feature {in_colname_feature}" if in_colname_feature else ""
        if args.in_colname_x:
            pmpoint_arg += f" --colname-x {args.in_colname_x}"
        if args.in_colname_y:
            pmpoint_arg += f" --colname-y {args.in_colname_y}"
        col_rename_arg = ""
        if args.col_rename is not None and len(args.col_rename) > 0:
            for col_rename in args.col_rename:
                col_rename_arg += f" --col-rename {col_rename}"

        # The gene->bin assignment is produced by 'assign-feature2bin' and consumed by
        # 'split-mol2bin'. When --in-bin-json is given (e.g. a shared assignment from
        # run_cartload2_multi), reuse it so gene-to-bin assignment is identical across
        # samples; otherwise derive it from this dataset's own feature counts.
        bin_json = f"{args.out_prefix}_bin_counts.json"
        if args.in_bin_json is not None:
            logger.info(f"Reusing gene->bin assignment from {args.in_bin_json}")
            if os.path.abspath(args.in_bin_json) != os.path.abspath(bin_json):
                shutil.copyfile(args.in_bin_json, bin_json)
        else:
            assign_cmd = f"""'{args.spatula}' assign-feature2bin \\
                    --feature-tsv '{args.in_features}' \\
                    --out-json '{bin_json}' \\
                    --bin-count {args.bin_count} \\
                    --in-feature-tsv-delim '{args.in_features_delim}'
            """
            print(assign_cmd)
            result = subprocess.run(assign_cmd, shell=True)
            if result.returncode != 0:
                logger.error("Error in assigning features to bins (assign-feature2bin)")
                sys.exit(1)

        cmd = f"""'{args.spatula}' split-mol2bin \\
                --mol-tsv '{args.in_molecules}' \\
                --bin-json '{bin_json}' \\
                --out-prefix '{args.out_prefix}' \\
                --in-mol-tsv-delim '{args.in_molecules_delim}' \\
                --out-mol-tsv-delim '{args.out_molecules_delim}' \\
                --out-feature-tsv-delim '{args.out_features_delim}' \\
                --out-mol-suffix '{args.out_molecules_suffix}' \\
                --out-feature-suffix '{args.out_features_suffix}' {pmpoint_arg} {col_rename_arg} \\
                """ + ("--skip-original" if args.skip_original else "")

        print(cmd)
        result = subprocess.run(cmd, shell=True)

        if result.returncode != 0:
            logger.error("Error in splitting the input TSV file into per-bin files")
            sys.exit(1)

    # pmpoint reads the SPLIT output, whose header has already been rewritten by
    # --col-rename, so its X/Y column names are the renamed ones.
    def _renamed(name):
        for r in args.col_rename:
            old, sep, new = r.partition(":")
            if sep and old == name:
                return new
        return name
    pmpoint_colname_x = _renamed(args.in_colname_x or "X")
    pmpoint_colname_y = _renamed(args.in_colname_y or "Y")

    # 2. Perform conversion:
    if args.convert:
        ## open index file
        df = pd.read_csv(f"{args.out_prefix}_index.tsv", sep="\t")

        ## add targets for each bin
        for i, row in df.iterrows():
            bin_id = row["bin_id"]
            csv_path = out_dir + "/" + row["molecules_path"]
            if bin_id == "all":
                #pmtiles_path = args.out_prefix + "_all.pmtiles"
                pmtiles_prefix = args.out_prefix + "_all"
            else:
                #pmtiles_path = args.out_prefix + "_bin" + bin_id + ".pmtiles"
                pmtiles_prefix = args.out_prefix + "_bin" + bin_id 
            cmds = cmd_separator([], f"Converting bin {bin_id} to pmtiles")
            if args.use_pmpoint:
                cmds.append(f"mkdir -p {args.tmp_dir}/{bin_id}")
                cmds.append(f"'{args.pmpoint}' build-point-pmtiles --tmp-dir {args.tmp_dir}/{bin_id} --in {csv_path} --out {pmtiles_prefix}.z{args.max_zoom}.pmtiles --zoom {args.max_zoom} --colname-x {pmpoint_colname_x} --colname-y {pmpoint_colname_y} --delim ',' --threads {args.threads} --format {args.tile_format_pmpoint}")
                cmds.append(f"'{args.pmpoint}' build-pyramid-pmtiles --scale-factor-compression {args.pmpoint_compression_scale} --tmp-dir {args.tmp_dir}/{bin_id} --in {pmtiles_prefix}.z{args.max_zoom}.pmtiles --out {pmtiles_prefix}.pmtiles --min-zoom {args.min_zoom} --max-tile-bytes {args.max_tile_bytes} --max-tile-features {args.max_feature_counts} --threads {args.threads}")
                cmds.append(f"rm {pmtiles_prefix}.z{args.max_zoom}.pmtiles")
                cmds.append(f"rm -rf {args.tmp_dir}/{bin_id}")
            else:
                cmds.append(f"TIPPECANOE_MAX_THREADS={args.threads} '{args.tippecanoe}' -t {args.tmp_dir} -o {pmtiles_prefix}.pmtiles -Z {args.min_zoom} -z {args.max_zoom} --force -s EPSG:3857 -M {args.max_tile_bytes} -O {args.max_feature_counts} --drop-densest-as-needed --extend-zooms-if-still-dropping '--preserve-point-density-threshold={args.preserve_point_density_thres}' --no-duplication --no-clipping --buffer 0 {csv_path}")
            mm.add_target(f"{pmtiles_prefix}.pmtiles", [csv_path], cmds)

        if len(mm.targets) == 0:
            logging.error("There is no target to run. Please make sure that at least one run option was turned on")
            sys.exit(1)

        ## write makefile
        make_f = f"{args.out_prefix}.mk"
        mm.write_makefile(make_f)

        logger.info("Running makefile to convert the CSV files to pmtiles")

        execute_makefile(make_f, dry_run=False, restart=args.restart, n_jobs=args.n_jobs)


    # 3. clean the intermediate files and write new output files
    if not args.keep_intermediate_files:
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

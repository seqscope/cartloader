import sys, os, gzip, argparse, logging, warnings, shutil
from ficture.utils.minimake import minimake

def parse_arguments(_args):
    if len(_args) == 0:
        parser.print_help()
        return
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    # run params
    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default="Makefile", help='The name of the Makefile to generate')
    # basic params
    bas_params = parser.add_argument_group("Basic Parameters", "Basic parameters for the FICTURE pipeline.")
    bas_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope"], help='Platform of the raw input file to infer the format of the input file.')
    bas_params.add_argument('--dummy-genes', type=str, default="BLANK\|NegCon\|NegPrb", help='Name of the negative controls, could pass regex to match multiple name patterns. Required platform(s): cosmx_smi; 10x_xenium; vizgen_merscope. ')
    bas_params.add_argument('--units-per-um', type=float, default=1.00, help='Output Options. Coordinate unit per um (conversion factor) (default: 1.00). Required platform(s):') 
    bas_params.add_argument('--precision-um', type=int, default=2, help='Number of digits to store the transcript coordinates (only if --px_to_um is in use). Set it to 0 to round to integer. Required platform(s): cosmx_smi; vizgen_merscope. Default: 2')
    # input params
    input_params = parser.add_argument_group("Input Parameters", "Input parameters, including platform information and the input files.")
    input_params.add_argument('--in-sge', type=str, default=None, help='Path to input raw sge directory. Required platform(s): 10x_visium_hd')
    input_params.add_argument('--sge-bcd', type=str, default="barcodes.tsv.gz", help='(Optional) Barcode file name in SGE directory. Default: barcodes.tsv.gz. Required platform(s): 10x_visium_hd.')
    input_params.add_argument('--sge-ftr', type=str, default="features.tsv.gz", help='(Optional)Feature file name in SGE directory. Default: features.tsv.gz. Required platform(s): 10x_visium_hd.')
    input_params.add_argument('--sge-mtx', type=str, default="matrix.mtx.gz", help='(Optional) Matrix file name in SGE directory. Default: matrix.mtx.gz. Required platform(s): 10x_visium_hd.')
    input_params.add_argument('--in-parquet', type=str, default=None, help='Path to input raw parquet file for spatial coordinates in parquet format, e.g., tissue_positions.parquet. Required platform(s): 10x_visium_hd.')
    input_params.add_argument('--in-json', type=str, default=None, help='Path to input scale factor json file, e.g., scalefactors_json.json. Required platform(s): 10x_visium_hd.')
    input_params.add_argument('--in-csv', type=str, default=None, help='Path to input raw csv/tsv file. Required platform(s): 10x_xenium; cosmx_smi.')
    # output params
    output_params=parser.add_argument_group("Output Parameters", "Output parameters for the FICTURE pipeline.")
    output_params.add_argument('--out-dir', type=str, help='The output directory.')
    output_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='The output compressed transcript-indexed SGE file in TSV format. Default: transcripts.unsorted.tsv.gz')
    output_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='The output coordinate minmax TSV file. Default: coordinate_minmax.tsv')
    output_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='The output files for gene. Default: feature.clean.tsv.gz')
    # AUX input SGE params
    aux_in_sge_params = parser.add_argument_group("SGE Auxiliary Parameters", "Auxiliary Parameters for SGE input and additional barcode position input.")
    aux_in_sge_params.add_argument('--icols-mtx', type=str, default='1,2,3,4,5', help='Expected column index in SGE input. Comma-separated 1-based column indices use as the count (default: 1,2,3,4,5)')
    aux_in_sge_params.add_argument('--icol-bcd-barcode', type=int, default=1, help='Expected column index in SGE input. 1-based column index of barcode in the barcode file (default: 1)')
    aux_in_sge_params.add_argument('--icol-bcd-x', type=int, default=6, help='Expected column index in SGE input. 1-based column index of x coordinate in the barcode file (default: 6)')
    aux_in_sge_params.add_argument('--icol-bcd-y', type=int, default=7, help='Expected column index in SGE input. 1-based column index of y coordinate in the barcode file (default: 7)')
    aux_in_sge_params.add_argument('--icol-ftr-id', type=int, default=1, help='Expected column index in SGE input. 1-based column index of feature ID in the barcode file (default: 1)')
    aux_in_sge_params.add_argument('--icol-ftr-name', type=int, default=2, help='Expected column index in SGE input. 1-based column index of feature name in the barcode file (default: 2)')
    aux_in_sge_params.add_argument('--pos-colname-barcode', type=str, default='barcode', help='Additional Barcode Position File. Column name for barcode in the position file (default: barcode)')
    aux_in_sge_params.add_argument('--pos-colname-x', type=str, default='pxl_row_in_fullres', help='Additional Barcode Position File. Column name for X-axis in the position file (default: pxl_row_in_fullres)')
    aux_in_sge_params.add_argument('--pos-colname-y', type=str, default='pxl_col_in_fullres', help='Additional Barcode Position File. Column name for Y-axis in the position file (default: pxl_col_in_fullres)')
    aux_in_sge_params.add_argument('--pos-delim', type=str, default=',', help='Additional Barcode Position File. Delimiter for the position file (default: ",")')
    # AUX input csv params
    aux_in_csv_params = parser.add_argument_group("CSV Auxiliary Parameters", "Auxiliary Parameters TSV/CSV input.")
    aux_in_csv_params.add_argument('--csv-colname-x',  type=str, default='x_location', help='Column name for X-axis (default: x_location)')
    aux_in_csv_params.add_argument('--csv-colname-y',  type=str, default='y_location', help='Column name for Y-axis (default: y_location)')
    aux_in_csv_params.add_argument('--csv-colname-feature-name', type=str, default='feature_name', help='Column name for gene name(default: feature_name)')
    aux_in_csv_params.add_argument('--csv-colname-phredscore', type=str, default='qv', help='Column name for Phred-scaled quality value (Q-Score) estimating the probability of incorrect call(default: qv)')
    aux_in_csv_params.add_argument('--csv-colnames-others', nargs='+', default=[], help='Columns names to keep.e.g., cell_id, overlaps_nucleus')
    aux_in_csv_params.add_argument('--csv-colnames-count', type=str, default="MIDCounts", help='Column name for gene id (e.g.: transcript_id)')
    # AUX output params
    aux_output_params = parser.add_argument_group("Shared Output Auxiliary Parameters", "Auxiliary Shared Output column index parameters for the FICTURE pipeline.")
    aux_output_params.add_argument('--colname-x', type=str, default='X', help='Column name for X (default: X)')
    aux_output_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y (default: Y)')
    aux_output_params.add_argument('--colnames-count', type=str, default='gn', help='Comma-separate column names for Count (default: gn)')
    aux_output_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for feature/gene name (default: gene)')
    aux_output_params.add_argument('--colname-feature-id', type=str, default='gene_id', help='Column name for feature/gene ID (default: gene_id)')
    aux_output_params.add_argument('--colname-molecule-id', type=str, default='MoleculeID', help='Column name for MoleculeID (default: MoleculeID)')
    # Aux other params
    aux_params = parser.add_argument_group("Shared Auxiliary Parameters", "Auxiliary parameters.")
    aux_params.add_argument('--unique', action='store_true', default=False, help='Merge pixels with (almost?) identical coordinates. Applies to cosmx_smi only.')
    aux_params.add_argument('--include-feature-list', type=str, default=None, help='Input Filtering Options. A file containing a list of input genes to be included (feature name of IDs) (default: "")')
    aux_params.add_argument('--exclude-feature-list', type=str, default=None, help='Input Filtering Options. A file containing a list of input genes to be excluded (feature name of IDs) (default: "")')
    aux_params.add_argument('--include-feature-substr', type=str, default=None, help='Input Filtering Options. A substring of feature/gene names to be included (default: "")')
    aux_params.add_argument('--exclude-feature-substr', type=str, default=None, help='Input Filtering Options. A substring of feature/gene names to be excluded (default: "")')
    aux_params.add_argument('--include-feature-regex', type=str, default=None, help='Input Filtering Options. A regex pattern of feature/gene names to be included (default: "")')
    aux_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Input Filtering Options. A regex pattern of feature/gene names to be excluded (default: "")')
    aux_params.add_argument('--print-feature-id', action='store_true', help='Output Options. Print feature ID in output file (default: OFF)')
    aux_params.add_argument('--allow-duplicate-gene-names', action='store_true', help='Output Options. Allow duplicate gene names in the output file (default: OFF)')
    # feature_params = parser.add_argument_group("Filtering Feature Parameters", "Filtering feature paramter using a reference gene info file.")
    # feature_params.add_argument('--filter-feature', action='store_true', default=False, help='Filter the feature using the gene information file.')
    # feature_params.add_argument('--geneinfo', type=str, default=None, help='Gene information file. Required when --filter-feature is used.')
    # feature_params.add_argument('--kept-gene-type', type=str, default="protein_coding,lncRNA", help='Gene type to keep. Default: protein_coding,lncRNA')
    # feature_params.add_argument('--rm-gene-regex', type=str, default="^Gm\d+|^mt-|^MT-", help='Gene regex to remove. Default: ^Gm\d+|^mt-|^MT-')
    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to {args.gzip} binary. For faster processing, use "pigz -p 4".')
    env_params.add_argument('--sort', type=str, default="sort", help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G".')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary.')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools binary.')
    args= parser.parse_args()
    return args

def convert_in_by_platform(args):
    if args.platform == "10x_visium_hd":
        required_files = [os.path.join(args.in_sge, f) for f in [args.sge_bcd, args.sge_ftr, args.sge_mtx]] + [args.in_parquet, args.in_json]
    elif args.platform == "10x_xenium":
        required_files = [args.in_csv]
    elif args.platform == "bgi_stereoseq":
        required_files = [args.in_csv]
    elif args.platform == "cosmx_smi":
        required_files = [args.in_csv]
    elif args.platform == "vizgen_merscope":
        required_files = [args.in_csv]
    return required_files

def create_minmax(cmds, args):
    # removed the px-to-um conversion.
    minmax_cmd = f"""{args.gzip} -cd {args.out_dir}/{args.out_transcript} | awk 'BEGIN{{FS=OFS="\\t"}} NR==1{{for(i=1;i<=NF;i++){{if($i=="X")x=i;if($i=="Y")y=i}}print $x,$y;next}}{{print $x,$y}}' | awk -F'\\t' ' BEGIN {{ min1 = "undef"; max1 = "undef"; min2 = "undef"; max2 = "undef"; }} {{ if (NR == 2 || $1 < min1) min1 = $1; if (NR == 2 || $1 > max1) max1 = $1; if (NR == 2 || $2 < min2) min2 = $2; if (NR == 2 || $2 > max2) max2 = $2; }} END {{ print "xmin\\t", min1; print "xmax\\t", max1; print "ymin\\t", min2; print "ymax\\t", max2; }}' > {args.out_dir}/{args.out_minmax}"""
    cmds.append(minmax_cmd)
    return cmds

def add_param_to_cmd(cmd, args, aux_argset):
    aux_args = {k: v for k, v in vars(args).items() if k in aux_argset}
    for arg, value in aux_args.items():
        if value or isinstance(value, bool):
            arg_name = arg.replace('_', '-')
            if isinstance(value, bool) and value:
                cmd += f" --{arg_name}"
            elif isinstance(value, list):
                cmd += f" --{arg_name} {' '.join(value)}"
            elif not isinstance(value, bool):
                cmd += f" --{arg_name} {value}"
    return cmd

def convert_visiumhd(cmds, args):
    scheck_app(args.spatula)
    # input: in_sge, in_parquet, in_json
    # output: out_transcript, out_minmax, out_feature, (out_sge)
    tmp_parquet = f"{args.out_dir}/tissue_positions.csv.gz"
    # 1) convert parquet to csv
    cmds.append(f"{args.parquet_tools} csv {args.in_parquet} |  {args.gzip} -c > {tmp_parquet}")
    # 2) convert sge to tsv (output: out_transcript, out_minmax, out_feature, (optional) out_sge)
    #   - if out_sge is provided, spatula will generate sge files
    format_cmd=f"{args.spatula} convert-sge --in-sge {args.in_sge} --out-tsv {args.out_dir} --pos {tmp_parquet} --tsv-mtx {args.out_transcript} --tsv-ftr {args.out_feature} --tsv-minmax {args.out_minmax}"
    aux_argset = {
        'sge_bcd', 'sge_ftr', 'sge_mtx', 'out_sge',
        'icols_mtx', 'icol_bcd_barcode', 'icol_bcd_x', 'icol_bcd_y', 'icol_ftr_id', 'icol_ftr_name',
        'pos_colname_barcode', 'pos_colname_x', 'pos_colname_y', 'pos_delim',
        'include_feature_list', 'exclude_feature_list', 'include_feature_substr', 'exclude_feature_substr', 'include_feature_regex', 'exclude_feature_regex', 
        'units_per_um', 'precision_um',
        'print_feature_id', 'allow_duplicate_gene_names', 'colname_feature_name', 'colname_feature_id', 'colnames_count', 'colname_x', 'colname_y'
    }
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    cmds.append(f"rm {tmp_parquet}")
    return cmds

def convert_xenium(cmds, args):
    # input: in_csv
    # output: out_transcript, out_minmax, out_feature
    transcript_tsv = args.out_transcript.replace(".gz", "")
    format_cmd = f"python {cartloader}/scripts/format_xenium_custom.py --input {args.in_csv} --out-dir {args.out_dir} --out-transcript {transcript_tsv} --out-feature {args.out_feature} --out-minmax {args.out_minmax}"
    aux_argset = {
        'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name', 'csv_colname_phredscore', 'csv_colnames_others',
        'min_phred_score', 'dummy_genes', 
        'colname_x', 'colname_y','colname_feature_name', 'colnames_count'
    }
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    cmds.append(f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}")
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

def convert_merscope(cmds, args):
    # input: in_csv
    # output: out_transcript, out_minmax, out_feature
    transcript_tsv = args.out_transcript.replace(".gz", "")
    format_cmd=f"python {cartloader}/scripts/format_merscope_custom.py --input {args.in_csv} --out-dir {args.out_dir} --out-transcript {transcript_tsv} --out-feature {args.out_feature} --out-minmax {args.out_minmax}"
    aux_argset = {
        'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name', 'csv_colname_feature_id', 'csv_colnames_others',
        'dummy_genes', 'precision_um',
        'colname_x', 'colname_y', 'colname_feature_name', 'colname_feature_id', 'colnames_count', 'colname_molecule_id'
    }
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    cmds.append(f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}")
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

def uniq_transcript(transcript_tsv, args):
    if args.unique:
        uniq_cmd = f"awk 'BEGIN {{ OFS=\"\\t\"; print \"{args.colname_x}\", \"{args.colname_y}\", \"{args.colname_feature_name}\", \"{args.annotation}\", \"{args.colnames_count}\" }} NR > 1 {{ if ($1 == prevX && $2 == prevY) {{ sumCount += $5; }} else {{ if (NR > 2) {{ print prevX, prevY, prevGene, firstCellID, sumCount; }} prevX = $1; prevY = $2; prevGene = $3; firstCellID = $4; sumCount = $5; }} }} END {{ print prevX, prevY, prevGene, firstCellID, sumCount; }}' {args.out_dir}/{transcript_tsv} | {args.gzip} -c > {args.out_dir}/{args.out_transcript}"
    else:
        uniq_cmd = f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}"
    return uniq_cmd

def convert_smi(cmds, args):
    # input: in_csv
    # output: out_transcript, out_minmax, out_feature
    transcript_tsv = args.out_transcript.replace(".gz", "")
    format_cmd=f"python {cartloader}/scripts/format_smi_custom.py --input {args.in_csv} --out-dir {args.out_dir} --out-transcript {transcript_tsv} --out-feature {args.out_feature} --out-minmax {args.out_minmax}"
    aux_argset = {
        'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name', 'csv_colnames_others',
        'units_per_um', 'dummy_genes', 'precision_um',
        'colname_x', 'colname_y','colname_feature_name', 'colnames_count'
    }
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    uniq_cmd = uniq_transcript(transcript_tsv, args)
    cmds.append(uniq_cmd)
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

def convert_stereoseq(cmds, args):
    # # cmds.append(f"(echo -e \"X\tY\tgene\tgn\"; {args.gzip} -cd {args.in_csv} | tail -n +2 | perl -lane 'print join(\"\\t\",$F[1]*0.5,$F[2]*0.5,$F[0],$F[3])' | {args.gzip} -c > {args.out_dir}/{args.out_transcript}")
    # input: in_csv
    # output: out_transcript, out_minmax, out_feature
    transcript_tsv = args.out_transcript.replace(".gz", "")
    format_cmd=f"python {cartloader}/scripts/format_stereoseq_custom.py --input {args.in_csv} --out-dir {args.out_dir} --out-transcript {transcript_tsv} --out-feature {args.out_feature} --out-minmax {args.out_minmax}"
    aux_argset = {
        'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name', 'csv_colnames_count',
        'units_per_um', 'dummy_genes', 'precision_um',
        'colname_x', 'colname_y','colname_feature_name', 'colnames_count'
    }
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    cmds.append(f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}")
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

def main(_args):
    # args
    args=parse_arguments(_args)
    scheck_app(args.gzip)
    scheck_app(args.sort)

    # input
    in_raw_filelist=convert_in_by_platform(args)
    for f in in_raw_filelist:
        assert f is not None, f"Missing input file for {f}. Please provide --{f}"
    # output
    os.makedirs(args.out_dir, exist_ok=True)
    # start cmds
    cmds = cmd_separator([], f"Converting input for raw data from: {args.platform}...")
    if args.platform == "10x_visium_hd":
        cmds = convert_visiumhd(cmds, args)
    elif args.platform == "10x_xenium":
        cmds = convert_xenium(cmds, args)
    elif args.platform == "cosmx_smi":
        cmds = convert_smi(cmds, args)      # return transcript, minmax
    elif args.platform == "bgi_stereoseq":
        cmds = convert_stereoseq(cmds, args)
    elif args.platform == "vizgen_merscope":
        cmds = convert_merscope(cmds, args)
    # mm
    mm = minimake()
    mm.add_target(args.out_transcript, in_raw_filelist, cmds) 
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    mm.write_makefile(f"{args.out_dir}/{args.makefn}")
    ## run makefile
    if args.dry_run:
        os.system(f"make -f {args.out_dir}/{args.makefn} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")
    else:
        os.system(f"make -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")

if __name__ == "__main__":
    global cartloader
    cartloader=os.path.dirname(os.path.abspath(__file__))
    # cartloader="/nfs/turbo/sph-hmkang/index/data/weiqiuc/cartloader/cartloader/"
    sys.path.extend(dir_i for dir_i in [cartloader] if dir_i not in sys.path)
    from utils import cmd_separator, scheck_app

    main(sys.argv[1:])
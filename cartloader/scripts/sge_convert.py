import sys, os, argparse, logging,  inspect, json, subprocess
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd
from cartloader.utils.sge_helper import aux_sge_args, input_by_platform, update_csvformat_by_platform
from cartloader.scripts.feature_filtering import filter_feature_by_type


# get the path of the cu

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""
                                     Standardize Spatial Transcriptomics (ST) datasets into a transcript-indexed SGE in TSV format.  
                                     Platform Supports: 10X Visium HD, SeqScope, 10X Xenium, BGI Stereoseq, Cosmx SMI, Vizgen Merscope, Pixel-Seq, and Nova-ST.
                                     Outputs: A transcript-indexed SGE file, a coordinate minmax TSV file, and a feature file counting UMIs per gene. All output are in micro-meter precision.
                                     Options: filtering SGE by quality, gene, or density; coordinate conversion.
                                     """)
    #parser = argparse.ArgumentParser()
    # run params
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_convert.mk", help='Makefile name (default: sge_convert.mk)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    
    # Input/output/key params
    inout_params = parser.add_argument_group("Input/Output Parameters", "Parameters to specify platform, input, output, and units per um, precision, and density-filtering for output.")
    inout_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic"], required=True, help='Platform of the raw input file to infer the format of the input file')
    # - input
    inout_params.add_argument('--in-mex', type=str, default=os.getcwd(), help='(10x_visium_hd and seqscope only) Directory path to input files in Market Exchange (MEX) format. Defaults to the current working directory.') # 10x_visium_hd, seqscope 
    inout_params.add_argument('--in-parquet', type=str, default="tissue_positions.parquet", help='(10x_visium_hd only) Path to the input parquet file for spatial coordinates (default: tissue_positions.parquet)') # 10x_visium_hd
    inout_params.add_argument('--in-csv', type=str, default=None, help='(10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st only) Path to the input raw CSV/TSV file (default: None).') # 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st
    inout_params.add_argument('--units-per-um', type=float, default=1.00, help='Coordinate unit per um in the input files (default: 1.00). Alternatively, for 10x Visium HD, skip --units-per-um and use --scale-json to auto-compute.')  
    inout_params.add_argument('--scale-json', type=str, default=None, help="(10x_visium_hd only) Path to a scale json file for calculating --units-per-um (default: None; Typical naming convention: scalefactors_json.json)") # 10x_visium_hd
    # - output
    inout_params.add_argument('--out-dir', type=str, required=True, help='The output directory to host files from SGE format conversion, files from density-filtering, and the make file.')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='Output for SGE format conversion. The compressed transcript-indexed SGE file in TSV format (default: transcripts.unsorted.tsv.gz).')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='Output for SGE format conversion. The coordinate minmax TSV file (default: coordinate_minmax.tsv).')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='Output for SGE format conversion. The compressed UMI count per gene TSV file (default: feature.clean.tsv.gz).')
    # - density filtering
    inout_params.add_argument('--filter-by-density', action='store_true', default=False, help='Filter SGE from format conversion by density (default: False). If enabled, check the density-filtering auxiliary parameters.')
    inout_params.add_argument('--out-filtered-prefix', type=str, default="filtered", help='Output for density-filtering. If --filter-by-density, define the prefix for filtered SGE (default: filtered)')
    # - sge visualization
    inout_params.add_argument('--sge-visual', action='store_true', default=False, help='Visualize the output SGE. If --filter-by-density, both unfiltered and filtered SGE will be visualized (default: False)')

    # AUX input MEX params
    aux_in_mex_params = parser.add_argument_group( "IN-MEX Auxiliary Parameters", "(10x_visium_hd and seqscope only) Auxiliary parameters for input MEX and parquet files. Required if --in-mex is used." )
    aux_in_mex_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='Barcode file name (default: barcodes.tsv.gz)')
    aux_in_mex_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='Feature file name (default: features.tsv.gz)')
    aux_in_mex_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='Matrix file name (default: matrix.mtx.gz)')
    aux_in_mex_params.add_argument('--icols-mtx', type=str, default=1, help='Comma-separated 1-based indices for the target genomic features among the count columns in the input matrix file (default: 1). For example, with SeqScope, "--icols-mtx 1,2,3,4,5" would include the first five count columns.')
    aux_in_mex_params.add_argument('--icol-ftr-id', type=int, default=1, help='1-based column index of feature ID in the input feature file (default: 1)')
    aux_in_mex_params.add_argument('--icol-ftr-name', type=int, default=2, help='1-based column index of feature name in the input feature file (default: 2)')
    aux_in_mex_params.add_argument('--icol-bcd-barcode', type=int, default=1, help='(seqscope only) 1-based column index of barcode in the input barcode file (default: 1)')
    aux_in_mex_params.add_argument('--icol-bcd-x', type=int, default=6, help='(seqscope only) 1-based column index of x coordinate in the input barcode file (default: 6)')
    aux_in_mex_params.add_argument('--icol-bcd-y', type=int, default=7, help='(seqscope only) 1-based column index of y coordinate in the input barcode file (default: 7)')
    aux_in_mex_params.add_argument('--pos-colname-barcode', type=str, default='barcode', help='(10x_visium_hd only) Column name for barcode in the input parquet file (default: barcode)')
    aux_in_mex_params.add_argument('--pos-colname-x', type=str, default='pxl_row_in_fullres', help='(10x_visium_hd only) Column name for X-axis in the input parquet file (default: pxl_row_in_fullres)')
    aux_in_mex_params.add_argument('--pos-colname-y', type=str, default='pxl_col_in_fullres', help='(10x_visium_hd only) Column name for Y-axis in the input parquet file (default: pxl_col_in_fullres)')
    aux_in_mex_params.add_argument('--pos-delim', type=str, default=',', help='(10x_visium_hd only) Delimiter for the input parquet file (default: ",")')
    aux_in_mex_params.add_argument('--print-feature-id', action='store_true', help='Print feature ID in the output (default: False)')
    aux_in_mex_params.add_argument('--allow-duplicate-gene-names', action='store_true', help='Allow duplicate gene names in the output (default: False)')
    aux_in_mex_params.add_argument('--keep-mismatches', action='store_true', help='(seqscope only) For debugging purposes, keep mismatches in the output (default: False)')

    # AUX input csv params
    aux_in_csv_params = parser.add_argument_group( "IN-CSV Auxiliary Parameters", "(10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st only) Auxiliary parameters for input TSV/CSV files")
    aux_in_csv_params.add_argument('--csv-comment', action='store_true', help="If comment in included in the csv file, specify the comment character (default: False for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, and pixel_seq; True for nova_st)")
    aux_in_csv_params.add_argument('--csv-delim', type=str, default=None, help='Delimiter for the additional input tsv/csv file (default: "," for 10x_xenium, cosmx_smi, and vizgen_merscope; "\\t" for bgi_stereoseq, pixel_seq, and nova_st) ')
    aux_in_csv_params.add_argument('--csv-colname-x',  type=str, default=None, help='Column name for X-axis (default: x_location for 10x_xenium; x for bgi_stereoseq; x_local_px for cosmx_smi; global_x for vizgen_merscope; xcoord for pixel_seq; x for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-y',  type=str, default=None, help='Column name for Y-axis (default: y_location for 10x_xenium; y for bgi_stereoseq; y_local_px for cosmx_smi; global_y for vizgen_merscope; ycoord for pixel_seq; y for nova_st)')
    aux_in_csv_params.add_argument('--csv-colnames-count', type=str, default=None, help='Comma-separated Column name for expression count. If not provided, a count of 1 will be added for a feature in a pixel (default: MIDCounts for bgi_stereoseq; MIDCount for nova_st; None for the rest platforms).')
    aux_in_csv_params.add_argument('--csv-colname-feature-name', type=str, default=None, help='Column name for gene name (default: feature_name for 10x_xenium; geneID for bgi_stereoseq; target for cosmx_smi; gene for vizgen_merscope; geneName for pixel_seq; geneID for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-feature-id', type=str, default=None, help='Column name for gene id (default: None)')
    aux_in_csv_params.add_argument('--csv-colnames-others', nargs='*', default=[], help='Columns names to keep (e.g., cell_id, overlaps_nucleus) (default: None)')
    aux_in_csv_params.add_argument('--csv-colname-phredscore', type=str, default=None, help='Column name for Phred-scaled quality value (Q-Score) estimating the probability of incorrect call (default: qv for 10x_xenium and None for the rest platforms).') # qv
    aux_in_csv_params.add_argument('--min-phred-score', type=float, default=None, help='Specify the Phred-scaled quality score cutoff (default: 20 for 10x_xenium and None for the rest platforms).')
    aux_in_csv_params.add_argument('--add-molecule-id', action='store_true', default=False, help='If enabled, a column of "molecule_id" will be added to the output file to track the index of the original input will be stored in (default: False).')
    
    # AUX output params
    aux_out_params = parser.add_argument_group("Output Auxiliary Parameters", "Auxiliary parameters for the output files (Recommand to use the default values)")
    aux_out_params.add_argument('--precision-um', type=int, default=2, help='Precision for transcript coordinates. Set it to 0 to round to integer (default: 2)')
    aux_out_params.add_argument('--colname-x', type=str, default='X', help='Column name for X (default: X)')
    aux_out_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y (default: Y)')
    aux_out_params.add_argument('--colnames-count', type=str, default='count', help='Comma-separated column names for count (default: count)')
    aux_out_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for gene name (default: gene)')
    aux_out_params.add_argument('--colname-feature-id', type=str, default=None, help='Column name for gene ID. Required only when --csv-colname-feature-id or --print-feature-id is applied (default: None)') 

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Filtering Auxiliary Parameters", "Auxiliary parameters for filtering features by their name or type using an additional file, substring, or regex pattern")
    aux_ftrfilter_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-substr', type=str, default=None, help='A substring of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-substr', type=str, default=None, help='A substring of feature/gene names to be excluded (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None).') # (e.g. protein_coding|lncRNA)
    aux_ftrfilter_params.add_argument('--csv-colname-feature-type', type=str, default=None, help='The input column name in the input that corresponding to the gene type information, if your input file has gene type information(default: None)')
    aux_ftrfilter_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated reference file to provide gene type information for each each per row (default: None)')
    aux_ftrfilter_params.add_argument('--feature-type-ref-delim', type=str, default=None, help='Delimiter used in the reference file (default: tab).')
    aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='Column index for gene name in the reference file (default: None).')
    aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='Column index for gene type in the reference file (default: None).')
    aux_ftrfilter_params.add_argument('--print-removed-transcripts', action='store_true', default=False, help='Print the list of removed transcript with corresponding filtering criteria (default: False)')

    # AUX polygon-filtering params
    aux_polyfilter_params = parser.add_argument_group('Density/Polygon Filtering Auxiliary Parameters','Auxiliary parameters for filtering polygons based on the number of vertices. Required when --filter-by-density is enabled.')
    aux_polyfilter_params.add_argument('--genomic-feature', type=str, default=None, help='Column name of genomic feature for polygon-filtering (default to --colnames-count if --colnames-count only specifies one column)')
    aux_polyfilter_params.add_argument('--mu-scale', type=int, default=1, help='Scale factor for the polygon area calculation (default: 1.0)')
    aux_polyfilter_params.add_argument('--radius', type=int, default=15, help='Radius for the polygon area calculation (default: 15)')
    aux_polyfilter_params.add_argument('--quartile', type=int, default=2, help='Quartile for the polygon area calculation (default: 2)')
    aux_polyfilter_params.add_argument('--hex-n-move', type=int, default=1, help='Sliding step (default: 1)')
    aux_polyfilter_params.add_argument('--polygon-min-size', type=int, default=500, help='The minimum polygon size (default: 500)')
    # gene_header and count_header will be automatically based on the --colname-feature-name, --colname-feature-id and --colnames-count

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4".')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary. ')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools binary')
    
    # not in use 
    #aux_ftrfilter_params.add_argument('--unique', action='store_true', default=False, help='Merge pixels with (almost?) identical coordinates. Applies to cosmx_smi only.')
 
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

#================================================================================================
#
# in-mex general functions
#
#================================================================================================

def extract_unit2px_from_json(scale_json):
    # purpose: Extract the microns per pixel value from the scale json file and calculate the units per um.
    print(f"As --scale-json is provided, calculating units per um based on the microns per pixel value from {scale_json}")
    assert os.path.exists(scale_json), f"The scale json file ({scale_json}) does not exist. Please provide the correct path using --scale-json."
    with open(scale_json, 'r') as file:
        scale_data = json.load(file)
    # Extract the value for 'microns_per_pixel'
    microns_per_pixel = scale_data['microns_per_pixel']
    # microns_per_pixel cannot be NA or zero
    assert microns_per_pixel is not None, "The value for 'microns_per_pixel' is not found in the json file."
    assert microns_per_pixel != 0, "The value for 'microns_per_pixel' is 0. Please check the json file."
    print(f"    - Microns per pixel: {microns_per_pixel}")
    return 1/microns_per_pixel

mexarg_mapping = {
    "mex_bcd": "--sge-bcd",
    "mex_ftr": "--sge-ftr",
    "mex_mtx": "--sge-mtx"
}

def add_mexparam_to_cmd(format_cmd, args, arg_mapping):
    'Add mex file names to the command'
    for arg_name, cmd_option in arg_mapping.items():
        arg_value = getattr(args, arg_name, None)
        if arg_value is not None:
            format_cmd = f"{format_cmd} {cmd_option} {arg_value}"
    return format_cmd

def convert_visiumhd(cmds, args):
    # tools:
    # 1) spatula
    # if args.spatula is None:
    #     args.spatula = os.path.join(cartloader_repo, "submodules", "spatula", "bin", "spatula")
    scheck_app(args.spatula)
    # input: in_mex, in_parquet, scale_json
    # output: out_transcript, out_minmax, out_feature
    tmp_parquet = f"{args.out_dir}/tissue_positions.csv.gz"
    # * --in_parquet: convert parquet to csv
    cmds.append(f"{args.parquet_tools} csv {args.in_parquet} |  {args.gzip} -c > {tmp_parquet}")
    # * --scale_json: if applicable
    if args.scale_json is not None:
        args.units_per_um = extract_unit2px_from_json(args.scale_json)
    # * --included_feature_type_regex
    if args.include_feature_type_regex is not None:
        # purpose: Given the spatula doesn't support filtering based on gene type, this function prepares a feature list based on the gene type reference file.
        #           Note, when the --include-feature-type-regex is enabled, the --include-feature-list argument will be updated with the path to the prepared feature list.
        # since MEX will not be able to provide ftr type column.
        ftrlist_tsv = f"{args.out_dir}/feature_list.tsv"
        cmd = " ".join(["cartloader", "feature_filtering_by_type",
                        f"--include-feature-type-regex {args.include_feature_type_regex}", 
                        f"--feature-type-ref {args.feature_type_ref}",
                        f"--feature-type-ref-colidx-name {args.feature_type_ref_colidx_name}",
                        f"--feature-type-ref-colidx-type {args.feature_type_ref_colidx_type}",
                        f"--feature-type-ref-delim {args.feature_type_ref_delim}" if args.feature_type_ref_delim is not None else "",
                        f"--include-feature-list {args.include_feature_list}"  if args.include_feature_list is not None else "",
                        f"--exclude-feature-list {args.exclude_feature_list}"  if args.include_feature_list is not None else "",
                        f"--out-feature {ftrlist_tsv}",
                        "--log"])
        args.include_feature_list = ftrlist_tsv
        args.exclude_feature_list = None
        cmds.append(cmd)
    # * check --icols-mtx has the same number as --colnames-count
    args.icols_mtx=str(args.icols_mtx)
    if len(args.icols_mtx.split(",")) != len(args.colnames_count.split(",")):
        raise ValueError(f"The number of columns in --icols-mtx ({args.icols_mtx}) should be the same as the number of columns in --colnames-count ({args.colnames_count}).")
    # * convert sge to tsv (output: out_transcript, out_minmax, out_feature, (optional) out_sge)
    format_cmd=f"{args.spatula} convert-sge --in-sge {args.in_mex} --out-tsv {args.out_dir} --pos {tmp_parquet} --tsv-mtx {args.out_transcript} --tsv-ftr {args.out_feature} --tsv-minmax {args.out_minmax}"
    aux_argset = set(item for lst in [aux_sge_args["out"], aux_sge_args["inftr"], aux_sge_args["inmtx"], aux_sge_args["inpos"], aux_sge_args["spatula"], aux_sge_args["ftrname"]] for item in lst)
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    format_cmd = add_mexparam_to_cmd(format_cmd, args, mexarg_mapping)
    cmds.append(format_cmd)
    cmds.append(f"rm {tmp_parquet}")
    return cmds

def convert_seqscope(cmds, args):
    # tools:
    # 1) spatula
    # if args.spatula is None:
    #     args.spatula = os.path.join(cartloader_repo, "submodules", "spatula", "bin", "spatula")
    scheck_app(args.spatula)
    # input: in_mex
    # output: out_transcript, out_minmax, out_feature
    # * --included_feature_type_regex
    if args.include_feature_type_regex is not None:
        # purpose: Given the spatula doesn't support filtering based on gene type, this function prepares a feature list based on the gene type reference file.
        #           Note, when the --include-feature-type-regex is enabled, the --include-feature-list argument will be updated with the path to the prepared feature list.
        # since MEX will not be able to provide ftr type column.
        ftrlist_tsv = f"{args.out_dir}/feature_list.tsv"
        cmd = " ".join(["cartloader", "feature_filtering_by_type",
                        f"--include-feature-type-regex {args.include_feature_type_regex}", 
                        f"--feature-type-ref {args.feature_type_ref}",
                        f"--feature-type-ref-colidx-name {args.feature_type_ref_colidx_name}",
                        f"--feature-type-ref-colidx-type {args.feature_type_ref_colidx_type}",
                        f"--feature-type-ref-delim {args.feature_type_ref_delim}" if args.feature_type_ref_delim is not None else "",
                        f"--include-feature-list {args.include_feature_list}"  if args.include_feature_list is not None else "",
                        f"--exclude-feature-list {args.exclude_feature_list}"  if args.include_feature_list is not None else "",
                        f"--out-feature {ftrlist_tsv}",
                        "--log"])
        args.include_feature_list = ftrlist_tsv
        args.exclude_feature_list = None
        cmds.append(cmd)
    # * check --icols-mtx has the same number as --colnames-count
    args.icols_mtx=str(args.icols_mtx)
    if len(args.icols_mtx.split(",")) != len(args.colnames_count.split(",")):
        raise ValueError(f"The number of columns in --icols-mtx ({args.icols_mtx}) should be the same as the number of columns in --colnames-count ({args.colnames_count}).")
    # * convert sge to tsv (output: out_transcript, out_minmax, out_feature
    format_cmd=f"{args.spatula} convert-sge --in-sge {args.in_mex} --out-tsv {args.out_dir} --tsv-mtx {args.out_transcript} --tsv-ftr {args.out_feature} --tsv-minmax {args.out_minmax}"
    aux_argset = set(item for lst in [aux_sge_args["out"], aux_sge_args["inftr"], aux_sge_args["inbcd"], aux_sge_args["inmtx"], aux_sge_args["ftrname"]] for item in lst)
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    format_cmd = add_mexparam_to_cmd(format_cmd, args, mexarg_mapping)
    cmds.append(format_cmd)
    # * drop mismatches if --keep-mismatches is not enabled
    if not args.keep_mismatches:   
        drop_cmd = f"cartloader sge_drop_mismatches --in-dir {args.out_dir} --transcript {args.out_transcript} --feature {args.out_feature} --minmax {args.out_minmax} --gzip {args.gzip}"    
        cmds.append(drop_cmd)
    return cmds

#================================================================================================
#
# in-tsv/csv general functions
#
#================================================================================================
def convert_tsv(cmds, args):
    # input: in_csv
    # output: out_transcript, out_minmax, out_feature
    #  * update csv_colname_*  and csv_delim based on the platform
    args = update_csvformat_by_platform(args)
    #print(args)
    #  * 10x_xenium: update default value for the phred score filtering 
    if args.platform == "10x_xenium":
        if args.csv_colname_phredscore is None:
            args.csv_colname_phredscore = "qv"
        if args.min_phred_score is None:
            args.min_phred_score = 20
    # * no need to check --csv-colnames-count has the same number as --colnames-count (for spme platform, no need to provide --csv-colnames-count)
    # if len(args.csv_colnames_count.split(",")) != len(args.colnames_count.split(",")):
    #     raise ValueError(f"The number of columns in --csv-colnames-count ({args.csv_colnames_count}) should be the same as the number of columns in --colnames-count ({args.colnames_count}).")
    # main commands
    transcript_tsv = args.out_transcript.replace(".gz", "")
    format_cmd=f"cartloader format_generic --input {args.in_csv} --out-dir {args.out_dir} --out-transcript {transcript_tsv} --out-feature {args.out_feature} --out-minmax {args.out_minmax}"
    # aux args
    aux_argset = set(item for lst in [aux_sge_args["out"], aux_sge_args["incsv"], aux_sge_args["ftrname"], aux_sge_args["ftrtype"]] for item in lst)
    aux_argset.add('print_removed_transcripts')
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    # append to cmds
    cmds.append(format_cmd)
    cmds.append(f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}")
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

#================================================================================================
#
# density filtering functions
#
#================================================================================================

def sge_density_filtering(mm, args):
    gene_header=[args.colname_feature_name] if args.colname_feature_id is None else [args.colname_feature_name, args.colname_feature_id]
    count_header=args.colnames_count.split(",")

    if args.genomic_feature is None:
        if len(args.colnames_count.split(",")) > 1:
            logging.error("Missing --genomic-feature. Cannot use --colnames-count with multiple columns as --genomic-feature. Please provide one column name for density-filtering.")
            sys.exit(1)
        args.genomic_feature = args.colnames_count
        
    sge_convert_flag = os.path.join(args.out_dir, "sge_convert.done")
    filtered_transcript_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.transcripts.unsorted.tsv.gz")
    cmds=cmd_separator([], "Filtering the converted data by density...")
    cmd = " ".join([
            "ficture", "filter_by_density",
            f"--input {args.out_dir}/{args.out_transcript}",
            f"--feature {args.out_dir}/{args.out_feature}",
            f"--output {filtered_transcript_f}",
            f"--output_boundary {args.out_dir}/{args.out_filtered_prefix}",
            f"--filter_based_on {args.genomic_feature}",
            f"--mu_scale {args.mu_scale}",
            f"--radius {args.radius}",
            f"--quartile {args.quartile}",
            f"--hex_n_move {args.hex_n_move}",
            f"--remove_small_polygons {args.polygon_min_size}",
            f"--gene_header {' '.join(gene_header)}",
            f"--count_header {' '.join(count_header)}",
        ])
    cmds.append(cmd)
    mm.add_target(f"{filtered_transcript_f}", [sge_convert_flag], cmds)
    return mm

#================================================================================================
#
# main functions
#
#================================================================================================
def sge_visual(mm, args, transcript_f, xy_f, prereq):
    scheck_app(args.spatula)
    # draw xy plot for visualization
    cmds = cmd_separator([], f"Drawing XY plot for SGE: {transcript_f}")
    draw_cmd=f"{args.gzip} -dc {transcript_f} | tail -n +2 | cut -f 1,2 | {args.spatula} draw-xy --tsv /dev/stdin --out {xy_f}"
    cmds.append(draw_cmd)
    mm.add_target(xy_f, prereq, cmds)
    return mm

def sge_convert(_args):
    # args
    args=parse_arguments(_args)
    if args.exclude_feature_regex == "":
        args.exclude_feature_regex = None
    scheck_app(args.gzip)

    # input
    in_raw_filelist=input_by_platform(args)

    # output
    os.makedirs(args.out_dir, exist_ok=True)

    out_transcript_f = os.path.join(args.out_dir, args.out_transcript)
    out_xy_f = os.path.join(args.out_dir, "xy.png")

    filtered_transcript_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.transcripts.unsorted.tsv.gz")
    filtered_xy_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.xy.png")

    # mm
    mm = minimake()

    # sge_convert
    sge_convert_flag = os.path.join(args.out_dir, "sge_convert.done")

    cmds = cmd_separator([], f"Converting input for raw data from: {args.platform}...")
    if args.platform == "10x_visium_hd":
        cmds = convert_visiumhd(cmds, args)
    elif args.platform == "seqscope":
        cmds = convert_seqscope(cmds, args)
    elif args.platform in ["10x_xenium", "cosmx_smi", "bgi_stereoseq", "vizgen_merscope", "pixel_seq", "nova_st", "generic"]:
        cmds = convert_tsv(cmds, args)
    cmds.append(f"[ -f {out_transcript_f} ] && [ -f {os.path.join(args.out_dir, args.out_feature)} ] && [ -f {os.path.join(args.out_dir, args.out_minmax)} ] && touch {sge_convert_flag}")
    mm.add_target(sge_convert_flag, in_raw_filelist, cmds) 

    if args.sge_visual:
        mm = sge_visual(mm, args, 
                        out_transcript_f,
                        out_xy_f,
                        [sge_convert_flag])        
            
   
    if args.filter_by_density:
        mm = sge_density_filtering(mm, args)
    
    if args.filter_by_density and args.sge_visual:
        mm = sge_visual(mm, args, 
                        filtered_transcript_f,
                        filtered_xy_f,
                        [filtered_transcript_f])

    # write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)
    
    if args.dry_run:
        dry_cmd=f"make -f {make_f} -n {'-B' if args.restart else ''} "
        os.system(dry_cmd)
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}")
    else:
        exe_cmd=f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}"
        result = subprocess.run(exe_cmd, shell=True)
        if result.returncode != 0:
            print(f"Error in executing: {exe_cmd}")
            sys.exit(1)

if __name__ == "__main__":
    # get the cartloader path
    global cartloader_repo
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])


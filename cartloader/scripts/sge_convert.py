import sys, os, argparse, logging,  inspect, json, subprocess
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, write_dict_to_file, load_file_to_dict, execute_makefile
from cartloader.utils.sge_helper import aux_sge_args, input_by_platform, update_csvformat_by_platform
from cartloader.scripts.feature_filtering import filter_feature_by_type

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""
                                     Standardize Format of Spatial Transcriptomics (ST) datasets. Each returns a transcript-indexed SGE in TSV format.  
                                     Platform Supports: 10X Visium HD, SeqScope, 10X Xenium, BGI Stereoseq, Cosmx SMI, Vizgen Merscope, Pixel-Seq, and Nova-ST. For SGE from others platforms or from custom/preprocessed sources, sge_convert provides a generic option that accepts CSV/TSV files with basic required fields
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
    inout_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic"], required=True, help='Platform of the raw input file to infer the format of the input file. "generic" refers to SGE from platforms not yet explicitly supported by cartloader, or from custom/preprocessed sources')
    # - input
    inout_params.add_argument('--in-json', type=str, default=None, help='(Shortcut; currenly only support 10x_xenium, 10x_visium_hd) Path to a JSON file to provide paths to the input file. If provided, omit --in-parquet/--in-csv/--scale-json.')
    inout_params.add_argument('--in-mex', type=str, default=os.getcwd(), help='(10x_visium_hd and seqscope only) Directory path to input files in Market Exchange (MEX) format. Defaults to the current working directory.') # 10x_visium_hd, seqscope 
    inout_params.add_argument('--in-parquet', type=str, default=None, help='(10x_visium_hd and 10x_xenium only) For 10X Visium HD platform, specify to path to the input parquet file for spatial coordinates (default: tissue_positions.parquet). For 10X Xenium, if the input transcript file is in parquet format, specify its path here and skip --in-csv (default: None)') # 10x_visium_hd
    inout_params.add_argument('--in-csv', type=str, default=None, help='(10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st only) Path to the input raw CSV/TSV file if raw CSV/TSV file exists(default: None).') 
    inout_params.add_argument('--units-per-um', type=float, default=1.00, help='Coordinate unit per um in the input files (default: 1.00). Alternatively, for 10x Visium HD, skip --units-per-um and use --scale-json to auto-compute.')  
    inout_params.add_argument('--scale-json', type=str, default=None, help="(Shortcut; 10x_visium_hd only) Path to a scale json file for calculating --units-per-um (default: None; Typical naming convention: scalefactors_json.json)") 
    # - output
    inout_params.add_argument('--out-dir', type=str, required=True, help='The output directory to host files from SGE format conversion, files from density-filtering, and the make file.')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='Output for SGE format conversion. The compressed transcript-indexed SGE file in TSV format (default: transcripts.unsorted.tsv.gz).')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='Output for SGE format conversion. The coordinate minmax TSV file (default: coordinate_minmax.tsv).')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='Output for SGE format conversion. The compressed UMI count per gene TSV file (default: feature.clean.tsv.gz).')
    # - density filtering
    inout_params.add_argument('--filter-by-density', action='store_true', default=False, help='Filter SGE from format conversion by density (default: False). If enabled, check the density-filtering auxiliary parameters.')
    inout_params.add_argument('--out-filtered-prefix', type=str, default="filtered", help='Output for density-filtering. If --filter-by-density, define the prefix for filtered SGE and its visualization (default: filtered)')
    # - sge visualization
    inout_params.add_argument('--sge-visual', action='store_true', default=False, help='Visualize the output SGE. If --filter-by-density, both unfiltered and filtered SGE will be visualized (default: False)')
    inout_params.add_argument('--out-xy', type=str, default="xy.png", help='Output for SGE visualization image (default: xy.png)')
    # - sge json
    inout_params.add_argument('--out-json',  type=str, default=None, help='Output json summarizing SGE information. (default: out_dir/sge_assets.json).')

    # AUX input MEX params
    aux_in_mex_params = parser.add_argument_group( "IN-MEX Auxiliary Parameters", "(10x_visium_hd and seqscope only) Auxiliary parameters for input MEX and parquet files. Required if --in-mex is used." )
    aux_in_mex_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='Barcode file name (default: barcodes.tsv.gz)')
    aux_in_mex_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='Feature file name (default: features.tsv.gz)')
    aux_in_mex_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='Matrix file name (default: matrix.mtx.gz)')
    aux_in_mex_params.add_argument('--icols-mtx', type=str, default=1, help='1-based indices for the target genomic feature among the count columns in the input matrix file (default: 1). For example, when focusing on "gn" of SeqScope datasets, use --icols-mtx 1.')
    aux_in_mex_params.add_argument('--icol-ftr-id', type=int, default=1, help='1-based column index of feature ID in the input feature file (default: 1)')
    aux_in_mex_params.add_argument('--icol-ftr-name', type=int, default=2, help='1-based column index of feature name in the input feature file (default: 2)')
    aux_in_mex_params.add_argument('--icol-bcd-barcode', type=int, default=1, help='(seqscope only) 1-based column index of barcode in the input barcode file (default: 1)')
    aux_in_mex_params.add_argument('--icol-bcd-x', type=int, default=6, help='(seqscope only) 1-based column index of x coordinate in the input barcode file (default: 6)')
    aux_in_mex_params.add_argument('--icol-bcd-y', type=int, default=7, help='(seqscope only) 1-based column index of y coordinate in the input barcode file (default: 7)')
    aux_in_mex_params.add_argument('--pos-colname-barcode', type=str, default='barcode', help='(10x_visium_hd only) Column name for barcode in the input parquet file (default: barcode)')
    aux_in_mex_params.add_argument('--pos-colname-x', type=str, default='pxl_row_in_fullres', help='(10x_visium_hd only) Column name for X-axis in the input parquet file (default: pxl_row_in_fullres)')
    aux_in_mex_params.add_argument('--pos-colname-y', type=str, default='pxl_col_in_fullres', help='(10x_visium_hd only) Column name for Y-axis in the input parquet file (default: pxl_col_in_fullres)')
    aux_in_mex_params.add_argument('--pos-delim', type=str, default=',', help='(10x_visium_hd only) Delimiter for the input parquet file (default: ",")')
    #aux_in_mex_params.add_argument('--print-feature-id', action='store_true', help='Print feature ID in the output (default: False)')
    aux_in_mex_params.add_argument('--allow-duplicate-gene-names', action='store_true', help='Allow duplicate gene names in the output (default: False)')
    aux_in_mex_params.add_argument('--keep-mismatches', action='store_true', help='(seqscope only) For debugging purposes, keep mismatches in the output (default: False)')

    # AUX input csv params
    aux_in_csv_params = parser.add_argument_group( "IN-CSV Auxiliary Parameters", "(10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st only) Auxiliary parameters for input TSV/CSV files")
    aux_in_csv_params.add_argument('--csv-comment', action='store_true', help="If enabled, lines starts with # in the CSV file will be skipped (default: False for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, and pixel_seq; True for nova_st)")
    aux_in_csv_params.add_argument('--csv-delim', type=str, default=None, help='Delimiter for the input TSV/CSV file (default: "," for 10x_xenium, cosmx_smi, and vizgen_merscope; "\\t" for bgi_stereoseq, pixel_seq, and nova_st) ')
    aux_in_csv_params.add_argument('--csv-colname-x',  type=str, default=None, help='Column name for X-axis (default: x_location for 10x_xenium; x for bgi_stereoseq; x_local_px for cosmx_smi; global_x for vizgen_merscope; xcoord for pixel_seq; x for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-y',  type=str, default=None, help='Column name for Y-axis (default: y_location for 10x_xenium; y for bgi_stereoseq; y_local_px for cosmx_smi; global_y for vizgen_merscope; ycoord for pixel_seq; y for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-count', type=str, default=None, help='Comma-separated Column name for expression count. If not provided, a count of 1 will be added for a feature in a pixel (default: MIDCounts for bgi_stereoseq; MIDCount for nova_st; None for the rest platforms).')
    aux_in_csv_params.add_argument('--csv-colname-feature-name', type=str, default=None, help='Column name for gene name (default: feature_name for 10x_xenium; geneID for bgi_stereoseq; target for cosmx_smi; gene for vizgen_merscope; geneName for pixel_seq; geneID for nova_st)')
    # aux_in_csv_params.add_argument('--csv-colname-feature-id', type=str, default=None, help='Column name for gene id (default: None)')
    aux_in_csv_params.add_argument('--csv-colnames-others', nargs='*', default=[], help='Columns names to keep (e.g., cell_id, overlaps_nucleus) (default: None)')
    aux_in_csv_params.add_argument('--csv-colname-phredscore', type=str, default=None, help='Column name for Phred-scaled quality value (Q-Score) estimating the probability of incorrect call (default: qv for 10x_xenium and None for the rest platforms).') # qv
    aux_in_csv_params.add_argument('--min-phred-score', type=float, default=None, help='Specify the Phred-scaled quality score cutoff (default: 20 for 10x_xenium and None for the rest platforms).')
    #aux_in_csv_params.add_argument('--add-molecule-id', action='store_true', default=False, help='If enabled, a column of "molecule_id" will be added to the output file to track the index of the original input will be stored in (default: False).')
    
    # AUX output params
    aux_out_params = parser.add_argument_group("Output Auxiliary Parameters", "Auxiliary parameters for the output files (Recommand to use the default values)")
    aux_out_params.add_argument('--precision-um', type=int, default=2, help='Precision for transcript coordinates. Set it to 0 to round to integer (default: 2)')
    aux_out_params.add_argument('--colname-x', type=str, default='X', help='Column name for X (default: X)')
    aux_out_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y (default: Y)')
    aux_out_params.add_argument('--colname-count', type=str, default='count', help='Comma-separated column names for count (default: count)')
    aux_out_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for gene name (default: gene)')
    # aux_out_params.add_argument('--colname-feature-id', type=str, default=None, help='Column name for gene ID. Required only when --csv-colname-feature-id or --print-feature-id is applied (default: None)') 

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group( "Feature Filtering Auxiliary Parameters", "Auxiliary parameters for filtering features by their name or type using an additional file, or regex pattern")
    # aux_ftrfilter_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    # aux_ftrfilter_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: "^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)" for 10_xenium, None for the rest)')
    # aux_ftrfilter_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None).') # (e.g. protein_coding|lncRNA)
    # aux_ftrfilter_params.add_argument('--csv-colname-feature-type', type=str, default=None, help='The input column name in the input that corresponding to the gene type information, if your input file has gene type information(default: None)')
    # aux_ftrfilter_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated reference file to provide gene type information for each gene per row (default: None)')
    # aux_ftrfilter_params.add_argument('--feature-type-ref-delim', type=str, default=None, help='Delimiter used in the reference file (default: tab).')
    # aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-name', type=str, default=None, help='Column index for gene name in the reference file (default: None).')
    # aux_ftrfilter_params.add_argument('--feature-type-ref-colidx-type', type=str, default=None, help='Column index for gene type in the reference file (default: None).')
    # aux_ftrfilter_params.add_argument('--print-removed-transcripts', action='store_true', default=False, help='Print the list of removed transcript with corresponding filtering criteria (default: False)')

    # AUX polygon-filtering params
    aux_polyfilter_params = parser.add_argument_group('Density/Polygon Filtering Auxiliary Parameters','Auxiliary parameters for filtering polygons based on the number of vertices. Required when --filter-by-density is enabled.')
    # aux_polyfilter_params.add_argument('--genomic-feature', type=str, default=None, help='Column name of genomic feature for polygon-filtering (default to --colname-count)')
    aux_polyfilter_params.add_argument('--radius', type=int, default=15, help='Advanced parameter. Radius for the polygon area calculation (default: 15)')
    aux_polyfilter_params.add_argument('--quartile', type=int, default=2, help='Quartile for the polygon area calculation (default: 2)')
    aux_polyfilter_params.add_argument('--hex-n-move', type=int, default=1, help='Sliding step (default: 1)')
    aux_polyfilter_params.add_argument('--polygon-min-size', type=int, default=500, help='The minimum polygon size (default: 500)')
    # gene_header and count_header will be automatically based on the --colname-feature-name, --colname-feature-id and --colnames-count

    # AUX visualization params
    aux_visual_params = parser.add_argument_group("North-up Auxiliary Parameters", "Auxiliary parameters for visualizing the output SGE in a north-up orientation. It is optional if --sge-visual is enabled.")
    aux_visual_params.add_argument('--north-up', action='store_true', default=False, help='If enabled, the SGE will be visualized in a tif image north-up (default: False).')
    aux_visual_params.add_argument('--out-northup-tif', type=str, default="xy_northup.tif", help='Output for SGE visualization after north-up orientation (default: xy_northup.tif).')
    aux_visual_params.add_argument('--srs', type=str, default='EPSG:3857', help='If --north-up, define the spatial reference system (default: EPSG:3857)')
    aux_visual_params.add_argument('--resample', type=str, default='cubic', help='If --north-up, define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4" (default: gzip)')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary (default: spatula).')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='If --in-parquet is enabled, path to parquet-tools binary (default: parquet-tools)')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='If --north-up, provide path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdalwarp', type=str, default=f"gdalwarp", help='If --north-up, provide path to gdalwarp binary (default: gdalwarp)')

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
    ## tools:
    scheck_app(args.spatula)
    ## input: in_mex, in_parquet, scale_json
    ## output: out_transcript, out_minmax, out_feature
    tmp_parquet = f"{args.out_dir}/tissue_positions.csv.gz"
    # * --in_parquet: convert parquet to csv
    if args.parquet is None:
        args.parquet = "tissue_positions.parquet"
    cmds.append(f"{args.parquet_tools} csv {args.in_parquet} |  {args.gzip} -c > {tmp_parquet}")
    # * --scale_json: if applicable
    if args.scale_json is not None:
        args.units_per_um = extract_unit2px_from_json(args.scale_json)
    # * --included_feature_type_regex
    # if args.include_feature_type_regex is not None:
    #     # purpose: Given the spatula doesn't support filtering based on gene type, this function prepares a feature list based on the gene type reference file.
    #     #           Note, when the --include-feature-type-regex is enabled, the --include-feature-list argument will be updated with the path to the prepared feature list.
    #     # since MEX will not be able to provide ftr type column.
    #     ftrlist_tsv = f"{args.out_dir}/feature_list.tsv"
    #     cmd = " ".join(["cartloader", "feature_filtering_by_type",
    #                     f"--include-feature-type-regex {args.include_feature_type_regex}", 
    #                     f"--feature-type-ref {args.feature_type_ref}",
    #                     f"--feature-type-ref-colidx-name {args.feature_type_ref_colidx_name}",
    #                     f"--feature-type-ref-colidx-type {args.feature_type_ref_colidx_type}",
    #                     f"--feature-type-ref-delim {args.feature_type_ref_delim}" if args.feature_type_ref_delim is not None else "",
    #                     f"--include-feature-list {args.include_feature_list}"  if args.include_feature_list is not None else "",
    #                     f"--exclude-feature-list {args.exclude_feature_list}"  if args.include_feature_list is not None else "",
    #                     f"--out-feature {ftrlist_tsv}",
    #                     "--log"])
    #     args.include_feature_list = ftrlist_tsv
    #     args.exclude_feature_list = None
    #     cmds.append(cmd)
    # * check --icols-mtx has the same number as --colnames-count
    # args.icols_mtx=str(args.icols_mtx)
    # if len(args.icols_mtx.split(",")) != len(args.colnames_count.split(",")):
    #     raise ValueError(f"The number of columns in --icols-mtx ({args.icols_mtx}) should be the same as the number of columns in --colnames-count ({args.colnames_count}).")
    # * convert sge to tsv (output: out_transcript, out_minmax, out_feature, (optional) out_sge)
    cmd = " ".join([f"{args.spatula} convert-sge",
                    f"--in-sge {args.in_mex}",
                    f"--out-tsv {args.out_dir}",
                    f"--pos {tmp_parquet}",
                    f"--tsv-mtx {args.out_transcript}",
                    f"--tsv-ftr {args.out_feature}",
                    f"--tsv-minmax {args.out_minmax}",
                    f"--colnames-count {args.colname_count}" if args.colname_count else ""])
    aux_argset = set(item for lst in [aux_sge_args["out"], aux_sge_args["inftr"], aux_sge_args["inmtx"], aux_sge_args["inpos"], aux_sge_args["spatula"], aux_sge_args["ftrname"]] for item in lst)
    cmd = add_param_to_cmd(cmd, args, aux_argset)
    cmd = add_mexparam_to_cmd(cmd, args, mexarg_mapping)
    cmds.append(cmd)
    cmds.append(f"rm {tmp_parquet}")
    return cmds

def convert_seqscope(cmds, args):
    ## tools:
    scheck_app(args.spatula)
    ## input: in_mex
    ## output: out_transcript, out_minmax, out_feature
    # * --included_feature_type_regex
    # if args.include_feature_type_regex is not None:
    #     # purpose: Given the spatula doesn't support filtering based on gene type, this function prepares a feature list based on the gene type reference file.
    #     #           Note, when the --include-feature-type-regex is enabled, the --include-feature-list argument will be updated with the path to the prepared feature list.
    #     # since MEX will not be able to provide ftr type column.
    #     ftrlist_tsv = f"{args.out_dir}/feature_list.tsv"
    #     cmd = " ".join(["cartloader", "feature_filtering_by_type",
    #                     f"--include-feature-type-regex {args.include_feature_type_regex}", 
    #                     f"--feature-type-ref {args.feature_type_ref}",
    #                     f"--feature-type-ref-colidx-name {args.feature_type_ref_colidx_name}",
    #                     f"--feature-type-ref-colidx-type {args.feature_type_ref_colidx_type}",
    #                     f"--feature-type-ref-delim {args.feature_type_ref_delim}" if args.feature_type_ref_delim is not None else "",
    #                     f"--include-feature-list {args.include_feature_list}"  if args.include_feature_list is not None else "",
    #                     f"--exclude-feature-list {args.exclude_feature_list}"  if args.include_feature_list is not None else "",
    #                     f"--out-feature {ftrlist_tsv}",
    #                     "--log"])
    #     args.include_feature_list = ftrlist_tsv
    #     args.exclude_feature_list = None
    #     cmds.append(cmd)
    # * check --icols-mtx has the same number as --colnames-count
    # args.icols_mtx=str(args.icols_mtx)
    # if len(args.icols_mtx.split(",")) != len(args.colnames_count.split(",")):
    #     raise ValueError(f"The number of columns in --icols-mtx ({args.icols_mtx}) should be the same as the number of columns in --colnames-count ({args.colnames_count}).")
    # * convert sge to tsv (output: out_transcript, out_minmax, out_feature
    cmd = " ".join([f"{args.spatula} convert-sge",
                f"--in-sge {args.in_mex}",
                f"--out-tsv {args.out_dir}",
                f"--tsv-mtx {args.out_transcript}",
                f"--tsv-ftr {args.out_feature}",
                f"--tsv-minmax {args.out_minmax}",
                f"--colnames-count {args.colname_count}" if args.colname_count else ""])
    aux_argset = set(item for lst in [aux_sge_args["out"], aux_sge_args["inftr"], aux_sge_args["inbcd"], aux_sge_args["inmtx"], aux_sge_args["ftrname"]] for item in lst)
    cmd = add_param_to_cmd(cmd, args, aux_argset)
    cmd = add_mexparam_to_cmd(cmd, args, mexarg_mapping)
    cmds.append(cmd)
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
    #  * 10x_xenium: update default value for the phred score filtering 
    if args.platform == "10x_xenium":
        if args.csv_colname_phredscore is None:
            args.csv_colname_phredscore = "qv"
        if args.min_phred_score is None:
            args.min_phred_score = 20

    transcript_tsv = args.out_transcript.replace(".gz", "")
    cmd =  " ".join([f"cartloader format_generic",
                     f"--input {args.in_csv}",
                     f"--out-dir {args.out_dir}",
                     f"--out-transcript {transcript_tsv}",
                     f"--out-feature {args.out_feature}",
                     f"--out-minmax {args.out_minmax}",
                     f"--colname-count {args.colname_count}" if args.colname_count else ""])      
    # aux args
    aux_argset = set(item for lst in [aux_sge_args["out"], 
                                      aux_sge_args["incsv"], 
                                      aux_sge_args["ftrname"] #, aux_sge_args["ftrtype"]
                                      ] for item in lst)
    #aux_argset.add('print_removed_transcripts')
    cmd = add_param_to_cmd(cmd, args, aux_argset)
    # append to cmds
    cmds.append(cmd)
    cmds.append(f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}")
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

#================================================================================================
#
# density filtering functions
#
#================================================================================================
def sge_density_filtering(mm, sge_filtering_dict):
    genomic_feature=sge_filtering_dict["genomic_feature"]
    count_header=sge_filtering_dict["count_header"]
    # if genomic_feature is None:
    #     if len(count_header) > 1:
    #         logging.error("Missing --genomic-feature. Cannot use ---count with multiple columns as --genomic-feature. Please provide one column name for density-filtering.")
    #         sys.exit(1)
    #     genomic_feature = count_header[0]
    cmds=cmd_separator([], "Filtering the converted data by density...")
    cmd = " ".join([
            "ficture", "filter_by_density",
            f"--input {sge_filtering_dict['raw_transcript']}",
            f"--feature {sge_filtering_dict['raw_feature']}",
            f"--output", sge_filtering_dict["filtered_transcript"],
            f"--output_boundary  {sge_filtering_dict['filtered_prefix']}",
            f"--filter_based_on {genomic_feature}", 
            f"--mu_scale {sge_filtering_dict['mu_scale']}",
            f"--radius {sge_filtering_dict['radius']}", 
            f"--quartile {sge_filtering_dict['quartile']}",
            f"--hex_n_move {sge_filtering_dict['hex_n_move']}",
            f"--remove_small_polygons {sge_filtering_dict['polygon_min_size']}",
            f"--gene_header {' '.join(sge_filtering_dict['gene_header'])}",
            f"--count_header {' '.join(count_header)}",
        ])
    cmds.append(cmd)
    cmds.append(f"[ -f {sge_filtering_dict['filtered_transcript']} ] && [ -f {sge_filtering_dict['filtered_minmax']} ] && [ -f {sge_filtering_dict['filtered_feature']} ] && touch {sge_filtering_dict['flag']}")
    mm.add_target(sge_filtering_dict['flag'], sge_filtering_dict["prereq"], cmds)
    return mm

#================================================================================================
#
# visual functions
#
#================================================================================================

def sge_visual(mm, transcript_f, minmax_f, xy_f, prereq, spatula):
    scheck_app(spatula)
    # draw xy plot for visualization
    cmds = cmd_separator([], f"Drawing XY plot for SGE: {transcript_f}")
    #draw_cmd=f"{args.gzip} -dc {transcript_f} | tail -n +2 | cut -f 1,2 | {args.spatula} draw-xy --tsv /dev/stdin --out {xy_f}"
    cmds.append(f"XMIN=$(awk '/xmin/' {minmax_f}|cut -f 2) && \\")
    cmds.append(f"XMAX=$(awk '/xmax/' {minmax_f}|cut -f 2) && \\")
    cmds.append(f"YMIN=$(awk '/ymin/' {minmax_f}|cut -f 2) && \\")
    cmds.append(f"YMAX=$(awk '/ymax/' {minmax_f}|cut -f 2) && \\")
    draw_cmd = " ".join([
        f"'{spatula}'",
        "draw-xy",
        "--tsv", transcript_f,
        "--out", xy_f,
        "--icol-x 0", # str(icol_x),
        "--icol-y 1", # str(icol_y),
        "--icol-cnt -1", # str(icol_cnt) if icol_cnt is not None else "-1",
        "--ullr", "$XMIN,$YMIN,$XMAX,$YMAX",
        "--auto-adjust",
        "--skip-lines", "1",
    ])
    cmds.append(draw_cmd)
    mm.add_target(xy_f, prereq, cmds)
    return mm

def sge_visual_northup(mm, xy_f, xy_northup_f, minmax_f, prereq, srs="EPSG:3857", resample="cubic", gdalwarp="gdalwarp", gdal_translate="gdal_translate"):
    cmds = cmd_separator([], f"Create a north-up tif image for SGE: {xy_f}")
    temp_tif= xy_northup_f+".temp.tif"
    cmd = f"cartloader image_north_up --in-png {xy_f} --in-tif {temp_tif} --out-tif {xy_northup_f} --minmax {minmax_f} --srs {srs} --resample {resample} --gdalwarp {gdalwarp} --gdal_translate {gdal_translate}"
    cmds.append(cmd)
    cmds.append(f"rm {temp_tif}")
    mm.add_target(xy_northup_f, [xy_f]+prereq, cmds)
    return mm

#================================================================================================
#
# main functions
#
#================================================================================================

def sge_convert(_args):
    # args
    args=parse_arguments(_args)
    if args.exclude_feature_regex == "":
        args.exclude_feature_regex = None
    scheck_app(args.gzip)

    # input
    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"The input json file doesn't exist: {args.in_json}"

        print(f"Loading input files from input JSON {args.in_json}")        
        raw_data = load_file_to_dict(args.in_json)
        raw_tx = raw_data["TRANSCRIPT"]
        if raw_tx.endswith("parquet"):
            args.in_parquet = raw_tx
            print(f" * --in-parquet {args.in_parquet}")
        elif raw_tx.endswith("csv.gz") or raw_tx.endswith("tsv.gz") or raw_tx.endswith("tsv") or raw_tx.endswith("csv"):
            args.in_csv = raw_tx
            print(f" * --in-csv {args.in_csv}")
        args.scale_json = raw_data.get("SCALE", None)
    
    in_raw_filelist=input_by_platform(args)

    # output
    os.makedirs(args.out_dir, exist_ok=True)

    out_transcript_f = os.path.join(args.out_dir, args.out_transcript)
    out_minmax_f = os.path.join(args.out_dir, args.out_minmax)
    out_feature_f = os.path.join(args.out_dir, args.out_feature)
    out_xy_f = os.path.join(args.out_dir, args.out_xy)

    filtered_transcript_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.transcripts.unsorted.tsv.gz")
    filtered_feature_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.feature.lenient.tsv.gz")
    filtered_minmax_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.coordinate_minmax.tsv")
    filtered_xy_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.{args.out_xy}")

    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, "sge_assets.json")

    # params
    if args.platform == "10x_xenium":
        if args.exclude_feature_regex is None:
            # Negative probe patterns:
            #   BLANK_*
            #   DeprecatedCodeword_*
            #   NegControlCodeword_*
            #   NegControlProbe_*
            #   UnassignedCodeword_*
            args.exclude_feature_regex = "^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)"
            print(f"Update --exclude-feature-regex: {args.exclude_feature_regex }")
    
    # mm
    mm = minimake()

    # sge_convert
    sge_convert_flag = os.path.join(args.out_dir, "sge_convert.done")

    cmds = cmd_separator([], f"Converting input for raw data from: {args.platform}...")
    if args.platform == "10x_visium_hd":
        cmds = convert_visiumhd(cmds, args)
    elif args.platform == "seqscope":
        cmds = convert_seqscope(cmds, args)
    elif args.platform == "10x_xenium":
        #  * 10x_xenium: convert parquet to csv
        if args.in_parquet is not None and args.in_csv is not None:
            raise ValueError("For 10X Xenium, only one of the two input options (--in-parquet or --in-csv) can be specified for transcript data. Providing both will result in an error.")
        if args.in_parquet is not None:
            cmds = cmd_separator([], f"Converting input parquet into a csv file : (platform: {args.platform})...")
            args.in_csv = f"{args.out_dir}/transcripts.parquet.csv.gz"
            cmds.append(f"{args.parquet_tools} csv {args.in_parquet} |  {args.gzip} -c > {args.in_csv}")
            mm.add_target(args.in_csv, [args.in_parquet], cmds) 
            # update the prereq for convert
            in_raw_filelist=[args.in_csv] 
        cmds = convert_tsv(cmds, args)
    elif args.platform in ["cosmx_smi", "bgi_stereoseq", "vizgen_merscope", "pixel_seq", "nova_st", "generic"]:
        cmds = convert_tsv(cmds, args)
    cmds.append(f"[ -f {out_transcript_f} ] && [ -f {out_feature_f} ] && [ -f {out_minmax_f} ] && touch {sge_convert_flag}")
    mm.add_target(sge_convert_flag, in_raw_filelist, cmds)

    sge_assets={
        "transcript": out_transcript_f,
        "feature": out_feature_f,
        "minmax": out_minmax_f,
        "density_filtering": False
    }

    if args.sge_visual:
        mm = sge_visual(mm, 
                        out_transcript_f,
                        out_minmax_f,
                        out_xy_f,
                        [sge_convert_flag], 
                        args.spatula)
        if args.north_up:
            out_xyn_f= os.path.join(args.out_dir, args.out_northup_tif)
            mm = sge_visual_northup(mm, out_xy_f, out_xyn_f, out_minmax_f, [sge_convert_flag],
                                  srs=args.srs, resample=args.resample, gdalwarp=args.gdalwarp, gdal_translate=args.gdal_translate)
        sge_assets["visual"]={
            "png": out_xyn_f if args.north_up else out_xy_f,
            "northup": args.north_up
        }

    # filtering
    if args.filter_by_density:
        sge_filtered_flag = os.path.join(args.out_dir, "sge_density_filtering.done")
        sge_filtering_dict={
            "raw_transcript": out_transcript_f,
            "raw_feature": out_feature_f,
            "prereq": [sge_convert_flag],
            "filtered_transcript": filtered_transcript_f,
            "filtered_minmax": filtered_minmax_f,
            "filtered_feature": filtered_feature_f,
            "filtered_xy": filtered_xy_f,
            "filtered_prefix": os.path.join(args.out_dir, args.out_filtered_prefix),
            "flag": sge_filtered_flag,
            "gene_header": [args.colname_feature_name],
            "count_header": [args.colname_count], #args.colnames_count.split(","), #list
            "genomic_feature": args.colname_count, #args.genomic_feature,
            "mu_scale": 1,
            "radius": args.radius,
            "quartile": args.quartile,
            "hex_n_move": args.hex_n_move,
            "polygon_min_size": args.polygon_min_size,
        }
        mm = sge_density_filtering(mm, sge_filtering_dict)
    
        sge_assets["transcript"] = filtered_transcript_f
        sge_assets["feature"] = filtered_feature_f
        sge_assets["minmax"] = filtered_minmax_f
        sge_assets["density_filtering"] = True

        if args.sge_visual:
            mm = sge_visual(mm, 
                            filtered_transcript_f,
                            filtered_minmax_f,
                            filtered_xy_f,
                            [sge_filtered_flag],
                            args.spatula)
            if args.north_up:
                filtered_xyn_f= os.path.join(args.out_dir, f"{args.out_filtered_prefix}.{args.out_northup_tif}")
                mm = sge_visual_northup(mm, filtered_xy_f, filtered_xyn_f, filtered_minmax_f, [sge_filtered_flag],
                                        srs=args.srs, resample=args.resample, gdalwarp=args.gdalwarp, gdal_translate=args.gdal_translate)

            sge_assets["visual"]={
                "png": filtered_xyn_f if args.north_up else filtered_xy_f,
                "northup": args.north_up
            } 

    # write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    # write down a json file when execute
    write_dict_to_file(sge_assets, args.out_json, check_equal=True)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)



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


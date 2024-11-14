import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd

# get the path of the cu

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""
                                     Standardize Spatial Transcriptomics (ST) datasets into a transcript-indexed SGE in TSV format for downstream analysis.  
                                     1) Platform Support: 10X Visium HD, SeqScope, 10X Xenium, BGI Stereoseq, Cosmx SMI, Vizgen Merscope, Pixel-Seq, and Nova-ST
                                     2) Main functions: it generates a Makefile for a transcript-indexed SGE file in TSV format, a coordinate minmax TSV file, and a feature file counting UMIs per gene. 
                                     When executed with --dry-run, it generates only the Makefile without running it. 
                                     3) Additional Features: Includes options for filtering based on Phred-scaled quality score, gene name, and gene type, as well as converting pixel coordinates to micrometers
                                     """)
    #parser = argparse.ArgumentParser()
    # run params
    run_params = parser.add_argument_group("Run Options", "Run options for converting SGE files into a tsv format.")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_convert.mk", help='The file name of the Makefile to generate (default: sge_convert.mk)')
    
    # Key params
    key_params = parser.add_argument_group(
        "Key Parameters", 
        """
        Parameters to specify platform, units per um, and precision for the output files.
        Two ways to specify the conversion from units to micrometers:
        1) Use --units-per-um directly.
        2) For 10x_visium_hd datasets, use --scale-json to provide a scale json valid file to calculate the units per um.
        """)
    key_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st"], required=True, help='Platform of the raw input file to infer the format of the input file.')
    key_params.add_argument('--units-per-um', type=float, default=1.00, help='Coordinate unit per um (conversion factor) (default: 1.00)') 
    key_params.add_argument('--scale-json', type=str, default=None, help="For 10x_visium_hd datasets, users could use --scale-json to provide the path to the scale json file for calculating units-per-um (default: None). Typical naming convention: scalefactors_json.json")
    key_params.add_argument('--precision-um', type=int, default=2, help='Number of digits to store the transcript coordinates (only if --px_to_um is in use). Set it to 0 to round to integer (default: 2)')
    key_params.add_argument('--filter-by-density', action='store_true', default=False, help='(Optional) Filter the transcript-indexed SGE file by density (default: False)')

    # Output dir/file params
    output_params=parser.add_argument_group("Output Directory/File Parameters", "Output Parameters.")
    output_params.add_argument('--out-dir', type=str, required=True, help='The output directory, which will host the transcript-indexed SGE file/coordinate minmax TSV file/feature file, as well as the Makefile.')
    output_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='The converted compressed transcript-indexed SGE file in TSV format before filter_by_density (default: transcripts.unsorted.tsv.gz)')
    output_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='The converted coordinate minmax TSV file before filter_by_density (default: coordinate_minmax.tsv)')
    output_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='The converted file collects UMI counts on a per-gene basis before filter_by_density(default: feature.clean.tsv.gz)')
    output_params.add_argument('--out-filtered-transcript', type=str, default="filtered.transcripts.unsorted.tsv.gz", help='The filtered transcript-indexed SGE file in TSV format after filter_by_density')
    output_params.add_argument('--out-filtered-prefix', type=str, default="filtered", help='The prefix for filtered files (default: filtered)')

    # Input dir/file params
    input_params = parser.add_argument_group(
        "Input Directory/File Parameters",
        """
        Specify the input paths for the required files based on the platform you are using.
        1) For 10x_visium_hd: Specify the input paths for the SGE directory and parquet file. Additionally, provide the file names for barcode, feature, and matrix files in the SGE directory, if they differ from the defaults.
        2) For seqscope, specify the input path for the SGE directory.
        2) For 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, nova_st: Specify the path to the input raw CSV/TSV file.
        """
    )
    # - Arguments for 10x_visium_hd platform
    input_params.add_argument('--in-parquet', type=str, default="tissue_positions.parquet", help='Path to the input raw parquet file for spatial coordinates. Required for 10x_visium_hd platform (default: tissue_positions.parquet).') # naming convention: tissue_positions.parquet
    input_params.add_argument('--in-sge', type=str, default=os.getcwd(), help='Path to the input SGE directory. Required for 10x_visium_hd and seqscope platform. Defaults to the current working directory.')
    input_params.add_argument('--sge-bcd', type=str, default="barcodes.tsv.gz", help='Barcode file name in the SGE directory. Required for 10x_visium_hd platform if not using the default (barcodes.tsv.gz).')
    input_params.add_argument('--sge-ftr', type=str, default="features.tsv.gz", help='Feature file name in the SGE directory. Required for 10x_visium_hd platform if not using the default (features.tsv.gz).')
    input_params.add_argument('--sge-mtx', type=str, default="matrix.mtx.gz", help='Matrix file name in the SGE directory. Required for 10x_visium_hd platform if not using the default (matrix.mtx.gz).')
    # - Arguments for other platforms
    input_params.add_argument('--in-csv', type=str, default=None, help='Path to the input raw CSV/TSV file. Required for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq, and nova_st platforms (default: None).')

    # AUX input SGE params
    aux_in_sge_params = parser.add_argument_group(
        "IN-SGE Auxiliary Parameters",
        """
        Auxiliary parameters for processing input files when input files are SGE files. This applies to 10x_visium_hd and seqscope datasets.
        1) Use --icol-* parameters to specify the column indices in the input SGE files.
        2) Use --pos-* parameters to specify the column names and the delimiter for an additional barcode position file when required.
        3) Use --print-feature-id and --allow-duplicate-gene-names to customize the output.
        """)
    aux_in_sge_params.add_argument('--icols-mtx', type=str, default='1,2,3,4,5', help='Input column indices (comma-separated 1-based) in the SGE matrix file (default: 1,2,3,4,5).')
    aux_in_sge_params.add_argument('--icol-bcd-barcode', type=int, default=1, help='1-based column index of barcode in the SGE barcode file (default: 1).')
    aux_in_sge_params.add_argument('--icol-bcd-x', type=int, default=6, help='1-based column index of x coordinate in the SGE barcode file (default: 6).')
    aux_in_sge_params.add_argument('--icol-bcd-y', type=int, default=7, help='1-based column index of y coordinate in the SGE barcode file (default: 7).')
    aux_in_sge_params.add_argument('--icol-ftr-id', type=int, default=1, help='1-based column index of feature ID in the SGE feature file (default: 1).')
    aux_in_sge_params.add_argument('--icol-ftr-name', type=int, default=2, help='1-based column index of feature name in the SGE feature file (default: 2).')
    aux_in_sge_params.add_argument('--pos-colname-barcode', type=str, default='barcode', help='Column name for barcode in the additional barcode position file (default: barcode).')
    aux_in_sge_params.add_argument('--pos-colname-x', type=str, default='pxl_row_in_fullres', help='Column name for X-axis in the additional barcode position file (default: pxl_row_in_fullres).')
    aux_in_sge_params.add_argument('--pos-colname-y', type=str, default='pxl_col_in_fullres', help='Column name for Y-axis in the additional barcode position file (default: pxl_col_in_fullres).')
    aux_in_sge_params.add_argument('--pos-delim', type=str, default=',', help='Delimiter for the additional barcode position file (default: ",").')
    aux_in_sge_params.add_argument('--print-feature-id', action='store_true', help='Print feature ID in the output file (default: False).')
    aux_in_sge_params.add_argument('--allow-duplicate-gene-names', action='store_true', help='Allow duplicate gene names in the output file (default: False).')
    aux_in_sge_params.add_argument('--keep-mismatches', action='store_true', help='Keep mismatches in the output file (default: False). This applies to seqscope datasets only.')

    # AUX input csv params
    aux_in_csv_params = parser.add_argument_group(
        "IN-CSV Auxiliary Parameters", 
        """
        Auxiliary parameters for processing input files when input files are in TSV/CSV format. This applies to datasets from 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, pixel_seq or nova_st.
        1) Use --csv-colname-* parameters to specify the column names in the input CSV/TSV files.
        2) Use --csv-colname-phredscore with --min-phred-score to filter the input data based on the Phred-scaled quality score.
        Please note --csv-columns-x, --csv-columns-y, --csv-columns-feature-name are mandatory when the input file is in CSV format.
        """)
    aux_in_csv_params.add_argument('--csv-delim', type=str, default=None, help='Delimiter for the additional input tsv/csv file. Required if not using the default: 1) "," for 10x_xenium, cosmx_smi, and vizgen_merscope; 2)"\\t" for bgi_stereoseq, pixel_seq, and nova_st.')
    aux_in_csv_params.add_argument('--csv-colname-x',  type=str, default=None, help='Column name for X-axis (default: x_location for 10x_xenium; x for bgi_stereoseq; x_local_px for cosmx_smi; global_x for vizgen_merscope; xcoord for pixel_seq; x for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-y',  type=str, default=None, help='Column name for Y-axis (default: y_location for 10x_xenium; y for bgi_stereoseq; y_local_px for cosmx_smi; global_y for vizgen_merscope; ycoord for pixel_seq; y for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-feature-name', type=str, default=None, help='Column name for gene name (default: feature_name for 10x_xenium; geneID for bgi_stereoseq; target for cosmx_smi; gene for vizgen_merscope; geneName for pixel_seq; geneID for nova_st)')
    aux_in_csv_params.add_argument('--csv-colnames-count', type=str, default=None, help='Column name for feature expression count. If not provided, a count of 1 will be added for a feature in a pixel (default: MIDCounts for bgi_stereoseq; MIDCount for nova_st; None for the rest platforms).')
    aux_in_csv_params.add_argument('--csv-colname-feature-id', type=str, default=None, help='Column name for gene id (default: None)')
    aux_in_csv_params.add_argument('--csv-colnames-others', nargs='+', default=[], help='Columns names to keep (e.g., cell_id, overlaps_nucleus) (default: None)')
    aux_in_csv_params.add_argument('--csv-colname-phredscore', type=str, default=None, help='Column name for Phred-scaled quality value (Q-Score) estimating the probability of incorrect call (default: qv for 10x_xenium and None for the rest platforms).') # qv
    aux_in_csv_params.add_argument('--min-phred-score', type=float, default=None, help='Specify the Phred-scaled quality score cutoff (default: 20 for 10x_xenium and None for the rest platforms).')
    aux_in_csv_params.add_argument('--add-molecule-id', action='store_true', default=False, help='If enabled, a column of "molecule_id" will be added to the output file to track the index of the original input will be stored in (default: False).')
    
    # AUX output params
    aux_output_params = parser.add_argument_group(
        "Output Column Auxiliary Parameters",  
        "Auxiliary Output column parameters for the output files.")
    aux_output_params.add_argument('--colname-x', type=str, default='X', help='Column name for X (default: X)')
    aux_output_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y (default: Y)')
    aux_output_params.add_argument('--colnames-count', type=str, default='gn', help='Comma-separate column names for Count (default: gn)')
    aux_output_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for feature/gene name (default: gene)')
    aux_output_params.add_argument('--colname-feature-id', type=str, default=None, help='Column name for feature/gene ID. This is only required when --csv-colname-feature-id or --print-feature-id is applied (default: None)') 

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group(
        "Feature Filtering Auxiliary Parameters", 
        """
        Auxiliary parameters for filtering features by their name or type using an additional file, substring, or regex pattern:
        1) Use the ---feature-list, ---feature-substr, and --*-feature-regex parameters to exclude or include input data based on feature names.
        2) Use the --*-feature-type-regex parameter in conjunction with --csv-colname-feature-type or --feature-type-ref to filter input data based on feature type.
        """)
    aux_ftrfilter_params.add_argument('--include-feature-list', type=str, default=None, help='A file containing a list of input genes to be included (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-list', type=str, default=None, help='A file containing a list of input genes to be excluded (feature name of IDs) (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-substr', type=str, default=None, help='A substring of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-substr', type=str, default=None, help='A substring of feature/gene names to be excluded (default: None)')
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included (default: None)')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default="^(BLANK|Blank-|NegCon|NegPrb)", help='A regex pattern of feature/gene names to be excluded (default: "^(BLANK|Blank-|NegCon|NegPrb)"). To avoid filtering features by name using regex, use --exclude-feature-regex "". ')
    aux_ftrfilter_params.add_argument('--include-feature-type-regex', type=str, default=None, help='A regex pattern of feature/gene type to be included (default: None). To enable filtering genes by gene types, users must specify --csv-colname-feature-type or --feature-type-ref to provide gene type information') # (e.g. protein_coding|lncRNA)
    aux_ftrfilter_params.add_argument('--csv-colname-feature-type', type=str, default=None, help='The input column name in the input that corresponding to the gene type information, if your input file has gene type information (default: None)')
    aux_ftrfilter_params.add_argument('--feature-type-ref', type=str, default=None, help='Specify the path to a tab-separated gene information reference file to provide gene type information. The format should be: chrom, start position, end position, gene id, gene name, gene type (default: None)')
    aux_ftrfilter_params.add_argument('--print-removed-transcripts', action='store_true', default=False, help='For debugging purposes, print the list of features removed based on the filtering criteria (default: False). Currently it only works for datasets from 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope or pixel_seq.')

    # AUX polygon-filtering params
    aux_polyfilter_params = parser.add_argument_group('Polygon Filtering Auxiliary Parameters','Auxiliary parameters for filtering polygons based on the number of vertices. Only applicable when --filter-by-density is enabled.')
    aux_polyfilter_params.add_argument('--genomic-feature', type=str, default='gn', help='The column name of genomic feature for polygon-filtering.')
    aux_polyfilter_params.add_argument('--mu-scale', type=int, default=1, help='The scale factor for the polygon area calculation (default: 1.0)')
    aux_polyfilter_params.add_argument('--radius', type=int, default=15, help='The radius for the polygon area calculation (default: 15)')
    aux_polyfilter_params.add_argument('--quartile', type=int, default=2, help='The quartile for the polygon area calculation (default: 2)')
    aux_polyfilter_params.add_argument('--hex-n-move', type=int, default=1, help='The sliding step for polygon-filtering (default: 1)')
    aux_polyfilter_params.add_argument('--polygon-min-size', type=int, default=500, help='The minimum polygon size to be included in the output file (default: 500)')
    aux_polyfilter_params.add_argument('--gene-header', type=str, default='gene', help='The column names of gene name in the input file for polygon-filtering.')
    aux_polyfilter_params.add_argument('--count-header', type=str, default='gn', help='The column names of count in the input file for polygon-filtering.')
   
    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4".')
    env_params.add_argument('--spatula', type=str, default=None, help='Path to spatula binary. When not provided, it will use the spatula from the submodules.')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools binary.')
    
    # not in use 
    #aux_ftrfilter_params.add_argument('--unique', action='store_true', default=False, help='Merge pixels with (almost?) identical coordinates. Applies to cosmx_smi only.')
 
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def convert_in_by_platform(args):
    if args.platform == "10x_visium_hd":
        required_files = [os.path.join(args.in_sge, f) for f in [args.sge_bcd, args.sge_ftr, args.sge_mtx]] + [args.in_parquet]
    elif args.platform == "seqscope":
        required_files = [os.path.join(args.in_sge, f) for f in [args.sge_bcd, args.sge_ftr, args.sge_mtx]]
    elif args.platform in ["10x_xenium", "cosmx_smi", "bgi_stereoseq", "vizgen_merscope", "pixel_seq", "nova_st"]:
        required_files = [args.in_csv]
    return required_files

def create_minmax(cmds, args):
    # removed the px-to-um conversion.
    minmax_cmd = f"""{args.gzip} -cd {args.out_dir}/{args.out_transcript} | awk 'BEGIN{{FS=OFS="\\t"}} NR==1{{for(i=1;i<=NF;i++){{if($i=="X")x=i;if($i=="Y")y=i}}print $x,$y;next}}{{print $x,$y}}' | awk -F'\\t' ' BEGIN {{ min1 = "undef"; max1 = "undef"; min2 = "undef"; max2 = "undef"; }} {{ if (NR == 2 || $1 < min1) min1 = $1; if (NR == 2 || $1 > max1) max1 = $1; if (NR == 2 || $2 < min2) min2 = $2; if (NR == 2 || $2 > max2) max2 = $2; }} END {{ print "xmin\\t", min1; print "xmax\\t", max1; print "ymin\\t", min2; print "ymax\\t", max2; }}' > {args.out_dir}/{args.out_minmax}"""
    cmds.append(minmax_cmd)
    return cmds

# arg names
ans_out      = ['units_per_um', 'precision_um'] + ['colname_x', 'colname_y', 'colnames_count', 'colname_feature_name', 'colname_feature_id']
ans_insge    = ['sge_bcd', 'sge_ftr', 'sge_mtx']+ ['icols_mtx', 'icol_bcd_barcode', 'icol_bcd_x', 'icol_bcd_y', 'icol_ftr_id', 'icol_ftr_name',]
ans_inpos    = ['pos_colname_barcode', 'pos_colname_x', 'pos_colname_y', 'pos_delim']
ans_spatula  = ['print_feature_id', 'allow_duplicate_gene_names']
ans_incsv    = ['csv_delim', 'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name', 'csv_colnames_count', 'csv_colname_feature_id', 'csv_colnames_others', 'csv_colname_phredscore', 'min_phred_score', 'add_molecule_id']
ans_ftrname  = ['include_feature_list', 'exclude_feature_list', 'include_feature_substr', 'exclude_feature_substr', 'include_feature_regex', 'exclude_feature_regex'] 
ans_ftrtype  = ['include_feature_type_regex',  'csv_colname_feature_type', 'feature_type_ref']

#================================================================================================
#
# 10x_visium_hd and seqscope
#
#================================================================================================

def extract_unit2px_from_json(scale_json):
    # purpose: Extract the microns per pixel value from the scale json file and calculate the units per um.
    print(f"As --units-per-um-from-json is enabled, calculating units per um based on the microns per pixel value...")
    assert os.path.exists(scale_json), f"The scale json file ({scale_json}) does not exist. Please provide the correct path using --scale-json."
    with open(scale_json, 'r') as file:
        scale_data = json.load(file)
    # Extract the value for 'microns_per_pixel'
    microns_per_pixel = scale_data['microns_per_pixel']
    # microns_per_pixel cannot be NA or zero
    assert microns_per_pixel is not None, "The value for 'microns_per_pixel' is not found in the json file."
    assert microns_per_pixel != 0, "The value for 'microns_per_pixel' is zero. Please check the json file."
    print(f"    - Microns per pixel: {microns_per_pixel}")
    return 1/microns_per_pixel

def write_ftrlist_from_ftrtype(args):
    # purpose: Given the spatula doesn't support filtering based on gene type, this function prepares a feature list based on the gene type reference file.
    #           Note, when the --include-feature-type-regex is enabled, the --include-feature-list argument will be updated with the path to the prepared feature list.
    # Validate arguments
    assert args.feature_type_ref is not None, "The argument --include-feature-type-regex is enabled. Please provide the gene type reference file by --feature-type-ref."
    assert args.csv_colname_feature_type is None, "The argument --include-feature-type-regex is enabled. For 10x_visium_hd datasets, please use --feature-type-ref instead of --csv-colname-feature-type."
    # Read the gene type reference file
    df_ftrtype = pd.read_csv(args.feature_type_ref, sep="\t", header=None, names=["chrom", "start", "end", "gene_id", "gene_name", "gene_type"])
    # Filter df_ftrtype based on the gene type regex
    df_ftrtype = df_ftrtype[df_ftrtype["gene_type"].str.contains(args.include_feature_type_regex)]
    # If the include_feature_list argument is provided, merge the gene names from the feature list and the gene type reference file
    if args.include_feature_list is not None:
        df_ftrlist_raw = pd.read_csv(args.include_feature_list, header=None, names=["gene_name"])
        df_ftrlist     = pd.concat([df_ftrlist_raw, df_ftrtype["gene_name"]], ignore_index=True)
        df_ftrlist     = df_ftrlist.drop_duplicates()
    else:
        df_ftrlist     = df_ftrtype["gene_name"].drop_duplicates()
    # Save the feature list to a TSV file
    ftrlist_tsv = f"{args.out_dir}/include_feature_list.tsv"
    df_ftrlist.to_csv(ftrlist_tsv, sep="\t", header=False, index=False)
    # Update the include_feature_list argument with the path to the prepared feature list
    return ftrlist_tsv

def convert_visiumhd(cmds, args):
    # tools:
    # 1) spatula
    if args.spatula is None:
        args.spatula = os.path.join(cartloader_repo, "submodules", "spatula", "bin", "spatula")
        print(f"Given spatula is not provided, using the spatula from submodules: {args.spatula}.")
    scheck_app(args.spatula)
    # input: in_sge, in_parquet, scale_json
    # output: out_transcript, out_minmax, out_feature
    tmp_parquet = f"{args.out_dir}/tissue_positions.csv.gz"
    # 1) --in_parquet: convert parquet to csv
    cmds.append(f"{args.parquet_tools} csv {args.in_parquet} |  {args.gzip} -c > {tmp_parquet}")
    # 2) --scale_json: if applicable
    if args.scale_json is not None:
        args.units_per_um = extract_unit2px_from_json(args.scale_json)
    # 3) --included_feature_type_regex
    if args.include_feature_type_regex is not None:
        args.include_feature_list = write_ftrlist_from_ftrtype(args)
    # 4) convert sge to tsv (output: out_transcript, out_minmax, out_feature, (optional) out_sge)
    format_cmd=f"{args.spatula} convert-sge --in-sge {args.in_sge} --out-tsv {args.out_dir} --pos {tmp_parquet} --tsv-mtx {args.out_transcript} --tsv-ftr {args.out_feature} --tsv-minmax {args.out_minmax}"
    aux_argset = set(item for lst in [ans_out, ans_insge, ans_inpos, ans_spatula, ans_ftrname] for item in lst)
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    cmds.append(f"rm {tmp_parquet}")
    return cmds

def convert_seqscope(cmds, args):
    # tools:
    # 1) spatula
    if args.spatula is None:
        args.spatula = os.path.join(cartloader_repo, "submodules", "spatula", "bin", "spatula")
        print(f"Given spatula is not provided, using the spatula from submodules: {args.spatula}.")
    scheck_app(args.spatula)
    # input: in_sge
    # output: out_transcript, out_minmax, out_feature
    # 1) --included_feature_type_regex
    if args.include_feature_type_regex is not None:
        args.include_feature_list = write_ftrlist_from_ftrtype(args)
    # 2) convert sge to tsv (output: out_transcript, out_minmax, out_feature
    format_cmd=f"{args.spatula} convert-sge --in-sge {args.in_sge} --out-tsv {args.out_dir} --tsv-mtx {args.out_transcript} --tsv-ftr {args.out_feature} --tsv-minmax {args.out_minmax}"
    aux_argset = set(item for lst in [ans_out, ans_insge, ans_ftrname] for item in lst)
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    cmds.append(format_cmd)
    # 3) drop mismatches if --keep-mismatches is not enabled
    if not args.keep_mismatches:   
        drop_cmd = f"cartloader sge_drop_mismatches --in-dir {args.out_dir} --transcript {args.out_transcript} --feature {args.out_feature} --minmax {args.out_minmax} --gzip {args.gzip}"    
        cmds.append(drop_cmd)
    return cmds

#================================================================================================
#
# in-tsv/csv general functions
#
#================================================================================================

def update_csvformat_by_platform(args):
    platform_colnames_mapping={
        "x":{
            "10x_xenium": "x_location",
            "bgi_stereoseq": "x",
            "cosmx_smi": "x_global_px",
            "vizgen_merscope": "global_x",
            "pixel_seq": "xcoord",
            "nova_st": "x"
        },
        "y":{
            "10x_xenium": "y_location",
            "bgi_stereoseq": "y",
            "cosmx_smi": "y_global_px",
            "vizgen_merscope": "global_y",
            "pixel_seq": "ycoord",
            "nova_st": "y"
        },
        "feature_name":{
            "10x_xenium": "feature_name",
            "bgi_stereoseq": "geneID",
            "cosmx_smi": "target",
            "vizgen_merscope": "gene",
            "pixel_seq": "geneName",
            "nova_st": "geneID"
        },
        "count":{
            "10x_xenium": None,
            "bgi_stereoseq": "MIDCounts",
            "cosmx_smi": None,
            "vizgen_merscope": None,
            "pixel_seq": None,
            "nova_st": "MIDCount"
        }
    }
    platform_csvdelim_mapping={
        # If None, format_generic will use the default delimiter \t
        "10x_xenium": ",",
        "bgi_stereoseq": None,
        "cosmx_smi": ",",
        "vizgen_merscope":  ",",
        "pixel_seq": None,
        "nova_st": None
    }
    if args.csv_colname_x is None:
        args.csv_colname_x = platform_colnames_mapping["x"][args.platform]
    if args.csv_colname_y is None:
        args.csv_colname_y = platform_colnames_mapping["y"][args.platform]
    if args.csv_colname_feature_name is None:
        args.csv_colname_feature_name = platform_colnames_mapping["feature_name"][args.platform]
    if args.csv_colnames_count is None:
        args.csv_colnames_count = platform_colnames_mapping["count"][args.platform]
    if args.csv_delim is None:
        args.csv_delim = platform_csvdelim_mapping[args.platform]
    return args

def convert_tsv(cmds, args):
    # input: in_csv
    # output: out_transcript, out_minmax, out_feature
    # update csv_colname_*  and csv_delim based on the platform
    args = update_csvformat_by_platform(args)
    #  - 10x_xenium: update default value for the phred score filtering 
    if args.platform == "10x_xenium":
        if args.csv_colname_phredscore is None:
            args.csv_colname_phredscore = "qv"
        if args.min_phred_score is None:
            args.min_phred_score = 20
    # main commands
    transcript_tsv = args.out_transcript.replace(".gz", "")
    format_cmd=f"cartloader format_generic --input {args.in_csv} --out-dir {args.out_dir} --out-transcript {transcript_tsv} --out-feature {args.out_feature} --out-minmax {args.out_minmax}"
    # aux args
    aux_argset = set(item for lst in [ans_out, ans_incsv, ans_ftrname, ans_ftrtype] for item in lst)
    aux_argset.add('print_removed_transcripts')
    format_cmd = add_param_to_cmd(format_cmd, args, aux_argset)
    # append to cmds
    cmds.append(format_cmd)
    cmds.append(f"{args.gzip} -c {args.out_dir}/{transcript_tsv} > {args.out_dir}/{args.out_transcript}")
    cmds.append(f"rm {args.out_dir}/{transcript_tsv}")
    return cmds

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
    in_raw_filelist=convert_in_by_platform(args)
    for f in in_raw_filelist:
        assert f is not None, f"Missing input file for {f}. Please provide --{f}"
    # output
    os.makedirs(args.out_dir, exist_ok=True)

    # mm
    mm = minimake()

    # sge_convert
    cmds = cmd_separator([], f"Converting input for raw data from: {args.platform}...")
    if args.platform == "10x_visium_hd":
        cmds = convert_visiumhd(cmds, args)
    elif args.platform == "seqscope":
        cmds = convert_seqscope(cmds, args)
    elif args.platform in ["10x_xenium", "cosmx_smi", "bgi_stereoseq", "vizgen_merscope", "pixel_seq", "nova_st"]:
        cmds = convert_tsv(cmds, args)
    # add done if out_transcript, out_minmax, out_feature are generated
    sge_convert_flag = os.path.join(args.out_dir, "sge_convert.done")
    cmds.append(f"[ -f {os.path.join(args.out_dir, args.out_transcript)} ] && [ -f {os.path.join(args.out_dir, args.out_feature)} ] && [ -f {os.path.join(args.out_dir, args.out_minmax)} ] && touch {sge_convert_flag}")
    mm.add_target(sge_convert_flag, in_raw_filelist, cmds) 

    if args.filter_by_density:
        cmds=cmd_separator([], "Filtering the converted data by density...")
        cmd = " ".join([
                "ficture", "filter_by_density",
                f"--input {args.out_dir}/{args.out_transcript}",
                f"--feature {args.out_dir}/{args.out_feature}",
                f"--output {args.out_dir}/{args.out_filtered_transcript}",
                f"--output_boundary {args.out_dir}/{args.out_filtered_prefix}",
                f"--filter_based_on {args.genomic_feature}",
                f"--mu_scale {args.mu_scale}",
                f"--radius {args.radius}",
                f"--quartile {args.quartile}",
                f"--hex_n_move {args.hex_n_move}",
                f"--remove_small_polygons {args.polygon_min_size}",
                f"--gene_header \"{args.gene_header}\"",
                f"--count_header \"{args.count_header}\""
            ])
        cmds.append(cmd)
        mm.add_target(f"{args.out_dir}/{args.out_filtered_transcript}", [sge_convert_flag], cmds)


    # write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    mm.write_makefile(f"{args.out_dir}/{args.makefn}")
    # if args.platform == "seqscope":
    #     mm.write_makefile(f"{args.out_dir}/{args.makefn}") #, use_bash=True)
    # else:
    #     mm.write_makefile(f"{args.out_dir}/{args.makefn}")

    # run makefile
    if args.dry_run:
        os.system(f"make -f {args.out_dir}/{args.makefn} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")
    else:
        os.system(f"make -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")

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
import sys, os, argparse, logging,  inspect, json, subprocess
import pandas as pd
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, read_minmax, write_dict_to_file, load_file_to_dict, execute_makefile
from cartloader.utils.sge_helper import aux_sge_args, input_by_platform, update_csvformat_by_platform
from cartloader.scripts.feature_filtering import filter_feature_by_type

def parse_arguments(_args):

    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description=(
            "Standardize spatial transcriptomics inputs into SGE TSV outputs. "
            "Platform Supports: 10x Visium HD, SeqScope, 10x Xenium, BGI Stereo-seq, CosMx SMI, Vizgen MERSCOPE, Pixel-seq, Nova-ST, and a generic CSV/TSV option. "
            "Outputs: A transcript-indexed SGE file, a coordinate minmax TSV file, and a feature file counting UMIs per gene. All output are in micrometer precision."
            "Options: filtering SGE by gene, density, or quality score; SGE visualization"
        )
    )
    #parser = argparse.ArgumentParser()
    # run params
    run_params = parser.add_argument_group("Run Options", "Run options.")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate Makefile and print commands without executing')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--makefn', type=str, default="sge_convert.mk", help='File name of Makefile to write (default: sge_convert.mk)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    
    # Input/output/key params
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output paths and core settings.")
    inout_params.add_argument('--platform', type=str, choices=["10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic"], required=True, help='Input platform. Use "generic" for CSV/TSV from unsupported or custom sources')
    # - input
    inout_params.add_argument('--in-json', type=str, default=None, help='Path to input manifest JSON. If set, omits --in-parquet/--in-csv/--pos-parquet/--scale-json (platform: 10x_xenium, 10x_visium_hd)')
    inout_params.add_argument('--in-mex', type=str, default=os.getcwd(), help='Path to input MEX directory (platform: 10x Visium HD, SeqScope; default: current working directory)') # 10x_visium_hd, seqscope 
    inout_params.add_argument('--in-csv', type=str, default=None, help='Path to input CSV/TSV (platform: 10x Xenium, BGI Stereo-seq, CosMx SMI, Vizgen MERSCOPE, Pixel-seq, Nova-ST)') 
    inout_params.add_argument('--in-parquet', type=str, default=None, help='Path to input transcript parquet (platform: 10x Xenium)') 
    # - additional pos
    inout_params.add_argument('--pos-parquet', type=str, default=None, help='Path to input position parquet providing spatial coordinates (platform: 10x Visium HD; typical: tissue_positions.parquet)') # 10x_visium_hd
    # - scaling
    inout_params.add_argument('--scale-json', type=str, default=None, help='Path to input scale JSON. If set, defaults --units-per-um from microns_per_pixel in this JSON file (platform: 10x Visium HD; typical: scalefactors_json.json)') 
    inout_params.add_argument('--units-per-um', type=float, default=1.00, help='Coordinate units per µm in inputs (default: 1.00). For 10x Visium HD, prefer --scale-json')  
    # - output
    inout_params.add_argument('--out-dir', type=str, required=True, help='Path to output directory for converted SGE, filtered SGE, visualizations, and Makefile')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='File name of output transcript-indexed SGE TSV under --out-dir (default: transcripts.unsorted.tsv.gz)')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='File name of output compressed per-gene UMI count TSV under --out-dir (default: feature.clean.tsv.gz)')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='File name of output coordinate min/max TSV under --out-dir (default: coordinate_minmax.tsv)')
    inout_params.add_argument('--out-json',  type=str, default=None, help='Path to output JSON manifest of SGE paths (default: <out-dir>/sge_assets.json)')

    # AUX input MEX params
    aux_in_mex_params = parser.add_argument_group("IN-MEX Auxiliary Parameters", "10x Visium HD and SeqScope: auxiliary inputs")
    aux_in_mex_params.add_argument('--mex-bcd', type=str, default="barcodes.tsv.gz", help='File name of barcode in --in-mex (default: barcodes.tsv.gz)')
    aux_in_mex_params.add_argument('--mex-ftr', type=str, default="features.tsv.gz", help='File name of feature in --in-mex (default: features.tsv.gz)')
    aux_in_mex_params.add_argument('--mex-mtx', type=str, default="matrix.mtx.gz", help='File name of matrix in --in-mex (default: matrix.mtx.gz)')
    aux_in_mex_params.add_argument('--icols-mtx', type=int, default=1, help='1-based indices for the target genomic feature among count columns in --mex-mtx (default: 1). For example, when focusing on "gn" of SeqScope datasets, use --icols-mtx 1.')
    aux_in_mex_params.add_argument('--icol-ftr-id', type=int, default=1, help='1-based column index of feature ID in --mex-ftr (default: 1)')
    aux_in_mex_params.add_argument('--icol-ftr-name', type=int, default=2, help='1-based column index of feature name in --mex-ftr (default: 2)')
    aux_in_mex_params.add_argument('--icol-bcd-barcode', type=int, default=1, help='1-based column index of barcode in --mex-bcd (platform: SeqScope; default: 1)')
    aux_in_mex_params.add_argument('--icol-bcd-x', type=int, default=6, help='1-based column index of x coordinate in --mex-bcd (platform: SeqScope; default: 6)')
    aux_in_mex_params.add_argument('--icol-bcd-y', type=int, default=7, help='1-based column index of y coordinate in --mex-bcd (platform: SeqScope; default: 7)')
    aux_in_mex_params.add_argument('--pos-colname-barcode', type=str, default='barcode', help='Column name for barcode in --pos-parquet (platform: 10x Visium HD; default: barcode)')
    aux_in_mex_params.add_argument('--pos-colname-x', type=str, default='pxl_col_in_fullres', help='Column name for X coordinates in --pos-parquet (platform: 10x Visium HD; default: pxl_row_in_fullres)')
    aux_in_mex_params.add_argument('--pos-colname-y', type=str, default='pxl_row_in_fullres', help='Column name for Y coordinates in --pos-parquet (platform: 10x Visium HD; default: pxl_col_in_fullres)')
    aux_in_mex_params.add_argument('--pos-delim', type=str, default=',', help='Delimiter in --pos-parquet (platform: 10x Visium HD; default: ",")')
    #aux_in_mex_params.add_argument('--print-feature-id', action='store_true', help='Print feature ID in the output')
    aux_in_mex_params.add_argument('--allow-duplicate-gene-names', action='store_true', help='Allow duplicate gene names in the output')
    aux_in_mex_params.add_argument('--keep-mismatches', action='store_true', help='For debugging purposes, keep mismatches in the output (platform: SeqScope)')

    # AUX input csv params
    aux_in_csv_params = parser.add_argument_group("IN-CSV Auxiliary Parameters", "CSV/TSV inputs for Xenium/Stereo-seq/CosMx/MERSCOPE/Pixel-seq/Nova-ST")
    aux_in_csv_params.add_argument('--csv-comment', action='store_true', help="Enable to skip lines starts with # in the CSV file (default: False for 10x_xenium, bgi_stereoseq, cosmx_smi, vizgen_merscope, and pixel_seq, True for nova_st)")
    aux_in_csv_params.add_argument('--csv-delim', type=str, default=None, help='Delimiter in --in-csv (default: "," for 10x_xenium, cosmx_smi, and vizgen_merscope, "\\t" for bgi_stereoseq, pixel_seq, and nova_st) ')
    aux_in_csv_params.add_argument('--csv-colname-x', type=str, default=None, help='Column name for X coordinates in --in-csv (default: x_location for 10x_xenium, x for bgi_stereoseq, x_local_px for cosmx_smi, global_x for vizgen_merscope, xcoord for pixel_seq, x for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-y', type=str, default=None, help='Column name for Y coordinates in --in-csv (default: y_location for 10x_xenium, y for bgi_stereoseq, y_local_px for cosmx_smi, global_y for vizgen_merscope, ycoord for pixel_seq, y for nova_st)')
    aux_in_csv_params.add_argument('--csv-colname-count', type=str, default=None, help='Column name for expression count in --in-csv. If not provided, a count of 1 will be added for a feature in a pixel (default: MIDCounts for bgi_stereoseq, MIDCount for nova_st, None for the rest platforms).')
    aux_in_csv_params.add_argument('--csv-colname-feature-name', type=str, default=None, help='Column name for gene name in --in-csv (default: feature_name for 10x_xenium, geneID for bgi_stereoseq, target for cosmx_smi, gene for vizgen_merscope, geneName for pixel_seq, geneID for nova_st)')
    # aux_in_csv_params.add_argument('--csv-colname-feature-id', type=str, default=None, help='Column name for gene id')
    aux_in_csv_params.add_argument('--csv-colnames-others', nargs='*', default=[], help='Columns names to keep in --in-csv (e.g., cell_id, overlaps_nucleus)')
    aux_in_csv_params.add_argument('--csv-colname-phredscore', type=str, default=None, help='Column name for Phred-scaled quality value in --in-csv. This is also named as Q-Score, which estimates the probability of incorrect call (default: qv for 10x_xenium and None for the rest platforms)') # qv
    aux_in_csv_params.add_argument('--min-phred-score', type=float, default=None, help='Phred-scaled quality score cutoff (default: 20 for 10x_xenium and None for the rest platforms).')
    #aux_in_csv_params.add_argument('--add-molecule-id', action='store_true', default=False, help='If enabled, a column of "molecule_id" will be added to the output file to track the index of the original input will be stored in.')
    
    # AUX output params
    aux_out_params = parser.add_argument_group("Output Auxiliary Parameters")
    aux_out_params.add_argument('--precision-um', type=int, default=2, help='Precision for transcript coordinates in the output. Set it to 0 to round to integer (default: 2)')
    aux_out_params.add_argument('--colname-x', type=str, default='X', help='Column name for X in the output (default: X)')
    aux_out_params.add_argument('--colname-y', type=str, default='Y', help='Column name for Y in the output(default: Y)')
    aux_out_params.add_argument('--colname-count', type=str, default='count', help='Comma-separated column names for count in the output (default: count)')
    aux_out_params.add_argument('--colname-feature-name', type=str, default='gene', help='Column name for gene name in the output (default: gene)')
    # aux_out_params.add_argument('--colname-feature-id', type=str, default=None, help='Column name for gene ID. Required only when --csv-colname-feature-id or --print-feature-id is applied') 

    # AUX gene-filtering params
    aux_ftrfilter_params = parser.add_argument_group("Feature Filtering Parameters")
    aux_ftrfilter_params.add_argument('--include-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be included')
    aux_ftrfilter_params.add_argument('--exclude-feature-regex', type=str, default=None, help='A regex pattern of feature/gene names to be excluded (default: "^(BLANK_|Blank-|DeprecatedCodeword_|NegCon|UnassignedCodeword_)" for 10_xenium, None for the rest)')

    # AUX polygon-filtering params
    aux_polyfilter_params = parser.add_argument_group('Density/Polygon Filtering Parameters')
    # genomic feature, gene_header and count_header will be automatically based on the --colname-feature-name and --colnames-count
    aux_polyfilter_params.add_argument('--filter-by-density', action='store_true', default=False, help='Enable density-based filtering')
    aux_polyfilter_params.add_argument('--out-filtered-prefix', type=str, default="filtered", help='Prefix for filtered outputs and images under --out-dir (default: filtered)')
    aux_polyfilter_params.add_argument('--radius', type=int, default=15, help='Radius for the polygon area calculation (default: 15)')
    aux_polyfilter_params.add_argument('--quartile', type=int, default=2, help='Quartile for the polygon area calculation (default: 2)')
    aux_polyfilter_params.add_argument('--hex-n-move', type=int, default=1, help='Sliding step (default: 1)')
    aux_polyfilter_params.add_argument('--polygon-min-size', type=int, default=500, help='Minimum polygon size (default: 500)')

    # AUX visualization params
    aux_visual_params = parser.add_argument_group("Visualization Parameters")
    aux_visual_params.add_argument('--sge-visual', action='store_true', default=False, help='Generate SGE visualization; if filtering is enabled, visualize both. See "North-up Auxiliary Parameters"')
    aux_visual_params.add_argument('--north-up', action='store_true', default=False, help='Enable to visualize SGE in a tif image north-up.')
    aux_visual_params.add_argument('--out-xy', type=str, default="xy.png", help='File name of output visualization image under --out-dir (default: xy.png). Filtered image is prefixed by --out-filtered-prefix')
    aux_visual_params.add_argument('--out-northup-tif', type=str, default="xy_northup.tif", help='File name of output SGE visualization after north-up orientation under --out-dir (default: xy_northup.tif).')
    aux_visual_params.add_argument('--srs', type=str, default='EPSG:3857', help='Spatial reference system (used with --north-up; default: EPSG:3857)')
    aux_visual_params.add_argument('--resample', type=str, default='cubic', help='Resampling method (used with --north-up; options: near, bilinear, cubic, etc; default: cubic).')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Paths to external tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary (default: gzip). For speed, consider "pigz -p 4"')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary (default: spatula)')
    env_params.add_argument('--parquet-tools', type=str, default="parquet-tools", help='Path to parquet-tools (used with --in-parquet or --pos-parquet; default: parquet-tools)')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate (used with --north-up; default: gdal_translate)')
    env_params.add_argument('--gdalwarp', type=str, default=f"gdalwarp", help='Path to gdalwarp (used with --north-up; default: gdalwarp)')
 
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
    # purpose: Extract microns_per_pixel from scale JSON and derive units per µm.
    print(f"Deriving --units-per-um from --scale-json: {scale_json}")
    with open(scale_json, 'r') as file:
        scale_data = json.load(file)
    # Extract the value for 'microns_per_pixel'
    microns_per_pixel = scale_data['microns_per_pixel']
    # microns_per_pixel cannot be NA or zero
    assert microns_per_pixel is not None, f"Key not found: 'microns_per_pixel'. Check your scale JSON {scale_json} (--scale-json)"
    assert microns_per_pixel != 0, f"Invalid value: 'microns_per_pixel' == 0. Check your scale JSON {scale_json} (--scale-json)"
    print(f"microns_per_pixel: {microns_per_pixel}")
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
    ## input: in_mex, pos_parquet, scale_json
    ## output: out_transcript, out_minmax, out_feature
    # * --pos-parquet: convert parquet to csv
    tmp_parquet = f"{args.out_dir}/tissue_positions.csv.gz"
    cmds.append(f"{args.parquet_tools} csv {args.pos_parquet} |  {args.gzip} -c > {tmp_parquet}")
    # * --scale_json: if applicable
    if args.scale_json is not None:
        assert os.path.exists(args.scale_json), f"File not found: {args.scale_json} (--scale-json)"
        args.units_per_um = extract_unit2px_from_json(args.scale_json)
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
    cmd =  " ".join([f"cartloader sge_format_generic",
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
        assert os.path.exists(args.in_json), f"File not found: {args.in_json} (--in-json)"

        print(f"Reading inputs from manifest: {args.in_json}")        
        raw_data = load_file_to_dict(args.in_json)
        raw_sge = raw_data.get("SGE", raw_data) # use raw_data as default to support the flat dict build in the old scripts
        if args.platform == "10x_xenium":
            raw_tx = raw_sge["TRANSCRIPT"]
            if raw_tx.endswith("parquet"):
                args.in_parquet = raw_tx
                print(f"Resolved transcript source: --in-parquet {args.in_parquet}")
            elif raw_tx.endswith("csv.gz") or raw_tx.endswith("tsv.gz") or raw_tx.endswith("tsv") or raw_tx.endswith("csv"):
                args.in_csv = raw_tx
                print(f"Resolved transcript source: --in-csv {args.in_csv}")
        elif args.platform == "10x_visium_hd":
            args.in_mex = raw_sge.get("TRANSCRIPT_MEX", None)
            print(f"Resolved input: --in-mex {args.in_mex}")
            args.pos_parquet = raw_sge.get("POSITION", None)
            print(f"Resolved input: --pos-parquet {args.pos_parquet}")
            args.scale_json = raw_sge.get("SCALE", None)
            print(f"Resolved input: --scale-json {args.scale_json}")
    
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
            args.exclude_feature_regex = "^(BLANK_|Blank-|DeprecatedCodeword_|NegCon|UnassignedCodeword_)"
            print(f"Using --exclude-feature-regex: {args.exclude_feature_regex }")
    
    # mm
    mm = minimake()

    #  * 10x_xenium: convert parquet to csv
    if args.platform == "10x_xenium":
        if args.in_parquet is not None and args.in_csv is not None:
            raise ValueError("For 10X Xenium, only one of the two input options (--in-parquet or --in-csv) can be specified for transcript data. Providing both will result in an error.")
        if args.in_parquet is not None:
            cmds = cmd_separator([], f"Converting input parquet into a csv file : (platform: {args.platform})...")
            args.in_csv = f"{args.out_dir}/transcripts.parquet.csv.gz"
            cmds.append(f"{args.parquet_tools} csv {args.in_parquet} |  {args.gzip} -c > {args.in_csv}")
            mm.add_target(args.in_csv, [args.in_parquet], cmds) 
            # update the prereq for convert
            in_raw_filelist=[args.in_csv] 

    # sge_convert (from csv/mex)
    sge_convert_flag = os.path.join(args.out_dir, "sge_convert.done")

    cmds = cmd_separator([], f"Converting input for raw data from: {args.platform}...")
    if args.platform == "10x_visium_hd":
        cmds = convert_visiumhd(cmds, args)
    elif args.platform == "seqscope":
        cmds = convert_seqscope(cmds, args)
    elif args.platform in ["cosmx_smi", "bgi_stereoseq", "vizgen_merscope", "pixel_seq", "nova_st", "generic", "10x_xenium"]:
        cmds = convert_tsv(cmds, args)
    cmds.append(f"[ -f {out_transcript_f} ] && [ -f {out_feature_f} ] && [ -f {out_minmax_f} ] && touch {sge_convert_flag}")
    mm.add_target(sge_convert_flag, in_raw_filelist, cmds)

    # original assets (will be updated if density filtering is enabled)
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
            "count_header": [args.colname_count],
            "genomic_feature": args.colname_count,
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
        raise ValueError("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    # write down a json file when execute
    write_dict_to_file(sge_assets, args.out_json, check_equal=True)
    
    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

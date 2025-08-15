import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess

from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, log_dataframe, write_dict_to_file, execute_makefile
from cartloader.utils.minimake import minimake
from cartloader.scripts.sge_convert import sge_visual, sge_density_filtering, sge_visual_northup
from cartloader.utils.image_helper import update_orient

def sge_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description=""""
        Stitching multiple tile SGEs into one SGE.
        Outputs: A transcript-indexed SGE file, a coordinate minmax TSV file, and a feature file counting UMIs per feature. All output are in micro-meter precision.
        """
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_stitch.mk", help='Makefile name (default: sge_stitch.mk)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Parameters to specify platform, input, output, and units per um, precision, and density-filtering for output.")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of the input SGEs tiles in a format of <transcript>,<feature>,<minmax>,<row>,<col>,<rotate>,<vertical_flip>,<horizontal_flip>. These are referred to as SGE tiles to distinguish them from the combined output SGE. ")
    inout_params.add_argument("--out-dir", type=str, help="Output directory.")
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='Output for SGE stitch. File name of the compressed transcript-indexed SGE file in TSV format (default: transcripts.unsorted.tsv.gz).')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='Output for SGE stitch. File name of the coordinate minmax TSV file (default: coordinate_minmax.tsv).')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='Output for SGE stitch. File name of the compressed UMI count per feature TSV file (default: feature.clean.tsv.gz).')
    inout_params.add_argument('--out-tile-minmax', type=str, default="coordinate_minmax_per_tile.tsv", help='Output for SGE stitch. File name of the coordinate minmax per tile TSV file (default: coordinate_minmax_per_tile.tsv).')
    inout_params.add_argument('--out-json',  type=str, default="sge_assets.json", help='Output json summarizing SGE information. (default: sge_assets.json).')
    inout_params.add_argument('--feature-distribution', action='store_true', default=False, help='Action. Create a file to summarize the feature distribution across the input tiles')
    inout_params.add_argument('--out-feature-distribution', type=str, default="feature.distribution.tsv.gz", help='Output for --feature-distribution. File name of the output summarizing the feature distribution across the input tiles (default: feature.distribution.tsv.gz)')
    inout_params.add_argument('--filter-by-density', action='store_true', default=False, help='Action. Enable density-filtering for the output SGE(default: False). If enabled, check the density-filtering auxiliary parameters.')
    inout_params.add_argument('--out-filtered-prefix', type=str, default="filtered", help='Output for --filter-by-density. If --filter-by-density, define the prefix for filtered SGE (default: filtered)')

    key_params = parser.add_argument_group("Key Parameters", "Key Parameters")
    key_params.add_argument("--colnames-count", type=str, default="count", help="Comma-separated column names for count (default: count)")
    key_params.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    key_params.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    key_params.add_argument('--colname-x', type=str, default="X", help='X column name (default: X)')
    key_params.add_argument('--colname-y', type=str, default="Y", help='Y column name (default: Y)')
    key_params.add_argument('--colnames-others', nargs='*', default=[], help='Columns names to keep (e.g., cell_id, overlaps_nucleus) (default: None)')
    key_params.add_argument('--units-per-um', type=float, default=1.0, help='Units per um in the input transcript tsv files (default: 1.0)')
    key_params.add_argument("--precision", type=int, default=2, help="Precision for the output minmax and transcript files (default: 2)")

    # polygon-filtering params
    polyfilter_params = parser.add_argument_group('Density/Polygon Filtering Parameters','Parameters for filtering polygons based on the number of vertices.')
    polyfilter_params.add_argument('--genomic-feature', type=str, default=None, help='Column name of genomic feature for polygon-filtering (default to --colnames-count if --colnames-count only specifies one column)')
    polyfilter_params.add_argument('--mu-scale', type=int, default=1, help='Scale factor for the polygon area calculation (default: 1.0)')
    polyfilter_params.add_argument('--radius', type=int, default=15, help='Radius for the polygon area calculation (default: 15)')
    polyfilter_params.add_argument('--quartile', type=int, default=2, help='Quartile for the polygon area calculation (default: 2)')
    polyfilter_params.add_argument('--hex-n-move', type=int, default=1, help='Sliding step (default: 1)')
    polyfilter_params.add_argument('--polygon-min-size', type=int, default=500, help='The minimum polygon size (default: 500)')
    # gene_header and count_header will be automatically based on the --colname-feature-name, --colname-feature-id and --colnames-count
    
    # AUX visualization params
    visual_params = parser.add_argument_group("SGE Visualization Parameters", "Parameters for visualizing the output SGE in a north-up orientation.")
    visual_params.add_argument('--sge-visual', action='store_true', default=False, help='(Optional) Plot the SGE in a greyscale PNG file. (default: False)')
    visual_params.add_argument('--out-xy', type=str, default="xy.png", help='Output for SGE visualization image (default: xy.png)')
    visual_params.add_argument('--north-up', action='store_true', default=False, help='If enabled, the sge will be visualized in a tif image north-up (default: False).')
    visual_params.add_argument('--out-northup-tif', type=str, default="xy_northup.tif", help='Output for SGE visualization. The prefix for the output north-up image (default: north_up.tif)')
    visual_params.add_argument('--srs', type=str, default='EPSG:3857', help='If --north-up, define the spatial reference system (default: EPSG:3857)')
    visual_params.add_argument('--resample', type=str, default='cubic', help='Define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4".')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary.')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='If --sge-visual with --north-up, provide path to gdal_translate binary')
    env_params.add_argument('--gdalwarp', type=str, default=f"gdalwarp", help='If --sge-visual with --north-up, provide path to gdalwarp binary')

    args = parser.parse_args(_args)

    mm = minimake()

    updated_tiles = []
    updated_ftrs = []
    updated_minmax = []
    prerequisities = []

    os.makedirs(args.out_dir, exist_ok=True)
    
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, "sge_assets.json")

    # sge adds on
    tile_indices=[]
    for in_tile in args.in_tiles:
        transcript, feature, minmax, row, col, rotate, vflip, hflip = in_tile.split(",")

        tile_idx=f"r{row}c{col}"
        if tile_idx not in tile_indices:
            tile_indices.append(tile_idx)
        else:
            raise ValueError(f"Found duplicate layout index (row {row}, col {col}).")
        
        if not os.path.exists(transcript):
            raise FileNotFoundError(f"Transcript file {transcript} does not exist.")
        
        # missing feature or minmax
        missing_ftr = feature is None or feature == "None" or feature == ""
        missing_minmax = minmax is None or minmax == "None" or minmax == ""
        if missing_ftr or missing_minmax:
            in_dir = os.path.dirname(transcript)
            if os.path.basename(transcript) == "transcripts.unsorted.tsv.gz":
                in_id = None
            else:
                in_id = os.path.basename(transcript).replace('.tsv.gz', '').replace(".transcripts", "").replace(".transcript", "")
        if missing_minmax: 
            cmds = cmd_separator([], f"Creating missing minmax file for {transcript}")
            out_minmax = os.path.join(in_dir, f"{in_id}.minmax.tsv") if in_id else os.path.join(in_dir, f"coordinate_minmax.tsv")
            add_minmax_cmd =" ".join([f"cartloader", "sge_adds_on",
                                    f"--in-transcript {transcript}",
                                    "--add-minmax",
                                    f"--out-minmax {out_minmax}"
                                ])
            cmds.append(add_minmax_cmd)
            mm.add_target(out_minmax, [transcript], cmds)
            minmax = out_minmax
        if missing_ftr: 
            cmds = cmd_separator([], f"Creating missing feature file for {transcript}")
            out_ftr = os.path.join(in_dir, f"{in_id}.feature.tsv.gz") if in_id else os.path.join(in_dir, f"feature.clean.tsv.gz")
            add_ftr_cmd =" ".join([f"cartloader", "sge_adds_on",
                                    f"--in-transcript {transcript}",
                                    "--add-feature",
                                    f"--colname-feature-name {args.colname_feature_name}",
                                    f"--colname-feature-id {args.colname_feature_id}" if args.colname_feature_id else "",
                                    f"--colname-count {args.colname_count}",
                                    f"--out-feature {out_ftr}"
                                ])
            cmds.append(add_ftr_cmd)
            mm.add_target(out_ftr, [transcript], cmds)
            feature = out_ftr

        # orientation
        rotate = str(rotate) if rotate else None
        vflip = str(vflip).lower() == "true"
        hflip = str(hflip).lower() == "true"
        rotate, vflip, hflip = update_orient(rotate, vflip, hflip, f"tile:{tile_idx}")
        
        if rotate is not None or vflip or hflip:
            ori_abbr = [tile_idx]
            if rotate is not None:
                ori_abbr.append(f"rot{rotate}")
            if vflip:
                ori_abbr.append("vfip")
            if hflip:
                ori_abbr.append("hfip")
            ori_prefix = os.path.join(args.out_dir, ".".join(ori_abbr))

            ori_transcript = f"{ori_prefix}.transcripts.unsorted.tsv.gz"
            ori_feature = f"{ori_prefix}.feature.clean.tsv.gz"
            ori_minmax = f"{ori_prefix}.coordinate_minmax.tsv"
            ori_flag = f"{ori_prefix}.done"
            cmds = cmd_separator([], f"Orientate tile {tile_idx} with rotation {rotate}, vertical flip {vflip}, horizontal flip {hflip}")
            orientate_cmd=" ".join([f"cartloader", "sge_orientate",
                                    f"--in-transcript {transcript}",
                                    f"--out-transcript {ori_transcript}",
                                    f"--in-minmax {minmax}",
                                    f"--out-minmax {ori_minmax}",
                                    f"--in-feature {feature}",
                                    f"--out-feature {ori_feature}",
                                    f"--rotate {rotate}",
                                    f"--flip-horizontal" if hflip else "",
                                    f"--flip-vertical" if vflip else "",
                                ])
            cmds.append(orientate_cmd)
            cmds.append(f'[ -f {ori_transcript} ] && [ -f {ori_feature} ] && [ -f {ori_minmax} ] && touch {ori_flag}')
            mm.add_target(ori_flag, [transcript, feature, minmax], cmds)
            
            transcript = ori_transcript
            feature = ori_feature
            minmax = ori_minmax
            prerequisities.append(ori_flag)
        else:
            prerequisities.extend([transcript, feature, minmax])
        
        updated_tiles.append(",".join([transcript, feature, minmax, row, col]))
        updated_ftrs.append(",".join([feature, row, col]))
        updated_minmax.append(",".join([minmax, row, col]))

    # coordinate transform from local to global and from unit to um 
    tile_minmax = os.path.join(args.out_dir, args.out_tile_minmax)
    cmds = cmd_separator([], f"Transforming local coordinates to global coordinates, and converting to um")
    transform_coord_cmd=" ".join([f"cartloader", "sge_tile_coord_transform",
                            f"--in-tiles {' '.join(updated_minmax)}",
                            f"--output {tile_minmax}",
                            f"--units-per-um {args.units_per_um}",
                        ])
    cmds.append(transform_coord_cmd)
    mm.add_target(tile_minmax, prerequisities, cmds)

    # sge stitch 
    out_transcript_f = os.path.join(args.out_dir, args.out_transcript)
    out_minmax_f = os.path.join(args.out_dir, args.out_minmax)
    out_feature_f = os.path.join(args.out_dir, args.out_feature)
    out_xy_f = os.path.join(args.out_dir, args.out_xy)

    # combine sges
    cmds = cmd_separator([], f"Combining SGEs")
    combine_cmd=" ".join([f"cartloader", "sge_combine_tiles",
                            f"--in-tiles {' '.join(updated_tiles)}",
                            f"--in-tile-minmax {tile_minmax}",
                            f"--out-dir {args.out_dir}",
                            f"--out-transcript {args.out_transcript}",
                            f"--out-minmax {args.out_minmax}",
                            f"--out-feature {args.out_feature}",
                            f"--colname-count {args.colname_count}",
                            f"--colname-feature-name {args.colname_feature_name}",
                            f"--colname-feature-id {args.colname_feature_id}" if args.colname_feature_id else "",
                            f"--colname-x {args.colname_x}",
                            f"--colname-y {args.colname_y}",
                            f"--units-per-um {args.units_per_um}",
                            f"--precision {args.precision}",
                            f"--colnames-others {' '.join(args.colnames_others)}"
                        ])
    cmds.append(combine_cmd)
    sge_stitch_flag = os.path.join(args.out_dir, "sge_stitch.done")
    cmds.append(f'[ -f {out_transcript_f} ] && [ -f {out_feature_f} ] && [ -f {out_minmax_f} ] && touch {sge_stitch_flag}')
    mm.add_target(sge_stitch_flag, prerequisities+[tile_minmax], cmds)

    sge_assets={
            "transcript": out_transcript_f,
            "feature": out_feature_f,
            "minmax": out_minmax_f,
            "density_filtering": False
        }

    # filter by density
    if args.filter_by_density:
        filtered_transcript_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.transcripts.unsorted.tsv.gz")
        filtered_minmax_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.coordinate_minmax.tsv")
        filtered_feature_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.feature.lenient.tsv.gz")
        filtered_xy_f = os.path.join(args.out_dir, f"{args.out_filtered_prefix}.xy.png")
        
        sge_filtered_flag = os.path.join(args.out_dir, "sge_density_filtering.done")
        sge_filtering_dict={
            "raw_transcript": out_transcript_f,
            "raw_feature": out_feature_f,
            "prereq": [sge_stitch_flag],
            "filtered_transcript": filtered_transcript_f,
            "filtered_minmax": filtered_minmax_f,
            "filtered_feature": filtered_feature_f,
            "filtered_xy": filtered_xy_f,
            "filtered_prefix": os.path.join(args.out_dir, args.out_filtered_prefix),
            "flag": sge_filtered_flag,
            "gene_header": [args.colname_feature_name] if args.colname_feature_id is None else [args.colname_feature_name, args.colname_feature_id], # list
            "count_header": [args.colname_count],
            "genomic_feature": args.colname_count,
            "mu_scale": args.mu_scale,
            "radius": args.radius,
            "quartile": args.quartile,
            "hex_n_move": args.hex_n_move,
            "polygon_min_size": args.polygon_min_size,
        }
        mm = sge_density_filtering(mm, sge_filtering_dict)

        # update sge_assets
        sge_assets["transcript"] = filtered_transcript_f
        sge_assets["feature"] = filtered_feature_f
        sge_assets["minmax"] = filtered_minmax_f
        sge_assets["density_filtering"] = True

    # generate feature distribution file
    out_dist_f = os.path.join(args.out_dir, "feature.distribution.tsv.gz")    
    cmds = cmd_separator([], f"Creating feature distribution file: {out_dist_f}")
    overlap_cmd=" ".join([f"cartloader", "feature_distribution",
                        f"--in-tiles {' '.join(updated_ftrs)}",
                        f"--colname-feature-name {args.colname_feature_name}",
                        f"--colname-count {args.colname_count}",
                        f"--output {out_dist_f}"
                    ])
    cmds.append(overlap_cmd)
    mm.add_target(out_dist_f, prerequisities, cmds)

    # draw xy plot for visualization
    if args.sge_visual:
        mm = sge_visual(mm, out_transcript_f, out_minmax_f, out_xy_f, [sge_stitch_flag], args.spatula)
        if args.north_up:
            out_xyn_f = os.path.join(args.out_dir, args.out_northup_tif)
            mm = sge_visual_northup(mm, out_xy_f, out_xyn_f, out_minmax_f, [sge_stitch_flag], 
                                    srs=args.srs, resample=args.resample, gdalwarp=args.gdalwarp, gdal_translate=args.gdal_translate)
        sge_assets["visual"]={
            "png": out_xyn_f if args.north_up else out_xy_f,
            "northup": args.north_up
        }

        if args.filter_by_density:
            mm = sge_visual(mm, filtered_transcript_f, filtered_minmax_f, filtered_xy_f, [sge_filtered_flag], args.spatula)
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
    
    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

    # write down a json file when execute
    write_dict_to_file(sge_assets, args.out_json, check_equal=True)


if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])

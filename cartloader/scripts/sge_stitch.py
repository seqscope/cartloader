import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess

from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, log_dataframe, write_dict_to_file, execute_makefile
from cartloader.utils.minimake import minimake
from cartloader.scripts.sge_convert import sge_visual, sge_density_filtering, sge_visual_northup
from cartloader.utils.orient_helper import update_orient

def sge_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description=""""
        Stitch multiple SGE tiles into one SGE.
        Outputs: transcript-indexed SGE, coordinate min/max TSV, and per-feature UMI counts; all in micrometer precision.
        """
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate Makefile and print commands without executing')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_stitch.mk", help='File name of Makefile to write (default: sge_stitch.mk)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads per job (default: 1)')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output paths and core settings.")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help='List of input SGE tiles as "<transcript_path>,<feature_path>,<minmax_path>,<row>,<col>,<rotate>,<vertical_flip>,<horizontal_flip>"')
    inout_params.add_argument("--out-dir", type=str, help='Path to output directory for stitched SGE, filtered outputs, visualizations, and Makefile')
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='File name of output transcript-indexed SGE TSV under --out-dir (default: transcripts.unsorted.tsv.gz)')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='File name of output per-gene UMI count TSV under --out-dir (default: feature.clean.tsv.gz)')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='File name of output coordinate min/max TSV under --out-dir (default: coordinate_minmax.tsv)')
    inout_params.add_argument('--out-tile-minmax', type=str, default="coordinate_minmax_per_tile.tsv", help='File name of output per-tile coordinate min/max TSV under --out-dir (default: coordinate_minmax_per_tile.tsv)')
    inout_params.add_argument('--out-json',  type=str, default=None, help='Path to output JSON summarizing SGE assets (default: <out-dir>/sge_assets.json)')

    key_params = parser.add_argument_group("Key Parameters", "Key output column names and numeric settings.")
    key_params.add_argument("--colname-count", type=str, default="count", help='Comma-separated output column name(s) for UMI count (default: count)')
    key_params.add_argument('--colname-feature-name', type=str, default='gene', help='Output column name for feature (gene) name (default: gene)')
    key_params.add_argument('--colname-feature-id', type=str, default=None, help='Output column name for feature (gene) ID')
    key_params.add_argument('--colname-x', type=str, default="X", help='Output column name for X (default: X)')
    key_params.add_argument('--colname-y', type=str, default="Y", help='Output column name for Y (default: Y)')
    key_params.add_argument('--colnames-others', nargs='*', default=[], help='Output column names to preserve (e.g., cell_id, overlaps_nucleus) (default: [])')
    key_params.add_argument('--units-per-um', type=float, default=1.0, help='Coordinate units per µm in inputs (default: 1.0)')
    key_params.add_argument("--precision", type=int, default=2, help='Precision for coordinates in stitched outputs (default: 2)')

    # polygon-filtering params
    polyfilter_params = parser.add_argument_group('Density/Polygon Filtering Parameters')
    polyfilter_params.add_argument('--filter-by-density', action='store_true', default=False, help='Enable density-based filtering (default: False). See "Density/Polygon Filtering Parameters"')
    polyfilter_params.add_argument('--out-filtered-prefix', type=str, default="filtered", help='Prefix for filtered outputs under --out-dir (default: filtered)')
    polyfilter_params.add_argument('--radius', type=int, default=15, help='Radius for the polygon area calculation (default: 15)')
    polyfilter_params.add_argument('--quartile', type=int, default=2, help='Quartile for the polygon area calculation (default: 2)')
    polyfilter_params.add_argument('--hex-n-move', type=int, default=1, help='Sliding step (default: 1)')
    polyfilter_params.add_argument('--polygon-min-size', type=int, default=500, help='The minimum polygon size (default: 500)')
    # gene_header and count_header will be automatically based on the --colname-feature-name, --colname-feature-id and --colnames-count
    
    # AUX visualization params
    visual_params = parser.add_argument_group("SGE Visualization Parameters")
    visual_params.add_argument('--sge-visual', action='store_true', default=False, help='Enable visualization; if filtering is enabled, visualizes both unfiltered and filtered SGEs')
    visual_params.add_argument('--out-xy', type=str, default="xy.png", help='File name of output visualization image under --out-dir (default: xy.png). Filtered image is prefixed with --out-filtered-prefix')
    visual_params.add_argument('--north-up', action='store_true', default=False, help='If enabled, produce north-up GeoTIFF visualization')
    visual_params.add_argument('--out-northup-tif', type=str, default="xy_northup.tif", help='File name of north-up GeoTIFF under --out-dir (default: xy_northup.tif)')
    visual_params.add_argument('--srs', type=str, default='EPSG:3857', help='If --north-up, spatial reference system (default: EPSG:3857)')
    visual_params.add_argument('--resample', type=str, default='cubic', help='Resampling method (default: cubic). Options: near, bilinear, cubic, etc.')

    ftrdis_params = parser.add_argument_group("Feature Distribution Parameters")
    ftrdis_params.add_argument('--feature-distribution', action='store_true', default=False, help='Enable feature distribution summary across input tiles')
    ftrdis_params.add_argument('--out-feature-distribution', type=str, default="feature.distribution.tsv.gz", help='File name of output feature distribution TSV under --out-dir (default: feature.distribution.tsv.gz)')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for external tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster compression, use "pigz -p 4" (default: gzip)')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary (default: spatula)')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='If --sge-visual with --north-up, path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdalwarp', type=str, default=f"gdalwarp", help='If --sge-visual with --north-up, path to gdalwarp binary (default: gdalwarp)')

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
            raise FileNotFoundError(f"File not found: {transcript} (--in-tiles)")
        
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
            add_minmax_cmd =" ".join([f"cartloader", "sge_supp",
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
            add_ftr_cmd =" ".join([f"cartloader", "sge_supp",
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
            "mu_scale": 1,
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

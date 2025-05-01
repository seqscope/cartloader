import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess

from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, log_dataframe
from cartloader.utils.minimake import minimake

def sge_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description=""""
        Stitching multiple tile SGEs into one SGE.
        Outputs: A transcript-indexed SGE file, a coordinate minmax TSV file, and a feature file counting UMIs per gene. All output are in micro-meter precision.

        """
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_stitch.mk", help='Makefile name (default: sge_stitch.mk)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output Parameters")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of the input tiles in a format of <transcript>,<feature>,<minmax>,<row>,<col>. ")
    inout_params.add_argument("--out-dir", type=str, help="Output directory.")
    inout_params.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='Output for SGE stitch. The compressed transcript-indexed SGE file in TSV format (default: transcripts.unsorted.tsv.gz).')
    inout_params.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='Output for SGE stitch. The coordinate minmax TSV file (default: coordinate_minmax.tsv).')
    inout_params.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='Output for SGE stitch. The compressed UMI count per gene TSV file (default: feature.clean.tsv.gz).')
    inout_params.add_argument('--out-tile-minmax', type=str, default="coordinate_minmax_per_tile.tsv", help='Output for SGE stitch. The coordinate minmax per tile TSV file (default: coordinate_minmax_per_tile.tsv).')
    inout_params.add_argument("--colnames-count", type=str, default="count", help="Comma-separated column names for count (default: count)")
    inout_params.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    inout_params.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    inout_params.add_argument('--colname-x', type=str, default="X", help='X column name (default: X)')
    inout_params.add_argument('--colname-y', type=str, default="Y", help='Y column name (default: Y)')
    inout_params.add_argument('--colnames-others', nargs='*', default=[], help='Columns names to keep (e.g., cell_id, overlaps_nucleus) (default: None)')
    inout_params.add_argument('--units-per-um', type=float, default=1.0, help='Units per um in the input transcript tsv files (default: 1.0)')
    inout_params.add_argument("--precision", type=int, default=2, help="Precision for the output minmax and transcript files (default: 2)")
    inout_params.add_argument('--sge-visual', action='store_true', default=False, help='Plot the SGE in a PNG file. (default: False)')
    inout_params.add_argument('--generate-tile-minmax-only', action='store_true', default=False, help='Only generate the tile minmax file without stitching the SGE. (default: False)')

    # env params
    env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4".')
    env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary.')
    args = parser.parse_args(_args)

    mm = minimake()

    updated_tiles = []
    prerequisities = []

    os.makedirs(args.out_dir, exist_ok=True)

    # sge adds on
    layout_indices=[]
    for in_tile in args.in_tiles:
        transcript, feature, minmax, row, col = in_tile.split(",")
        
        layout_idx=f"{row},{col}"
        if layout_idx not in layout_indices:
            layout_indices.append(layout_idx)
        else:
            raise ValueError(f"Found duplicate layout index (row {row}, col {col}).")
        
        if not os.path.exists(transcript):
            raise FileNotFoundError(f"Transcript file {transcript} does not exist.")
        
        # missing feature or minmax
        missing_ftr = feature is None or feature == "None"
        missing_minmax = minmax is None or minmax == "None"
        if missing_ftr or missing_minmax:
            in_dir = os.path.dirname(transcript)
            in_id = os.path.basename(transcript).replace('.tsv.gz', '').replace(".transcripts", "").replace(".transcript", "")
        if missing_ftr and not args.generate_tile_minmax_only:
            cmds = cmd_separator([], f"Creating missing feature file for {transcript}")
            out_ftr = os.path.join(in_dir, f"{in_id}.feature.tsv.gz")
            add_ftr_cmd =" ".join([f"cartloader", "sge_adds_on",
                                    f"--in-transcript {transcript}",
                                    "--add-feature",
                                    f"--colname-feature-name {args.colname_feature_name}",
                                    f"--colname-feature-id {args.colname_feature_id}" if args.colname_feature_id else "",
                                    f"--colnames-count {args.colnames_count}",
                                    f"--out-feature {out_ftr}"
                                ])
            cmds.append(add_ftr_cmd)
            mm.add_target(out_ftr, [transcript], cmds)
            feature = out_ftr
        if missing_minmax: # should always be true 
            cmds = cmd_separator([], f"Creating missing minmax file for {transcript}")
            out_minmax = os.path.join(in_dir, f"{in_id}.minmax.tsv")
            add_minmax_cmd =" ".join([f"cartloader", "sge_adds_on",
                                    f"--in-transcript {transcript}",
                                    "--add-minmax",
                                    f"--out-minmax {out_minmax}"
                                ])
            cmds.append(add_minmax_cmd)
            mm.add_target(out_minmax, [transcript], cmds)
            minmax = out_minmax
        updated_tiles.append(",".join([transcript, feature, minmax, row, col]))
        prerequisities.extend([transcript, feature, minmax])

    # coordinate transform from local to global and from unit to um 
    tile_minmax = os.path.join(args.out_dir, args.out_tile_minmax)
    cmds = cmd_separator([], f"Transforming local coordinates to global coordinates, and converting to um")
    transform_coord_cmd=" ".join([f"cartloader", "sge_tile_coord_transform",
                            f"--in-tiles {' '.join(updated_tiles)}",
                            f"--output {tile_minmax}",
                            f"--units-per-um {args.units_per_um}",
                        ])
    cmds.append(transform_coord_cmd)
    mm.add_target(tile_minmax, prerequisities, cmds)
    
   
    if not args.generate_tile_minmax_only:
         # combine sges
        cmds = cmd_separator([], f"Combining SGEs")
        if "," in args.colnames_count:
            colnames_count = args.colnames_count.split(",")
        else:
            colnames_count = [args.colnames_count]
        combine_cmd=" ".join([f"cartloader", "sge_combine_tiles",
                                f"--in-tiles {' '.join(updated_tiles)}",
                                f"--in-tile-minmax {tile_minmax}",
                                f"--out-dir {args.out_dir}",
                                f"--out-transcript {args.out_transcript}",
                                f"--out-minmax {args.out_minmax}",
                                f"--out-feature {args.out_feature}",
                                f"--colnames-count {' '.join(colnames_count)}",
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
        cmds.append(f'[ -f {os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")} ] && [ -f {os.path.join(args.out_dir, "feature.clean.tsv.gz")} ] && [ -f {os.path.join(args.out_dir, "coordinate_minmax.tsv")} ] && touch {sge_stitch_flag}')
        mm.add_target(sge_stitch_flag, prerequisities+[tile_minmax], cmds)

        # draw xy plot for visualization
        if args.sge_visual:
            cmds = cmd_separator([], f"Drawing XY plot")
            out_transcript=os.path.join(args.out_dir, args.out_transcript)
            out_xypng=os.path.join(args.out_dir, "xy.png")
            draw_cmd=f"{args.gzip} -dc {out_transcript} | tail -n +2 | cut -f 1,2 | {args.spatula} draw-xy --tsv /dev/stdin --out {out_xypng}"
            cmds.append(draw_cmd)
            mm.add_target(out_xypng, [sge_stitch_flag], cmds)

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
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])

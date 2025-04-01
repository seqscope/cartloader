import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess

from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, log_dataframe
from cartloader.utils.minimake import minimake

def hist_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Stitching tiles into one tif."
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_stitch.mk", help='Makefile name (default: sge_stitch.mk)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output Parameters")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of the input tiles in a specific format: <path>,<row>,<col>.")
    inout_params.add_argument("--output", type=str, help="Output path to store the merged file.")
    inout_params.add_argument("--in-offsets", type=str, help="Path to the input offsets file.")
    inout_params.add_argument("--colname-x-offset", type=str, default="x_offset", help="Column name for x offset in the input offsets file.")
    inout_params.add_argument("--colname-y-offset", type=str, default="y_offset", help="Column name for y offset in the input offsets file.")
    inout_params.add_argument('--units-per-um', type=float, default=1.0, help='Units per um in the input transcript tsv files (default: 1.0)')
    inout_params.add_argument('--convert-to-um', action='store_true', help='Convert output to um.')

    # env params
    # env_params = parser.add_argument_group("ENV Parameters", "Environment parameters for the tools")
    # env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p 4".')
    # env_params.add_argument('--spatula', type=str, default="spatula", help='Path to spatula binary.')
    args = parser.parse_args(_args)

    mm = minimake()

    updated_tiles = []
    prerequisities = []

    out_dir=os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    # read in_offsets
    offsets = pd.read_csv(args.in_offsets, sep="\t")
    offsets["row"] = offsets["row"].astype(int)
    offsets["col"] = offsets["col"].astype(int)
    # create dict between (row, col) to x_offset, y_offset
    idx2offsets = {}
    for i, row in offsets.iterrows():
        idx2offsets[(row['row'], row['col'])] = (row[args.colname_x_offset], row[args.colname_y_offset])

    # hist add offsets
    for in_tile in args.in_tiles:
        path, row, col = in_tile.split(",")
        x_offset, y_offset = idx2offsets[(int(row), int(col))]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Transcript file {path} does not exist.")
        # missing feature or minmax
        missing_ftr = feature is None or feature == "None"
        missing_minmax = minmax is None or minmax == "None"
        if missing_ftr or missing_minmax:
            in_dir = os.path.dirname(transcript)
            in_id = os.path.basename(transcript).replace('.tsv.gz', '').replace(".transcripts", "").replace(".transcript", "")
        if missing_ftr:
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
        if missing_minmax:
            cmds = cmd_separator([], f"Creating missing minmax file for {transcript}")
            out_minmax = os.path.join(in_dir, f"{in_id}.minmax.tsv")
            add_minmax_cmd =" ".join([f"cartloader", "sge_adds_on",
                                    f"--in-transcript {transcript}",
                                    "--add-minmax",
                                    f"--out-minmax {out_minmax}",
                                    f"--mu-scale {args.units_per_um}" if args.minmax_in_um else "",
                                ])
            cmds.append(add_minmax_cmd)
            mm.add_target(out_minmax, [transcript], cmds)
            minmax = out_minmax
        updated_tiles.append(",".join([transcript, feature, minmax, row, col]))
        prerequisities.extend([transcript, feature, minmax])

    
    # combine sge
    cmds = cmd_separator([], f"Combining SGEs")
    if "," in args.colnames_count:
        colnames_count = args.colnames_count.split(",")
    else:
        colnames_count = [args.colnames_count]
    combine_cmd=" ".join([f"cartloader", "combine_sges_by_layout",
                            f"--in-tiles {' '.join(updated_tiles)}",
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
                            '--minmax-in-um' if args.minmax_in_um else '',
                            '--convert-to-um' if args.convert_to_um else '',
                        ])
    cmds.append(combine_cmd)
    sge_stitch_flag = os.path.join(args.out_dir, "sge_stitch.done")
    cmds.append(f'[ -f {os.path.join(args.out_dir, "transcripts.unsorted.tsv.gz")} ] && [ -f {os.path.join(args.out_dir, "feature.clean.tsv.gz")} ] && [ -f {os.path.join(args.out_dir, "coordinate_minmax.tsv")} ] && touch {sge_stitch_flag}')
    mm.add_target(sge_stitch_flag, prerequisities, cmds)

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

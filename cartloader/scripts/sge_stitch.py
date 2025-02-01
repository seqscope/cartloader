import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess

from cartloader.utils.utils import cmd_separator, scheck_app, add_param_to_cmd, log_dataframe
from cartloader.utils.minimake import minimake

def sge_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Stitching tiles into one SGE."
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default="sge_stitch.mk", help='Makefile name (default: sge_stitch.mk)')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output Parameters")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of the input tiles in a specific format.")
    inout_params.add_argument("--out-dir", type=str, help="Output directory.")
    inout_params.add_argument("--colnames-count", type=str, nargs='*', help="Columns to sum (default: count).", default=['count'])
    inout_params.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    inout_params.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    inout_params.add_argument('--colname-x', type=str, default="X", help='X column name (default: X)')
    inout_params.add_argument('--colname-y', type=str, default="Y", help='Y column name (default: Y)')
    inout_params.add_argument('--units-per-um', type=float, default=1.0, help='Units per um in the input transcript tsv files (default: 1.0)')
    inout_params.add_argument('--minmax-in-um', action='store_true', help='Input minmax is in um while input transcript is based on unit.')
    inout_params.add_argument('--convert-to-um', action='store_true', help='Convert output to um.')
    args = parser.parse_args(_args)

    mm = minimake()

    updated_tiles = []
    prerequisities = []
    colnames_ftr=[
        args.colname_feature_name,
        args.colname_feature_id if args.colname_feature_id else "",
    ]
    
    # sge adds on
    for in_tile in args.in_tiles:
        transcript, feature, minmax, row, col = in_tile.split(",")
        if not os.path.exists(transcript):
            raise FileNotFoundError(f"Transcript file {transcript} does not exist.")
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
                                    f"--index-col {' '.join(colnames_ftr)}",
                                    f"--count-col {' '.join(args.colnames_count)}",
                                    f"--out-feature {out_ftr}",
                                    f"--mu-scale {args.units_per_um}" if args.minmax_in_um else "",
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
                                    f"--index-col {' '.join(colnames_ftr)}",
                                    f"--count-col {' '.join(args.colnames_count)}",
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
    combine_cmd=" ".join([f"cartloader", "combine_sges_by_layout",
                            f"--in-tiles {' '.join(updated_tiles)}",
                            f"--out-dir {args.out_dir}",
                            f"--colnames-count {' '.join(args.colnames_count)}",
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

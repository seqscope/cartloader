import sys
import os
import argparse
import inspect
import copy

from cartloader.image import register_georeference_stage, register_orientation_stage, configure_color_mode
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import scheck_actions, execute_makefile
from cartloader.utils.orient_helper import update_orient


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Apply georeferencing and/or orientation (rotate/flip) to raster images using minimake-backed stages.",
    )

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile but do not execute it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of parallel make jobs to run')
    run_params.add_argument('--makefn', type=str, default=None, help='Custom Makefile name (default: {out-prefix}.geoorient.mk)')

    action_params = parser.add_argument_group("Actions")
    action_params.add_argument('--georeference', action='store_true', default=False, help='Create a georeferenced GeoTIFF from the input image using supplied bounds')
    action_params.add_argument('--orient', action='store_true', default=False, help='Re-orient the image via rotation/flips')

    orient_params = parser.add_argument_group("Orientation Options")
    orient_params.add_argument('--rotate', type=str, choices=["90", "180", "270"], help='Rotate clockwise by 90/180/270 degrees (applied before flips)')
    orient_params.add_argument('--flip-vertical', action='store_true', default=False, help='Flip vertically (around X axis); applied after rotation')
    orient_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Flip horizontally (around Y axis); applied after rotation')

    inout_params = parser.add_argument_group("Input/Output")
    inout_params.add_argument('--in-img', type=str, required=True, help='Input raster (PNG/TIF)')
    inout_params.add_argument('--out-prefix', type=str, required=True, help='Prefix for outputs written by each stage')

    georef_params = parser.add_argument_group("Georeference Parameters")
    georef_params.add_argument('--georef-pixel-tsv', type=str, default=None, help='Bounds from *.pixel.sorted.tsv.gz (see run_ficture)')
    georef_params.add_argument('--georef-bounds-tsv', type=str, default=None, help='Bounds TSV with single line: <ulx>,<uly>,<lrx>,<lry>')
    georef_params.add_argument('--georef-bounds', type=str, default=None, help='Bounds string: "<ulx>,<uly>,<lrx>,<lry>"')
    georef_params.add_argument('--georef-detect', type=str, default=None, help='Auto-detect bounds (currently supports "ome")')
    georef_params.add_argument('--srs', type=str, default='EPSG:3857', help='Spatial reference system for georeferencing (default: EPSG:3857)')

    tool_params = parser.add_argument_group("Tool Paths")
    tool_params.add_argument('--gdal_translate', type=str, default='gdal_translate', help='Path to gdal_translate')
    tool_params.add_argument('--gdalinfo', type=str, default='gdalinfo', help='Path to gdalinfo (used for dimension extraction)')
    tool_params.add_argument('--gdalwarp', type=str, default='gdalwarp', help='Path to gdalwarp (used for reorientation)')

    color_params = parser.add_argument_group("Color Options")
    color_params.add_argument('--mono', action='store_true', default=False, help='Treat the input as single-band (mono) data')
    color_params.add_argument('--rgba', action='store_true', default=False, help='Treat the input as 4-band RGBA data')
    color_params.add_argument('--color-mode-record', type=str, default=None, help='Read color mode from a file (mutually exclusive with --mono/--rgba)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # infer orient action if rotation/flip flags are supplied
    if args.rotate or args.flip_vertical or args.flip_horizontal:
        args.orient = True

    scheck_actions(args, ['--georeference', '--orient', '--rotate', '--flip-vertical', '--flip-horizontal'], context='actions')

    if not args.georeference and not args.orient:
        parser.error('Nothing to do: enable at least --georeference or --orient (or provide rotate/flip flags)')

    return args


def image_geoorient(_args):
    args = parse_arguments(_args)

    assert os.path.exists(args.in_img), f"Input image not found: {args.in_img}"

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    mm = minimake()

    stage_args = copy.deepcopy(args)

    if stage_args.color_mode_record:
        print(f"--color-mode-record is provided. Read color mode from {stage_args.color_mode_record}")
    configure_color_mode(stage_args)

    current_img = stage_args.in_img

    if stage_args.georeference:
        current_img = register_georeference_stage(
            mm,
            stage_args,
            in_img=current_img,
            out_prefix=stage_args.out_prefix,
        )

    if stage_args.orient:
        stage_args.rotate, stage_args.flip_vertical, stage_args.flip_horizontal = update_orient(
            stage_args.rotate, stage_args.flip_vertical, stage_args.flip_horizontal
        )
        current_img = register_orientation_stage(
            mm,
            stage_args,
            src_tif=current_img,
            out_prefix=stage_args.out_prefix,
        )

    make_f = (
        os.path.join(out_dir, stage_args.makefn)
        if stage_args.makefn is not None and out_dir
        else f"{stage_args.out_prefix}.geoorient.mk"
    )
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=stage_args.dry_run, restart=stage_args.restart, n_jobs=stage_args.n_jobs)


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name)
    func(sys.argv[1:])

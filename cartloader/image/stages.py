"""Reusable building blocks for image conversion pipelines."""

from __future__ import annotations

import gzip
import os
from dataclasses import dataclass
from typing import Dict, Optional

import shlex

import tifffile

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app
from cartloader.utils.orient_helper import (
    get_orientation_suffix,
    orient2axisorder,
)

_REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_GDAL_GET_SIZE_SCRIPT = os.path.join(_REPO_DIR, "utils", "gdal_get_size.sh")


@dataclass
class Png2PmtilesResult:
    """Outputs produced when registering the PNG→PMTiles workflow."""

    georef_tif: str
    oriented_tif: str
    final_tif: str
    mbtile_flag: Optional[str]
    mbtile_path: Optional[str]
    pmtiles_path: Optional[str]


def configure_color_mode(args) -> None:
    """Normalise mono/rgba flags using a recorded color-mode file when needed."""

    color_mode_record = getattr(args, "color_mode_record", None)
    needs_color_mode = (
        getattr(args, "flip_vertical", False)
        or getattr(args, "flip_horizontal", False)
        or getattr(args, "rotate", None) is not None
        or getattr(args, "geotif2mbtiles", False)
    )

    setattr(args, "color_mode_pending", False)

    if not needs_color_mode or not color_mode_record:
        return

    if getattr(args, "mono", False) or getattr(args, "rgba", False):
        raise ValueError("--color-mode-record cannot be combined with --mono or --rgba")

    allow_missing = getattr(args, "allow_missing_color_mode_record", False)

    if not os.path.exists(color_mode_record):
        if allow_missing:
            args.color_mode_pending = True
            return
        raise FileNotFoundError(f"File not found: {color_mode_record} (--color-mode-record)")

    with open(color_mode_record, "r", encoding="utf-8") as handle:
        color_mode = handle.readline().strip().lower()

    if color_mode not in {"rgb", "rgba", "mono"}:
        raise ValueError(f"Invalid color mode '{color_mode}' in file: {color_mode_record}")

    if color_mode == "rgba":
        args.rgba = True
        args.mono = False
    elif color_mode == "mono":
        args.rgba = False
        args.mono = True
    else:  # rgb
        args.rgba = False
        args.mono = False


def _resolve_bounds_from_args(args, *, in_img: str) -> Optional[Dict[str, float]]:

    # only one should be provided and indicate that current georef_detect only supports ome.
    georef_bounds=getattr(args, "georef_bounds", None)
    georef_bounds_tsv=getattr(args, "georef_bounds_tsv", None)
    georef_pixel_tsv=getattr(args, "georef_pixel_tsv", None)
    georef_detect=getattr(args, "georef_detect", None)

    georef_inputs = [georef_bounds is not None, georef_bounds_tsv is not None, georef_pixel_tsv is not None, georef_detect is not None]
    assert sum(georef_inputs) == 1, (
        f"georeferencing bounds not found or ambiguous. Provide exactly one of --georef-bounds, --georef-bound-tsv, --georef-pixel-tsv, or georef_detect"
    )

    if georef_bounds:
        ulx, uly, lrx, lry = georef_bounds.split(",")
        return {
            "ulx": float(ulx),
            "uly": float(uly),
            "lrx": float(lrx),
            "lry": float(lry),
        }

    if georef_bounds_tsv:
        with open(georef_bounds_tsv, "r", encoding="utf-8") as handle:
            line = handle.readline()
        ulx, uly, lrx, lry = line.strip().lower().split(",")
        return {
            "ulx": float(ulx),
            "uly": float(uly),
            "lrx": float(lrx),
            "lry": float(lry),
        }

    if georef_pixel_tsv:
        with gzip.open(georef_pixel_tsv, "rt", encoding="utf-8") as handle:
            for _ in range(3):
                line = handle.readline()
        ann2val = {
            token.split("=")[0]: token.split("=")[1]
            for token in line.strip().replace("##", "").split(";")
        }
        ulx = float(ann2val["OFFSET_X"])
        uly = float(ann2val["OFFSET_Y"])
        lrx = float(ann2val["SIZE_X"]) + 1 + ulx
        lry = float(ann2val["SIZE_Y"]) + 1 + uly
        return {"ulx": ulx, "uly": uly, "lrx": lrx, "lry": lry}

    if georef_detect:
        detect_mode = georef_detect.lower()
        if detect_mode != "ome":
            raise ValueError(f"Unsupported --georef-detect mode: {georef_detect}")
        with tifffile.TiffFile(in_img) as tif:
            meta = tifffile.xml2dict(tif.ome_metadata)["OME"]["Image"]["Pixels"]
            physical_size_x = meta["PhysicalSizeX"]
            physical_size_y = meta["PhysicalSizeY"]
            size_x = meta["SizeX"]
            size_y = meta["SizeY"]
            px_size_unit = meta.get("PhysicalSizeXUnit", "um")
            if px_size_unit not in {"um", "µm"}:
                raise ValueError(
                    f"Physical size unit is not supported for --georef-detect=ome: {px_size_unit}"
                )
            ulx = float(meta.get("OffsetX", 0))
            uly = float(meta.get("OffsetY", 0))
            lrx = ulx + float(physical_size_x) * int(size_x)
            lry = uly + float(physical_size_y) * int(size_y)
            return {"ulx": ulx, "uly": uly, "lrx": lrx, "lry": lry}

    return None


def register_georeference_stage(
    mm: minimake,
    args,
    *,
    in_img: str,
    out_prefix: str,
) -> str:
    scheck_app(args.gdal_translate)

    bounds = _resolve_bounds_from_args(args, in_img=in_img)
    if bounds is None:
        raise ValueError(
            "Georeferencing requested but no bounds provided via --georef-*, or --georef-detect"
        )

    georef_f = f"{out_prefix}.georef.tif"
    cmds = cmd_separator([], f"Geo-referencing {in_img} to {georef_f}")
    ullr = "{ulx} {uly} {lrx} {lry}".format(**bounds)
    cmds.append(
        " ".join(
            [
                args.gdal_translate,
                "-of GTiff",
                f"-a_srs {args.srs}",
                f"-a_ullr {ullr}",
                in_img,
                georef_f,
            ]
        )
    )
    mm.add_target(georef_f, [in_img], cmds)
    return georef_f


def register_dimension_stage(mm: minimake, *, src_tif: str, gdalinfo: str) -> str:
    dim_f = os.path.splitext(src_tif)[0] + ".dim.tsv"
    cmds = cmd_separator([], f"Extract dimensions from: {src_tif}")
    cmds.append(f"{_GDAL_GET_SIZE_SCRIPT} {src_tif} {dim_f} {gdalinfo}")
    mm.add_target(dim_f, [src_tif], cmds)
    return dim_f


def register_orientation_stage(
    mm: minimake,
    args,
    *,
    src_tif: str,
    out_prefix: str,
    dim_f: Optional[str] = None,
) -> str:
    if not (
        getattr(args, "flip_vertical", False)
        or getattr(args, "flip_horizontal", False)
        or getattr(args, "rotate", None) is not None
    ):
        return src_tif

    axis_order = orient2axisorder.get(
        (args.rotate, args.flip_vertical, args.flip_horizontal)
    )
    if axis_order is None:
        raise ValueError("Invalid combination of rotation and flip options.")

    dim_f = dim_f or register_dimension_stage(mm, src_tif=src_tif, gdalinfo=args.gdalinfo)

    ort_suffix = get_orientation_suffix(args.rotate, args.flip_vertical, args.flip_horizontal)
    ort_f = f"{out_prefix}.{ort_suffix}.tif"

    if axis_order.startswith("1") or axis_order.startswith("-1"):
        out_dim = "$WIDTH $HEIGHT"
    else:
        out_dim = "$HEIGHT $WIDTH"

    msg = " ".join(
        [
            "Orientate",
            src_tif,
            "(vertical flip)" if args.flip_vertical else "",
            "(horizontal flip)" if args.flip_horizontal else "",
            f"rotate {args.rotate} deg" if args.rotate else "",
        ]
    ).strip()

    cmds = cmd_separator([], msg)
    cmds.append(f"WIDTH=$(awk '/WIDTH/' {dim_f}|cut -f 2) && \\")
    cmds.append(f"HEIGHT=$(awk '/HEIGHT/' {dim_f}|cut -f 2) && \\")
    band_args = "-b 1"
    if getattr(args, "rgba", False):
        band_args = "-b 1 -b 2 -b 3 -b 4"
    elif not getattr(args, "mono", False):
        band_args = "-b 1 -b 2 -b 3"

    gdalwarp_bin = getattr(args, "gdalwarp", "gdalwarp")

    cmd = " ".join(
        [
            gdalwarp_bin,
            f'"{src_tif}"',
            f'"{ort_f}"',
            band_args,
            f"-ct \"+proj=pipeline +step +proj=axisswap +order={axis_order}\"",
            "-overwrite",
            "-ts",
            out_dim,
        ]
    )
    cmds.append(cmd)
    mm.add_target(ort_f, [src_tif, dim_f], cmds)
    return ort_f


def create_mbtile_flag(mbtile_flag: str, mbtile_f: str, partial_db: str, journal_db: str) -> None:
    os.makedirs(os.path.dirname(mbtile_flag), exist_ok=True)
    conditions_met = (
        not os.path.exists(mbtile_flag)
        and os.path.exists(mbtile_f)
        and not os.path.exists(partial_db)
        and not os.path.exists(journal_db)
    )
    if conditions_met:
        with open(mbtile_flag, "a", encoding="utf-8"):
            pass


def register_geotif2mbtiles_stage(
    mm: minimake,
    args,
    *,
    src_tif: str,
    out_prefix: str,
) -> Optional[Dict[str, str]]:
    if not getattr(args, "geotif2mbtiles", False):
        return None

    scheck_app(args.gdal_translate)

    mbtile_f = f"{out_prefix}.pmtiles.mbtiles"
    mbtile_flag = f"{mbtile_f}.done"
    partial_db = mbtile_f.replace(".mbtiles", ".partial_tiles.db")
    journal_db = f"{mbtile_f}-journal"

    create_mbtile_flag(mbtile_flag, mbtile_f, partial_db, journal_db)

    cmds = cmd_separator([], f"Converting from geotif to mbtiles: {src_tif}")
    cleanup_cmd = (
        f"if [ -f {journal_db} ] || [ -f {partial_db} ] ; then echo 'Warning: Cleaning up incomplete previous conversion...' ; rm -f {mbtile_f} {journal_db} {partial_db} ; fi"
    )
    cmds.append(cleanup_cmd)

    color_mode_record = getattr(args, "color_mode_record", None)
    resolve_at_runtime = bool(color_mode_record) and getattr(args, "color_mode_pending", False)

    if resolve_at_runtime:
        quoted_record = shlex.quote(color_mode_record)
        translate_cmd = "; ".join(
            [
                f"COLOR_MODE=$(head -n 1 {quoted_record} | tr '[:upper:]' '[:lower:]')",
                "case \"$COLOR_MODE\" in",
                "rgba) BAND_ARGS=\"-b 1 -b 2 -b 3 -b 4\"; SCALE_FLAG=\"\" ;;",
                "mono) BAND_ARGS=\"-b 1\"; SCALE_FLAG=\"-scale\" ;;",
                "rgb) BAND_ARGS=\"-b 1 -b 2 -b 3\"; SCALE_FLAG=\"\" ;;",
                "*) echo \"Invalid color mode '$COLOR_MODE' in {color_mode_record}\" >&2; exit 1 ;;",
                "esac",
                args.gdal_translate,
                "$BAND_ARGS",
                "-strict",
                f"-co \"ZOOM_LEVEL_STRATEGY=UPPER\"",
                f"-co \"RESAMPLING={args.resample}\"",
                f"-co \"BLOCKSIZE={args.blocksize}\"",
                "-ot Byte",
                "$SCALE_FLAG",
                "-of mbtiles",
                f"-a_srs {args.srs}",
                src_tif,
                mbtile_f,
            ]
        )
    else:
        band_args = "-b 1"
        if getattr(args, "rgba", False):
            band_args = "-b 1 -b 2 -b 3 -b 4"
        elif not getattr(args, "mono", False):
            band_args = "-b 1 -b 2 -b 3"

        translate_cmd = " ".join(
            [
                args.gdal_translate,
                band_args,
                "-strict",
                f"-co \"ZOOM_LEVEL_STRATEGY=UPPER\"",
                f"-co \"RESAMPLING={args.resample}\"",
                f"-co \"BLOCKSIZE={args.blocksize}\"",
                "-ot Byte",
                "-scale" if getattr(args, "mono", False) else "",
                "-of mbtiles",
                f"-a_srs {args.srs}",
                src_tif,
                mbtile_f,
            ]
        )

    cmds.append(translate_cmd)
    validation_cmd = (
        f" [ -f {mbtile_f} ]  && [ ! -f {journal_db} ] && [ ! -f {partial_db} ] && touch {mbtile_flag}"
    )
    cmds.append(validation_cmd)
    prereqs = [src_tif]
    if resolve_at_runtime:
        prereqs.append(color_mode_record)
    mm.add_target(mbtile_flag, prereqs, cmds)

    return {
        "mbtile_flag": mbtile_flag,
        "mbtile_f": mbtile_f,
        "partial_db": partial_db,
        "journal_db": journal_db,
        "mbtile_resampled": f"{out_prefix}.pmtiles.{args.resample}.mbtiles",
    }


def register_mbtiles2pmtiles_stage(
    mm: minimake,
    args,
    *,
    mbtile_flag: str,
    mbtile_f: str,
    mbtile_resampled: str,
    out_prefix: str,
) -> Optional[str]:
    if not getattr(args, "mbtiles2pmtiles", False):
        return None

    scheck_app(args.pmtiles)
    scheck_app(args.gdaladdo)

    pmtiles_f = f"{out_prefix}.pmtiles"

    cmds = cmd_separator([], f"Resampling mbtiles and converting to pmtiles: {mbtile_f}")
    cmds.append(f"cp {mbtile_f} {mbtile_resampled}")
    cmds.append(
        f"'{args.gdaladdo}' {mbtile_resampled} -r {args.resample} 2 4 8 16 32 64 128 256 512 1024 2048 4096"
    )
    cmds.append(f"'{args.pmtiles}' convert --force {mbtile_resampled} {pmtiles_f}")
    cmds.append(f" [ -f {pmtiles_f} ] && rm {mbtile_resampled}")
    mm.add_target(pmtiles_f, [mbtile_flag], cmds)

    return pmtiles_f


def register_png2pmtiles_pipeline(
    mm: minimake,
    args,
    *,
    in_img: Optional[str] = None,
    out_prefix: Optional[str] = None,
) -> Png2PmtilesResult:
    src_img = in_img if in_img is not None else args.in_img
    prefix = out_prefix if out_prefix is not None else args.out_prefix

    georef_f = src_img
    if getattr(args, "georeference", False):
        georef_f = register_georeference_stage(mm, args, in_img=src_img, out_prefix=prefix)

    oriented_f = register_orientation_stage(mm, args, src_tif=georef_f, out_prefix=prefix)

    mbtile_info = register_geotif2mbtiles_stage(mm, args, src_tif=oriented_f, out_prefix=prefix)

    pmtiles_f = None
    mbtile_flag = None
    mbtile_path = None
    if mbtile_info:
        mbtile_flag = mbtile_info["mbtile_flag"]
        mbtile_path = mbtile_info["mbtile_f"]
        pmtiles_f = register_mbtiles2pmtiles_stage(
            mm,
            args,
            mbtile_flag=mbtile_flag,
            mbtile_f=mbtile_path,
            mbtile_resampled=mbtile_info["mbtile_resampled"],
            out_prefix=prefix,
        )

    final_tif = oriented_f if oriented_f else georef_f

    return Png2PmtilesResult(
        georef_tif=georef_f,
        oriented_tif=oriented_f,
        final_tif=final_tif,
        mbtile_flag=mbtile_flag,
        mbtile_path=mbtile_path,
        pmtiles_path=pmtiles_f,
    )

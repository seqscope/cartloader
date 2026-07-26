"""Reusable building blocks for image conversion pipelines."""

from __future__ import annotations

import gzip
import os
from dataclasses import dataclass
from typing import Dict, Optional
import json
import subprocess

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

import json
import subprocess

def _needs_rgb_expansion_cli(image_path: str, args) -> bool:
    """Checks if an image needs '-expand rgb' using the gdalinfo CLI tool."""
    scheck_app(args.gdalinfo)

    try:
        # Run gdalinfo and capture the JSON output
        result = subprocess.run(
            [args.gdalinfo, "-json", image_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        info_json = json.loads(result.stdout)
        bands = info_json.get('bands', [])
        
        if bands:
            return bands[0].get('colorInterpretation', '') == 'Palette'
            
    except (subprocess.CalledProcessError, json.JSONDecodeError, KeyError, IndexError):
        # Handle cases where the file doesn't exist or gdalinfo fails
        pass
        
    return False

# --- Example Usage ---
# file_path = "cartostore/visiumhd-public-dataset-collection/.../rep1.t12_f48_pixel.png"
# if needs_rgb_expansion(file_path):
#     print("Use: -expand rgb")
# else:
#     print("Do not use: -expand rgb")
def _read_png_size(image_path: str) -> tuple[int, int]:
    """Return (width, height) of a PNG by parsing its IHDR chunk (no gdal/PIL)."""
    with open(image_path, "rb") as handle:
        header = handle.read(24)
    # 8-byte PNG signature, then 4-byte length + 'IHDR' type, then width/height (big-endian uint32)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError(f"Not a valid PNG file: {image_path}")
    width = int.from_bytes(header[16:20], "big")
    height = int.from_bytes(header[20:24], "big")
    return width, height


def _get_image_size(image_path: str) -> tuple[int, int]:
    """Return (width, height) of a TIF or PNG image without relying on gdal."""
    ext = os.path.splitext(image_path)[1].lower()
    if ext in {".tif", ".tiff"}:
        with tifffile.TiffFile(image_path) as tif:
            page = tif.pages[0]
            return int(page.imagewidth), int(page.imagelength)
    if ext == ".png":
        return _read_png_size(image_path)
    raise ValueError(
        f"Cannot determine image size for --georef-plain (unsupported extension '{ext}'): {image_path}"
    )

# Conversion factors from an OME PhysicalSize unit to microns (um).
_OME_LENGTH_TO_UM = {
    "nm": 1e-3,
    "µm": 1.0,
    "um": 1.0,
    "micron": 1.0,
    "microns": 1.0,
    "mm": 1e3,
    "cm": 1e4,
    "m": 1e6,
}


def _ome_physical_size_um(pixels: dict, axis: str) -> float:
    """Return PhysicalSize<axis> from an OME Pixels dict, converted to microns."""
    size = pixels.get(f"PhysicalSize{axis}")
    if size is None:
        raise ValueError(
            f"OME metadata has no PhysicalSize{axis}; cannot use --georef-detect=ome"
        )
    unit = pixels.get(f"PhysicalSize{axis}Unit", "µm")
    scale = _OME_LENGTH_TO_UM.get(unit)
    if scale is None:
        raise ValueError(
            f"Unsupported OME PhysicalSize{axis}Unit '{unit}' for --georef-detect=ome "
            f"(supported: {', '.join(sorted(_OME_LENGTH_TO_UM))})"
        )
    return float(size) * scale


def _select_ome_pixels(ome: dict, *, in_img: str) -> dict:
    """Pick the OME Pixels block matching the input image among one or more Images.

    An OME-TIFF may embed multiple <Image> elements (e.g. a full-resolution
    pyramid plus derived masks/downsamples), so ``xml2dict`` yields a list. Match
    the Image whose pixel dimensions equal the actual raster, falling back to the
    first Image when none matches.
    """
    images = ome["Image"]
    if isinstance(images, dict):
        images = [images]

    actual_x, actual_y = _get_image_size(in_img)
    for image in images:
        pixels = image["Pixels"]
        if int(pixels["SizeX"]) == actual_x and int(pixels["SizeY"]) == actual_y:
            return pixels
    return images[0]["Pixels"]


def _resolve_bounds_from_args(args, *, in_img: str) -> Optional[Dict[str, float]]:

    # only one should be provided and indicate that current georef_detect only supports ome.
    georef_bounds=getattr(args, "georef_bounds", None)
    georef_bounds_tsv=getattr(args, "georef_bounds_tsv", None)
    georef_pixel_tsv=getattr(args, "georef_pixel_tsv", None)
    georef_detect=getattr(args, "georef_detect", None)
    georef_plain=getattr(args, "georef_plain", False)

    georef_inputs = [georef_bounds is not None, georef_bounds_tsv is not None, georef_pixel_tsv is not None, georef_detect is not None, georef_plain]
    assert sum(georef_inputs) == 1, (
        f"georeferencing bounds not found or ambiguous. Provide exactly one of --georef-bounds, --georef-bound-tsv, --georef-pixel-tsv, --georef-detect, or --georef-plain"
    )

    if georef_bounds:
        ulx, uly, lrx, lry = georef_bounds.split(",")
        return {
            "ulx": float(ulx) + args.georef_offset_x,
            "uly": float(uly) + args.georef_offset_y,
            "lrx": float(lrx) + args.georef_offset_x,
            "lry": float(lry) + args.georef_offset_y,
        }

    if georef_bounds_tsv:
        with open(georef_bounds_tsv, "r", encoding="utf-8") as handle:
            line = handle.readline()
        ulx, uly, lrx, lry = line.strip().lower().split(",")
        return {
            "ulx": float(ulx) + args.georef_offset_x,
            "uly": float(uly) + args.georef_offset_y,
            "lrx": float(lrx) + args.georef_offset_x,
            "lry": float(lry) + args.georef_offset_y,
        }

    if georef_pixel_tsv:
        with gzip.open(georef_pixel_tsv, "rt", encoding="utf-8") as handle:
            for _ in range(3):
                line = handle.readline()
        ann2val = {
            token.split("=")[0]: token.split("=")[1]
            for token in line.strip().replace("##", "").split(";")
        }
        ulx = float(ann2val["OFFSET_X"]) + args.georef_offset_x
        uly = float(ann2val["OFFSET_Y"]) + args.georef_offset_y
        lrx = float(ann2val["SIZE_X"]) + 1 + ulx + args.georef_offset_x
        lry = float(ann2val["SIZE_Y"]) + 1 + uly + args.georef_offset_y
        return {"ulx": ulx, "uly": uly, "lrx": lrx, "lry": lry}

    if georef_detect:
        detect_mode = georef_detect.lower()
        if detect_mode != "ome":
            raise ValueError(f"Unsupported --georef-detect mode: {georef_detect}")
        with tifffile.TiffFile(in_img) as tif:
            ome = tifffile.xml2dict(tif.ome_metadata)["OME"]
            meta = _select_ome_pixels(ome, in_img=in_img)
            size_x = int(meta["SizeX"])
            size_y = int(meta["SizeY"])
            um_per_x = _ome_physical_size_um(meta, "X")
            um_per_y = _ome_physical_size_um(meta, "Y")
            ulx = float(meta.get("OffsetX", 0)) + args.georef_offset_x
            uly = float(meta.get("OffsetY", 0)) + args.georef_offset_y
            lrx = ulx + um_per_x * size_x
            lry = uly + um_per_y * size_y
            return {"ulx": ulx, "uly": uly, "lrx": lrx, "lry": lry}

    if georef_plain:
        um_per_pixel = float(getattr(args, "um_per_pixel", 1.0))
        size_x, size_y = _get_image_size(in_img)
        ulx = args.georef_offset_x
        uly = args.georef_offset_y
        lrx = ulx + size_x * um_per_pixel
        lry = uly + size_y * um_per_pixel
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

    #print(f"args.mono = {args.mono}, args.rgba = {args.rgba}, {in_img.endswith('.png')} {getattr(args, 'mono', False)}")

    georef_f = f"{out_prefix}.georef.tif"

    # --georef-detect gtiff: the input already carries the geotransform we want (e.g. a
    # Seq-Scope H&E TIF registered upstream, whose corner coordinates are already the
    # transcript um extent), but no CRS. That combination is not "already georeferenced"
    # as far as the tilers are concerned: geotiff2pmtiles resolves the max zoom from the
    # CRS, so without one it settles on 0, writes an empty PMTiles and still exits 0.
    # Stamp --srs on and leave the geotransform alone -- no -a_ullr, so the corner
    # coordinates are not round-tripped through a decimal rendering of themselves.
    if str(getattr(args, "georef_detect", "") or "").lower() == "gtiff":
        cmds = cmd_separator([], f"Assigning {args.srs} to the existing geotransform of {in_img}")
        cmds.append(
            " ".join([args.gdal_translate, "-of GTiff", f"-a_srs {args.srs}", in_img, georef_f])
        )
        mm.add_target(georef_f, [in_img], cmds)
        return georef_f

    bounds = _resolve_bounds_from_args(args, in_img=in_img)
    if bounds is None:
        raise ValueError(
            "Georeferencing requested but no bounds provided via --georef-*, or --georef-detect"
        )

    cmds = cmd_separator([], f"Geo-referencing {in_img} to {georef_f}")
    ullr = "{ulx} {uly} {lrx} {lry}".format(**bounds)
    ## check if rgb expansion is needed
    if in_img.endswith(".png") and _needs_rgb_expansion_cli(in_img, args):
        if getattr(args, "rgba", False):
            expand_str = "-expand rgba"
        elif getattr(args, "mono", False):
            expand_str = "-expand gray"
        else:
            expand_str = "-expand rgb"
    else:
        expand_str = ""
    cmds.append(
        " ".join(
            [
                args.gdal_translate,
                "-of GTiff",
                f"-a_srs {args.srs}",
                f"-a_ullr {ullr}",
                expand_str,
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

def register_gdalwarp_stage(
    mm: minimake,
    args,
    *,
    src_tif: str,
    out_prefix: str,
) -> Optional[str]:
    
    gdalwarp_bin = getattr(args, "gdalwarp", "gdalwarp")
    warped_f = f"{out_prefix}.warped.tif"

    scheck_app(gdalwarp_bin)

    cmds = cmd_separator([], f"Performing pixel-level operations on GeoTIFF: {src_tif}")
    cmds.append(f"'{gdalwarp_bin}' -r bilinear -of GTiff {src_tif} {warped_f}")
    mm.add_target(warped_f, [src_tif], cmds)

    return warped_f

def register_geotiff2pmtiles_stage(
    mm: minimake,
    args,
    *,
    src_tif: str,
    out_prefix: str,
) -> Optional[str]:
    if not getattr(args, "geotiff2pmtiles", False):
        return None

    scheck_app(args.geotiff2pmtiles)

    pmtiles_f = f"{out_prefix}.pmtiles"

    # Rescale controls for 16-bit imagery: geotiff2pmtiles errors on 16-bit input
    # unless given an explicit --rescale-range (e.g. some Stereo-seq H&E TIFs).
    rescale = ""
    if getattr(args, "rescale", None):
        rescale += f"--rescale {args.rescale} "
    if getattr(args, "rescale_range", None):
        rescale += f"--rescale-range {args.rescale_range} "

    cmds = cmd_separator([], f"Converting from geotiff to pmtiles: {src_tif}")
    # geotiff2pmtiles derives its zoom range from the CRS, and on an input that has none it
    # resolves the max zoom to 0 -- below --min-zoom -- so it writes an empty PMTiles and
    # still exits 0. Catch that here rather than letting an empty layer reach the catalog.
    cmds.append(
        f"if ! '{getattr(args, 'gdalinfo', 'gdalinfo')}' {src_tif} | grep -q 'Coordinate System is'; then "
        f"echo 'ERROR: {src_tif} has no CRS, so geotiff2pmtiles would write an empty PMTiles. "
        f"Re-run with --georeference plus a --georef-* mode (use --georef-detect gtiff to keep "
        f"an existing geotransform and only assign the CRS).' >&2; exit 1; fi"
    )
    cmds.append(f"'{args.geotiff2pmtiles}' --format {args.tile_format} --min-zoom {args.min_zoom} " + (f"--max-zoom {args.max_zoom} " if args.max_zoom is not None else "") + rescale + f"{src_tif} {pmtiles_f}")
    mm.add_target(pmtiles_f, [src_tif], cmds)

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

    if args.method == "gdal": ## use gdal-based pipeline for performing pmtiles conversion
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
    elif args.method == "geotiff2pmtiles": ## use geotiff2pmtiles for direct conversion to pmtiles without mbtiles intermediate
        georef_f = src_img
        if getattr(args, "georeference", False):
            georef_f = register_georeference_stage(mm, args, in_img=src_img, out_prefix=prefix)

        gdalwarp_f = register_gdalwarp_stage(mm, args, src_tif=georef_f, out_prefix=prefix)

        pmtiles_f = register_geotiff2pmtiles_stage(mm, args, src_tif=gdalwarp_f, out_prefix=prefix)

        return Png2PmtilesResult(
            georef_tif=georef_f,
            oriented_tif=gdalwarp_f,
            final_tif=gdalwarp_f,
            mbtile_flag=None,
            mbtile_path=None,
            pmtiles_path=pmtiles_f
        )
    # elif args.method == "ficture2": ## use ficture2 for direct conversion to pmtiles without mbtiles intermediate
    #     pmtiles_f = register_ficture2_stage(mm, args, src_img=src_img, out_prefix=prefix)

    #     return Png2PmtilesResult(
    #         georef_tif=None,
    #         oriented_tif=None,
    #         final_tif=None,
    #         mbtile_flag=None,
    #         mbtile_path=None,
    #         pmtiles_path=pmtiles_f
    #     )
    else:
        raise ValueError(f"Unsupported method: {args.method}")

import sys, os, argparse, inspect, json, gzip, io, mmap, logging
from multiprocessing import Pool

from cartloader.utils.utils import create_custom_logger

# Raster tile types this tool can flip (pmtiles TileType values -> Pillow format name).
RASTER_FORMATS = {2: "PNG", 3: "JPEG", 4: "WEBP"}


def parse_arguments(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
                                     description="""
                                     Flip the sign of the y coordinate of a raster PMTiles (y -> -y), i.e. mirror it
                                     about y = 0. Use it for an image built with the GeoTIFF 'north-up' convention
                                     (pixel row r at y = -r), so that it lines up with cartloader layers, which put
                                     row r at y = +r. Tiles move from row y to row 2^z-1-y and each tile is flipped
                                     vertically; Web Mercator is symmetric about y = 0, so the mirror is exact. The
                                     header and metadata bounds are negated to match. Each tile is decoded and
                                     re-encoded in its original format.
                                     """)
    parser.add_argument('--in-pmtiles', type=str, required=True, help='Input raster PMTiles (PNG, JPEG or WebP tiles)')
    parser.add_argument('--out-pmtiles', type=str, required=True, help='Output PMTiles with the y coordinate negated')
    parser.add_argument('--quality', type=int, default=90, help='Encoding quality for WebP/JPEG tiles (default: 90; the input tiles are usually already lossy, so keep it at or above theirs)')
    parser.add_argument('--lossless', action='store_true', default=False, help='Re-encode WebP tiles losslessly (larger output; ignored for PNG, which is always lossless, and JPEG)')
    parser.add_argument('--threads', type=int, default=4, help='Number of processes for decoding and re-encoding tiles (default: 4)')
    parser.add_argument('--log', action='store_true', default=False, help='Write a log file next to the output')
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def walk_entries(buf, header):
    """Yield (tile_id, data_offset, data_length) for every addressed tile, without reading tile data."""
    from pmtiles.tile import deserialize_directory

    def walk(offset, length):
        for e in deserialize_directory(buf[offset:offset + length]):
            if e.run_length > 0:
                for i in range(e.run_length):
                    yield e.tile_id + i, header["tile_data_offset"] + e.offset, e.length
            else:
                yield from walk(header["leaf_directory_offset"] + e.offset, e.length)

    yield from walk(header["root_offset"], header["root_length"])


# Per-worker encoding settings, set once by the pool initializer.
_ENC = {}


def _init_worker(fmt, tile_gzip, quality, lossless):
    _ENC.update(fmt=fmt, tile_gzip=tile_gzip, quality=quality, lossless=lossless)


def flip_tile(data):
    from PIL import Image
    if _ENC["tile_gzip"]:
        data = gzip.decompress(data)
    im = Image.open(io.BytesIO(data))
    im.load()
    im = im.transpose(Image.Transpose.FLIP_TOP_BOTTOM)
    out = io.BytesIO()
    fmt = _ENC["fmt"]
    if fmt == "WEBP":
        if _ENC["lossless"]:
            im.save(out, format="WEBP", lossless=True)
        else:
            im.save(out, format="WEBP", quality=_ENC["quality"])
    elif fmt == "JPEG":
        im.save(out, format="JPEG", quality=_ENC["quality"])
    else:
        im.save(out, format="PNG", optimize=False)
    data = out.getvalue()
    return gzip.compress(data) if _ENC["tile_gzip"] else data


def negate_bounds_str(s):
    """'minlon,minlat,maxlon,maxlat' -> same with the latitudes negated (and swapped)."""
    v = [float(x) for x in str(s).split(",")]
    return f"{v[0]:.6f},{-v[3] + 0.0:.6f},{v[2]:.6f},{-v[1] + 0.0:.6f}"   # + 0.0 turns -0.0 into 0.0


def negate_center_str(s):
    v = str(s).split(",")
    v[1] = f"{-float(v[1]) + 0.0:.6f}"
    return ",".join(v)


def image_pmtiles_flip_y(_args):
    args = parse_arguments(_args)
    try:
        from pmtiles.tile import deserialize_header, zxy_to_tileid, tileid_to_zxy, Compression
        from pmtiles.writer import Writer
    except ImportError:
        sys.exit("ERROR: image_pmtiles_flip_y needs the 'pmtiles' Python package (pip install pmtiles).")

    out_dir = os.path.dirname(os.path.abspath(args.out_pmtiles))
    os.makedirs(out_dir, exist_ok=True)
    logger = create_custom_logger(__name__, args.out_pmtiles + ".log" if args.log else None)
    if os.path.abspath(args.in_pmtiles) == os.path.abspath(args.out_pmtiles):
        sys.exit("ERROR: --out-pmtiles must differ from --in-pmtiles.")

    with open(args.in_pmtiles, "rb") as fin:
        buf = mmap.mmap(fin.fileno(), 0, access=mmap.ACCESS_READ)
        header = deserialize_header(buf[0:127])
        ttype = header["tile_type"].value if hasattr(header["tile_type"], "value") else int(header["tile_type"])
        if ttype not in RASTER_FORMATS:
            sys.exit(f"ERROR: {args.in_pmtiles} holds tile type {header['tile_type']}; only raster PNG/JPEG/WebP "
                     f"tiles can be flipped (vector tiles would need their geometry rewritten).")
        tcomp = header["tile_compression"]
        if tcomp not in (Compression.NONE, Compression.UNKNOWN, Compression.GZIP):
            sys.exit(f"ERROR: unsupported tile compression {tcomp} in {args.in_pmtiles}.")
        meta_raw = buf[header["metadata_offset"]:header["metadata_offset"] + header["metadata_length"]]
        if header["internal_compression"] == Compression.GZIP:
            meta_raw = gzip.decompress(meta_raw)
        metadata = json.loads(meta_raw) if meta_raw else {}

        # New tile id of every addressed tile: row y -> 2^z-1-y. Writing in new-id order keeps
        # the archive clustered.
        entries = []
        for tid, off, ln in walk_entries(buf, header):
            z, x, y = tileid_to_zxy(tid)
            entries.append((zxy_to_tileid(z, x, (1 << z) - 1 - y), off, ln))
        entries.sort()
        logger.info(f"{args.in_pmtiles}: {len(entries)} tiles, zoom {header['min_zoom']}-{header['max_zoom']}, "
                    f"{RASTER_FORMATS[ttype]} tiles")

        # Flip each distinct source tile once (run-length/deduplicated tiles share a data offset).
        uniq = sorted({(off, ln) for _, off, ln in entries})
        flipped = {}
        with Pool(processes=max(1, args.threads), initializer=_init_worker,
                  initargs=(RASTER_FORMATS[ttype], tcomp == Compression.GZIP, args.quality, args.lossless)) as pool:
            for key, data in zip(uniq, pool.imap(flip_tile, (bytes(buf[o:o + l]) for o, l in uniq), chunksize=64)):
                flipped[key] = data

        tmp_out = args.out_pmtiles + ".tmp"
        with open(tmp_out, "wb") as fout:
            writer = Writer(fout)
            for new_id, off, ln in entries:
                writer.write_tile(new_id, flipped[(off, ln)])
            # Negate the latitude bounds (and swap min/max) and the center latitude.
            header["min_lat_e7"], header["max_lat_e7"] = -header["max_lat_e7"], -header["min_lat_e7"]
            header["center_lat_e7"] = -header["center_lat_e7"]
            if "bounds" in metadata:
                metadata["bounds"] = negate_bounds_str(metadata["bounds"])
            if "center" in metadata:
                metadata["center"] = negate_center_str(metadata["center"])
            note = "y coordinate negated (mirrored about y = 0) by cartloader image_pmtiles_flip_y"
            metadata["description"] = (metadata.get("description", "") + "\n\n" + note).strip()
            writer.finalize(header, metadata)
        buf.close()
    os.replace(tmp_out, args.out_pmtiles)
    logger.info(f"Wrote {args.out_pmtiles} (lat bounds {header['min_lat_e7'] / 1e7}..{header['max_lat_e7'] / 1e7})")


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

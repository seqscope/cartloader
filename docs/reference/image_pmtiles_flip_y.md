# PMTiles Y Flip

## Overview

`image_pmtiles_flip_y` negates the y coordinate of a raster PMTiles (y → −y), mirroring it about y = 0.

CartLoader places pixel row *r* of an image at y = +*r* µm, the same frame as the transcripts. An image tiled from a standard north-up GeoTIFF whose top-left corner is at (0, 0) puts row *r* at y = −*r* instead. Its y bounds are negative, so it appears just below the data and upside down relative to it. This tool moves such an image into CartLoader's frame without rebuilding it from the source image.

`run_together` runs it for a prebuilt image layer with `flip_y=true` (see [image fields](./run_together_images.md#fields)), for example:

```bash
--image type=hne,kind=prebuilt,flip_y=true,source=/path/to/hne.pmtiles
```

## How it works

- Each tile moves from row *y* to row 2<sup>*z*</sup> − 1 − *y* at the same zoom *z* and column, and its pixels are flipped top to bottom. Web Mercator is symmetric about y = 0, so the mirror is exact; no tile is resampled.
- The header and metadata bounds and the center latitude are negated, and a note is added to the metadata description.
- Each tile is decoded and re-encoded in its original format (PNG, JPEG or WebP). PNG stays lossless. WebP and JPEG are re-encoded at `--quality`; use `--lossless` to keep WebP tiles exact at the cost of a larger file.
- Only raster tiles are supported; vector tiles would need their geometry rewritten.

It needs the `pmtiles` Python package (`pip install pmtiles`).

## Usage

```bash
cartloader image_pmtiles_flip_y \
    --in-pmtiles hne.northup.pmtiles \
    --out-pmtiles hne.pmtiles \
    --threads 8
```

## Parameters

- `--in-pmtiles` (str, required): input raster PMTiles (PNG, JPEG or WebP tiles).
- `--out-pmtiles` (str, required): output PMTiles; must differ from the input.
- `--quality` (int): WebP/JPEG re-encoding quality (default: 90). The input tiles are usually already lossy, so keep it at or above their quality.
- `--lossless` (flag): re-encode WebP tiles losslessly.
- `--threads` (int): processes used to decode and re-encode tiles (default: 4).
- `--log` (flag): also write `<out-pmtiles>.log`.

#!/usr/bin/env python3
"""
annotate_spatial_factors_with_ai.py

Combines three capabilities:
  1. render_pmtiles   — renders a PMTiles file into a B&W background PNG
  2. plot_factor_overlay — overlays per-factor expression points onto the background
  3. annotate_bulk_de_with_ai — calls a generative AI API to annotate each factor
                                using marker genes + the per-factor spatial image

For each factor the script will:
  - Render the background once (reused for all factors)
  - Generate a per-factor overlay PNG showing only that factor's expression
  - Include that image alongside the marker-gene prompt sent to the AI
  - Write the final alias table to --out

Usage:
    python3 annotate_spatial_factors_with_ai.py \\
        --pmtiles   /path/to/sge-mono-dark.pmtiles \\
        --zoom      15 \\
        --tsv       factors.tsv.gz \\
        --de        bulk_de.tsv \\
        --tissue    "mouse brain" \\
        --organism  mouse \\
        --api-type  claude \\
        --out       annotations.tsv \\
        --tmp-dir   /tmp/spatial_annot

Dependencies:
    pip install pmtiles Pillow numpy pandas requests google-genai

Environment variables (set before running):
    ANTHROPIC_API_KEY   (for --api-type claude)
    OPENAI_API_KEY      (for --api-type openai)
    GEMINI_API_KEY      (for --api-type google)
"""

import argparse
import base64
import colorsys
import gzip
import inspect
import io
import json
import logging
import math
import os
import re
import sys
import time
import zlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests
from PIL import Image, ImageDraw

# ── optional Google genai (only needed for --api-type google) ─────────────────
try:
    from google import genai
    from google.genai import types as genai_types
    _HAVE_GENAI = True
except ImportError:
    _HAVE_GENAI = False

# =============================================================================
# Constants / defaults
# =============================================================================

OPENAI_API_KEY_ENV    = "OPENAI_API_KEY"
GEMINI_API_KEY_ENV    = "GEMINI_API_KEY"
ANTHROPIC_API_KEY_ENV = "ANTHROPIC_API_KEY"

OPENAI_MODEL_ENV    = "OPENAI_MODEL"
GOOGLE_MODEL_ENV    = "GOOGLE_MODEL"
ANTHROPIC_MODEL_ENV = "ANTHROPIC_MODEL"

DEFAULT_OPENAI_MODEL    = "gpt-5-mini"
DEFAULT_GOOGLE_MODEL    = "gemini-3-flash-preview"
DEFAULT_ANTHROPIC_MODEL = "claude-opus-4-6"

# EPSG:3857 / PMTiles geometry constants — see pmpoint/pmt_utils.{h,cpp}.
# Keep these in sync with the C++ side: `epsg3857totilecoord` assumes a
# 256-px tile, and the half-bound is the Web-Mercator extent in meters.
EPSG_3857_BOUND   = 20037508.3428         # = 6378137.0 * pi
PMTILES_TILE_PX   = 256.0                 # native tile size used by pmt_utils
EARTH_RADIUS_M    = 6378137.0


def _wgs84_to_epsg3857(lon: float, lat: float) -> Tuple[float, float]:
    """WGS84 lon/lat (degrees) -> EPSG:3857 (meters)."""
    x = math.radians(lon) * EARTH_RADIUS_M
    y = EARTH_RADIUS_M * math.log(math.tan(math.pi / 4.0 + math.radians(lat) / 2.0))
    return x, y


def _epsg3857_to_world_pixel(ix, iy, zoom: int):
    """
    Mirror pmt_utils::epsg3857totilecoord, returning a single global pixel
    coordinate on the PMTiles 256-px tile grid: world_pixel = tile * 256 +
    local_pixel.  Works on scalars or numpy arrays.
    """
    world_size_px = (1 << zoom) * PMTILES_TILE_PX
    px_per_meter  = world_size_px / (2.0 * EPSG_3857_BOUND)
    return (ix + EPSG_3857_BOUND) * px_per_meter, \
           (EPSG_3857_BOUND - iy) * px_per_meter


# =============================================================================
# ── SECTION 1: PMTiles rendering (from render_pmtiles.py) ────────────────────
# =============================================================================

def _decompress_tile(data: bytes, compression) -> bytes:
    from pmtiles.tile import Compression
    if compression == Compression.GZIP:
        return zlib.decompress(data, 16 + zlib.MAX_WBITS)
    elif compression == Compression.ZSTD:
        try:
            import zstandard as zstd
            return zstd.ZstdDecompressor().decompress(data)
        except ImportError:
            sys.exit("Install zstandard: pip install zstandard")
    elif compression == Compression.BROTLI:
        try:
            import brotli
            return brotli.decompress(data)
        except ImportError:
            sys.exit("Install brotli: pip install brotli")
    return data


def _tiles_at_zoom(source, header, zoom: int):
    from pmtiles.reader import traverse
    for (z, x, y), data in traverse(source, header,
                                     header['root_offset'],
                                     header['root_length']):
        if z == zoom:
            yield z, x, y, data


def render_background(pmtiles_path: str, zoom: int) -> Tuple[Image.Image, dict]:
    """
    Render a PMTiles file at the given zoom level and return
    (PIL Image, transform dict).  The transform dict contains the
    parameters needed to map TSV coordinates to image pixels.
    """
    try:
        from pmtiles.reader import MmapSource, Reader
        from pmtiles.tile import TileType
    except ImportError:
        sys.exit("Install pmtiles: pip install pmtiles")

    f      = open(pmtiles_path, 'rb')
    source = MmapSource(f)
    reader = Reader(source)
    header = reader.header()

    tt = header['tile_type']
    raster_types = {TileType.PNG, TileType.JPEG, TileType.WEBP, TileType.AVIF}

    tiles = list(_tiles_at_zoom(source, header, zoom))
    if not tiles:
        sys.exit(f"No tiles at zoom {zoom} in {pmtiles_path}. "
                 f"Available: {header['min_zoom']}–{header['max_zoom']}")

    xs = [t[1] for t in tiles]
    ys = [t[2] for t in tiles]
    tx_min, tx_max = min(xs), max(xs)
    ty_min, ty_max = min(ys), max(ys)

    sample_data = _decompress_tile(tiles[0][3], header['tile_compression'])
    sample_img  = Image.open(io.BytesIO(sample_data))
    tile_w, tile_h = sample_img.size

    cols = tx_max - tx_min + 1
    rows = ty_max - ty_min + 1

    if tt in raster_types:
        canvas = Image.new('RGBA', (cols * tile_w, rows * tile_h), (0, 0, 0, 255))
        for _, tx, ty, raw in tiles:
            data = _decompress_tile(raw, header['tile_compression'])
            img  = Image.open(io.BytesIO(data)).convert('RGBA')
            canvas.paste(img, ((tx - tx_min) * tile_w, (ty - ty_min) * tile_h))
    elif tt.name == 'MVT':
        canvas = _render_mvt_tiles(tiles, header, tx_min, ty_min, tile_w, tile_h)
    else:
        sys.exit(f"Unsupported tile_type: {tt}")

    f.close()

    # Build transform: pixel coords of the bbox corners. Project lon/lat
    # bounds through EPSG:3857 -> world-pixel -> canvas-pixel using the
    # same formulas as pmt_utils::epsg3857totilecoord; that way the bbox
    # and the per-point overlay share one projection.
    lon_min = header['min_lon_e7'] / 1e7
    lon_max = header['max_lon_e7'] / 1e7
    lat_min = header['min_lat_e7'] / 1e7
    lat_max = header['max_lat_e7'] / 1e7

    canvas_scale_x = tile_w / PMTILES_TILE_PX
    canvas_scale_y = tile_h / PMTILES_TILE_PX

    def lon_to_px(lon):
        ex, _ = _wgs84_to_epsg3857(lon, 0.0)
        wpx, _ = _epsg3857_to_world_pixel(ex, 0.0, zoom)
        return wpx * canvas_scale_x - tx_min * tile_w

    def lat_to_py(lat):
        _, ey = _wgs84_to_epsg3857(0.0, lat)
        _, wpy = _epsg3857_to_world_pixel(0.0, ey, zoom)
        return wpy * canvas_scale_y - ty_min * tile_h

    transform = {
        "zoom": zoom, "tx_min": tx_min, "ty_min": ty_min,
        "tile_w": tile_w, "tile_h": tile_h,
        "lon_min": lon_min, "lon_max": lon_max,
        "lat_min": lat_min, "lat_max": lat_max,
        "bbox_px": {
            "left":   lon_to_px(lon_min),
            "right":  lon_to_px(lon_max),
            "top":    lat_to_py(lat_max),
            "bottom": lat_to_py(lat_min),
        }
    }

    return canvas, transform


def _render_mvt_tiles(tiles, header, tx_min, ty_min, tile_w, tile_h) -> Image.Image:
    try:
        import mapbox_vector_tile
    except ImportError:
        sys.exit("Install mapbox-vector-tile: pip install mapbox-vector-tile")
    xs = [t[1] for t in tiles]; ys = [t[2] for t in tiles]
    cols = max(xs) - tx_min + 1; rows = max(ys) - ty_min + 1
    mvt_extent = 4096; scale = tile_w / mvt_extent
    canvas = Image.new('RGBA', (cols * tile_w, rows * tile_h), (0, 0, 0, 255))
    draw   = ImageDraw.Draw(canvas)
    for _, tx, ty, raw in tiles:
        data = _decompress_tile(raw, header['tile_compression'])
        tile = mapbox_vector_tile.decode(data)
        ox = (tx - tx_min) * tile_w; oy = (ty - ty_min) * tile_h
        for layer in tile.values():
            for feature in layer['features']:
                geom = feature['geometry']
                cl = geom['coordinates']
                if geom['type'] == 'Polygon': cl = [cl]
                elif geom['type'] != 'MultiPolygon': continue
                for polygon in cl:
                    for ring in polygon:
                        pts = [(ox + c[0]*scale, oy + (mvt_extent-c[1])*scale) for c in ring]
                        if len(pts) >= 3:
                            draw.polygon(pts, fill=(200,200,200,255), outline=(120,120,120,180))
    return canvas


# =============================================================================
# ── SECTION 2: Factor overlay (from plot_factor_overlay.py) ──────────────────
# =============================================================================

def _make_palette(n: int) -> List[Tuple[int, int, int]]:
    colours = []
    for i in range(n):
        h = i / n
        l = 0.55 if i % 2 == 0 else 0.65
        r, g, b = colorsys.hls_to_rgb(h, l, 0.85)
        colours.append((int(r*255), int(g*255), int(b*255)))
    return colours


# Heatmap colour scale: dark blue (low prob) → cyan → green → yellow → red (high prob).
# This is a perceptually intuitive "cool-to-warm" scale that reads clearly on a
# dark greyscale background.
_HEATMAP_STOPS: List[Tuple[float, Tuple[int, int, int]]] = [
    (0.00, (  0,   0, 180)),   # dark blue   — very low probability
    (0.20, (  0, 180, 220)),   # cyan        — low-medium
    (0.40, ( 50, 200,  50)),   # green       — medium
    (0.60, (240, 220,   0)),   # yellow      — medium-high
    (1.00, (220,  30,  30)),   # red         — high probability
]

def _prob_to_color(p: float) -> Tuple[int, int, int]:
    """Map a probability in [0, 1] to an RGB colour via the heatmap stops."""
    p = max(0.0, min(1.0, float(p)))
    stops = _HEATMAP_STOPS
    # Find the two bracketing stops
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if p <= t1:
            frac = (p - t0) / (t1 - t0) if t1 > t0 else 0.0
            return tuple(int(c0[j] + frac * (c1[j] - c0[j])) for j in range(3))
    return stops[-1][1]


def _deduplicate(df: pd.DataFrame, radius: float) -> pd.DataFrame:
    if radius <= 0 or len(df) == 0:
        return df
    sort_col = 'topP' if 'topP' in df.columns else df.columns[0]
    df = df.sort_values(sort_col, ascending=False).reset_index(drop=True)
    occupied: set = set()
    keep = []
    for _, row in df.iterrows():
        key = (int(row['X'] / radius), int(row['Y'] / radius))
        if key not in occupied:
            occupied.add(key); keep.append(True)
        else:
            keep.append(False)
    return df[keep].reset_index(drop=True)


def _make_coord_transform(df: pd.DataFrame, transform: dict):
    """
    Returns a vectorised function (xs, ys) -> (px_array, py_array).

    Treats input (xs, ys) as EPSG:3857 (Web-Mercator meters) — the same
    coordinate space pmpoint feeds to pmt_utils::epsg3857totilecoord when
    building the PMTiles file — and projects directly onto the rendered
    canvas. The C++ function works on a 256-px tile; here we rescale the
    local pixel offset to the actual rendered tile_w / tile_h so the same
    formula handles raster and MVT canvases that use different tile sizes.
    """
    del df  # data-driven scaling no longer required
    zoom   = transform['zoom']
    tx_min = transform['tx_min']
    ty_min = transform['ty_min']
    tile_w = transform['tile_w']
    tile_h = transform['tile_h']

    px_per_meter   = (1 << zoom) * PMTILES_TILE_PX / (2.0 * EPSG_3857_BOUND)
    canvas_scale_x = tile_w / PMTILES_TILE_PX
    canvas_scale_y = tile_h / PMTILES_TILE_PX
    canvas_origin_x = tx_min * tile_w
    canvas_origin_y = ty_min * tile_h

    def fn(xs: np.ndarray, ys: np.ndarray):
        wpx = (xs + EPSG_3857_BOUND) * px_per_meter
        wpy = (EPSG_3857_BOUND - ys) * px_per_meter
        px = wpx * canvas_scale_x - canvas_origin_x
        py = wpy * canvas_scale_y - canvas_origin_y
        return px.astype(int), py.astype(int)

    return fn


def render_factor_overlay(
    bg: Image.Image,
    df: pd.DataFrame,
    transform: dict,
    factor_id: int,
    colour: Tuple[int, int, int],
    point_size: int = 3,
    alpha: float = 0.85,
    use_heatmap: bool = False,
) -> Image.Image:
    """
    Composite a single factor's points over the background image.

    When use_heatmap=True the topP column must be present: each point is
    coloured by its probability using a blue→cyan→green→yellow→red scale
    (low → high).  The scale is fixed to [0, 1] so it is consistent across
    all factors.

    When use_heatmap=False every point is drawn in the single `colour` passed
    in (uniform mode).

    Returns a new RGBA PIL Image.
    """
    coord_fn = _make_coord_transform(df, transform)
    img_w, img_h = bg.size
    overlay   = Image.new('RGBA', (img_w, img_h), (0, 0, 0, 0))
    draw      = ImageDraw.Draw(overlay)
    alpha_int = max(0, min(255, int(alpha * 255)))
    r         = max(1, point_size)

    subset = df[df['topK'] == factor_id]
    if len(subset) > 0:
        px, py = coord_fn(subset['X'].values.astype(float),
                          subset['Y'].values.astype(float))
        if use_heatmap:
            probs = subset['topP'].values
            for x, y, p in zip(px.tolist(), py.tolist(), probs.tolist()):
                c = _prob_to_color(p)
                draw.ellipse([x-r, y-r, x+r, y+r], fill=(*c, alpha_int))
        else:
            rgba = (*colour, alpha_int)
            for x, y in zip(px.tolist(), py.tolist()):
                draw.ellipse([x-r, y-r, x+r, y+r], fill=rgba)

    return Image.alpha_composite(bg, overlay)


# =============================================================================
# ── SECTION 3: AI annotation helpers (from annotate_bulk_de_with_ai.py) ──────
# =============================================================================

def _normalize_alias(raw: str) -> str:
    if raw is None:
        return "Unknown"
    s = raw.strip()
    if s.startswith("{") and s.endswith("}"):
        try:
            obj = json.loads(s)
            if isinstance(obj, dict) and "alias" in obj:
                s = str(obj["alias"]).strip()
        except Exception:
            pass
    s = s.strip("`\"'")
    s = s.splitlines()[0].strip()
    s = re.split(r"[.;]", s, maxsplit=1)[0].strip()
    if " " in s:
        parts = re.split(r"\s+", s)
        s = "".join(p[:1].upper() + p[1:] for p in parts if p)
    s = re.sub(r"[^A-Za-z0-9+\-/]", "", s)
    return s if s else "Unknown"


def _extract_json_alias(text: str) -> str:
    if not text:
        return "Unknown"
    t = text.strip()
    if t.startswith("{") and t.endswith("}"):
        try:
            obj = json.loads(t)
            if isinstance(obj, dict) and "alias" in obj:
                return _normalize_alias(str(obj["alias"]))
        except Exception:
            pass
    m = re.search(r"\{.*\}", t, flags=re.DOTALL)
    if m:
        try:
            obj = json.loads(m.group(0))
            if isinstance(obj, dict) and "alias" in obj:
                return _normalize_alias(str(obj["alias"]))
        except Exception:
            pass
    return _normalize_alias(t)


def _make_prompt(tissue: str, organism: str, genes: List[str],
                 factor_id: int, has_image: bool,
                 use_heatmap: bool = False) -> str:
    gene_list = ", ".join(genes)
    image_block = ""
    if has_image:
        if use_heatmap:
            colour_desc = (
                "Points are colored by assignment probability using a heatmap scale: "
                "dark blue = low probability, cyan = low-medium, green = medium, "
                "yellow = medium-high, red = high probability. "
                "The scale is fixed to [0, 1] and is identical across all factors. "
                "Regions with predominantly red/yellow points indicate high-confidence "
                "assignments; blue/cyan regions are lower confidence."
            )
        else:
            colour_desc = (
                "Points are shown in a single uniform colour; all cells assigned to "
                "this factor are coloured identically regardless of confidence."
            )
        image_block = (
            "\nAdditional input:\n"
            f"- A spatial expression image is attached. The background is a black-and-white "
            f"image of the {tissue} tissue section. Coloured points overlaid on it show the "
            f"distribution of factor {factor_id}. {colour_desc} "
            f"Use the spatial pattern ONLY as supporting evidence to disambiguate "
            f"the cell type suggested by the marker genes.\n"
            "- Prioritize marker genes first; use spatial pattern second.\n"
        )
    return (
        "You are annotating latent factors from bulk differential expression.\n"
        "Task: Identify the single most likely cell type represented by these "
        "ordered, top marker genes.\n\n"
        f"Organism: {organism}\n"
        f"Tissue: {tissue}\n"
        f"Factor index: {factor_id}\n"
        f"Top marker genes (comma-separated): {gene_list}\n"
        f"{image_block}\n"
        "Return ONLY a JSON object with this exact schema:\n"
        "{\"alias\": \"UpperCamelCaseCellType\"}\n\n"
        "Rules:\n"
        "- alias must be terse and singular.\n"
        "- Use informative shorthand when appropriate "
        "(e.g., CD4+T, CD8+T, NKCell, BCell, PlasmaCell).\n"
        "- Avoid long phrases, parentheses, or multi-sentence explanations.\n"
        "- Do not include any extra keys besides 'alias'.\n"
    )


def _image_to_b64(img: Image.Image, max_dim: int = 1568) -> Tuple[str, str]:
    # Downsample to fit within max_dim on the longest side before encoding.
    # Anthropic recommends images no larger than 1568px on the long edge.
    w, h = img.size
    if max(w, h) > max_dim:
        scale = max_dim / max(w, h)
        img = img.resize((int(w * scale), int(h * scale)), Image.LANCZOS)
    buf = io.BytesIO()
    img.convert('RGB').save(buf, format='PNG')
    return "image/png", base64.b64encode(buf.getvalue()).decode('utf-8')


def _load_thumbnail_b64(path: Optional[str]) -> Optional[Tuple[str, str]]:
    if not path:
        return None
    if not os.path.exists(path):
        raise FileNotFoundError(f"--thumbnail not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    mime = "image/png" if ext == ".png" else "image/jpeg"
    with open(path, "rb") as f:
        return mime, base64.b64encode(f.read()).decode("utf-8")


def call_openai(prompt, thumbnail, model_name, request_timeout, max_retries):
    api_key = os.environ.get(OPENAI_API_KEY_ENV, "")
    model   = os.environ.get(OPENAI_MODEL_ENV, DEFAULT_OPENAI_MODEL) \
              if model_name is None else model_name
    url     = "https://api.openai.com/v1/responses"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    content = [{"type": "input_text", "text": prompt}]
    if thumbnail:
        mime, b64 = thumbnail
        content.append({"type": "input_image",
                        "image_url": f"data:{mime};base64,{b64}"})
    payload = {"model": model,
               "input": [{"role": "user", "content": content}],
               "reasoning": {"effort": "medium"}}
    data = {}
    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
            r.raise_for_status(); data = r.json(); break
        except Exception as e:
            if attempt == max_retries - 1: raise
            time.sleep(2 ** attempt + 5)
    if isinstance(data, dict) and data.get("output_text"):
        return str(data["output_text"])
    try:
        texts = []
        for item in data.get("output", []):
            for c in item.get("content", []):
                if c.get("type") == "output_text":
                    texts.append(c["text"])
        return "\n".join(texts).strip()
    except Exception:
        return ""


def call_google(prompt, thumbnail, model_name, request_timeout, max_retries):
    if not _HAVE_GENAI:
        sys.exit("Install google-genai: pip install google-genai")
    client   = genai.Client()
    model_id = os.environ.get(GOOGLE_MODEL_ENV, DEFAULT_GOOGLE_MODEL) \
               if model_name is None else model_name
    config   = genai_types.GenerateContentConfig(
        thinking_config=genai_types.ThinkingConfig(
            thinking_level=genai_types.ThinkingLevel.HIGH))
    if thumbnail is None:
        contents = prompt
    else:
        mime, b64 = thumbnail
        contents = [prompt,
                    genai_types.Part.from_bytes(data=base64.b64decode(b64),
                                                mime_type=mime)]
    for attempt in range(max_retries):
        try:
            resp = client.models.generate_content(model=model_id,
                                                  contents=contents,
                                                  config=config)
            return resp.text.strip()
        except Exception as e:
            if attempt == max_retries - 1: raise
            time.sleep(2 ** attempt + 5)
    return ""


def call_claude(prompt, thumbnail, model_name, request_timeout, max_retries):
    api_key = os.environ.get(ANTHROPIC_API_KEY_ENV, "")
    model   = os.environ.get(ANTHROPIC_MODEL_ENV, DEFAULT_ANTHROPIC_MODEL) \
              if model_name is None else model_name
    url     = "https://api.anthropic.com/v1/messages"
    headers = {"x-api-key": api_key, "anthropic-version": "2023-06-01",
               "content-type": "application/json"}
    content_blocks = []
    if thumbnail:
        mime, b64 = thumbnail
        content_blocks.append({"type": "image",
                                "source": {"type": "base64",
                                           "media_type": mime,
                                           "data": b64}})
    content_blocks.append({"type": "text", "text": prompt})
    payload = {"model": model, "max_tokens": 256, "temperature": 0.2,
               "messages": [{"role": "user", "content": content_blocks}]}
    data = {}
    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
            r.raise_for_status(); data = r.json(); break
        except Exception as e:
            if attempt == max_retries - 1: raise
            time.sleep(2 ** attempt + 5)
    try:
        blocks = data.get("content", [])
        return "\n".join(b.get("text", "") for b in blocks
                         if b.get("type") == "text").strip()
    except Exception:
        return ""


# =============================================================================
# ── SECTION 4: DE reading / gene ranking (from annotate_bulk_de_with_ai.py) ──
# =============================================================================

def read_bulk_de(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"gene", "factor"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input missing required columns: {sorted(required)}")
    return df


def top_genes_per_factor(df, primary_rank, secondary_rank, top_n):
    for col in (primary_rank, secondary_rank):
        if col not in df.columns:
            raise ValueError(f"Ranking column '{col}' not found. "
                             f"Available: {list(df.columns)}")
    try:
        factors = sorted(df["factor"].astype(int).unique().tolist())
        df = df.copy(); df["factor"] = df["factor"].astype(int)
    except Exception:
        factors = sorted(df["factor"].unique().tolist())

    out: Dict[int, List[str]] = {}
    for f in factors:
        sub  = df[df["factor"] == f].copy()
        prim = sub.sort_values(primary_rank, ascending=False)["gene"] \
                  .astype(str).head(top_n).tolist()
        sec  = sub.sort_values(secondary_rank, ascending=False)["gene"] \
                  .astype(str).head(top_n).tolist()
        seen: set = set(); merged: List[str] = []
        for i in range(top_n):
            for g in ([prim[i]] if i < len(prim) else []) + \
                     ([sec[i]]  if i < len(sec)  else []):
                if g not in seen:
                    merged.append(g); seen.add(g)
        out[int(f)] = merged
    return out


# =============================================================================
# ── SECTION 5: Main orchestration ────────────────────────────────────────────
# =============================================================================

def annotate_spatial_factors_with_ai(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── PMTiles / spatial rendering ───────────────────────────────────────
    spatial = parser.add_argument_group(
        "Spatial rendering",
        "Parameters for rendering the background PMTiles image and factor overlay."
    )
    spatial.add_argument('--pmtiles', type=str, default=None,
                         help='PMTiles file for background spatial image. '
                              'If omitted, no spatial image is sent to the AI.')
    spatial.add_argument('--zoom', type=int, default=15,
                         help='Zoom level for PMTiles rendering (default: 15)')
    spatial.add_argument('--tsv', type=str, default=None,
                         help='Factor TSV (or .tsv.gz) with columns X, Y, topK, topP. '
                              'Required when --pmtiles is provided.')
    spatial.add_argument('--point-size', type=int, default=3,
                         help='Marker radius in pixels for factor overlay (default: 3)')
    spatial.add_argument('--alpha', type=float, default=0.85,
                         help='Marker opacity 0–1 (default: 0.85)')
    spatial.add_argument('--dedup-radius', type=float, default=0.0,
                         help='Dedup radius in TSV units to remove tile-boundary '
                              'duplicates (default: 0.0)')
    spatial.add_argument('--save-overlays', action='store_true',
                         help='Save per-factor overlay PNGs to --tmp-dir '
                              '(useful for inspection)')
    spatial.add_argument('--tmp-dir', type=str, default=None,
                         help='Directory for temporary overlay PNGs. '
                              'Defaults to a sibling of --out.')

    # ── DE / gene annotation ──────────────────────────────────────────────
    inout = parser.add_argument_group("Input/Output", "DE file and output.")
    inout.add_argument('--de', required=True, type=str,
                       help='TSV with bulk DE test results (columns: gene, factor, …)')
    inout.add_argument('--out', required=True, type=str,
                       help='Output TSV with columns: index, alias')
    inout.add_argument('--tissue', required=True, type=str,
                       help='Tissue name used in the prompt')
    inout.add_argument('--organism', default="human",
                       help='Organism (default: human)')
    inout.add_argument('--api-type', required=True,
                       choices=['openai', 'google', 'claude'],
                       help='Generative AI backend')
    inout.add_argument('--thumbnail', type=str, default=None,
                       help='Optional static thumbnail for ALL factors '
                            '(overridden per-factor when --pmtiles is given)')

    aux = parser.add_argument_group("Auxiliary", "Ranking / model options.")
    aux.add_argument('--primary-rank',   default="Chi2")
    aux.add_argument('--secondary-rank', default="FoldChange")
    aux.add_argument('--top-n',          type=int, default=10)
    aux.add_argument('--model-name',     type=str, default=None)
    aux.add_argument('--request-timeout', type=int, default=60)
    aux.add_argument('--max-retries',    type=int, default=3)

    args = parser.parse_args(_args)

    # ── logging ───────────────────────────────────────────────────────────
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s - %(levelname)s - %(message)s]",
        datefmt="[%Y-%m-%d %H:%M:%S]",
    )
    logger = logging.getLogger(__name__)

    # ── validate spatial args ─────────────────────────────────────────────
    use_spatial = args.pmtiles is not None
    if use_spatial and not args.tsv:
        parser.error("--tsv is required when --pmtiles is provided")

    # ── render background once ────────────────────────────────────────────
    bg_img   = None
    transform = None
    df_tsv   = None
    palette  = None

    if use_spatial:
        logger.info(f"Rendering background from {args.pmtiles} at zoom {args.zoom}…")
        bg_img, transform = render_background(args.pmtiles, args.zoom)
        logger.info(f"  Background size: {bg_img.width}×{bg_img.height} px")

        logger.info(f"Loading factor TSV: {args.tsv}")
        opener = gzip.open if args.tsv.endswith('.gz') else open
        with opener(args.tsv, 'rt') as fh:
            df_tsv = pd.read_csv(fh, sep='\t')
        logger.info(f"  {len(df_tsv):,} points")

        df_tsv = _deduplicate(df_tsv, args.dedup_radius)
        logger.info(f"  After dedup: {len(df_tsv):,} points")

        # Detect whether per-point probabilities are available
        use_heatmap = 'topP' in df_tsv.columns
        if use_heatmap:
            logger.info("  topP column detected — using heatmap colour scale "
                        "(blue=low → red=high probability)")
        else:
            logger.info("  No topP column — using uniform colour for all points")

        # Uniform fallback colour (red)
        uniform_colour = (220, 30, 30)

        # tmp dir for overlays
        tmp_dir = Path(args.tmp_dir) if args.tmp_dir \
                  else Path(args.out).parent / "spatial_overlays"
        tmp_dir.mkdir(parents=True, exist_ok=True)

    # ── read DE + top genes ───────────────────────────────────────────────
    logger.info(f"Reading DE results from {args.de}")
    df_de = read_bulk_de(args.de)
    logger.info(f"Collecting top {args.top_n} genes per factor…")
    factor2genes = top_genes_per_factor(
        df_de,
        primary_rank=args.primary_rank,
        secondary_rank=args.secondary_rank,
        top_n=args.top_n,
    )
    logger.info(f"Found {len(factor2genes)} factors")

    # ── fallback static thumbnail ─────────────────────────────────────────
    static_thumbnail = _load_thumbnail_b64(args.thumbnail)

    # ── annotate factor by factor ─────────────────────────────────────────
    results: List[List] = []
    alias2cnts: Dict[str, int] = {}

    for idx in sorted(factor2genes.keys()):
        genes = factor2genes[idx]

        # Build per-factor spatial image if PMTiles provided
        if use_spatial and df_tsv is not None:
            overlay_img = render_factor_overlay(
                bg=bg_img,
                df=df_tsv,
                transform=transform,
                factor_id=idx,
                colour=uniform_colour,
                point_size=args.point_size,
                alpha=args.alpha,
                use_heatmap=use_heatmap,
            )
            thumbnail = _image_to_b64(overlay_img)

            if args.save_overlays:
                overlay_path = tmp_dir / f"factor_{idx:03d}.png"
                overlay_img.save(str(overlay_path))
                logger.info(f"  Saved overlay: {overlay_path}")
        else:
            thumbnail = static_thumbnail

        has_image = thumbnail is not None
        prompt = _make_prompt(
            tissue=args.tissue,
            organism=args.organism,
            genes=genes,
            factor_id=idx,
            has_image=has_image,
            use_heatmap=use_heatmap if use_spatial else False,
        )

        logger.info(f"Annotating factor {idx} ({len(genes)} genes, "
                    f"image={'yes' if has_image else 'no'}) with {args.api_type}…")

        dispatch = {"openai": call_openai,
                    "google": call_google,
                    "claude": call_claude}
        text  = dispatch[args.api_type](
            prompt, thumbnail, args.model_name,
            args.request_timeout, args.max_retries,
        )
        alias = _extract_json_alias(text)
        logger.info(f"  Factor {idx} -> {alias}")
        results.append([idx, alias])
        alias2cnts[alias] = alias2cnts.get(alias, 0) + 1

    # ── resolve duplicate aliases ─────────────────────────────────────────
    logger.info("Resolving duplicate aliases…")
    alias2iter: Dict[str, int] = {}
    for i in range(len(results)):
        idx, alias = results[i]
        if alias2cnts[alias] > 1:
            alias2iter[alias] = alias2iter.get(alias, 0) + 1
            results[i][1] = f"{alias}_{alias2iter[alias]}"

    # ── write output ──────────────────────────────────────────────────────
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_out = pd.DataFrame(results, columns=["index", "alias"]).sort_values("index")
    df_out.to_csv(str(out_path), sep="\t", index=False)
    logger.info(f"Written {len(df_out)} annotations to {out_path}")


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name, None)
    if func is None:
        func = annotate_spatial_factors_with_ai
    func(sys.argv[1:])

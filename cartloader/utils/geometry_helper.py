"""Geometry helpers shared across cell-segmentation import and cell decode.

Default 10x Visium HD segmentation output carries cell polygons but no cell
centroids, so the centroid must be derived from the boundary geometry. This
module is the single source of truth for that derivation; both
``import_visiumhd_cell`` and ``run_ficture2_multi_cells`` call it.
"""
import json
import math

from shapely.geometry import shape
from shapely.affinity import scale as shapely_scale


def rescale_geometry(geom, units_per_um):
    """Rescale geometry into microns when needed.

    ``units_per_um`` is the number of coordinate units per micron in the input
    (e.g. pixels-per-um, i.e. ``1 / microns_per_pixel`` from a Visium HD scale
    JSON). A value of 1 (or None) means the coordinates are already in microns
    and the geometry is returned unchanged.
    """
    if units_per_um is None or math.isclose(units_per_um, 1.0):
        return geom
    if units_per_um == 0:
        raise ValueError("units_per_um must be non-zero")
    scale_factor = 1.0 / units_per_um
    return shapely_scale(geom, xfact=scale_factor, yfact=scale_factor, origin=(0, 0))


def iter_geojson_cell_centroids(cells_geojson, units_per_um=1.0,
                                cell_id_format="cellid_{:09d}-1",
                                cell_id_prop="cell_id"):
    """Yield ``(cell_id, x, y)`` polygon centroids from a cell-boundary GeoJSON.

    Each feature's geometry is rescaled by ``units_per_um`` (into microns) and
    reduced to its shapely centroid. The raw feature cell id (read from
    ``cell_id_prop``) is formatted with ``cell_id_format`` so the emitted id
    matches the barcode convention used downstream (e.g. Visium HD segmented
    barcodes ``cellid_{:09d}-1``); pass ``cell_id_format=None`` to keep the raw
    id verbatim.
    """
    with open(cells_geojson) as f:
        data = json.load(f)
    for feature in data.get("features", []):
        raw_id = feature["properties"][cell_id_prop]
        cid = cell_id_format.format(raw_id) if cell_id_format else str(raw_id)
        centroid = rescale_geometry(shape(feature["geometry"]), units_per_um).centroid
        yield cid, centroid.x, centroid.y


# Boundary file formats from which per-cell centroids can be derived. Extend this
# as new parsers are added (e.g. a Xenium-style cell_id/vertex_x/vertex_y CSV).
CENTROID_SUPPORTED_FORMATS = ("geojson",)


def write_cell_centroids_tsv(boundary_file, out_tsv, fmt="geojson", units_per_um=1.0,
                             cell_id_format="cellid_{:09d}-1", cell_id_prop="cell_id",
                             colname_cell_id="cell_id", colname_x="X", colname_y="Y",
                             opener=open):
    """Write a ``cell_id/X/Y`` TSV of centroids derived from ``boundary_file``.

    ``opener`` lets callers pass a gzip-aware opener (e.g. ``flexopen``). Returns
    the number of cells written. Raises ``ValueError`` for an unsupported format.
    """
    if fmt not in CENTROID_SUPPORTED_FORMATS:
        raise ValueError(
            f"Cannot derive centroids from boundary format {fmt!r}; "
            f"supported: {', '.join(CENTROID_SUPPORTED_FORMATS)}")
    n = 0
    with opener(out_tsv, "wt") as wf:
        wf.write(f"{colname_cell_id}\t{colname_x}\t{colname_y}\n")
        for cid, x, y in iter_geojson_cell_centroids(
                boundary_file, units_per_um, cell_id_format, cell_id_prop):
            wf.write(f"{cid}\t{x}\t{y}\n")
            n += 1
    return n

# Platform: BGI Stereo-seq (`--platform stereoseq`)

Built-in profile for BGI Stereo-seq / SAW output. Unlike every other platform, the inputs are **binary GEF files that only SAW can read**, so this profile shells out to the SAW binary (`--saw`) to expand them into text GEMs before the ordinary ingest runs.

Because a SAW sample is a set of files sharing a **name prefix** rather than a directory layout, this platform is addressed with `--in-prefix` instead of `--in-dir`.

## Expected inputs

All resolved by suffixing `--in-prefix`; each is optional except the tissue GEF.

| Role | File (auto-detected from `--in-prefix`) | Column | Purpose |
|------|------------------------------------------|--------|---------|
| Tissue GEF | `{prefix}.tissue.gef` | `gef` | **Required.** Pixel-level transcripts |
| Cell-bin GEF | `{prefix}.cellbin.gef` | `cellbin_gef` | Cell segmentation; enables cell-level analysis |
| H&E | `{prefix}_HE_regist.tif` | `--image` | Registered RGB histology |
| ssDNA | `{prefix}_ssDNA_regist.tif` | `--image` | Registered single-channel histology |
| DAPI | `{prefix}_DAPI_regist.tif` | `--image` | Registered single-channel histology |

Missing files are simply skipped. Naming is not always consistent across SAW versions — when a file does not match, give it explicitly (`--image type=hne,source=/p/whatever.tif`, or a `gef` / `cellbin_gef` sample-sheet column), which always wins over auto-detection.

**Defaults:** FICTURE `width=12`, `n_factor=12,24,48`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`; histology at 0.5 µm/pixel.

---
## 1. Single sample

```bash
cartloader run_together --platform stereoseq \
    --saw /path/to/saw \
    --in-prefix /data/C04687E314 \
    --out-dir OUT --width 12 --n-factor 24
```

The sample id defaults to the prefix basename (`C04687E314`); override with `--id`.

---
## 2. Multi-sample (sample sheet)

```
id            in_prefix                  cellbin_gef
C04687E314    /data/C04687E314
Y40320NA      /data/chipB/Y40320NA       /data/alt/Y40320NA.cellbin.gef
```

```bash
cartloader run_together --platform stereoseq --saw /path/to/saw \
    --samples samples.tsv --out-dir OUT -j 8
```

---
## Counts and gene names

Both ingests read **`ExonCount`** (the exonic subset of the MIDs), not `MIDCount`, and **rows with a zero count are dropped** — a large share of rows, since many MIDs have no exonic overlap.

!!! note "The count column is resolved by name, not position"
    `ExonCount` is the last column of the bin1 GEM but **not** of the cell-bin GEM, where `CellID` follows it. A positional filter (`grep -v -w 0$`) is therefore wrong on cellbin in both directions: it keeps zero-count rows whose `CellID` is non-zero, and discards real counts whose `CellID` is `0`. Both ingests read the header and select the named column, so layout differences between GEM variants are handled. A `CellID` of `0` is preserved and treated as unassigned later by `pixel2sptsv --ignore-ids`.

Features are read from **`geneName`** (the gene symbol), not `geneID` (Ensembl id). Symbols are what CartoScope displays, and they are what the exclude-feature regexes are written against — `Gm[0-9]`, `mt-`, `Rps`, `Rpl` match nothing at all against Ensembl ids, so using `geneID` would silently disable that filtering.

Override both together if a GEM variant differs:

```json
{ "ingest": { "csv_colnames": { "feature": "geneID", "count": "MIDCount" } } }
```

---
## Coordinates

Stereo-seq GEM coordinates are on a **0.5 µm grid**, so ingest passes `--units-per-um 2` to convert them to microns. The registered histology is at the same 0.5 µm/pixel, so it is georeferenced with `--px-per-um-x/y 2` (single-channel) or `--um-per-pixel 0.5` (RGB) and lands in the same frame with no further transform.

The GEM `#OffsetX` / `#OffsetY` headers are **not** applied. A non-zero offset would shift the transcripts relative to the histology, so the converter fails rather than silently producing a misaligned sample; override with `--allow-offset` only if the offset is already accounted for.

---
## Ingest: what SAW is run for

Two conversions, each deleted as soon as it has been read — the GEMs are tens of GB of text:

```bash
# pixel-level transcripts -> the FICTURE transcript
saw convert gef2gem --bin-size 1 --gef {prefix}.tissue.gef --gem bin1.gem

# cell segmentation -> the cell-level counts
saw convert gef2gem --cellbin-gef {prefix}.cellbin.gef \
                    --gef {prefix}.tissue.gef --cellbin-gem cellbin.gem
```

!!! warning "Re-running ingest re-runs SAW"
    The GEMs are removed after ingest, so `--restart` (or deleting `mk/sge.*.done`) pays the SAW conversion again. Use `--skip ingest` to resume from a completed ingest.

---
## Cell-level analysis

Stereo-seq is the one platform where cell identity **cannot be attached to the pixel-level transcript**: the cell-bin GEM's rows cannot be mapped back onto the bin1 GEM. So instead of a `cell_id` column on the transcript, the cell-bin GEM becomes a standalone pixel TSV:

```
cartloader convert_stereoseq_cellbin   # cellbin.gem -> X, Y, gene, count, cell_id
```

which the cells stage passes to `run_ficture2_multi_cells --tsv-list`. Cells are clustered from that file (`spatula pixel2sptsv`), and the resulting clusters are projected back onto **all** bin1 pixels by the pixel decode. Cell centroids come from the same step, so no separate cell-metadata file is needed. There are no cell polygons, so no cell boundaries are rendered.

Cell analysis runs only when `{prefix}.cellbin.gef` is present; without it the run is pixel-level (FICTURE → packaging).

!!! danger "Feature naming must match between the two GEMs"
    The cell clusters are projected onto a model trained from the bin1 features. If the two GEMs name genes differently (`geneID` vs `geneName`, Ensembl id vs symbol), the projection yields **near-empty cells rather than an error**. The converter therefore checks its features against the pixel run's feature file and fails the run if fewer than 50% overlap. Both ingests read `geneName` by default; to switch, set both together:

    ```json
    { "ingest": { "csv_colnames": { "feature": "geneID" } } }
    ```

---
## Images

The registered TIFs are plain (non-OME) TIFFs with no pixel-size metadata, so the scale is stated rather than detected. The profile handles the three standard suffixes automatically; declare anything else explicitly:

```bash
cartloader run_together --platform stereoseq --saw /path/to/saw \
  --in-prefix /data/C04687E314 --out-dir OUT \
  --image type=hne,source=/data/odd_name_HE.tif,um_per_pixel=0.5,georef_plain=true
```

Override the scale per image with `um_per_pixel` when a TIF is at a different resolution (e.g. `um_per_pixel=1.0` for a 2× downsampled export).

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)

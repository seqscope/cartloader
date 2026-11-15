# End-to-End Example Scripts

The repository's [`examples/`](https://github.com/seqscope/cartloader/tree/main/examples) directory holds complete shell scripts that take a public Xenium dataset from download to a CartoScope-ready catalog. Unlike the [quick starts](./intro.md), which run everything with one [`run_together`](../../reference/run_together.md) command, these scripts call each stage separately, so they double as a template for the [step-by-step](../step_by_step.md) route.

## Available scripts

| Dataset | Docker | Local |
|---------|--------|-------|
| Xenium human kidney, IgA nephropathy ([GSE264334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264334); Bull et al. 2024) | [`xenium_human_kidney_igan_bull2024_end_to_end.docker.sh`](https://github.com/seqscope/cartloader/blob/main/examples/xenium_human_kidney_igan_bull2024_end_to_end.docker.sh) | [`xenium_human_kidney_igan_bull2024_end_to_end.local.sh`](https://github.com/seqscope/cartloader/blob/main/examples/xenium_human_kidney_igan_bull2024_end_to_end.local.sh) |
| Xenium mouse brain, coronal subset (10x Genomics) | [`xenium_v1_mouse_brain_coronal_subset_end_to_end.docker.sh`](https://github.com/seqscope/cartloader/blob/main/examples/xenium_v1_mouse_brain_coronal_subset_end_to_end.docker.sh) | [`xenium_v1_mouse_brain_coronal_subset_end_to_end.local.sh`](https://github.com/seqscope/cartloader/blob/main/examples/xenium_v1_mouse_brain_coronal_subset_end_to_end.local.sh) |

The **Docker** scripts run every command in the `hyunminkang/cartloader` image (pulled at the start); the **Local** scripts call an installed `cartloader` (see [Installation](../../installation.md)).

## What the scripts do

1. Download and unpack the dataset.
2. Convert the Xenium transcripts into CartLoader's TSV layout ([`sge_convert`](../../reference/sge_convert.md)).
3. Train FICTURE models with 12 µm hexagons and 12, 24, and 48 factors, with single-molecule pixel decoding ([`run_ficture2_multi`](../../reference/run_ficture2_multi.md)).
4. Run two cell analyses that project the largest model ([`run_ficture2_multi_cells`](../../reference/run_ficture2_multi_cells.md)): `cartloader` (LDA-based clustering of the segmented cells) and `xeniumranger` (the Xenium Ranger clusters, centroids, and boundaries).
5. Package everything into PMTiles and a `catalog.yaml` ([`run_cartload2`](../../reference/run_cartload2.md), with `--pmtiles-pmpoint`).
6. Import the morphology image(s) ([`import_image`](../../reference/import_image.md)) and add them to the catalog.

The deployable output is `out/cartl/samples/rep1/`; a commented-out line at the end of each script shows how to upload it to S3.

## Running the kidney example

```bash
# Docker: writes data/ and out/ next to the script
bash examples/xenium_human_kidney_igan_bull2024_end_to_end.docker.sh

# Local: writes ./out (input unpacked into ./out/outs)
bash examples/xenium_human_kidney_igan_bull2024_end_to_end.local.sh
```

- The input (about 1.5 GB, including a 1.3 GB morphology image) is downloaded from the public `cartostore` S3 bucket once. If `transcripts.csv.gz` is already in the input directory, the download is skipped, so you can point the script at a copy you already have: `REAL_INDIR=/path/to/gse264334-reformatted` (Docker) or `INDIR=/path/to/gse264334-reformatted` (Local). `REAL_OUTDIR` / `OUTDIR` set the output directory the same way.
- `THREADS` and `JOBS` near the top of each script control parallelism (default: 4 threads, 2 jobs).
- The kidney panel has 377 genes, so the scripts use 50 gene bins (`BIN_COUNT=50`) for the point layers.
- On a laptop (4 threads, 2 jobs), the local kidney script finishes in about 6 minutes, excluding the download and the image import.

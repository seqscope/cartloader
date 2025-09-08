# Reference Overview

This section documents the `CartLoader` CLI by task. The Reference panel organizes task‑oriented docs for the CartLoaderCLI.

## SGE Preparation
<div class="grid cards forref" markdown>
-   <ul>
      <li>
        <h4>SGE Format Conversion</h4>
        <p><code>sge_convert</code>: Standardize raw platform outputs into the unified format; optional visualization; optional density-based filtering</p>
    </li>
    </ul>
    [Read](./sge_convert.md){ .md-button .md-button--primary .button-tight-small }
-   <ul>
    <li>
        <h4>SGE Additional Utilities</h4>
        <p><code>sge_stitch</code>: Stitch >=2 SGEs into one.</p>
        <p><code>sge_orientation</code>: Orientate SGE by rotation and flips. </p>
    </li>
    </ul>
    [Read](./sge_utilities.md){ .md-button .md-button--primary .button-tight-small }
</div>

## FICTURE Analysis
<div class="grid cards forref" markdown>
-   <ul>
      <li>
        <h4>Single-sample Analysis</h4>
        <p><code>run_ficture2</code>: Train spatial factors and decode pixels; writes a manifest for packaging</p>
    </li>
    </ul>
    [Read](./run_ficture2.md){ .md-button .md-button--primary .button-tight-small }
-   <ul>
    <li>
        <h4>Multi-sample Analysis</h4>
        <p><code>run_ficture2_multi</code>: Train spatial factors across >=2 samples and perform per‑sample pixel‑level decoding, writing individual manifests.</p>
    </li>
    </ul>
    [Read](./run_ficture2_multi.md){ .md-button .md-button--primary .button-tight-small }
</div>

## Asset Packaging
<div class="grid cards forref" markdown>
-   <ul>
      <li>
        <h4>Asset Packaging</h4>
        <p><code>run_cartload2</code>: Package SGE and (optionally) FICTURE outputs into PMTiles and write a catalog of layers.</p>
    </li>
    </ul>
    [Read](./run_cartload2.md){ .md-button .md-button--primary .button-tight-small }
-   <ul>
    <li>
        <h4>Multi-sample Packaging</h4>
        <p><code>run_cartload2_multi</code>: For each sample, package SGE and (optionally) FICTURE outputs into PMTiles and write a catalog of layers.</p>
    </li>
    </ul>
    [Read](./run_cartload2_multi.md){ .md-button .md-button--primary .button-tight-small }
-   <ul>
    <li>
        <h4>Import Images</h4>
        <p><code>import_image</code>: Convert an image into a georeferenced PMTiles layer (optionally reorient) and register it in the catalog.</p>
    </li>
    </ul>
    [Read](./import_image.md){ .md-button .md-button--primary .button-tight-small }
-   <ul>
    <li>
        <h4>Import Cell Analysis</h4>
        <p><code>import_xenium_cell</code>, <code>import_visiumhd_cell</code>: Import platform‑specific cell analysis outputs and add them to the catalog.</p>
    </li>
    </ul>
    [Read](./import_cell.md){ .md-button .md-button--primary .button-tight-small }
</div>


## Upload
<div class="grid cards forref" markdown>
-   <ul>
      <li>
        <h4>Data Repository Upload</h4>
        <p><code>upload_aws</code>, <code>upload_zenodo</code>: Publish outputs to S3 or Zenodo for sharing or web visualization.</p>
      </li>
    </ul>
    [Read](./upload_repository.md){ .md-button .md-button--primary .button-tight-small }
</div>


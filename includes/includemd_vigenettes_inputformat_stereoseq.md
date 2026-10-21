Stereo‑seq SGE includes one tab‑delimited “Bin1” gene expression matrix (GEM) file, `input.tsv.gz` in the starter data:

!!! info "TSV File Format"

    ```text
    geneID	        x	    y	    MIDCounts
    0610005C13Rik	6632	9074	1
    0610005C13Rik	8651	8935	1
    0610005C13Rik	7228	12814	2
    ```

    * `geneID`: Gene symbol.
    * `x`: X coordinate of each DNB on the capture chip.
    * `y`: Y coordinate of each DNB on the capture chip.
    * `MIDCounts`: UMI count for each gene at each DNB.

!!! note "Column names differ between SAW versions"
    This MOSTA file names its columns `geneID` (gene symbol) and `MIDCounts`. GEM files from current SAW versions use `geneName` and `ExonCount` instead, which are the defaults of the `bgi_stereoseq` preset, so this tutorial names the starter file's columns explicitly.

Stereo‑seq SGE includes one tab‑delimited “Bin1” gene expression matrix file:

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

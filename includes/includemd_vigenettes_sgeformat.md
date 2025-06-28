
??? "`transcripts.unsorted.tsv.gz`: transcript-indexed SGE in TSV"

    ```text
    X        Y        gene     count
    295.29   1422.35  Myo3a    0
    1745.54  1110.72  Med14    1
    1745.54  1110.72  Ntpcr    1
    ```

    * `X`: X coordinates in um
    * `Y`: Y coordinates in um
    * `gene`: gene symbols
    * `count`: expression count per pixel per gene

??? "`feature.clean.tsv.gz`: UMI counts on a per-gene basis in TSV"

    ```text
    gene           gene_id             count
    Gm29155        ENSMUSG00000100764  1
    Pcmtd1         ENSMUSG00000051285  431
    Gm26901        ENSMUSG00000097797  1
    ```
    * `gene`: gene symbols
    * `gene_id`: gene IDs
    * `count`: expression count per gene

??? "`coordinate_minmax.tsv`: X Y min/max includemd_vigenettes_sgeformat.mdcoordinates"

    ```text
    xmin	0.14
    xmax	2359.90
    ymin	0.23
    ymax	1439.95
    ```

    * `xmin` `xmax`: min and max X coordinates in um
    * `ymin` `ymax`: min and max Y coordinates in um

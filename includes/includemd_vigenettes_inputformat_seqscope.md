SeqScope provides SGE with three files:

??? "`barcodes.tsv.gz` – spatial barcode metadata"

    ```text
    AAAACAAAAACCTTCTTCGGACACTGGTCT	1	20	1	1	295288	1422349	0,1,0,0,0
    AAAACAAAAATCCTGTTATACATGCCATGG	2	45	1	1	1745544	1110720	2,2,1,0,1
    AAAACAAAACACGGGAAAAAACTATAGGTG	3	58	1	1	887244	250820	7,7,5,0,1
    ```

    * Column 1: Sorted spatial barcodes
    * Column 2: 1-based integer index of spatial barcodes, used in `matrix.mtx.gz`
    * Column 3: 1-based integer index from the full barcode that is in the STARsolo output
    * Column 4: Lane ID (fixed as `1`)  
    * Column 5: Tile ID (fixed as `1`)  
    * Column 6: X-coordinates
    * Column 7: Y-coordinates
    * Column 8: Five comma-separated numbers denote the count per spatial barcode for "Gene", "GeneFull", "Spliced", "Unspliced", and "Ambiguous".

??? "`features.tsv.gz` – feature metadata "

    ```text
    ENSMUSG00000100764	Gm29155	1	1,1,1,0,0
    ENSMUSG00000100635	Gm29157	2	0,0,0,0,0
    ENSMUSG00000100480	Gm29156	3	0,0,0,0,0
    ```

    * Column 1: Feature ID
    * Column 2: Feature symbol
    * Column 3: 1-based integer index of genes, used in `matrix.mtx.gz`
    * Column 4: Five comma-separated numbers denote the count per gene "Gene", "GeneFull", "Spliced", "Unspliced", and "Ambiguous".


??? "`matrix.mtx.gz` – expression count matrix "

    ```text
    %%MatrixMarket matrix coordinate integer general
    %
    33989 2928173 5404336
    2487 1 0 1 0 0 0
    5104 2 1 1 0 0 1
    ```
    
    * `Header`: Initial lines form the header, declaring the matrix's adherence to the [Market Matrix (MTX) format](https://math.nist.gov/MatrixMarket/formats.html), outlining its traits. This may include comments (lines beginning with `%`) for extra metadata, all marked by a “%”.
    * `Dimensions`: Following the header, the first line details the matrix dimensions: the count of rows (features), columns (barcodes), and non-zero entries.
    * `Data Entries`: Post-dimensions, subsequent lines enumerate non-zero entries in seven columns: row index (feature index), column index (barcode index), and five values (expression levels) corresponds to "Gene", "GeneFull", "Spliced", "Unspliced", and "Ambiguous".       
        * "Gene": represents unique, confidently mapped transcript count ("gene name"-based);
        * "GeneFull": denotes total transcript count assigned to gene (includes ambiguities).

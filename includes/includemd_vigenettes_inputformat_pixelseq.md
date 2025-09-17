
The Pixel-Seq SGE includes one tab-delimited text file, where each row represents a unique RNA molecule detected within a defined field of view (FOV), with associated genomic and spatial metadata.

!!! info "TSV file format"

    ```text
    FOVx  FOVy  xcoord     ycoord     UMIs     SpatialBarcode            MapStrand  Chrom  Start      STARmapping        Counts  geneID              geneName  bioType                 intronRatio
    017   005   26691.5    5786.5     TAACGAA  AAGGTTCATACCTACGACTGTTAA  16         1      24613729   150M               1       ENSMUSG00000101111  Gm28437   unprocessed_pseudogene  0.00
    016   007   27590.25   4639.0815  TAATATA  AATGGCGCATTTTGCTGTTTAGGC  16         2      39001628   138M2341N12M       1       ENSMUSG00000062997  Rpl35     protein_coding          0.00
    018   006   25099.945  5621.8335  AGTTGTA  CTGCATATGTGTCACCTAGGTAGC  16         1      24615767   150M               1       ENSMUSG00000101249  Gm29216   unprocessed_pseudogene  0.00
    ```

	* `FOVx`, `FOVy`: Field-of-view indices indicating the imaging tile coordinates in the x and y directions.
	* `xcoord`, `ycoord`: Spatial coordinates (in microns or pixels).
	* `UMIs`: Unique molecular identifier (UMI) sequence.
	* `SpatialBarcode`: Spatial barcode capturing the location and identity.
	* `MapStrand`: Indicates the strand orientation of the mapped read.
	* `Chrom`, `Start`: Chromosome number and start position of the mapped read on the genome.
	* `STARmapping`: Alignment pattern (CIGAR string) from the STAR aligner indicating how the transcript maps to the genome.
    * `Counts`: Number of times the UMI/gene combination was observed.
	* `geneID`, `geneName`: Ensembl gene ID and gene symbol.
	* `bioType`: Gene biotype.
    * `intronRatio`: Fraction of UMI counts assigned to intronic regions.

Visium HD slides use a 2×2 µm grid of barcoded squares (`square_002um`) for high-resolution spatial gene mapping. Use `filtered_feature_bc_matrix` to only process tissue-associated signals. 

!!! warning "ATTENTION"
    The file‑format examples below use a sample Visium HD dataset. Paths, IDs, and values are illustrative and may not match the dataset used in this tutorial.

SGE comprises several key files as below: 

??? "`filtered_feature_bc_matrix/barcodes.tsv.gz` – spatial barcode for tissue locations "

    ```text
    s_002um_00639_00600-1
    s_002um_00923_01639-1
    s_002um_01050_01530-1
    ```

    * Column 1: Spatial barcodes corresponding to specific locations on the tissue section.


??? "`filtered_feature_bc_matrix/features.tsv.gz` – feature metadata "

    ```text
    ENSMUSG00000051951	Xkr4	Gene Expression
    ENSMUSG00000025900	Rp1	    Gene Expression
    ENSMUSG00000025902	Sox17	Gene Expression
    ```

    * Column 1: Feature ID
    * Column 2: Feature symbol
    * Column 3: Feature type

??? "`filtered_feature_bc_matrix/matrix.mtx.gz` – expression count matrix "

    ```text
    %%MatrixMarket matrix coordinate integer general
    %
    19059 869411 11376563
    3606 1 1
    8957 1 1
    9733 1 1
    ```
    
    * `Header`: Initial lines form the header, declaring the matrix's adherence to the [Market Matrix (MTX) format](https://math.nist.gov/MatrixMarket/formats.html), outlining its traits. This may include comments (lines beginning with `%`) for extra metadata, all marked by a “%”.
    * `Dimensions`: Following the header, the first line details the matrix dimensions: the count of rows (features), columns (barcodes), and non-zero entries.
    * `Data Entries`: Post-dimensions, subsequent lines enumerate non-zero entries in seven columns: row index (feature index), column index (barcode index), and one value presenting the expression count per barcode per feature.

??? "`tissue_positions.parquet` – spatial barcode metadata "

    ```text
    barcode                 in_tissue   array_row   array_col   pxl_row_in_fullres  pxl_col_in_fullres
    s_002um_00434_01637-1   1           434         1637        3396.371014         9125.919898
    ```
    
    * `barcode`: Unique spatial barcode associated with each capture spot.
    * `in_tissue`: Binary flag (1 = in tissue, 0 = background) indicating whether the spot falls within the tissue boundary.
    * `array_row`, `array_col`: Integer indices representing the position of the spot on the capture array grid.
    * `pxl_row_in_fullres`, `pxl_col_in_fullres`: Floating point coordinates locating the spot in full-resolution tissue image pixels.

??? "`scalefactors_json.json` – pixel-to-micrometer scaling factors "

    ```json
    {
        "spot_diameter_fullres": 7.303953797779634,
        "bin_size_um": 2.0,
        "microns_per_pixel": 0.2738242950835738,
        "regist_target_img_scalef": 0.2505533,
        "tissue_lowres_scalef": 0.02505533,
        "fiducial_diameter_fullres": 1205.1523766336395,
        "tissue_hires_scalef": 0.2505533
    }
    ```

    * `spot_diameter_fullres`: Estimated diameter of a barcoded spot in full-resolution pixels.
    * `bin_size_um`: Physical size (in micrometers) of the smallest bin, typically 2.0 µm for Visium HD.
    * `microns_per_pixel`: Resolution of the full-res image, used to convert pixel distances to micrometers.
    * `regist_target_img_scalef`: Scaling factor applied during image registration to the target image.
    * `tissue_lowres_scalef`: Downscaling factor from full-res to low-resolution tissue image.
    * `fiducial_diameter_fullres`: Diameter of fiducial markers in full-resolution pixels, useful for alignment.
    * `tissue_hires_scalef`: Downscaling factor from full-res to high-resolution tissue image.

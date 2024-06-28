

## 10x_visium_hd 
### data and data format:

1. Data:
    - raw_feature_bc_matrix: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
    - only one gene expression
    - tissue_positions.parquet: spatial coordinates in parquet format
    - scalefactors_json.json: the metadata of their scale
2. Spatula does not support filtering dummy genes.
3. Commands:
    ```bash
    cartloader=/nfs/turbo/sph-hmkang/index/data/weiqiuc/cartloader/cartloader

    python $cartloader/sge_convert.py --platform 10x_visium_hd \
        --in-sge /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/10x_visium_hd/raw/raw_feature_bc_matrix \
        --in-parquet /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/10x_visium_hd/raw/tissue_positions.parquet \
        --in-json /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/10x_visium_hd/raw/scalefactors_json.json \
        --icols-mtx 1 \
        --units-per-um 3.6519768988898167 \
        --out-dir /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/10x_visium_hd/out \
        --print-feature-id \
        --spatula /nfs/turbo/sph-hmkang/tools/dev/spatula/bin/spatula 
    ```

## 10x_xenium
1. Data and data format:
- *.csv.gz
    - header:
    ```
    "transcript_id","cell_id","overlaps_nucleus","feature_name","x_location","y_location","z_location","qv"
    281474976710664,67490,1,"Bhlhe40",4843.046,6427.73,19.068869,40.0
    ```
    - columns:
    ```
    transcript_id:	Unique ID of the transcript
    cell_id:	Unique ID of the cell, consisting of a cell prefix and dataset suffix
    overlaps_nucleus:	Binary value to indicate if the transcript falls within the segmented nucleus of the cell (1) or not (0)
    feature_name:	Gene or control name
    x_location:	X location in µm
    y_location:	Y location in µm
    z_location:	Z location in µm
    qv:	Phred-scaled quality value (Q-Score) estimating the probability of incorrect call
    fov_name:	Name of the field of view (FOV) that the transcript falls within
    nucleus_distance:	The distance between the transcript and the cell centroid in µm based on segmentation mask boundaries. Transcripts localized within the nucleus have a distance of 0.0 µm.
    codeword_index:	An integer index for each codeword used to decode transcripts (same value as codewords in the gene_panel.json file).
    codeword_category:	Codeword category
    is_gene:	Boolean value to indicate whether transcript feature is "Gene Expression" or not.
    ```
2. Dummy genes: "BLANK_*", "NegControlCodeword_*", "NegControlProbe_*"
3. Commands:
    ```bash
    cartloader=/nfs/turbo/sph-hmkang/index/data/weiqiuc/cartloader/cartloader
    python $cartloader/sge_convert.py \
        --platform 10x_xenium \
        --in-csv /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/10x_xenium/raw/transcripts.csv.gz \
        --out-dir /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/10x_xenium/out \
        --csv-colname-x x_location \
        --csv-colname-y y_location \
        --csv-colname-feature-name feature_name \
        --csv-colname-phredscore qv \
        --csv-colnames-others cell_id overlaps_nucleus \
        --precision-um 2
    ```

## cosmX_smi
1. Data and data format:
- *_tx_file.csv.gz
    - format:
    ```
    "fov","cell_ID","x_global_px","y_global_px","x_local_px","y_local_px","z","target","CellComp"
    1,0,298943.990047619,19493.2809095238,896.371,4433.7571,0,"Snap25","None"
    ```
2. ctrl probes: NegPrb*
3. Commands:
    ```bash
    cartloader=/nfs/turbo/sph-hmkang/index/data/weiqiuc/cartloader/cartloader
    python $cartloader/sge_convert.py \
        --platform cosmx_smi \
        --in-csv /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/cosmx_smi/raw/Run5642_S3_Quarter_tx_file.csv.gz \
        --out-dir /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/cosmx_smi/out \
        --csv-colname-x x_global_px \
        --csv-colname-y y_global_px \
        --csv-colname-feature-name target \
        --csv-colnames-others cell_ID CellComp \
        --precision-um 2 
    ```

## stereoseq 

1. Data and data format:

- *_bin1.tsv.gz
    - format:
    ```
    geneID	x	y	MIDCounts
    0610005C13Rik	11864	11748	3
    ```
    - columns:
    ```
    geneID:	Gene ID
    x/y:	"X" and "Y" correspond to the coordinates of each DNB on the captured chip
    MIDCounts:	the number of UMI for each gene at each DNB.
    ```
2. No negative/ctrl probe was found.
3. Commands:
    ```bash
    cartloader=/nfs/turbo/sph-hmkang/index/data/weiqiuc/cartloader/cartloader
    python $cartloader/sge_convert.py \
        --platform bgi_stereoseq \
        --in-csv /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/bgi_stereoseq/raw/Mouse_brain_Adult_GEM_bin1.tsv.gz \
        --out-dir /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/bgi_stereoseq/out \
        --csv-colname-x x \
        --csv-colname-y y \
        --csv-colname-feature-name geneID \
        --csv-colnames-count MIDCounts \
        --precision-um 2 \
        --dry-run
    ```

## merscope
1. Data and data format
    - detected_transcripts.csv.gz
    - format:
        ```
        ,barcode_id,global_x,global_y,global_z,x,y,fov,gene,transcript_id
        747,0,2123.7725,156.28409,1.0,349.0,1891.4949,0,PDK4,ENST00000005178
        ```
    - columns:
        ```
        [BLANK] A numeric index that uniquely identifies a transcript within a field of view. The index is non-consecutive and ascending within each field of view. 
        barcode_id The row index of the identified transcript “barcode” in the codebook file (zero indexed). 
        fov and barcode_id are a composite primary key for the detected_ transcripts.csv table. 
        global_x The x-coordinate of the transcript (µm), relative to the space of the experimental region. global_x may be negative in some circumstances depending on the alignment between fields of view. 
        global_y The y-coordinate of the transcript (µm), relative to the space of the experimental region. global_y may be negative in some circumstances depending on the alignment between fields of view. 
        global_z The index of the z-position. The position is a zero-indexed integer. global_z can be translated into microns using the entry in the first row of the zPos column of the dataorganization.csv file sorted in ascending order. 
        x The x-coordinate of the transcript (µm), within the coordinate space of the field of view in which it was imaged. 
        y The y-coordinate of the transcript (µm), within the coordinate space of the field of view in which it was imaged. 
        fov The index of the field of view in which the transcript was imaged (zero indexed). fov and barcode_id are a composite primary key for the detected_transcripts. csv table. 
        gene The human readable name of the gene this transcript is associated with. Gene is derived from the “name” column of the codebook file. 
        transcript_id A unique identifier of the gene that this transcript is associated with. transcript_ id is derived from the “id” column of the codebook file. 
        cell_id IF cell segmentation was performed: The numeric index of the cell that contains this transcript, if any. If this transcript is not associated with a cell, cell_id will be -1. cell_id maps to the EntityID field found in the cell_boundaries.parquet and cell_metadata.csv. 
        ```
3. Commands

        'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name', 'csv_colname_feature_id', 'csv_colnames_others',
        'dummy_genes', 'precision_um',
        'colname_x', 'colname_y', 'colname_feature_name', 'colname_feature_id', 'colnames_count', 'colname_molecule_id'
    ```bash
    cartloader=/nfs/turbo/sph-hmkang/index/data/weiqiuc/cartloader/cartloader
    python $cartloader/sge_convert.py \
        --platform  vizgen_merscope \
        --in-csv /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/vizgen_merscope/raw/datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_detected_transcripts_S1R1.csv.gz\
        --out-dir /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/vizgen_merscope/out \
        --csv-colname-x global_x \
        --csv-colname-y global_y \
        --csv-colname-feature-name gene \
        --precision-um 2 
    ```

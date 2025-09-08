
The MERSCOPE input SGE includes one comma‑delimited text file in the following format:

!!! info "CSV file format"

    ```text
    ,barcode_id,global_x,global_y,global_z,x,y,fov,gene
    0,22,56.930107,3851.851,5.0,147.80061,1711.9067,0,Adgre1
    1,22,183.60107,3874.0085,5.0,1320.6799,1917.0692,0,Adgre1
    2,22,59.750736,3666.5576,5.0,132.66754,1844.2372,1,Adgre1
    ```

    * Column 1: Unique numeric index for each transcript within a field of view (non-consecutive, ascending).  
    * `barcode_id`: Zero-based index of the transcript barcode in the codebook; forms a composite key with `fov`.
    * `global_x`: Transcript x coordinates (µm) in the experimental region; may be negative due to alignment.  
    * `global_y`: Transcript y coordinates (µm) in the experimental region; may be negative due to alignment.  
    * `global_z`: Zero‑based z‑position index.
    * `x`: The x-coordinate of the transcript (µm), within the coordinate space of the field of view.
    * `y`: The y-coordinate of the transcript (µm), within the coordinate space of the field of view.
    * `fov`: Zero-based field of view index; forms a composite key with `barcode_id`.  
    * `gene`: Gene name associated with the transcript.

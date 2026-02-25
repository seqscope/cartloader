NanoString CosMx SMI produces single‑molecule spatial transcriptomics data as a comma‑separated values (CSV) table.

!!! info "CSV File Format"

    ```text
    "fov","cell_ID","x_global_px","y_global_px","x_local_px","y_local_px","z","target","CellComp"
    64,0,-473043,7954.533,4015.3,4246.2,1,"Gfap","None"
    64,0,-473022.9,7902.723,4035.48,4194.39,1,"Fth1","None"
    64,0,-473132,7836.476,3926.34,4128.143,1,"Ptn","None"
    ```

    * `fov`: The field of view (FOV) number.
    * `cell_ID`: Unique identifier for a single cell within a given FOV; 0 if background or unassigned molecules.
    * `x_global_px`, `y_global_px`: Global pixel coordinates relative to the tissue.
    * `x_local_px`, `y_local_px`: The x or y position (in pixels) relative to the given FOV.
    * `z`: Z-plane index representing the depth (optical section) where the transcript was detected.
    * `target`: Target name.
    * `CellComp`: Subcellular location of the target.

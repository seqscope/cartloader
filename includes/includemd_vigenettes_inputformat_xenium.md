10X Genomics Xenium platform outputs Spatial Gene Expression (SGE) data in a comma-separated values (CSV) format. 

!!! info "CSV File Format"

    ```text
    "transcript_id","cell_id","overlaps_nucleus","feature_name","x_location","y_location","z_location","qv"
    281827164036151,133793,1,"Sox10",2350.0232,4153.6846,16.592316,40.0
    281827164036152,133793,1,"Sox10",2350.2585,4154.5225,17.237207,10.514394
    281827164036164,151216,0,"Sox10",2350.5874,4277.699,14.285685,40.0
    ```

    * "`transcript_id`": Unique identifier for each detected transcript molecule.  
    * "`cell_id`": ID of the segmented cell associated with the transcript.
    * "`overlaps_nucleus`": 1 if the transcript overlaps the nucleus mask, 0 otherwise.  
    * "`feature_name`": Gene or other features name corresponding to the transcript.  
    * "`x_location`": X-coordinate of the transcript.  
    * "`y_location`": Y-coordinate of the transcript.  
    * "`z_location`": Z-coordinate (depth) of the transcript.  
    * "`qv`": Quality value indicating confidence in transcript detection.

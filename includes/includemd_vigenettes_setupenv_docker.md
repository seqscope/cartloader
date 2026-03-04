Define paths to all required binaries and resources in the docker.

```bash
# ====
# Replace each placeholder with the actual path on your system.  
# ====
work_dir=/path/to/work/directory                        # path to work directory that contains the downloaded input data
cd $work_dir

# Following are fixed path in the docker. DONT REVISE 
spatula=/app/cartloader/submodules/spatula/bin/spatula  # path to spatula executable
punkst=/app/cartloader/submodules/punkst                # path to FICTURE2 (punkst) executable
tippecanoe=/usr/local/bin/tippecanoe                    # path to tippecanoe executable
pmtiles=/usr/local/bin/pmtiles                          # path to pmtiles executable
aws=/usr/local/bin/aws                                  # path to AWS CLI binary

# (Optional) Define path to color map. 
cmap=/app/cartloader/assets/fixed_color_map_256.tsv     # Path to fixed color map. `CartLoader` includes one at cartloader/assets/fixed_color_map_256.tsv.

# Number of jobs
n_jobs=10                                               # If not specified, the number of jobs defaults to 1.

# Docker tag 
docker_tag=20260304a
```
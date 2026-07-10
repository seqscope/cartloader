!!! warning "Fixed paths in the Docker Image"
    Tools and dependencies have **fixed** paths in the Docker image (for example, `/usr/local/bin/pmtiles`).
    
    DO **NOT** modify paths of tools and dependencies manually.

```bash
# ====
# Replace user-specific placeholders with actual paths on your system.
# ====
work_dir=/path/to/work/directory                        # path to work directory that contains the downloaded input data
cd $work_dir

# The following paths are fixed inside Docker. Do not modify them.
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
docker_tag=20260306
```


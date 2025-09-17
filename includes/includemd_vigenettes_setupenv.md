
!!! warning "Pre-installed tools"

    Please ensure you have installed all required tools (See [Installation](../docs/installation.md)).

Define paths to all required binaries and resources. Optionally, specify a fixed color map for consistent rendering.

```bash
# ====
# Replace each placeholder with the actual path on your system.  
# ====

work_dir=/path/to/work/directory        # path to work directory that contains the downloaded input data
cd $work_dir

# Define paths to required binaries and resources
spatula=/path/to/spatula/binary         # path to spatula executable
punkst=/path/to/punkst/binary           # path to FICTURE2 (punkst) executable
tippecanoe=/path/to/tippecanoe/binary   # path to tippecanoe executable
pmtiles=/path/to/pmtiles/binary         # path to pmtiles executable
aws=/path/to/aws/cli/binary             # path to AWS CLI binary

# (Optional) Define path to color map. 
cmap=/path/to/color/map                 # Path to fixed color map. `CartLoader` includes one at cartloader/assets/fixed_color_map_256.tsv.

# Number of jobs
n_jobs=10                               # If not specified, the number of jobs defaults to 1.

# Activate the bioconda environment
conda activate ENV_NAME                 # replace ENV_NAME with your conda environment name
```

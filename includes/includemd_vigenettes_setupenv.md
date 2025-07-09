Define paths to all required binaries and resources, and target AWS S3 bucket. Optionally, specify a fixed color map for consistent rendering.

```bash
# ====
# Replace each placeholder with the actual path on your system.  
# ====
work_dir=/path/to/work/directory        # path to work directory that contains the downloaded input data
cd $work_dir

# Define paths to required binaries and resources
spatula=/path/to/spatula/binary         # path to spatula executable
punkst=/path/to/punkst/binary           # path to FICTURE2/punkst executable
tippecanoe=/path/to/tippecanoe/binary   # path to tippecanoe executable
pmtiles=/path/to/pmtiles/binary         # path to pmtiles executable
aws=/path/to/aws/cli/binary             # path to AWS CLI binary

# (Optional) Define path to color map. 
cmap=/path/to/color/map                 # Path to the fixed color map for rendering. cartloader provides a fixed color map at cartloader/assets/fixed_color_map_256.tsv.

# AWS S3 target location for cartostore
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace EXAMPLE_AWS_BUCKET with your actual S3 bucket name

# Number of jobs
n_jobs=10

# Activate the bioconda environment
conda activate ENV_NAME              # replace BIOENV_NAME with your bioconda environment name
```
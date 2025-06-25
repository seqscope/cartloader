Define paths to all required binaries and resources, and target AWS S3 bucket. Optionally, specify a fixed color map for consistent rendering.

```bash
# ====
# Replace each placeholder with the actual path on your system.  
# ====

# Define paths to required binaries and resources
spatula=/path/to/spatula/binary         # path to spatula executable
punkst=/path/to/punkst/binary           # path to FICTURE2/punkst executable
tippecanoe=/path/to/tippecanoe/binary   # path to tippecanoe executable
pmtiles=/path/to/pmtiles/binary         # path to pmtiles executable
aws=/path/to/aws/cli/binary             # path to AWS CLI binary

# (Optional) Define path to color map. 
cmap=/path/to/color/map                 # Path to the fixed color map for rendering. cartloader provides a fixed color map at cartloader/assets/fixed_color_map_256.tsv.

# AWS S3 target location for cartostore
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace this with your actual S3 bucket name
```
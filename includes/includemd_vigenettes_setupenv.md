```bash
# ====
# Replace the placeholders with values for your system.
# ====
work_dir=/path/to/work/directory        # directory holding the input data; outputs are written here too
cd $work_dir

n_jobs=10                               # number of parallel jobs (also used as threads per job below)

# Activate the Python environment in which CartLoader is installed
conda activate ENV_NAME                 # replace ENV_NAME (or: source /path/to/venv/bin/activate)
```

!!! info "No tool paths needed"
    Once `submodules/build.sh` has run (see [Installation](../docs/installation.md)), every `CartLoader` command finds `spatula`, `punkst`, `tippecanoe`, `pmtiles`, and `pmpoint` in the `submodules/` directory of the repository, and the default color map at `assets/default_color_map.tsv`. Pass `--spatula`, `--ficture2` (the `punkst` repository directory, which contains `bin/punkst`), `--tippecanoe`, `--pmtiles`, or `--cmap-file` only to use a different build or file. `upload_aws` runs the `aws` CLI from your `PATH` (`--aws` to override).

import rasterio
import os 

# return the number of bands in a raster image
def get_band_count(input_path, max_bands=4):
    """Return the number of bands in a raster image."""
    with rasterio.open(input_path) as src:
        Nbands= src.count
    if Nbands > max_bands:
        raise ValueError(f"Error: {input_path} has {Nbands} bands, which is more than the defined maximum of {max_bands}.")
    return Nbands

def check_north_up(filename):
    with rasterio.open(filename) as src:
        t = src.transform
        return not (t.e > 0 or t.b != 0 or t.d != 0)

def check_georef(file_path):
    try:
        with rasterio.open(file_path) as src:
            has_transform = src.transform is not None and src.transform != rasterio.Affine.identity()
            has_crs = src.crs is not None
            return has_transform and has_crs
    except Exception as e:
        raise ValueError(f"Error checking georeference for {file_path}: {e}")


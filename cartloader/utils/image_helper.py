import rasterio

def scheck_georeference(tif_path):
    """
    Check if a TIFF file is georeferenced (i.e., a GeoTIFF).
    """
    with rasterio.open(tif_path) as src:
        has_transform = src.transform != rasterio.Affine.identity()
        has_crs = src.crs is not None
        return has_transform and has_crs
    

orient2axisorder = {
    # rotate, flip_vertical, flip_horizontal: order_arg
    # no rotation,
    (None, False, False): "1,2",     # No rotation
    (None, False, True): "-1,2",     # Flip X-axis
    (None, True, False): "1,-2",     # Flip Y-axis
    (None, True, True): "-1,-2",     # Flip both axes

    # Rotate 90° clockwise
    ("90", False, False): "2,-1",    # Rotate 90° clockwise
    ("90", False, True): "-2,-1",    # Rotate 90°, then flip X
    ("90", True, False): "2,1",      # Rotate 90°, then flip Y
    ("90", True, True): "-2,1",      # Rotate 90°, then flip both

    # Rotate 180° clockwise
    ("180", False, False): "-1,-2",  # Rotate 180° clockwise
    ("180", False, True): "1,-2",    # Rotate 180°, then flip X
    ("180", True, False): "-1,2",    # Rotate 180°, then flip Y
    ("180", True, True): "1,2",      # Rotate 180°, then flip both (cancels out)

    # Rotate 270° clockwise
    ("270", False, False): "-2,1",   # Rotate 270° clockwise
    ("270", False, True): "2,1",     # Rotate 270°, then flip X
    ("270", True, False): "-2,-1",   # Rotate 270°, then flip Y
    ("270", True, True): "2,-1",     # Rotate 270°, then flip both
}


# Simplified map of equivalent transformations
def update_orient(rotation, flip_vertical, flip_horizontal, image_f):
    rotation = str(rotation) if rotation is not None else None
    orient_map = {
        # No transformation
        (None, False, False): (None, False, False),

        # Flip X-axis
        (None, False, True): (None, False, True),
        ("180", True, False): (None, False, True),  # Equivalent to flip X-axis

        # Flip Y-axis
        (None, True, False): (None, True, False),
        ("180", False, True): (None, True, False),  # Equivalent to flip Y-axis

        # Flip both axes (equivalent to 180° rotation)
        (None, True, True): (None, True, True),
        ("180", False, False): (None, True, True),  # Equivalent to flipping both axes

        # Rotate 90° clockwise
        ("90", False, False): ("90", False, False),
        ("270", True, True): ("90", False, False),  # Equivalent to Rotate 90°

        # Rotate 90° clockwise, then flip X
        ("90", False, True): ("90", False, True),
        ("270", True, False): ("90", False, True),  # Equivalent to Rotate 90° + flip X

        # Rotate 90° clockwise, then flip Y
        ("90", True, False): ("90", True, False),
        ("270", False, True): ("90", True, False),  # Equivalent to Rotate 90° + flip Y

        # Rotate 90° clockwise, then flip both
        ("90", True, True): ("90", True, True),
        ("270", False, False): ("90", True, True),  # Equivalent to Rotate 90° + flip both
    }

    if rotation is not None:
        if isinstance(rotation, int):
            rotation = str(rotation)
            if rotation not in ["90", "180", "270"]:
                raise ValueError(f"Error: Invalid rotate value ({rotation}) for {image_f}. Rotation must be None, 90, 180, or 270.")

    new_rot, new_vflip, new_hflip = orient_map.get(
        (rotation, flip_vertical, flip_horizontal)
    )

    return new_rot, new_vflip, new_hflip

# Update the orientation of a histology image
def update_orient_in_histology(histology):
    hist_path = histology["path"]

    rotation = histology.get("rotate", None)
    
    flip = histology.get("flip", None)
    flip_vertical = flip in [ "vertical", "both"]
    flip_horizontal = flip in ["horizontal", "both"]

    new_rot, new_vflip, new_hflip = update_orient(rotation, flip_vertical, flip_horizontal, hist_path)

    # Update the histology dictionary with the best solution
    histology["rotate"] = new_rot
    if new_vflip and new_hflip:
        histology["flip"] = "both"
    elif new_vflip:
        histology["flip"] = "vertical"
    elif new_hflip:
        histology["flip"] = "horizontal"
    else:
        histology["flip"] = None
    return histology
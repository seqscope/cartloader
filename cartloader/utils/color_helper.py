

def hex_to_rgb01(hex_code):
    """Convert RGB HEX code to normalized RGB values (0-1)"""
    hex_code = hex_code.lstrip('#')
    r, g, b = [int(hex_code[i:i+2], 16) / 255.0 for i in (0, 2, 4)]
    return r, g, b

def hex_to_rgb255(hex_code):
    """Parse an RGB HEX code (e.g., "#FFA500") to RGB values(0-255)"""
    hex_code = hex_code.lstrip('#')  # Remove the '#' character if present
    r = int(hex_code[0:2], 16)  # First two characters -> Red
    g = int(hex_code[2:4], 16)  # Next two characters -> Green
    b = int(hex_code[4:6], 16)  # Last two characters -> Blue
    return r, g, b


def normalize_rgb(r, g, b):
    """Return tuple in 0-1 range. Accepts 0-1 floats or 0-255 integers/floats."""
    # If any channel > 1.0, assume 0-255 scale
    if max(r, g, b) > 1.0:
        return (float(r) / 255.0, float(g) / 255.0, float(b) / 255.0)
    return (float(r), float(g), float(b))
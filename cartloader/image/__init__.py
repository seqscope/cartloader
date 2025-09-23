"""Shared helpers for image processing pipelines."""

from .stages import (
    configure_color_mode,
    register_georeference_stage,
    register_orientation_stage,
    register_png2pmtiles_pipeline,
)

__all__ = [
    "configure_color_mode",
    "register_georeference_stage",
    "register_orientation_stage",
    "register_png2pmtiles_pipeline",
]

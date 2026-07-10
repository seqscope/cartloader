# How do I define scaling (`--units-per-um`) for `sge_convert`?

`--units-per-um` means: **how many coordinate units in your input equal 1 micrometer (um)**.

In other words:

`units_per_um = input_coordinate_units / um`

## Quick rules

- If input coordinates are already in **um**: use `--units-per-um 1`.
- If input coordinates are in **pixels**: use  
  `--units-per-um = pixels_per_um = 1 / microns_per_pixel`.
- If input coordinates are in **nanometers (nm)**: use `--units-per-um 1000`.

## Platform guidance

- **10x Visium HD**: prefer `--scale-json` (from `scalefactors_json.json`) so scaling is computed automatically.
- **Other platforms / generic CSV**: set `--units-per-um` based on your file's coordinate unit metadata.

## Common mistakes

- Using `1` when coordinates are actually in pixels.
- Mixing units across samples in one project.
- Rounding scale too aggressively (keep enough precision).

## Sanity checks after conversion

- Tissue width/height should be in a realistic micrometer range for your assay.
- Spot/transcript spacing should look plausible in generated visualization.
- If geometry looks too small/large by a constant factor, re-check `--units-per-um`.

## Examples

```bash
# Coordinates already in um
cartloader sge_convert ... --units-per-um 1

# Coordinates in pixels, and microns_per_pixel = 0.2125
# units_per_um = 1 / 0.2125 = 4.705882...
cartloader sge_convert ... --units-per-um 4.705882

# Visium HD: use scale JSON instead of manual factor
cartloader sge_convert ... --scale-json /path/to/scalefactors_json.json
```

See also:

- [SGE Format Conversion reference](../reference/sge_convert.md)

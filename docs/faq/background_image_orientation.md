# How do I decide orientation for background image import?

Use the SGE visualization from `sge_convert` as your orientation reference. Compare your raw background image against it to determine whether the image needs rotation and/or flips before it is imported.

!!! warning "Use the north-up reference"
    `sge_convert --sge-visual --north-up` writes `xy_northup.tif` (next to `xy.png`), the SGE rendered in the north-up view that CartoScope displays. Compare against `xy_northup.tif`, not `xy.png`.

## Recommended workflow

1. Generate SGE reference image:
   - Run `sge_convert` with `--sge-visual --north-up`.
   - Locate `xy_northup.tif` in the SGE output directory.
2. Open your background image and `xy_northup.tif` side by side.
3. Match major landmarks (tissue outline, empty corners, distinctive structures).
4. Infer the transform:
   - Rotation: `--rotate 90`, `--rotate 180`, or `--rotate 270`
   - Flips: `--flip-horizontal` and/or `--flip-vertical`
5. Run a quick trial import and verify overlap in the viewer.
6. If needed, adjust transform and rerun.

## Applying the transform

- **[`import_image`](../reference/import_image.md)**: `--rotate 90|180|270`, `--flip-horizontal`, `--flip-vertical`.
- **`run_together`**: the same transform as the `rotate`, `flip_horizontal`, and `flip_vertical` keys of the image's entry (`--image`, the `--images` sheet, or the config). See [Image Modalities](../reference/run_together_images.md) for which image kinds honor them.

## Practical tips

- Start with rotation only; add flips only if left/right or up/down are mirrored.
- Rotation is applied before flips.
- Keep one transform per image source and reuse it for consistency.

See also:

- [Background Image Import reference](../reference/import_image.md)
- [SGE Format Conversion reference](../reference/sge_convert.md)

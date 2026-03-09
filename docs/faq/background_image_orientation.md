# How do I decide orientation for background image import?

Use `xy.png` from `sge_convert` as your orientation reference. Compare your raw background image against this `xy.png` to determine whether your image needs rotation and/or flips before `import_image`.

!!! warning "Reference `xy.png`"
      Ensure `xy.png` is already georeferenced in a north-up view (see `--north-up` in `sge_convert` )

## Recommended workflow

1. Generate SGE reference image:
   - Run `sge_convert` with `--sge-visual --north-up`.
   - Locate `xy.png` in the SGE output directory.
2. Open your background image and `xy.png` side by side.
3. Match major landmarks (tissue outline, empty corners, distinctive structures).
4. Infer transform for `import_image`:
   - Rotation: `--rotate 90`, `--rotate 180`, or `--rotate 270`
   - Flips: `--flip-horizontal` and/or `--flip-vertical`
5. Run a quick trial import and verify overlap in the viewer.
6. If needed, adjust transform and rerun.

## Practical tips

- Start with rotation only; add flips only if left/right or up/down are mirrored.
- Rotation is applied before flips in `import_image`.
- Keep one transform per image source and reuse it for consistency.

See also:

- [Background Image Import reference](../reference/import_image.md)
- [SGE Format Conversion reference](../reference/sge_convert.md)

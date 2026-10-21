# How can I speed up runs safely?

The main controls, available in `run_together` and most individual commands:

- `-j` / `--n-jobs`: number of parallel jobs (pipeline-level parallelism)
- `--threads`: threads per job (tool-level parallelism)

The two multiply: `-j 4 --threads 4` can keep about 16 cores busy.

## Tuning guide

- Start with `-j` × `--threads` close to the number of CPU cores.
- Increase gradually while monitoring RAM and disk I/O.
- If the machine becomes unstable (heavy swapping, stalled I/O), reduce `-j` first.
- Use higher `--threads` when the work is dominated by a few compute-heavy steps (`FICTURE` training and decoding), and higher `-j` when there are many small steps (packaging many gene bins, many samples).

## Rules of thumb

- `run_together` on a multi-sample run: `-j` parallelizes across samples and stages; per-sample packaging in a joint run also uses `-j`.
- Fewer gene bins (`--bin-count`) make packaging faster, at the cost of coarser grouping in the point layers.
- Shared servers and HPC login nodes: be conservative and follow local policy.

## Example

```bash
cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT -j 4 --threads 4
```

See also:

- [`run_together` reference](../reference/run_together.md)
- [Asset Packaging reference](../reference/run_cartload2.md)

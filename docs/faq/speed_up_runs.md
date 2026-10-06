# How can I speed up runs safely?

The main controls are:

- `--n-jobs`: number of parallel jobs (pipeline-level parallelism)
- `--threads`: threads per job (tool-level parallelism; available in some modules)

## Tuning guide

- Start with `--n-jobs` equal to available CPU cores divided by 2.
- Increase gradually while monitoring RAM and disk I/O.
- If the machine becomes unstable (heavy swapping, stalled I/O), reduce `--n-jobs` first.
- Use higher `--threads` only when the underlying step benefits from multithreading.

## Rules of thumb

- Many small files/steps: increase `--n-jobs`.
- Few compute-heavy steps: balance `--n-jobs` and `--threads`.
- Shared servers/HPC login nodes: be conservative and follow local policy.

## Example

```bash
cartloader run_cartload2 ... --n-jobs 8 --threads 4
```

See also:

- [FICTURE Analysis reference](../reference/run_ficture2.md)
- [Asset Packaging reference](../reference/run_cartload2.md)

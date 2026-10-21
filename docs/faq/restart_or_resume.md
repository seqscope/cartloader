# How do I resume a run, or rerun everything from scratch?

`run_together` and the individual `CartLoader` commands plan their work as a Makefile, so:

- **Resume:** rerun the same command. Steps whose outputs already exist are skipped, so an interrupted or failed run continues from where it stopped.
- **Rerun everything:** add `--restart`, which rebuilds every step (`make -B`).
- **Preview:** add `--dry-run` to write the Makefile and print the planned commands without executing them.

## With `run_together`

```bash
# Resume after an interruption: the same command again
cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT --id my_sample

# Run or skip particular stages (stages outside the selection are assumed done)
cartloader run_together ... --only ingest,ficture
cartloader run_together ... --skip images

# Rebuild everything
cartloader run_together ... --restart
```

Changing a setting such as `--n-factor` or `--width` adds new models alongside the existing ones on the next run. `--restart` is only needed to redo steps whose outputs already exist.

## With the individual commands

```bash
cartloader run_ficture2 ... --dry-run     # inspect commands only
cartloader run_ficture2 ... --restart     # force a full rerun
```

See also:

- [`run_together` reference](../reference/run_together.md)
- [FICTURE Analysis reference](../reference/run_ficture2.md)

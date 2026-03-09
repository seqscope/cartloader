# How do I resume a run, or rerun everything from scratch?

Most `CartLoader` modules use Makefile-based execution behavior:

- If outputs already exist, completed steps are typically skipped.
- Use `--restart` to ignore existing outputs and rerun all steps.
- Use `--dry-run` to generate the Makefile and inspect planned commands without executing.

## Practical workflow

1. Run normally first.
2. If a run was interrupted, rerun the same command to continue from missing targets.
3. If inputs/parameters changed and you want a clean rerun, add `--restart`.
4. If you want to inspect what will run, add `--dry-run`.

## Example

```bash
# Inspect commands only
cartloader run_ficture2 ... --dry-run

# Force full rerun
cartloader run_ficture2 ... --restart
```

See also:

- [FICTURE Analysis reference](../reference/run_ficture2.md)
- [Asset Packaging reference](../reference/run_cartload2.md)
- [SGE Conversion reference](../reference/sge_convert.md)

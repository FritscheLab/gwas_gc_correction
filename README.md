# GC Correction Tools

This repository provides an R-only utility to apply Genomic Control (GC) to GWAS or SAIGE step‑2 summary statistics.

- `scripts/gc_adjust_general.R`: R script using optparse + data.table to construct a chi‑square statistic from common columns (CHISQ, Tstat/var, BETA/SE, p.value.NA, p.value), estimate one lambda with optional MAC/MAF filter, optionally use a fixed prespecified lambda, and output GC‑adjusted results.
- `data/test_data.csv`: small synthetic dataset to try the tool.

## Quick start

### Run GC adjustment (R)

You need R with packages `data.table` and `optparse`.

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_r
```

This will write `out/test_gc_r.tsv.gz` and `out/test_gc_r.lambda.tsv`.

Use a fixed lambda (bypass estimation):

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_r_fixed \
  --lambda_fixed 1.12
```

### Multithreading

Set threads for data.table with `--cpus` (defaults to `parallel::detectCores()`):

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_r_threads \
  --cpus 8
```

## Usage details

See `docs/USAGE.md` for CLI options and examples. `docs/METHODS.md` summarizes the chi‑square fallbacks and lambda estimation.

## Notes

- Lambda is clamped to >= 1 unless you pass `--allow_deflation`.
- When using SAIGE binary outputs, the primary `p.value` is SPA-based; the scripts warn if they must derive chi‑square from that.
- Optional filtering for lambda estimation supports MAC or MAF. If both are provided, MAC is preferred.

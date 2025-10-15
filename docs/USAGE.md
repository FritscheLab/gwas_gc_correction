# Usage (R)

`scripts/gc_adjust_general.R` computes GC adjustment for summary statistics (tsv/tsv.gz).

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_r
```

With MAC filter for lambda estimation:

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_mac20 \
  --mac_thr 20
```

Use a fixed lambda (skip estimation and ignore filters):

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_fixed \
  --lambda_fixed 1.10
```

Outputs: `<out>.tsv.gz` and `<out>.lambda.tsv`.

## Multithreading

Control data.table threads with `--cpus` (defaults to `parallel::detectCores()`):

```bash
Rscript scripts/gc_adjust_general.R \
  --input data/test_data.csv \
  --out out/test_gc_threads \
  --cpus 8
```

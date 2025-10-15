# Methods overview

We compute a single genomic-control lambda and apply it to GWAS/SAIGE summary statistics. Lambda is defined as

$\lambda = \mathrm{median}(\chi^2) / F_{\chi^2_1}^{-1}(0.5)$,

where $F_{\chi^2_1}^{-1}(0.5)$ is the median of the $\chi^2$ distribution with 1 degree of freedom.

## Chi-square construction fallbacks

Priority order:

1. Use `CHISQ` if present.
2. Else use `Tstat` with a variance column (`var`, `varT`, `varTstar`, `var1`), via `Tstat^2 / var`.
3. Else use `BETA` and `SE` to get `Z^2`.
4. Else use the normal-approx p value (`p.value.NA`).
5. Else use the primary p value (`p.value`). For SAIGE binary, this is SPA-based, not a direct normal-approx chi-square.

## Optional MAC/MAF filter

Optionally apply a single MAC or MAF threshold when estimating lambda (not for output inclusion). If both thresholds are provided, MAC is preferred. Minor MAC can be derived from AC and either `N` or `AF` if requested.

## Adjustments

- `CHISQ_GC = CHISQ / lambda`
- `PVALUE_GC = 1 - F_{\chi^2_1}(CHISQ_GC)`
- If `SE` exists, `SE_GC = SE * sqrt(lambda)` so that `Z = BETA/SE` is rescaled.

We clamp lambda to 1 unless deflation is explicitly allowed.

## Fixed lambda option

You can bypass lambda estimation entirely by supplying a fixed value with `--lambda_fixed`. In this mode:

- No MAC/MAF filters are applied for estimation.
- The metadata file will indicate `fixed lambda (no estimation)` and `n_used_for_lambda = 0`.

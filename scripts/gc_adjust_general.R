#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---------------- CLI ----------------
opt_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "GWAS summary stats (tsv, tsv.gz)"),
  make_option(c("-o", "--out"),
    type = "character", default = "gc_adjusted",
    help = "Output prefix [default: %default]"
  ),

  # Generic columns
  make_option(c("--p_col"),
    type = "character", default = "p.value",
    help = "Primary p value column name [default: %default]"
  ),
  make_option(c("--pna_col"),
    type = "character", default = "p.value.NA",
    help = "Normal-approx p value column (if present) [default: %default]"
  ),
  make_option(c("--beta_col"),
    type = "character", default = "BETA",
    help = "BETA column [default: %default]"
  ),
  make_option(c("--se_col"),
    type = "character", default = "SE",
    help = "SE column [default: %default]"
  ),
  make_option(c("--chisq_col"),
    type = "character", default = "CHISQ",
    help = "CHISQ column name if present [default: %default]"
  ),

  # SAIGE-friendly names (script will also try smart fallbacks)
  make_option(c("--tstat_col"),
    type = "character", default = "Tstat",
    help = "Score test statistic column [default: %default]"
  ),
  make_option(c("--var_col"),
    type = "character", default = "var",
    help = "Variance of score statistic column; fallbacks varT, varTstar, var1 are auto-tried [default: %default]"
  ),

  # Optional single filter used ONLY for lambda estimation
  make_option(c("--mac_thr"),
    type = "double", default = NA,
    help = "MAC threshold used only to estimate lambda, example 20 [default: NA]"
  ),
  make_option(c("--mac_col"),
    type = "character", default = "AC_Allele2",
    help = "Column to treat as allele count for MAC filtering [default: %default]"
  ),
  make_option(c("--maf_thr"),
    type = "double", default = NA,
    help = "MAF threshold used only to estimate lambda [default: NA]"
  ),
  make_option(c("--maf_col"),
    type = "character", default = "AF_Allele2",
    help = "MAF column name [default: %default]"
  ),
  make_option(c("--n_col"),
    type = "character", default = "N",
    help = "Sample size column, used to derive minor MAC when available [default: %default]"
  ),
  make_option(c("--derive_minor_mac"),
    action = "store_true", default = TRUE,
    help = "If TRUE and both allele count and AF or N exist, derive minor MAC = min(AC, 2N-AC) for filtering [default: %default]"
  ),

  # Behavior
  make_option(c("--allow_deflation"),
    action = "store_true", default = FALSE,
    help = "Allow lambda < 1 to deflate statistics [default: off]"
  ),
  make_option(c("--min_n"),
    type = "integer", default = 1000,
    help = "Warn if fewer than this number of variants are used for lambda [default: %default]"
  ),
  make_option(c("--lambda_fixed"),
    type = "double", default = NA,
    help = "Use this fixed lambda value and skip estimation (overrides filters)"
  ),
  make_option(c("--cpus"),
    type = "integer", default = NA,
    help = "Threads for data.table; default uses parallel::detectCores()"
  )
)

opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input)) stop("Provide --input")

# ---------------- Threads ----------------
cpus_detected <- tryCatch(
  {
    parallel::detectCores(logical = TRUE)
  },
  error = function(e) NA_integer_
)
cpus_use <- opt$cpus
if (is.na(cpus_use) || cpus_use <= 0) {
  cpus_use <- if (is.na(cpus_detected)) 1L else cpus_detected
}
if (!is.na(cpus_detected)) {
  cpus_use <- max(1L, min(as.integer(cpus_use), as.integer(cpus_detected)))
} else {
  cpus_use <- max(1L, as.integer(cpus_use))
}
data.table::setDTthreads(cpus_use)

# ---------------- IO ----------------
dt <- fread(opt$input, showProgress = FALSE)
n_total <- nrow(dt)
if (n_total == 0) stop("Input has zero rows")

has_col <- function(x) !is.na(x) && x %in% names(dt)

pick_first_col <- function(cands) {
  cands <- cands[!is.na(cands)]
  cands[cands %in% names(dt)][1]
}

# ---------------- Build chi-square with fallbacks ----------------
# Order:
#   1) CHISQ column
#   2) Tstat + var (try var, then varT, varTstar, var1)
#   3) BETA + SE
#   4) p.value.NA  -> chi
#   5) p.value     -> chi  (note: SPA for SAIGE binary)
tiny <- 1e-300
qmed <- qchisq(0.5, df = 1) # 0.4549364

chi <- rep(NA_real_, n_total)
chi_source <- NA_character_
var_col_used <- NA_character_

# 1) CHISQ
if (has_col(opt$chisq_col)) {
  chi <- as.numeric(dt[[opt$chisq_col]])
  chi_source <- paste0("CHISQ from ", opt$chisq_col)
}

# 2) Tstat + var (if needed)
if (all(!is.finite(chi))) {
  tcol <- pick_first_col(c(opt$tstat_col, "TSTAT", "T_STAT", "STAT"))
  vcol <- pick_first_col(c(opt$var_col, "varT", "varTstar", "var1", "VAR", "Var"))
  if (!is.na(tcol) && !is.na(vcol)) {
    tstat <- as.numeric(dt[[tcol]])
    v <- as.numeric(dt[[vcol]])
    good <- is.finite(tstat) & is.finite(v) & v > 0
    chi[good] <- (tstat[good]^2) / v[good]
    chi_source <- paste0("Tstat^2/", vcol)
    var_col_used <- vcol
  }
}

# 3) BETA + SE
if (all(!is.finite(chi)) && has_col(opt$beta_col) && has_col(opt$se_col)) {
  b <- as.numeric(dt[[opt$beta_col]])
  se <- as.numeric(dt[[opt$se_col]])
  good <- is.finite(b) & is.finite(se) & se > 0
  chi[good] <- (b[good] / se[good])^2
  chi_source <- "Z^2 from BETA/SE"
}

# 4) p.value.NA
if (all(!is.finite(chi)) && has_col(opt$pna_col)) {
  pna <- as.numeric(dt[[opt$pna_col]])
  pna[!is.finite(pna) | pna <= 0] <- tiny
  pna[pna > 1] <- 1
  chi <- qchisq(pna, df = 1, lower.tail = FALSE)
  chi_source <- paste0("qchisq(", opt$pna_col, ")")
}

# 5) p.value
if (all(!is.finite(chi)) && has_col(opt$p_col)) {
  p <- as.numeric(dt[[opt$p_col]])
  p[!is.finite(p) | p <= 0] <- tiny
  p[p > 1] <- 1
  chi <- qchisq(p, df = 1, lower.tail = FALSE)
  chi_source <- paste0("qchisq(", opt$p_col, ")")
  warning("Using p.value to build chi. For SAIGE binary traits this is SPA-based and not the normal-approx chi-square.")
}

if (all(!is.finite(chi))) stop("Could not construct chi-square from any input columns")

# ---------------- Optional single MAC or MAF filter for lambda ----------------
use_fixed_lambda <- !is.na(opt$lambda_fixed)

use_mac <- !is.na(opt$mac_thr)
use_maf <- !is.na(opt$maf_thr)
if (use_mac && use_maf) {
  warning("Both --mac_thr and --maf_thr were supplied. Using MAC and ignoring MAF.")
  use_maf <- FALSE
}

keep_lambda <- is.finite(chi)
filt_label <- if (use_fixed_lambda) "fixed lambda (no estimation)" else "no MAC/MAF filter"

if (!use_fixed_lambda && use_mac) {
  if (!has_col(opt$mac_col)) stop("Requested MAC filter but --mac_col not found in data")
  ac <- as.numeric(dt[[opt$mac_col]])
  mac_used <- ac

  if (opt$derive_minor_mac && has_col(opt$n_col)) {
    N <- as.numeric(dt[[opt$n_col]])
    twoN <- 2 * N
    mac_minor <- ifelse(is.finite(twoN), pmin(ac, twoN - ac), NA_real_)
    mac_used <- ifelse(is.finite(mac_minor) & mac_minor >= 0, mac_minor, mac_used)
  } else if (opt$derive_minor_mac && has_col(opt$maf_col)) {
    af <- as.numeric(dt[[opt$maf_col]])
    twoN_est <- ifelse(is.finite(af) & af > 0, round(ac / af), NA_real_)
    mac_minor <- pmin(ac, pmax(twoN_est - ac, NA_real_), na.rm = FALSE)
    mac_used <- ifelse(is.finite(mac_minor) & mac_minor >= 0, mac_minor, mac_used)
  }

  keep_lambda <- keep_lambda & is.finite(mac_used) & mac_used >= opt$mac_thr
  filt_label <- sprintf("MAC >= %s on %s", format(opt$mac_thr, scientific = FALSE), opt$mac_col)
}

if (!use_fixed_lambda && use_maf) {
  if (!has_col(opt$maf_col)) stop("Requested MAF filter but --maf_col not found in data")
  maf <- as.numeric(dt[[opt$maf_col]])
  keep_lambda <- keep_lambda & is.finite(maf) & maf >= opt$maf_thr
  filt_label <- sprintf("MAF >= %s on %s", format(opt$maf_thr, scientific = FALSE), opt$maf_col)
}

n_used <- if (use_fixed_lambda) 0L else sum(keep_lambda, na.rm = TRUE)
if (!use_fixed_lambda && n_used < opt$min_n) {
  warning(sprintf("Only %d variants available for lambda estimation after filter. You set --min_n=%d", n_used, opt$min_n))
}

# ---------------- Lambda and GC adjustments ----------------
lambda_raw <- if (use_fixed_lambda) opt$lambda_fixed else median(chi[keep_lambda], na.rm = TRUE) / qmed
lambda_use <- if (use_fixed_lambda) opt$lambda_fixed else if (opt$allow_deflation) lambda_raw else max(1.0, lambda_raw)

chi_gc <- chi / lambda_use
p_gc <- pchisq(chi_gc, df = 1, lower.tail = FALSE)

dt[, CHISQ_FROM := chi_source]
dt[, CHISQ := chi]
dt[, CHISQ_GC := chi_gc]
dt[, PVALUE_GC := p_gc]
dt[, LAMBDA_GC := lambda_use]

if (has_col(opt$se_col)) {
  dt[, SE_GC := as.numeric(dt[[opt$se_col]]) * sqrt(lambda_use)]
}

# ---------------- Write outputs ----------------
out_main <- paste0(opt$out, ".tsv.gz")
out_meta <- paste0(opt$out, ".lambda.tsv")
fwrite(dt, out_main, sep = "\t")

meta <- data.table(
  n_total = n_total,
  n_used_for_lambda = n_used,
  filter_for_lambda = filt_label,
  chisq_source = chi_source,
  var_col_used = ifelse(is.na(var_col_used), "", var_col_used),
  lambda_raw = lambda_raw,
  lambda = lambda_use
)
fwrite(meta, out_meta, sep = "\t")

message("Done. Wrote: ", out_main, " and ", out_meta)

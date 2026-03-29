#!/usr/bin/env Rscript

# Scenario-driven benchmark runner for stglmnet vs fmrilss.
# Usage:
#   Rscript bench/run_stglmnet_vs_fmrilss_suite.R

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}
source("bench/stglmnet_vs_fmrilss_harness.R")
has_hrfals <- requireNamespace("hrfals", quietly = TRUE)
if (!has_hrfals) {
  cat("Package 'hrfals' not available: running without hrfals methods.\n")
}

scenarios <- default_stglmnet_vs_fmrilss_scenarios("wide")

suite <- benchmark_stglmnet_vs_fmrilss_suite(
  scenarios = scenarios,
  base_args = list(
    n_reps = 20,
    n_tp = c(120, 120),
    n_runs = 2,
    n_trials = 48,
    n_vox = 16,
    isi = 1.2,
    cv_folds = 5,
    include_oasis = TRUE,
    include_hrfals = has_hrfals,
    hrfals_modes = c("r", "cpp"),
    hrfals_fit_methods = c("ls_svd_only", "ls_svd_1als", "cf_als"),
    hrfals_basis = "spmg2",
    trace = FALSE
  ),
  base_seed = 20260307L,
  trace = TRUE
)

cat("\n=== Pairwise Summary (primary st vs fmrilss_lss) ===\n")
print(suite$pairwise_summary[order(suite$pairwise_summary$st_better_corr_rate, decreasing = TRUE), ])

cat("\n=== Top Method by Correlation Per Scenario ===\n")
print(suite$top_by_cor[order(suite$top_by_cor$mean_cor, decreasing = TRUE), ])

cat("\n=== Scenario Method Summary (first 20 rows) ===\n")
print(utils::head(suite$scenario_method_summary, 20))

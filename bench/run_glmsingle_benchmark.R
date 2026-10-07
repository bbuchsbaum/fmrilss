# End-to-end timing of glmsingle() against pinned Python GLMsingle.
#
# 1. OPENBLAS_NUM_THREADS=1 python -I bench/glmsingle/bench_python.py OUT 8 200 20000
# 2. OPENBLAS_NUM_THREADS=1 Rscript bench/run_glmsingle_benchmark.R OUT
#
# Uses the same simulated inputs and fixed ON-OFF thresholds on both sides.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
suppressMessages(if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else library(fmrilss))
meta <- jsonlite::fromJSON(file.path(out, "python_result.json"))
rd <- function(f, n) readBin(f, "double", n = n, size = 4L, endian = "little")
Y <- lapply(seq_len(meta$n_runs), function(r) {
  matrix(rd(file.path(out, "inputs", sprintf("data_run%d.f32", r)), meta$n_vox * meta$n_time),
         meta$n_time, meta$n_vox)  # file is voxels x time in C order
})
design <- lapply(seq_len(meta$n_runs), function(r) {
  matrix(rd(file.path(out, "inputs", sprintf("design_run%d.f32", r)), meta$n_time * meta$n_cond),
         meta$n_time, meta$n_cond, byrow = TRUE)
})
gc()
t0 <- proc.time()[["elapsed"]]
fit <- glmsingle(Y, design, tr = meta$tr, stimdur = meta$stimdur,
                 brain_r2 = meta$brainR2, pc_r2_cutoff = meta$pcR2cutoff,
                 extras_in_denoise = "with_pcs", zero_sd_cv = "python", verbose = FALSE)
elapsed <- proc.time()[["elapsed"]] - t0
pyb <- matrix(rd(file.path(out, "python_betas_typed.f32"), meta$n_vox * meta$n_trials),
              meta$n_trials, meta$n_vox)
cor_med <- stats::median(vapply(seq_len(meta$n_vox), function(v) {
  suppressWarnings(stats::cor(fit$typed$betasmd[, v], pyb[, v]))
}, numeric(1)), na.rm = TRUE)
res <- list(
  dims = sprintf("%d runs x %d TRs x %d voxels, %d trials", meta$n_runs, meta$n_time, meta$n_vox, meta$n_trials),
  blas = sessionInfo()$BLAS,
  python_seconds = meta$python_seconds, r_seconds = elapsed,
  speedup = meta$python_seconds / elapsed,
  pcnum = c(python = meta$pcnum, r = fit$typed$pcnum),
  median_beta_correlation = cor_med,
  max_rel_diff = max(abs(fit$typed$betasmd - pyb), na.rm = TRUE) / max(abs(pyb), na.rm = TRUE),
  stage_seconds = fit$timing
)
print(res)
saveRDS(res, file.path(out, "r_result.rds"))

# Compare fmrilss with Python LSS implementations on a shared simulation.
#
# Usage: Rscript run_comparison.R <sim_dir> [reps]
#
# Expects simulate.R and python_lss.py to have been run on <sim_dir>. Writes
# fmrilss timings and the accuracy/timing table to <sim_dir>/comparison.csv.

suppressPackageStartupMessages(library(fmrilss))
args <- commandArgs(trailingOnly = TRUE)
sim_dir <- if (length(args) >= 1) args[[1]] else "sim"
reps <- if (length(args) >= 2) as.integer(args[[2]]) else 3L

man <- jsonlite::fromJSON(file.path(sim_dir, "manifest.json"))
rd <- function(nm, nc) {
  matrix(readBin(file.path(sim_dir, paste0(nm, ".bin")), "double",
                 n = man$n_scans * nc, endian = "little"), man$n_scans, nc)
}
Y <- rd("Y", man$n_vox)
X <- rd("X", man$n_trials)
Z <- rd("Z", man$n_z)
motion <- rd("motion", man$n_motion)
cond <- readLines(file.path(sim_dir, "condition.txt"))
beta_true <- matrix(readBin(file.path(sim_dir, "beta_true.bin"), "double",
                            n = man$n_trials * man$n_vox, endian = "little"),
                    man$n_trials, man$n_vox)

fm <- list(
  fmrilss_lss = function() lss(Y, X, Z, motion),
  fmrilss_lss_n = function() lss(Y, X, Z, motion, trial_groups = cond),
  fmrilss_lss_n_cpp = function() lss(Y, X, Z, motion, trial_groups = cond,
                                     method = "cpp_optimized"),
  fmrilss_lss_n_ar1 = function() lss(Y, X, Z, motion, trial_groups = cond,
                                     prewhiten = list(method = "ar", p = 1)),
  fmrilss_lss_n_ar1_voxel = function() {
    lss(Y, X, Z, motion, trial_groups = cond,
        prewhiten = list(method = "ar", p = 1, pooling = "voxel"))
  },
  fmrilss_lss_n_ridge = function() {
    lss(Y, X, Z, motion, trial_groups = cond, ridge = c(1, 0))
  },
  fmrilss_lss_n_ridge_ar1_voxel = function() {
    lss(Y, X, Z, motion, trial_groups = cond, ridge = c(1, 0),
        prewhiten = list(method = "ar", p = 1, pooling = "voxel"))
  }
)

time_it <- function(f) {
  times <- numeric(reps)
  for (r in seq_len(reps)) times[r] <- system.time(out <- f())[["elapsed"]]
  list(beta = unname(as.matrix(out)), time = stats::median(times))
}

results <- lapply(fm, time_it)
py_times <- jsonlite::fromJSON(file.path(sim_dir, "py_timings.json"))
for (nm in names(py_times)) {
  path <- file.path(sim_dir, paste0("py_", nm, ".bin"))
  if (!file.exists(path)) next
  B <- matrix(readBin(path, "double", n = man$n_trials * man$n_vox, endian = "little"),
              man$n_trials, man$n_vox)
  results[[paste0("python_", nm)]] <- list(beta = B, time = py_times[[nm]])
}

# Accuracy: trial-wise correlation with the truth per voxel (overall and
# within condition, which isolates trial-to-trial variability), and RMSE.
within_r <- function(B) {
  mean(vapply(unique(cond), function(cc) {
    idx <- cond == cc
    mean(vapply(seq_len(ncol(B)), function(v) {
      stats::cor(B[idx, v], beta_true[idx, v])
    }, numeric(1)))
  }, numeric(1)))
}
tab <- do.call(rbind, lapply(names(results), function(nm) {
  B <- results[[nm]]$beta
  data.frame(
    method = nm,
    seconds = results[[nm]]$time,
    r_trial = mean(vapply(seq_len(ncol(B)), function(v) stats::cor(B[, v], beta_true[, v]), 1)),
    r_within_condition = within_r(B),
    rmse = sqrt(mean((B - beta_true)^2)),
    stringsAsFactors = FALSE
  )
}))

# Numerical agreement between algebraically identical estimators
agree <- function(a, b) {
  if (is.null(results[[a]]) || is.null(results[[b]])) return(NA_real_)
  max(abs(results[[a]]$beta - results[[b]]$beta))
}
tab$max_abs_diff_vs_equivalent <- NA_real_
pairs <- c(python_nilearn_ols = "fmrilss_lss_n", python_numpy_loop_n = "fmrilss_lss_n",
           python_numpy_loop = "fmrilss_lss", python_numpy_closed = "fmrilss_lss")
for (a in names(pairs)) tab$max_abs_diff_vs_equivalent[tab$method == a] <- agree(a, pairs[[a]])

tab$n_vox <- man$n_vox
tab$n_trials <- man$n_trials
tab$n_scans <- man$n_scans
utils::write.csv(tab, file.path(sim_dir, "comparison.csv"), row.names = FALSE)
print(tab, digits = 4, row.names = FALSE)

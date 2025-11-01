#!/usr/bin/env Rscript
# Debug script to trace sign issues in SBHM pipeline

devtools::load_all()
library(fmrihrf)

set.seed(3301)

# Small synthetic setup (same as failing test)
Tlen <- 120; V <- 3; r <- 4
sframe <- sampling_frame(blocklens = Tlen, TR = 1)

# Simple library matrix; HRF-like decays
H <- cbind(
  exp(-seq(0, 30, length.out = Tlen)/4),
  exp(-seq(0, 30, length.out = Tlen)/6),
  exp(-seq(0, 30, length.out = Tlen)/8),
  exp(-seq(0, 30, length.out = Tlen)/10)
)
sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

# Trials and design
onsets <- seq(10, 90, by = 10)
ntrials <- length(onsets)
design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

# Build per-trial SBHM-basis regressors for simulation
regs <- vector("list", ntrials)
for (t in seq_len(ntrials)) {
  rr_t <- fmrihrf::regressor(onsets = onsets[t], hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
  regs[[t]] <- fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
}

# Assign each voxel a library HRF coordinate and trial amplitudes
idx <- sample(seq_len(ncol(sbhm$A)), V, replace = TRUE)
alpha <- sbhm$A[, idx, drop = FALSE]
amps_true <- matrix(rnorm(ntrials * V, mean = 2, sd = 0.4), ntrials, V)

cat("=== TRUE SIMULATION PARAMETERS ===\n")
cat("Library indices used:", idx, "\n")
cat("True alpha signs (first element):", sign(alpha[1,]), "\n")
cat("True amplitude range:", range(amps_true), "\n")
cat("Mean true amplitude:", mean(amps_true), "\n\n")

# Simulate Y = sum_t amps[t]*X_t * alpha + noise
Y <- matrix(rnorm(Tlen * V, sd = 0.3), Tlen, V)
for (v in seq_len(V)) {
  for (t in seq_len(ntrials)) {
    Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(regs[[t]] %*% alpha[, v])
  }
}

# Run SBHM with diagnostics
cat("=== RUNNING SBHM PIPELINE ===\n")

# 1) Prepass
pre <- fmrilss:::sbhm_prepass(
  Y = Y,
  sbhm = sbhm,
  design_spec = design_spec,
  ridge = list(mode = "fractional", lambda = 0.01)
)
cat("Prepass beta_bar signs (first element):", sign(pre$beta_bar[1,]), "\n")
cat("Prepass beta_bar vs true alpha correlation:",
    sapply(1:V, function(v) cor(pre$beta_bar[,v], alpha[,v])), "\n\n")

# 2) Matching (before any modification in lss_sbhm)
m <- fmrilss:::sbhm_match(
  beta_bar = pre$beta_bar,
  S = sbhm$S,
  A = sbhm$A,
  shrink = list(tau = 0, ref = sbhm$ref$alpha_ref, snr = NULL),
  topK = 1,
  whiten = TRUE,
  orient_ref = TRUE
)
cat("Matched library indices:", m$idx, "(true:", idx, ")\n")
cat("Matched alpha_hat signs (first element):", sign(m$alpha_hat[1,]), "\n")
cat("Matched alpha_hat vs true alpha correlation:",
    sapply(1:V, function(v) cor(m$alpha_hat[,v], alpha[,v])), "\n\n")

# 3) OASIS trial-wise fit
res <- lss_sbhm(Y, sbhm, design_spec,
                prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
                return = "both")

cat("=== FINAL RESULTS ===\n")
cat("Final alpha_coords signs (first element):", sign(res$alpha_coords[1,]), "\n")
cat("Final alpha_coords vs true alpha correlation:",
    sapply(1:V, function(v) cor(res$alpha_coords[,v], alpha[,v])), "\n")
cat("Recovered amplitude range:", range(res$amplitude), "\n")
cat("Mean recovered amplitude:", mean(res$amplitude), "\n")
cat("Overall amplitude correlation:", cor(as.vector(res$amplitude), as.vector(amps_true)), "\n")
cat("Per-voxel amplitude correlations:",
    sapply(1:V, function(v) cor(res$amplitude[,v], amps_true[,v])), "\n")

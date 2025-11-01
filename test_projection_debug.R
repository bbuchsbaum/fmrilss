#!/usr/bin/env Rscript
# Debug projection step

devtools::load_all()
library(fmrihrf)

set.seed(3301)

Tlen <- 120; V <- 3; r <- 4
sframe <- sampling_frame(blocklens = Tlen, TR = 1)

H <- cbind(
  exp(-seq(0, 30, length.out = Tlen)/4),
  exp(-seq(0, 30, length.out = Tlen)/6),
  exp(-seq(0, 30, length.out = Tlen)/8),
  exp(-seq(0, 30, length.out = Tlen)/10)
)
sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

onsets <- seq(10, 90, by = 10)
ntrials <- length(onsets)
design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

regs <- vector("list", ntrials)
for (t in seq_len(ntrials)) {
  rr_t <- fmrihrf::regressor(onsets = onsets[t], hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
  regs[[t]] <- fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
}

idx <- sample(seq_len(ncol(sbhm$A)), V, replace = TRUE)
alpha_true <- sbhm$A[, idx, drop = FALSE]
amps_true <- matrix(rnorm(ntrials * V, mean = 2, sd = 0.4), ntrials, V)

Y <- matrix(rnorm(Tlen * V, sd = 0.3), Tlen, V)
for (v in seq_len(V)) {
  for (t in seq_len(ntrials)) {
    Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(regs[[t]] %*% alpha_true[, v])
  }
}

# Run SBHM
res <- lss_sbhm(Y, sbhm, design_spec,
                prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
                return = "both")

# Check projection manually
beta_rt <- res$coeffs_r  # r x ntrials x V
alpha_hat <- res$alpha_coords  # r x V

cat("=== PROJECTION ANALYSIS ===\n")
cat("True amplitudes (first 3 trials, voxel 1):", amps_true[1:3, 1], "\n")
cat("Recovered amplitudes (first 3 trials, voxel 1):", res$amplitude[1:3, 1], "\n\n")

# Manual projection for voxel 1
manual_amps_v1 <- numeric(ntrials)
for (t in 1:ntrials) {
  manual_amps_v1[t] <- sum(beta_rt[, t, 1] * alpha_hat[, 1])
}
cat("Manual projection (before normalization):", manual_amps_v1[1:3], "\n")

# With normalization
den <- sum(alpha_hat[, 1]^2)
manual_amps_v1_norm <- manual_amps_v1 / den
cat("Manual projection (after normalization):", manual_amps_v1_norm[1:3], "\n")
cat("Matches sbhm_project output?", all.equal(manual_amps_v1_norm, res$amplitude[, 1]), "\n\n")

# Check if beta_rt actually recovers the true signal structure
cat("=== BETA_RT STRUCTURE ===\n")
cat("Expected: beta_rt[,t,v] should be proportional to amps_true[t,v] * alpha_true[,v]\n\n")

# For voxel 1, check if the trial-wise coefficients follow the pattern
for (t in 1:3) {
  expected <- amps_true[t, 1] * alpha_true[, 1]
  actual <- beta_rt[, t, 1]
  cat(sprintf("Trial %d, voxel 1:\n", t))
  cat("  Expected beta:", expected, "\n")
  cat("  Actual beta:  ", actual, "\n")
  cat("  Correlation:", cor(expected, actual), "\n\n")
}

# The key question: why is the projection giving wrong amplitudes?
# Let's check if alpha_hat is the problem or beta_rt
cat("=== DIAGNOSIS ===\n")
cat("alpha_hat[,1] vs alpha_true[,1] correlation:", cor(alpha_hat[,1], alpha_true[,1]), "\n")
cat("alpha_hat[,1] / ||alpha_hat||^2:", alpha_hat[,1] / sum(alpha_hat[,1]^2), "\n")
cat("alpha_true[,1] / ||alpha_true||^2:", alpha_true[,1] / sum(alpha_true[,1]^2), "\n\n")

# Check projection formula
# Under model: y_t = X_t(alpha) * amp_t + noise
# We fit: beta_t in basis space
# Should have: beta_t ≈ amp_t * alpha (in coefficient space)
# So: amp_t ≈ <beta_t, alpha> / ||alpha||^2
#
# But if beta_t is off by a scaling, this fails!
cat("Mean ||beta_rt|| per trial for voxel 1:", mean(apply(beta_rt[,,1], 2, function(x) sqrt(sum(x^2)))), "\n")
cat("Expected scale (mean amp * ||alpha||):", mean(amps_true[,1]) * sqrt(sum(alpha_true[,1]^2)), "\n")

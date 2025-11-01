#!/usr/bin/env Rscript
# Debug: check trial-wise coefficient signs

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
alpha <- sbhm$A[, idx, drop = FALSE]
amps_true <- matrix(rnorm(ntrials * V, mean = 2, sd = 0.4), ntrials, V)

Y <- matrix(rnorm(Tlen * V, sd = 0.3), Tlen, V)
for (v in seq_len(V)) {
  for (t in seq_len(ntrials)) {
    Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(regs[[t]] %*% alpha[, v])
  }
}

# Get matched alpha_hat
pre <- fmrilss:::sbhm_prepass(Y, sbhm, design_spec, ridge = list(mode = "fractional", lambda = 0.01))
m <- fmrilss:::sbhm_match(pre$beta_bar, sbhm$S, sbhm$A,
                          shrink = list(tau = 0, ref = sbhm$ref$alpha_ref, snr = NULL),
                          topK = 1, whiten = TRUE, orient_ref = TRUE)
alpha_hat <- m$alpha_hat

# Run OASIS to get trial-wise coefficients
hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
spec <- design_spec
spec$cond$hrf <- hrf_B
BetaMat <- lss(Y = Y, X = NULL, Z = NULL, Nuisance = NULL,
               method = "oasis", oasis = list(design_spec = spec, K = r),
               prewhiten = NULL)
beta_rt <- array(BetaMat, dim = c(r, ntrials, ncol(BetaMat)))
beta_mean <- apply(beta_rt, c(1, 3), mean)  # r × V

cat("=== ANALYSIS ===\n")
cat("True alpha (from library):\n")
print(alpha)
cat("\nMatched alpha_hat (after sbhm_match):\n")
print(alpha_hat)
cat("\nTrial-wise beta_mean (from OASIS LSS):\n")
print(beta_mean)

cat("\n=== DOT PRODUCTS ===\n")
dots <- colSums(beta_mean * alpha_hat)
cat("beta_mean · alpha_hat:", dots, "\n")
cat("Would flip:", which(dots < 0), "\n\n")

# The key insight: what should we compare to?
# Expected: amplitudes * alpha for each trial, averaged
# So beta_mean should be proportional to mean(amplitudes) * alpha
expected_mean <- colMeans(amps_true) * alpha
cat("Expected beta_mean (mean(amps) * alpha):\n")
print(expected_mean)
cat("\nCorrelation beta_mean vs expected:\n")
print(diag(cor(beta_mean, expected_mean)))

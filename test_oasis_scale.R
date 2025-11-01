#!/usr/bin/env Rscript
# Check OASIS scaling

devtools::load_all()
library(fmrihrf)

set.seed(3301)

Tlen <- 120; r <- 4
sframe <- sampling_frame(blocklens = Tlen, TR = 1)

H <- cbind(
  exp(-seq(0, 30, length.out = Tlen)/4),
  exp(-seq(0, 30, length.out = Tlen)/6),
  exp(-seq(0, 30, length.out = Tlen)/8),
  exp(-seq(0, 30, length.out = Tlen)/10)
)
sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

# Simple test: 1 voxel, 1 trial, known amplitude
onsets <- c(30)
amp_true <- 2.5
alpha_true <- sbhm$A[, 2]  # Pick library member 2

cat("=== SIMPLE OASIS TEST ===\n")
cat("True amplitude:", amp_true, "\n")
cat("True alpha:", alpha_true, "\n")
cat("||alpha||^2:", sum(alpha_true^2), "\n\n")

# Build regressor for this trial
rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
X_trial <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")

# Simulate Y = amp * X * alpha + noise
Y <- matrix(rnorm(Tlen, sd = 0.1), Tlen, 1)
Y[, 1] <- Y[, 1] + amp_true * as.numeric(X_trial %*% alpha_true)

cat("Signal strength in Y:", sd(as.numeric(X_trial %*% alpha_true)), "\n")
cat("Noise strength:", 0.1, "\n\n")

# Run OASIS LSS
design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
BetaMat <- lss(Y = Y, X = NULL, Z = NULL, Nuisance = NULL,
               method = "oasis",
               oasis = list(design_spec = design_spec, K = r),
               prewhiten = NULL)

cat("=== OASIS OUTPUT ===\n")
cat("BetaMat dimensions:", dim(BetaMat), "\n")
cat("BetaMat values:", BetaMat[, 1], "\n")
cat("||BetaMat||:", sqrt(sum(BetaMat[,1]^2)), "\n\n")

# Projection
proj <- sum(BetaMat[, 1] * alpha_true)
proj_norm <- proj / sum(alpha_true^2)

cat("Projection <beta, alpha>:", proj, "\n")
cat("Normalized projection:", proj_norm, "\n")
cat("Expected amplitude:", amp_true, "\n")
cat("Error:", proj_norm - amp_true, "\n")
cat("Ratio (recovered / true):", proj_norm / amp_true, "\n\n")

# Compare with direct GLM on the summate regressor
rr_agg <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = TRUE)
X_agg <- fmrihrf::evaluate(rr_agg, grid = sbhm$tgrid, precision = 0.1, method = "conv")

fit_agg <- lm(Y ~ X_agg - 1)
cat("=== AGGREGATE GLM (for comparison) ===\n")
cat("Coefficients from aggregate GLM:", coef(fit_agg), "\n")
cat("Projection <coef, alpha> / ||alpha||^2:", sum(coef(fit_agg) * alpha_true) / sum(alpha_true^2), "\n")

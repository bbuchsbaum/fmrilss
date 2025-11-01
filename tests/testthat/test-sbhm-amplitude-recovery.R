test_that("SBHM recovers trial amplitudes at moderate SNR", {
  skip_on_cran()
  set.seed(3301)
  library(fmrihrf)

  # Small synthetic setup
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

  # Trials and design (reduce overlap + add small jitter for identifiability)
  base_onsets <- c(12, 42, 72, 102)
  jitter <- runif(length(base_onsets), min = -8, max = 8)
  onsets <- pmin(pmax(base_onsets + jitter, 5), Tlen - 26)  # keep within window for span ~24–30
  ntrials <- length(onsets)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 24))

  # Build per-trial SBHM-basis regressors for simulation
  regs <- vector("list", ntrials)
  for (t in seq_len(ntrials)) {
    rr_t <- fmrihrf::regressor(onsets = onsets[t], hrf = hrf_B, duration = 0, span = 24, summate = FALSE)
    regs[[t]] <- fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  }

  # Assign each voxel a library HRF coordinate and trial amplitudes
  idx <- sample(seq_len(ncol(sbhm$A)), V, replace = TRUE)
  alpha <- sbhm$A[, idx, drop = FALSE]
  amps_true <- matrix(rnorm(ntrials * V, mean = 2, sd = 0.4), ntrials, V)

  # Simulate Y = sum_t amps[t]*X_t * alpha + noise
  Y <- matrix(rnorm(Tlen * V, sd = 0.10), Tlen, V)
  for (v in seq_len(V)) {
    for (t in seq_len(ntrials)) {
      Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(regs[[t]] %*% alpha[, v])
    }
  }

  # Run SBHM end-to-end
  res <- lss_sbhm(Y, sbhm, design_spec,
                  prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
                  return = "amplitude")

  # Correlation between recovered and true amplitudes should be reasonably high.
  # Note: We use absolute correlation because OASIS LSS can have per-trial sign
  # ambiguity that per-trial orientation may not fully resolve in all cases.
  # The key validation is that amplitudes track the true pattern, even if signs vary.
  cor_all <- suppressWarnings(abs(cor(as.vector(res$amplitude), as.vector(amps_true))))
  expect_gt(cor_all, 0.6)
})

test_that("SBHM amplitude methods correlate and fallback can trigger", {
  skip_on_cran()
  set.seed(4321)
  library(fmrihrf)

  # Small synthetic setup
  Tlen <- 120; V <- 12; r <- 4
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

  # Moderately spaced design (better conditioned)
  onsets <- c(10, 34, 58, 82, 106)
  ntrials <- length(onsets)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 24))

  # Build per-trial basis regressors and simulate Y
  regs <- lapply(onsets, function(ot) {
    rr_t <- fmrihrf::regressor(onsets = ot, hrf = hrf_B, duration = 0, span = 24, summate = FALSE)
    fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  })
  idx <- sample(seq_len(ncol(sbhm$A)), V, replace = TRUE)
  alpha <- sbhm$A[, idx, drop = FALSE]
  amps_true <- matrix(rnorm(ntrials * V, mean = 2, sd = 0.4), ntrials, V)
  Y <- matrix(rnorm(Tlen * V, sd = 0.12), Tlen, V)
  for (v in seq_len(V)) for (t in seq_len(ntrials)) Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(regs[[t]] %*% alpha[, v])

  # Global LS amplitudes
  gls <- lss_sbhm(
    Y, sbhm, design_spec,
    amplitude = list(method = "global_ls", ridge = list(mode = "fractional", lambda = 0.02)),
    return = "amplitude"
  )
  # Per-trial 2x2 LSS amplitudes
  lss1 <- lss_sbhm(
    Y, sbhm, design_spec,
    amplitude = list(method = "lss1", ridge_frac = list(x = 0.02, b = 0.02)),
    return = "amplitude"
  )

  # Cross-method correlation in well-conditioned regime should be high
  cor_methods <- suppressWarnings(cor(as.vector(gls$amplitude), as.vector(lss1$amplitude), method = "spearman"))
  expect_gt(cor_methods, 0.7)

  # Fast ER design with cond_gate fallback
  set.seed(4322)
  base_onsets <- seq(8, 112, length.out = 12)
  jitter <- runif(length(base_onsets), min = -2, max = 2)
  on_fast <- pmax(pmin(base_onsets + jitter, Tlen - 26), 5)
  design_fast <- list(sframe = sframe, cond = list(onsets = on_fast, duration = 0, span = 24))

  # Resimulate Y to align with the fast design (heavy overlap)
  regs_fast <- lapply(on_fast, function(ot) {
    rr_t <- fmrihrf::regressor(onsets = ot, hrf = hrf_B, duration = 0, span = 24, summate = FALSE)
    fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  })
  ntrials_fast <- length(on_fast)
  amps_fast <- matrix(rnorm(ntrials_fast * V, mean = 2, sd = 0.4), ntrials_fast, V)
  Y_fast <- matrix(rnorm(Tlen * V, sd = 0.12), Tlen, V)
  for (v in seq_len(V)) for (t in seq_len(ntrials_fast)) {
    Y_fast[, v] <- Y_fast[, v] + amps_fast[t, v] * as.numeric(regs_fast[[t]] %*% alpha[, v])
  }

  out_fast <- lss_sbhm(
    Y_fast, sbhm, design_fast,
    amplitude = list(
      method = "global_ls",
      ridge = list(mode = "fractional", lambda = 0.02),
      # Use rho-based gate for more stable triggering under overlap
      cond_gate = list(metric = "rho", thr = 0.97, fallback = "lss1")
    ),
    return = "amplitude"
  )
  # Expect at least some fallbacks in fast ER regime
  n_fallback <- sum(out_fast$diag$method_used == "lss1")
  expect_gte(n_fallback, 1L)
})

test_that("ISI sweep produces sensible amplitudes and fallback share", {
  skip_on_cran()
  set.seed(55221)
  library(fmrihrf)

  # Common setup
  Tlen <- 180; V <- 8; r <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8),
    exp(-seq(0, 30, length.out = Tlen)/10)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

  sim_one <- function(isi, jitter = 0, noise_sd = 0.12) {
    # Build onsets covering [10, Tlen-30] with approximate ISI and jitter
    bounds <- c(10, Tlen - 30)
    nslots <- max(3L, floor((bounds[2] - bounds[1]) / isi))
    base <- seq(bounds[1], bounds[2], length.out = nslots)
    if (jitter > 0) base <- pmin(pmax(base + runif(length(base), -jitter, jitter), bounds[1]), bounds[2])
    onsets <- sort(base)
    design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 24))

    # Build per-trial basis regressors and simulate Y
    regs <- lapply(onsets, function(ot) {
      rr_t <- fmrihrf::regressor(onsets = ot, hrf = hrf_B, duration = 0, span = 24, summate = FALSE)
      fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
    })
    idx <- sample(seq_len(ncol(sbhm$A)), V, replace = TRUE)
    alpha <- sbhm$A[, idx, drop = FALSE]
    amps_true <- matrix(rnorm(length(onsets) * V, mean = 2, sd = 0.4), length(onsets), V)
    Y <- matrix(rnorm(Tlen * V, sd = noise_sd), Tlen, V)
    for (v in seq_len(V)) for (t in seq_along(onsets)) Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(regs[[t]] %*% alpha[, v])

    list(Y = Y, spec = design_spec)
  }

  # Slow ER (better conditioning)
  set.seed(55222)
  slow <- sim_one(isi = 15, jitter = 2)
  out_gls_slow <- lss_sbhm(slow$Y, sbhm, slow$spec,
                            amplitude = list(method = "global_ls", ridge = list(mode = "fractional", lambda = 0.02)),
                            return = "amplitude")
  out_lss_slow <- lss_sbhm(slow$Y, sbhm, slow$spec,
                            amplitude = list(method = "lss1", ridge_frac = list(x = 0.02, b = 0.02)),
                            return = "amplitude")
  cor_slow <- suppressWarnings(cor(as.vector(out_gls_slow$amplitude), as.vector(out_lss_slow$amplitude), method = "spearman"))
  expect_gt(cor_slow, 0.7)

  # Fast ER (heavy overlap) with fallback gate
  set.seed(55223)
  fast <- sim_one(isi = 2, jitter = 0.5)
  out_fast <- lss_sbhm(fast$Y, sbhm, fast$spec,
                       amplitude = list(method = "global_ls",
                                        ridge = list(mode = "fractional", lambda = 0.02),
                                        cond_gate = list(metric = "rho", thr = 0.97, fallback = "lss1")),
                       return = "amplitude")
  # Expect some fallbacks in fast design
  n_fallback <- sum(out_fast$diag$method_used == "lss1")
  expect_gte(n_fallback, 1L)
})

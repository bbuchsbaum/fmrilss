test_that("sbhm_prepass recovers coefficients with small ridge", {
  set.seed(1001)
  library(fmrihrf)

  Tlen <- 120; V <- 1; r <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  # Simple library matrix; columns roughly HRF-like decays
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8),
    exp(-seq(0, 30, length.out = Tlen)/10)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)

  onsets <- seq(5, 95, by = 10)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

  # Build aggregate per-basis regressor for signal simulation
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = TRUE)
  Xagg <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")

  alpha_true <- rnorm(r)
  Y <- matrix(rnorm(Tlen*V, sd = 0.3), Tlen, V)
  Y[,1] <- Y[,1] + Xagg %*% alpha_true

  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
  pre <- sbhm_prepass(Y, sbhm, design_spec, ridge = list(mode = "absolute", lambda = 0.01))

  expect_equal(nrow(pre$beta_bar), r)
  expect_equal(ncol(pre$beta_bar), V)
  # Correlation should be reasonably high with modest ridge
  expect_gt(cor(pre$beta_bar[,1], alpha_true), 0.5)
})

test_that("Aggregated per-basis design equals fmrihrf summate=TRUE", {
  set.seed(1002)
  library(fmrihrf)

  Tlen <- 90; r <- 3
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  onsets <- seq(6, 72, by = 12)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  spec  <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30, hrf = hrf_B))

  built <- fmrilss:::.oasis_build_X_from_events(spec)
  X_trials <- built$X_trials
  idx_by_basis <- lapply(seq_len(r), function(k) seq.int(k, ncol(X_trials), by = r))
  A <- do.call(cbind, lapply(idx_by_basis, function(idx) rowSums(X_trials[, idx, drop = FALSE])))

  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = TRUE)
  Xagg <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")

  expect_equal(dim(A), dim(Xagg))
  expect_lt(max(abs(A - Xagg)), 1e-10)
})

test_that("Factorized path matches dense path (no prewhitening)", {
  set.seed(1003)
  library(fmrihrf)

  Tlen <- 80; V <- 5; r <- 3; q <- 6
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 24, length.out = Tlen)/4),
    exp(-seq(0, 24, length.out = Tlen)/6),
    exp(-seq(0, 24, length.out = Tlen)/8)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)

  onsets <- seq(5, 65, by = 10)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 24))

  # Create exactly factorized data Y = Scores %*% Loadings
  Scores <- matrix(rnorm(Tlen*q), Tlen, q)
  Load   <- matrix(rnorm(q*V), q, V)
  Y <- Scores %*% Load

  # Dense path
  pre_dense <- sbhm_prepass(Y, sbhm, design_spec, prewhiten = list(method = "none"))
  # Factorized path
  pre_fac   <- sbhm_prepass(Y, sbhm, design_spec, data_fac = list(scores = Scores, loadings = Load))

  expect_equal(dim(pre_dense$beta_bar), dim(pre_fac$beta_bar))
  expect_equal(pre_dense$diag$K, pre_fac$diag$K)
  expect_lt(max(abs(pre_dense$beta_bar - pre_fac$beta_bar)), 1e-8)
})

test_that("Including design_spec$others as nuisances reduces bias", {
  set.seed(1004)
  library(fmrihrf)

  Tlen <- 120; V <- 1; r <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8),
    exp(-seq(0, 30, length.out = Tlen)/10)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

  on_main <- seq(8, 110, by = 14)
  on_oth  <- seq(4, 104, by = 20)

  # Build A (aggregated basis) for main condition
  spec_main <- list(sframe = sframe, cond = list(onsets = on_main, duration = 0, span = 30, hrf = hrf_B))
  built_main <- fmrilss:::.oasis_build_X_from_events(spec_main)
  X_trials <- built_main$X_trials
  idx_by_basis <- lapply(seq_len(r), function(k) seq.int(k, ncol(X_trials), by = r))
  A <- do.call(cbind, lapply(idx_by_basis, function(idx) rowSums(X_trials[, idx, drop = FALSE])))

  # Build X_other for nuisance condition (aggregated to single column)
  spec_both <- list(
    sframe = sframe,
    cond   = list(onsets = on_main, duration = 0, span = 30, hrf = hrf_B),
    others = list(list(onsets = on_oth, duration = 0, span = 30, hrf = hrf_B))
  )
  built_both <- fmrilss:::.oasis_build_X_from_events(spec_both)
  X_other <- built_both$X_other

  # Simulate signal: main basis coefficients + nuisance condition leak + small noise
  alpha_true <- rnorm(r)
  gamma <- 2.0
  Y <- matrix(rnorm(Tlen*V, sd = 0.15), Tlen, V)
  Y[,1] <- Y[,1] + A %*% alpha_true + as.numeric(X_other) * gamma

  # Prepass without others as nuisance
  pre_no_oth <- sbhm_prepass(Y, sbhm, design_spec = spec_main, ridge = list(mode = "absolute", lambda = 0.1))
  # Prepass with others included as nuisance
  pre_with_oth <- sbhm_prepass(Y, sbhm, design_spec = spec_both, ridge = list(mode = "absolute", lambda = 0.1))

  rmse <- function(b) sqrt(mean((b - alpha_true)^2))
  err_no <- rmse(pre_no_oth$beta_bar[,1])
  err_w  <- rmse(pre_with_oth$beta_bar[,1])

  expect_lt(err_w, err_no)
})

test_that("Prewhitening path runs and sets flag", {
  set.seed(1005)
  library(fmrihrf)

  Tlen <- 100; V <- 2; r <- 3
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 25, length.out = Tlen)/4),
    exp(-seq(0, 25, length.out = Tlen)/6),
    exp(-seq(0, 25, length.out = Tlen)/8)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  onsets <- seq(6, 80, by = 12)
  spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 25, hrf = hrf_B))

  # Simulate AR(1) noise + tiny signal
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 25, summate = TRUE)
  Xagg <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  alpha_true <- rnorm(r, sd = 0.2)

  make_ar1 <- function(n, phi = 0.7, sd = 0.5) {
    e <- rnorm(n, sd = sd)
    y <- numeric(n)
    for (t in 2:n) y[t] <- phi * y[t-1] + e[t]
    y
  }

  Y <- cbind(Xagg %*% alpha_true + make_ar1(Tlen),
             Xagg %*% rnorm(r, sd = 0.2) + make_ar1(Tlen))

  pre_w <- sbhm_prepass(Y, sbhm, spec, prewhiten = list(method = "ar", p = 1L, pooling = "global", exact_first = "ar1"))

  expect_true(isTRUE(pre_w$diag$used_prewhiten))
  expect_equal(nrow(pre_w$beta_bar), r)
  expect_equal(ncol(pre_w$beta_bar), ncol(Y))
})

test_that("Small ridge reduces error vs no ridge in ill-conditioned design", {
  set.seed(1006)
  library(fmrihrf)

  Tlen <- 120; V <- 1; r <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8),
    exp(-seq(0, 30, length.out = Tlen)/10)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

  onsets <- seq(5, 95, by = 10)
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = TRUE)
  Xagg <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")

  alpha_true <- rnorm(r)
  Y <- matrix(rnorm(Tlen*V, sd = 0.3), Tlen, V)
  Y[,1] <- Y[,1] + Xagg %*% alpha_true

  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
  pre0 <- sbhm_prepass(Y, sbhm, design_spec, ridge = list(mode = "absolute", lambda = 0))
  pre1 <- sbhm_prepass(Y, sbhm, design_spec, ridge = list(mode = "absolute", lambda = 0.01))

  rmse <- function(b) sqrt(mean((b - alpha_true)^2))
  expect_lt(rmse(pre1$beta_bar[,1]), rmse(pre0$beta_bar[,1]))
})

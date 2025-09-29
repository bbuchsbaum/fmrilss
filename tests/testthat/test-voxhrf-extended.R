test_that("VOXHRF with AR1 whitening returns finite betas", {
  skip_on_cran()
  set.seed(101)
  n_time <- 140
  n_vox  <- 16
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  onsets <- seq(8, n_time - 24, by = 12)

  beta <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30)
      ),
      hrf_mode = "voxel_ridge",
      whiten = "ar1",
      lambda_shape = 0
    )
  )

  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(length(onsets), n_vox))
  expect_true(all(is.finite(beta)))
})

test_that("VOXHRF runs with design_spec$others and nuisance", {
  skip_on_cran()
  set.seed(102)
  n_time <- 160
  n_vox  <- 12
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  onsets_main <- seq(10, n_time - 30, by = 15)
  onsets_other <- seq(5, n_time - 35, by = 20)
  Nuis <- cbind(1, scale(seq_len(n_time)))

  beta <- lss(
    Y = Y,
    X = NULL,
    Z = NULL,
    Nuisance = Nuis,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(onsets = onsets_main, hrf = fmrihrf::HRF_SPMG3, span = 30),
        others = list(list(onsets = onsets_other))
      ),
      hrf_mode = "voxel_ridge",
      lambda_shape = 1
    )
  )

  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(length(onsets_main), n_vox))
  expect_true(all(is.finite(beta)))
})

test_that("estimate_voxel_hrf_fast returns KxV and unit-energy shapes", {
  skip_on_cran()
  set.seed(103)
  n_time <- 120
  n_vox  <- 10
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)
  onsets <- seq(8, n_time - 24, by = 10)

  spec <- list(
    sframe = sframe,
    cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30),
    precision = 0.1, method = "conv"
  )
  built <- fmrilss:::`.oasis_build_X_from_events`(spec)
  X_trials <- built$X_trials
  K <- built$K

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  vhrf <- fmrilss:::`.estimate_voxel_hrf_fast`(
    Y = Y, X_trials = X_trials, design_spec = spec, N_nuis = NULL, K = NULL,
    lambda_shape = 0, mu_rough = 0, ref_hrf = NULL, shrink_global = 0, orient_ref = TRUE
  )

  expect_s3_class(vhrf, "VoxelHRF")
  expect_equal(dim(vhrf$coefficients), c(K, n_vox))

  # Check unit-energy normalization on the time grid
  r <- fmrihrf::regressor(onsets = 0, hrf = spec$cond$hrf, duration = 0, span = 30)
  B_time <- fmrihrf::evaluate(r, grid = times, precision = 0.1, method = "conv")
  if (inherits(B_time, "Matrix")) B_time <- as.matrix(B_time)
  if (is.vector(B_time)) B_time <- cbind(B_time)
  H <- B_time %*% vhrf$coefficients
  energy <- sqrt(colSums(H * H))
  expect_true(max(abs(energy - 1)) < 1e-6)
})

test_that("Roughness penalty reduces HRF curvature energy", {
  skip_on_cran()
  set.seed(104)
  n_time <- 150
  n_vox  <- 8
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)
  onsets <- seq(12, n_time - 24, by = 12)

  spec <- list(
    sframe = sframe,
    cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30),
    precision = 0.1, method = "conv"
  )
  built <- fmrilss:::`.oasis_build_X_from_events`(spec)
  X_trials <- built$X_trials

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  v0 <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, X_trials, spec, NULL, NULL, 0, 0, NULL, 0, TRUE)
  v1 <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, X_trials, spec, NULL, NULL, 0, 0.5, NULL, 0, TRUE)

  r <- fmrihrf::regressor(onsets = 0, hrf = spec$cond$hrf, duration = 0, span = 30)
  B_time <- fmrihrf::evaluate(r, grid = times, precision = 0.1, method = "conv")
  if (inherits(B_time, "Matrix")) B_time <- as.matrix(B_time)
  if (is.vector(B_time)) B_time <- cbind(B_time)
  H0 <- B_time %*% v0$coefficients
  H1 <- B_time %*% v1$coefficients

  # Discrete 2nd diff curvature energy
  D2 <- diff(diag(nrow(H0)), differences = 2)
  curv0 <- colSums((D2 %*% H0)^2)
  curv1 <- colSums((D2 %*% H1)^2)
  expect_true(mean(curv1) <= mean(curv0) + 1e-8)
})

test_that("Global shrinkage lowers across-voxel weight variance", {
  skip_on_cran()
  set.seed(105)
  n_time <- 120
  n_vox  <- 20
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  onsets <- seq(8, n_time - 24, by = 10)
  spec <- list(sframe = sframe, cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30))
  built <- fmrilss:::`.oasis_build_X_from_events`(spec)

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  v0 <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, built$X_trials, spec, NULL, NULL, 0, 0, NULL, 0, TRUE)
  v1 <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, built$X_trials, spec, NULL, NULL, 0, 0, NULL, 0.1, TRUE)

  # Row-wise variance across voxels should decrease with shrinkage
  rowvar <- function(M) apply(M, 1, var)
  expect_true(mean(rowvar(v1$coefficients)) < mean(rowvar(v0$coefficients)))
})

test_that("Orientation flip aligns shapes with reference", {
  skip_on_cran()
  set.seed(106)
  n_time <- 100
  n_vox  <- 10
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)
  onsets <- seq(8, n_time - 24, by = 10)
  spec <- list(sframe = sframe, cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG1, span = 25))
  built <- fmrilss:::`.oasis_build_X_from_events`(spec)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Make ref negative of canonical to exercise flip
  r <- fmrihrf::regressor(onsets = 0, hrf = spec$cond$hrf, duration = 0, span = 25)
  B_time <- fmrihrf::evaluate(r, grid = times, precision = 0.1, method = "conv")
  if (is.vector(B_time)) B_time <- cbind(B_time)
  ref <- -as.numeric(B_time[, 1])

  v <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, built$X_trials, spec, NULL, NULL, 0, 0, ref, 0, TRUE)
  H <- B_time %*% v$coefficients
  dots <- colSums(H * ref)
  expect_true(all(dots >= 0))
})

test_that("Custom ref_hrf increases alignment as lambda_shape increases", {
  skip_on_cran()
  set.seed(107)
  n_time <- 140
  n_vox  <- 12
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)
  onsets <- seq(12, n_time - 24, by = 12)
  spec <- list(sframe = sframe, cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30))
  built <- fmrilss:::`.oasis_build_X_from_events`(spec)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Create a shifted reference shape on the time grid
  r <- fmrihrf::regressor(onsets = 0, hrf = fmrihrf::HRF_SPMG1, duration = 0, span = 30)
  href <- as.numeric(fmrihrf::evaluate(r, grid = times, precision = 0.1, method = "conv"))
  href_shift <- c(rep(0, 1), href[-length(href)])

  v0 <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, built$X_trials, spec, NULL, NULL, 0, 0, href_shift, 0, TRUE)
  v1 <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, built$X_trials, spec, NULL, NULL, 5, 0, href_shift, 0, TRUE)

  # Compare mean shape across voxels
  B3 <- {
    rr <- fmrihrf::regressor(onsets = 0, hrf = spec$cond$hrf, duration = 0, span = 30)
    Bt <- fmrihrf::evaluate(rr, grid = times, precision = 0.1, method = "conv")
    if (inherits(Bt, "Matrix")) Bt <- as.matrix(Bt)
    if (is.vector(Bt)) Bt <- cbind(Bt)
    Bt
  }
  H0 <- B3 %*% v0$coefficients
  H1 <- B3 %*% v1$coefficients
  m0 <- rowMeans(H0); m1 <- rowMeans(H1)
  cor0 <- abs(cor(m0, href_shift))
  cor1 <- abs(cor(m1, href_shift))
  expect_true(cor1 >= cor0 - 1e-8)
})

test_that("K autodetect matches HRF nbasis for voxel HRF fit", {
  skip_on_cran()
  set.seed(108)
  n_time <- 100
  n_vox  <- 6
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  onsets <- seq(10, n_time - 20, by = 15)
  spec <- list(sframe = sframe, cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 25))
  built <- fmrilss:::`.oasis_build_X_from_events`(spec)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  v <- fmrilss:::`.estimate_voxel_hrf_fast`(Y, built$X_trials, spec, NULL, NULL, 0, 0, NULL, 0, TRUE)
  expect_equal(nrow(v$coefficients), fmrihrf::nbasis(spec$cond$hrf))
})

test_that("VOXHRF works with advanced prewhitening (fmriAR) when available", {
  skip_on_cran()
  skip_if_not_installed("fmriAR")
  set.seed(109)
  n_time <- 120
  n_vox  <- 10
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  onsets <- seq(8, n_time - 24, by = 10)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  beta <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    prewhiten = list(method = "arma", p = 1, q = 0),
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30)
      ),
      hrf_mode = "voxel_ridge",
      lambda_shape = 0
    )
  )

  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(length(onsets), n_vox))
  expect_true(all(is.finite(beta)))
})

test_that("VOXHRF handles degenerate trial designs without crashing", {
  skip_on_cran()
  set.seed(404)
  n_time <- 80
  n_vox  <- 3
  onsets <- c(10, 10, 40, 40)  # repeated events -> collinear Xi/Xother
  durations <- c(3, 3, 2, 2)
  events <- data.frame(onset = onsets, duration = durations, condition = "cond")

  basis <- fmrihrf::HRF_SPMG1
  coeffs <- matrix(1, nrow = 1L, ncol = n_vox)
  vhrf <- list(coefficients = coeffs, basis = basis, conditions = "cond")
  class(vhrf) <- "VoxelHRF"

  Y <- matrix(rnorm(n_time * n_vox, sd = 0.1), n_time, n_vox)

  expect_error(
    {
      out <- lss_with_hrf(
        Y = Y,
        events = events,
        hrf_estimates = vhrf,
        nuisance_regs = NULL,
        engine = "C++",
        verbose = FALSE
      )
      mat <- as.matrix(out$betas[])
      expect_equal(nrow(mat), length(onsets))
      expect_equal(ncol(mat), n_vox)
    },
    NA
  )
})

test_that("Manual single AR1 whiten == internal VOXHRF AR1 whiten", {
  skip_on_cran()
  set.seed(110)
  n_time <- 150
  n_vox  <- 18
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  onsets <- seq(10, n_time - 25, by = 15)
  spec <- list(sframe = sframe, cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30))

  built <- fmrilss:::`.oasis_build_X_from_events`(spec)
  X_trials <- built$X_trials
  K <- fmrihrf::nbasis(spec$cond$hrf)

  # Random data
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Manual AR1 whiten once for both Y and X
  # Internal path uses default intercept in nuisance; emulate that here
  Z <- matrix(1, n_time, 1)
  w <- fmrilss:::`.oasis_ar1_whitener`(Y, X_nuis = Z)
  Yw <- w$Wy
  Xw <- w$W_apply(X_trials)
  Zw <- w$W_apply(Z)

  # Internal path: whiten inside with oasis$whiten = "ar1"
  beta_internal <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = spec,
      hrf_mode = "voxel_ridge",
      whiten = "ar1",
      K = K,
      lambda_shape = 0
    )
  )

  # Manual path: pass whitened Y and X; disable whitening inside
  beta_manual <- lss(
    Y = Yw,
    X = Xw,
    Z = Zw,
    method = "oasis",
    oasis = list(
      design_spec = spec,
      hrf_mode = "voxel_ridge",
      whiten = "none",
      K = K,
      lambda_shape = 0
    )
  )

  # Compare within tight tolerance
  expect_equal(beta_manual, beta_internal, tolerance = 1e-8)
})

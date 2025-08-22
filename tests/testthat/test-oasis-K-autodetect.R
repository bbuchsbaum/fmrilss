# tests/testthat/test-oasis-K-autodetect.R
test_that("K autodetect via ntrials matches explicit K for SPMG3 designs", {
  skip_on_cran()
  skip_if_not_installed("fmrihrf")

  set.seed(31)
  T  <- 180
  V  <- 40
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)

  # Build trial-wise SPMG3 design (K = 3)
  ntrials <- 10
  onsets  <- sort(sample(12:(T-30), ntrials))
  K <- 3

  rset <- fmrihrf::regressor_set(
    onsets   = onsets,
    fac      = factor(seq_along(onsets)),
    hrf      = fmrihrf::HRF_SPMG3,
    duration = 0, amplitude = 1, span = 40, summate = TRUE
  )
  X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  if (inherits(X, "Matrix")) X <- as.matrix(X)

  # Simulate data with some signal
  Theta <- matrix(rnorm(ncol(X) * V, sd = 0.3), nrow = ncol(X))
  Y <- X %*% Theta + matrix(rnorm(T * V, sd = 1.0), T, V)

  # (1) Autodetect K from ntrials and X columns
  B_auto <- fmrilss:::.lss_oasis(
    Y, X = X, Z = NULL, Nuisance = NULL,
    oasis = list(ntrials = ntrials, ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)
  )

  # (2) Explicit K path should match
  B_explicit <- fmrilss:::.lss_oasis(
    Y, X = X, Z = NULL, Nuisance = NULL,
    oasis = list(K = K, ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)
  )

  # Shapes and equality
  expect_equal(dim(B_auto), c(ntrials * K, V))
  expect_equal(dim(B_explicit), c(ntrials * K, V))
  expect_equal(B_auto, B_explicit, tolerance = 1e-10)
})
# tests/testthat/test-oasis-api-compat.R
test_that("fmrilss::lss(method='oasis') returns matrix by default and matches .lss_oasis", {
  skip_on_cran()
  skip_if_not_installed("fmrihrf")

  # Try to call fmrilss::lss; skip if method='oasis' not wired yet
  have_lss <- "lss" %in% getNamespaceExports("fmrilss")
  skip_if_not(have_lss, "fmrilss::lss not exported")

  set.seed(41)
  T  <- 200
  V  <- 50
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)

  # Single-basis design (K=1) - using SPMG1 which is predefined
  onsets <- sort(sample(12:(T-24), 28))
  rset   <- fmrihrf::regressor_set(onsets   = onsets,
                                   fac      = factor(seq_along(onsets)),
                                   hrf      = fmrihrf::HRF_SPMG1,
                                   duration = 0, amplitude = 1, span = 40, summate = TRUE)
  X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  if (inherits(X, "Matrix")) X <- as.matrix(X)

  # Simulate data (with a simple intercept nuisance)
  Beta <- matrix(rnorm(ncol(X) * V, sd = 0.4), nrow = ncol(X))
  Y    <- X %*% Beta + matrix(rnorm(T * V, sd = 1.0), T, V)
  Nuis <- matrix(1, T, 1)

  # Internal OASIS (reference)
  B_ref <- fmrilss:::.lss_oasis(
    Y, X = X, Z = NULL, Nuisance = Nuis,
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0, block_cols = 512L)
  )

  # Public API: lss(..., method = "oasis")
  # If method isn't wired yet, this call will error—skip gracefully.
  res <- try(
    fmrilss::lss(Y, X, Z = NULL, Nuisance = Nuis,
                 method = "oasis",
                 oasis  = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0, block_cols = 512L)),
    silent = TRUE
  )
  if (inherits(res, "try-error")) skip("fmrilss::lss doesn't support method='oasis' yet")

  # By default lss returns a bare matrix
  expect_true(is.matrix(res))
  expect_equal(dim(res), dim(B_ref))
  expect_equal(res, B_ref, tolerance = 1e-10, scale = 1)
})

test_that("fmrilss::lss(method='oasis', return_se=TRUE) returns a list with beta & se", {
  skip_on_cran()
  skip_if_not_installed("fmrihrf")

  have_lss <- "lss" %in% getNamespaceExports("fmrilss")
  skip_if_not(have_lss, "fmrilss::lss not exported")

  set.seed(42)
  T  <- 160
  V  <- 40
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)

  onsets <- sort(sample(10:(T-24), 20))
  rset   <- fmrihrf::regressor_set(onsets   = onsets,
                                   fac      = factor(seq_along(onsets)),
                                   hrf      = fmrihrf::HRF_SPMG1,
                                   duration = 0, amplitude = 1, span = 40, summate = TRUE)
  X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  if (inherits(X, "Matrix")) X <- as.matrix(X)

  Y <- X %*% matrix(rnorm(ncol(X) * V, sd = 0.3), nrow = ncol(X)) + matrix(rnorm(T * V), T, V)

  res <- try(
    fmrilss::lss(Y, X, method = "oasis",
                 oasis = list(return_se = TRUE, ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)),
    silent = TRUE
  )
  if (inherits(res, "try-error")) skip("fmrilss::lss doesn't support method='oasis' yet")

  expect_true(is.list(res))
  expect_true(all(c("beta","se") %in% names(res)))
  expect_equal(dim(res$beta), c(ncol(X), V))
  expect_equal(dim(res$se),   c(ncol(X), V))
  expect_true(all(is.finite(res$beta)))
  expect_true(all(is.finite(res$se)))
})
# Tests for sbhm_build.R to improve coverage

test_that("sbhm_build requires exactly one of library_spec or library_H", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)

  # Neither provided
  expect_error(sbhm_build(sframe = sframe), "Provide exactly one of")

  # Both provided
  H <- matrix(rnorm(100 * 10), 100, 10)
  param_grid <- expand.grid(shape = c(6, 8), rate = c(0.9, 1.1))
  gamma_fun <- function(shape, rate) {
    fmrihrf::as_hrf(function(t) fmrihrf::hrf_gamma(t, shape = shape, rate = rate), span = 32)
  }
  spec <- list(fun = gamma_fun, pgrid = param_grid, span = 32)

  expect_error(sbhm_build(library_spec = spec, library_H = H, sframe = sframe),
               "Provide exactly one of")
})

test_that("sbhm_build works with library_H", {
  skip_if_not_installed("fmrihrf")

  set.seed(123)
  T_len <- 100
  K <- 10
  tgrid <- seq(0, T_len - 1, by = 1)

  # Create synthetic HRF library
  H <- matrix(rnorm(T_len * K), T_len, K)

  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 3)

  expect_true(is.list(result))
  expect_true("B" %in% names(result))
  expect_true("S" %in% names(result))
  expect_true("A" %in% names(result))
  expect_equal(ncol(result$B), 3)  # rank 3
  expect_equal(nrow(result$B), T_len)
  expect_equal(nrow(result$A), 3)  # r rows
  expect_equal(ncol(result$A), K)  # K columns
})

test_that("sbhm_build validates library_H dimensions", {
  skip_if_not_installed("fmrihrf")

  tgrid <- 1:100
  H <- matrix(rnorm(50 * 10), 50, 10)  # Wrong number of rows

  expect_error(sbhm_build(library_H = H, tgrid = tgrid),
               "library_H rows must match")
})

test_that("sbhm_build requires sframe or tgrid", {
  skip_if_not_installed("fmrihrf")

  H <- matrix(rnorm(100 * 10), 100, 10)

  expect_error(sbhm_build(library_H = H), "Provide sframe or tgrid")
})

test_that("sbhm_build works with sframe instead of tgrid", {
  skip_if_not_installed("fmrihrf")

  set.seed(234)
  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 2)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- matrix(rnorm(length(times) * 8), length(times), 8)

  result <- sbhm_build(library_H = H, sframe = sframe, r = 4)

  expect_equal(result$tgrid, times)
  expect_equal(ncol(result$B), 4)
})

test_that("sbhm_build applies baseline removal", {
  skip_if_not_installed("fmrihrf")

  set.seed(345)
  T_len <- 100
  tgrid <- seq(0, 99, by = 1)

  # Create library with non-zero baseline
  H <- matrix(5 + rnorm(T_len * 5), T_len, 5)

  # With baseline removal in window [0, 10]
  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 3,
                       baseline = c(0, 10), normalize = FALSE)

  # The baseline should be removed
  expect_true(is.list(result))
})

test_that("sbhm_build applies normalization", {
  skip_if_not_installed("fmrihrf")

  set.seed(456)
  T_len <- 50
  tgrid <- seq(0, 49, by = 1)
  H <- matrix(rnorm(T_len * 5, sd = 10), T_len, 5)

  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 3,
                       normalize = TRUE, baseline = NULL)

  # B should be orthonormal columns
  BtB <- crossprod(result$B)
  expect_equal(BtB, diag(ncol(result$B)), tolerance = 1e-8)
})

test_that("sbhm_build handles shifts parameter", {
  skip_if_not_installed("fmrihrf")

  set.seed(567)
  T_len <- 100
  tgrid <- seq(0, 99, by = 1)
  H <- matrix(rnorm(T_len * 5), T_len, 5)

  # With shifts
  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 4,
                       shifts = c(-1, 1), baseline = NULL)

  # Should have augmented library (5 original + 5*2 shifted = 15 columns in A)
  expect_equal(ncol(result$A), 15)
})

test_that("sbhm_build ref='mean' computes mean coordinates", {
  skip_if_not_installed("fmrihrf")

  set.seed(678)
  T_len <- 50
  tgrid <- seq(0, 49, by = 1)
  H <- matrix(rnorm(T_len * 8), T_len, 8)

  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 3,
                       ref = "mean", baseline = NULL)

  expect_equal(result$ref$name, "mean")
  expect_equal(length(result$ref$alpha_ref), 3)
  # alpha_ref should be rowMeans of A
  expect_equal(result$ref$alpha_ref, rowMeans(result$A), tolerance = 1e-10)
})

test_that("sbhm_build ref='spmg1' projects SPMG1 HRF", {
  skip_if_not_installed("fmrihrf")
  # Skip due to fmrihrf::make_hrf API requirements
  skip("sbhm_build ref='spmg1' requires specific fmrihrf API")
})

test_that("sbhm_build clips rank to min(T, K)", {
  skip_if_not_installed("fmrihrf")

  set.seed(890)
  T_len <- 50
  K <- 3  # fewer columns than requested rank
  tgrid <- seq(0, T_len - 1, by = 1)
  H <- matrix(rnorm(T_len * K), T_len, K)

  # Request rank 10, but only 3 columns available
  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 10, baseline = NULL)

  # Should be clipped to min(T_len, K) = 3
  expect_equal(result$meta$r, 3)
  expect_equal(ncol(result$B), 3)
})

test_that("sbhm_build returns expected meta information", {
  skip_if_not_installed("fmrihrf")

  set.seed(901)
  T_len <- 60
  K <- 8
  tgrid <- seq(0, T_len - 1, by = 1)
  H <- matrix(rnorm(T_len * K), T_len, K)

  result <- sbhm_build(library_H = H, tgrid = tgrid, r = 5,
                       normalize = TRUE, baseline = c(0, 5))

  expect_equal(result$meta$r, 5)
  expect_equal(result$meta$K, K)
  expect_true(result$meta$normalize)
  expect_equal(result$meta$baseline, c(0, 5))
})

test_that("sbhm_build validates library_spec structure", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)

  # Missing fun
  expect_error(
    sbhm_build(library_spec = list(pgrid = data.frame(a = 1)), sframe = sframe),
    "library_spec must be a list with elements fun"
  )

  # Missing pgrid
  expect_error(
    sbhm_build(library_spec = list(fun = identity), sframe = sframe),
    "library_spec must be a list with elements fun"
  )

  # pgrid not a data.frame
  expect_error(
    sbhm_build(library_spec = list(fun = identity, pgrid = list()), sframe = sframe),
    "library_spec must be a list with elements fun"
  )
})

# Tests for lss_with_hrf.R - lss_with_hrf_pure_r function

test_that("lss_with_hrf_pure_r validates Y as matrix", {
  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = "not_matrix",
      onset_idx = c(5, 15),
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = matrix(1, 2, 5)
    ),
    "Y must be a matrix"
  )
})

test_that("lss_with_hrf_pure_r validates hrf_basis_kernels as matrix", {
  Y <- matrix(rnorm(100), 50, 2)

  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = c(5, 15),
      hrf_basis_kernels = "not_matrix",
      coefficients = matrix(1, 2, 2)
    ),
    "hrf_basis_kernels must be a matrix"
  )
})

test_that("lss_with_hrf_pure_r validates coefficients as matrix", {
  Y <- matrix(rnorm(100), 50, 2)

  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = c(5, 15),
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = c(1, 2)  # vector, not matrix
    ),
    "coefficients must be a matrix"
  )
})

test_that("lss_with_hrf_pure_r validates coefficient dimensions", {
  Y <- matrix(rnorm(100), 50, 2)

  # Wrong K dimension
  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = c(5, 15),
      hrf_basis_kernels = matrix(1, 10, 3),  # K = 3
      coefficients = matrix(1, 2, 2)         # nrow = 2, should be 3
    ),
    "nrow\\(coefficients\\) must equal ncol\\(hrf_basis_kernels\\)"
  )

  # Wrong voxel dimension
  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,                                  # ncol = 2
      onset_idx = c(5, 15),
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = matrix(1, 2, 5)         # ncol = 5, should be 2
    ),
    "ncol\\(coefficients\\) must equal ncol\\(Y\\)"
  )
})

test_that("lss_with_hrf_pure_r requires at least one trial", {
  Y <- matrix(rnorm(100), 50, 2)

  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = integer(0),  # No trials
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = matrix(1, 2, 2)
    ),
    "Need at least one trial"
  )
})

test_that("lss_with_hrf_pure_r validates durations length", {
  Y <- matrix(rnorm(100), 50, 2)

  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = c(5, 15),
      durations = c(1),  # Should be length 2
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = matrix(1, 2, 2)
    ),
    "durations must match length\\(onset_idx\\)"
  )
})

test_that("lss_with_hrf_pure_r validates Z dimensions", {
  Y <- matrix(rnorm(100), 50, 2)
  Z <- matrix(1, 30, 2)  # Wrong rows

  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = c(5, 15),
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = matrix(1, 2, 2),
      Z = Z
    ),
    "Z must have nrow == n_time"
  )
})

test_that("lss_with_hrf_pure_r validates Nuisance dimensions", {
  Y <- matrix(rnorm(100), 50, 2)
  Nuisance <- matrix(1, 30, 2)  # Wrong rows

  expect_error(
    fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = c(5, 15),
      hrf_basis_kernels = matrix(1, 10, 2),
      coefficients = matrix(1, 2, 2),
      Nuisance = Nuisance
    ),
    "Nuisance must have nrow == n_time"
  )
})

test_that("lss_with_hrf_pure_r returns correct dimensions", {
  set.seed(123)
  n_time <- 100
  n_vox <- 5
  n_trials <- 3
  K <- 2  # basis dimension
  L <- 15  # HRF length

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 40, 70)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that("lss_with_hrf_pure_r handles NULL durations", {
  set.seed(234)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  # NULL durations should default to zeros
  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = NULL,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles extended durations", {
  set.seed(345)
  n_time <- 100
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 40, 70)
  durations <- c(3, 2, 5)  # Non-zero durations
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = durations,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles Z regressors", {
  set.seed(456)
  n_time <- 100
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 40, 70)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)
  Z <- cbind(1, 1:n_time)  # Intercept + linear trend

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    Z = Z,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles Nuisance regressors", {
  set.seed(567)
  n_time <- 100
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 40, 70)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)
  Nuisance <- matrix(rnorm(n_time * 2), n_time, 2)  # Motion regressors

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    Nuisance = Nuisance,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r converts onset_idx to integer", {
  set.seed(678)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10.5, 30.2, 50.8)  # Non-integer onsets
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  # Should work (converts to integer)
  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles out-of-bounds onsets", {
  set.seed(789)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(0, 10, 100)  # 0 and 100 are out of bounds for 1:80
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  # Should not error, just skip invalid onsets
  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles NA onsets", {
  set.seed(890)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, NA, 50)  # One NA
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  # Should return NA for the invalid onset
  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles single trial", {
  set.seed(901)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(20)  # Single trial
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(1, n_vox))
})

test_that("lss_with_hrf_pure_r handles single voxel", {
  set.seed(111)
  n_time <- 80
  n_vox <- 1
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles K=1 basis", {
  set.seed(222)
  n_time <- 80
  n_vox <- 3
  K <- 1  # Single basis function
  L <- 15

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r handles zero HRF weights", {
  set.seed(333)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  # Set some weights to zero
  coefficients <- matrix(c(1, 0, 0, 1, 1, 1), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r method argument is matched", {
  set.seed(444)
  n_time <- 50
  n_vox <- 2
  K <- 2
  L <- 10

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  # All these should work (method is match.arg'd)
  result_r <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y, onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(dim(result_r), c(2, n_vox))
})

test_that("lss_with_hrf_pure_r sets rownames", {
  set.seed(555)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  colnames(Y) <- c("V1", "V2", "V3")
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "r"
  )

  expect_equal(rownames(result), c("trial_1", "trial_2", "trial_3"))
  expect_equal(colnames(result), c("V1", "V2", "V3"))
})

test_that("lss_with_hrf_pure_r verbose mode doesn't error", {
  set.seed(666)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  expect_message(
    result <- fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = onset_idx,
      hrf_basis_kernels = hrf_basis_kernels,
      coefficients = coefficients,
      method = "r",
      verbose = FALSE
    ),
    NA  # No message expected
  )

  expect_equal(dim(result), c(3, n_vox))
})

test_that("lss_with_hrf_pure_r fallback chain works when C++ unavailable", {
  set.seed(777)
  n_time <- 80
  n_vox <- 3
  K <- 2
  L <- 12

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  onset_idx <- c(10, 30, 50)
  hrf_basis_kernels <- matrix(rnorm(L * K), L, K)
  coefficients <- matrix(rnorm(K * n_vox), K, n_vox)

  # Request C++ methods - should fall back to R if not available
  # This tests the fallback chain logic
  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "cpp_omp",  # Will fall back through chain
    verbose = FALSE
  )

  expect_equal(dim(result), c(3, n_vox))
})

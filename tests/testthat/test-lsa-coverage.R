# Tests for lsa.R (Least Squares All)

test_that("lsa returns correct dimensions with R method", {
  set.seed(123)
  n_time <- 50
  n_trials <- 5
  n_vox <- 10

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  result <- lsa(Y, X, method = "r")

  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that("lsa returns correct dimensions with cpp method", {
  set.seed(234)
  n_time <- 50
  n_trials <- 5
  n_vox <- 10

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # cpp method falls back to R currently
  result <- lsa(Y, X, method = "cpp")

  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that("lsa works with Z nuisance regressors", {
  set.seed(345)
  n_time <- 60
  n_trials <- 4
  n_vox <- 8

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  Z <- cbind(1, 1:n_time)  # Intercept and linear trend

  result <- lsa(Y, X, Z = Z, method = "r")

  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that("lsa accepts Nuisance as alias for Z", {
  set.seed(456)
  n_time <- 50
  n_trials <- 4
  n_vox <- 8

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  nuis <- cbind(1, 1:n_time)

  result_Z <- lsa(Y, X, Z = nuis, method = "r")
  result_N <- lsa(Y, X, Nuisance = nuis, method = "r")

  expect_equal(result_Z, result_N)
})

test_that("lsa Z takes precedence over Nuisance", {
  set.seed(567)
  n_time <- 50
  n_trials <- 4
  n_vox <- 8

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  Z1 <- cbind(1, 1:n_time)
  Z2 <- cbind(1, (1:n_time)^2)

  result_Z <- lsa(Y, X, Z = Z1, method = "r")
  result_both <- lsa(Y, X, Z = Z1, Nuisance = Z2, method = "r")

  # Z takes precedence, so results should match
  expect_equal(result_Z, result_both)
})

test_that("lsa coerces non-matrix inputs", {
  set.seed(678)
  n_time <- 40
  n_trials <- 3
  n_vox <- 5

  # Create data.frames instead of matrices
  X <- data.frame(matrix(rnorm(n_time * n_trials), n_time, n_trials))
  Y <- data.frame(matrix(rnorm(n_time * n_vox), n_time, n_vox))

  result <- lsa(Y, X, method = "r")

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that("lsa errors on dimension mismatch between Y and X", {
  X <- matrix(rnorm(50 * 5), 50, 5)
  Y <- matrix(rnorm(40 * 10), 40, 10)  # Different rows

  expect_error(lsa(Y, X), "same number of rows")
})

test_that("lsa errors on dimension mismatch between Y and Z", {
  X <- matrix(rnorm(50 * 5), 50, 5)
  Y <- matrix(rnorm(50 * 10), 50, 10)
  Z <- matrix(rnorm(40 * 2), 40, 2)  # Different rows

  expect_error(lsa(Y, X, Z = Z), "same number of rows")
})

test_that("lsa handles missing values with warning", {
  set.seed(789)
  n_time <- 50
  n_trials <- 4
  n_vox <- 6

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Add some NAs
  Y[5, 2] <- NA
  Y[10, 3] <- NA

  expect_warning(result <- lsa(Y, X, method = "r"), "Removed.*rows with missing values")

  # Should still return valid result
  expect_true(is.matrix(result))
})

test_that("lsa preserves dimension names", {
  set.seed(890)
  n_time <- 40
  n_trials <- 3
  n_vox <- 4

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  colnames(X) <- c("trial_A", "trial_B", "trial_C")
  colnames(Y) <- c("V1", "V2", "V3", "V4")

  result <- lsa(Y, X, method = "r")

  expect_equal(rownames(result), c("trial_A", "trial_B", "trial_C"))
  expect_equal(colnames(result), c("V1", "V2", "V3", "V4"))
})

test_that("lsa method argument is validated", {
  X <- matrix(rnorm(50 * 5), 50, 5)
  Y <- matrix(rnorm(50 * 10), 50, 10)

  expect_error(lsa(Y, X, method = "invalid"), "should be one of")
})

test_that("lsa handles single voxel (Y is n_time x 1)", {
  set.seed(901)
  n_time <- 50
  n_trials <- 5

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * 1), n_time, 1)  # Single voxel

  result <- lsa(Y, X, method = "r")

  expect_equal(dim(result), c(n_trials, 1))
})

test_that("lsa R and cpp methods give same results", {
  set.seed(111)
  n_time <- 50
  n_trials <- 5
  n_vox <- 10

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  result_r <- lsa(Y, X, method = "r")
  result_cpp <- lsa(Y, X, method = "cpp")

  # cpp currently falls back to R, so they should be identical
  expect_equal(result_r, result_cpp)
})

test_that("lsa with Z recovers correct betas", {
  set.seed(222)
  n_time <- 100
  n_trials <- 3
  n_vox <- 5

  # Create known design
  X <- matrix(0, n_time, n_trials)
  X[10:20, 1] <- 1
  X[40:50, 2] <- 1
  X[70:80, 3] <- 1

  # Known betas
  true_betas <- matrix(c(2, -1, 3, 1, 0, -2, 0.5, 1.5, -0.5, 2, -1, 0.5, 1, 0, -1), n_trials, n_vox)

  # Create Y with signal + noise
  Y <- X %*% true_betas + matrix(rnorm(n_time * n_vox, sd = 0.1), n_time, n_vox)

  result <- lsa(Y, X, method = "r")

  # Betas should be close to true values (ignore dimnames)
  expect_equal(unname(result), unname(true_betas), tolerance = 0.5)
})

test_that(".lsa_r handles single response correctly", {
  set.seed(333)
  n_time <- 50
  n_trials <- 4

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time), n_time, 1)  # Single column
  Z <- matrix(1, n_time, 1)  # Just intercept

  result <- fmrilss:::.lsa_r(Y, X, Z)

  expect_equal(dim(result), c(n_trials, 1))
})

test_that(".lsa_cpp calls .lsa_r (current fallback)", {
  set.seed(444)
  n_time <- 40
  n_trials <- 3
  n_vox <- 5

  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  result_r <- fmrilss:::.lsa_r(Y, X, NULL)
  result_cpp <- fmrilss:::.lsa_cpp(Y, X, NULL)

  expect_equal(result_r, result_cpp)
})

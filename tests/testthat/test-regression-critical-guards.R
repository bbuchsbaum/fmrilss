library(testthat)
library(fmrilss)

test_that("lss_beta_cpp does not produce NaNs when bt2 is zero", {
  set.seed(1)

  # Construct C so that for trial 1, the "other trials" regressor is identically 0.
  # This yields bt2 == 0 for that trial and previously produced 0/0 -> NaNs.
  C <- cbind(c(1, 0, 0, 0), rep(0, 4))
  Y <- matrix(rnorm(4 * 3), 4, 3)

  out <- lss_beta_cpp(C, Y)
  expect_false(any(!is.finite(out)))
})

test_that("mixed_solve errors cleanly when X saturates the data", {
  set.seed(2)
  y <- rnorm(3)
  X <- diag(3) # n_filtered == p -> no residual degrees of freedom

  expect_error(
    mixed_solve(Y = y, X = X),
    "Need more non-NA observations than columns in X"
  )
})

test_that("oasisk_betas errors cleanly when Gram matrix is not SPD", {
  K <- 1
  N <- 1
  V <- 2

  D <- array(-1, dim = c(K, K, N))
  C <- array(0, dim = c(K, K, N))
  E <- array(-1, dim = c(K, K, N))
  N1 <- matrix(0, nrow = K * N, ncol = V)
  SY <- matrix(0, nrow = K, ncol = V)

  expect_error(
    oasisk_betas(D, C, E, N1, SY, diag_eps = 1e-10),
    "Cholesky failed"
  )
})


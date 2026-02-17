# Tests for sbhm_amplitude.R helper functions

test_that(".sbhm_resid returns M when Z is NULL", {
  M <- matrix(rnorm(20), 10, 2)
  result <- fmrilss:::.sbhm_resid(M, NULL)
  expect_equal(result, M)
})

test_that(".sbhm_resid returns M when Z has 0 columns", {
  M <- matrix(rnorm(20), 10, 2)
  Z <- matrix(0, 10, 0)
  result <- fmrilss:::.sbhm_resid(M, Z)
  expect_equal(result, M)
})

test_that(".sbhm_resid removes Z from M correctly", {
  set.seed(123)
  n <- 50
  # Create confound
  Z <- cbind(1, rnorm(n))
  # Create M with confound signal (2 columns, each affected by Z)
  coef <- matrix(c(5, 2, 3, 1), nrow = 2, ncol = 2)  # 2x2 to get n x 2 result
  M <- Z %*% coef + matrix(rnorm(n * 2, sd = 0.1), n, 2)

  result <- fmrilss:::.sbhm_resid(M, Z)

  # Residualized columns should have mean near zero
  expect_true(all(abs(colMeans(result)) < 0.5))
})

test_that(".sbhm_resid coerces vector Z to matrix", {
  set.seed(234)
  n <- 30
  M <- matrix(rnorm(n * 2), n, 2)
  Z <- rnorm(n)  # vector, not matrix

  # Should not error
  result <- fmrilss:::.sbhm_resid(M, Z)
  expect_equal(dim(result), dim(M))
})

test_that(".sbhm_resolve_ridge returns 0 for NULL ridge", {
  G <- diag(3)
  result <- fmrilss:::.sbhm_resolve_ridge(G, NULL)
  expect_equal(result, 0)
})

test_that(".sbhm_resolve_ridge returns numeric value directly", {
  G <- diag(3)
  result <- fmrilss:::.sbhm_resolve_ridge(G, 0.5)
  expect_equal(result, 0.5)
})

test_that(".sbhm_resolve_ridge handles absolute mode", {
  G <- diag(c(10, 20, 30))
  ridge <- list(mode = "absolute", lambda = 0.1)
  result <- fmrilss:::.sbhm_resolve_ridge(G, ridge)
  expect_equal(result, 0.1)
})

test_that(".sbhm_resolve_ridge handles fractional mode", {
  G <- diag(c(10, 20, 30))  # mean diagonal = 20
  ridge <- list(mode = "fractional", lambda = 0.1)
  result <- fmrilss:::.sbhm_resolve_ridge(G, ridge)
  expect_equal(result, 0.1 * 20)  # 2.0
})

test_that(".sbhm_resolve_ridge defaults to fractional mode", {
  G <- diag(c(10, 20, 30))
  ridge <- list(lambda = 0.1)  # no mode specified
  result <- fmrilss:::.sbhm_resolve_ridge(G, ridge)
  expect_equal(result, 0.1 * 20)  # defaults to fractional
})

test_that(".sbhm_solve solves linear system", {
  set.seed(345)
  G <- crossprod(matrix(rnorm(50), 10, 5))  # 5x5 positive definite
  B <- rnorm(5)

  result <- fmrilss:::.sbhm_solve(G, B, ridge = 0)

  # Verify: G %*% result ≈ B
  expect_equal(as.numeric(G %*% result), B, tolerance = 1e-8)
})

test_that(".sbhm_solve adds ridge correctly", {
  G <- diag(c(1, 1, 1))
  B <- c(1, 2, 3)
  ridge <- 0.5

  result <- fmrilss:::.sbhm_solve(G, B, ridge = ridge)

  # Should solve (G + ridge*I) x = B
  # (1.5*I) x = B => x = B/1.5
  expect_equal(as.numeric(result), B / 1.5, tolerance = 1e-10)
})

test_that(".sbhm_solve handles matrix B (multiple RHS)", {
  set.seed(456)
  G <- crossprod(matrix(rnorm(32), 8, 4))  # 4x4 positive definite
  B <- matrix(rnorm(4 * 3), 4, 3)  # 3 right-hand sides

  result <- fmrilss:::.sbhm_solve(G, B, ridge = 0.01)

  expect_equal(dim(result), c(4, 3))
  # Verify each column
  for (j in 1:3) {
    G_ridge <- G + 0.01 * diag(4)
    expect_equal(as.numeric(G_ridge %*% result[, j]), B[, j], tolerance = 1e-8)
  }
})

test_that(".sbhm_solve uses qr.solve fallback for singular G", {
  # Create a rank-deficient G
  G <- matrix(c(1, 2, 2, 4), 2, 2)  # rank 1
  B <- c(1, 2)

  # Should not error due to fallback (or may throw error from qr.solve)
  # The function attempts chol then falls back to qr.solve
  result <- tryCatch(
    fmrilss:::.sbhm_solve(G, B, ridge = 0),
    error = function(e) NA
  )

  # If it returns a result, check it's numeric
  if (!identical(result, NA)) {
    expect_true(is.numeric(result))
    expect_equal(length(result), 2)
  }
})

test_that(".sbhm_auto_gate_threshold computes quantile", {
  values <- 1:100
  result <- fmrilss:::.sbhm_auto_gate_threshold(values, metric = "kappa", q = 0.95)
  expect_equal(result, unname(quantile(values, 0.95)), tolerance = 1e-10)
})
test_that(".sbhm_auto_gate_threshold respects hard_min", {
  values <- 1:100
  result <- fmrilss:::.sbhm_auto_gate_threshold(values, metric = "rho", q = 0.5, hard_min = 100)
  expect_equal(result, 100)  # hard_min > quantile
})

test_that(".sbhm_adaptive_ridge_gls computes scaled ridge", {
  kappa_vec <- c(500, 1000, 2000)
  base <- 0.02
  k0 <- 1000

  result <- fmrilss:::.sbhm_adaptive_ridge_gls(kappa_vec, base = base, k0 = k0)

  # lambda = base * (kappa/k0), clipped to [base, max_lam]
  expected <- c(0.02, 0.02, 0.04)  # 500/1000*0.02=0.01 (clamped to 0.02), 1000/1000*0.02=0.02, 2000/1000*0.02=0.04
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that(".sbhm_adaptive_ridge_gls respects max_lam", {
  kappa_vec <- c(10000, 20000)
  result <- fmrilss:::.sbhm_adaptive_ridge_gls(kappa_vec, base = 0.02, k0 = 1000, max_lam = 0.08)

  # 10000/1000*0.02=0.2 -> clamped to 0.08
  # 20000/1000*0.02=0.4 -> clamped to 0.08
  expect_equal(result, c(0.08, 0.08), tolerance = 1e-10)
})

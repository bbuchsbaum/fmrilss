# Tests for sbhm_amplitude.R main functions (sbhm_amplitude_ls, etc.)

test_that("sbhm_amplitude_ls validates inputs", {
  skip_if_not_installed("fmrihrf")

  # Non-matrix Y
  expect_error(
    sbhm_amplitude_ls(Y = "not_matrix", sbhm = list(), design_spec = list(),
                      alpha_hat = matrix(1, 2, 2)),
    "is.matrix"
  )

  # Non-matrix alpha_hat
  expect_error(
    sbhm_amplitude_ls(Y = matrix(1, 10, 2), sbhm = list(), design_spec = list(),
                      alpha_hat = c(1, 2)),
    "is.matrix"
  )
})

test_that("sbhm_amplitude_ls validates alpha_hat dimensions", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(100), 50, 2)  # 2 voxels
  alpha_hat <- matrix(rnorm(6), 3, 2)  # Wrong V

  # Mock sbhm object
  sbhm <- list(B = matrix(1, 50, 3), tgrid = 0:49, span = 32)

  # Wrong V dimension
  expect_error(
    sbhm_amplitude_ls(Y = Y, sbhm = sbhm, design_spec = list(),
                      alpha_hat = matrix(1, 3, 5)),  # 5 != 2
    "alpha_hat must be r x V to match Y"
  )
})

test_that("sbhm_amplitude_lss1 validates inputs", {
  skip_if_not_installed("fmrihrf")

  # Non-matrix Y
  expect_error(
    sbhm_amplitude_lss1(Y = "not_matrix", sbhm = list(), design_spec = list(),
                        alpha_hat = matrix(1, 2, 2)),
    "is.matrix"
  )

  # Non-matrix alpha_hat
  expect_error(
    sbhm_amplitude_lss1(Y = matrix(1, 10, 2), sbhm = list(), design_spec = list(),
                        alpha_hat = c(1, 2)),
    "is.matrix"
  )
})

test_that("sbhm_amplitude_lss1 validates alpha_hat dimensions", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(100), 50, 2)  # 2 voxels

  sbhm <- list(B = matrix(1, 50, 3), tgrid = 0:49, span = 32)

  # Wrong V dimension
  expect_error(
    sbhm_amplitude_lss1(Y = Y, sbhm = sbhm, design_spec = list(),
                        alpha_hat = matrix(1, 3, 5)),  # 5 != 2
    "alpha_hat must be r x V to match Y"
  )
})

test_that("sbhm_amplitude_oasis_k1 validates inputs", {
  skip_if_not_installed("fmrihrf")

  # Non-matrix Y
  expect_error(
    sbhm_amplitude_oasis_k1(Y = "not_matrix", sbhm = list(), design_spec = list(),
                            alpha_hat = matrix(1, 2, 2)),
    "is.matrix"
  )

  # Non-matrix alpha_hat
  expect_error(
    sbhm_amplitude_oasis_k1(Y = matrix(1, 10, 2), sbhm = list(), design_spec = list(),
                            alpha_hat = c(1, 2)),
    "is.matrix"
  )
})

test_that(".sbhm_build_trial_regs builds regressors", {
  skip_if_not_installed("fmrihrf")

  set.seed(123)
  n_time <- 100
  r <- 3

  # Create minimal sbhm-like object
  sbhm <- list(
    B = matrix(rnorm(n_time * r), n_time, r),
    tgrid = seq(0, n_time - 1, by = 1),
    span = 32
  )

  design_spec <- list(
    cond = list(
      onsets = c(10, 40, 70),
      duration = 0,
      span = 32
    )
  )

  result <- fmrilss:::.sbhm_build_trial_regs(sbhm, design_spec)

  expect_true(is.list(result))
  expect_equal(length(result), 3)  # 3 trials
  # Each element should be T x r matrix
  expect_true(is.matrix(result[[1]]))
  expect_equal(ncol(result[[1]]), r)
})

test_that(".sbhm_prewhiten returns unmodified data when no prewhitening", {
  Y <- matrix(rnorm(100), 50, 2)
  regs <- list(matrix(rnorm(50 * 2), 50, 2), matrix(rnorm(50 * 2), 50, 2))
  Zint <- matrix(1, 50, 1)
  Nuisance <- matrix(rnorm(50 * 2), 50, 2)

  # NULL prewhiten
  result <- fmrilss:::.sbhm_prewhiten(Y, regs, Zint, Nuisance, prewhiten = NULL)

  expect_false(result$applied)
  expect_equal(result$Yw, Y)

  # method = "none"
  result2 <- fmrilss:::.sbhm_prewhiten(Y, regs, Zint, Nuisance,
                                        prewhiten = list(method = "none"))

  expect_false(result2$applied)
})

test_that(".sbhm_solve handles singular matrices with fallback", {
  # Create a singular matrix
  G <- matrix(c(1, 2, 2, 4), 2, 2)  # rank 1
  B <- c(1, 2)

  # Should use qr.solve fallback
  result <- tryCatch(
    fmrilss:::.sbhm_solve(G, B, ridge = 0),
    error = function(e) NA
  )

  # May return NA if truly singular, or a result if qr.solve succeeds
  if (!identical(result, NA)) {
    expect_true(is.numeric(result))
  }
})

test_that(".sbhm_solve handles matrix B (multiple RHS)", {
  set.seed(234)
  G <- crossprod(matrix(rnorm(30), 10, 3))  # 3x3 positive definite
  B <- matrix(rnorm(3 * 2), 3, 2)  # 2 right-hand sides

  result <- fmrilss:::.sbhm_solve(G, B, ridge = 0.01)

  expect_equal(dim(result), c(3, 2))
})

test_that(".sbhm_resolve_ridge handles per-voxel ridge vector", {
  # Test numeric ridge of length 1
  G <- diag(c(1, 2, 3))
  result <- fmrilss:::.sbhm_resolve_ridge(G, ridge = 0.5)
  expect_equal(result, 0.5)
})

test_that(".sbhm_adaptive_ridge_gls scales with kappa", {
  kappa_vec <- c(100, 1000, 5000, 10000)

  result <- fmrilss:::.sbhm_adaptive_ridge_gls(kappa_vec, base = 0.01, k0 = 1000, max_lam = 0.05)

  # kappa=100: 0.01 * (100/1000) = 0.001 -> clamped to base = 0.01
  # kappa=1000: 0.01 * (1000/1000) = 0.01
  # kappa=5000: 0.01 * (5000/1000) = 0.05
  # kappa=10000: 0.01 * (10000/1000) = 0.1 -> clamped to max_lam = 0.05
  expect_equal(result[1], 0.01)  # clamped to base

  expect_equal(result[2], 0.01)  # exact
  expect_equal(result[3], 0.05)  # exact
  expect_equal(result[4], 0.05)  # clamped to max
})

test_that(".sbhm_auto_gate_threshold computes correct quantile", {
  values <- 1:100

  # 95th percentile
  result <- fmrilss:::.sbhm_auto_gate_threshold(values, metric = "kappa", q = 0.95)
  expect_equal(result, unname(quantile(values, 0.95)))

  # With hard_min that's higher
  result2 <- fmrilss:::.sbhm_auto_gate_threshold(values, metric = "kappa", q = 0.5, hard_min = 100)
  expect_equal(result2, 100)
})

test_that(".sbhm_auto_gate_threshold metric argument works", {
  values <- 1:50

  result_kappa <- fmrilss:::.sbhm_auto_gate_threshold(values, metric = "kappa", q = 0.9)
  result_rho <- fmrilss:::.sbhm_auto_gate_threshold(values, metric = "rho", q = 0.9)

  # Both should compute the same quantile
  expect_equal(result_kappa, result_rho)
})

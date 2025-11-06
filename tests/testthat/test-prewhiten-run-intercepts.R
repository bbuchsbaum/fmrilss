test_that("prewhiten is robust with per-run intercept dummies in Nuisance", {
  set.seed(7)
  n_per_run <- 30L; R <- 3L
  n <- n_per_run * R
  V <- 6L

  # Data
  Y <- matrix(rnorm(n * V), n, V)

  # Random trial design (not singular)
  X <- matrix(rnorm(n * 8), n, 8)

  # Build per-run intercept dummies (one constant column per run)
  runs <- rep(seq_len(R), each = n_per_run)
  N_run <- stats::model.matrix(~ 0 + factor(runs))

  # Case 1: only run dummies (no global intercept)
  beta1_r <- lss(Y, X, Z = NULL,
                 Nuisance = N_run,
                 method = "r_optimized",
                 prewhiten = list(method = "ar", p = 1))
  expect_equal(dim(beta1_r), c(ncol(X), V))

  beta1_o <- lss(Y, X, Z = NULL,
                 Nuisance = N_run,
                 method = "oasis",
                 prewhiten = list(method = "ar", p = 1))
  expect_equal(dim(beta1_o), c(ncol(X), V))

  # Case 2: run dummies + explicit global intercept (perfect collinearity)
  N_collinear <- cbind(1, N_run)
  beta2_r <- lss(Y, X, Z = NULL,
                 Nuisance = N_collinear,
                 method = "r_optimized",
                 prewhiten = list(method = "ar", p = 1))
  expect_equal(dim(beta2_r), c(ncol(X), V))

  beta2_o <- lss(Y, X, Z = NULL,
                 Nuisance = N_collinear,
                 method = "oasis",
                 prewhiten = list(method = "ar", p = 1))
  expect_equal(dim(beta2_o), c(ncol(X), V))
})


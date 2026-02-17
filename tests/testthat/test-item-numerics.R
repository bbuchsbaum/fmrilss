test_that("item_compute_u falls back under collinearity and returns finite output", {
  set.seed(2001)

  X_base <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  X_t <- cbind(X_base, X_base[, 1, drop = FALSE])

  expect_warning(
    U <- item_compute_u(X_t, ridge = 0, method = "chol"),
    "requested 'chol' but used"
  )

  expect_true(all(is.finite(U)))

  diag_u <- attr(U, "item_diagnostics")
  expect_true(diag_u$solver_path %in% c("svd", "pinv"))
  expect_true(is.finite(diag_u$condition_number) || is.infinite(diag_u$condition_number))
})

test_that("item_fit falls back on singular GLS system and returns finite W_hat", {
  set.seed(2002)

  n_trials <- 24
  Gamma <- matrix(rnorm(n_trials * 4), nrow = n_trials, ncol = 4)
  Gamma <- cbind(Gamma, Gamma[, 1, drop = FALSE])

  T_train <- matrix(rnorm(n_trials * 2), nrow = n_trials, ncol = 2)
  U_train <- diag(n_trials)

  expect_warning(
    W_hat <- item_fit(Gamma, T_train, U_train, ridge = 0, method = "chol"),
    "requested 'chol' but used"
  )

  expect_equal(dim(W_hat), c(ncol(Gamma), ncol(T_train)))
  expect_true(all(is.finite(W_hat)))

  diag_w <- attr(W_hat, "item_diagnostics")
  expect_true(diag_w$solver_path %in% c("svd", "pinv"))
})

test_that("item_compute_u accepts sparse Matrix input when available", {
  skip_if_not_installed("Matrix")
  set.seed(2003)

  X_t <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
  V_prec <- Matrix::Diagonal(n = 30, x = 2)

  U <- item_compute_u(X_t, V = V_prec, v_type = "precision", ridge = 1e-6)

  expect_equal(dim(U), c(5, 5))
  expect_equal(U, t(U), tolerance = 1e-8)
  expect_true(all(is.finite(U)))
})

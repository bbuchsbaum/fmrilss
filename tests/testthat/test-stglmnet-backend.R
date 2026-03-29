test_that("lss(method='stglmnet') returns trial-by-voxel beta matrix", {
  skip_if_not_installed("glmnet")

  set.seed(4101)
  n <- 90
  p <- 12
  v <- 3
  X <- matrix(rnorm(n * p), n, p)
  Z <- cbind(1, scale(seq_len(n)))
  B <- matrix(rnorm(p * v), p, v)
  Y <- X %*% B + matrix(rnorm(n * v, sd = 0.35), n, v)

  beta_hat <- lss(
    Y = Y,
    X = X,
    Z = Z,
    method = "stglmnet",
    stglmnet = stglmnet_options(
      mode = "cv",
      lambda = exp(seq(log(0.3), log(0.03), length.out = 10)),
      cv_folds = 3,
      cv_type.measure = "mse"
    )
  )

  expect_true(is.matrix(beta_hat))
  expect_equal(dim(beta_hat), c(p, v))
  expect_true(all(is.finite(beta_hat)))
})

test_that("stglmnet return_fit exposes selected lambda and beta matrix", {
  skip_if_not_installed("glmnet")

  set.seed(4102)
  n <- 80
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(X %*% rnorm(p) + rnorm(n, sd = 0.4), ncol = 1L)

  fit <- lss(
    Y = Y,
    X = X,
    method = "stglmnet",
    stglmnet = stglmnet_options(
      mode = "cv",
      lambda = exp(seq(log(0.25), log(0.025), length.out = 8)),
      cv_folds = 4,
      cv_type.measure = "correlation",
      return_fit = TRUE
    )
  )

  expect_true(is.list(fit))
  expect_true(is.matrix(fit$beta))
  expect_true(is.finite(fit$lambda))
  expect_true(inherits(fit$fit, "fmrilss_stglmnet_fit"))
  expect_true(inherits(fit$cv, "fmrilss_cv_stglmnet"))
})

test_that("pool_to_mean preserves ridge solutions under isotropic penalties", {
  skip_if_not_installed("glmnet")

  set.seed(4103)
  n <- 96
  p <- 12
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(X %*% rnorm(p) + rnorm(n, sd = 0.5), ncol = 1L)
  lam <- 0.08

  beta_std <- lss(
    Y = Y,
    X = X,
    method = "stglmnet",
    stglmnet = stglmnet_options(
      mode = "fixed",
      alpha = 0,
      lambda = lam,
      overlap_strategy = "none"
    )
  )

  beta_pool <- lss(
    Y = Y,
    X = X,
    method = "stglmnet",
    stglmnet = stglmnet_options(
      mode = "fixed",
      alpha = 0,
      lambda = lam,
      overlap_strategy = "none",
      pool_to_mean = TRUE,
      pool_strength = 1,
      pool_mean_penalty = 1,
      pool_scale_by_overlap = FALSE
    )
  )

  expect_equal(as.numeric(beta_pool), as.numeric(beta_std), tolerance = 3e-4)
})

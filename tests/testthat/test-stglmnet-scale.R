test_that("automatic stglmnet paths preserve response units and CV scores", {
  set.seed(1717)
  n <- 96L
  p <- 8L
  X <- matrix(rnorm(n * p), n, p)
  Z <- cbind(1, seq(-1, 1, length.out = n))
  Y <- X %*% matrix(rnorm(p * 3), p, 3) +
    Z %*% rbind(c(100, 400, 900), c(1, 3, -2)) +
    matrix(rnorm(n * 3, sd = 0.2), n, 3)
  folds <- rep(1:4, length.out = n)
  for (nv in c(1L, 3L)) {
    opts <- list(cv_foldid = folds, cv_type.measure = "mse", return_fit = TRUE)
    a <- lss(Y[, seq_len(nv), drop = FALSE], X, Z, method = "stglmnet", stglmnet = opts)
    b <- lss(5000 * Y[, seq_len(nv), drop = FALSE], X, Z, method = "stglmnet", stglmnet = opts)
    expect_equal(b$beta / 5000, a$beta, tolerance = 1e-8)
    expect_equal(b$fit$fit$lambda, a$fit$fit$lambda, tolerance = 1e-10)
    expect_equal(b$lambda, a$lambda, tolerance = 1e-10)
    expect_equal(b$cv$cvm / 5000^2, a$cv$cvm, tolerance = 1e-8)
    expect_equal(b$fit$response_scale / 5000, a$fit$response_scale, tolerance = 1e-12)
    expect_identical(a$fit$lambda_scale, "normalized_response")
    expect_equal(
      predict(b$fit$fit, newx = b$fit$x_fit, s = b$lambda) / 5000,
      predict(a$fit$fit, newx = a$fit$x_fit, s = a$lambda), tolerance = 1e-8
    )
    expect_equal(b$fit$fit$nulldev / 5000^2, a$fit$fit$nulldev, tolerance = 1e-8)
  }
})

test_that("explicit stglmnet penalties retain raw glmnet semantics", {
  set.seed(1718)
  X <- matrix(rnorm(120 * 6), 120, 6)
  Y <- 100 * matrix(rnorm(120 * 3), 120, 3)
  for (alpha in c(0, 0.2, 1)) {
    for (lambda in list(0.2, c(0.3, 0.2, 0.1))) {
      actual <- fmrilss:::.stg_fit_core(Y, X, alpha = alpha, lambda = lambda)
      reference <- glmnet::glmnet(X, Y, family = "mgaussian", alpha = alpha,
                                  lambda = lambda, standardize = FALSE, intercept = FALSE)
      expect_equal(actual$fit$beta, reference$beta, tolerance = 1e-12)
      expect_identical(actual$response_scale, 1)
      expect_identical(actual$lambda_scale, "original_response")
    }
  }
})

test_that("multivariate ridge agrees with a direct penalized least-squares solve", {
  set.seed(1719)
  X <- matrix(rnorm(150 * 5), 150, 5)
  Z <- matrix(rnorm(150), 150, 1)
  Y <- 50 * matrix(rnorm(150 * 3), 150, 3)
  fit <- fmrilss:::.stg_fit_core(Y, X, Z, alpha = 0, lambda = 0.2,
                                overlap_strategy = "additive", overlap_strength = 2)
  D <- cbind(X, Z)
  penalty <- fit$penalty.factor * ncol(D) / sum(fit$penalty.factor)
  oracle <- solve(crossprod(D) / nrow(D) + diag(0.2 * penalty),
                  crossprod(D, Y) / nrow(D))
  actual <- do.call(cbind, lapply(fit$fit$beta, as.matrix))
  expect_equal(unname(actual), unname(oracle), tolerance = 2e-4)
})

test_that("response scaling covers pooled effects and zero-response boundaries", {
  set.seed(1720)
  X <- matrix(rnorm(80 * 6), 80, 6)
  Y <- cbind(X %*% rnorm(6) + rnorm(80), rep(0, 80))
  for (pooling in list(list(pool_to_mean = TRUE), list(graph_pool = TRUE))) {
    opts <- c(list(mode = "fixed", return_fit = TRUE), pooling)
    a <- lss(Y, X, method = "stglmnet", stglmnet = opts)
    b <- lss(5000 * Y, X, method = "stglmnet", stglmnet = opts)
    expect_true(all(is.finite(a$beta)))
    expect_equal(b$beta / 5000, a$beta, tolerance = 1e-8)
    expect_equal(unname(a$beta[, 2]), rep(0, ncol(X)), tolerance = 1e-12)
  }
  expect_error(lss(Y * 0, X, method = "stglmnet"),
               "non-zero finite response after nuisance projection")
})

test_that("structured run drifts do not collapse automatic paths at BOLD scales", {
  set.seed(1721)
  n <- 180L
  onsets <- seq(5, n - 15, length.out = 18)
  X <- vapply(onsets, function(onset) {
    t <- pmax(0, seq_len(n) - onset)
    h <- t^2 * exp(-t / 1.2)
    h / max(h)
  }, numeric(n))
  runs <- rep(1:3, each = 60)
  poly <- stats::poly(seq_len(60), degree = 6)
  Z <- cbind(1, do.call(cbind, lapply(1:3, function(r) {
    out <- matrix(0, n, 6)
    out[runs == r, ] <- poly
    out
  })))
  expect_equal(qr(Z)$rank, 19L)
  Y <- X %*% matrix(rnorm(18 * 4, 1, .5), 18, 4) +
    Z %*% rbind(rep(100, 4), matrix(rnorm(18 * 4), 18, 4)) +
    matrix(rnorm(n * 4, sd = .5), n, 4)
  opts <- list(cv_foldid = runs, cv_type.measure = "mse", return_fit = TRUE)
  reference <- lss(Y, X, Z, method = "stglmnet", stglmnet = opts)
  expect_true(all(is.finite(reference$beta)))
  for (scale in c(10, 100, 1000, 5000)) {
    result <- lss(Y * scale, X, Z, method = "stglmnet", stglmnet = opts)
    expect_equal(result$beta / scale, reference$beta, tolerance = 1e-7)
    expect_equal(length(result$fit$fit$lambda), length(reference$fit$fit$lambda))
  }
})


test_that("glmnet family objects retain their existing explicit-lambda behavior", {
  set.seed(1722)
  X <- matrix(rnorm(120 * 6), 120, 6)
  Y <- matrix(rnorm(120), 120, 1)
  family <- stats::gaussian()
  fit <- fmrilss:::.stg_fit_core(Y, X, family = family, lambda = .1)
  direct <- glmnet::glmnet(X, Y, family = family, alpha = .2, lambda = .1,
                           standardize = FALSE, intercept = FALSE)
  expect_equal(fit$fit$beta, direct$beta)
  expect_identical(fit$response_scale, 1)
  expect_identical(fit$lambda_scale, "original_response")
})

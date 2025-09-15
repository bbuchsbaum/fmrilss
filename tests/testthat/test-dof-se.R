skip_on_cran()

test_that("SEs respond to nuisance rank (QR residualization)", {
  set.seed(42)
  n  <- 120
  V  <- 20
  Tt <- 15

  # Trial design
  X <- matrix(rnorm(n*Tt), n, Tt)
  # True betas + noise
  Btrue <- matrix(rnorm(Tt*V, sd=0.3), Tt, V)
  E     <- matrix(rnorm(n*V, sd=0.7), n, V)
  Y     <- X %*% Btrue + E

  # Baseline confounds: intercept only
  Z0 <- matrix(1, n, 1)

  # Add 5 more nuisance columns (orthonormal), increasing rank by 5
  N_extra <- qr.Q(qr(matrix(rnorm(n*5), n, 5)))

  out0 <- lss(Y = Y, X = X, Z = Z0,
              method = "oasis",
              oasis  = list(return_se = TRUE, ridge_mode = "absolute",
                            ridge_x = 1e-8, ridge_b = 1e-8))

  out1 <- lss(Y = Y, X = X, Z = cbind(Z0, N_extra),
              method = "oasis",
              oasis  = list(return_se = TRUE, ridge_mode = "absolute",
                            ridge_x = 1e-8, ridge_b = 1e-8))

  # Same dimensions
  expect_equal(dim(out0$beta), c(Tt, V))
  expect_equal(dim(out1$beta), c(Tt, V))
  expect_equal(dim(out0$se),   c(Tt, V))
  expect_equal(dim(out1$se),   c(Tt, V))

  # With fewer residual degrees of freedom we expect larger SEs on average
  mean_se0 <- mean(out0$se, na.rm = TRUE)
  mean_se1 <- mean(out1$se, na.rm = TRUE)
  expect_gt(mean_se1, mean_se0)
})


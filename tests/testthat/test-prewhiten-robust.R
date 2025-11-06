test_that("prewhiten handles rank-deficient or wide designs without error", {
  set.seed(42)
  n <- 80; V <- 5
  Y <- matrix(rnorm(n * V), n, V)

  # Build X with near rank-deficiency (duplicate columns)
  base <- matrix(rnorm(n * 10), n, 10)
  dup  <- base[, 1:5]
  X <- cbind(base, dup)

  # Nuisance with run dummies + an explicit intercept to create collinearity
  runs <- rep(1:2, each = n/2)
  Nuis <- cbind(1, model.matrix(~ 0 + factor(runs)))

  # Should not error and should return the right shape
  beta <- lss(Y, X, Z = NULL,
              Nuisance = Nuis,
              method = "r_optimized",
              prewhiten = list(method = "ar", p = 1))
  expect_equal(dim(beta), c(ncol(X), V))
})


skip_on_cran()

test_that("OASIS adds an intercept by default and matches explicit intercept", {
  n  <- 80
  V  <- 5
  Tt <- 12

  set.seed(1)
  # Simple single-basis design: staggered impulses, lightly smoothed
  X <- matrix(0, n, Tt)
  for (i in seq_len(Tt)) {
    idx <- ( (i-1)*5 + 1 )
    if (idx <= n) X[idx, i] <- 1
  }
  X <- apply(X, 2, function(col) stats::filter(col, rep(1/3,3), sides=1))
  X[is.na(X)] <- 0

  beta_true <- matrix(rnorm(Tt*V, sd=0.4), Tt, V)
  Y <- X %*% beta_true + matrix(rnorm(n*V, mean = 0.2), n, V) # nonzero mean

  # 1) Default (adds intercept via OASIS)
  B1 <- lss(Y = Y, X = X, method = "oasis", oasis = list())

  # 2) Explicit intercept but disable default
  Z  <- matrix(1, n, 1)
  B2 <- lss(Y = Y, X = X, Z = Z, method = "oasis", oasis = list(add_intercept = FALSE))

  expect_equal(dim(B1), c(Tt, V))
  expect_equal(dim(B2), c(Tt, V))
  expect_lt(max(abs(B1 - B2)), 1e-7)
})

test_that("OASIS accepts Matrix inputs and preserves Y colnames", {
  skip_if_not_installed("Matrix")

  n <- 60; V <- 4; Tt <- 8
  set.seed(2)
  X <- matrix(rnorm(n*Tt), n, Tt)
  Y <- matrix(rnorm(n*V),  n, V)

  colnames(Y) <- paste0("vox", seq_len(V))

  # Convert to Matrix
  X_M <- Matrix::Matrix(X, sparse = FALSE)
  Z_M <- Matrix::Matrix(matrix(1, n, 1), sparse = FALSE)
  N_M <- Matrix::Matrix(matrix(rnorm(n*2), n, 2), sparse = FALSE)

  B <- lss(Y = Y, X = X_M, Z = Z_M, Nuisance = N_M, method = "oasis", oasis = list())
  expect_equal(dim(B), c(Tt, V))
  expect_identical(colnames(B), colnames(Y))
  expect_true(all(is.finite(B)))
})

test_that("OASIS AR(1) whitening path runs and returns finite betas", {
  n <- 80; V <- 3; Tt <- 9
  set.seed(4)
  X <- matrix(rnorm(n*Tt), n, Tt)

  # Build AR(1) noise
  ar <- 0.5
  E  <- matrix(0, n, V)
  E[1,] <- rnorm(V)
  for (t in 2:n) E[t,] <- ar*E[t-1,] + rnorm(V, sd=sqrt(1-ar^2))
  Y <- X %*% matrix(rnorm(Tt*V, sd=0.4), Tt, V) + E

  # Use new prewhiten API
  B <- lss(Y = Y, X = X, method = "oasis", prewhiten = list(method = "ar", p = 1))
  expect_equal(dim(B), c(Tt, V))
  expect_true(all(is.finite(B)))
})


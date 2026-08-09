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

test_that("mixed_solve handles matrix responses column by column", {
  set.seed(22)
  n <- 30
  X <- cbind(intercept = 1, trend = seq_len(n))
  Z <- matrix(rnorm(n * 3), nrow = n)
  Y <- cbind(voxel_a = rnorm(n), voxel_b = rnorm(n))

  multi <- mixed_solve(Y = Y, X = X, Z = Z)
  single <- lapply(seq_len(ncol(Y)), function(j) {
    mixed_solve(Y = Y[, j], X = X, Z = Z)
  })

  expect_equal(unname(multi$beta), t(vapply(single, `[[`, numeric(ncol(X)), "beta")))
  expect_equal(unname(multi$u), t(vapply(single, `[[`, numeric(ncol(Z)), "u")))
  expect_equal(unname(multi$Vu), vapply(single, `[[`, numeric(1), "Vu"))
  expect_equal(unname(multi$Ve), vapply(single, `[[`, numeric(1), "Ve"))
  expect_equal(rownames(multi$beta), colnames(Y))
})

test_that("all LSS backends tolerate rank-deficient confounds", {
  set.seed(23)
  n <- 40
  Y <- matrix(rnorm(n * 3), nrow = n)
  X <- matrix(0, nrow = n, ncol = 4)
  X[cbind(c(4, 12, 20, 28), seq_len(4))] <- 1
  Z <- cbind(intercept = rep(1, n), duplicate_intercept = rep(1, n))
  methods <- c("r_optimized", "cpp_optimized", "r_vectorized", "cpp", "naive")

  fits <- lapply(methods, function(method) {
    lss(Y, X, Z = Z, method = method)
  })

  expect_true(all(vapply(fits, function(x) all(is.finite(x)), logical(1))))
  expect_equal(fits[[2]], fits[[1]], tolerance = 1e-7)
  expect_equal(fits[[4]], fits[[1]], tolerance = 1e-7)
  expect_equal(fits[[5]], fits[[3]], tolerance = 1e-7)
})

test_that("block sizes are validated before blocked loops", {
  set.seed(24)
  Y <- matrix(rnorm(20), nrow = 10)
  X <- cbind(c(0, 1, rep(0, 8)), c(rep(0, 5), 1, rep(0, 4)))

  for (bad in list(0, -1, 1.5, NA_real_, Inf)) {
    expect_error(
      lss(Y, X, method = "cpp_optimized", block_size = bad),
      "block_size must be a positive integer"
    )
    expect_error(
      lss(Y, X, method = "oasis", oasis = list(block_cols = bad)),
      "block_cols must be a positive integer"
    )
  }

  expect_error(
    lss_fused_optim_cpp(matrix(1, 10, 1), Y, X, block_size = 0L),
    "block_size must be positive"
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

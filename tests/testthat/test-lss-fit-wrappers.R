test_that("lss_naive_fit produces correct results", {
  set.seed(123)
  n_timepoints <- 60
  n_trials <- 10
  n_voxels <- 20

  # Create trial design matrix
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  # Fixed effects
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])

  # Generate synthetic data
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test that wrapper produces same result as direct call
  beta_wrapper <- lss_naive_fit(Y, X, Z)
  beta_direct <- lss(Y, X, Z, method = "naive")

  expect_equal(beta_wrapper, beta_direct, tolerance = 1e-12)
})

test_that("lss_naive_fit handles NULL Z and Nuisance", {
  set.seed(234)
  n_timepoints <- 50
  n_trials <- 8
  n_voxels <- 15

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    X[onset:(onset + 2), i] <- 1
  }

  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Should work with NULL Z
  expect_no_error(beta <- lss_naive_fit(Y, X))
  expect_equal(dim(beta), c(n_trials, n_voxels))
})

test_that("lss_optimized_fit with cpp engine produces correct results", {
  set.seed(345)
  n_timepoints <- 60
  n_trials <- 10
  n_voxels <- 20

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test cpp engine
  beta_wrapper <- lss_optimized_fit(Y, X, Z, engine = "cpp")
  beta_direct <- lss(Y, X, Z, method = "cpp_optimized")

  expect_equal(beta_wrapper, beta_direct, tolerance = 1e-12)
})

test_that("lss_optimized_fit with r engine produces correct results", {
  set.seed(456)
  n_timepoints <- 60
  n_trials <- 10
  n_voxels <- 20

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test r engine
  beta_wrapper <- lss_optimized_fit(Y, X, Z, engine = "r")
  beta_direct <- lss(Y, X, Z, method = "r_optimized")

  expect_equal(beta_wrapper, beta_direct, tolerance = 1e-12)
})

test_that("lss_optimized_fit passes block_size correctly", {
  set.seed(567)
  n_timepoints <- 80
  n_trials <- 15
  n_voxels <- 30

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints))
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Different block sizes should produce same results (within numerical tolerance)
  beta_block32 <- lss_optimized_fit(Y, X, Z, engine = "cpp", block_size = 32)
  beta_block128 <- lss_optimized_fit(Y, X, Z, engine = "cpp", block_size = 128)

  expect_equal(beta_block32, beta_block128, tolerance = 1e-8)
})

test_that("lss_optimized_fit validates engine parameter", {
  set.seed(678)
  n_timepoints <- 40
  n_trials <- 5
  n_voxels <- 10

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 7 + 5
    X[onset:(onset + 2), i] <- 1
  }

  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  expect_error(lss_optimized_fit(Y, X, engine = "invalid"), "'arg' should be one of")
})

test_that("lss_optimized_fit handles Nuisance regressors", {
  set.seed(789)
  n_timepoints <- 70
  n_trials <- 12
  n_voxels <- 25

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints))
  Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # With and without nuisance should give different results
  beta_no_nuisance <- lss_optimized_fit(Y, X, Z, engine = "cpp")
  beta_with_nuisance <- lss_optimized_fit(Y, X, Z, Nuisance = Nuisance, engine = "cpp")

  expect_false(isTRUE(all.equal(beta_no_nuisance, beta_with_nuisance)))
})

test_that("lss_naive_fit and lss_optimized_fit produce equivalent results", {
  set.seed(890)
  n_timepoints <- 60
  n_trials <- 10
  n_voxels <- 20

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  beta_naive <- lss_naive_fit(Y, X, Z)
  beta_optimized_cpp <- lss_optimized_fit(Y, X, Z, engine = "cpp")
  beta_optimized_r <- lss_optimized_fit(Y, X, Z, engine = "r")

  # All implementations should be equivalent
  expect_equal(beta_naive, beta_optimized_cpp, tolerance = 1e-6)
  expect_equal(beta_naive, beta_optimized_r, tolerance = 1e-6)
  expect_equal(beta_optimized_cpp, beta_optimized_r, tolerance = 1e-6)
})

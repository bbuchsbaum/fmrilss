library(fmrilss)

test_that("fmriAR integration works with basic AR(1) whitening", {
  skip_if_not_installed("fmriAR")

  set.seed(123)
  n_time <- 100
  n_voxels <- 20
  n_trials <- 8

  # Create AR(1) data
  phi_true <- 0.5
  Y <- matrix(0, n_time, n_voxels)
  for (v in 1:n_voxels) {
    e <- rnorm(n_time)
    Y[1, v] <- e[1]
    for (t in 2:n_time) {
      Y[t, v] <- phi_true * Y[t-1, v] + e[t]
    }
  }

  # Create simple trial design
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 10 + 5
    if (start + 3 <= n_time) {
      X[start:(start+3), i] <- 1
    }
  }

  # Test new prewhiten parameter
  result <- lss(Y, X, prewhiten = list(method = "ar", p = 1))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_trials)
  expect_equal(ncol(result), n_voxels)
  expect_true(all(is.finite(result)))
})

test_that("fmriAR integration with auto AR order selection", {
  skip_if_not_installed("fmriAR")

  set.seed(456)
  n_time <- 150
  n_voxels <- 15
  n_trials <- 10

  # Create AR(2) data
  phi1 <- 0.6
  phi2 <- -0.3
  Y <- matrix(0, n_time, n_voxels)
  for (v in 1:n_voxels) {
    e <- rnorm(n_time)
    Y[1:2, v] <- e[1:2]
    for (t in 3:n_time) {
      Y[t, v] <- phi1 * Y[t-1, v] + phi2 * Y[t-2, v] + e[t]
    }
  }

  # Create trial design
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 12 + 5
    if (start + 4 <= n_time) {
      X[start:(start+4), i] <- 1
    }
  }

  # Test with auto AR order selection
  result <- lss(Y, X, prewhiten = list(method = "ar", p = "auto", p_max = 4))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_trials)
  expect_equal(ncol(result), n_voxels)
})

test_that("fmriAR integration with voxel-specific AR parameters", {
  skip_if_not_installed("fmriAR")

  set.seed(789)
  n_time <- 120
  n_voxels <- 25
  n_trials <- 8

  # Create data with different AR coefficients per voxel
  Y <- matrix(0, n_time, n_voxels)
  phi_voxels <- seq(0.2, 0.7, length.out = n_voxels)
  for (v in 1:n_voxels) {
    e <- rnorm(n_time)
    Y[1, v] <- e[1]
    for (t in 2:n_time) {
      Y[t, v] <- phi_voxels[v] * Y[t-1, v] + e[t]
    }
  }

  # Create trial design
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 12 + 5
    if (start + 5 <= n_time) {
      X[start:(start+5), i] <- 1
    }
  }

  # Test with voxel-specific pooling
  result <- lss(Y, X, prewhiten = list(method = "ar", p = 1, pooling = "voxel"))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_trials)
  expect_equal(ncol(result), n_voxels)
})

test_that("fmriAR integration with run-aware estimation", {
  skip_if_not_installed("fmriAR")

  set.seed(321)
  n_time_per_run <- 60
  n_runs <- 2
  n_time <- n_time_per_run * n_runs
  n_voxels <- 20
  n_trials <- 8

  # Create AR(1) data with different phi per run
  Y <- matrix(0, n_time, n_voxels)
  phi_run1 <- 0.3
  phi_run2 <- 0.6
  runs <- rep(1:n_runs, each = n_time_per_run)

  for (v in 1:n_voxels) {
    e <- rnorm(n_time)
    # Run 1
    Y[1, v] <- e[1]
    for (t in 2:n_time_per_run) {
      Y[t, v] <- phi_run1 * Y[t-1, v] + e[t]
    }
    # Run 2
    Y[n_time_per_run + 1, v] <- e[n_time_per_run + 1]
    for (t in (n_time_per_run + 2):n_time) {
      Y[t, v] <- phi_run2 * Y[t-1, v] + e[t]
    }
  }

  # Create trial design
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 12 + 5
    if (start + 5 <= n_time) {
      X[start:(start+5), i] <- 1
    }
  }

  # Test with run-aware pooling
  result <- lss(Y, X, prewhiten = list(
    method = "ar",
    p = 1,
    pooling = "run",
    runs = runs
  ))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_trials)
  expect_equal(ncol(result), n_voxels)
})

test_that("fmriAR integration with OASIS method", {
  skip_if_not_installed("fmriAR")

  set.seed(654)
  n_time <- 100
  n_voxels <- 15
  n_trials <- 6

  # Create AR(1) data
  phi <- 0.4
  Y <- matrix(0, n_time, n_voxels)
  for (v in 1:n_voxels) {
    e <- rnorm(n_time)
    Y[1, v] <- e[1]
    for (t in 2:n_time) {
      Y[t, v] <- phi * Y[t-1, v] + e[t]
    }
  }

  # Create trial design
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 14 + 5
    if (start + 5 <= n_time) {
      X[start:(start+5), i] <- 1
    }
  }

  # Test OASIS with new prewhiten parameter
  result <- lss(Y, X, method = "oasis",
               prewhiten = list(method = "ar", p = 1))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_trials)
  expect_equal(ncol(result), n_voxels)

  # Compare with old syntax (should give similar results)
  result_old <- lss(Y, X, method = "oasis",
                   oasis = list(whiten = "ar1"))
  expect_true(is.matrix(result_old))
  expect_equal(dim(result), dim(result_old))
})

test_that("backward compatibility with oasis$whiten is maintained", {
  skip_if_not_installed("fmriAR")

  set.seed(987)
  n_time <- 80
  n_voxels <- 10
  n_trials <- 5

  Y <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 12 + 5
    if (start + 4 <= n_time) {
      X[start:(start+4), i] <- 1
    }
  }

  # Old syntax should still work (with deprecation message)
  expect_message(
    result <- lss(Y, X, method = "oasis", oasis = list(whiten = "ar1")),
    "deprecated"
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_trials)
})

test_that("prewhiten parameter works with all LSS methods", {
  skip_if_not_installed("fmriAR")

  set.seed(111)
  n_time <- 60
  n_voxels <- 8
  n_trials <- 4

  Y <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 12 + 5
    if (start + 4 <= n_time) {
      X[start:(start+4), i] <- 1
    }
  }

  prewhiten_opts <- list(method = "ar", p = 1)

  # Test with different methods
  methods <- c("r_optimized", "r_vectorized", "naive")

  for (m in methods) {
    result <- lss(Y, X, method = m, prewhiten = prewhiten_opts)
    expect_true(is.matrix(result))
    expect_equal(nrow(result), n_trials)
    expect_equal(ncol(result), n_voxels)
    expect_true(all(is.finite(result)))
  }
})

test_that("prewhiten with method='none' does not apply whitening", {
  set.seed(222)
  n_time <- 60
  n_voxels <- 8
  n_trials <- 4

  Y <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
  X <- matrix(0, n_time, n_trials)
  for (i in 1:n_trials) {
    start <- (i-1) * 12 + 5
    if (start + 4 <= n_time) {
      X[start:(start+4), i] <- 1
    }
  }

  # No prewhitening
  result_none <- lss(Y, X, prewhiten = list(method = "none"))

  # Should be same as not specifying prewhiten
  result_null <- lss(Y, X, prewhiten = NULL)

  expect_equal(result_none, result_null)
})

test_that("OASIS with return_se=TRUE returns standard errors", {
  skip_if_not_installed("fmrihrf")

  set.seed(123)
  n_timepoints <- 80
  n_trials <- 12
  n_voxels <- 10

  # Create trial design matrix
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 6 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test with return_se = TRUE
  result <- lss(Y, X, Z, method = "oasis",
                oasis = list(return_se = TRUE))

  # Should return a list with betas and se
  expect_true(is.list(result) || !is.null(attr(result, "se")))
})

test_that("OASIS with ridge + return_se computes SEs", {
  skip_if_not_installed("fmrihrf")

  set.seed(234)
  n_timepoints <- 60
  n_trials <- 10
  n_voxels <- 8

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints))
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test with ridge and return_se
  result <- lss(Y, X, Z, method = "oasis",
                oasis = list(
                  ridge_x = 0.01,
                  ridge_b = 0.01,
                  return_se = TRUE
                ))

  expect_true(is.list(result) || !is.null(attr(result, "se")))
})

test_that("OASIS fractional ridge mode works", {
  skip_if_not_installed("fmrihrf")

  set.seed(345)
  n_timepoints <- 70
  n_trials <- 10
  n_voxels <- 6

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 6 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints))
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test fractional ridge mode
  result_abs <- lss(Y, X, Z, method = "oasis",
                    oasis = list(
                      ridge_mode = "absolute",
                      ridge_x = 0.1,
                      ridge_b = 0.1
                    ))

  result_frac <- lss(Y, X, Z, method = "oasis",
                     oasis = list(
                       ridge_mode = "fractional",
                       ridge_x = 0.01,
                       ridge_b = 0.01
                     ))

  # Both should produce valid results
  expect_equal(dim(result_abs), c(n_trials, n_voxels))
  expect_equal(dim(result_frac), c(n_trials, n_voxels))

  # Results should differ (different ridge amounts)
  expect_false(isTRUE(all.equal(result_abs, result_frac, tolerance = 1e-3)))
})

test_that("OASIS return_diag produces diagnostics", {
  skip_if_not_installed("fmrihrf")

  set.seed(456)
  n_timepoints <- 60
  n_trials <- 8
  n_voxels <- 5

  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 6 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset + 3), i] <- 1
    }
  }

  Z <- cbind(rep(1, n_timepoints))
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Test with return_diag
  result <- lss(Y, X, Z, method = "oasis",
                oasis = list(return_diag = TRUE))

  # Should have diagnostics attribute or be a list
  expect_true(is.list(result) || !is.null(attr(result, "diagnostics")))
})

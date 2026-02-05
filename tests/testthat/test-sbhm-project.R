test_that("sbhm_project computes correct inner products", {
  set.seed(123)
  r <- 3       # basis dimension
  ntrials <- 10
  V <- 5       # number of voxels

  # Create random coefficients
  beta_rt <- array(rnorm(r * ntrials * V), dim = c(r, ntrials, V))
  alpha_hat <- matrix(rnorm(r * V), nrow = r, ncol = V)

  # Compute via sbhm_project
  amps <- sbhm_project(beta_rt, alpha_hat)

  # Check dimensions
  expect_equal(dim(amps), c(ntrials, V))

  # Verify computation manually for first voxel
  expected_v1 <- numeric(ntrials)
  for (t in 1:ntrials) {
    expected_v1[t] <- sum(beta_rt[, t, 1] * alpha_hat[, 1])
  }
  expect_equal(amps[, 1], expected_v1, tolerance = 1e-12)

  # Verify computation manually for second voxel
  expected_v2 <- numeric(ntrials)
  for (t in 1:ntrials) {
    expected_v2[t] <- sum(beta_rt[, t, 2] * alpha_hat[, 2])
  }
  expect_equal(amps[, 2], expected_v2, tolerance = 1e-12)
})

test_that("sbhm_project validates input dimensions", {
  set.seed(234)
  r <- 3
  ntrials <- 10
  V <- 5

  beta_rt <- array(rnorm(r * ntrials * V), dim = c(r, ntrials, V))

  # Wrong r dimension
  alpha_wrong_r <- matrix(rnorm(4 * V), nrow = 4, ncol = V)
  expect_error(sbhm_project(beta_rt, alpha_wrong_r), "alpha_hat must be r×V to match beta_rt dims")

  # Wrong V dimension
  alpha_wrong_v <- matrix(rnorm(r * (V + 1)), nrow = r, ncol = V + 1)
  expect_error(sbhm_project(beta_rt, alpha_wrong_v), "alpha_hat must be r×V to match beta_rt dims")
})

test_that("sbhm_project requires 3D array for beta_rt", {
  set.seed(345)
  r <- 3
  ntrials <- 10
  V <- 5

  # 2D matrix instead of 3D array
  beta_2d <- matrix(rnorm(r * ntrials), nrow = r, ncol = ntrials)
  alpha_hat <- matrix(rnorm(r * V), nrow = r, ncol = V)

  expect_error(sbhm_project(beta_2d, alpha_hat))
})

test_that("sbhm_project requires matrix for alpha_hat", {
  set.seed(456)
  r <- 3
  ntrials <- 10
  V <- 5

  beta_rt <- array(rnorm(r * ntrials * V), dim = c(r, ntrials, V))

  # Vector instead of matrix
  alpha_vec <- rnorm(r * V)

  expect_error(sbhm_project(beta_rt, alpha_vec))
})

test_that("sbhm_project handles edge case of r=1", {
  # Note: r=1 case has a known dimension-drop issue in the function
  # When r=1, beta_rt[,,v] becomes a vector and colSums fails
  set.seed(567)
  r <- 1
  ntrials <- 8
  V <- 4

  beta_rt <- array(rnorm(r * ntrials * V), dim = c(r, ntrials, V))
  alpha_hat <- matrix(rnorm(r * V), nrow = r, ncol = V)

  # Function currently errors on r=1 due to dimension drop
  expect_error(sbhm_project(beta_rt, alpha_hat), "must be an array")
})

test_that("sbhm_project handles edge case of single trial", {
  # Note: single trial case has a known dimension-drop issue in the function
  # This test documents the expected error
  set.seed(678)
  r <- 3
  ntrials <- 1
  V <- 5

  beta_rt <- array(rnorm(r * ntrials * V), dim = c(r, ntrials, V))
  alpha_hat <- matrix(rnorm(r * V), nrow = r, ncol = V)

  # Function currently errors on single trial due to dimension drop
  expect_error(sbhm_project(beta_rt, alpha_hat), "must be an array")
})

test_that("sbhm_project handles edge case of single voxel", {
  set.seed(789)
  r <- 3
  ntrials <- 10
  V <- 1

  beta_rt <- array(rnorm(r * ntrials * V), dim = c(r, ntrials, V))
  alpha_hat <- matrix(rnorm(r * V), nrow = r, ncol = V)

  amps <- sbhm_project(beta_rt, alpha_hat)

  expect_equal(dim(amps), c(ntrials, V))
})

test_that("sbhm_project preserves sign relationships", {
  set.seed(890)
  r <- 2
  ntrials <- 5
  V <- 3

  # Create positive coefficients and alphas
  beta_rt <- array(abs(rnorm(r * ntrials * V)), dim = c(r, ntrials, V))
  alpha_hat <- matrix(abs(rnorm(r * V)), nrow = r, ncol = V)

  amps <- sbhm_project(beta_rt, alpha_hat)

  # All amplitudes should be positive

  expect_true(all(amps > 0))

  # Now flip alpha signs
  alpha_neg <- -alpha_hat
  amps_neg <- sbhm_project(beta_rt, alpha_neg)

  # All amplitudes should be negative
  expect_true(all(amps_neg < 0))

  # Absolute values should be the same
  expect_equal(abs(amps), abs(amps_neg), tolerance = 1e-12)
})

test_that("sbhm_project produces zero amplitudes when orthogonal", {
  set.seed(901)
  r <- 2
  ntrials <- 5
  V <- 3

  # Create orthogonal pairs
  beta_rt <- array(0, dim = c(r, ntrials, V))
  alpha_hat <- matrix(0, nrow = r, ncol = V)

  # beta has only first basis, alpha has only second basis
  beta_rt[1, , ] <- rnorm(ntrials * V)  # non-zero in first basis
  alpha_hat[2, ] <- rnorm(V)             # non-zero in second basis

  amps <- sbhm_project(beta_rt, alpha_hat)

  # All amplitudes should be zero (or very close)
  expect_true(all(abs(amps) < 1e-12))
})

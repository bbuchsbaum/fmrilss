# Extended tests for oasis_hrf_recovery.R functions

test_that("generate_rapid_design returns valid onsets", {
  result <- generate_rapid_design(n_events = 20, total_time = 200,
                                  min_isi = 2, max_isi = 4, seed = 123)

  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
  expect_true(all(result >= 0))
  expect_true(all(result < 200 - 20))  # Within bounds
})

test_that("generate_rapid_design respects ISI bounds", {
  result <- generate_rapid_design(n_events = 50, total_time = 300,
                                  min_isi = 3, max_isi = 5, seed = 234)

  # Calculate ISIs
  isis <- diff(result)

  # All ISIs should be within bounds (with small tolerance for cumsum effects)
  expect_true(all(isis >= 2.9))  # Slightly below min_isi due to numerical issues
  expect_true(all(isis <= 5.1))  # Slightly above max_isi
})

test_that("generate_rapid_design is reproducible with seed", {
  result1 <- generate_rapid_design(n_events = 15, seed = 456)
  result2 <- generate_rapid_design(n_events = 15, seed = 456)

  expect_equal(result1, result2)
})

test_that("generate_rapid_design handles edge cases", {
  # Very short total time
  result_short <- generate_rapid_design(n_events = 5, total_time = 50,
                                        min_isi = 2, max_isi = 3, seed = 567)
  expect_true(length(result_short) > 0)

  # Single event
  result_single <- generate_rapid_design(n_events = 1, total_time = 100,
                                         min_isi = 2, max_isi = 4, seed = 678)
  expect_equal(length(result_single), 1)
  expect_equal(result_single, 5)  # Starts at 5s
})

test_that("generate_lwu_data requires fmrihrf", {
  skip_if_not_installed("fmrihrf")

  onsets <- c(5, 15, 25, 35)
  result <- generate_lwu_data(onsets, tau = 6, sigma = 2.5, rho = 0.35,
                              TR = 1, total_time = 100, n_voxels = 5, seed = 123)

  expect_true(is.list(result))
  expect_true("Y" %in% names(result))
  expect_true("X" %in% names(result))
  expect_true("true_hrf" %in% names(result))
  expect_true("true_betas" %in% names(result))
  expect_true("sframe" %in% names(result))
})

test_that("generate_lwu_data returns correct dimensions", {
  skip_if_not_installed("fmrihrf")

  onsets <- c(5, 15, 25)
  n_voxels <- 3
  total_time <- 50
  TR <- 1

  result <- generate_lwu_data(onsets, total_time = total_time, TR = TR,
                              n_voxels = n_voxels, seed = 234)

  n_time_expected <- length(seq(0, total_time, by = TR))

  expect_equal(nrow(result$Y), n_time_expected)
  expect_equal(ncol(result$Y), n_voxels)
  expect_equal(nrow(result$X), n_time_expected)
  expect_equal(ncol(result$X), length(onsets))
  expect_equal(nrow(result$true_betas), length(onsets))
  expect_equal(ncol(result$true_betas), n_voxels)
})

test_that("generate_lwu_data uses provided amplitudes", {
  skip_if_not_installed("fmrihrf")

  onsets <- c(5, 15)
  amplitudes <- c(2.0, 0.5)

  result <- generate_lwu_data(onsets, amplitudes = amplitudes,
                              total_time = 50, n_voxels = 2, seed = 345)

  expect_equal(result$amplitudes, amplitudes)
})

test_that("generate_lwu_data scalar amplitude is expanded", {
  skip_if_not_installed("fmrihrf")

  onsets <- c(5, 15, 25)
  amplitude <- 1.5

  result <- generate_lwu_data(onsets, amplitudes = amplitude,
                              total_time = 60, n_voxels = 2, seed = 456)

  expect_equal(result$amplitudes, rep(amplitude, 3))
})

test_that("generate_lwu_data stores hrf_params", {
  skip_if_not_installed("fmrihrf")

  onsets <- c(5, 15)
  tau <- 7
  sigma <- 3
  rho <- 0.4

  result <- generate_lwu_data(onsets, tau = tau, sigma = sigma, rho = rho,
                              total_time = 50, n_voxels = 2, seed = 567)

  expect_equal(result$hrf_params$tau, tau)
  expect_equal(result$hrf_params$sigma, sigma)
  expect_equal(result$hrf_params$rho, rho)
})

test_that("create_lwu_grid creates grid with correct dimensions", {
  result <- create_lwu_grid(tau_range = c(4, 8),
                            sigma_range = c(1.5, 3.5),
                            rho_range = c(0.1, 0.5),
                            n_tau = 3, n_sigma = 2, n_rho = 2)

  expect_true(is.list(result))
  expect_true("hrfs" %in% names(result))
  expect_true("parameters" %in% names(result))

  # Grid should have n_tau * n_sigma * n_rho rows
  expected_n <- 3 * 2 * 2  # 12
  expect_equal(length(result$hrfs), expected_n)
  expect_equal(nrow(result$parameters), expected_n)
})

test_that("create_lwu_grid HRFs have correct attributes", {
  result <- create_lwu_grid(tau_range = c(5, 7), n_tau = 3,
                            sigma_range = c(2, 3), n_sigma = 2,
                            rho_range = c(0.2, 0.4), n_rho = 2)

  # Check first HRF has expected attributes
  hrf1 <- result$hrfs[[1]]
  expect_true(is.function(hrf1))
  expect_true(!is.null(attr(hrf1, "tau")))
  expect_true(!is.null(attr(hrf1, "sigma")))
  expect_true(!is.null(attr(hrf1, "rho")))
  expect_equal(attr(hrf1, "span"), 30)
})

test_that("create_lwu_grid parameters match HRF attributes", {
  result <- create_lwu_grid(tau_range = c(4, 6), n_tau = 2,
                            sigma_range = c(2, 3), n_sigma = 2,
                            rho_range = c(0.3, 0.4), n_rho = 2)

  # Check that each HRF's attributes match the parameter grid
  for (i in seq_len(nrow(result$parameters))) {
    hrf <- result$hrfs[[i]]
    expect_equal(attr(hrf, "tau"), result$parameters$tau[i])
    expect_equal(attr(hrf, "sigma"), result$parameters$sigma[i])
    expect_equal(attr(hrf, "rho"), result$parameters$rho[i])
  }
})

test_that("create_lwu_grid HRFs are callable", {
  skip_if_not_installed("fmrihrf")

  result <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)

  # Try calling each HRF
  times <- seq(0, 20, by = 1)
  for (i in seq_along(result$hrfs)) {
    hrf_values <- result$hrfs[[i]](times)
    expect_true(is.numeric(hrf_values))
    expect_equal(length(hrf_values), length(times))
  }
})

test_that("fit_oasis_grid validates inputs", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(100 * 5), 100, 5)
  onsets <- c(5, 25, 45, 65)
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)
  hrf_grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)

  # Should not error with valid inputs
  expect_no_error({
    # Just run - may take a moment but should work
    result <- fit_oasis_grid(Y, onsets, sframe, hrf_grid,
                             ridge_x = 0.01, ridge_b = 0.01)
  })
})

test_that("fit_oasis_grid returns expected structure", {
  skip_if_not_installed("fmrihrf")

  set.seed(123)
  Y <- matrix(rnorm(100 * 3), 100, 3)
  onsets <- c(10, 30, 50, 70)
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)
  hrf_grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)

  result <- fit_oasis_grid(Y, onsets, sframe, hrf_grid)

  expect_true(is.list(result))
  expect_true("best_idx" %in% names(result))
  expect_true("best_params" %in% names(result))
  expect_true("best_hrf" %in% names(result))
  expect_true("beta" %in% names(result))
  expect_true("scores" %in% names(result))
  expect_true("grid" %in% names(result))

  expect_true(result$best_idx >= 1)
  expect_true(result$best_idx <= length(hrf_grid$hrfs))
})

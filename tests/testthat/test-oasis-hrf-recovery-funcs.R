test_that("generate_rapid_design creates valid event onsets", {
  onsets <- generate_rapid_design(n_events = 20, total_time = 200, min_isi = 2, max_isi = 4, seed = 123)

  # Check we got some events (may be fewer than requested due to time constraints)
  expect_gt(length(onsets), 0)

  # Check onsets are sorted
  expect_equal(onsets, sort(onsets))

  # Check onsets respect total_time (leaving room for HRF tail)
  expect_true(all(onsets < (200 - 20)))

  # Check first onset is at 5s
  expect_equal(onsets[1], 5)
})

test_that("generate_rapid_design respects ISI constraints", {
  set.seed(456)
  onsets <- generate_rapid_design(n_events = 15, total_time = 150, min_isi = 3, max_isi = 5, seed = 456)

  # Calculate actual ISIs
  isis <- diff(onsets)

  # All ISIs should be within specified range
  expect_true(all(isis >= 3))
  expect_true(all(isis <= 5))
})

test_that("generate_rapid_design uses seed for reproducibility", {
  onsets1 <- generate_rapid_design(n_events = 10, seed = 999)
  onsets2 <- generate_rapid_design(n_events = 10, seed = 999)

  expect_equal(onsets1, onsets2)
})

test_that("generate_lwu_data creates valid synthetic data", {
  skip_if_not_installed("fmrihrf")

  onsets <- generate_rapid_design(n_events = 10, total_time = 100, seed = 123)
  data <- generate_lwu_data(
    onsets,
    tau = 6, sigma = 2.5, rho = 0.35,
    TR = 2.0, total_time = 100, n_voxels = 5,
    noise_sd = 0.1, seed = 456
  )

  # Check output structure
  expect_true(is.list(data))
  expect_true("Y" %in% names(data))
  expect_true("X" %in% names(data))
  expect_true("true_hrf" %in% names(data))
  expect_true("true_betas" %in% names(data))
  expect_true("onsets" %in% names(data))
  expect_true("sframe" %in% names(data))

  # Check dimensions
  n_time <- length(data$time_points)
  n_trials <- length(data$onsets)
  n_voxels <- 5
  expect_equal(dim(data$Y), c(n_time, n_voxels))
  expect_equal(dim(data$X), c(n_time, n_trials))
  expect_equal(dim(data$true_betas), c(n_trials, n_voxels))
})

test_that("generate_lwu_data accepts custom amplitudes", {
  skip_if_not_installed("fmrihrf")

  onsets <- generate_rapid_design(n_events = 5, total_time = 60, seed = 111)

  # Scalar amplitude
  data1 <- generate_lwu_data(onsets, total_time = 60, n_voxels = 3, amplitudes = 2.0, seed = 222)
  expect_equal(data1$amplitudes, rep(2.0, length(onsets)))

  # Vector amplitude
  amps <- c(1, 2, 3, 4, 5)
  data2 <- generate_lwu_data(onsets, total_time = 60, n_voxels = 3, amplitudes = amps, seed = 333)
  expect_equal(data2$amplitudes, amps)
})

test_that("create_lwu_grid generates parameter grid", {
  grid <- create_lwu_grid(
    tau_range = c(4, 8), n_tau = 3,
    sigma_range = c(1.5, 3.5), n_sigma = 2,
    rho_range = c(0.1, 0.5), n_rho = 2
  )

  # Check structure
  expect_true(is.list(grid))
  expect_true("hrfs" %in% names(grid))
  expect_true("parameters" %in% names(grid))

  # Check grid size
  expected_size <- 3 * 2 * 2  # n_tau * n_sigma * n_rho
  expect_equal(length(grid$hrfs), expected_size)
  expect_equal(nrow(grid$parameters), expected_size)

  # Check parameter ranges
  expect_true(all(grid$parameters$tau >= 4 & grid$parameters$tau <= 8))
  expect_true(all(grid$parameters$sigma >= 1.5 & grid$parameters$sigma <= 3.5))
  expect_true(all(grid$parameters$rho >= 0.1 & grid$parameters$rho <= 0.5))

  # Check HRFs are functions with attributes
  expect_true(is.function(grid$hrfs[[1]]))
  expect_true(!is.null(attr(grid$hrfs[[1]], "tau")))
  expect_true(!is.null(attr(grid$hrfs[[1]], "sigma")))
  expect_true(!is.null(attr(grid$hrfs[[1]], "rho")))
})

test_that("create_lwu_grid produces callable HRF functions", {
  skip_if_not_installed("fmrihrf")

  grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)

  # Each HRF should be evaluable
  times <- seq(0, 30, by = 1)
  for (i in 1:length(grid$hrfs)) {
    hrf_vals <- grid$hrfs[[i]](times)
    expect_true(is.numeric(hrf_vals))
    expect_equal(length(hrf_vals), length(times))
    expect_false(any(is.na(hrf_vals)))
  }
})

test_that("fit_oasis_grid runs and returns expected structure", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()  # Skip on CRAN due to computational cost

  # Create small test data
  onsets <- generate_rapid_design(n_events = 5, total_time = 40, seed = 123)
  data <- generate_lwu_data(
    onsets, tau = 6, sigma = 2.5, rho = 0.3,
    TR = 2.0, total_time = 40, n_voxels = 2,
    noise_sd = 0.5, seed = 456
  )

  # Create small grid
  grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 1)

  # Fit with grid search
  result <- fit_oasis_grid(
    Y = data$Y,
    onsets = data$onsets,
    sframe = data$sframe,
    hrf_grid = grid,
    ridge_x = 0.1,
    ridge_b = 0.1
  )

  # Check structure
  expect_true(is.list(result))
  expect_true("best_idx" %in% names(result))
  expect_true("best_params" %in% names(result))
  expect_true("beta" %in% names(result))
  expect_true("scores" %in% names(result))

  # Check best_idx is valid
  expect_true(result$best_idx >= 1 && result$best_idx <= length(grid$hrfs))

  # Check scores length matches grid size
  expect_equal(length(result$scores), length(grid$hrfs))
})

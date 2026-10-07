test_that("LWU simulation preserves fractional physical event times", {
  skip_if_not_installed("fmrihrf")
  onsets <- c(5.23, 29.67, 62.41)
  amps <- c(.7, 1.3, 2)
  sim <- generate_lwu_data(onsets, TR = 1.7, total_time = 110,
                           amplitudes = amps, noise_sd = 0, n_voxels = 2, seed = 91)
  # Independent causal two-Gaussian LWU formula, not the regressor builder.
  raw <- function(t) {
    z <- exp(-(t - 6)^2 / (2 * 2.5^2)) -
      .35 * exp(-(t - 11)^2 / (2 * 4^2))
    z[t < 0 | t > 30] <- 0
    z
  }
  height <- max(abs(raw(seq(0, 30, by = .1))))
  oracle <- vapply(seq_along(onsets), function(j) {
    amps[j] * raw(sim$time_points - onsets[j]) / height
  }, numeric(length(sim$time_points)))
  # Onset discretization and interpolation at the production 0.1s precision.
  expect_lt(max(abs(sim$X - oracle)), .004)
  expect_equal(unname(sim$signal), unname(sim$X %*% sim$true_betas))
  expect_equal(sim$Y, sim$signal)
  expect_equal(sim$hrf_times, seq(0, 30, by = 1.7))
})

test_that("HRF selection uses joint profile residuals independently of LSS ridge", {
  skip_if_not_installed("fmrihrf")
  onsets <- c(5.23, 34.67, 66.41, 99.18)
  sim <- generate_lwu_data(onsets, tau = 6, sigma = 2.5, rho = .35,
                           TR = 1.7, total_time = 140, amplitudes = c(.7, 2, 1, 1.4),
                           noise_sd = 0, n_voxels = 3, seed = 92)
  grid <- create_lwu_grid(c(5, 7), c(2, 3), c(.2, .5), 3, 3, 3)
  fit <- fit_oasis_grid(sim$Y, onsets, sim$sframe, grid)
  expect_equal(unname(unlist(fit$best_params)), c(6, 2.5, .35))
  expect_equal(max(fit$scores), 1, tolerance = 1e-12)
  shifted <- sweep(sim$Y * 5000, 2, c(100, -300, 10000), "+")
  other <- fit_oasis_grid(shifted, onsets, sim$sframe, grid,
                          ridge_x = .3, ridge_b = .4)
  expect_equal(other$scores, fit$scores, tolerance = 1e-12)
  expect_equal(other$best_idx, fit$best_idx)
  perm <- fit_oasis_grid(sim$Y, rev(onsets), sim$sframe, grid)
  expect_equal(perm$scores, fit$scores, tolerance = 1e-12)
  expect_equal(perm$best_idx, fit$best_idx)

  # Duplicate event columns: independent SVD oracle handles deficient rank.
  dup_onsets <- c(onsets, onsets[1])
  dup <- fit_oasis_grid(sim$Y, dup_onsets, sim$sframe, grid)
  oracle <- vapply(seq_along(grid$hrfs), function(i) {
    X <- fmrilss:::.oasis_build_X_from_events(list(sframe = sim$sframe,
      cond = list(onsets = dup_onsets, hrf = structure(grid$hrfs[[i]],
        class = c("hrf", "function")), span = 30)))$X_trials
    s <- svd(cbind(1, X))
    U <- s$u[, s$d > max(s$d) * 1e-10, drop = FALSE]
    rss <- sum((sim$Y - U %*% crossprod(U, sim$Y))^2)
    1 - rss / sum(sweep(sim$Y, 2, colMeans(sim$Y))^2)
  }, numeric(1))
  expect_equal(dup$scores, oracle, tolerance = 1e-12)
  expect_error(fit_oasis_grid(sim$Y * 0, onsets, sim$sframe, grid),
               "within-voxel variation")
  invalid <- grid
  invalid$parameters$sigma[] <- 0
  expect_error(fit_oasis_grid(sim$Y, onsets, sim$sframe, invalid),
               "No HRF candidate")
})

test_that("recovery metrics use supplied waveform, time units, and effective betas", {
  skip_if_not_installed("fmrihrf")
  sim <- generate_lwu_data(c(5.23, 34.67, 66.41), total_time = 110, TR = 1.7,
                           amplitudes = c(.7, 2, 1), noise_sd = 0, n_voxels = 2, seed = 93)
  grid <- create_lwu_grid(c(6, 6), c(2.5, 2.5), c(.35, .35), 1, 1, 1)
  res <- suppressMessages(compare_hrf_recovery(sim, grid))
  expect_equal(res$true_betas, sim$true_betas * sim$amplitudes)
  met <- calculate_recovery_metrics(res, sim$true_hrf)
  expect_equal(met$mse[1], 0, tolerance = 1e-12)
  expect_true(is.na(met$beta_correlation[3]))
  # Supply the canonical curve as truth to prove the argument is honored.
  truth <- fmrihrf::evaluate(fmrihrf::HRF_SPMG1, sim$hrf_times)
  canonical <- calculate_recovery_metrics(res, truth)
  expect_equal(canonical$mse[2], 0, tolerance = 1e-12)
  expect_equal(canonical$peak_time_error[2], 0)
  expect_equal(canonical$mse[2], canonical$mse[3], tolerance = 1e-12)
})

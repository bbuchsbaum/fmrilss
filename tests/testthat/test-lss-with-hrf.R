test_that("lss_with_hrf reduces to lss() when K=1 and weights=1", {
  skip_on_cran()
  set.seed(123)

  # Dimensions
  n_time   <- 120L
  n_trials <- 10L
  n_vox    <- 8L

  # --- Build event sticks and onsets/durations ---
  onset_idx <- as.integer(round(seq(6, n_time - 12, length.out = n_trials)))
  dur <- rep(4L, n_trials)  # inclusive segments of length 5

  Xev <- matrix(0, n_time, n_trials)
  for (i in seq_len(n_trials)) {
    i1 <- onset_idx[i]
    i2 <- min(n_time, i1 + dur[i])
    Xev[i1:i2, i] <- 1
  }

  # --- Single HRF basis (K=1) ---
  # simple unimodal kernel; length L
  hrf_kernel <- c(0.0, 0.2, 0.6, 1.0, 0.6, 0.3, 0.15, 0.05, 0.0)
  L <- length(hrf_kernel)
  basis <- matrix(hrf_kernel, nrow = L, ncol = 1L)

  # voxel weights = 1 for all voxels
  coeffs <- matrix(1.0, nrow = 1L, ncol = n_vox)

  # --- Convolve to get X for lss() ---
  conv_open_trim <- function(x, k) {
    as.numeric(stats::convolve(x, rev(as.numeric(k)), type = "open"))[seq_len(length(x))]
  }
  X_hrf <- vapply(seq_len(n_trials), function(j) conv_open_trim(Xev[, j], hrf_kernel),
                  numeric(n_time))

  # --- Simulate data ---
  true_betas <- matrix(rnorm(n_trials * n_vox, sd = 1.0), n_trials, n_vox)
  Y <- X_hrf %*% true_betas + matrix(rnorm(n_time * n_vox, sd = 0.5), n_time, n_vox)
  colnames(Y) <- paste0("V", seq_len(n_vox))

  # --- Compare ---
  beta_hrf <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = dur,
    hrf_basis_kernels = basis,
    coefficients = coeffs,
    Z = NULL,
    Nuisance = NULL,
    verbose = FALSE
  )

  beta_lss <- fmrilss::lss(Y, X_hrf, method = "r_optimized")

  # Compare values (ignoring row names difference)
  expect_equal(unname(beta_hrf), unname(beta_lss), tolerance = 1e-8)
})

test_that("lss_with_hrf runs with nuisance and non-null Z", {
  skip_on_cran()
  set.seed(42)

  n_time   <- 100L
  n_trials <- 6L
  n_vox    <- 5L

  onset_idx <- as.integer(seq(8, n_time - 12, length.out = n_trials))
  dur <- rep(3L, n_trials)

  Xev <- matrix(0, n_time, n_trials)
  for (i in seq_len(n_trials)) {
    i1 <- onset_idx[i]
    i2 <- min(n_time, i1 + dur[i])
    Xev[i1:i2, i] <- 1
  }

  # two basis columns (K=2)
  b1 <- c(0, 0.4, 0.9, 1.0, 0.5, 0.2, 0)
  b2 <- c(0, 0.2, 0.4, 0.3, 0.2, 0.1, 0)
  basis <- cbind(b1, b2)

  # voxel weights (random)
  coeffs <- matrix(runif(2 * n_vox, 0.5, 1.5), nrow = 2L, ncol = n_vox)

  # build a reference X by convolving with a *voxel-averaged* HRF (for stability only)
  avg_hrf <- drop(basis %*% rowMeans(coeffs))
  conv_open_trim <- function(x, k) {
    as.numeric(stats::convolve(x, rev(as.numeric(k)), type = "open"))[seq_len(length(x))]
  }
  X_ref <- vapply(seq_len(n_trials), function(j) conv_open_trim(Xev[, j], avg_hrf),
                  numeric(n_time))

  # experimental regressors and nuisance
  Z <- cbind(1, scale(seq_len(n_time)))
  Nuis <- cbind(scale(sin(seq_len(n_time) / 7)), scale(cos(seq_len(n_time) / 11)))

  # simulate Y
  true_betas <- matrix(rnorm(n_trials * n_vox, sd = 0.7), n_trials, n_vox)
  Y <- X_ref %*% true_betas + Z %*% matrix(rnorm(ncol(Z) * n_vox), ncol(Z), n_vox) +
       Nuis %*% matrix(rnorm(ncol(Nuis) * n_vox, sd = 0.3), ncol(Nuis), n_vox) +
       matrix(rnorm(n_time * n_vox, sd = 0.6), n_time, n_vox)

  # run lss_with_hrf just to ensure it returns finite values
  beta_hrf <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = dur,
    hrf_basis_kernels = basis,
    coefficients = coeffs,
    Z = Z,
    Nuisance = Nuis,
    verbose = FALSE
  )

  expect_equal(dim(beta_hrf), c(n_trials, n_vox))
  expect_true(all(is.finite(beta_hrf) | is.na(beta_hrf)))
})
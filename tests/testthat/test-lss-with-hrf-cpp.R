test_that("cpp backend matches R backend on a small problem", {
  skip_on_cran()
  # If the C++ symbol isn't built yet, skip
  sym_avail <- FALSE
  try({
    get("lss_engine_vox_hrf_cpp", envir = asNamespace("fmrilss"))
    sym_avail <- TRUE
  }, silent = TRUE)
  if (!sym_avail) skip("C++ backend not available in this build")

  set.seed(2025)
  n_time <- 90L; n_trials <- 8L; n_vox <- 6L
  onset_idx <- as.integer(seq(6, n_time - 12, length.out = n_trials))
  dur <- rep(2L, n_trials)

  # Event matrix
  Xev <- matrix(0, n_time, n_trials)
  for (i in seq_len(n_trials)) {
    i1 <- onset_idx[i]; i2 <- min(n_time, i1 + dur[i])
    Xev[i1:i2, i] <- 1
  }

  # Basis (K=3) and random voxel weights
  b1 <- c(0, 0.3, 0.9, 1.0, 0.6, 0.2, 0)
  b2 <- c(0, 0.1, 0.4, 0.5, 0.3, 0.1, 0)
  b3 <- c(0, 0.2, 0.5, 0.4, 0.25, 0.1, 0)
  basis <- cbind(b1, b2, b3)
  coeffs <- matrix(runif(3 * n_vox, 0.5, 1.5), nrow = 3L, ncol = n_vox)

  # Build a reference X using average HRF to simulate Y
  avg_hrf <- drop(basis %*% rowMeans(coeffs))
  conv_open_trim <- function(x, k) {
    as.numeric(stats::convolve(x, rev(as.numeric(k)), type = "open"))[seq_len(length(x))]
  }
  X_ref <- vapply(seq_len(n_trials), function(j) conv_open_trim(Xev[, j], avg_hrf),
                  numeric(n_time))

  # Regressors and nuisance
  Z <- cbind(1, scale(seq_len(n_time)))
  Nuis <- cbind(scale(sin(seq_len(n_time) / 9)), scale(cos(seq_len(n_time) / 10)))

  # Simulate data
  true_betas <- matrix(rnorm(n_trials * n_vox, sd = 0.6), n_trials, n_vox)
  Y <- X_ref %*% true_betas + Z %*% matrix(rnorm(ncol(Z) * n_vox), ncol(Z), n_vox) +
       matrix(rnorm(n_time * n_vox, sd = 0.5), n_time, n_vox) +
       Nuis %*% matrix(rnorm(ncol(Nuis) * n_vox, sd = 0.2), ncol(Nuis), n_vox)

  # Run both backends
  beta_r <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = dur,
    hrf_basis_kernels = basis,
    coefficients = coeffs,
    Z = Z,
    Nuisance = Nuis,
    method = "r",
    verbose = FALSE
  )

  beta_cpp <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = dur,
    hrf_basis_kernels = basis,
    coefficients = coeffs,
    Z = Z,
    Nuisance = Nuis,
    method = "cpp",
    verbose = FALSE
  )

  expect_equal(beta_cpp, beta_r, tolerance = 1e-8)
})
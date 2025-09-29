test_that("lss_with_hrf C++ honors durations and matches R backend closely", {
  skip_on_cran()

  # If the C++ symbol isn't built yet, skip
  sym_avail <- FALSE
  try({
    get("lss_engine_vox_hrf", envir = asNamespace("fmrilss"))
    sym_avail <- TRUE
  }, silent = TRUE)
  if (!sym_avail) skip("C++ engine not available in this build")

  set.seed(2121)
  n_time <- 120L
  n_vox  <- 5L
  onsets <- c(10, 35, 60, 85)
  durations <- c(3, 2, 4, 1)  # non-zero durations

  # Events data frame (seconds, TR=1)
  events <- data.frame(
    onset = onsets,
    duration = durations,
    condition = "cond"
  )

  # HRF basis (K=1) and coefficients (all ones)
  basis <- fmrihrf::HRF_SPMG1
  coeffs <- matrix(1, nrow = 1L, ncol = n_vox)
  vhrf <- list(coefficients = coeffs, basis = basis, conditions = "cond")
  class(vhrf) <- "VoxelHRF"

  # Build synthetic data using the TR-level design (durations respected)
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  times  <- fmrihrf::samples(sframe, global = TRUE)
  rset   <- fmrihrf::regressor_set(onsets = onsets,
                                   fac = factor(seq_along(onsets)),
                                   hrf = basis,
                                   duration = durations,
                                   span = 30,
                                   summate = FALSE)
  X_trials <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  if (inherits(X_trials, "Matrix")) X_trials <- as.matrix(X_trials)

  set.seed(313)
  betas_true <- matrix(rnorm(length(onsets) * n_vox, sd = 0.5), length(onsets), n_vox)
  Y_noise <- matrix(rnorm(n_time * n_vox, sd = 0.05), n_time, n_vox)
  Y <- X_trials %*% betas_true + Y_noise

  # Compute betas via R engine
  beta_r <- lss_with_hrf(
    Y = Y, events = events, hrf_estimates = vhrf,
    nuisance_regs = NULL, engine = "R", verbose = FALSE
  )

  # Compute betas via C++ engine
  beta_cpp_obj <- lss_with_hrf(
    Y = Y, events = events, hrf_estimates = vhrf,
    nuisance_regs = NULL, engine = "C++", verbose = FALSE
  )
  beta_cpp <- as.matrix(beta_cpp_obj$betas[])

  # The two backends should be very close for integer onsets/durations
  expect_equal(dim(beta_cpp), dim(beta_r))
  expect_lt(max(abs(beta_cpp - beta_r)), 5e-3)
})

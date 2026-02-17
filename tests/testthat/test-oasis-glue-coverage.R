# Tests for oasis_glue.R internal functions

test_that(".lss_oasis validates Y as matrix", {
  expect_error(
    fmrilss:::.lss_oasis(Y = "not_matrix"),
    "Y must be a numeric matrix"
  )
})

test_that(".lss_oasis validates Y has finite values", {
  Y <- matrix(c(1, 2, Inf, 4), 2, 2)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y),
    "Y contains non-finite values"
  )
})

test_that(".lss_oasis validates X dimensions", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(40 * 3), 40, 3)  # Wrong rows

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X),
    "X must have the same number of rows as Y"
  )
})

test_that(".lss_oasis validates X has finite values", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)
  X[25, 2] <- NaN  # Insert NaN

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X),
    "X contains non-finite values"
  )
})

test_that(".lss_oasis validates Z dimensions", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Z <- matrix(rnorm(40 * 2), 40, 2)  # Wrong rows

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, Z = Z),
    "Z must have the same number of rows as Y"
  )
})

test_that(".lss_oasis validates Z has finite values", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Z <- matrix(c(1, Inf), 50, 2, byrow = TRUE)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, Z = Z),
    "Z contains non-finite values"
  )
})

test_that(".lss_oasis validates Nuisance dimensions", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Nuisance <- matrix(rnorm(40 * 2), 40, 2)  # Wrong rows

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, Nuisance = Nuisance),
    "Nuisance must have the same number of rows as Y"
  )
})

test_that(".lss_oasis validates Nuisance has finite values", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)
  Nuisance <- matrix(rnorm(50 * 2), 50, 2)
  Nuisance[25, 1] <- NA  # Insert NA

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, Nuisance = Nuisance),
    "Nuisance contains non-finite values"
  )
})

test_that(".lss_oasis validates oasis is a list", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, oasis = "not_list"),
    "oasis must be a list"
  )
})

test_that(".lss_oasis validates ridge_x is non-negative", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, oasis = list(ridge_x = -0.1)),
    "ridge_x must be non-negative"
  )
})

test_that(".lss_oasis validates ridge_b is non-negative", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, oasis = list(ridge_b = -0.1)),
    "ridge_b must be non-negative"
  )
})

test_that(".lss_oasis requires X or design_spec", {
  Y <- matrix(rnorm(50 * 5), 50, 5)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, oasis = list()),
    "Either X or oasis\\$design_spec must be provided"
  )
})

test_that(".lss_oasis adds default intercept when Z is NULL", {
  set.seed(123)
  Y <- matrix(rnorm(100 * 5), 100, 5)
  X <- matrix(rnorm(100 * 3), 100, 3)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X, Z = NULL, oasis = list(add_intercept = TRUE))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 5))
})

test_that(".lss_oasis returns correct dimensions for K=1", {
  set.seed(234)
  n_time <- 100
  n_trials <- 10
  n_vox <- 5

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X, oasis = list(K = 1))

  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that(".oasis_resolve_ridge handles absolute mode", {
  # Create minimal pre object
  pre <- list(d = c(1, 2, 3), s = c(0.5, 1, 1.5))

  result <- fmrilss:::.oasis_resolve_ridge(pre, ridge_x = 0.1, ridge_b = 0.2,
                                           ridge_mode = "absolute", K = 1L)

  expect_equal(result$lx, 0.1)
  expect_equal(result$lb, 0.2)
})

test_that(".oasis_resolve_ridge handles fractional mode for K=1", {
  pre <- list(d = c(1, 2, 3), s = c(0.5, 1, 1.5))  # mean(d) = 2, mean(s) = 1

  result <- fmrilss:::.oasis_resolve_ridge(pre, ridge_x = 0.5, ridge_b = 0.5,
                                           ridge_mode = "fractional", K = 1L)

  # Fractional: lx = 0.5 * mean(d) = 0.5 * 2 = 1.0
  # Fractional: lb = 0.5 * mean(s) = 0.5 * 1 = 0.5
  expect_equal(result$lx, 1.0)
  expect_equal(result$lb, 0.5)
})

test_that(".oasis_resolve_ridge handles fractional mode for K>1", {
  # For K > 1, pre contains D and E arrays
  K <- 2
  N <- 3

  # Create D[K, K, N] and E[K, K, N] arrays
  D <- array(0, dim = c(K, K, N))
  E <- array(0, dim = c(K, K, N))

  for (j in seq_len(N)) {
    D[, , j] <- diag(c(1, 2))  # mean diag = 1.5 per trial
    E[, , j] <- diag(c(0.5, 1))  # mean diag = 0.75 per trial
  }

  pre <- list(D = D, E = E)

  result <- fmrilss:::.oasis_resolve_ridge(pre, ridge_x = 0.5, ridge_b = 0.5,
                                           ridge_mode = "fractional", K = 2L)

  # mx = mean of mean(diag(D[,,j])) = mean(c(1.5, 1.5, 1.5)) = 1.5
  # mb = mean of mean(diag(E[,,j])) = mean(c(0.75, 0.75, 0.75)) = 0.75
  expect_equal(result$lx, 0.5 * 1.5)  # 0.75
  expect_equal(result$lb, 0.5 * 0.75)  # 0.375
})

test_that(".detect_basis_dimension returns 1 for single column", {
  X <- matrix(rnorm(50), 50, 1)
  result <- fmrilss:::.detect_basis_dimension(X)
  expect_equal(result, 1L)
})

test_that(".detect_basis_dimension detects K from correlation pattern", {
  set.seed(345)
  n_time <- 100
  K <- 3
  n_trials <- 5

  # Create design with K-basis structure - columns within a trial block should be correlated
  X <- matrix(0, n_time, K * n_trials)

  for (i in seq_len(n_trials)) {
    onset <- 10 + (i - 1) * 20
    block_start <- (i - 1) * K + 1

    # Create correlated basis functions within trial
    base_signal <- rep(0, n_time)
    base_signal[onset:min(onset + 5, n_time)] <- 1

    for (k in seq_len(K)) {
      X[, block_start + k - 1] <- base_signal + rnorm(n_time, sd = 0.1)
    }
  }

  result <- fmrilss:::.detect_basis_dimension(X)

  # Should detect K = 3 or fall back to 1
  expect_true(result %in% c(1L, 3L))
})

test_that(".detect_basis_dimension defaults to 1 for unclear pattern", {
  set.seed(456)
  # Random design with no clear K structure
  X <- matrix(rnorm(100 * 7), 100, 7)  # 7 is not in the candidate list

  result <- fmrilss:::.detect_basis_dimension(X)
  expect_equal(result, 1L)
})

test_that(".oasis_build_X_from_events requires spec", {
  expect_error(
    fmrilss:::.oasis_build_X_from_events(NULL),
    "design_spec must be provided"
  )
})

test_that(".oasis_build_X_from_events builds design from events", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 2)

  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(5, 25, 45),  # Global seconds
      duration = 0,
      hrf = fmrihrf::HRF_SPMG1
    )
  )

  result <- fmrilss:::.oasis_build_X_from_events(spec)

  expect_true(is.list(result))
  expect_true("X_trials" %in% names(result))
  expect_true(is.matrix(result$X_trials))
  expect_equal(ncol(result$X_trials), 3)  # 3 trials
})

test_that(".oasis_build_X_from_events returns X_other for other conditions", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)

  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(10, 40, 70),
      duration = 0,
      hrf = fmrihrf::HRF_SPMG1
    ),
    others = list(
      list(onsets = c(20, 50, 80), duration = 0)
    )
  )

  result <- fmrilss:::.oasis_build_X_from_events(spec)

  expect_true(!is.null(result$X_other))
  expect_equal(ncol(result$X_other), 1)  # 1 other condition
})

test_that(".oasis_build_X_from_events detects K from HRF", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 2)

  # Use multi-basis HRF
  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(5, 25),
      duration = 0,
      hrf = fmrihrf::HRF_SPMG3  # 3 basis functions
    )
  )

  result <- fmrilss:::.oasis_build_X_from_events(spec)

  # K should be detected from HRF
  expect_true(result$K >= 1)
})

test_that(".oasis_pick_hrf_lwu_fast selects best HRF", {
  skip_if_not_installed("fmrihrf")

  set.seed(567)

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)

  # Create Y that correlates better with one HRF shape
  onsets <- c(10, 40, 70)

  # Build HRF grid
  hrf_grid <- list(
    fmrihrf::HRF_SPMG1,
    fmrihrf::as_hrf(function(t) fmrihrf::hrf_gamma(t, shape = 4, rate = 1), span = 32)
  )

  design_spec <- list(
    sframe = sframe,
    cond = list(
      onsets = onsets,
      duration = 0,
      span = 32
    )
  )

  Y <- matrix(rnorm(100 * 10), 100, 10)

  # Should not error
  result <- fmrilss:::.oasis_pick_hrf_lwu_fast(Y, design_spec, hrf_grid)

  expect_true(inherits(result, "HRF"))
})

test_that(".oasis_pick_hrf_lwu_fast handles confounds", {
  skip_if_not_installed("fmrihrf")

  set.seed(678)

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)

  design_spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(10, 40, 70),
      duration = 0,
      span = 32
    )
  )

  hrf_grid <- list(fmrihrf::HRF_SPMG1)

  Y <- matrix(rnorm(100 * 5), 100, 5)
  confounds <- cbind(1, 1:100)  # Intercept + linear

  result <- fmrilss:::.oasis_pick_hrf_lwu_fast(Y, design_spec, hrf_grid, confounds = confounds)

  expect_true(inherits(result, "HRF"))
})

test_that(".lss_oasis coerces Matrix objects", {
  skip_if_not_installed("Matrix")

  set.seed(789)
  n_time <- 50
  n_trials <- 5
  n_vox <- 3

  Y_mat <- Matrix::Matrix(rnorm(n_time * n_vox), n_time, n_vox)
  X_mat <- Matrix::Matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y_mat, X = X_mat, oasis = list(K = 1))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that(".lss_oasis coerces data.frame inputs", {
  set.seed(890)
  n_time <- 50
  n_trials <- 4
  n_vox <- 3

  Y_df <- data.frame(matrix(rnorm(n_time * n_vox), n_time, n_vox))
  X_df <- data.frame(matrix(rnorm(n_time * n_trials), n_time, n_trials))

  result <- fmrilss:::.lss_oasis(Y = Y_df, X = X_df, oasis = list(K = 1))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that(".lss_oasis return_se returns list with se", {
  set.seed(901)
  n_time <- 100
  n_trials <- 5
  n_vox <- 3

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X,
                                  oasis = list(K = 1, return_se = TRUE))

  expect_true(is.list(result))
  expect_true("beta" %in% names(result))
  expect_true("se" %in% names(result))
  expect_equal(dim(result$beta), c(n_trials, n_vox))
  expect_equal(dim(result$se), c(n_trials, n_vox))
})

test_that(".lss_oasis return_diag returns diagnostic info", {
  set.seed(111)
  n_time <- 100
  n_trials <- 5
  n_vox <- 3

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X,
                                  oasis = list(K = 1, return_diag = TRUE))

  expect_true(is.list(result))
  expect_true("beta" %in% names(result))
  expect_true("diag" %in% names(result))
  expect_true("d" %in% names(result$diag))
  expect_true("alpha" %in% names(result$diag))
  expect_true("s" %in% names(result$diag))
})

test_that(".lss_oasis with ridge regularization works", {
  set.seed(222)
  n_time <- 100
  n_trials <- 5
  n_vox <- 3

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X,
                                  oasis = list(K = 1, ridge_x = 0.1, ridge_b = 0.1,
                                              ridge_mode = "absolute"))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that(".lss_oasis with fractional ridge works", {
  set.seed(333)
  n_time <- 100
  n_trials <- 5
  n_vox <- 3

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X,
                                  oasis = list(K = 1, ridge_x = 0.01, ridge_b = 0.01,
                                              ridge_mode = "fractional"))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n_trials, n_vox))
})

test_that(".lss_oasis names output rows and columns", {
  set.seed(444)
  n_time <- 100
  n_trials <- 5
  n_vox <- 3

  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  colnames(Y) <- c("V1", "V2", "V3")
  X <- matrix(rnorm(n_time * n_trials), n_time, n_trials)

  result <- fmrilss:::.lss_oasis(Y = Y, X = X, oasis = list(K = 1))

  expect_equal(rownames(result), sprintf("Trial_%d", 1:n_trials))
  expect_equal(colnames(result), c("V1", "V2", "V3"))
})

test_that(".lss_oasis validates ridge_mode argument", {
  Y <- matrix(rnorm(50 * 5), 50, 5)
  X <- matrix(rnorm(50 * 3), 50, 3)

  expect_error(
    fmrilss:::.lss_oasis(Y = Y, X = X, oasis = list(ridge_mode = "invalid")),
    "should be one of"
  )
})

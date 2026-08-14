contract_tol <- function(reference) {
  1e-10 + 1e-8 * max(abs(reference), na.rm = TRUE)
}

test_that("project_confounds uses only the rank-revealed nuisance span", {
  set.seed(9187)
  n <- 80
  nuisance <- cbind(a = rnorm(n), b = rnorm(n))
  redundant <- cbind(nuisance, duplicate_a = 2 * nuisance[, "a"])
  Q <- project_confounds(redundant)
  expected_basis <- qr.Q(qr(redundant))[, seq_len(qr(redundant)$rank), drop = FALSE]
  expected <- diag(n) - tcrossprod(expected_basis)
  expect_lte(max(abs(Q - expected)), contract_tol(expected))
  expect_lte(max(abs(Q %*% redundant)), contract_tol(redundant))
  expect_equal(sum(diag(Q)), n - qr(redundant)$rank, tolerance = 1e-8)

  X <- matrix(rnorm(n * 5), n, 5)
  Y <- matrix(rnorm(n * 3), n, 3)
  Z <- cbind(Intercept = 1, Trend = seq_len(n))
  full_fit <- lss(Y, X, Z = Z, Nuisance = redundant)
  reduced_fit <- lss(Y, X, Z = Z, Nuisance = nuisance)
  expect_lte(max(abs(full_fit - reduced_fit)), contract_tol(reduced_fit))

  cached_fit <- lss(
    Q %*% Y, Q %*% X, Z = Q %*% Z, method = "cpp_optimized"
  )
  expect_lte(max(abs(full_fit - cached_fit)), contract_tol(full_fit))
})

test_that("raw multi-basis OASIS uses an explicit, permutation-safe identity map", {
  set.seed(9188)
  T <- 90
  K <- 3
  ntrials <- 4
  X <- matrix(rnorm(T * K * ntrials), T)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  Y <- matrix(rnorm(T * 2), T, 2)
  map <- data.frame(
    column = colnames(X), trial = rep(seq_len(ntrials), each = K),
    basis = rep(seq_len(K), times = ntrials)
  )
  opts <- list(
    K = K, ntrials = ntrials, trial_basis_map = map,
    ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
  )
  reference <- lss(Y, X, method = "oasis", oasis = opts)

  basis_major <- unlist(lapply(seq_len(K), function(k) seq.int(k, ncol(X), by = K)))
  X_permuted <- X[, basis_major, drop = FALSE]
  permuted_map <- map[match(colnames(X_permuted), map$column), , drop = FALSE]
  got <- lss(
    Y, X_permuted, method = "oasis",
    oasis = utils::modifyList(opts, list(trial_basis_map = permuted_map))
  )
  expect_lte(max(abs(got - reference)), contract_tol(reference))
  expect_identical(rownames(got), rownames(reference))
  expect_identical(rownames(got)[1], "Trial_1:Basis_1")
  expect_identical(
    attr(got, "trial_basis_map")[, c("trial", "basis")],
    attr(reference, "trial_basis_map")[, c("trial", "basis")]
  )

  expect_error(
    lss(Y, X, method = "oasis", oasis = opts[names(opts) != "trial_basis_map"]),
    "trial_basis_map"
  )
  expect_error(
    lss(Y, unname(X), method = "oasis", oasis = opts),
    "unique, non-empty column names"
  )
  fractional_map <- map
  fractional_map$trial <- fractional_map$trial + 0.9
  expect_error(
    lss(
      Y, X, method = "oasis",
      oasis = utils::modifyList(opts, list(trial_basis_map = fractional_map))
    ),
    "exact integers"
  )
})

test_that("fmridesign mapping is invariant to upstream column permutation", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)
  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg3", nbasis = 3),
    data = trials, block = ~run, sampling_frame = sframe
  )
  dm <- fmridesign::design_matrix(emod)
  reference <- fmrilss:::.prepare_fmridesign_lss(dm, emod)
  permutation <- rev(seq_len(ncol(dm)))
  dm_permuted <- as.matrix(dm)[, permutation, drop = FALSE]
  metadata <- attr(dm, "col_metadata")
  attr(dm_permuted, "col_metadata") <- metadata[permutation, , drop = FALSE]
  got <- fmrilss:::.prepare_fmridesign_lss(dm_permuted, emod)
  expect_lte(max(abs(got$X - reference$X)), contract_tol(reference$X))
  expect_identical(got$map[, c("trial", "trial_name", "basis")],
                   reference$map[, c("trial", "trial_name", "basis")])
  renamed <- as.matrix(dm)
  colnames(renamed) <- paste0("opaque_", seq_len(ncol(renamed)))
  attr(renamed, "col_metadata") <- metadata
  expect_error(
    fmrilss:::.prepare_fmridesign_lss(renamed, emod),
    "lack basis identity"
  )
  forged <- as.matrix(dm)
  forged_names <- colnames(forged)
  target <- metadata$term_tag == names(emod$terms)[1L]
  forged_names[target] <- sub(
    paste0("^", names(emod$terms)[1L], "_[^_]+"),
    paste0(names(emod$terms)[1L], "_forged"), forged_names[target]
  )
  colnames(forged) <- forged_names
  attr(forged, "col_metadata") <- metadata
  expect_error(
    fmrilss:::.prepare_fmridesign_lss(forged, emod),
    "condition metadata"
  )
})

test_that("rank-deficient OASIS inference fails closed", {
  set.seed(9189)
  T <- 60
  K <- 2
  x <- matrix(rnorm(T * K), T, K)
  X <- cbind(x, x)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  map <- data.frame(
    column = colnames(X), trial = rep(1:2, each = K),
    basis = rep(seq_len(K), times = 2)
  )
  expect_error(
    lss(
      matrix(rnorm(T * 2), T, 2), X, method = "oasis",
      oasis = list(
        K = K, ntrials = 2, trial_basis_map = map, return_se = TRUE,
        ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
      )
    ),
    "rank-deficient trial model"
  )
  x1 <- matrix(rnorm(T), T, 1)
  expect_error(
    lss(
      matrix(rnorm(T * 2), T, 2), cbind(x1, x1), method = "oasis",
      oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)
    ),
    "Unpenalized OASIS is undefined"
  )
})

test_that("design-spec OASIS rejects K and trial counts that contradict events", {
  skip_if_not_installed("fmrihrf")
  T <- 90
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = c(10, 30, 50, 70), hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  Y <- matrix(rnorm(T * 2), T, 2)
  for (wrong_K in c(2, 4)) {
    expect_error(
      lss(
        Y, X = NULL, method = "oasis",
        oasis = list(
          design_spec = spec, K = wrong_K, ntrials = 4,
          ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
        )
      ),
      "disagrees with design_spec HRF basis dimension"
    )
  }
  expect_error(
    lss(
      Y, X = NULL, method = "oasis",
      oasis = list(design_spec = spec, K = 3, ntrials = 3)
    ),
    "disagrees with design_spec event count"
  )
})

wilson_interval <- function(successes, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

test_that("unpenalized OASIS standard errors calibrate under iid Gaussian noise", {
  set.seed(9190)
  T <- 80
  reps <- 4000
  Z <- matrix(1, T, 1)
  for (K in c(1L, 3L)) {
    ntrials <- 4L
    X <- matrix(rnorm(T * ntrials * K), T)
    colnames(X) <- paste0("x", seq_len(ncol(X)))
    opts <- list(
      return_se = TRUE, ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
    )
    if (K > 1L) {
      opts$K <- K
      opts$ntrials <- ntrials
      opts$trial_basis_map <- data.frame(
        column = colnames(X), trial = rep(seq_len(ntrials), each = K),
        basis = rep(seq_len(K), times = ntrials)
      )
    }
    fit <- lss(matrix(rnorm(T * reps), T, reps), X, Z = Z,
               method = "oasis", oasis = opts)
    model_df <- T - 1L - 2L * K
    t_value <- fit$beta[1, ] / fit$se[1, ]
    covered <- abs(t_value) <= stats::qt(0.975, df = model_df)
    coverage <- mean(covered)
    interval <- wilson_interval(sum(covered), length(covered))
    expect_gte(coverage, 0.935)
    expect_lte(coverage, 0.965)
    expect_lte(interval[1], 0.95)
    expect_gte(interval[2], 0.95)
    type1 <- 1 - coverage
    type1_interval <- 1 - rev(interval)
    expect_gte(type1, 0.035)
    expect_lte(type1, 0.065)
    expect_lte(type1_interval[1], 0.05)
    expect_gte(type1_interval[2], 0.05)
  }
})

test_that("estimated-prewhitening uncertainty fails closed", {
  skip_if_not_installed("fmriAR")
  X <- matrix(rnorm(60 * 4), 60, 4)
  Y <- matrix(rnorm(60 * 2), 60, 2)
  expect_error(
    lss(
      Y, X, method = "oasis",
      oasis = list(
        return_se = TRUE, ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
      ),
      prewhiten = list(method = "ar", p = 1, pooling = "global")
    ),
    "not calibrated for estimated prewhitening"
  )
})

test_that("global and run whitening match explicit whitened full GLMs", {
  skip_if_not_installed("fmriAR")
  set.seed(9191)
  T <- 100
  runs <- rep(1:2, each = T / 2)
  X <- matrix(rnorm(T * 4), T, 4)
  Z <- cbind(run_1 = as.numeric(runs == 1), run_2 = as.numeric(runs == 2))
  Nuisance <- cbind(trend = rep(seq(-1, 1, length.out = T / 2), 2))
  Y <- matrix(rnorm(T * 3), T, 3)

  for (pooling in c("global", "run")) {
    whitening <- list(method = "ar", p = 1, pooling = pooling, runs = runs)
    got <- lss(
      Y, X, Z = Z, Nuisance = Nuisance, method = "oasis",
      oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
      prewhiten = whitening
    )
    pw <- fmrilss:::.prewhiten_data(Y, X, Z, Nuisance, whitening)
    reference <- matrix(NA_real_, ncol(X), ncol(Y))
    for (j in seq_len(ncol(X))) {
      design <- cbind(
        pw$Z_whitened, pw$Nuisance_whitened,
        Xj = pw$X_whitened[, j],
        Xother = rowSums(pw$X_whitened) - pw$X_whitened[, j]
      )
      reference[j, ] <- lm.fit(design, pw$Y_whitened)$coefficients[ncol(design) - 1L, ]
    }
    expect_lte(max(abs(got - reference)), contract_tol(reference))
  }
})

test_that("run whitening does not filter across run boundaries", {
  skip_if_not_installed("fmriAR")
  set.seed(9192)
  T <- 80
  runs <- rep(1:2, each = T / 2)
  Y <- matrix(rnorm(T * 3), T, 3)
  X <- matrix(0, T, 1)
  X[T / 2, 1] <- 1
  whitening <- list(
    method = "ar", p = 1, pooling = "run", runs = runs,
    compute_residuals = FALSE
  )
  pw <- fmrilss:::.prewhiten_data(Y, X, prewhiten = whitening)
  expect_equal(max(abs(pw$X_whitened[runs == 2, 1])), 0)
})

test_that("voxel-ridge OASIS solves the final LSS model in whitened space", {
  skip_if_not_installed("fmriAR")
  skip_if_not_installed("fmrihrf")
  set.seed(9194)
  T <- 140
  V <- 4
  K <- 3
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = seq(10, 110, by = 20), hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  built <- fmrilss:::.oasis_build_X_from_events(spec)
  X <- built$X_trials
  Z <- matrix(1, T, 1)
  innovations <- matrix(rnorm(T * V), T, V)
  Y <- apply(innovations, 2, stats::filter, filter = 0.7, method = "recursive")
  Y <- as.matrix(Y)
  whitening <- list(method = "ar", p = 1, pooling = "global")

  got <- lss(
    Y, X = NULL, method = "oasis",
    oasis = list(
      design_spec = spec, hrf_mode = "voxel_ridge", K = K,
      ntrials = length(spec$cond$onsets), lambda_shape = 0
    ),
    prewhiten = whitening
  )

  pw <- fmrilss:::.prewhiten_data(Y, X, Z, NULL, whitening)
  qr_nuis <- qr(pw$Z_whitened)
  nuisance_basis <- qr.Q(qr_nuis)[, seq_len(qr_nuis$rank), drop = FALSE]
  vhrf <- fmrilss:::.estimate_voxel_hrf_fast(
    pw$Y_whitened, pw$X_whitened, spec,
    N_nuis = nuisance_basis, K = K, lambda_shape = 0
  )
  basis_convolved <- lapply(seq_len(K), function(k) {
    qr.resid(qr_nuis, pw$X_whitened[, seq.int(k, ncol(X), by = K), drop = FALSE])
  })
  Y_res <- qr.resid(qr_nuis, pw$Y_whitened)
  coeff <- vhrf$coefficients * vhrf$ref_norm
  reference <- matrix(NA_real_, length(spec$cond$onsets), V)
  for (v in seq_len(V)) {
    Xv <- Reduce(`+`, lapply(seq_len(K), function(k) basis_convolved[[k]] * coeff[k, v]))
    for (j in seq_len(ncol(Xv))) {
      design <- cbind(
        nuisance_basis, Xj = Xv[, j],
        Xother = rowSums(Xv) - Xv[, j]
      )
      reference[j, v] <- lm.fit(design, Y_res[, v])$coefficients[ncol(design) - 1L]
    }
  }
  expect_lte(max(abs(got - reference)), contract_tol(reference))
})

test_that("voxel-ridge OASIS honors an explicit no-intercept model", {
  skip_if_not_installed("fmrihrf")
  set.seed(9196)
  T <- 120
  V <- 3
  K <- 3
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = seq(10, 90, by = 20), hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  X <- fmrilss:::.oasis_build_X_from_events(spec)$X_trials
  Y <- matrix(4 + rnorm(T * V), T, V)
  got <- lss(
    Y, X = NULL, method = "oasis",
    oasis = list(
      design_spec = spec, hrf_mode = "voxel_ridge", K = K,
      ntrials = length(spec$cond$onsets), lambda_shape = 0,
      add_intercept = FALSE
    )
  )

  vhrf <- fmrilss:::.estimate_voxel_hrf_fast(Y, X, spec, K = K, lambda_shape = 0)
  coeff <- vhrf$coefficients * vhrf$ref_norm
  basis_convolved <- lapply(seq_len(K), function(k) {
    X[, seq.int(k, ncol(X), by = K), drop = FALSE]
  })
  reference <- matrix(NA_real_, length(spec$cond$onsets), V)
  for (v in seq_len(V)) {
    Xv <- Reduce(`+`, lapply(seq_len(K), function(k) basis_convolved[[k]] * coeff[k, v]))
    for (j in seq_len(ncol(Xv))) {
      design <- cbind(Xj = Xv[, j], Xother = rowSums(Xv) - Xv[, j])
      reference[j, v] <- lm.fit(design, Y[, v])$coefficients[1L]
    }
  }
  expect_lte(max(abs(got - reference)), contract_tol(reference))
})

test_that("voxel-ridge global shrinkage restores unit-energy shapes", {
  skip_if_not_installed("fmrihrf")
  set.seed(9195)
  T <- 120
  V <- 12
  K <- 3
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = seq(10, 90, by = 20), hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  X <- fmrilss:::.oasis_build_X_from_events(spec)$X_trials
  fit <- fmrilss:::.estimate_voxel_hrf_fast(
    matrix(rnorm(T * V), T, V), X, spec,
    K = K, shrink_global = 0.7
  )
  times <- fmrihrf::samples(sframe, global = TRUE)
  B_time <- fmrihrf::evaluate(
    fmrihrf::regressor(onsets = 0, hrf = spec$cond$hrf, duration = 0, span = 30),
    grid = times, precision = 0.1, method = "conv"
  )
  if (inherits(B_time, "Matrix")) B_time <- as.matrix(B_time)
  energy <- sqrt(colSums((B_time %*% fit$coefficients)^2))
  expect_lte(max(abs(energy - 1)), 1e-8)
})

test_that("voxel-ridge handles zero-energy shapes explicitly and quietly", {
  skip_if_not_installed("fmrihrf")
  T <- 100
  V <- 2
  K <- 3
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = seq(10, 70, by = 20), hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  X <- fmrilss:::.oasis_build_X_from_events(spec)$X_trials
  zero_fit <- fmrilss:::.estimate_voxel_hrf_fast(matrix(0, T, V), X, spec, K = K)
  near_fit <- fmrilss:::.estimate_voxel_hrf_fast(matrix(1e-20, T, V), X, spec, K = K)
  times <- fmrihrf::samples(sframe, global = TRUE)
  B_time <- fmrihrf::evaluate(
    fmrihrf::regressor(onsets = 0, hrf = spec$cond$hrf, duration = 0, span = 30),
    grid = times, precision = 0.1, method = "conv"
  )
  if (inherits(B_time, "Matrix")) B_time <- as.matrix(B_time)
  expect_lte(max(abs(sqrt(colSums((B_time %*% zero_fit$coefficients)^2)) - 1)), 1e-8)
  expect_lte(max(abs(sqrt(colSums((B_time %*% near_fit$coefficients)^2)) - 1)), 1e-8)
  expect_identical(zero_fit$fallback_voxels, seq_len(V))
  expect_warning(
    got <- lss(
      matrix(0, T, V), X = NULL, method = "oasis",
      oasis = list(
        design_spec = spec, hrf_mode = "voxel_ridge", K = K,
        ntrials = length(spec$cond$onsets)
      )
    ),
    NA
  )
  expect_true(all(is.finite(got)))
  expect_equal(unname(got), matrix(0, nrow(got), ncol(got)), tolerance = 0)
})

test_that("voxel-ridge shape controls reject values outside their domains", {
  skip_if_not_installed("fmrihrf")
  T <- 60
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = c(10, 30), hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  X <- fmrilss:::.oasis_build_X_from_events(spec)$X_trials
  Y <- matrix(rnorm(T), T, 1)
  expect_error(
    fmrilss:::.estimate_voxel_hrf_fast(Y, X, spec, K = 3, shrink_global = 1.5),
    "between 0 and 1"
  )
  expect_error(
    fmrilss:::.estimate_voxel_hrf_fast(Y, X, spec, K = 3, shrink_global = -0.1),
    "non-negative"
  )
  expect_error(
    fmrilss:::.estimate_voxel_hrf_fast(Y, X, spec, K = 3, lambda_shape = Inf),
    "non-negative"
  )
  expect_error(
    fmrilss:::.estimate_voxel_hrf_fast(Y, X, spec, K = 3, mu_rough = -1),
    "non-negative"
  )
  expect_error(
    lss(
      Y, X = NULL, method = "oasis",
      oasis = list(
        design_spec = spec, hrf_mode = "voxel_ridge", K = 3,
        ntrials = length(spec$cond$onsets), return_se = TRUE
      )
    ),
    "uncertainty is not calibrated"
  )
  expect_error(
    lss(
      Y, X = NULL, method = "oasis",
      oasis = list(
        design_spec = spec, hrf_mode = "voxel_ridge", K = 3,
        ntrials = length(spec$cond$onsets), return_diag = TRUE
      )
    ),
    "diagnostics are not defined"
  )
})

test_that("voxel-HRF designs preserve run boundaries and event amplitudes", {
  skip_if_not_installed("fmrihrf")
  sframe <- fmrihrf::sampling_frame(blocklens = c(20, 20), TR = 1)
  events <- data.frame(
    onset = c(18, 2), duration = c(0, 1), condition = "A", run = c(1, 2)
  )
  built <- fmrilss:::.voxhrf_trial_basis(events, fmrihrf::HRF_SPMG1, sframe)
  expect_equal(max(abs(built$X[21:40, 1])), 0)
  expect_equal(max(abs(built$X[1:20, 2])), 0)

  scaled <- events
  scaled$amplitude <- c(3, 1)
  built_scaled <- fmrilss:::.voxhrf_trial_basis(scaled, fmrihrf::HRF_SPMG1, sframe)
  expect_lte(max(abs(built_scaled$X[, 1] - 3 * built$X[, 1])),
             contract_tol(built$X[, 1]))
  expect_error(
    fmrilss:::.voxhrf_trial_basis(events[, names(events) != "run"],
                                  fmrihrf::HRF_SPMG1, sframe),
    "explicit run column"
  )
  malformed <- events
  malformed$run <- c(1.9, 2.1)
  expect_error(
    fmrilss:::.voxhrf_trial_basis(malformed, fmrihrf::HRF_SPMG1, sframe),
    "exact integer run indices"
  )
})

test_that("LSSBeta matrix extraction preserves units and sampling metadata", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")
  T <- 60
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  events <- data.frame(onset = c(5, 25, 45), duration = 0, condition = "A")
  Y <- matrix(rnorm(T * 2), T, 2)
  est <- estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG1, sframe = sframe)
  r_fit <- lss_with_hrf(Y, events, est, engine = "R", verbose = FALSE)
  cpp_fit <- as.matrix(lss_with_hrf(Y, events, est, engine = "C++", verbose = FALSE))
  expect_identical(attr(r_fit, "normalization"), attr(cpp_fit, "normalization"))
  expect_identical(attr(r_fit, "units"), "peak-response amplitude")
  expect_identical(attr(cpp_fit, "sframe"), sframe)
})

test_that("unnormalized VoxelHRF fits do not claim peak-response units", {
  skip_if_not_installed("fmrihrf")
  T <- 50
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  events <- data.frame(onset = c(5, 25), duration = 0, condition = "A")
  Y <- matrix(rnorm(T * 2), T, 2)
  est <- estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG1, sframe = sframe)
  est$normalization <- NULL
  got <- lss_with_hrf(Y, events, est, engine = "R", verbose = FALSE)
  expect_identical(attr(got, "normalization"), "unspecified")
  expect_identical(
    attr(got, "units"),
    "coefficient on supplied event design relative to supplied HRF shape"
  )
})

test_that("voxel-HRF estimates bind shapes to voxel identity", {
  skip_if_not_installed("fmrihrf")
  T <- 70
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  events <- data.frame(onset = c(5, 25, 45), duration = 0, condition = "A")
  Y <- matrix(rnorm(T * 2), T, 2, dimnames = list(NULL, c("vA", "vB")))
  est <- estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG1, sframe = sframe)
  expect_identical(colnames(est$coefficients), colnames(Y))
  expect_identical(names(est$amplitude_scale), colnames(Y))
  reference <- lss_with_hrf(Y, events, est, engine = "R", verbose = FALSE)

  permuted <- est
  permuted$coefficients <- permuted$coefficients[, 2:1, drop = FALSE]
  got <- lss_with_hrf(Y, events, permuted, engine = "R", verbose = FALSE)
  expect_lte(max(abs(got - reference)), contract_tol(reference))

  duplicated <- est
  colnames(duplicated$coefficients) <- c("vA", "vA")
  expect_error(lss_with_hrf(Y, events, duplicated), "same unique voxel names")
  unnamed <- est
  colnames(unnamed$coefficients) <- NULL
  expect_error(lss_with_hrf(Y, events, unnamed), "both have voxel names")
})

test_that("internal voxel-HRF Armadillo helper accepts explicit multi-run identity", {
  skip_if_not_installed("fmrihrf")
  T <- 80
  sframe <- fmrihrf::sampling_frame(blocklens = c(40, 40), TR = 1)
  events <- data.frame(
    onset = c(5, 25, 5, 25), duration = 0, condition = "A",
    run = c(1, 1, 2, 2)
  )
  Y <- matrix(rnorm(T * 2), T, 2)
  est <- estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG1, sframe = sframe)
  got <- fmrilss:::.voxhrf_betas_cpp_arma(
    Y, events$onset, events$duration, est,
    sframe = sframe, runs = events$run
  )
  expect_equal(dim(got), c(nrow(events), ncol(Y)))
  expect_true(all(is.finite(got)))
})

test_that("internal voxel-HRF Armadillo helper binds shapes to named voxels", {
  skip_if_not_installed("fmrihrf")
  T <- 70
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  events <- data.frame(onset = c(5, 25, 45), duration = 0, condition = "A")
  Y <- matrix(rnorm(T * 2), T, 2, dimnames = list(NULL, c("vA", "vB")))
  est <- estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG3, sframe = sframe)
  reference <- fmrilss:::.voxhrf_betas_cpp_arma(
    Y, events$onset, events$duration, est, sframe = sframe
  )

  permuted <- est
  permuted$coefficients <- permuted$coefficients[, 2:1, drop = FALSE]
  got <- fmrilss:::.voxhrf_betas_cpp_arma(
    Y, events$onset, events$duration, permuted, sframe = sframe
  )
  expect_lte(max(abs(got - reference)), contract_tol(reference))

  duplicated <- est
  colnames(duplicated$coefficients) <- c("vA", "vA")
  expect_error(
    fmrilss:::.voxhrf_betas_cpp_arma(
      Y, events$onset, events$duration, duplicated, sframe = sframe
    ),
    "same unique voxel names"
  )
  unnamed <- est
  colnames(unnamed$coefficients) <- NULL
  expect_error(
    fmrilss:::.voxhrf_betas_cpp_arma(
      Y, events$onset, events$duration, unnamed, sframe = sframe
    ),
    "both have voxel names"
  )
})

test_that("lss_with_hrf validates chunk_size as a positive integer", {
  skip_if_not_installed("fmrihrf")
  T <- 30
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  events <- data.frame(onset = 5, duration = 0, condition = "A")
  Y <- matrix(rnorm(T), T, 1)
  est <- estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG1, sframe = sframe)
  expect_error(lss_with_hrf(Y, events, est, chunk_size = 1.5), "positive integer")
  expect_error(lss_with_hrf(Y, events, est, chunk_size = Inf), "positive integer")
})

test_that("multi-basis OASIS uncertainty is estimator and voxel-set invariant", {
  skip_if_not_installed("fmrihrf")
  set.seed(9101)
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(10, 30, 50, 70),
      hrf = fmrihrf::HRF_SPMG3,
      span = 30
    )
  )
  y <- rnorm(100)
  fit <- function(V, return_se) {
    lss(
      matrix(rep(y, V), 100, V), X = NULL, method = "oasis",
      oasis = list(
        design_spec = spec,
        return_se = return_se,
        ridge_mode = "absolute",
        ridge_x = 0,
        ridge_b = 0
      )
    )
  }
  with_se_8 <- fit(8, TRUE)
  with_se_40 <- fit(40, TRUE)
  plain <- fit(8, FALSE)

  expect_lte(max(abs(with_se_8$beta - plain)), contract_tol(plain))
  expect_lte(
    max(abs(with_se_8$se[, 1] - with_se_40$se[, 1])),
    contract_tol(with_se_8$se[, 1])
  )
  expect_identical(dimnames(with_se_8$beta), dimnames(with_se_8$se))
})

test_that("unpenalized multi-basis OASIS SE agrees with a direct GLM", {
  skip_if_not_installed("fmrihrf")
  set.seed(9102)
  T <- 120
  V <- 2
  K <- 3
  onsets <- c(10, 35, 65, 90)
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30)
  )
  built <- fmrilss:::.oasis_build_X_from_events(spec)
  X <- built$X_trials
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  identity_map <- data.frame(
    column = colnames(X), trial = rep(seq_along(onsets), each = K),
    basis = rep(seq_len(K), times = length(onsets))
  )
  Z <- cbind(Intercept = 1, trend = seq(-1, 1, length.out = T))
  Y <- matrix(rnorm(T * V), T, V)
  got <- lss(
    Y, X, Z = Z, method = "oasis",
    oasis = list(
      K = K, ntrials = length(onsets), trial_basis_map = identity_map,
      return_se = TRUE,
      ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
    )
  )

  ref_beta <- matrix(NA_real_, ncol(X), V)
  ref_se <- matrix(NA_real_, ncol(X), V)
  for (j in seq_along(onsets)) {
    idx <- ((j - 1L) * K + 1L):(j * K)
    Aj <- X[, idx, drop = FALSE]
    # Preserve the complete K-dimensional other-trial span.
    Bj <- Reduce(`+`, lapply(setdiff(seq_along(onsets), j), function(h) {
      X[, ((h - 1L) * K + 1L):(h * K), drop = FALSE]
    }))
    D <- cbind(Z, Aj, Bj)
    fit <- lm.fit(D, Y)
    df <- nrow(D) - fit$rank
    inv <- solve(crossprod(D))
    pos <- ncol(Z) + seq_len(K)
    ref_beta[idx, ] <- fit$coefficients[pos, , drop = FALSE]
    residual_var <- colSums(fit$residuals^2) / df
    ref_se[idx, ] <- sqrt(outer(diag(inv)[pos], residual_var))
  }
  expect_lte(max(abs(got$beta - ref_beta)), contract_tol(ref_beta))
  expect_lte(max(abs(got$se - ref_se)), contract_tol(ref_se))
})

test_that("ridge uncertainty fails closed", {
  set.seed(9103)
  X <- matrix(rnorm(80 * 5), 80, 5)
  Y <- matrix(rnorm(80 * 3), 80, 3)
  expect_error(
    lss(Y, X, method = "oasis", oasis = list(return_se = TRUE)),
    "only for unpenalized OASIS"
  )
})

test_that("multi-basis other conditions retain their complete span", {
  skip_if_not_installed("fmrihrf")
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)
  built <- fmrilss:::.oasis_build_X_from_events(list(
    sframe = sframe,
    cond = list(onsets = c(10, 30), hrf = fmrihrf::HRF_SPMG3),
    others = list(list(onsets = c(20, 50), hrf = fmrihrf::HRF_SPMG3))
  ))
  expect_equal(ncol(built$X_other), 3L)
  expect_equal(qr(built$X_other)$rank, 3L)
})

test_that("SBHM amplitude adapters preserve event fields and nuisance spans", {
  skip_if_not_installed("fmrihrf")
  set.seed(9193)
  T <- 120
  r <- 3
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = T) / 4),
    exp(-seq(0, 30, length.out = T) / 7),
    exp(-seq(0, 30, length.out = T) / 10)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(10, 45, 80), duration = c(0, 3, 1),
      amplitude = c(2, -4, 1), span = 30
    ),
    others = list(list(onsets = c(25, 65), duration = c(2, 0), span = 30))
  )
  regs <- fmrilss:::.sbhm_build_trial_regs(sbhm, spec)
  expect_length(regs, 3L)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  direct_second <- fmrihrf::evaluate(
    fmrihrf::regressor(
      onsets = spec$cond$onsets[2], hrf = hrf_B,
      duration = spec$cond$duration[2], amplitude = spec$cond$amplitude[2],
      span = 30, summate = FALSE
    ),
    grid = sbhm$tgrid, precision = 0.1, method = "conv"
  )
  expect_lte(max(abs(regs[[2]] - direct_second)), contract_tol(direct_second))

  alpha <- matrix(c(0.7, -0.2, 0.4), r, 1)
  X_main <- do.call(cbind, lapply(regs, function(Xt) as.numeric(Xt %*% alpha)))
  X_other <- fmrilss:::.sbhm_build_design(sbhm, spec)$X_other
  beta_truth <- c(1.5, -0.5, 2)
  gamma <- seq(0.4, 1.2, length.out = ncol(X_other))
  Y <- matrix(3 + X_main %*% beta_truth + X_other %*% gamma, T, 1)
  got <- fmrilss:::sbhm_amplitude_ls(
    Y, sbhm, spec, alpha, ridge = 0, return_se = FALSE
  )
  expect_lte(max(abs(got[, 1] - beta_truth)), contract_tol(beta_truth))

  out <- lss_sbhm(
    Y, sbhm, spec,
    match = list(topK = 1, soft_blend = FALSE, alpha_source = "prepass"),
    amplitude = list(method = "global_ls", ridge = 0),
    return = "amplitude"
  )
  X_diag <- do.call(cbind, lapply(
    regs, function(Xt) as.numeric(Xt %*% out$alpha_coords[, 1])
  ))
  X_diag <- fmrilss:::.sbhm_resid(
    X_diag, cbind(fmrilss:::.sbhm_build_design(sbhm, spec)$intercepts, X_other)
  )
  G <- crossprod(X_diag)
  scale <- sqrt(diag(G))
  cosine <- G / outer(scale, scale)
  diag(cosine) <- 0
  rho_reference <- max(abs(cosine[upper.tri(cosine)]))
  expect_lte(abs(out$diag$rho_max[1] - rho_reference), contract_tol(rho_reference))

  gated <- lss_sbhm(
    Y, sbhm, spec,
    match = list(topK = 1, soft_blend = FALSE, alpha_source = "prepass"),
    amplitude = list(
      method = "global_ls", ridge = 0,
      cond_gate = list(metric = "rho", thr = rho_reference - 1e-8, fallback = "lss1")
    ),
    return = "amplitude"
  )
  expect_identical(unname(gated$diag$method_used), "lss1")

  projected_with <- fmrilss:::.sbhm_alpha_from_trial_projections(Y, sbhm, spec)
  spec_without <- spec
  spec_without$others <- NULL
  projected_without <- fmrilss:::.sbhm_alpha_from_trial_projections(Y, sbhm, spec_without)
  expect_gt(max(abs(projected_with$alpha - projected_without$alpha)), 1e-6)
})

test_that("fmridesign mapping is trial-major and validation invariant", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  set.seed(9104)
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)
  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg3", nbasis = 3),
    data = trials, block = ~run, sampling_frame = sframe
  )
  Y <- matrix(rnorm(100 * 3), 100, 3)
  opts <- list(K = 3, ntrials = nrow(trials), ridge_mode = "absolute",
               ridge_x = 0, ridge_b = 0)
  expect_warning(
  checked <- suppressMessages(lss_design(
      Y, emod, oasis = opts, validate = TRUE
    )),
    "collinearity.*ridge",
    ignore.case = TRUE
  )
  unchecked <- lss_design(Y, emod, oasis = opts, validate = FALSE)
  expect_lte(max(abs(checked - unchecked)), contract_tol(checked))

  X_basis_major <- as.matrix(fmridesign::design_matrix(emod))
  trial_levels <- as.character(emod$terms$trial$condition_levels)
  expected_names <- unlist(lapply(
    paste0("trial_", trial_levels),
    function(trial_name) paste0(trial_name, sprintf("_b%02d", 1:3))
  ), use.names = FALSE)
  idx <- match(expected_names, colnames(X_basis_major))
  expect_false(anyNA(idx))
  X_direct <- X_basis_major[, idx, drop = FALSE]
  opts$trial_basis_map <- data.frame(
    column = colnames(X_direct), trial = rep(seq_len(nrow(trials)), each = 3),
    basis = rep(seq_len(3), times = nrow(trials))
  )
  direct <- lss(Y, X_direct, method = "oasis", oasis = opts)
  expect_lte(max(abs(checked - direct)), contract_tol(direct))
  expect_match(rownames(checked)[1], "basis_1$")

  inferred <- lss(
    Y,
    fmridesign::design_matrix(emod),
    method = "oasis",
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)
  )
  expect_lte(max(abs(checked - inferred)), contract_tol(inferred))
  expect_identical(rownames(inferred), rownames(checked))
  expect_identical(
    attr(inferred, "trial_basis_map")[, c("trial", "basis", "output_name")],
    attr(checked, "trial_basis_map")[, c("trial", "basis", "output_name")]
  )
})

test_that("fmridesign non-trial event terms are fixed rather than targets", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  sframe <- fmrihrf::sampling_frame(blocklens = 120, TR = 2)
  events <- data.frame(
    onset = seq(10, 60, by = 10), run = 1,
    RT = seq(0.4, 0.9, length.out = 6)
  )
  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1") + fmridesign::hrf(RT),
    data = events, block = ~run, sampling_frame = sframe
  )
  dm <- as.matrix(fmridesign::design_matrix(emod))
  trial_names <- paste0("trial_", emod$terms$trial$condition_levels)
  trial_columns <- match(trial_names, colnames(dm))
  fixed_columns <- setdiff(seq_len(ncol(dm)), trial_columns)
  expect_false(anyNA(trial_columns))
  expect_length(fixed_columns, 1L)

  X <- dm[, trial_columns, drop = FALSE]
  fixed <- dm[, fixed_columns, drop = FALSE]
  set.seed(91041)
  Y <- 4 * fixed %*% matrix(c(1, -0.5), 1L, 2L) +
    matrix(rnorm(120 * 2, sd = 0.05), 120, 2)
  result <- lss_design(
    Y, emod, validate = FALSE,
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)
  )
  expect_equal(nrow(result), nrow(events))

  direct <- without_fixed <- matrix(NA_real_, nrow(events), ncol(Y))
  for (j in seq_len(nrow(events))) {
    target <- X[, j]
    other <- rowSums(X[, -j, drop = FALSE])
    direct[j, ] <- lm.fit(cbind(target, other, 1, fixed), Y)$coefficients[1L, ]
    without_fixed[j, ] <- lm.fit(cbind(target, other, 1), Y)$coefficients[1L, ]
  }
  expect_lte(max(abs(unname(result) - direct)), contract_tol(direct))
  expect_gt(max(abs(without_fixed - direct)), 1e-3)
  expect_identical(
    attr(result, "trial_basis_map")$source_column,
    trial_columns
  )
})

test_that("voxel and parcel whitening fail closed for a shared design", {
  skip_if_not_installed("fmriAR")
  Y <- matrix(rnorm(80 * 4), 80, 4)
  X <- matrix(rnorm(80 * 5), 80, 5)
  for (pooling in c("voxel", "parcel")) {
    opts <- list(method = "ar", p = 1, pooling = pooling)
    if (pooling == "parcel") opts$parcels <- c(1, 1, 2, 2)
    expect_error(
      lss(Y, X, prewhiten = opts),
      "cannot be applied to a shared design matrix"
    )
  }
})

test_that("voxel-HRF amplitude, time, nuisance, and engine contracts hold", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")
  set.seed(9105)
  for (TR in c(0.8, 1, 2)) {
    T <- 140
    V <- 2
    sframe <- fmrihrf::sampling_frame(blocklens = T, TR = TR)
    events <- data.frame(
      onset = c(0, 16.3, 38.7, 65.1, 92.4),
      duration = c(0, 1.2, 0, 2.1, 0.6),
      condition = "A"
    )
    basis <- fmrihrf::HRF_SPMG1
    built <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)
    shape <- fmrilss:::.normalize_voxel_hrf_coefficients(
      matrix(1, 1, 1), basis, span = 30
    )$coefficients[1, 1]
    truth <- cbind(rep(2, nrow(events)), rep(4, nrow(events)))
    fixed <- cbind(trend = seq(-1, 1, length.out = T))
    nuisance <- cbind(sine = sin(seq_len(T) / 11), cosine = cos(seq_len(T) / 17))
    common <- cbind(1, fixed, nuisance)
    Y <- (built$X * shape) %*% truth + common %*% matrix(rnorm(ncol(common) * V), ncol(common), V)

    estimated <- estimate_voxel_hrf(
      Y, events, basis, nuisance_regs = nuisance,
      sframe = sframe, fixed_regs = fixed
    )
    got_r <- lss_with_hrf(
      Y, events, estimated, nuisance_regs = nuisance,
      fixed_regs = fixed, engine = "R", verbose = FALSE
    )
    got_cpp <- lss_with_hrf(
      Y, events, estimated, nuisance_regs = nuisance,
      fixed_regs = fixed, engine = "C++", verbose = FALSE
    )
    expect_lte(max(abs(got_r - truth)), 5e-4)
    expect_lte(max(abs(as.matrix(got_cpp) - got_r)), contract_tol(got_r))

    X <- built$X * estimated$coefficients[1, 1]
    ref <- matrix(NA_real_, nrow(events), V)
    for (j in seq_len(nrow(events))) {
      D <- cbind(1, fixed, nuisance, X[, j], rowSums(X) - X[, j])
      ref[j, ] <- lm.fit(D, Y)$coefficients[ncol(D) - 1L, ]
    }
    expect_lte(max(abs(got_r - ref)), contract_tol(ref))
  }
})

test_that("voxel-HRF public APIs require physical-time metadata", {
  skip_if_not_installed("fmrihrf")
  Y <- matrix(rnorm(40), 20, 2)
  events <- data.frame(onset = 0, duration = 0, condition = "A")
  expect_error(
    estimate_voxel_hrf(Y, events, fmrihrf::HRF_SPMG1),
    "explicit fmrihrf sampling_frame"
  )
})

test_that("unpenalized one-basis OASIS remains scale equivariant", {
  set.seed(9106)
  T <- 80
  N <- 4
  V <- 3
  X0 <- matrix(rnorm(T * N), T, N)
  Y <- matrix(rnorm(T * V), T, V)
  Z <- cbind(intercept = 1, trend = seq(-1, 1, length.out = T))
  opts <- oasis_options(
    ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
  )

  for (design_scale in c(1e-6, 1e-5, 1, 1e5, 1e6)) {
    X <- X0 * design_scale
    got <- lss(Y, X, Z = Z, method = "oasis", oasis = opts)
    ref <- matrix(NA_real_, N, V)
    for (j in seq_len(N)) {
      D <- cbind(X[, j], rowSums(X[, -j, drop = FALSE]), Z)
      ref[j, ] <- qr.coef(qr(D), Y)[1L, ]
    }
    relative_error <- max(abs(got - ref) / pmax(abs(ref), 1e-12))
    expect_lte(relative_error, 1e-8)
  }
})

test_that("one-trial OASIS reduces to the identifiable target model", {
  set.seed(9107)
  T <- 70
  V <- 3
  Z <- cbind(intercept = 1, trend = seq(-1, 1, length.out = T))
  Y <- matrix(rnorm(T * V), T, V)

  X1 <- matrix(rnorm(T), T, 1, dimnames = list(NULL, "trial_1"))
  Qz <- qr.Q(qr(Z))
  RX1 <- X1 - Qz %*% crossprod(Qz, X1)
  RY <- Y - Qz %*% crossprod(Qz, Y)
  d1 <- drop(crossprod(RX1))
  expected_default <- crossprod(RX1, RY) / (d1 + 0.05 * d1)
  got_default <- lss(Y, X1, Z = Z, method = "oasis", oasis = oasis_options())
  expect_equal(unname(got_default), unname(expected_default), tolerance = 1e-10)

  got1 <- lss(
    Y, X1, Z = Z, method = "oasis",
    oasis = oasis_options(
      ridge_mode = "absolute", ridge_x = 0, ridge_b = 0, return_se = TRUE
    )
  )
  D1 <- cbind(X1, Z)
  coef1 <- solve(crossprod(D1), crossprod(D1, Y))
  sigma1 <- colSums((Y - D1 %*% coef1)^2) / (T - ncol(D1))
  se1 <- sqrt(solve(crossprod(D1))[1L, 1L] * sigma1)
  expect_equal(unname(got1$beta), unname(coef1[1L, , drop = FALSE]), tolerance = 1e-10)
  expect_equal(unname(got1$se), matrix(se1, 1L), tolerance = 1e-10)

  K <- 3L
  X3 <- matrix(rnorm(T * K), T, K)
  colnames(X3) <- sprintf("trial_1_basis_%d", seq_len(K))
  map3 <- data.frame(column = colnames(X3), trial = 1L, basis = seq_len(K))
  got3 <- lss(
    Y, X3, Z = Z, method = "oasis",
    oasis = oasis_options(
      K = K, ntrials = 1L, trial_basis_map = map3,
      ridge_mode = "absolute", ridge_x = 0, ridge_b = 0, return_se = TRUE
    )
  )
  D3 <- cbind(X3, Z)
  coef3 <- solve(crossprod(D3), crossprod(D3, Y))
  sigma3 <- colSums((Y - D3 %*% coef3)^2) / (T - ncol(D3))
  se3 <- sqrt(diag(solve(crossprod(D3)))[seq_len(K)] %o% sigma3)
  expect_equal(
    matrix(as.numeric(got3$beta), K, V),
    unname(coef3[seq_len(K), , drop = FALSE]),
    tolerance = 1e-10
  )
  expect_equal(
    matrix(as.numeric(got3$se), K, V),
    unname(se3),
    tolerance = 1e-10
  )
})

test_that("OASIS inference fails closed outside its calibrated domain", {
  set.seed(9108)
  X <- matrix(rnorm(3 * 2), 3, 2)
  Y <- matrix(rnorm(3 * 2), 3, 2)
  expect_error(
    lss(
      Y, X, Z = matrix(1, 3, 1), method = "oasis",
      oasis = oasis_options(
        ridge_mode = "absolute", ridge_x = 0, ridge_b = 0, return_se = TRUE
      )
    ),
    "positive residual degrees of freedom"
  )

  expect_error(
    lss(
      matrix(rnorm(40), 20, 2), X = NULL, method = "oasis",
      oasis = list(
        design_spec = list(hrf_grid = list("candidate")),
        ridge_mode = "absolute", ridge_x = 0, ridge_b = 0,
        return_se = TRUE
      )
    ),
    "not calibrated after data-adaptive HRF grid selection"
  )
})

test_that("HRF grid selection is invariant to redundant confounds", {
  skip_if_not_installed("fmrihrf")
  T <- 100
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  onsets <- c(10, 25, 40, 55, 70)
  grid <- list(fmrihrf::HRF_SPMG1, fmrihrf::HRF_GAMMA)
  times <- fmrihrf::samples(sframe, global = TRUE)
  candidate_designs <- lapply(grid, function(hrf) {
    reg <- fmrihrf::regressor(
      onsets = onsets, hrf = hrf, duration = 0,
      amplitude = 1, span = 30, summate = TRUE
    )
    evaluated <- fmrihrf::evaluate(
      reg, times, precision = 0.1, method = "conv"
    )
    if (is.matrix(evaluated)) rowSums(evaluated) else as.numeric(evaluated)
  })

  set.seed(3)
  confound <- rnorm(T)
  full_q <- qr.Q(qr(cbind(confound, 2 * confound)))
  completion_direction <- full_q[, 2L]
  Y <- matrix(100 * completion_direction + candidate_designs[[2L]], T, 1L)
  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = onsets, duration = 0, amplitude = 1, span = 30
    ),
    precision = 0.1,
    method = "conv"
  )

  reduced <- fmrilss:::.oasis_pick_hrf_lwu_fast(
    Y, spec, grid, confounds = matrix(confound, T, 1L)
  )
  redundant <- fmrilss:::.oasis_pick_hrf_lwu_fast(
    Y, spec, grid, confounds = cbind(confound, 2 * confound)
  )

  expect_identical(reduced, grid[[1L]])
  expect_identical(redundant, reduced)
})

test_that("HRF grid selection conditions on modeled other conditions", {
  skip_if_not_installed("fmrihrf")
  T <- 180L
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  target_onsets <- c(10, 40, 70, 100, 130)
  other_onsets <- c(20, 50, 80, 110, 140)
  grid <- list(fmrihrf::HRF_SPMG1, fmrihrf::HRF_GAUSSIAN)

  aggregate <- function(onsets, hrf) {
    reg <- fmrihrf::regressor(
      onsets = onsets, hrf = hrf, duration = 0,
      amplitude = 1, span = 30, summate = TRUE
    )
    evaluated <- fmrihrf::evaluate(
      reg, times, precision = 0.1, method = "conv"
    )
    if (is.matrix(evaluated)) rowSums(evaluated) else as.numeric(evaluated)
  }

  target_signal <- aggregate(target_onsets, fmrihrf::HRF_GAUSSIAN)
  other_signal <- aggregate(other_onsets, fmrihrf::HRF_SPMG1)
  Y <- matrix(target_signal + 8 * other_signal, T, 1L)
  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = target_onsets, duration = 0, amplitude = 1, span = 30
    ),
    others = list(list(
      onsets = other_onsets, hrf = fmrihrf::HRF_SPMG1,
      duration = 0, amplitude = 1, span = 30
    )),
    precision = 0.1,
    method = "conv"
  )

  selected <- fmrilss:::.oasis_pick_hrf_lwu_fast(
    Y, spec, grid, confounds = matrix(1, T, 1L)
  )
  selected_redundant <- fmrilss:::.oasis_pick_hrf_lwu_fast(
    Y, spec, grid, confounds = cbind(rep(1, T), rep(2, T))
  )
  expect_identical(selected, fmrihrf::HRF_GAUSSIAN)
  expect_identical(selected_redundant, selected)

  grid_spec <- spec
  grid_spec$hrf_grid <- grid
  explicit_spec <- spec
  explicit_spec$cond$hrf <- fmrihrf::HRF_GAUSSIAN
  opts <- list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0)
  fit_selected <- lss(
    Y, X = NULL, Z = matrix(1, T, 1L), method = "oasis",
    oasis = utils::modifyList(opts, list(design_spec = grid_spec))
  )
  fit_explicit <- lss(
    Y, X = NULL, Z = matrix(1, T, 1L), method = "oasis",
    oasis = utils::modifyList(opts, list(design_spec = explicit_spec))
  )
  expect_equal(fit_selected, fit_explicit, tolerance = 1e-12)
})

test_that("HRF grid selection rejects a candidate in the common span", {
  skip_if_not_installed("fmrihrf")
  T <- 100L
  sframe <- fmrihrf::sampling_frame(blocklens = T, TR = 1)
  onsets <- c(10, 30, 50, 70)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = onsets, duration = 0, amplitude = 1, span = 30),
    others = list(list(
      onsets = onsets, hrf = fmrihrf::HRF_SPMG1,
      duration = 0, amplitude = 1, span = 30
    )),
    precision = 0.1,
    method = "conv"
  )
  expect_error(
    fmrilss:::.oasis_pick_hrf_lwu_fast(
      matrix(rnorm(T), T, 1L), spec,
      list(fmrihrf::HRF_SPMG1, fmrihrf::HRF_GAUSSIAN),
      confounds = matrix(1, T, 1L)
    ),
    "zero residual norm"
  )
})

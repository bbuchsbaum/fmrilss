test_that("SBHM library and matching inputs fail closed", {
  skip_if_not_installed("fmrihrf")
  tgrid <- 0:39
  H <- cbind(a = dgamma(tgrid, 5, 1), b = dgamma(tgrid, 8, 1))

  expect_error(sbhm_build(library_H = H, tgrid = tgrid, r = 1.5), "positive integer")
  expect_error(
    sbhm_build(library_H = cbind(H, zero = 0), tgrid = tgrid, r = 2),
    "non-zero finite energy"
  )
  expect_error(
    sbhm_build(
      library_H = cbind(H, zero = 0), tgrid = tgrid, r = 2,
      normalize = FALSE
    ),
    "non-zero finite energy"
  )
  H_dup_names <- H
  colnames(H_dup_names) <- c("same", "same")
  expect_error(
    sbhm_build(library_H = H_dup_names, tgrid = tgrid, r = 2),
    "complete and unique"
  )

  expect_warning(
    rank_one <- sbhm_build(
      library_H = cbind(a = H[, 1], b = 2 * H[, 1]),
      tgrid = tgrid, r = 2, baseline = NULL
    ),
    "numerical rank 1"
  )
  expect_identical(rank_one$meta$r, 1L)
  expect_equal(rank_one$meta$retained_energy_fraction, 1, tolerance = 1e-12)

  built <- sbhm_build(library_H = H, tgrid = tgrid, r = 2, baseline = NULL)
  beta <- matrix(
    c(1, 0.2), 2, 1,
    dimnames = list(rownames(built$A), "voxel_a")
  )
  expect_error(sbhm_match(beta, built$S, built$A, topK = 1.5), "topK")
  expect_error(sbhm_match(beta, built$S, built$A, whiten_power = 1.2), "whiten_power")
  expect_error(sbhm_match(beta, built$S, built$A, shrink = list(tau = 1.2)), "shrink\\$tau")
  expect_error(sbhm_match(beta * 0, built$S, built$A), "zero-norm voxel")
  expect_error(sbhm_match(beta + NA_real_, built$S, built$A), "finite")

  beta_named <- built$A[, 1, drop = FALSE]
  colnames(beta_named) <- "voxel_a"
  matched_reference <- sbhm_match(beta_named, built$S, built$A, whiten = FALSE)
  matched_permuted <- sbhm_match(
    beta_named[rev(rownames(beta_named)), , drop = FALSE],
    built$S[rev(names(built$S))],
    built$A,
    whiten = FALSE
  )
  expect_identical(matched_permuted$idx, matched_reference$idx)
  expect_error(
    sbhm_match(unname(beta_named), built$S, built$A, whiten = FALSE),
    "basis names must be supplied"
  )
})

test_that("factorized SBHM prepass is dimensionally honest for q smaller than V", {
  skip_if_not_installed("fmrihrf")
  set.seed(81101)
  Tlen <- 90L
  q <- 3L
  V <- 8L
  sframe <- fmrihrf::sampling_frame(blocklens = Tlen, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(
    fast = dgamma(times, 5, 1),
    medium = dgamma(times, 7, 1),
    slow = dgamma(times, 9, 1)
  )
  sbhm <- sbhm_build(library_H = H, tgrid = times, r = 3, baseline = NULL)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = c(8, 24, 40, 56, 72), duration = 0, span = 30)
  )
  Scores <- matrix(rnorm(Tlen * q), Tlen, q)
  Load <- matrix(rnorm(q * V), q, V)
  colnames(Scores) <- paste0("factor_", seq_len(q))
  rownames(Load) <- colnames(Scores)
  Y <- Scores %*% Load
  colnames(Y) <- paste0("v", seq_len(V))
  colnames(Load) <- colnames(Y)

  dense <- sbhm_prepass(Y, sbhm, spec, ridge = list(mode = "absolute", lambda = 0))
  factorized <- sbhm_prepass(
    Y, sbhm, spec, ridge = list(mode = "absolute", lambda = 0),
    data_fac = list(scores = Scores, loadings = Load)
  )
  expect_equal(factorized$beta_bar, dense$beta_bar, tolerance = 1e-10)
  reversed <- sbhm_prepass(
    Y, sbhm, spec, ridge = list(mode = "absolute", lambda = 0),
    data_fac = list(scores = Scores, loadings = Load[, V:1, drop = FALSE])
  )
  expect_equal(reversed$beta_bar, dense$beta_bar, tolerance = 1e-10)
  factor_reversed <- sbhm_prepass(
    Y, sbhm, spec, ridge = list(mode = "absolute", lambda = 0),
    data_fac = list(scores = Scores, loadings = Load[q:1, , drop = FALSE])
  )
  expect_equal(factor_reversed$beta_bar, dense$beta_bar, tolerance = 1e-10)
  expect_error(
    sbhm_prepass(Y, sbhm, spec, data_fac = list(scores = Scores, loadings = t(Load))),
    "scores.*loadings"
  )
  expect_error(
    sbhm_prepass(
      Y, sbhm, spec,
      prewhiten = list(method = "ar", p = 1L),
      data_fac = list(scores = Scores, loadings = Load)
    ),
    "does not support active prewhitening"
  )

  out <- lss_sbhm(
    Y, sbhm, spec,
    prepass = list(
      ridge = list(mode = "absolute", lambda = 0),
      data_fac = list(scores = Scores, loadings = Load)
    ),
    match = list(topK = 1, soft_blend = FALSE),
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
    amplitude = list(method = "global_ls", ridge = 0),
    return = "both"
  )
  expect_equal(dim(out$alpha_coords), c(3L, V))
  expect_equal(dim(out$amplitude), c(5L, V))
  expect_equal(dim(out$coeffs_r), c(3L, 5L, V))
  expect_identical(colnames(out$alpha_coords), colnames(Y))
})

test_that("SBHM designs are run-safe and retain fmridesign event roles", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("fmridesign")
  sframe <- fmrihrf::sampling_frame(blocklens = c(50, 60), TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(short = dgamma(times, 5, 1), long = dgamma(times, 8, 1))
  sbhm <- sbhm_build(library_H = H, tgrid = times, r = 2, baseline = NULL)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = c(48.5, 8), duration = c(2, 4), run = c(1, 2), span = 30)
  )
  built <- fmrilss:::.sbhm_build_design(sbhm, spec)
  expect_equal(sum(abs(built$regs[[1]][51:110, ])), 0, tolerance = 0)
  expect_equal(sum(abs(built$regs[[2]][1:50, ])), 0, tolerance = 0)
  expect_equal(built$events$duration, c(2, 4))
  expect_equal(unname(colSums(built$intercepts)), c(50, 60))
  run_safe_pw <- fmrilss:::.sbhm_run_safe_prewhiten(
    list(method = "ar", p = 1L), sframe
  )
  expect_identical(run_safe_pw$runs, fmrihrf::blockids(sframe))
  set.seed(81105)
  pw <- fmrilss:::.sbhm_prewhiten(
    matrix(rnorm(110 * 2), 110, 2), built$regs, built$intercepts,
    Nuisance = NULL, prewhiten = run_safe_pw, X_other = NULL
  )
  expect_equal(sum(abs(pw$regs_w[[1]][51:110, ])), 0, tolerance = 0)
  expect_equal(sum(abs(pw$regs_w[[2]][1:50, ])), 0, tolerance = 0)
  Y_pw <- matrix(rnorm(110 * 2), 110, 2)
  master <- fmrilss:::.prewhiten_data(
    Y_pw, built$X_trials, NULL, built$intercepts, run_safe_pw
  )
  reused_opts <- run_safe_pw
  reused_opts$.whiten_plan <- master$whiten_plan
  aggregate <- fmrilss:::.prewhiten_data(
    Y_pw, built$A_agg, NULL, built$intercepts, reused_opts
  )
  aggregate_from_trials <- do.call(cbind, lapply(seq_len(built$K), function(k) {
    rowSums(master$X_whitened[, seq.int(
      k, ncol(master$X_whitened), by = built$K
    ), drop = FALSE])
  }))
  expect_equal(aggregate$Y_whitened, master$Y_whitened, tolerance = 0)
  expect_equal(
    unname(aggregate$X_whitened), unname(aggregate_from_trials), tolerance = 1e-12
  )
  expect_error(
    fmrilss:::.sbhm_run_safe_prewhiten(
      list(method = "ar", p = 1L, runs = rep(1L, 110L)), sframe
    ),
    "must match the sampling-frame run boundaries"
  )
  zero_amplitude <- spec
  zero_amplitude$cond$amplitude <- c(1, 0)
  expect_error(
    fmrilss:::.sbhm_build_design(sbhm, zero_amplitude),
    "amplitude values must be non-zero"
  )

  events <- data.frame(
    onset = c(10, 28, 12, 32),
    duration = c(1, 2, 3, 4),
    RT = c(0.4, 1.2, 0.7, 1.5),
    run = c(1, 1, 2, 2)
  )
  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1") + fmridesign::hrf(RT),
    data = events, block = ~run, sampling_frame = sframe,
    durations = events$duration
  )
  prepared <- fmrilss:::.prepare_fmridesign_lss(
    fmridesign::design_matrix(emod), emod
  )
  mapped <- fmrilss:::.sbhm_design_spec_from_event_model(emod, sbhm, prepared)
  expect_equal(mapped$cond$onsets, events$onset)
  expect_equal(mapped$cond$duration, events$duration)
  expect_equal(mapped$cond$run, events$run)
  expect_equal(ncol(prepared$fixed), 1L)

  set.seed(81102)
  Y <- matrix(rnorm(110 * 2), 110, 2)
  Y <- Y + 4 * prepared$fixed[, 1]
  with_fixed <- lss_sbhm_design(
    Y, sbhm, emod,
    match = list(topK = 1, soft_blend = FALSE),
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
    amplitude = list(method = "global_ls", ridge = 0),
    return = "both"
  )
  emod_no_fixed <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = events, block = ~run, sampling_frame = sframe,
    durations = events$duration
  )
  without_fixed <- lss_sbhm_design(
    Y, sbhm, emod_no_fixed,
    match = list(topK = 1, soft_blend = FALSE),
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
    amplitude = list(method = "global_ls", ridge = 0),
    return = "both"
  )
  expect_gt(max(abs(with_fixed$amplitude - without_fixed$amplitude)), 1e-4)
  expect_identical(attr(with_fixed, "trial_basis_map"), with_fixed$trial_basis_map)
})

test_that("SBHM fallback coordinates fail closed when the shape is undefined", {
  skip_if_not_installed("fmrihrf")
  set.seed(81104)
  Tlen <- 70L
  sframe <- fmrihrf::sampling_frame(blocklens = Tlen, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(a = dgamma(times, 5, 1), b = dgamma(times, 8, 1))
  sbhm <- sbhm_build(library_H = H, tgrid = times, r = 2, baseline = NULL)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = c(8, 24, 40, 56), duration = 0, span = 30)
  )
  Y <- matrix(rnorm(Tlen * 2L), Tlen, 2L)
  base_args <- list(
    Y = Y, sbhm = sbhm, design_spec = spec,
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
    amplitude = list(method = "global_ls", ridge = 0)
  )
  expect_error(
    do.call(lss_sbhm, c(base_args, list(
      match = list(min_margin = 2, fallback_ref = c(0, 0))
    ))),
    "finite non-zero vector"
  )
  expect_error(
    do.call(lss_sbhm, c(base_args, list(
      match = list(min_margin = 2, fallback_ref = c(NA_real_, 1))
    ))),
    "finite non-zero vector"
  )
})

test_that("SBHM override lists reject misspelled option names", {
  skip_if_not_installed("fmrihrf")
  Tlen <- 60L
  sframe <- fmrihrf::sampling_frame(Tlen, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(a = dgamma(times, 5, 1), b = dgamma(times, 8, 1))
  sbhm <- sbhm_build(library_H = H, tgrid = times, r = 2, baseline = NULL)
  spec <- list(sframe = sframe, cond = list(onsets = c(8, 24, 40), span = 30))
  Y <- matrix(rnorm(Tlen * 2), Tlen, 2)
  expect_error(lss_sbhm(Y, sbhm, spec, match = list(soft_blnd = FALSE)),
               "unknown option: soft_blnd")
  expect_error(lss_sbhm(Y, sbhm, spec, amplitude = list(methd = "global_ls")),
               "unknown option: methd")
  expect_error(lss_sbhm(Y, sbhm, spec, prepass = list(rdig = list())),
               "unknown option: rdig")
  expect_error(lss_sbhm(Y, sbhm, spec, oasis = list(ridg_x = 0)),
               "unknown option: ridg_x")
  expect_error(lss_sbhm(Y, sbhm, spec, match = list(shrink = list(tua = 0))),
               "unknown option: tua")
  expect_error(lss_sbhm(Y, sbhm, spec, match = list(min_margin = Inf)),
               "documented domain")
  expect_error(lss_sbhm(Y, sbhm, spec, match = list(whiten = 1)),
               "match\\$whiten must be TRUE or FALSE")
  expect_error(lss_sbhm(Y, sbhm, spec, match = list(orient_ref = NA)),
               "match\\$orient_ref must be TRUE or FALSE")
  expect_error(
    lss_sbhm(Y, sbhm, spec, amplitude = list(return_se = 1)),
    "amplitude\\$return_se must be TRUE or FALSE"
  )
  expect_error(
    lss_sbhm(Y, sbhm, spec, amplitude = list(adaptive = list(enable = "yes"))),
    "amplitude\\$adaptive\\$enable must be TRUE or FALSE"
  )
  expect_error(lss_sbhm(Y, sbhm, spec, oasis = list(return_se = 1)),
               "oasis\\$return_se must be TRUE or FALSE")
  expect_error(lss_sbhm(Y, sbhm, spec, oasis = list(return_se = TRUE)),
               "does not support OASIS standard errors")
  expect_error(lss_sbhm(Y, sbhm, spec, oasis = list(return_diag = TRUE)),
               "does not expose OASIS design diagnostics")
})

test_that("SBHM output preserves candidate, trial, and voxel identity", {
  skip_if_not_installed("fmrihrf")
  set.seed(81103)
  Tlen <- 80L
  sframe <- fmrihrf::sampling_frame(Tlen, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(short = dgamma(times, 5, 1), long = dgamma(times, 8, 1))
  sbhm <- sbhm_build(library_H = H, tgrid = times, r = 2, baseline = NULL)
  spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(10, 30, 50), duration = 0, span = 30,
      trial_names = c("event_a", "event_b", "event_c")
    )
  )
  Y <- matrix(rnorm(Tlen * 2), Tlen, 2, dimnames = list(NULL, c("left", "right")))
  out <- lss_sbhm(
    Y, sbhm, spec,
    match = list(topK = 1, soft_blend = FALSE),
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
    amplitude = list(method = "global_ls", ridge = 0),
    return = "both"
  )
  expect_identical(names(out$matched_idx), colnames(Y))
  expect_true(all(out$matched_name %in% colnames(H)))
  expect_identical(colnames(out$alpha_coords), colnames(Y))
  expect_identical(dimnames(out$amplitude), list(spec$cond$trial_names, colnames(Y)))
  expect_identical(dimnames(out$coeffs_r)[2:3], list(spec$cond$trial_names, colnames(Y)))
  expect_identical(out$trial_basis_map$output_name, out$trial_basis_map$column)
  expect_identical(out$trial_basis_map$run, rep.int(1L, nrow(out$trial_basis_map)))
  expect_identical(unname(out$event_amplitude), rep(1, 3))
  expect_identical(unname(out$event_duration), rep(0, 3))
  expect_identical(unname(out$event_run), rep(1L, 3))
  expect_identical(names(out$prepass_fallback), colnames(Y))

  partial <- lss_sbhm(
    Y, sbhm, spec,
    match = list(alpha_source = "prepass"),
    oasis = list(ridge_mode = "absolute", ridge_x = 0, ridge_b = 0),
    amplitude = list(method = "global_ls", ridge = 0)
  )
  expect_equal(nrow(partial$topK_idx), 2L)
  expect_error(
    lss_sbhm(
      Y, sbhm, spec,
      match = list(topK = 1, soft_blend = FALSE),
      amplitude = list(
        method = "global_ls", ridge = 0,
        cond_gate = list(metric = "not_a_metric", thr = 0, fallback = "lss1")
      )
    ),
    "rho.*kappa"
  )
  gated <- lss_sbhm(
    Y, sbhm, spec,
    match = list(topK = 1, soft_blend = FALSE),
    amplitude = list(
      method = "lss1",
      cond_gate = list(metric = "rho", thr = -1, fallback = "global_ls")
    )
  )
  expect_true(all(gated$diag$method_used == "global_ls"))
})

test_that("SBHM scalar output is a coefficient on the supplied event design", {
  skip_if_not_installed("fmrihrf")
  Tlen <- 100L
  sframe <- fmrihrf::sampling_frame(Tlen, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(short = dgamma(times, 5, 1), long = dgamma(times, 8, 1))
  sbhm <- sbhm_build(library_H = H, tgrid = times, r = 2, baseline = NULL)
  alpha <- sbhm$A[, 1, drop = FALSE]
  spec_one <- list(
    sframe = sframe,
    cond = list(
      onsets = c(12, 42, 72), duration = c(0, 2, 4),
      amplitude = 1, span = 30
    )
  )
  built_one <- fmrilss:::.sbhm_build_design(sbhm, spec_one)
  X_one <- do.call(cbind, lapply(
    built_one$regs, function(Xt) as.numeric(Xt %*% alpha[, 1])
  ))
  truth <- c(3, 3, 3)
  Y <- matrix(2 + X_one %*% truth, Tlen, 1)
  beta_one <- fmrilss:::sbhm_amplitude_ls(
    Y, sbhm, spec_one, alpha, ridge = 0
  )[, 1]
  spec_two <- spec_one
  spec_two$cond$amplitude <- 2
  beta_two <- fmrilss:::sbhm_amplitude_ls(
    Y, sbhm, spec_two, alpha, ridge = 0
  )[, 1]
  expect_equal(beta_one, truth, tolerance = 1e-10)
  expect_equal(beta_two, truth / 2, tolerance = 1e-10)
  expect_equal(beta_one, 2 * beta_two, tolerance = 1e-10)
})

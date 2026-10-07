test_that("lss_sbhm_design basic functionality (single run)", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  # Temporal setup
  sframe <- fmrihrf::sampling_frame(blocklens = 120, TR = 1)
  trials <- data.frame(onset = seq(6, 90, by = 12), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  # Minimal SBHM basis (r = 2) built from a simple time library
  Tlen <- length(fmrihrf::samples(sframe, global = TRUE))
  H <- cbind(
    exp(-seq(0, 3, length.out = Tlen)),
    exp(-seq(0, 6, length.out = Tlen))
  )
  sbhm <- sbhm_build(library_H = H, r = 2, sframe = sframe, normalize = TRUE)

  set.seed(1)
  Y <- matrix(rnorm(Tlen * 8), Tlen, 8)

  out <- lss_sbhm_design(Y, sbhm, emod, return = "coefficients", validate = TRUE)
  expect_true(is.list(out))
  expect_true(!is.null(out$coeffs_r))
  # coeffs_r should be r × ntrials × V
  expect_equal(dim(out$coeffs_r)[1], 2)
  expect_equal(dim(out$coeffs_r)[3], ncol(Y))
  # Metadata
  expect_equal(attr(out, "method"), "lss_sbhm_design")
  expect_true(!is.null(attr(out, "event_model")))
  expect_true(!is.null(attr(out, "sampling_frame")))
})

test_that("lss_sbhm_design rejects malformed validation and reserved dots", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
  trials <- data.frame(onset = c(5, 20, 35), run = 1)
  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )
  H <- cbind(exp(-seq(0, 3, length.out = 60)), exp(-seq(0, 6, length.out = 60)))
  sbhm <- sbhm_build(library_H = H, r = 2, sframe = sframe, normalize = TRUE)
  Y <- matrix(rnorm(60 * 2), 60, 2)

  expect_error(
    lss_sbhm_design(Y, sbhm, emod, validate = NA),
    "validate must be TRUE or FALSE"
  )
  expect_error(
    lss_sbhm_design(Y, sbhm, emod, typo = 1),
    "unsupported argument"
  )
})

test_that("lss_sbhm_design with baseline_model (multi-run)", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = c(80, 80), TR = 1.5)
  trials <- data.frame(
    onset = c(6, 18, 30, 6, 18, 30),
    run = rep(1:2, each = 3)
  )

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  bmodel <- fmridesign::baseline_model(
    basis = "poly", degree = 3, sframe = sframe, intercept = "runwise"
  )

  Tlen <- length(fmrihrf::samples(sframe, global = TRUE))
  H <- cbind(
    exp(-seq(0, 3, length.out = Tlen)),
    exp(-seq(0, 6, length.out = Tlen))
  )
  sbhm <- sbhm_build(library_H = H, r = 2, sframe = sframe, normalize = TRUE)

  set.seed(2)
  Y <- matrix(rnorm(Tlen * 6), Tlen, 6)

  out <- lss_sbhm_design(Y, sbhm, emod, baseline_model = bmodel,
                         return = "coefficients", validate = TRUE)

  expect_true(is.list(out))
  expect_true(!is.null(out$coeffs_r))
  # Multi-run: ntrials is 6
  expect_equal(dim(out$coeffs_r)[2], 6)
  # Metadata is attached when validate=TRUE
  expect_equal(attr(out, "method"), "lss_sbhm_design")
  expect_true(!is.null(attr(out, "event_model")))
  expect_true(!is.null(attr(out, "baseline_model")))
  expect_true(!is.null(attr(out, "sampling_frame")))
})

test_that("SBHM records and applies one run-aware whitening plan", {
  skip_if_not_installed("fmridesign")
  set.seed(1818)
  sframe <- fmrihrf::sampling_frame(c(90L, 90L), TR = 1)
  events <- data.frame(onset = rep(c(5, 24, 45, 65), 2), run = rep(1:2, each = 4))
  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"), data = events,
    block = ~run, sampling_frame = sframe
  )
  time <- seq(0, 30, by = 0.5)
  H <- cbind(dgamma(time, 5, 1), dgamma(time, 7, 1))
  sbhm <- sbhm_build(library_H = H, tgrid = time, span = 30, r = 2, baseline = NULL)
  spec <- fmrilss:::.sbhm_design_spec_from_event_model(emod, sbhm)
  built <- fmrilss:::.sbhm_build_design(sbhm, spec)
  Y <- matrix(0, 180, 3)
  for (v in 1:3) {
    signal <- rowSums(do.call(cbind, lapply(built$regs, function(x) x %*% sbhm$A[, 1])))
    for (rows in list(1:90, 91:180)) {
      Y[rows, v] <- signal[rows] + as.numeric(stats::filter(rnorm(90), 0.75, "recursive"))
    }
  }
  args <- list(Y = Y, sbhm = sbhm, event_model = emod, validate = FALSE,
               oasis = list(ridge_x = 0, ridge_b = 0),
               amplitude = list(method = "global_ls", ridge = list(mode = "absolute", lambda = 0),
                                adaptive = list(enable = FALSE)), return = "both")
  out <- do.call(lss_sbhm_design, c(args, list(prewhiten = list(method = "ar", p = 1, pooling = "run"))))
  plan <- attr(out, "whiten_plan")
  expect_s3_class(plan, "fmriAR_plan")
  expect_identical(attr(out$amplitude, "whiten_plan"), plan)
  expect_identical(attr(out$coeffs_r, "whiten_plan"), plan)
  expect_true(out$diag$prewhitening$requested)
  expect_true(out$diag$prewhitening$applied)
  expect_true(out$diag$prewhitening$response_changed)
  expect_identical(out$diag$prewhitening$residual_model, "aggregate")
  for (v in seq_len(ncol(Y))) {
    X <- do.call(cbind, lapply(built$regs, function(x) x %*% out$alpha_coords[, v]))
    whitened <- fmriAR::whiten_apply(plan, X = cbind(built$intercepts, X),
                                    Y = Y[, v, drop = FALSE], parcels = 1L)
    direct <- as.numeric(stats::lm.fit(whitened$X, whitened$Y)$coefficients)[-seq_len(ncol(built$intercepts))]
    expect_equal(unname(out$amplitude[, v]), unname(direct), tolerance = 1e-8)
  }
  none <- do.call(lss_sbhm_design, c(args, list(prewhiten = list(method = "none"))))
  expect_null(attr(none, "whiten_plan"))
  expect_false(none$diag$prewhitening$requested)
  expect_false(none$diag$prewhitening$applied)
  expect_gt(max(abs(as.numeric(out$amplitude) - as.numeric(none$amplitude))), 1e-5)
  pre <- sbhm_prepass(Y, sbhm, spec, prewhiten = list(method = "none"))
  expect_false(pre$diag$used_prewhiten)
  expect_null(attr(pre, "whiten_plan"))
})

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

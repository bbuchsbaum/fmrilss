context("lss_design integration with fmridesign")

test_that("lss_design requires fmridesign package", {
  skip_if_installed("fmridesign")

  Y <- matrix(rnorm(100 * 10), 100, 10)

  expect_error(
    lss_design(Y, NULL),
    "fmridesign.*required"
  )
})

test_that("lss_design validates event_model input", {
  skip_if_not_installed("fmridesign")

  Y <- matrix(rnorm(100 * 10), 100, 10)

  # Missing event_model
  expect_error(
    lss_design(Y),
    "event_model is required"
  )

  # Wrong class
  expect_error(
    lss_design(Y, list(foo = "bar")),
    "event_model.*object from fmridesign"
  )
})

test_that("lss_design basic functionality", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  # Setup
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Test basic call
  beta <- lss_design(Y, emod, method = "oasis")

  expect_equal(dim(beta), c(5, 50))  # 5 trials, 50 voxels
  expect_true(!any(is.na(beta)))
  expect_true(is.matrix(beta))
})

test_that("lss_design with baseline_model", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  bmodel <- fmridesign::baseline_model(
    basis = "bs",
    degree = 5,
    sframe = sframe,
    intercept = "runwise"
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)
  beta <- lss_design(Y, emod, bmodel, method = "oasis")

  expect_equal(dim(beta), c(5, 50))
  expect_true(!any(is.na(beta)))
})

test_that("lss_design with nuisance regressors", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  # Add motion parameters
  motion <- matrix(rnorm(100 * 6), 100, 6)
  bmodel <- fmridesign::baseline_model(
    basis = "bs",
    degree = 5,
    sframe = sframe,
    nuisance_list = list(motion)
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)
  beta <- lss_design(Y, emod, bmodel, method = "oasis")

  expect_equal(dim(beta), c(5, 50))
})

test_that("lss_design multi-run handling", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2)
  trials <- data.frame(
    onset = c(10, 30, 50, 10, 30, 50),
    run = rep(1:2, each = 3)
  )

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(200 * 50), 200, 50)
  beta <- lss_design(Y, emod, method = "oasis")

  expect_equal(nrow(beta), 6)  # 6 trials across 2 runs
  expect_equal(ncol(beta), 50)
})

test_that("lss_design multi-run with baseline", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2)
  trials <- data.frame(
    onset = c(10, 30, 50, 10, 30, 50),
    run = rep(1:2, each = 3)
  )

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  # Baseline with runwise intercepts and motion per run
  motion_r1 <- matrix(rnorm(100 * 6), 100, 6)
  motion_r2 <- matrix(rnorm(100 * 6), 100, 6)

  bmodel <- fmridesign::baseline_model(
    basis = "poly",
    degree = 3,
    sframe = sframe,
    intercept = "runwise",
    nuisance_list = list(motion_r1, motion_r2)
  )

  Y <- matrix(rnorm(200 * 50), 200, 50)
  beta <- lss_design(Y, emod, bmodel, method = "oasis")

  expect_equal(nrow(beta), 6)
  expect_equal(ncol(beta), 50)
})

test_that("lss_design validates temporal alignment", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  # Wrong number of timepoints
  Y_wrong <- matrix(rnorm(150 * 50), 150, 50)

  expect_error(
    lss_design(Y_wrong, emod, method = "oasis"),
    "150 rows but sampling_frame expects 100"
  )
})

test_that("lss_design validates sampling_frame consistency", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe1 <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  sframe2 <- fmrihrf::sampling_frame(blocklens = 100, TR = 1.5)  # Different TR

  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe1
  )

  bmodel <- fmridesign::baseline_model(
    basis = "poly",
    degree = 3,
    sframe = sframe2  # Different sampling frame!
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  expect_error(
    lss_design(Y, emod, bmodel, method = "oasis"),
    "same sampling_frame"
  )
})

test_that("lss_design detects multi-basis K", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg3", nbasis = 3),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Should auto-detect K = 3
  expect_message(
    beta <- lss_design(Y, emod, method = "oasis"),
    "Detected multi-basis HRF with K = 3"
  )

  expect_equal(nrow(beta), 15)  # 5 trials × 3 basis = 15
})

test_that("lss_design manual K specification", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg3", nbasis = 3),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Manual K specification should work
  beta <- lss_design(Y, emod, method = "oasis", oasis = list(K = 3))

  expect_equal(nrow(beta), 15)
})

test_that("lss_design attaches metadata", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  bmodel <- fmridesign::baseline_model(
    basis = "bs",
    degree = 5,
    sframe = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)
  beta <- lss_design(Y, emod, bmodel, method = "oasis")

  # Check metadata
  expect_true(!is.null(attr(beta, "event_model")))
  expect_true(!is.null(attr(beta, "baseline_model")))
  expect_true(!is.null(attr(beta, "sampling_frame")))
  expect_equal(attr(beta, "method"), "lss_design")
})

test_that("lss_design equivalent to manual lss() call", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  # Setup
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  bmodel <- fmridesign::baseline_model(
    basis = "bs",
    degree = 5,
    sframe = sframe,
    intercept = "runwise"
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Method 1: lss_design
  beta1 <- lss_design(Y, emod, bmodel, method = "oasis", validate = FALSE)

  # Method 2: Manual extraction + lss (using term_matrices)
  X <- as.matrix(fmridesign::design_matrix(emod))
  bm_term_mats <- fmridesign::term_matrices(bmodel)
  Z <- as.matrix(cbind(bm_term_mats$drift, bm_term_mats$block))
  beta2 <- lss(Y, X, Z, method = "oasis")

  expect_equal(beta1, beta2, tolerance = 1e-10)
})

test_that("lss_design with ridge regularization", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # With ridge regularization
  beta <- lss_design(
    Y, emod,
    method = "oasis",
    oasis = list(
      ridge_mode = "fractional",
      ridge_x = 0.02,
      ridge_b = 0.02
    )
  )

  expect_equal(dim(beta), c(5, 50))
  expect_true(!any(is.na(beta)))
})

test_that("lss_design can skip validation", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Should run without validation checks
  beta <- lss_design(Y, emod, method = "oasis", validate = FALSE)

  expect_equal(dim(beta), c(5, 50))
})

test_that("lss_design warns about collinearity", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 2)

  # Create events that are very close together (high collinearity)
  trials <- data.frame(onset = seq(1, 9, by = 2), run = 1)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(50 * 10), 50, 10)

  # May warn about collinearity
  expect_warning(
    beta <- lss_design(Y, emod, method = "oasis"),
    "collinearity.*ridge",
    ignore.case = TRUE
  )
})

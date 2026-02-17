# Tests for voxel_hrf.R to improve coverage

test_that("estimate_voxel_hrf validates Y as numeric matrix", {
  skip_if_not_installed("fmrihrf")

  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")
  basis <- fmrihrf::HRF_SPMG1

  expect_error(estimate_voxel_hrf("not_matrix", events, basis),
               "Y must be a numeric matrix")
  expect_error(estimate_voxel_hrf(list(1, 2), events, basis),
               "Y must be a numeric matrix")
})

test_that("estimate_voxel_hrf validates events data.frame", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 3), 50, 3)
  basis <- fmrihrf::HRF_SPMG1

  # Missing columns
  expect_error(estimate_voxel_hrf(Y, data.frame(onset = 5), basis),
               "events must be a data.frame with columns")

  # Not a data.frame
  expect_error(estimate_voxel_hrf(Y, list(onset = 5), basis),
               "events must be a data.frame with columns")
})

test_that("estimate_voxel_hrf validates basis type", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 3), 50, 3)
  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")

  expect_error(estimate_voxel_hrf(Y, events, "not_hrf"),
               "basis must be an 'HRF' object")
})

test_that("estimate_voxel_hrf validates nuisance_regs", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 3), 50, 3)
  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")
  basis <- fmrihrf::HRF_SPMG1

  # Not a matrix
  expect_error(estimate_voxel_hrf(Y, events, basis, nuisance_regs = "not_matrix"),
               "nuisance_regs must be a numeric matrix")

  # Wrong number of rows
  expect_error(estimate_voxel_hrf(Y, events, basis, nuisance_regs = matrix(1:10, 10, 1)),
               "nuisance_regs must have same number of rows as Y")
})

test_that("estimate_voxel_hrf returns VoxelHRF object", {
  skip_if_not_installed("fmrihrf")

  set.seed(123)
  n_time <- 50
  n_vox <- 3
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  events <- data.frame(onset = c(5, 15, 30), duration = 1, condition = "A")
  basis <- fmrihrf::HRF_SPMG1

  result <- estimate_voxel_hrf(Y, events, basis)

  expect_s3_class(result, "VoxelHRF")
  expect_true("coefficients" %in% names(result))
  expect_true("basis" %in% names(result))
  expect_true("conditions" %in% names(result))
  expect_true(is.matrix(result$coefficients))
})

test_that("estimate_voxel_hrf works with nuisance regressors", {
  skip_if_not_installed("fmrihrf")

  set.seed(234)
  n_time <- 60
  n_vox <- 4
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  events <- data.frame(onset = c(5, 20, 40), duration = 1, condition = "A")
  basis <- fmrihrf::HRF_SPMG1
  nuisance <- cbind(1:n_time, (1:n_time)^2)  # Linear and quadratic trends

  result <- estimate_voxel_hrf(Y, events, basis, nuisance_regs = nuisance)

  expect_s3_class(result, "VoxelHRF")
})

test_that("lss_with_hrf validates Y", {
  skip_if_not_installed("fmrihrf")

  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")

  # Create mock VoxelHRF object
  mock_hrf <- list(
    coefficients = matrix(1, 3, 2),
    basis = fmrihrf::HRF_SPMG1,
    conditions = "A"
  )
  class(mock_hrf) <- "VoxelHRF"

  expect_error(lss_with_hrf("not_matrix", events, mock_hrf),
               "Y must be a numeric matrix")
})

test_that("lss_with_hrf validates events", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 2), 50, 2)

  mock_hrf <- list(
    coefficients = matrix(1, 3, 2),
    basis = fmrihrf::HRF_SPMG1,
    conditions = "A"
  )
  class(mock_hrf) <- "VoxelHRF"

  expect_error(lss_with_hrf(Y, data.frame(onset = 5), mock_hrf),
               "events must be a data.frame with columns")
})

test_that("lss_with_hrf validates hrf_estimates", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 2), 50, 2)
  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")

  expect_error(lss_with_hrf(Y, events, "not_voxelhrf"),
               "hrf_estimates must be a VoxelHRF object")

  # Missing coefficients
  bad_hrf <- list(basis = fmrihrf::HRF_SPMG1)
  class(bad_hrf) <- "VoxelHRF"
  expect_error(lss_with_hrf(Y, events, bad_hrf),
               "hrf_estimates must be a VoxelHRF object")
})

test_that("lss_with_hrf validates nuisance_regs", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 2), 50, 2)
  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")

  mock_hrf <- list(
    coefficients = matrix(1, 2, 2),
    basis = fmrihrf::HRF_SPMG1,
    conditions = "A"
  )
  class(mock_hrf) <- "VoxelHRF"

  expect_error(lss_with_hrf(Y, events, mock_hrf, nuisance_regs = "not_matrix"),
               "nuisance_regs must be a numeric matrix")
})

test_that("lss_with_hrf validates chunk_size", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(50 * 2), 50, 2)
  events <- data.frame(onset = c(5, 15), duration = 1, condition = "A")

  mock_hrf <- list(
    coefficients = matrix(1, 2, 2),
    basis = fmrihrf::HRF_SPMG1,
    conditions = "A"
  )
  class(mock_hrf) <- "VoxelHRF"

  expect_error(lss_with_hrf(Y, events, mock_hrf, chunk_size = -1),
               "chunk_size must be a positive integer")
  expect_error(lss_with_hrf(Y, events, mock_hrf, chunk_size = "abc"),
               "chunk_size must be a positive integer")
})

test_that("lss_with_hrf R engine works", {
  # Skip: requires specific coordinate between estimate_voxel_hrf output
  # and lss_with_hrf expectations that are tested elsewhere
  skip("lss_with_hrf R engine tested in test-lss-with-hrf.R")
})

test_that("lss_with_hrf respects verbose parameter", {
  # Skip: requires specific coordinate between estimate_voxel_hrf output
  # and lss_with_hrf expectations that are tested elsewhere
  skip("lss_with_hrf verbose tested in test-lss-with-hrf.R")
})

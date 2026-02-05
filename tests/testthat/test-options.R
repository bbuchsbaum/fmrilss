test_that("oasis_options creates valid options with defaults", {
  opts <- oasis_options()

  expect_s3_class(opts, "fmrilss_oasis_options")
  expect_s3_class(opts, "list")

  # Check defaults
  expect_null(opts$design_spec)
  expect_null(opts$K)
  expect_equal(opts$ridge_mode, "absolute")
  expect_equal(opts$ridge_x, 0)
  expect_equal(opts$ridge_b, 0)
  expect_equal(opts$block_cols, 4096L)
  expect_false(opts$return_se)
  expect_false(opts$return_diag)
  expect_true(opts$add_intercept)
  expect_null(opts$hrf_mode)
})

test_that("oasis_options accepts custom values", {
  opts <- oasis_options(
    K = 3,
    ridge_mode = "fractional",
    ridge_x = 0.01,
    ridge_b = 0.05,
    block_cols = 2048L,
    return_se = TRUE,
    return_diag = TRUE,
    add_intercept = FALSE
  )

  expect_equal(opts$K, 3L)
  expect_equal(opts$ridge_mode, "fractional")
  expect_equal(opts$ridge_x, 0.01)
  expect_equal(opts$ridge_b, 0.05)
  expect_equal(opts$block_cols, 2048L)
  expect_true(opts$return_se)
  expect_true(opts$return_diag)
  expect_false(opts$add_intercept)
})

test_that("oasis_options validates ridge_x", {
  expect_error(oasis_options(ridge_x = -1), "ridge_x must be a non-negative scalar")
  expect_error(oasis_options(ridge_x = "abc"), "ridge_x must be a non-negative scalar")
  expect_error(oasis_options(ridge_x = c(0.1, 0.2)), "ridge_x must be a non-negative scalar")
})

test_that("oasis_options validates ridge_b", {
  expect_error(oasis_options(ridge_b = -0.5), "ridge_b must be a non-negative scalar")
  expect_error(oasis_options(ridge_b = "abc"), "ridge_b must be a non-negative scalar")
  expect_error(oasis_options(ridge_b = c(0.1, 0.2)), "ridge_b must be a non-negative scalar")
})

test_that("oasis_options validates block_cols", {
  expect_error(oasis_options(block_cols = 0), "block_cols must be a positive integer")
  expect_error(oasis_options(block_cols = -1), "block_cols must be a positive integer")
  expect_error(oasis_options(block_cols = NA), "block_cols must be a positive integer")
})

test_that("oasis_options validates ridge_mode", {
  expect_error(oasis_options(ridge_mode = "invalid"), "'arg' should be one of")
})

test_that("oasis_options allows extra arguments via ...", {
  opts <- oasis_options(custom_field = "test_value", another = 123)

  expect_equal(opts$custom_field, "test_value")
  expect_equal(opts$another, 123)
})

test_that("prewhiten_options creates valid options with defaults", {
  opts <- prewhiten_options()

  expect_s3_class(opts, "fmrilss_prewhiten_options")
  expect_s3_class(opts, "list")

  # Check defaults
  expect_equal(opts$method, "none")
  expect_equal(opts$p, "auto")
  expect_equal(opts$q, 0L)
  expect_equal(opts$p_max, 6L)
  expect_equal(opts$pooling, "global")
  expect_null(opts$runs)
  expect_null(opts$parcels)
  expect_equal(opts$exact_first, "ar1")
  expect_true(opts$compute_residuals)
})

test_that("prewhiten_options accepts custom values", {
  opts <- prewhiten_options(
    method = "ar",
    p = 2,
    q = 1L,
    p_max = 8L,
    pooling = "voxel",
    runs = c(1, 1, 2, 2),
    parcels = c(1, 2, 1, 2),
    exact_first = "none",
    compute_residuals = FALSE
  )

  expect_equal(opts$method, "ar")
  expect_equal(opts$p, 2)
  expect_equal(opts$q, 1L)
  expect_equal(opts$p_max, 8L)
  expect_equal(opts$pooling, "voxel")
  expect_equal(opts$runs, c(1, 1, 2, 2))
  expect_equal(opts$parcels, c(1, 2, 1, 2))
  expect_equal(opts$exact_first, "none")
  expect_false(opts$compute_residuals)
})

test_that("prewhiten_options validates method", {
  expect_error(prewhiten_options(method = "invalid"), "'arg' should be one of")

  # Valid methods should work
  expect_no_error(prewhiten_options(method = "none"))
  expect_no_error(prewhiten_options(method = "ar"))
  expect_no_error(prewhiten_options(method = "arma"))
})

test_that("prewhiten_options validates pooling", {
  expect_error(prewhiten_options(pooling = "invalid"), "'arg' should be one of")

  # Valid pooling options should work
  expect_no_error(prewhiten_options(pooling = "global"))
  expect_no_error(prewhiten_options(pooling = "voxel"))
  expect_no_error(prewhiten_options(pooling = "run"))
  expect_no_error(prewhiten_options(pooling = "parcel"))
})

test_that("prewhiten_options validates exact_first", {
  expect_error(prewhiten_options(exact_first = "invalid"), "'arg' should be one of")

  # Valid exact_first options should work
  expect_no_error(prewhiten_options(exact_first = "ar1"))
  expect_no_error(prewhiten_options(exact_first = "none"))
})

test_that("prewhiten_options coerces q and p_max to integer", {
  opts <- prewhiten_options(q = 2.5, p_max = 5.9)

  expect_type(opts$q, "integer")
  expect_type(opts$p_max, "integer")
  expect_equal(opts$q, 2L)
  expect_equal(opts$p_max, 5L)
})

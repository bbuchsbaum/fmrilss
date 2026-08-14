test_that("oasis_options creates valid options with defaults", {
  opts <- oasis_options()

  expect_s3_class(opts, "fmrilss_oasis_options")
  expect_s3_class(opts, "list")

  # Check defaults
  expect_null(opts$design_spec)
  expect_null(opts$K)
  expect_null(opts$ntrials)
  expect_null(opts$trial_basis_map)
  expect_equal(opts$ridge_mode, "fractional")
  expect_equal(opts$ridge_x, 0.05)
  expect_equal(opts$ridge_b, 0.05)
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
    ridge_x = 0,
    ridge_b = 0,
    block_cols = 2048L,
    return_se = TRUE,
    return_diag = TRUE,
    add_intercept = FALSE
  )

  expect_equal(opts$K, 3L)
  expect_equal(opts$ridge_mode, "fractional")
  expect_equal(opts$ridge_x, 0)
  expect_equal(opts$ridge_b, 0)
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
  expect_error(oasis_options(block_cols = 1.5), "block_cols must be a positive integer")
  expect_error(oasis_options(block_cols = Inf), "block_cols must be a positive integer")
})

test_that("oasis_options validates K", {
  expect_error(oasis_options(K = 0), "K must be a positive integer")
  expect_error(oasis_options(K = 1.5), "K must be a positive integer")
  expect_error(oasis_options(K = NA), "K must be a positive integer")
  expect_equal(oasis_options(K = 3)$K, 3L)
})

test_that("oasis_options rejects penalized standard errors", {
  expect_error(
    oasis_options(return_se = TRUE),
    "return_se requires ridge_x = ridge_b = 0"
  )
})

test_that("oasis_options validates ridge_mode", {
  expect_error(oasis_options(ridge_mode = "invalid"), "'arg' should be one of")
})

test_that("oasis_options rejects unknown and malformed controls", {
  expect_error(oasis_options(rdig_x = 0), "unknown option")
  expect_error(oasis_options(return_se = 1), "return_se must be TRUE or FALSE")
  expect_error(oasis_options(return_diag = NA), "return_diag must be TRUE or FALSE")
  expect_error(oasis_options(add_intercept = "yes"), "add_intercept must be TRUE or FALSE")
  expect_error(oasis_options(orient_ref = 1), "orient_ref must be TRUE or FALSE")
})

test_that("prewhiten_options creates valid options with defaults", {
  opts <- prewhiten_options()

  expect_s3_class(opts, "fmrilss_prewhiten_options")
  expect_s3_class(opts, "list")
  expect_named(opts, fmrilss:::.prewhiten_option_names(internal = FALSE))

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
  expect_null(opts$design)
  expect_null(opts$acvf_correction)
  expect_equal(opts$correction_max_lag, 25L)
})

test_that("prewhiten_options accepts custom values", {
  opts <- prewhiten_options(
    method = "ar",
    p = 2,
    q = 0L,
    p_max = 8L,
    pooling = "voxel",
    runs = c(1, 1, 2, 2),
    parcels = c(1, 2, 1, 2),
    exact_first = "none",
    compute_residuals = FALSE
  )

  expect_equal(opts$method, "ar")
  expect_equal(opts$p, 2)
  expect_equal(opts$q, 0L)
  expect_equal(opts$p_max, 8L)
  expect_equal(opts$pooling, "voxel")
  expect_equal(opts$runs, c(1, 1, 2, 2))
  expect_equal(opts$parcels, c(1, 2, 1, 2))
  expect_equal(opts$exact_first, "none")
  expect_false(opts$compute_residuals)
})

test_that("prewhiten_options accepts residual-bias correction controls", {
  opts <- prewhiten_options(
    method = "ar",
    design = diag(4),
    correction_max_lag = 12L
  )

  expect_equal(opts$design, diag(4))
  expect_null(opts$acvf_correction)
  expect_equal(opts$correction_max_lag, 12L)
})

test_that("prewhiten_options validates method", {
  expect_error(prewhiten_options(method = "invalid"), "'arg' should be one of")

  # Valid methods should work
  expect_no_error(prewhiten_options(method = "none"))
  expect_no_error(prewhiten_options(method = "ar"))
  expect_no_error(prewhiten_options(method = "arma", q = 1))
})

test_that("prewhiten_options validates pooling", {
  expect_error(prewhiten_options(pooling = "invalid"), "'arg' should be one of")

  # Valid pooling options should work
  expect_no_error(prewhiten_options(pooling = "global"))
  expect_no_error(prewhiten_options(pooling = "voxel"))
  expect_no_error(prewhiten_options(pooling = "run", runs = c(1, 1)))
  expect_no_error(prewhiten_options(pooling = "parcel", parcels = c(1, 2)))
})

test_that("prewhiten_options validates exact_first", {
  expect_error(prewhiten_options(exact_first = "invalid"), "'arg' should be one of")

  # Valid exact_first options should work
  expect_no_error(prewhiten_options(exact_first = "ar1"))
  expect_no_error(prewhiten_options(exact_first = "none"))
})

test_that("prewhiten_options rejects lossy or malformed controls", {
  expect_error(prewhiten_options(q = 2.5), "q must be a non-negative integer")
  expect_error(prewhiten_options(p_max = 5.9), "p_max must be a positive integer")
  expect_error(prewhiten_options(p = 1.5), "p must be a positive integer")
  expect_error(prewhiten_options(compute_residuals = 1), "must be TRUE or FALSE")
  expect_error(prewhiten_options(correction_max_lag = 5.5), "positive integer")
  expect_error(prewhiten_options(method = "arma", q = 0), "q must be a positive integer")
  expect_error(prewhiten_options(method = "ar", q = 1), "q must be 0")
  expect_error(prewhiten_options(pooling = "run"), "runs must be supplied")
  expect_error(prewhiten_options(pooling = "parcel"), "parcels must be supplied")
  expect_error(
    prewhiten_options(pooling = "parcel", parcels = c(1.9, 2.1)),
    "exact integer identifiers"
  )
  expect_error(
    prewhiten_options(method = "ar", design = diag(3), acvf_correction = diag(3)),
    "either design or acvf_correction"
  )
  expect_error(
    prewhiten_options(method = "ar", design = matrix("x", 2, 2)),
    "finite numeric matrix"
  )
  expect_error(
    prewhiten_options(method = "ar", acvf_correction = matrix(1, 2, 3)),
    "finite square numeric matrices"
  )
  expect_error(
    prewhiten_options(method = "arma", q = 1L, design = diag(3)),
    "requires method = 'ar'"
  )
  expect_error(
    prewhiten_options(
      method = "ar", pooling = "parcel", parcels = 1:3, design = diag(3)
    ),
    "requires pooling = 'global' or 'run'"
  )
})

test_that("execution paths revalidate direct option lists", {
  set.seed(42)
  Y <- matrix(rnorm(80), 40, 2)
  X <- matrix(rnorm(120), 40, 3)

  expect_error(
    lss(Y, X, method = "oasis", oasis = list(rdig_x = 0)),
    "unknown option"
  )
  expect_error(
    lss(Y, X, method = "oasis", oasis = list(
      return_se = 1, ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
    )),
    "oasis\\$return_se must be TRUE or FALSE"
  )
  expect_error(
    lss(Y, X, prewhiten = list(method = "none", typo = 1)),
    "unknown option"
  )
  expect_error(
    lss(Y, X, prewhiten = list(method = "ar", p = 1.5)),
    "prewhiten\\$p must be a positive integer"
  )
  expect_no_error(
    lss(
      Y, X, method = "oasis",
      prewhiten = list(
        method = "ar", p = 1L, design = cbind(1, X),
        correction_max_lag = 5L
      )
    )
  )
  expect_no_error(
    lss(
      Y, X, method = "oasis",
      prewhiten = prewhiten_options(method = "ar", p = 1L)
    )
  )
})

test_that("stglmnet options reject silent coercion and unknown fields", {
  expect_error(stglmnet_options(typo = 1), "unknown option")
  expect_error(stglmnet_options(pool_to_mean = 1), "must be TRUE or FALSE")
  expect_error(stglmnet_options(return_fit = NA), "must be TRUE or FALSE")
  expect_error(stglmnet_options(cv_folds = 3.9), "must be a positive integer")
  expect_error(stglmnet_options(alpha = c(0.1, 0.2)), "alpha must be one")
  expect_error(stglmnet_options(lambda = -0.1), "lambda must contain")
  expect_error(fmrilss:::.stg_resolve_options(list(typo = 1)), "unknown option")
  expect_error(
    fmrilss:::.stg_resolve_options(list(pool_to_mean = 1)),
    "must be TRUE or FALSE"
  )
})

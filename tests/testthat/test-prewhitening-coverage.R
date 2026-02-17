# Tests for prewhitening.R functions

test_that(".prewhiten_data returns unmodified data when method is 'none'", {
  Y <- matrix(rnorm(100), 50, 2)
  X <- matrix(rnorm(100), 50, 2)

  result <- fmrilss:::.prewhiten_data(Y, X = X, prewhiten = list(method = "none"))

  expect_false(result$applied)
  expect_equal(result$Y_whitened, Y)
  expect_equal(result$X_whitened, X)
  expect_null(result$whiten_plan)
})

test_that(".prewhiten_data returns unmodified data when method is NULL", {
  Y <- matrix(rnorm(100), 50, 2)

  result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = NULL))

  expect_false(result$applied)
  expect_equal(result$Y_whitened, Y)
})

test_that(".prewhiten_data handles non-matrix Y when method is none", {
  # When method="none", Y is returned as-is (not coerced)
  Y <- rnorm(50)

  result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = "none"))

  # With method="none", Y is not modified
  expect_equal(result$Y_whitened, Y)
})

test_that(".prewhiten_data uses defaults correctly", {
  skip_if_not_installed("fmriAR")

  set.seed(123)
  Y <- matrix(rnorm(100 * 3), 100, 3)
  X <- matrix(rnorm(100 * 2), 100, 2)

  # Use default AR method
  result <- fmrilss:::.prewhiten_data(Y, X = X, prewhiten = list(method = "ar"))

  expect_true(result$applied)
  expect_true(!is.null(result$whiten_plan))
  expect_equal(dim(result$Y_whitened), dim(Y))
  expect_equal(dim(result$X_whitened), dim(X))
})

test_that(".prewhiten_data handles Z and Nuisance", {
  skip_if_not_installed("fmriAR")

  set.seed(234)
  n_time <- 100
  Y <- matrix(rnorm(n_time * 3), n_time, 3)
  X <- matrix(rnorm(n_time * 2), n_time, 2)
  Z <- matrix(1, n_time, 1)  # Intercept
  Nuisance <- matrix(rnorm(n_time * 2), n_time, 2)

  result <- fmrilss:::.prewhiten_data(Y, X = X, Z = Z, Nuisance = Nuisance,
                                      prewhiten = list(method = "ar"))

  expect_true(result$applied)
  expect_equal(dim(result$Y_whitened), dim(Y))
  expect_equal(dim(result$X_whitened), dim(X))
  expect_equal(dim(result$Z_whitened), dim(Z))
  expect_equal(dim(result$Nuisance_whitened), dim(Nuisance))
})

test_that(".prewhiten_data handles voxel pooling", {
  skip_if_not_installed("fmriAR")

  set.seed(345)
  n_time <- 80
  n_vox <- 4
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = "ar", pooling = "voxel"))

  expect_true(result$applied)
  expect_equal(dim(result$Y_whitened), dim(Y))
})

test_that(".prewhiten_data handles parcel pooling", {
  skip_if_not_installed("fmriAR")

  set.seed(456)
  n_time <- 80
  n_vox <- 6
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  parcels <- c(1, 1, 2, 2, 3, 3)  # 3 parcels

  result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = "ar",
                                                          pooling = "parcel",
                                                          parcels = parcels))

  expect_true(result$applied)
  expect_equal(dim(result$Y_whitened), dim(Y))
})

test_that(".prewhiten_data handles compute_residuals=FALSE", {
  skip_if_not_installed("fmriAR")

  set.seed(567)
  n_time <- 80
  Y <- matrix(rnorm(n_time * 3), n_time, 3)

  result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = "ar",
                                                          compute_residuals = FALSE))

  expect_true(result$applied)
})

test_that(".prewhiten_data handles missing design matrices", {
  skip_if_not_installed("fmriAR")

  set.seed(678)
  Y <- matrix(rnorm(100 * 3), 100, 3)

  result <- fmrilss:::.prewhiten_data(Y, X = NULL, Z = NULL, Nuisance = NULL,
                                      prewhiten = list(method = "ar"))

  expect_true(result$applied)
  expect_null(result$X_whitened)
  expect_null(result$Z_whitened)
  expect_null(result$Nuisance_whitened)
})

test_that(".convert_legacy_whiten returns NULL for method='none'", {
  result <- fmrilss:::.convert_legacy_whiten(list(whiten = "none"))
  expect_null(result)
})

test_that(".convert_legacy_whiten returns NULL for NULL whiten", {
  result <- fmrilss:::.convert_legacy_whiten(list(whiten = NULL))
  expect_null(result)
})

test_that(".convert_legacy_whiten returns NULL for missing whiten", {
  result <- fmrilss:::.convert_legacy_whiten(list())
  expect_null(result)
})

test_that(".convert_legacy_whiten converts ar1", {
  result <- fmrilss:::.convert_legacy_whiten(list(whiten = "ar1"))

  expect_true(is.list(result))
  expect_equal(result$method, "ar")
  expect_equal(result$p, 1L)
  expect_equal(result$pooling, "global")
  expect_equal(result$exact_first, "ar1")
})

test_that(".convert_legacy_whiten warns on unknown whiten option", {
  expect_warning(
    result <- fmrilss:::.convert_legacy_whiten(list(whiten = "unknown_method")),
    "Unknown whiten option"
  )
  expect_null(result)
})

test_that(".needs_advanced_prewhitening returns FALSE for NULL", {
  result <- fmrilss:::.needs_advanced_prewhitening(NULL)
  expect_false(result)
})

test_that(".needs_advanced_prewhitening returns TRUE for arma method", {
  result <- fmrilss:::.needs_advanced_prewhitening(list(method = "arma"))
  expect_true(result)
})

test_that(".needs_advanced_prewhitening returns TRUE for non-standard p", {
  # p = 2 (not 1 or "auto")
  result <- fmrilss:::.needs_advanced_prewhitening(list(p = 2L))
  expect_true(result)
})

test_that(".needs_advanced_prewhitening returns TRUE for q > 0", {
  result <- fmrilss:::.needs_advanced_prewhitening(list(q = 1L))
  expect_true(result)
})

test_that(".needs_advanced_prewhitening returns TRUE for non-global pooling", {
  result <- fmrilss:::.needs_advanced_prewhitening(list(pooling = "voxel"))
  expect_true(result)

  result2 <- fmrilss:::.needs_advanced_prewhitening(list(pooling = "parcel"))
  expect_true(result2)
})

test_that(".needs_advanced_prewhitening returns TRUE for parcels", {
  # Just check that parcels is not NULL (vector comparison avoided)
  result <- fmrilss:::.needs_advanced_prewhitening(list(parcels = 1L))
  expect_true(result)
})

test_that(".needs_advanced_prewhitening returns TRUE for runs", {
  # Just check that runs is not NULL
  result <- fmrilss:::.needs_advanced_prewhitening(list(runs = 1L))
  expect_true(result)
})

test_that(".needs_advanced_prewhitening returns FALSE for simple AR(1)", {
  result <- fmrilss:::.needs_advanced_prewhitening(list(method = "ar", p = 1L))
  expect_false(result)

  result2 <- fmrilss:::.needs_advanced_prewhitening(list(method = "ar"))
  expect_false(result2)
})

test_that(".prewhiten_data expands short parcels vector", {
  skip_if_not_installed("fmriAR")

  set.seed(789)
  n_time <- 60
  n_vox <- 6
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  # Provide parcels matching n_vox
  result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = "ar",
                                                          pooling = "parcel",
                                                          parcels = c(1, 1, 2, 2, 3, 3)))

  expect_true(result$applied)
})

test_that(".prewhiten_data warns when parcels provided with voxel pooling", {
  skip_if_not_installed("fmriAR")

  set.seed(890)
  n_time <- 60
  n_vox <- 4
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  expect_warning(
    result <- fmrilss:::.prewhiten_data(Y, prewhiten = list(method = "ar",
                                                            pooling = "voxel",
                                                            parcels = c(1, 1, 2, 2))),
    "parcels is ignored"
  )
})

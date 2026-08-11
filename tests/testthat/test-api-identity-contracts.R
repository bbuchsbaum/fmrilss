test_that("lss and lsa enforce stable trial and voxel identities", {
  set.seed(901)
  Y <- matrix(rnorm(80), 40, 2)
  X <- matrix(rnorm(120), 40, 3)
  methods <- c("r_optimized", "cpp_optimized", "naive", "oasis")

  for (method in methods) {
    fit <- lss(Y, X, method = method)
    expect_identical(rownames(fit), paste0("Trial_", 1:3), info = method)
    expect_identical(colnames(fit), paste0("Voxel_", 1:2), info = method)

    X_bad <- X
    colnames(X_bad) <- c("same", "same", "third")
    expect_error(lss(Y, X_bad, method = method), "complete and unique", info = method)

    Y_bad <- Y
    colnames(Y_bad) <- c("vox", "vox")
    expect_error(lss(Y_bad, X, method = method), "complete and unique", info = method)
  }

  fit_lsa <- lsa(Y, X)
  expect_identical(rownames(fit_lsa), paste0("Trial_", 1:3))
  expect_identical(colnames(fit_lsa), paste0("Voxel_", 1:2))
  colnames(X) <- c("same", "same", "third")
  expect_error(lsa(Y, X), "complete and unique")
})

test_that("lss backends reject empty and non-finite scientific inputs", {
  set.seed(902)
  Y <- matrix(rnorm(80), 40, 2)
  X <- matrix(rnorm(120), 40, 3)
  methods <- c("r_optimized", "cpp_optimized", "naive", "oasis")

  for (method in methods) {
    Y_bad <- Y
    Y_bad[1, 1] <- NA_real_
    expect_error(lss(Y_bad, X, method = method), "Y contains non-finite", info = method)

    X_bad <- X
    X_bad[1, 1] <- Inf
    expect_error(lss(Y, X_bad, method = method), "X contains non-finite", info = method)

    expect_error(
      lss(matrix(numeric(), 0, 2), matrix(numeric(), 0, 3), method = method),
      "at least one timepoint", info = method
    )
    expect_error(
      lss(matrix(numeric(), 40, 0), matrix(numeric(), 40, 3), method = method),
      "at least one timepoint and one voxel", info = method
    )
    expect_error(
      lss(Y, matrix(numeric(), 40, 0), method = method),
      "at least one timepoint and one trial", info = method
    )
  }
})

test_that("identity follows trial and voxel permutations", {
  set.seed(903)
  Y <- matrix(rnorm(120), 40, 3, dimnames = list(NULL, c("v1", "v2", "v3")))
  X <- matrix(rnorm(120), 40, 3, dimnames = list(NULL, c("t1", "t2", "t3")))
  ref <- lss(Y, X, method = "oasis", oasis = list(
    ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
  ))
  perm <- lss(Y[, c(3, 1, 2)], X[, c(2, 3, 1)], method = "oasis", oasis = list(
    ridge_mode = "absolute", ridge_x = 0, ridge_b = 0
  ))
  expect_equal(perm[c("t1", "t2", "t3"), c("v1", "v2", "v3")], ref,
               tolerance = 1e-12)
})

test_that("OASIS return diagnostics and event-built maps have stable shapes", {
  skip_if_not_installed("fmrihrf")
  set.seed(904)
  n_time <- 100
  Y <- matrix(rnorm(n_time * 3), n_time, 3,
              dimnames = list(NULL, c("v1", "v2", "v3")))
  X <- matrix(rnorm(n_time * 4), n_time, 4,
              dimnames = list(NULL, paste0("t", 1:4)))
  one <- lss(Y, X, method = "oasis", oasis = list(return_diag = TRUE))
  expect_type(one$diag$d, "double")
  expect_null(dim(one$diag$d))
  expect_identical(names(one$diag$d), rownames(one$beta))
  expect_identical(names(one$diag$alpha), rownames(one$beta))
  expect_identical(names(one$diag$s), rownames(one$beta))

  spec <- list(
    sframe = fmrihrf::sampling_frame(blocklens = n_time, TR = 1),
    cond = list(
      onsets = c(10, 30, 50),
      hrf = fmrihrf::HRF_SPMG3,
      span = 25
    )
  )
  multi <- lss(Y, NULL, method = "oasis", oasis = list(
    design_spec = spec,
    ridge_mode = "absolute", ridge_x = 0, ridge_b = 0,
    return_se = TRUE, return_diag = TRUE
  ))
  beta_map <- attr(multi$beta, "trial_basis_map")
  se_map <- attr(multi$se, "trial_basis_map")
  expect_s3_class(beta_map, "data.frame")
  expect_identical(se_map, beta_map)
  expect_identical(beta_map$output_name, rownames(multi$beta))
  expect_identical(dim(multi$diag$D), c(3L, 3L, 3L))
  expect_identical(dim(multi$diag$C), c(3L, 3L, 3L))
  expect_identical(dim(multi$diag$E), c(3L, 3L, 3L))
})

test_that("HRF grids are accepted only inside design_spec", {
  skip_if_not_installed("fmrihrf")
  set.seed(905)
  n_time <- 80
  Y <- matrix(rnorm(n_time * 2), n_time, 2)
  grid <- list(fmrihrf::HRF_SPMG1, fmrihrf::HRF_GAUSSIAN)
  spec <- list(
    sframe = fmrihrf::sampling_frame(blocklens = n_time, TR = 1),
    cond = list(onsets = c(10, 30, 50), hrf = fmrihrf::HRF_SPMG1, span = 20)
  )
  expect_error(
    lss(Y, NULL, method = "oasis", oasis = list(design_spec = spec, hrf_grid = grid)),
    "unknown option"
  )
  spec$hrf_grid <- grid
  expect_no_error(lss(Y, NULL, method = "oasis", oasis = list(design_spec = spec)))
})

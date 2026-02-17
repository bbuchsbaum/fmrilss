test_that("item_compute_u returns symmetric matrix with expected dimensions", {
  set.seed(1001)

  X_t <- matrix(rnorm(30 * 6), nrow = 30, ncol = 6)
  U <- item_compute_u(X_t, ridge = 1e-8)

  expect_equal(dim(U), c(6, 6))
  expect_equal(U, t(U), tolerance = 1e-8)

  expected <- solve(crossprod(X_t) + diag(1e-8, 6))
  U_plain <- U
  attr(U_plain, "item_diagnostics") <- NULL
  expect_equal(unname(U_plain), unname(expected), tolerance = 1e-6)

  eigvals <- eigen(U, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(eigvals) > -1e-8)
})

test_that("item_compute_u supports run-block V list and matches dense equivalent", {
  set.seed(1002)

  X_t <- matrix(rnorm(24 * 5), nrow = 24, ncol = 5)

  V_blocks <- list(
    diag(12) * 2,
    diag(12) * 3
  )

  V_dense <- matrix(0, nrow = 24, ncol = 24)
  V_dense[1:12, 1:12] <- V_blocks[[1]]
  V_dense[13:24, 13:24] <- V_blocks[[2]]

  U_block <- item_compute_u(X_t, V = V_blocks, v_type = "cov", ridge = 1e-6)
  U_dense <- item_compute_u(X_t, V = V_dense, v_type = "cov", ridge = 1e-6)

  expect_equal(U_block, U_dense, tolerance = 1e-6)
})

test_that("item_slice_fold produces consistent train/test slicing", {
  set.seed(1003)

  n_time <- 60
  n_trials <- 12
  n_vox <- 4

  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  run_id <- rep(1:3, each = 4)

  bundle <- item_build_design(
    X_t = X_t,
    T_target = rep(c("A", "B"), length.out = n_trials),
    run_id = run_id
  )

  bundle$Gamma <- matrix(rnorm(n_trials * n_vox), nrow = n_trials, ncol = n_vox)
  bundle$U <- item_compute_u(X_t, ridge = 1e-6)

  fold <- item_slice_fold(bundle, test_run = 2)
  train_idx <- which(run_id != 2)
  test_idx <- which(run_id == 2)

  expect_equal(fold$train_idx, train_idx)
  expect_equal(fold$test_idx, test_idx)
  expect_equal(fold$Gamma_train, bundle$Gamma[train_idx, , drop = FALSE])
  expect_equal(fold$Gamma_test, bundle$Gamma[test_idx, , drop = FALSE])
  expect_equal(fold$U_train, bundle$U[train_idx, train_idx, drop = FALSE])
  expect_equal(fold$U_test, bundle$U[test_idx, test_idx, drop = FALSE])

  bundle_blocks <- bundle
  bundle_blocks$U <- NULL
  bundle_blocks$U_by_run <- item_compute_u(
    X_t,
    ridge = 1e-6,
    run_id = run_id,
    output = "by_run"
  )

  fold_blocks <- item_slice_fold(bundle_blocks, test_run = 2)
  expect_equal(fold_blocks$U_test, bundle_blocks$U_by_run[["2"]])
  expect_equal(nrow(fold_blocks$U_train), length(train_idx))
  expect_equal(ncol(fold_blocks$U_train), length(train_idx))
})

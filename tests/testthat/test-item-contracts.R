test_that("item_from_lsa returns aligned bundle and supports Z/Nuisance paths", {
  set.seed(4001)

  n_time <- 90
  n_trials <- 12
  n_vox <- 5

  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  Z <- cbind(1, scale(seq_len(n_time)))
  Z_alt <- cbind(1, scale((seq_len(n_time))^2))

  beta <- matrix(rnorm(n_trials * n_vox), nrow = n_trials, ncol = n_vox)
  beta_z <- matrix(rnorm(ncol(Z) * n_vox), nrow = ncol(Z), ncol = n_vox)
  Y <- X_t %*% beta + Z %*% beta_z + matrix(rnorm(n_time * n_vox, sd = 0.05), nrow = n_time)

  run_id <- rep(1:3, each = 4)
  T_target <- rep(c("A", "B"), length.out = n_trials)

  bundle_z <- item_from_lsa(
    Y = Y,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    Z = Z,
    lsa_method = "r",
    solver = "svd"
  )

  bundle_n <- item_from_lsa(
    Y = Y,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    Nuisance = Z,
    lsa_method = "r",
    solver = "svd"
  )

  bundle_both <- item_from_lsa(
    Y = Y,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    Z = Z,
    Nuisance = Z_alt,
    lsa_method = "r",
    solver = "svd"
  )

  expect_s3_class(bundle_z, "item_bundle")
  expect_equal(dim(bundle_z$Gamma), c(n_trials, n_vox))
  expect_equal(dim(bundle_z$U), c(n_trials, n_trials))
  expect_equal(bundle_z$run_id, run_id)

  expect_equal(bundle_z$Gamma, bundle_n$Gamma, tolerance = 1e-12)
  expect_equal(bundle_z$Gamma, bundle_both$Gamma, tolerance = 1e-12)
})

test_that("item_cv enforces trial hash when check_hash is TRUE", {
  set.seed(4002)

  n_trials <- 15
  n_features <- 4

  X_t <- matrix(rnorm(75 * n_trials), nrow = 75, ncol = n_trials)
  run_id <- rep(1:3, each = 5)
  trial_id <- sprintf("trial_%02d", seq_len(n_trials))
  trial_hash <- fmrilss:::.item_simple_hash(trial_id)

  bundle <- item_build_design(
    X_t = X_t,
    T_target = matrix(rnorm(n_trials), nrow = n_trials, ncol = 1),
    run_id = run_id,
    trial_id = trial_id,
    trial_hash = trial_hash
  )

  bundle$Gamma <- matrix(rnorm(n_trials * n_features), nrow = n_trials, ncol = n_features)
  bundle$U <- diag(n_trials)

  res <- item_cv(bundle, mode = "regression", metric = "correlation", check_hash = TRUE)
  expect_s3_class(res, "item_cv_result")

  bundle_bad <- bundle
  bundle_bad$trial_id[1] <- "tampered_trial"

  expect_error(
    item_cv(bundle_bad, mode = "regression", metric = "correlation", check_hash = TRUE),
    "Trial hash mismatch"
  )
})

test_that("item_cv is deterministic for fixed inputs and folds", {
  set.seed(4003)

  n_runs <- 4
  trials_per_run <- 20
  n_trials <- n_runs * trials_per_run
  n_features <- 6

  run_id <- rep(seq_len(n_runs), each = trials_per_run)
  Gamma <- matrix(rnorm(n_trials * n_features), nrow = n_trials, ncol = n_features)
  T_target <- matrix(rnorm(n_trials * 2), nrow = n_trials, ncol = 2)
  U <- diag(n_trials)

  res1 <- item_cv(
    Gamma = Gamma,
    T_target = T_target,
    U = U,
    run_id = run_id,
    mode = "regression",
    metric = "correlation",
    ridge = 1e-4,
    method = "svd"
  )

  res2 <- item_cv(
    Gamma = Gamma,
    T_target = T_target,
    U = U,
    run_id = run_id,
    mode = "regression",
    metric = "correlation",
    ridge = 1e-4,
    method = "svd"
  )

  expect_equal(res1$folds, res2$folds)
  expect_equal(res1$aggregate, res2$aggregate)
  expect_equal(res1$predictions$T_hat, res2$predictions$T_hat)
})

test_that("classification tie handling uses ties.method='first' and class order", {
  T_true <- rbind(
    c(0, 1, 0),
    c(1, 0, 0)
  )

  T_hat <- rbind(
    c(0.1, 0.9, 0.9),
    c(0.5, 0.5, 0.5)
  )

  scored <- fmrilss:::.item_score_fold(
    T_true = T_true,
    T_hat = T_hat,
    mode = "classification",
    metric = "accuracy",
    class_levels = c("A", "B", "C")
  )

  expect_equal(scored$pred_labels, c("B", "A"))
  expect_equal(scored$metric, 1)

  scored_reordered <- fmrilss:::.item_score_fold(
    T_true = T_true,
    T_hat = T_hat,
    mode = "classification",
    metric = "accuracy",
    class_levels = c("C", "B", "A")
  )

  expect_equal(scored_reordered$pred_labels, c("B", "C"))
})

test_that("ITEM helpers throw actionable errors on malformed dimensions/blocks", {
  set.seed(4004)

  Gamma <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  T_target <- matrix(rnorm(10), nrow = 10, ncol = 1)
  U <- diag(10)

  expect_error(
    item_cv(Gamma = Gamma, T_target = T_target, U = U, run_id = rep(1:3, length.out = 9), mode = "regression"),
    "run_id must have length 10"
  )

  X_t <- matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
  expect_error(
    item_compute_u(X_t = X_t, V = diag(19), v_type = "cov"),
    "V must be 20 x 20"
  )

  bundle <- item_build_design(
    X_t = matrix(rnorm(60 * 6), nrow = 60, ncol = 6),
    T_target = matrix(rnorm(6), nrow = 6, ncol = 1),
    run_id = rep(1:3, each = 2)
  )
  bundle$Gamma <- matrix(rnorm(6 * 2), nrow = 6, ncol = 2)
  bundle$U <- NULL

  bundle$U_by_run <- list("1" = diag(2), "2" = diag(2))
  expect_error(item_slice_fold(bundle, test_run = 1), "Missing U block for run '3'")

  bundle$U_by_run <- list("1" = diag(3), "2" = diag(2), "3" = diag(1))
  expect_error(item_slice_fold(bundle, test_run = 1), "must be 2 x 2")
})

test_that("balanced_accuracy and rmse metrics behave as expected", {
  T_true_cls <- rbind(
    c(1, 0, 0),
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )

  T_hat_cls <- rbind(
    c(0.9, 0.1, 0.0),
    c(0.8, 0.2, 0.0),
    c(0.7, 0.2, 0.1),
    c(0.6, 0.2, 0.2)
  )

  score_acc <- fmrilss:::.item_score_fold(
    T_true = T_true_cls,
    T_hat = T_hat_cls,
    mode = "classification",
    metric = "accuracy",
    class_levels = c("A", "B", "C")
  )

  score_bal <- fmrilss:::.item_score_fold(
    T_true = T_true_cls,
    T_hat = T_hat_cls,
    mode = "classification",
    metric = "balanced_accuracy",
    class_levels = c("A", "B", "C")
  )

  expect_equal(score_acc$metric, 0.5)
  expect_equal(score_bal$metric, 1 / 3, tolerance = 1e-12)

  T_true_reg <- matrix(c(0, 1, 2, 3), ncol = 1)
  T_hat_reg <- matrix(c(0, 2, 2, 4), ncol = 1)

  score_rmse <- fmrilss:::.item_score_fold(
    T_true = T_true_reg,
    T_hat = T_hat_reg,
    mode = "regression",
    metric = "rmse"
  )

  expect_equal(score_rmse$metric, sqrt(mean(c(0, 1, 0, 1))), tolerance = 1e-12)
})

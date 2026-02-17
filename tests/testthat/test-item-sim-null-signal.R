test_that("item_cv classification is near chance under null", {
  set.seed(3001)

  n_runs <- 4
  trials_per_run <- 45
  n_trials <- n_runs * trials_per_run
  n_features <- 8

  run_id <- rep(seq_len(n_runs), each = trials_per_run)
  Gamma <- matrix(rnorm(n_trials * n_features), nrow = n_trials, ncol = n_features)
  labels <- sample(c("A", "B", "C"), size = n_trials, replace = TRUE)
  U <- diag(n_trials)

  res <- item_cv(
    Gamma = Gamma,
    T_target = labels,
    U = U,
    run_id = run_id,
    mode = "classification",
    metric = "accuracy",
    ridge = 1e-4,
    method = "svd"
  )

  acc <- res$aggregate$mean
  chance <- 1 / 3

  expect_gt(acc, chance - 0.15)
  expect_lt(acc, chance + 0.15)
})

test_that("item_cv regression correlation is near zero under null", {
  set.seed(3002)

  n_runs <- 4
  trials_per_run <- 40
  n_trials <- n_runs * trials_per_run
  n_features <- 7

  run_id <- rep(seq_len(n_runs), each = trials_per_run)
  Gamma <- matrix(rnorm(n_trials * n_features), nrow = n_trials, ncol = n_features)
  T_target <- matrix(rnorm(n_trials * 2), nrow = n_trials, ncol = 2)
  U <- diag(n_trials)

  res <- item_cv(
    Gamma = Gamma,
    T_target = T_target,
    U = U,
    run_id = run_id,
    mode = "regression",
    metric = "correlation",
    ridge = 1e-4,
    method = "svd"
  )

  expect_lt(abs(res$aggregate$mean), 0.2)
})

test_that("item_cv regression performance improves with SNR", {
  set.seed(3003)

  n_runs <- 5
  trials_per_run <- 36
  n_trials <- n_runs * trials_per_run
  n_features <- 6
  n_targets <- 2

  run_id <- rep(seq_len(n_runs), each = trials_per_run)
  Gamma <- matrix(rnorm(n_trials * n_features), nrow = n_trials, ncol = n_features)
  W_true <- matrix(rnorm(n_features * n_targets), nrow = n_features, ncol = n_targets)

  signal <- Gamma %*% W_true
  T_low_snr <- signal + matrix(rnorm(n_trials * n_targets, sd = 2.0), nrow = n_trials)
  T_high_snr <- signal + matrix(rnorm(n_trials * n_targets, sd = 0.4), nrow = n_trials)

  U <- diag(n_trials)

  res_low <- item_cv(
    Gamma = Gamma,
    T_target = T_low_snr,
    U = U,
    run_id = run_id,
    mode = "regression",
    metric = "correlation",
    ridge = 1e-4,
    method = "svd"
  )

  res_high <- item_cv(
    Gamma = Gamma,
    T_target = T_high_snr,
    U = U,
    run_id = run_id,
    mode = "regression",
    metric = "correlation",
    ridge = 1e-4,
    method = "svd"
  )

  expect_gt(res_high$aggregate$mean, res_low$aggregate$mean)
})

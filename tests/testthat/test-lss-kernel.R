make_kernel_problem <- function(n = 160, n_trials = 24, V = 12, seed = 11) {
  set.seed(seed)
  onsets <- sort(sample(3:(n - 18), n_trials))
  h <- dgamma(0:15, 6, 1) - dgamma(0:15, 16, 1) / 6
  X <- vapply(onsets, function(o) {
    s <- numeric(n)
    s[o] <- 1
    out <- stats::filter(s, h, sides = 1)
    out[is.na(out)] <- 0
    as.numeric(out)
  }, numeric(n))
  Z <- cbind(1, seq_len(n) / n)
  N <- matrix(rnorm(n * 3), n, 3)
  noise <- apply(matrix(rnorm(n * V), n, V), 2, function(e) {
    as.numeric(stats::filter(e, 0.4, "recursive"))
  })
  B <- matrix(rnorm(n_trials * V), n_trials, V)
  list(Y = X %*% B + noise, X = X, Z = Z, N = N, n = n, T = n_trials, V = V)
}

test_that("weight-matrix kernels match the naive per-trial GLM", {
  p <- make_kernel_problem()
  ref <- lss(p$Y, p$X, p$Z, p$N, method = "naive")
  for (m in c("r_optimized", "cpp_optimized", "cpp")) {
    expect_equal(unname(lss(p$Y, p$X, p$Z, p$N, method = m)), unname(ref),
                 tolerance = 1e-10, info = m)
  }
  # Without separate nuisance regressors
  ref0 <- lss(p$Y, p$X, p$Z, method = "naive")
  for (m in c("r_optimized", "cpp_optimized", "cpp")) {
    expect_equal(unname(lss(p$Y, p$X, p$Z, method = m)), unname(ref0),
                 tolerance = 1e-10, info = m)
  }
})

test_that("cpp_optimized gives identical results with and without OpenMP blocks", {
  p <- make_kernel_problem(V = 50)
  Zb <- qr.Q(qr(p$Z))
  a <- lss_fused_optim_cpp(Zb, p$Y, p$X, block_size = 7L, use_omp = TRUE)
  b <- lss_fused_optim_cpp(Zb, p$Y, p$X, block_size = 7L, use_omp = FALSE)
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("LSS-N trial groups match the naive reference", {
  p <- make_kernel_problem()
  g <- rep(c("faces", "houses", "objects"), length.out = p$T)
  ref <- lss(p$Y, p$X, p$Z, p$N, method = "naive", trial_groups = g)
  for (m in c("r_optimized", "cpp_optimized", "cpp")) {
    expect_equal(unname(lss(p$Y, p$X, p$Z, p$N, method = m, trial_groups = g)),
                 unname(ref), tolerance = 1e-10, info = m)
  }
  # Factor and integer labels are equivalent to character labels
  expect_equal(lss(p$Y, p$X, p$Z, p$N, trial_groups = factor(g)),
               lss(p$Y, p$X, p$Z, p$N, trial_groups = g))
  expect_equal(lss(p$Y, p$X, p$Z, p$N, trial_groups = match(g, unique(g))),
               lss(p$Y, p$X, p$Z, p$N, trial_groups = g))
  # A single group reduces to classic LSS
  expect_equal(lss(p$Y, p$X, p$Z, p$N, trial_groups = rep("a", p$T)),
               lss(p$Y, p$X, p$Z, p$N))
  # Grouping changes the model when conditions are separated
  expect_gt(max(abs(ref - lss(p$Y, p$X, p$Z, p$N, method = "naive"))), 1e-6)
})

test_that("LSS-N handles a condition with a single trial", {
  p <- make_kernel_problem()
  g <- c("solo", rep("rest", p$T - 1))
  ref <- lss(p$Y, p$X, p$Z, method = "naive", trial_groups = g)
  expect_equal(unname(lss(p$Y, p$X, p$Z, trial_groups = g)), unname(ref),
               tolerance = 1e-10)
  expect_equal(unname(lss(p$Y, p$X, p$Z, method = "cpp_optimized", trial_groups = g)),
               unname(ref), tolerance = 1e-10)
})

test_that("trial_groups is validated", {
  p <- make_kernel_problem()
  expect_error(lss(p$Y, p$X, trial_groups = c("a", "b")), "one label per trial")
  expect_error(lss(p$Y, p$X, trial_groups = c(NA, rep("a", p$T - 1))), "must not contain NA")
  expect_error(lss(p$Y, p$X, method = "r_vectorized", trial_groups = rep("a", p$T)),
               "trial_groups")
})

test_that("non-finite data are rejected without allocating a full mask", {
  p <- make_kernel_problem()
  Y <- p$Y
  Y[3, 2] <- Inf
  expect_error(lss(Y, p$X), "non-finite")
  Y[3, 2] <- NA
  expect_error(lss(Y, p$X), "non-finite")
  expect_true(fmrilss:::.all_finite(p$Y))
  expect_false(fmrilss:::.all_finite(c(1, -Inf)))
})

test_that("prewhitened weight-matrix LSS matches the generic whitening path", {
  skip_if_not_installed("fmriAR")
  p <- make_kernel_problem()
  runs <- rep(1:2, each = p$n / 2)
  g <- rep(c("a", "b"), length.out = p$T)
  specs <- list(
    list(method = "ar", p = 1),
    list(method = "ar", p = 2, runs = runs),
    list(method = "ar", p = 1, pooling = "run", runs = runs)
  )
  for (pw in specs) {
    ref <- lss(p$Y, p$X, p$Z, p$N, method = "naive", prewhiten = pw)
    for (m in c("r_optimized", "cpp_optimized", "cpp")) {
      est <- lss(p$Y, p$X, p$Z, p$N, method = m, prewhiten = pw)
      expect_equal(unname(est), unname(ref), tolerance = 1e-10,
                   ignore_attr = TRUE, info = m)
      expect_s3_class(attr(est, "whiten_plan"), "fmriAR_plan")
    }
  }
  ref <- lss(p$Y, p$X, p$Z, p$N, method = "naive", prewhiten = specs[[1]],
             trial_groups = g)
  est <- lss(p$Y, p$X, p$Z, p$N, prewhiten = specs[[1]], trial_groups = g)
  expect_equal(unname(est), unname(ref), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("compressed noise fit reproduces the full-residual AR plan", {
  skip_if_not_installed("fmriAR")
  set.seed(5)
  n <- 90
  V <- 400
  runs <- rep(1:3, each = 30)
  D <- cbind(1, rnorm(n))
  Y <- apply(matrix(rnorm(n * V), n, V), 2, function(e) {
    as.numeric(stats::filter(e, c(0.4, 0.15), "recursive"))
  }) + 5
  specs <- list(
    list(method = "ar", p = 1),
    list(method = "ar", p = "auto"),
    list(method = "ar", p = 2, pooling = "run", runs = runs)
  )
  for (pw in specs) {
    opts <- fmrilss:::.resolve_prewhiten_options(pw)
    full <- fmrilss:::.fit_whiten_plan(fmrilss:::.noise_residuals(Y, D), opts)
    fast <- fmrilss:::.fit_noise_plan(Y, D, opts)
    expect_equal(unlist(fast$phi), unlist(full$phi), tolerance = 1e-10)
    expect_equal(unlist(fast$sigma2), unlist(full$sigma2), tolerance = 1e-10)
    expect_equal(fast$order, full$order)
  }
})

test_that("voxel pooling with one bin per voxel equals per-voxel GLS", {
  skip_if_not_installed("fmriAR")
  p <- make_kernel_problem(V = 8)
  est <- lss(p$Y, p$X, p$Z, p$N,
             prewhiten = list(method = "ar", p = 1, pooling = "voxel", voxel_bins = p$V))
  plan <- attr(est, "whiten_plan")
  expect_equal(length(unique(plan$parcels)), p$V)
  ref <- vapply(seq_len(p$V), function(j) {
    lss(p$Y[, j, drop = FALSE], p$X, p$Z, p$N, method = "naive",
        prewhiten = list(method = "ar", p = 1))[, 1]
  }, numeric(p$T))
  expect_equal(unname(est), unname(ref), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("voxel pooling bins voxels by autocorrelation", {
  skip_if_not_installed("fmriAR")
  set.seed(9)
  p <- make_kernel_problem(V = 60)
  phi <- rep(c(0, 0.6), each = 30)
  noise <- vapply(phi, function(f) {
    as.numeric(stats::filter(rnorm(p$n), f, "recursive"))
  }, numeric(p$n))
  Y <- p$X %*% matrix(rnorm(p$T * 60), p$T, 60) + noise
  est <- lss(Y, p$X, p$Z,
             prewhiten = list(method = "ar", p = 1, pooling = "voxel", voxel_bins = 2))
  plan <- attr(est, "whiten_plan")
  bins <- plan$parcels
  expect_equal(length(unique(bins)), 2L)
  # The two simulated populations land in different bins
  expect_gt(mean(bins[1:30] != bins[31:60][1]), 0.8)
  phis <- vapply(plan$phi_by_parcel, function(x) x[1], numeric(1))
  expect_gt(diff(range(phis)), 0.3)
  for (m in c("cpp_optimized", "cpp")) {
    expect_equal(lss(Y, p$X, p$Z, method = m,
                     prewhiten = list(method = "ar", p = 1, pooling = "voxel",
                                      voxel_bins = 2)),
                 est, tolerance = 1e-10, ignore_attr = TRUE, info = m)
  }
})

test_that("parcel pooling fits each parcel with its own filtered design", {
  skip_if_not_installed("fmriAR")
  p <- make_kernel_problem(V = 10)
  parcels <- rep(1:2, each = 5)
  pw <- list(method = "ar", p = 1, pooling = "parcel", parcels = parcels)
  est <- lss(p$Y, p$X, p$Z, p$N, prewhiten = pw)
  plan <- attr(est, "whiten_plan")
  K <- cbind(p$Z, p$N)
  wh <- fmriAR::whiten_apply(plan, X = cbind(p$X, K), Y = p$Y, parcels = parcels)
  for (key in names(wh$X_by)) {
    cols <- which(parcels == as.integer(key))
    D <- wh$X_by[[key]]
    ref <- lss(wh$Y[, cols, drop = FALSE], D[, seq_len(p$T)],
               Z = D[, -seq_len(p$T)], method = "naive")
    expect_equal(unname(est[, cols]), unname(ref), tolerance = 1e-10)
  }
})

test_that("voxel and parcel pooling are rejected by generic-whitening methods", {
  skip_if_not_installed("fmriAR")
  p <- make_kernel_problem()
  expect_error(
    lss(p$Y, p$X, method = "naive",
        prewhiten = list(method = "ar", p = 1, pooling = "voxel")),
    "cannot be applied to a shared design matrix"
  )
})

test_that("aggregate noise residuals avoid the many-trial AR bias", {
  skip_if_not_installed("fmriAR")
  set.seed(21)
  n <- 240
  n_trials <- 90
  V <- 200
  onsets <- sort(sample(3:(n - 18), n_trials))
  h <- dgamma(0:15, 6, 1) - dgamma(0:15, 16, 1) / 6
  X <- vapply(onsets, function(o) {
    s <- numeric(n)
    s[o] <- 1
    out <- stats::filter(s, h, sides = 1)
    out[is.na(out)] <- 0
    as.numeric(out)
  }, numeric(n))
  noise <- apply(matrix(rnorm(n * V), n, V), 2, function(e) {
    as.numeric(stats::filter(e, 0.4, "recursive"))
  })
  Y <- X %*% matrix(rnorm(n_trials * V, 1, 0.3), n_trials, V) + noise

  agg <- lss(Y, X, prewhiten = list(method = "ar", p = 1))
  full <- lss(Y, X, prewhiten = list(method = "ar", p = 1, residual_model = "full"))
  phi_agg <- attr(agg, "whiten_plan")$phi[[1]]
  phi_full <- attr(full, "whiten_plan")$phi[[1]]
  expect_lt(abs(phi_agg - 0.4), 0.1)
  # 90 trial columns in 240 scans bias the full-model residual AR estimate
  expect_lt(phi_full, phi_agg - 0.15)

  # The generic whitening path uses the same noise model
  ref <- lss(Y, X, method = "naive", prewhiten = list(method = "ar", p = 1))
  expect_equal(attr(ref, "whiten_plan")$phi[[1]], phi_agg, tolerance = 1e-10)
})

test_that("residual_model is validated against the bias correction", {
  expect_identical(
    fmrilss:::.resolve_prewhiten_options(list(method = "ar"))$residual_model,
    "aggregate"
  )
  expect_identical(
    fmrilss:::.resolve_prewhiten_options(list(method = "ar", design = diag(3)))$residual_model,
    "full"
  )
  expect_error(
    fmrilss:::.resolve_prewhiten_options(
      list(method = "ar", design = diag(3), residual_model = "aggregate")
    ),
    "residual_model = 'full'"
  )
  expect_error(prewhiten_options(method = "ar", residual_model = "bogus"))
})

test_that("fractional ridge matches the naive reference and OASIS", {
  p <- make_kernel_problem()
  g <- rep(c("a", "b", "c"), length.out = p$T)
  for (rg in list(0.3, c(0.1, 0.5), c(0.5, 0))) {
    ref <- lss(p$Y, p$X, p$Z, p$N, method = "naive", ridge = rg)
    refg <- lss(p$Y, p$X, p$Z, p$N, method = "naive", ridge = rg, trial_groups = g)
    for (m in c("r_optimized", "cpp_optimized", "cpp")) {
      expect_equal(unname(lss(p$Y, p$X, p$Z, p$N, method = m, ridge = rg)),
                   unname(ref), tolerance = 1e-10, info = m)
      expect_equal(unname(lss(p$Y, p$X, p$Z, p$N, method = m, ridge = rg,
                              trial_groups = g)),
                   unname(refg), tolerance = 1e-10, info = m)
    }
    rg2 <- rep_len(rg, 2L)
    oasis <- lss(p$Y, p$X, p$Z, p$N, method = "oasis",
                 oasis = list(ridge_x = rg2[1], ridge_b = rg2[2],
                              ridge_mode = "fractional"))
    expect_equal(unname(oasis), unname(ref), tolerance = 1e-10)
  }
  expect_equal(lss(p$Y, p$X, p$Z, ridge = 0), lss(p$Y, p$X, p$Z))
})

test_that("ridge composes with prewhitening", {
  skip_if_not_installed("fmriAR")
  p <- make_kernel_problem()
  pw <- list(method = "ar", p = 1)
  ref <- lss(p$Y, p$X, p$Z, p$N, method = "naive", prewhiten = pw, ridge = c(0.4, 0.1))
  for (m in c("r_optimized", "cpp_optimized", "cpp")) {
    expect_equal(unname(lss(p$Y, p$X, p$Z, p$N, method = m, prewhiten = pw,
                            ridge = c(0.4, 0.1))),
                 unname(ref), tolerance = 1e-10, ignore_attr = TRUE, info = m)
  }
})

test_that("ridge reduces error for overlapping rapid-design trials", {
  set.seed(31)
  n <- 400
  n_trials <- 100
  onsets <- cumsum(c(20, sample(2:4, n_trials - 1, replace = TRUE)))
  h <- dgamma(0:15, 6, 1) - dgamma(0:15, 16, 1) / 6
  X <- vapply(onsets, function(o) {
    out <- numeric(n)
    idx <- o + seq_along(h) - 1L
    keep <- idx <= n
    out[idx[keep]] <- h[keep]
    out
  }, numeric(n))
  expect_gt(min(colSums(X^2)), 0.9 * sum(h^2))
  V <- 50
  B <- matrix(rnorm(n_trials * V, 1, 0.5), n_trials, V)
  Y <- X %*% B + matrix(rnorm(n * V), n, V)
  rmse <- function(est) sqrt(mean((est - B)^2))
  expect_lt(rmse(lss(Y, X, ridge = c(0.5, 0))), rmse(lss(Y, X)))
})

test_that("ridge is validated", {
  p <- make_kernel_problem()
  expect_error(lss(p$Y, p$X, ridge = -1), "nonnegative")
  expect_error(lss(p$Y, p$X, ridge = c(1, 2, 3)), "one or two")
  expect_error(lss(p$Y, p$X, method = "oasis", ridge = 1), "oasis\\$ridge_x")
})

test_that("residual_model = 'corrected' equals a manually supplied correction design", {
  skip_if_not_installed("fmriAR")
  set.seed(4)
  n <- 160
  V <- 400
  X <- matrix(rnorm(n * 40), n, 40)
  Z <- cbind(1, seq_len(n) / n)
  Y <- apply(matrix(rnorm(n * V), n, V), 2, function(e) {
    as.numeric(stats::filter(e, 0.4, "recursive"))
  })
  auto <- lss(Y, X, Z, prewhiten = list(method = "ar", p = 1,
                                        residual_model = "corrected"))
  manual <- lss(Y, X, Z, prewhiten = list(method = "ar", p = 1,
                                          design = cbind(Z, X)))
  expect_equal(attr(auto, "whiten_plan")$phi, attr(manual, "whiten_plan")$phi,
               tolerance = 1e-10)
  expect_equal(unname(auto), unname(manual), tolerance = 1e-10, ignore_attr = TRUE)
  naive <- lss(Y[, 1:50], X, Z, method = "naive",
               prewhiten = list(method = "ar", p = 1, residual_model = "corrected"))
  fast <- lss(Y[, 1:50], X, Z,
              prewhiten = list(method = "ar", p = 1, residual_model = "corrected"))
  expect_equal(unname(fast), unname(naive), tolerance = 1e-10, ignore_attr = TRUE)
  expect_error(
    lss(Y, X, prewhiten = list(method = "ar", p = 1, residual_model = "corrected",
                               pooling = "voxel")),
    "global' or 'run"
  )
  expect_error(
    prewhiten_options(method = "ar", residual_model = "corrected", design = cbind(1, X)),
    "do not also supply"
  )
})

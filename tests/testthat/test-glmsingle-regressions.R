test_that("fraction interpolation handles zero and nonfinite coefficient norms", {
  grid <- c(0, 1, 3)
  norms <- cbind(c(2, 1, .5), c(0, 0, 0), c(NaN, NaN, NaN), c(Inf, 1, 0), c(2, NaN, 0))
  for (nt in c(1L, 2L)) {
    alpha <- glms_frac_alpha_grid(norms, grid, c(1, .5, .25), nt)
    expect_equal(alpha[, 1], grid, tolerance = 1e-14)
    expect_true(all(is.na(alpha[, -1])))
    expect_error(glms_frac_alpha_grid(matrix(0, 0, 1), numeric(), .5, nt), "nonempty")
  }
  # Equal interpolation knots follow numpy.interp's last-match rule.
  expect_equal(drop(glms_frac_alpha_grid(matrix(c(2, 1, 1), 3), grid, .5)), 1)
})

test_that("zero eigenvalues do not contaminate fractional ridge penalties", {
  sp <- list(a_all = matrix(c(2, 0), 2), s2_all = c(1, 0))
  for (method in c("grid", "exact")) {
    alpha <- .glms_frac_alphas(sp, c(1, .5), method)
    expect_equal(alpha[, 1], c(0, 1), tolerance = 1e-10)
    zero <- list(a_all = matrix(0, 2, 1), s2_all = c(0, 0))
    expect_true(all(is.na(.glms_frac_alphas(zero, c(1, .5), method))))
  }
})

test_that("singular public fits still shrink by the requested ridge fraction", {
  sim <- sim_glms(seed = 4, n_runs = 2, n_time = 60, n_vox = 5)
  extras <- lapply(sim$design, function(D) {
    on <- which(rowSums(D) > 0)
    E <- matrix(0, nrow(D), 2)
    E[on[1], 1] <- E[on[2], 2] <- 1
    E
  })
  for (method in c("grid", "exact")) {
    expect_warning(fit <- glmsingle(
      sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
      hrf_library = matrix(c(1, 0, 0), ncol = 1), nuisance = extras,
      n_pcs = 0, noise_pool_r2 = 2, pc_voxel_r2 = 2,
      singular = "pinv", ridge_fracs = .5, ridge_alpha = method,
      ridge_rescale = FALSE, percent_signal = FALSE, verbose = FALSE
    ), "minimum-norm")
    norm_ratio <- sqrt(colSums(fit$typed$betasmd^2) / colSums(fit$typec$betasmd^2))
    expect_equal(norm_ratio, rep(.5, 5), tolerance = if (method == "exact") 1e-10 else .005)
    expect_equal(fit$typed$FRACvalue, rep(.5, 5))
  }
})

test_that("event run IDs follow the order of the data, regardless of event row order", {
  sim <- sim_glms(seed = 21, n_runs = 3, n_time = 100, n_vox = 6)
  args <- list(tr = sim$tr, stim_dur = sim$stim_dur, max_pcs = 2,
               noise_pool_r2 = 2, pc_voxel_r2 = 2, verbose = FALSE)
  ref <- do.call(glmsingle, c(list(Y = sim$Y, design = sim$design), args))
  ids <- c(20L, 10L, 30L)
  events <- do.call(rbind, lapply(seq_along(sim$design), function(r) {
    w <- which(sim$design[[r]] == 1, arr.ind = TRUE)
    data.frame(run = ids[r], onset = (w[, 1] - 1) * sim$tr, condition = w[, 2])
  }))
  events <- events[rev(seq_len(nrow(events))), ]
  combined <- list(Y = do.call(rbind, sim$Y), design = events,
                   runs = rep(ids, vapply(sim$Y, nrow, integer(1))))
  fit <- do.call(glmsingle, c(combined, args))
  expect_identical(fit$design$onsets, ref$design$onsets)
  for (type in c("typeb", "typec", "typed")) {
    expect_equal(fit[[type]]$betasmd, ref[[type]]$betasmd, tolerance = 1e-10)
  }
  expect_identical(fit$typed$FRACvalue, ref$typed$FRACvalue)
  combined$design$run[1] <- 999L
  expect_error(do.call(glmsingle, c(combined, args)), "run IDs must match")
})

test_that("unused and degenerate automatic thresholds are well defined", {
  expect_equal(.glms_tail_threshold(7), 7)
  expect_equal(.glms_tail_threshold(rep(7, 5)), 7)
  expect_equal(.glms_tail_threshold(c(NA, NaN, Inf)), 0)
  sim <- sim_glms(seed = 21, n_runs = 3, n_time = 100, n_vox = 6)
  for (columns in list(1L, rep(1L, 4))) {
    Y <- lapply(sim$Y, function(y) y[, columns, drop = FALSE])
    fit <- glmsingle(Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
                     denoise = FALSE, ridge = FALSE, verbose = FALSE)
    expect_null(fit$settings$tail_threshold)
    expect_true(all(is.finite(fit$typeb$betasmd)))
    # Exercise the default threshold through a denoising call too.
    expect_warning(dn <- glmsingle(Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
                                   ridge = FALSE, verbose = FALSE), "noise-pool rank")
    expect_equal(dn$settings$tail_threshold, dn$typea$onoffR2[1])
    expect_identical(dn$typec$pcnum, 0L)
    expect_equal(dn$typec$betasmd, dn$typeb$betasmd, tolerance = 1e-10)
  }
})

test_that("empty noise pools cannot inject arbitrary regressors", {
  sim <- sim_glms(seed = 21, n_runs = 3, n_time = 100, n_vox = 6)
  for (pc_args in list(list(n_pcs = 2), list(pc_stop = 1.05))) {
    expect_warning(fit <- do.call(glmsingle, c(list(
      sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
      ridge = FALSE, max_pcs = 2,
      noise_pool_mask = rep(FALSE, 6), noise_pool_r2 = 2, pc_voxel_r2 = 2,
      verbose = FALSE
    ), pc_args)), "Limiting the number of noise PCs to 0")
    expect_false(any(fit$typec$noisepool))
    expect_identical(fit$typec$pcnum, 0L)
    expect_true(all(vapply(fit$typec$pcregressors, ncol, integer(1)) == 0L))
    expect_equal(fit$typec$betasmd, fit$typeb$betasmd, tolerance = 1e-10)
  }
})

test_that("noise PCs are restricted to the residual noise-pool span", {
  set.seed(101)
  n <- 60L
  Q <- .glms_poly_basis(n, 1L)
  signal <- .glms_resid(matrix(rnorm(n), n), Q)
  Y <- cbind(100 + signal, 100 + 2 * signal, 100)
  pcs <- .glms_noise_pcs(Y, rep(TRUE, 3), Q, 10L, 2L)
  expect_equal(ncol(pcs), 1L)
  expect_equal(abs(cor(pcs[, 1], signal[, 1])), 1, tolerance = 1e-12)
  expect_equal(ncol(.glms_noise_pcs(matrix(100, n, 2), c(TRUE, TRUE), Q, 10L, 2L)), 0L)

  # The common PC count is bounded by the weakest run, including forced PCs.
  sim <- sim_glms(seed = 21, n_runs = 3, n_time = 100, n_vox = 6)
  expect_warning(fit <- glmsingle(
    sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
    noise_pool_mask = c(TRUE, rep(FALSE, 5)), noise_pool_r2 = Inf,
    pc_voxel_r2 = 2, n_pcs = 3,
    ridge = FALSE, verbose = FALSE
  ), "Limiting the number of noise PCs to 1")
  expect_identical(fit$typec$pcnum, 1L)
  expect_true(all(vapply(fit$typec$pcregressors, ncol, integer(1)) == 1L))
})

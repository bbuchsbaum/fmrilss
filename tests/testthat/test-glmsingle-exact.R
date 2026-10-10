# Float64 exactness: the fast glmsingle() pipeline against dense references
# that rebuild GLMsingle's stacked single-trial computations literally.

glms_sim_fit <- local({
  cache <- list()
  function(key, ...) {
    if (is.null(cache[[key]])) {
      sim <- sim_glms(...)
      fit <- glmsingle(sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
                       nuisance = sim$extras, max_pcs = 3, verbose = FALSE,
                       noise_pool_r2 = 2, pc_voxel_r2 = 2, chunk_size = 7)
      cache[[key]] <<- list(sim = sim, fit = fit)
    }
    cache[[key]]
  }
})

test_that("type A and B match dense OLS fits", {
  x <- glms_sim_fit("plain", seed = 21)
  fit <- x$fit; sim <- x$sim
  deg <- fit$settings$poly_degree
  for (h in unique(fit$typeb$HRFindex)) {
    v <- which(fit$typeb$HRFindex == h)
    ref <- ref_fit_assume(lapply(sim$Y, function(y) y[, v, drop = FALSE]),
                          fit$design$onsets, fit$hrf_library[, h], deg)
    pb <- 100 / abs(fit$meanvol[v])
    expect_equal(unname(fit$typeb$betasmd[, v, drop = FALSE]), sweep(ref$beta, 2L, pb, "*"), tolerance = 1e-9)
    expect_equal(fit$typeb$R2[v], ref$R2, tolerance = 1e-9)
    expect_equal(fit$typeb$R2run[v, , drop = FALSE], ref$R2run, tolerance = 1e-9)
  }
  # the selected HRF maximises R^2 over the library
  expect_equal(fit$typeb$R2, apply(fit$typeb$FitHRFR2, 1, max))
})

test_that("PC-count cross-validation matches dense fits and calcbadness", {
  x <- glms_sim_fit("extras", seed = 22, n_extras = 2)
  fit <- x$fit; sim <- x$sim
  deg <- fit$settings$poly_degree
  g <- fit$design
  stimix <- lapply(g$validcolumns, function(cc) g$stimorder[cc])
  ix <- which(fit$typec$pcvoxels)
  for (h in unique(fit$typeb$HRFindex[ix])) {
    v <- ix[fit$typeb$HRFindex[ix] == h]
    res <- lapply(0:3, function(k) {
      ex <- lapply(seq_along(sim$Y), function(r) {
        cbind(sim$extras[[r]], fit$typec$pcregressors[[r]][, seq_len(k), drop = FALSE])
      })
      t(ref_fit_assume(lapply(sim$Y, function(y) y[, v, drop = FALSE]), g$onsets,
                       fit$hrf_library[, h], deg, ex)$beta)
    })
    ref <- ref_calcbadness(as.list(seq_along(sim$Y)), g$validcolumns, stimix, res,
                           rep(1, length(sim$Y)), python = FALSE)
    expect_equal(fit$typec$glmbadness[v, , drop = FALSE], ref, tolerance = 1e-8)
  }
})

test_that("type C and D match dense fracridge, calcbadness and autoscale", {
  x <- glms_sim_fit("extras", seed = 22, n_extras = 2)
  fit <- x$fit; sim <- x$sim
  deg <- fit$settings$poly_degree
  g <- fit$design
  stimix <- lapply(g$validcolumns, function(cc) g$stimorder[cc])
  pcnum <- fit$typec$pcnum
  ex <- lapply(seq_along(sim$Y), function(r) {
    cbind(sim$extras[[r]], fit$typec$pcregressors[[r]][, seq_len(pcnum), drop = FALSE])
  })
  fracs <- fit$settings$ridge_fracs
  for (h in unique(fit$typeb$HRFindex)) {
    v <- which(fit$typeb$HRFindex == h)
    Yv <- lapply(sim$Y, function(y) y[, v, drop = FALSE])
    pb <- 100 / abs(fit$meanvol[v])
    fits <- lapply(fracs, function(f) ref_fit_assume(Yv, g$onsets, fit$hrf_library[, h], deg, ex, frac = f))
    expect_equal(unname(fit$typec$betasmd[, v, drop = FALSE]), sweep(fits[[1]]$beta, 2L, pb, "*"), tolerance = 1e-8)
    expect_equal(fit$typec$R2[v], fits[[1]]$R2, tolerance = 1e-8)
    bad <- ref_calcbadness(as.list(seq_along(sim$Y)), g$validcolumns, stimix,
                           lapply(fits, function(f) t(f$beta)), rep(1, length(sim$Y)), python = FALSE)
    expect_equal(fit$typed$rrbadness[v, , drop = FALSE], bad, tolerance = 1e-8)
    idx <- apply(bad, 1, which.min)
    expect_equal(fit$typed$FRACvalue[v], fracs[idx])
    for (j in seq_along(v)) {
      b <- fits[[idx[j]]]$beta[, j]
      X <- cbind(b, 1)
      hh <- drop(ref_olsmatrix(X) %*% fits[[1]]$beta[, j])
      if (hh[1] < 0) hh <- c(1, 0)
      expect_equal(unname(fit$typed$betasmd[, v[j]]), drop(X %*% hh) * pb[j], tolerance = 1e-8)
      expect_equal(fit$typed$R2[v[j]], fits[[idx[j]]]$R2[j], tolerance = 1e-8)
    }
  }
})

test_that("results do not depend on chunking or input layout", {
  x <- glms_sim_fit("plain", seed = 21)
  sim <- x$sim
  runs <- rep(seq_along(sim$Y), vapply(sim$Y, nrow, integer(1)))
  ev <- do.call(rbind, lapply(seq_along(sim$design), function(r) {
    w <- which(sim$design[[r]] == 1, arr.ind = TRUE)
    data.frame(run = r, onset = (w[, 1] - 1) * sim$tr, condition = w[, 2])
  }))
  fit2 <- glmsingle(do.call(rbind, sim$Y), ev, tr = sim$tr, stim_dur = sim$stim_dur,
                    runs = runs, max_pcs = 3, verbose = FALSE, noise_pool_r2 = 2,
                    pc_voxel_r2 = 2, chunk_size = 1000)
  for (tp in c("typeb", "typec", "typed")) {
    expect_equal(unname(fit2[[tp]]$betasmd), unname(x$fit[[tp]]$betasmd), tolerance = 1e-10)
  }
  expect_equal(fit2$typed$FRACvalue, x$fit$typed$FRACvalue)
})

test_that("restricting PC cross-validation to the selection voxels is exact", {
  x <- glms_sim_fit("plain", seed = 21)
  sim <- x$sim
  full <- glmsingle(sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur, max_pcs = 3,
                    verbose = FALSE, noise_pool_r2 = 2, pc_voxel_r2 = 2, pc_cv_all_voxels = TRUE)
  expect_equal(full$typec$xvaltrend, x$fit$typec$xvaltrend)
  expect_identical(full$typec$pcnum, x$fit$typec$pcnum)
  expect_equal(full$typed$betasmd, x$fit$typed$betasmd)
  expect_true(all(is.finite(full$typec$glmbadness)))
  sel <- which(x$fit$typec$pcvoxels)
  expect_equal(full$typec$glmbadness[sel, ], x$fit$typec$glmbadness[sel, ])
  expect_true(all(is.na(x$fit$typec$glmbadness[-sel, ])))
})

test_that("extras policy only matters when extra regressors are supplied", {
  x <- glms_sim_fit("plain", seed = 21)
  sim <- x$sim
  alt <- glmsingle(sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur, max_pcs = 3,
                   verbose = FALSE, noise_pool_r2 = 2, pc_voxel_r2 = 2,
                   nuisance_in_denoise = "with_pcs", chunk_size = 7)
  expect_equal(alt$typed$betasmd, x$fit$typed$betasmd)
})

test_that("singular trial designs error by default and fall back with pinv", {
  sim <- sim_glms(seed = 4, n_runs = 2, n_time = 60, n_vox = 5)
  impulse <- matrix(c(1, 0, 0), ncol = 1)
  ex <- lapply(sim$design, function(D) {
    on <- which(rowSums(D) > 0)
    E <- matrix(0, nrow(D), 2); E[on[1], 1] <- 1; E[on[2], 2] <- 1
    E
  })
  args <- list(Y = sim$Y, design = sim$design, tr = sim$tr, stim_dur = sim$stim_dur,
               hrf_library = impulse, nuisance = ex, verbose = FALSE,
               denoise = FALSE, noise_pool_r2 = 2, pc_voxel_r2 = 2)
  expect_error(do.call(glmsingle, args), "Singular")
  w <- testthat::capture_warnings(fit <- do.call(glmsingle, c(args, list(singular = "pinv"))))
  expect_length(w, 1L)
  expect_match(w, "minimum-norm")
  expect_true(all(is.finite(fit$typeb$betasmd)))
})

test_that("input validation catches malformed designs", {
  sim <- sim_glms(seed = 4, n_runs = 2, n_time = 60, n_vox = 5)
  bad <- sim$design; bad[[1]][10, 1:2] <- 1
  expect_error(glmsingle(sim$Y, bad, tr = 1, stim_dur = 3, verbose = FALSE), "same trial onset")
  expect_error(glmsingle(sim$Y, sim$design[1], tr = 1, stim_dur = 3, verbose = FALSE), "one time x condition")
  ev <- data.frame(run = 1, onset = 2.5, condition = 1)
  expect_error(glmsingle(sim$Y[1], ev, tr = 1, stim_dur = 3, verbose = FALSE), "TR grid")
  expect_error(glmsingle(sim$Y, sim$design, tr = 1, stim_dur = 3, ridge_fracs = 1.5, verbose = FALSE), "ridge_fracs")
})

test_that("fmridesign front end matches the matrix interface", {
  skip_if_not_installed("fmridesign")
  sim <- sim_glms(seed = 31, n_runs = 3, n_time = 80, n_vox = 10)
  ev <- do.call(rbind, lapply(seq_along(sim$design), function(r) {
    w <- which(sim$design[[r]] == 1, arr.ind = TRUE)
    w <- w[order(w[, 1]), , drop = FALSE]
    data.frame(run = r, onset = (w[, 1] - 1) * sim$tr, stim = factor(w[, 2], levels = 1:8),
               duration = sim$stim_dur)
  }))
  sf <- fmrihrf::sampling_frame(blocklens = vapply(sim$Y, nrow, integer(1)), TR = sim$tr)
  em <- fmridesign::event_model(onset ~ fmridesign::hrf(stim), data = ev, block = ~run,
                                sampling_frame = sf, durations = ev$duration)
  a <- glmsingle_design(do.call(rbind, sim$Y), em, max_pcs = 2, verbose = FALSE,
                        noise_pool_r2 = 2, pc_voxel_r2 = 2)
  b <- glmsingle(sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur, max_pcs = 2,
                 verbose = FALSE, noise_pool_r2 = 2, pc_voxel_r2 = 2)
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
  expect_output(print(a), "glmsingle_fit")
  expect_s3_class(summary(a), "summary.glmsingle_fit")
  expect_length(coef(a, "a"), 10L)
})

test_that("results do not depend on the thread count", {
  x <- glms_sim_fit("plain", seed = 21)
  sim <- x$sim
  one <- glmsingle(sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur, max_pcs = 3,
                   verbose = FALSE, noise_pool_r2 = 2, pc_voxel_r2 = 2, n_threads = 1)
  many <- glmsingle(sim$Y, sim$design, tr = sim$tr, stim_dur = sim$stim_dur, max_pcs = 3,
                    verbose = FALSE, noise_pool_r2 = 2, pc_voxel_r2 = 2, n_threads = 2)
  expect_identical(one$typed$betasmd, many$typed$betasmd)
  expect_identical(one$typed$FRACvalue, many$typed$FRACvalue)
})

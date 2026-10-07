skip_if_not_installed("fmrihrf")

make_rank1_problem <- function(n = 160, n_trials = 22, V = 4, sd = 0.3,
                               seed = 3, basis = fmrihrf::HRF_SPMG3,
                               h = c(1, 0.5, -0.3), mean_amp = 1) {
  set.seed(seed)
  sframe <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  events <- data.frame(onset = sort(runif(n_trials, 5, n - 25)), duration = 0,
                       condition = rep(c("A", "B"), length.out = n_trials))
  X <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)$X
  K <- length(h)
  B <- matrix(rnorm(n_trials * V, mean_amp, 0.4), n_trials, V)
  N <- matrix(rnorm(n * 2), n, 2)
  Y <- X %*% kronecker(B, h) + N %*% matrix(rnorm(2 * V), 2, V) +
    matrix(rnorm(n * V, sd = sd), n, V)
  list(Y = Y, X = X, B = B, N = N, events = events, sframe = sframe,
       basis = basis, K = K, h = h, n = n, T = n_trials, V = V)
}

# Brute-force alternating least squares with explicit n-dimensional fits.
brute_rank1 <- function(p, model, iters, groups = NULL) {
  K <- p$K
  Q <- qr.Q(qr(cbind(1, p$N)))
  Xr <- p$X - Q %*% crossprod(Q, p$X)
  Yr <- p$Y - Q %*% crossprod(Q, p$Y)
  blk <- function(i) Xr[, ((i - 1) * K + 1):(i * K)]
  codes <- if (is.null(groups)) rep(1L, p$T) else match(groups, unique(groups))
  A <- lapply(seq_len(max(codes)), function(g) Reduce(`+`, lapply(which(codes == g), blk)))
  H <- matrix(NA, K, p$V)
  out <- matrix(NA, p$T, p$V)
  for (v in seq_len(p$V)) {
    y <- Yr[, v]
    h <- qr.coef(qr(Reduce(`+`, A)), y)
    h <- h / sqrt(sum(h^2))
    trial_design <- function(i, h) {
      ci <- blk(i) %*% h
      others <- vapply(seq_along(A), function(g) {
        drop(A[[g]] %*% h - (codes[i] == g) * ci)
      }, numeric(p$n))
      list(design = cbind(ci, others), keep = c(TRUE, colSums(others^2) > 1e-12))
    }
    bstep <- function(h) {
      if (model == "joint") {
        b <- qr.coef(qr(vapply(seq_len(p$T), function(i) drop(blk(i) %*% h), numeric(p$n))), y)
        return(cbind(b, matrix(0, p$T, length(A))))
      }
      t(vapply(seq_len(p$T), function(i) {
        td <- trial_design(i, h)
        coef <- numeric(length(td$keep))
        coef[td$keep] <- qr.coef(qr(td$design[, td$keep, drop = FALSE]), y)
        coef
      }, numeric(1 + length(A))))
    }
    for (it in seq_len(iters)) {
      br <- bstep(h)
      h <- if (model == "joint") {
        qr.coef(qr(Reduce(`+`, lapply(seq_len(p$T), function(i) br[i, 1] * blk(i)))), y)
      } else {
        D <- do.call(rbind, lapply(seq_len(p$T), function(i) {
          r <- br[i, -1]
          (br[i, 1] - r[codes[i]]) * blk(i) + Reduce(`+`, Map(`*`, r, A))
        }))
        qr.coef(qr(D), rep(y, p$T))
      }
      h <- h / sqrt(sum(h^2))
    }
    br <- bstep(h)
    H[, v] <- h
    out[, v] <- br[, 1]
  }
  list(H = H, beta = out)
}

products <- function(H, B) {
  vapply(seq_len(ncol(B)), function(v) kronecker(B[, v], H[, v]), numeric(nrow(B) * nrow(H)))
}

test_that("ALS kernels match brute-force alternating least squares", {
  p <- make_rank1_problem()
  iters <- 8
  for (model in c("separate", "joint")) {
    ref <- brute_rank1(p, model, iters)
    fit <- lss_rank1(p$Y, p$events, p$basis, p$sframe, nuisance_regs = p$N,
                     model = model, max_iter = iters, tol = 0)
    # h and beta are identified only up to scale; compare h (x) beta
    expect_equal(products(fit$hrf$coefficients, fit$beta), products(ref$H, ref$beta),
                 tolerance = 1e-6, info = model)
  }
  g <- p$events$condition
  ref <- brute_rank1(p, "separate", iters, groups = g)
  fit <- lss_rank1(p$Y, p$events, p$basis, p$sframe, nuisance_regs = p$N,
                   trial_groups = g, max_iter = iters, tol = 0)
  expect_equal(products(fit$hrf$coefficients, fit$beta), products(ref$H, ref$beta),
               tolerance = 1e-6)
})

test_that("the objective never increases across iterations", {
  p <- make_rank1_problem(sd = 1)
  for (args in list(list(model = "separate"), list(model = "joint"),
                    list(model = "separate", trial_groups = p$events$condition))) {
    objs <- vapply(1:6, function(it) {
      sum(do.call(lss_rank1, c(list(p$Y, p$events, p$basis, p$sframe,
                                    init = "reference", max_iter = it, tol = 0),
                               args))$objective)
    }, numeric(1))
    expect_true(all(diff(objs) <= 1e-8 * abs(objs[-1])), info = paste(unlist(args[1]), collapse = ""))
  }
})

test_that("the joint model is exact without noise", {
  p <- make_rank1_problem(sd = 1e-6, V = 3)
  wf <- fmrilss:::.rank1_waveforms(p$basis, fmrihrf::HRF_SPMG1, 24)
  truth <- wf$H %*% p$h
  fit <- lss_rank1(p$Y, p$events, p$basis, p$sframe, nuisance_regs = p$N,
                   model = "joint")
  est <- wf$H %*% fit$hrf$coefficients
  expect_true(all(fit$converged))
  expect_gt(min(apply(est, 2, cor, truth)), 0.99999)
  expect_gt(min(vapply(1:3, function(v) cor(fit$beta[, v], p$B[, v]), 1)), 0.99999)
  # unit positive peak, positively oriented
  expect_equal(unname(apply(est, 2, max)), rep(1, 3), tolerance = 1e-8)
})

test_that("the separate model recovers the HRF in a slow design", {
  # LSS pools all other trials into one regressor, which is exact only when
  # trial responses do not overlap; in rapid designs it is biased even
  # without noise.
  set.seed(5)
  n <- 400
  sframe <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  events <- data.frame(onset = seq(10, 370, by = 15), duration = 0, condition = "A")
  basis <- fmrihrf::HRF_SPMG3
  h <- c(1, 0.5, -0.3)
  X <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)$X
  B <- matrix(rnorm(nrow(events) * 3, 1, 0.4), nrow(events), 3)
  Y <- X %*% kronecker(B, h) + matrix(rnorm(n * 3, sd = 1e-4), n, 3)
  fit <- lss_rank1(Y, events, basis, sframe)
  wf <- fmrilss:::.rank1_waveforms(basis, fmrihrf::HRF_SPMG1, 24)
  expect_gt(min(apply(wf$H %*% fit$hrf$coefficients, 2, cor, wf$H %*% h)), 0.999)
  # the HRF undershoot still overlaps the next trial slightly
  expect_gt(min(vapply(1:3, function(v) cor(fit$beta[, v], B[, v]), 1)), 0.99)
})

test_that("the separate-model amplitudes are LSS estimates under the learned HRF", {
  p <- make_rank1_problem()
  fit <- lss_rank1(p$Y, p$events, p$basis, p$sframe, nuisance_regs = p$N)
  expect_s3_class(fit$hrf, "VoxelHRF")
  refit <- lss_with_hrf(p$Y, p$events, fit$hrf, nuisance_regs = p$N, verbose = FALSE)
  expect_equal(unname(unclass(refit)[seq_len(p$T), ]), unname(fit$beta),
               tolerance = 1e-8)
})

test_that("L-BFGS-B and ALS reach the same optimum", {
  p <- make_rank1_problem(V = 3)
  als <- lss_rank1(p$Y, p$events, p$basis, p$sframe, tol = 1e-12, max_iter = 1000)
  lbfgs <- lss_rank1(p$Y, p$events, p$basis, p$sframe, solver = "lbfgs",
                     tol = 10 * .Machine$double.eps)
  expect_equal(lbfgs$objective, als$objective, tolerance = 1e-8)
  expect_equal(products(lbfgs$hrf$coefficients, lbfgs$beta),
               products(als$hrf$coefficients, als$beta), tolerance = 1e-4)
})

test_that("trial_groups returns per-group other-trial amplitudes", {
  p <- make_rank1_problem()
  g <- p$events$condition
  fit <- lss_rank1(p$Y, p$events, p$basis, p$sframe, trial_groups = g)
  expect_equal(dim(fit$other), c(p$T, 2L, p$V))
  expect_identical(dimnames(fit$other)[[2]], c("A", "B"))
  ungrouped <- lss_rank1(p$Y, p$events, p$basis, p$sframe)
  expect_equal(dim(ungrouped$other), c(p$T, p$V))
  expect_equal(lss_rank1(p$Y, p$events, p$basis, p$sframe, trial_groups = rep("x", p$T)),
               ungrouped)
})

test_that("trial_groups recovers the HRF when conditions have opposite signs", {
  set.seed(8)
  p <- make_rank1_problem(n = 300, n_trials = 50, V = 40, sd = 0.3, mean_amp = 0)
  amp <- ifelse(p$events$condition == "A", 1, -1)
  p$B <- matrix(rnorm(p$T * p$V, amp, 0.3), p$T, p$V)
  p$Y <- p$X %*% kronecker(p$B, p$h) + matrix(rnorm(p$n * p$V, sd = 0.3), p$n, p$V)
  wf <- fmrilss:::.rank1_waveforms(p$basis, fmrihrf::HRF_SPMG1, 24)
  truth <- wf$H %*% p$h
  hrf_r <- function(fit) mean(apply(wf$H %*% fit$hrf$coefficients, 2, cor, truth))
  pooled <- lss_rank1(p$Y, p$events, p$basis, p$sframe)
  grouped <- lss_rank1(p$Y, p$events, p$basis, p$sframe,
                       trial_groups = p$events$condition)
  expect_gt(hrf_r(grouped), 0.95)
  expect_gt(hrf_r(grouped), hrf_r(pooled) + 0.1)
})

test_that("prewhitening is applied and recorded", {
  skip_if_not_installed("fmriAR")
  p <- make_rank1_problem(V = 30)
  fit <- lss_rank1(p$Y, p$events, p$basis, p$sframe,
                   prewhiten = list(method = "ar", p = 1))
  expect_s3_class(attr(fit, "whiten_plan"), "fmriAR_plan")
  expect_true(all(is.finite(fit$beta)))
  expect_error(
    lss_rank1(p$Y, p$events, p$basis, p$sframe,
              prewhiten = list(method = "ar", p = 1, pooling = "voxel")),
    "global' or 'run"
  )
})

test_that("lss_rank1 validates its arguments", {
  p <- make_rank1_problem()
  expect_error(lss_rank1(p$Y, p$events, p$basis, p$sframe, model = "joint",
                         solver = "lbfgs"), "separate")
  expect_error(lss_rank1(p$Y, p$events, p$basis, p$sframe, model = "joint",
                         trial_groups = p$events$condition), "separate")
  expect_error(lss_rank1(p$Y, p$events, p$basis, p$sframe, init = matrix(1, 2, 2)),
               "K x V")
  expect_error(lss_rank1(p$Y, p$events, "spmg3", p$sframe), "HRF")
  init <- matrix(c(1, 0, 0), 3, p$V)
  expect_equal(lss_rank1(p$Y, p$events, p$basis, p$sframe, init = init)$beta,
               lss_rank1(p$Y, p$events, p$basis, p$sframe, init = "reference")$beta,
               tolerance = 1e-6)
})

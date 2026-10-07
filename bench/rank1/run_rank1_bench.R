# Rank-1 GLM (lss_rank1) benchmarks.
#
# Usage: Rscript bench/rank1/run_rank1_bench.R [out_dir]
#
# 1. Solver comparison: exact alternating least squares (solver = "als") vs
#    joint L-BFGS-B (solver = "lbfgs", the approach of Pedregosa et al.,
#    2015) on the same objective and starting point, single-threaded.
# 2. Estimator comparison across amplitude regimes: canonical-HRF LSS,
#    estimate_voxel_hrf() + lss_with_hrf(), lss_rank1() variants, and
#    oracle LSS / LSS-N with each voxel's true HRF.

suppressPackageStartupMessages(library(fmrilss))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else "bench/rank1/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- simulation -------------------------------------------------------------
# 400 scans at TR 1 s; 10 groups of voxels whose double-gamma HRFs peak
# between 4 and 7.5 s; AR(1) noise.
simulate <- function(V, iti, amplitudes, noise_sd, seed = 11) {
  set.seed(seed)
  n <- 400
  sframe <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  on <- cumsum(c(12, runif(400, iti[1], iti[2])))
  on <- on[on < n - 30]
  Tn <- length(on)
  condition <- rep(c("A", "B"), length.out = Tn)
  events <- data.frame(onset = on, duration = 0, condition = condition)
  dt <- 0.1
  tt <- seq(0, 30, by = dt)
  ttp <- seq(4, 7.5, length.out = 10)
  shape_of <- function(p) {
    a <- p / 0.9 + 1
    h <- dgamma(tt, a, scale = 0.9) - dgamma(tt, a + 10, scale = 0.9) / 6
    h / max(h)
  }
  grp <- sample(1:10, V, replace = TRUE)
  th <- seq(0, n - 1)
  thi <- seq(0, n + 30, by = dt)
  Xg <- lapply(1:10, function(g) {
    h <- shape_of(ttp[g])
    pad <- length(h) - 1
    sapply(on, function(o) {
      u <- numeric(length(thi))
      u[which.min(abs(thi - o))] <- 1
      cv <- stats::filter(c(numeric(pad), u), h, sides = 1)[-seq_len(pad)]
      stats::approx(thi, cv, xout = th)$y
    })
  })
  B <- switch(amplitudes,
    positive = matrix(rnorm(Tn * V, 1, 0.5), Tn, V),
    zero_mean = matrix(rnorm(Tn * V, 0, 1), Tn, V),
    opposite = matrix(rnorm(Tn * V, ifelse(condition == "A", 1, -1), 0.5), Tn, V)
  )
  Y <- matrix(0, n, V)
  for (g in 1:10) {
    idx <- which(grp == g)
    Y[, idx] <- Xg[[g]] %*% B[, idx]
  }
  noise <- apply(matrix(rnorm(n * V), n, V), 2, function(e) {
    as.numeric(stats::filter(e, 0.3, "recursive"))
  })
  truth <- sapply(1:10, function(g) {
    stats::approx(tt, shape_of(ttp[g]), xout = seq(0, 24, by = 0.1))$y
  })
  list(Y = Y + noise_sd * noise + 50, B = B, events = events, sframe = sframe,
       Xg = Xg, grp = grp, truth = truth, V = V, Tn = Tn)
}

beta_r <- function(est, d) mean(vapply(seq_len(d$V), function(v) cor(est[, v], d$B[, v]), 1))
hrf_r <- function(fit, d) {
  wf <- fmrilss:::.rank1_waveforms(fit$hrf$basis, fmrihrf::HRF_SPMG1, 24, precision = 0.1)
  W <- wf$H %*% fit$hrf$coefficients
  mean(vapply(seq_len(d$V), function(v) cor(W[, v], d$truth[, d$grp[v]]), 1))
}
bases <- list(SPMG3 = fmrihrf::HRF_SPMG3,
              FIR12 = fmrihrf::hrf_fir_generator(nbasis = 12, span = 24))
scenarios <- list(
  positive = list(iti = c(3, 8), amplitudes = "positive", noise_sd = 0.8),
  positive_rapid = list(iti = c(1.5, 4), amplitudes = "positive", noise_sd = 1.5),
  zero_mean = list(iti = c(3, 8), amplitudes = "zero_mean", noise_sd = 0.8),
  opposite = list(iti = c(3, 8), amplitudes = "opposite", noise_sd = 0.8)
)

# --- 1. solver comparison ---------------------------------------------------
Sys.setenv(OMP_NUM_THREADS = 1)
solver_rows <- list()
for (sc in c("positive", "positive_rapid")) {
  d <- do.call(simulate, c(list(V = 1000), scenarios[[sc]]))
  for (bn in names(bases)) for (init in c("aggregate", "reference")) {
    fits <- list()
    settings <- list(
      list(label = "ALS, tol 1e-7", solver = "als", tol = 1e-7),
      list(label = "ALS, tol 1e-10", solver = "als", tol = 1e-10),
      list(label = "L-BFGS-B, factr 1e7", solver = "lbfgs", tol = 1e7 * .Machine$double.eps),
      list(label = "L-BFGS-B, factr 10", solver = "lbfgs", tol = 10 * .Machine$double.eps)
    )
    for (st in settings) {
      sec <- system.time(f <- lss_rank1(d$Y, d$events, bases[[bn]], d$sframe,
                                         init = init, solver = st$solver,
                                         tol = st$tol, max_iter = 2000))[["elapsed"]]
      fits[[st$label]] <- f
      solver_rows[[length(solver_rows) + 1L]] <- data.frame(
        scenario = sc, basis = bn, init = init, solver = st$label, seconds = sec,
        median_iterations = stats::median(f$iterations),
        converged = mean(f$converged), beta_r = beta_r(f$beta, d))
    }
    best <- do.call(pmin, lapply(fits, `[[`, "objective"))
    k <- length(solver_rows)
    for (j in seq_along(fits)) {
      gap <- (fits[[j]]$objective - best) / best
      solver_rows[[k - length(fits) + j]]$median_rel_gap <- stats::median(gap)
      solver_rows[[k - length(fits) + j]]$max_rel_gap <- max(gap)
    }
  }
}
Sys.unsetenv("OMP_NUM_THREADS")
solver_tab <- do.call(rbind, solver_rows)
utils::write.csv(solver_tab, file.path(out_dir, "solver_comparison.csv"), row.names = FALSE)
print(solver_tab, digits = 3, row.names = FALSE)

# --- 2. estimator comparison -------------------------------------------------
est_rows <- list()
add <- function(sc, method, sec, est, d, hr = NA_real_) {
  est_rows[[length(est_rows) + 1L]] <<- data.frame(
    scenario = sc, method = method, seconds = sec, beta_r = beta_r(est, d), hrf_r = hr)
}
for (sc in names(scenarios)) {
  d <- do.call(simulate, c(list(V = 2000), scenarios[[sc]]))
  g <- d$events$condition
  orc <- orc_n <- matrix(0, d$Tn, d$V)
  for (k in 1:10) {
    idx <- which(d$grp == k)
    orc[, idx] <- lss(d$Y[, idx, drop = FALSE], d$Xg[[k]])
    orc_n[, idx] <- lss(d$Y[, idx, drop = FALSE], d$Xg[[k]], trial_groups = g)
  }
  add(sc, "oracle LSS (true HRF)", NA, orc, d)
  add(sc, "oracle LSS-N (true HRF)", NA, orc_n, d)
  Xc <- fmrilss:::.voxhrf_trial_basis(d$events, fmrihrf::HRF_SPMG1, d$sframe)$X
  add(sc, "LSS, canonical HRF", system.time(e <- lss(d$Y, Xc))[["elapsed"]], e, d)
  add(sc, "LSS-N, canonical HRF",
      system.time(e <- lss(d$Y, Xc, trial_groups = g))[["elapsed"]], e, d)
  for (bn in names(bases)) {
    b <- bases[[bn]]
    sec <- system.time(e <- tryCatch({
      vh <- estimate_voxel_hrf(d$Y, d$events, b, sframe = d$sframe)
      unclass(lss_with_hrf(d$Y, d$events, vh, verbose = FALSE))[seq_len(d$Tn), , drop = FALSE]
    }, error = function(err) NULL))[["elapsed"]]
    if (is.null(e)) {
      est_rows[[length(est_rows) + 1L]] <- data.frame(
        scenario = sc, method = paste("estimate_voxel_hrf + lss_with_hrf,", bn, "(error)"),
        seconds = NA, beta_r = NA, hrf_r = NA)
    } else {
      add(sc, paste("estimate_voxel_hrf + lss_with_hrf,", bn), sec, e, d)
    }
    sec <- system.time(f <- lss_rank1(d$Y, d$events, b, d$sframe))[["elapsed"]]
    add(sc, paste("lss_rank1 separate,", bn), sec, f$beta, d, hrf_r(f, d))
    sec <- system.time(f <- lss_rank1(d$Y, d$events, b, d$sframe, trial_groups = g))[["elapsed"]]
    add(sc, paste("lss_rank1 separate + trial_groups,", bn), sec, f$beta, d, hrf_r(f, d))
    sec <- system.time(f <- lss_rank1(d$Y, d$events, b, d$sframe, model = "joint"))[["elapsed"]]
    add(sc, paste("lss_rank1 joint,", bn), sec, f$beta, d, hrf_r(f, d))
  }
}
est_tab <- do.call(rbind, est_rows)
utils::write.csv(est_tab, file.path(out_dir, "estimator_comparison.csv"), row.names = FALSE)
print(est_tab, digits = 3, row.names = FALSE)

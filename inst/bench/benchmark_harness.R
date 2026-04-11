# Benchmark harness under active development.
# Installed copy for package tests and local benchmark scripts.

# Utility: light null coalesce.
.bmk_or <- function(x, y) {
  if (is.null(x)) y else x
}

# Utility: clamp numeric values to [lo, hi].
.bmk_clip <- function(x, lo, hi) {
  pmax(lo, pmin(hi, x))
}

# Utility: normalize prewhiten aliases used across code paths.
.bmk_normalize_prewhiten <- function(prewhiten) {
  if (is.null(prewhiten) || !is.list(prewhiten)) return(prewhiten)
  out <- prewhiten
  if (!is.null(out$method)) {
    method <- tolower(as.character(out$method)[1L])
    if (identical(method, "ar1")) method <- "ar"
    out$method <- method
  }
  out
}

# Utility: repeatable split of total trial count across runs.
.bmk_split_trials <- function(n_trials, n_runs) {
  base <- rep(as.integer(n_trials %/% n_runs), as.integer(n_runs))
  rem <- as.integer(n_trials %% n_runs)
  if (rem > 0L) base[seq_len(rem)] <- base[seq_len(rem)] + 1L
  base
}

# Utility: sample ISI values under multiple distributions.
.bmk_draw_isi <- function(
  n,
  isi,
  isi_model = c("spread", "fixed", "uniform", "exponential", "gamma", "lognormal"),
  isi_params = list()
) {
  n <- as.integer(n)
  if (n <= 0L) return(numeric(0))

  isi <- as.numeric(isi)
  if (!is.finite(isi) || isi <= 0) {
    stop("isi must be a positive finite value.", call. = FALSE)
  }

  isi_model <- match.arg(isi_model)
  if (isi_model %in% c("spread", "fixed")) {
    return(rep(isi, n))
  }

  draws <- switch(
    isi_model,
    uniform = {
      width <- as.numeric(.bmk_or(isi_params$width, 0.5))
      width <- max(0, min(0.99, width))
      lo <- isi * (1 - width)
      hi <- isi * (1 + width)
      stats::runif(n, min = max(1e-4, lo), max = max(1e-4, hi))
    },
    exponential = {
      stats::rexp(n, rate = 1 / isi)
    },
    gamma = {
      shape <- as.numeric(.bmk_or(isi_params$shape, 2))
      if (!is.finite(shape) || shape <= 0) shape <- 2
      stats::rgamma(n, shape = shape, scale = isi / shape)
    },
    lognormal = {
      sdlog <- as.numeric(.bmk_or(isi_params$sdlog, 0.6))
      if (!is.finite(sdlog) || sdlog <= 0) sdlog <- 0.6
      meanlog <- log(isi) - 0.5 * sdlog^2
      stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog)
    },
    rep(isi, n)
  )

  pmax(1e-4, as.numeric(draws))
}

# Utility: build run-aware event table with optional onset jitter.
.bmk_make_events <- function(
  n_trials,
  n_runs,
  n_tp,
  tr,
  isi,
  isi_model = c("spread", "fixed", "uniform", "exponential", "gamma", "lognormal"),
  isi_params = list(),
  jitter_sd = 0,
  start_offset = 6,
  end_margin = NULL
) {
  n_runs <- as.integer(n_runs)
  n_tp <- as.integer(n_tp)
  tr <- as.numeric(tr)
  if (is.null(end_margin)) end_margin <- 2 * tr
  end_margin <- as.numeric(end_margin)
  isi_model <- match.arg(isi_model)

  trials_per_run <- .bmk_split_trials(as.integer(n_trials), n_runs)
  out <- vector("list", n_runs)
  trial_cursor <- 1L

  for (r in seq_len(n_runs)) {
    ntr <- trials_per_run[r]
    if (ntr < 1L) next

    run_dur <- n_tp[r] * tr
    start <- min(as.numeric(start_offset), max(0, run_dur - end_margin))
    end_cap <- max(start, run_dur - end_margin)

    if (ntr == 1L) {
      onset <- start
    } else if (isi_model == "spread") {
      onset <- seq(start, to = end_cap, length.out = ntr)
    } else {
      d_isi <- .bmk_draw_isi(ntr - 1L, isi = isi, isi_model = isi_model, isi_params = isi_params)
      onset <- start + c(0, cumsum(d_isi))
      max_onset <- max(onset)
      if (is.finite(max_onset) && max_onset > end_cap && max_onset > start) {
        onset <- start + (onset - start) * ((end_cap - start) / (max_onset - start))
      }
    }

    if (isTRUE(jitter_sd > 0)) {
      onset <- onset + stats::rnorm(ntr, sd = as.numeric(jitter_sd))
      onset <- pmax(0, pmin(end_cap, onset))
    }

    trial_idx <- seq.int(trial_cursor, length.out = ntr)
    trial_cursor <- trial_cursor + ntr
    out[[r]] <- data.frame(onset = as.numeric(onset), run = r, trial_id = trial_idx)
  }

  events <- do.call(rbind, out)
  rownames(events) <- NULL
  events
}

# Utility: apply true-onset shifts to induce model mismatch.
.bmk_apply_true_onset_shift <- function(events, n_tp, tr, shift_sd = 0, end_margin = NULL) {
  if (!isTRUE(shift_sd > 0)) return(events)

  tr <- as.numeric(tr)
  if (is.null(end_margin)) end_margin <- 2 * tr
  end_margin <- as.numeric(end_margin)

  out <- events
  runs <- sort(unique(as.integer(events$run)))
  for (r in runs) {
    idx <- which(events$run == r)
    if (length(idx) < 1L) next
    run_dur <- as.numeric(n_tp[r]) * tr
    end_cap <- max(0, run_dur - end_margin)
    shifted <- events$onset[idx] + stats::rnorm(length(idx), sd = as.numeric(shift_sd))
    out$onset[idx] <- pmax(0, pmin(end_cap, shifted))
  }
  out
}

# Utility: build trial-wise fmri design matrix via fmridesign.
.bmk_trial_design <- function(events, n_tp, tr, basis = "spmg1", precision = 0.3, return_model = FALSE) {
  basis_arg <- paste0("\"", gsub("\"", "\\\\\"", as.character(basis)), "\"")
  rhs <- sprintf(
    "fmridesign::trialwise(basis=%s, add_sum=FALSE)",
    basis_arg
  )
  form <- stats::as.formula(paste("onset ~", rhs))
  sframe <- fmrihrf::sampling_frame(blocklens = as.integer(n_tp), TR = as.numeric(tr))
  emod <- fmridesign::event_model(
    formula_or_list = form,
    data = as.data.frame(events),
    block = ~run,
    sampling_frame = sframe,
    durations = 0,
    precision = precision
  )

  X <- as.matrix(fmridesign::design_matrix(emod))
  if (ncol(X) != nrow(events)) {
    stop(
      sprintf(
        "Trial design mismatch: ncol(X)=%d but n_events=%d.",
        ncol(X),
        nrow(events)
      ),
      call. = FALSE
    )
  }

  out <- list(
    X = X,
    blockids = fmrihrf::blockids(sframe)
  )
  if (isTRUE(return_model)) {
    out$event_model <- emod
    out$sframe <- sframe
  }
  out
}

# Utility: shift trial regressors in time (seconds).
.bmk_shift_design <- function(X, tr, lag_sec = 0) {
  X <- as.matrix(X)
  lag_sec <- as.numeric(lag_sec)
  if (!is.finite(lag_sec) || abs(lag_sec) < 1e-12) return(X)

  tt <- seq(0, by = as.numeric(tr), length.out = nrow(X))
  out <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (j in seq_len(ncol(X))) {
    out[, j] <- stats::approx(tt, X[, j], xout = tt - lag_sec, rule = 2, ties = "ordered")$y
  }
  out
}

# Utility: broaden trial regressors by gaussian temporal smoothing.
.bmk_widen_design <- function(X, tr, widen_sec = 0) {
  X <- as.matrix(X)
  widen_sec <- as.numeric(widen_sec)
  if (!is.finite(widen_sec) || widen_sec <= 0) return(X)

  sd_samp <- widen_sec / as.numeric(tr)
  if (!is.finite(sd_samp) || sd_samp <= 0) return(X)

  radius <- max(1L, as.integer(ceiling(4 * sd_samp)))
  grid <- seq.int(-radius, radius)
  ker <- stats::dnorm(grid, mean = 0, sd = sd_samp)
  ker <- ker / sum(ker)

  out <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (j in seq_len(ncol(X))) {
    sm <- stats::filter(X[, j], filter = ker, sides = 2)
    sm <- as.numeric(sm)
    sm[!is.finite(sm)] <- 0
    out[, j] <- sm
  }
  out
}

# Utility: apply true-HRF mismatch transforms to generated design.
.bmk_apply_true_hrf_mismatch <- function(X, tr, lag_sec = 0, widen_sec = 0) {
  out <- .bmk_shift_design(X, tr = tr, lag_sec = lag_sec)
  out <- .bmk_widen_design(out, tr = tr, widen_sec = widen_sec)
  out
}

# Utility: apply per-trial mismatch vectors (lag/widen) to design columns.
.bmk_apply_true_hrf_mismatch_cols <- function(X, tr, lag_sec = 0, widen_sec = 0) {
  X <- as.matrix(X)
  p <- ncol(X)
  if (p < 1L) return(X)

  lag_vec <- if (length(lag_sec) == 1L) rep(as.numeric(lag_sec), p) else as.numeric(lag_sec)
  widen_vec <- if (length(widen_sec) == 1L) rep(as.numeric(widen_sec), p) else as.numeric(widen_sec)
  if (length(lag_vec) != p) stop("lag_sec must have length 1 or ncol(X).", call. = FALSE)
  if (length(widen_vec) != p) stop("widen_sec must have length 1 or ncol(X).", call. = FALSE)

  if (!any((abs(lag_vec) > 0) | (widen_vec > 0), na.rm = TRUE)) return(X)

  out <- matrix(0, nrow = nrow(X), ncol = p)
  for (j in seq_len(p)) {
    xj <- X[, j, drop = FALSE]
    lag_j <- lag_vec[j]
    widen_j <- widen_vec[j]
    if (is.finite(lag_j) && abs(lag_j) > 0) {
      xj <- .bmk_shift_design(xj, tr = tr, lag_sec = lag_j)
    }
    if (is.finite(widen_j) && widen_j > 0) {
      xj <- .bmk_widen_design(xj, tr = tr, widen_sec = widen_j)
    }
    out[, j] <- xj[, 1L]
  }
  out
}

# Utility: build true signal under voxel/trial HRF heterogeneity.
.bmk_make_true_signal <- function(
  X_true,
  beta_true,
  tr,
  true_hrf_lag_vox_sd = 0,
  true_hrf_widen_vox_sd = 0,
  true_hrf_lag_trial_sd = 0,
  true_hrf_widen_trial_sd = 0
) {
  X_true <- as.matrix(X_true)
  beta_true <- as.matrix(beta_true)

  n_trial <- ncol(X_true)
  n_vox <- ncol(beta_true)
  if (nrow(beta_true) != n_trial) {
    stop("nrow(beta_true) must equal ncol(X_true).", call. = FALSE)
  }

  lag_vox_sd <- max(0, as.numeric(true_hrf_lag_vox_sd))
  widen_vox_sd <- max(0, as.numeric(true_hrf_widen_vox_sd))
  lag_trial_sd <- max(0, as.numeric(true_hrf_lag_trial_sd))
  widen_trial_sd <- max(0, as.numeric(true_hrf_widen_trial_sd))

  hetero_on <- any(c(lag_vox_sd, widen_vox_sd, lag_trial_sd, widen_trial_sd) > 0)
  if (!hetero_on) {
    signal <- X_true %*% beta_true
    attr(signal, "hrf_heterogeneity_summary") <- c(
      mean_abs_lag = 0,
      max_abs_lag = 0,
      mean_widen = 0,
      max_widen = 0
    )
    return(signal)
  }

  lag_vox <- if (lag_vox_sd > 0) stats::rnorm(n_vox, mean = 0, sd = lag_vox_sd) else rep(0, n_vox)
  widen_vox <- if (widen_vox_sd > 0) stats::rnorm(n_vox, mean = 0, sd = widen_vox_sd) else rep(0, n_vox)

  signal <- matrix(0, nrow = nrow(X_true), ncol = n_vox)
  abs_lag_mean <- numeric(n_vox)
  abs_lag_max <- numeric(n_vox)
  widen_mean <- numeric(n_vox)
  widen_max <- numeric(n_vox)

  for (v in seq_len(n_vox)) {
    lag_trial <- if (lag_trial_sd > 0) stats::rnorm(n_trial, mean = 0, sd = lag_trial_sd) else rep(0, n_trial)
    widen_trial <- if (widen_trial_sd > 0) stats::rnorm(n_trial, mean = 0, sd = widen_trial_sd) else rep(0, n_trial)

    lag_vec <- lag_vox[v] + lag_trial
    widen_vec <- pmax(0, widen_vox[v] + widen_trial)

    Xv <- if (lag_trial_sd > 0 || widen_trial_sd > 0) {
      .bmk_apply_true_hrf_mismatch_cols(X_true, tr = tr, lag_sec = lag_vec, widen_sec = widen_vec)
    } else {
      .bmk_apply_true_hrf_mismatch(X_true, tr = tr, lag_sec = lag_vox[v], widen_sec = pmax(0, widen_vox[v]))
    }

    signal[, v] <- Xv %*% beta_true[, v]
    abs_lag_mean[v] <- mean(abs(lag_vec))
    abs_lag_max[v] <- max(abs(lag_vec))
    widen_mean[v] <- mean(widen_vec)
    widen_max[v] <- max(widen_vec)
  }

  attr(signal, "hrf_heterogeneity_summary") <- c(
    mean_abs_lag = mean(abs_lag_mean),
    max_abs_lag = max(abs_lag_max),
    mean_widen = mean(widen_mean),
    max_widen = max(widen_max)
  )
  signal
}

# Utility: identify exactly duplicated (or sign-flipped) design columns.
.bmk_find_duplicate_cols <- function(X, tol = 1e-10) {
  X <- as.matrix(X)
  p <- ncol(X)
  if (p < 2L) return(matrix(integer(0), ncol = 2L))

  pairs <- matrix(integer(0), ncol = 2L)
  for (i in seq_len(p - 1L)) {
    for (j in seq.int(i + 1L, p)) {
      d1 <- max(abs(X[, i] - X[, j]))
      d2 <- max(abs(X[, i] + X[, j]))
      if (is.finite(d1) && is.finite(d2) && (d1 <= tol || d2 <= tol)) {
        pairs <- rbind(pairs, c(i, j))
      }
    }
  }
  pairs
}

# Utility: rebuild design with tiny onset perturbations until no exact duplicates.
.bmk_resolve_duplicate_trials <- function(
  events,
  n_tp,
  tr,
  basis = "spmg1",
  precision = 0.3,
  return_model = FALSE,
  tol = 1e-10,
  max_iter = 12L,
  jitter_step = NULL
) {
  if (is.null(jitter_step)) jitter_step <- max(1e-3, precision * 0.75)
  jitter_step <- as.numeric(jitter_step)
  eps <- max(1e-4, jitter_step * 0.1)

  ev <- as.data.frame(events)
  n_tp <- as.integer(n_tp)
  tr <- as.numeric(tr)
  max_iter <- as.integer(max_iter)

  enforce_increasing <- function(df) {
    runs <- sort(unique(as.integer(df$run)))
    out <- df
    for (r in runs) {
      idx <- which(as.integer(out$run) == r)
      if (length(idx) < 1L) next
      ord <- order(out$onset[idx], out$trial_id[idx])
      id2 <- idx[ord]
      on <- as.numeric(out$onset[id2])
      if (length(on) > 1L) {
        for (k in 2:length(on)) on[k] <- max(on[k], on[k - 1L] + eps)
      }

      cap <- max(0, as.numeric(n_tp[r]) * tr - 2 * tr)
      if (on[length(on)] > cap) {
        on <- on - (on[length(on)] - cap)
      }
      if (on[1L] < 0) {
        on <- on - on[1L]
      }
      if (on[length(on)] > cap && length(on) > 1L) {
        span <- max(cap, eps * (length(on) - 1L))
        on <- seq(0, span, length.out = length(on))
      }
      if (length(on) > 1L) {
        for (k in 2:length(on)) on[k] <- max(on[k], on[k - 1L] + eps)
      }
      if (on[length(on)] > cap && on[length(on)] > 0) {
        on <- on * (cap / on[length(on)])
        if (length(on) > 1L) {
          for (k in 2:length(on)) on[k] <- max(on[k], on[k - 1L] + eps)
        }
      }

      out$onset[id2] <- as.numeric(on)
    }
    out <- out[order(out$run, out$onset, out$trial_id), , drop = FALSE]
    rownames(out) <- NULL
    out
  }

  for (it in seq_len(max_iter + 1L)) {
    ev <- enforce_increasing(ev)
    design <- .bmk_trial_design(
      events = ev,
      n_tp = n_tp,
      tr = tr,
      basis = basis,
      precision = precision,
      return_model = return_model
    )
    dup <- .bmk_find_duplicate_cols(design$X, tol = tol)
    if (nrow(dup) < 1L) {
      return(list(events = ev, design = design, n_fix = it - 1L, unresolved = FALSE))
    }
    if (it > max_iter) break

    cols_to_shift <- unique(dup[, 2L])
    for (k in seq_along(cols_to_shift)) {
      j <- cols_to_shift[k]
      run_j <- as.integer(ev$run[j])
      cap <- max(0, as.numeric(n_tp[run_j]) * tr - 2 * tr)
      step <- jitter_step * (1 + 0.05 * k + 0.05 * it)
      trial_new <- ev$onset[j] + step
      if (trial_new > cap) {
        trial_new <- max(0, ev$onset[j] - step)
      }
      ev$onset[j] <- as.numeric(trial_new)
    }
  }

  ev <- enforce_increasing(ev)
  design <- .bmk_trial_design(
    events = ev,
    n_tp = n_tp,
    tr = tr,
    basis = basis,
    precision = precision,
    return_model = return_model
  )
  list(events = ev, design = design, n_fix = max_iter, unresolved = TRUE)
}

# Utility: default SBHM library matrix from smooth decays on the TR grid.
.bmk_sbhm_library <- function(
  tgrid,
  lag_grid = seq(-3, 3, by = 1),
  width_scales = c(0.8, 1, 1.2, 1.5),
  gamma_shapes = c(4, 6, 8),
  gamma_rates = c(0.8, 1, 1.2)
) {
  tt <- as.numeric(tgrid) - min(as.numeric(tgrid))
  libs <- list()
  idx <- 1L

  for (ws in width_scales) {
    h0 <- fmrihrf::hrf_spmg1(tt / ws)
    h0 <- as.numeric(h0)
    h0[!is.finite(h0)] <- 0
    for (lag in lag_grid) {
      libs[[idx]] <- stats::approx(tt, h0, xout = tt - lag, rule = 2, ties = "ordered")$y
      idx <- idx + 1L
    }
  }

  for (sh in gamma_shapes) {
    for (rt in gamma_rates) {
      hg <- as.numeric(fmrihrf::hrf_gamma(tt, shape = sh, rate = rt))
      hg[!is.finite(hg)] <- 0
      libs[[idx]] <- hg
      idx <- idx + 1L
    }
  }

  H <- do.call(cbind, libs)
  H <- sweep(H, 2L, sqrt(pmax(colSums(H^2), .Machine$double.eps)), "/")
  H[!is.finite(H)] <- 0
  H
}

# Utility: scale design columns to unit L2 norm.
.bmk_scale_design <- function(X) {
  X <- as.matrix(X)
  scl <- sqrt(pmax(colSums(X^2), .Machine$double.eps))
  sweep(X, 2L, scl, "/")
}

# Utility: draw standardized innovation noise.
.bmk_draw_innov <- function(n, sd, dist = c("normal", "t"), df = 6) {
  dist <- match.arg(dist)
  sd <- as.numeric(sd)
  n <- as.integer(n)

  if (dist == "normal") {
    return(stats::rnorm(n, sd = sd))
  }

  df <- as.numeric(df)
  if (!is.finite(df) || df <= 2) {
    stop("noise_df must be > 2 when noise_innov_dist='t'.", call. = FALSE)
  }
  stats::rt(n, df = df) * sd * sqrt((df - 2) / df)
}

# Utility: simulate observation noise with optional AR/ARMA structure and spikes.
.bmk_sim_noise <- function(
  n_time,
  n_vox,
  run_id,
  noise_sd,
  noise_model = c("iid", "ar1", "ar2", "arma11"),
  ar_rho = 0.35,
  ar_rho_sd = 0,
  ar_rho_vox_sd = 0,
  ar_phi2 = 0.2,
  arma_theta = -0.3,
  noise_spike_prob = 0,
  noise_spike_scale = 0,
  noise_spike_shared = FALSE,
  noise_innov_dist = c("normal", "t"),
  noise_df = 6,
  spatial_rho = 0
) {
  noise_model <- match.arg(noise_model)
  noise_innov_dist <- match.arg(noise_innov_dist)

  n_time <- as.integer(n_time)
  n_vox <- as.integer(n_vox)

  simulate_arima <- function(innov, phi = numeric(0), theta = numeric(0)) {
    n <- length(innov)
    out <- numeric(n)
    p <- length(phi)
    q <- length(theta)
    if (n < 1L) return(out)
    for (t in seq_len(n)) {
      ar_part <- 0
      ma_part <- innov[t]
      if (p > 0L) {
        max_k <- min(p, t - 1L)
        if (max_k > 0L) {
          for (k in seq_len(max_k)) {
            ar_part <- ar_part + phi[k] * out[t - k]
          }
        }
      }
      if (q > 0L) {
        max_j <- min(q, t - 1L)
        if (max_j > 0L) {
          for (j in seq_len(max_j)) {
            ma_part <- ma_part + theta[j] * innov[t - j]
          }
        }
      }
      out[t] <- ar_part + ma_part
    }
    out
  }

  stabilize_ar2 <- function(phi1, phi2, margin = 0.98) {
    phi1 <- .bmk_clip(as.numeric(phi1), -margin, margin)
    phi2 <- .bmk_clip(as.numeric(phi2), -margin, margin)
    if (!is.finite(phi2) || phi2 <= -margin) phi2 <- -margin + 1e-4
    upper <- margin - phi2
    lower <- phi2 - margin
    phi1 <- max(lower, min(upper, phi1))
    c(phi1, phi2)
  }

  if (noise_model == "iid") {
    out <- matrix(
      .bmk_draw_innov(n_time * n_vox, sd = noise_sd, dist = noise_innov_dist, df = noise_df),
      nrow = n_time,
      ncol = n_vox
    )
  } else {
    base_rho <- .bmk_clip(as.numeric(ar_rho), -0.98, 0.98)
    rho_sd <- max(0, as.numeric(ar_rho_sd))
    rho_vox_sd <- max(0, as.numeric(ar_rho_vox_sd))
    phi2_base <- .bmk_clip(as.numeric(ar_phi2), -0.98, 0.98)
    theta_base <- .bmk_clip(as.numeric(arma_theta), -0.98, 0.98)

    runs <- as.integer(as.factor(run_id))
    out <- matrix(0, nrow = n_time, ncol = n_vox)

    for (r in sort(unique(runs))) {
      idx <- which(runs == r)
      nr <- length(idx)
      if (nr < 1L) next

      rho_r <- if (rho_sd > 0) stats::rnorm(1L, mean = base_rho, sd = rho_sd) else base_rho
      rho_r <- .bmk_clip(as.numeric(rho_r), -0.98, 0.98)

      for (v in seq_len(n_vox)) {
        rho_v <- if (rho_vox_sd > 0) rho_r + stats::rnorm(1L, mean = 0, sd = rho_vox_sd) else rho_r
        rho_v <- .bmk_clip(as.numeric(rho_v), -0.98, 0.98)

        innov <- .bmk_draw_innov(nr, sd = 1, dist = noise_innov_dist, df = noise_df)
        series <- switch(
          noise_model,
          ar1 = simulate_arima(innov, phi = c(rho_v), theta = numeric(0)),
          ar2 = {
            phi <- stabilize_ar2(rho_v, phi2_base)
            simulate_arima(innov, phi = phi, theta = numeric(0))
          },
          arma11 = simulate_arima(innov, phi = c(rho_v), theta = c(theta_base)),
          simulate_arima(innov, phi = c(rho_v), theta = numeric(0))
        )
        sd_series <- stats::sd(series)
        if (is.finite(sd_series) && sd_series > 1e-8) {
          series <- as.numeric(series) * as.numeric(noise_sd) / sd_series
        } else {
          series <- rep(0, nr)
        }
        out[idx, v] <- series
      }
    }
  }

  # Apply spatial (cross-voxel) correlation via exponential decay kernel.
  spatial_rho <- as.numeric(spatial_rho)
  if (is.finite(spatial_rho) && spatial_rho > 0 && n_vox > 1L) {
    dmat <- abs(outer(seq_len(n_vox), seq_len(n_vox), "-"))
    Sigma <- spatial_rho^dmat
    L <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (!is.null(L)) {
      out <- out %*% L
    }
  }

  spike_prob <- .bmk_clip(as.numeric(noise_spike_prob), 0, 1)
  spike_scale <- max(0, as.numeric(noise_spike_scale))
  if (isTRUE(spike_prob > 0) && isTRUE(spike_scale > 0)) {
    spike_sd <- as.numeric(noise_sd) * spike_scale
    if (isTRUE(noise_spike_shared)) {
      keep <- stats::runif(n_time) < spike_prob
      idx <- which(keep)
      if (length(idx) > 0L) {
        a <- .bmk_draw_innov(length(idx), sd = spike_sd, dist = noise_innov_dist, df = noise_df)
        out[idx, ] <- out[idx, ] + matrix(a, nrow = length(idx), ncol = n_vox)
      }
    } else {
      keep <- matrix(stats::runif(n_time * n_vox) < spike_prob, nrow = n_time, ncol = n_vox)
      nsp <- sum(keep)
      if (nsp > 0L) {
        out[keep] <- out[keep] + .bmk_draw_innov(nsp, sd = spike_sd, dist = noise_innov_dist, df = noise_df)
      }
    }
  }

  out
}

# Utility: sample per-trial amplitude modulation (mean ~= 1).
.bmk_trial_amplitudes <- function(
  n_trial,
  model = c("none", "lognormal", "gamma"),
  cv = 0,
  outlier_prob = 0,
  outlier_scale = 0
) {
  model <- match.arg(model)
  n_trial <- as.integer(n_trial)
  cv <- as.numeric(cv)

  if (n_trial < 1L || model == "none" || !is.finite(cv) || cv <= 0) {
    return(rep(1, max(0L, n_trial)))
  }

  amps <- switch(
    model,
    lognormal = {
      sdlog <- sqrt(log(cv^2 + 1))
      stats::rlnorm(n_trial, meanlog = -0.5 * sdlog^2, sdlog = sdlog)
    },
    gamma = {
      shape <- 1 / (cv^2)
      scale <- cv^2
      stats::rgamma(n_trial, shape = shape, scale = scale)
    },
    rep(1, n_trial)
  )

  prob <- .bmk_clip(as.numeric(outlier_prob), 0, 1)
  scale <- max(0, as.numeric(outlier_scale))
  if (isTRUE(prob > 0) && isTRUE(scale > 0)) {
    keep <- stats::runif(n_trial) < prob
    if (any(keep)) {
      amps[keep] <- amps[keep] * (1 + stats::rexp(sum(keep), rate = 1 / scale))
    }
  }

  amps / mean(amps)
}

# Utility: sample true trial betas under configurable distributions.
.bmk_sample_beta <- function(
  n_trial,
  n_vox,
  beta_dist = c("normal", "laplace", "t", "lognormal"),
  beta_scale = 1,
  beta_df = 4,
  beta_lognorm_sdlog = 0.6,
  beta_sparsity = 0
) {
  beta_dist <- match.arg(beta_dist)
  n <- as.integer(n_trial) * as.integer(n_vox)
  beta_scale <- as.numeric(beta_scale)

  if (!is.finite(beta_scale) || beta_scale <= 0) {
    stop("beta_scale must be a positive finite value.", call. = FALSE)
  }

  vals <- switch(
    beta_dist,
    normal = stats::rnorm(n, sd = beta_scale),
    laplace = {
      b <- beta_scale / sqrt(2)
      stats::rexp(n, rate = 1 / b) - stats::rexp(n, rate = 1 / b)
    },
    t = {
      beta_df <- as.numeric(beta_df)
      if (!is.finite(beta_df) || beta_df <= 2) beta_df <- 4
      stats::rt(n, df = beta_df) * beta_scale * sqrt((beta_df - 2) / beta_df)
    },
    lognormal = {
      sdlog <- as.numeric(beta_lognorm_sdlog)
      if (!is.finite(sdlog) || sdlog <= 0) sdlog <- 0.6
      draw <- stats::rlnorm(n, meanlog = -0.5 * sdlog^2, sdlog = sdlog)
      sign <- sample(c(-1, 1), size = n, replace = TRUE)
      draw * sign * beta_scale
    },
    stats::rnorm(n, sd = beta_scale)
  )

  out <- matrix(vals, nrow = as.integer(n_trial), ncol = as.integer(n_vox))

  sparsity <- as.numeric(beta_sparsity)
  if (is.finite(sparsity) && sparsity > 0) {
    sparsity <- max(0, min(1, sparsity))
    zero_mask <- matrix(stats::runif(length(out)) < sparsity, nrow = nrow(out), ncol = ncol(out))
    out[zero_mask] <- 0
  }

  out
}

# Utility: map pooled coefficients back to trial space when pooling is enabled.
.bmk_backtransform_trial <- function(theta, pool_fit = NULL) {
  theta <- as.matrix(theta)
  if (is.null(pool_fit) || !isTRUE(pool_fit$enabled)) return(theta)

  Tmat <- as.matrix(pool_fit$trial_to_original)
  if (nrow(theta) != ncol(Tmat)) {
    stop("Pooling back-transform dimension mismatch.", call. = FALSE)
  }
  as.matrix(Tmat %*% theta)
}

# Utility: extract trial coefficients at a specific lambda from glmnet fit.
# If pooling is enabled in stglmnet, coefficients are back-transformed to the
# original trial space before scoring.
.bmk_extract_trial_betas <- function(fit, ntrial, s, pool_fit = NULL) {
  b <- stats::coef(fit, s = s)
  if (is.list(b)) {
    theta <- matrix(NA_real_, nrow = ntrial, ncol = length(b))
    for (k in seq_along(b)) {
      theta[, k] <- as.numeric(b[[k]][1L + seq_len(ntrial), 1L, drop = TRUE])
    }
    .bmk_backtransform_trial(theta, pool_fit = pool_fit)
  } else {
    theta <- matrix(as.numeric(b[1L + seq_len(ntrial), 1L, drop = TRUE]), ncol = 1L)
    .bmk_backtransform_trial(theta, pool_fit = pool_fit)
  }
}

# Utility: correlation/RMSE and additional calibration diagnostics against known trial betas.
.bmk_score_betas <- function(beta_hat, beta_true) {
  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 1L) return(NA_real_)
    mean(x)
  }

  beta_hat <- as.matrix(beta_hat)
  beta_true <- as.matrix(beta_true)
  ncomp <- min(ncol(beta_hat), ncol(beta_true))
  if (ncomp < 1L) {
    return(list(
      cor = NA_real_,
      rmse = NA_real_,
      slope_bias = NA_real_,
      intercept_bias = NA_real_,
      sign_acc = NA_real_
    ))
  }

  bh <- beta_hat[, seq_len(ncomp), drop = FALSE]
  bt <- beta_true[, seq_len(ncomp), drop = FALSE]

  cors <- vapply(seq_len(ncomp), function(v) {
    suppressWarnings(stats::cor(bh[, v], bt[, v], use = "pairwise.complete.obs"))
  }, numeric(1))

  slopes <- rep(NA_real_, ncomp)
  intercepts <- rep(NA_real_, ncomp)
  sign_acc <- rep(NA_real_, ncomp)

  for (v in seq_len(ncomp)) {
    x <- as.numeric(bt[, v])
    y <- as.numeric(bh[, v])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2L) next

    x_ok <- x[ok]
    y_ok <- y[ok]
    vx <- stats::var(x_ok)
    if (is.finite(vx) && vx > 0) {
      slopes[v] <- stats::cov(x_ok, y_ok) / vx
      intercepts[v] <- mean(y_ok) - slopes[v] * mean(x_ok)
    }

    nz <- abs(x_ok) > 1e-8
    if (any(nz)) {
      sign_acc[v] <- mean(sign(y_ok[nz]) == sign(x_ok[nz]))
    }
  }

  list(
    cor = safe_mean(cors),
    rmse = sqrt(mean((bh - bt)^2)),
    slope_bias = safe_mean(slopes - 1),
    intercept_bias = safe_mean(intercepts),
    sign_acc = safe_mean(sign_acc),
    per_voxel_cor = as.numeric(cors)
  )
}

# Utility: generate condition-structured betas with discriminable condition means.
.bmk_sample_beta_conditions <- function(
  n_trial,
  n_vox,
  n_conditions = 3L,
  condition_effect = 1,
  beta_dist = "normal",
  beta_scale = 1,
  beta_df = 4,
  beta_lognorm_sdlog = 0.6,
  beta_sparsity = 0
) {
  n_conditions <- max(2L, as.integer(n_conditions))
  n_trial <- as.integer(n_trial)
  n_vox <- as.integer(n_vox)
  condition_labels <- rep(seq_len(n_conditions), length.out = n_trial)

  # Base trial variability
  beta_base <- .bmk_sample_beta(
    n_trial = n_trial,
    n_vox = n_vox,
    beta_dist = beta_dist,
    beta_scale = beta_scale,
    beta_df = beta_df,
    beta_lognorm_sdlog = beta_lognorm_sdlog,
    beta_sparsity = beta_sparsity
  )

  # Add condition-specific spatial patterns (each condition has a unique
  # voxel mean pattern, scaled by condition_effect)
  cond_patterns <- matrix(
    stats::rnorm(n_conditions * n_vox, sd = as.numeric(condition_effect)),
    nrow = n_conditions,
    ncol = n_vox
  )
  for (j in seq_len(n_trial)) {
    beta_base[j, ] <- beta_base[j, ] + cond_patterns[condition_labels[j], ]
  }
  list(beta = beta_base, condition_labels = condition_labels)
}

# Utility: leave-one-trial-out correlation-based classifier.
# For each held-out trial, compute mean pattern per condition from remaining
# trials, classify to nearest mean (by Pearson correlation), return accuracy.
.bmk_classif_accuracy <- function(beta_hat, condition_labels) {
  beta_hat <- as.matrix(beta_hat)
  n <- nrow(beta_hat)
  labels <- as.integer(condition_labels)
  if (n < 2L || length(unique(labels)) < 2L) return(NA_real_)
  if (length(labels) != n) return(NA_real_)

  correct <- 0L
  for (i in seq_len(n)) {
    train_idx <- seq_len(n)[-i]
    train_labels <- labels[train_idx]
    conds <- sort(unique(train_labels))
    means <- vapply(conds, function(c) {
      colMeans(beta_hat[train_idx[train_labels == c], , drop = FALSE])
    }, numeric(ncol(beta_hat)))
    # means is n_vox x n_conditions
    test_vec <- beta_hat[i, ]
    sims <- suppressWarnings(stats::cor(test_vec, means))
    if (any(is.finite(sims))) {
      pred <- conds[which.max(sims)]
      if (pred == labels[i]) correct <- correct + 1L
    }
  }
  correct / n
}

# Utility: summarize per-method benchmark results.
.bmk_method_summary <- function(metrics_long) {
  spl <- split(metrics_long, metrics_long$method)
  out <- do.call(rbind, lapply(names(spl), function(m) {
    d <- spl[[m]]
    data.frame(
      method = m,
      n = nrow(d),
      mean_cor = mean(d$cor, na.rm = TRUE),
      median_cor = stats::median(d$cor, na.rm = TRUE),
      sd_cor = stats::sd(d$cor, na.rm = TRUE),
      mean_rmse = mean(d$rmse, na.rm = TRUE),
      median_rmse = stats::median(d$rmse, na.rm = TRUE),
      sd_rmse = stats::sd(d$rmse, na.rm = TRUE),
      mean_slope_bias = mean(d$slope_bias, na.rm = TRUE),
      median_slope_bias = stats::median(d$slope_bias, na.rm = TRUE),
      sd_slope_bias = stats::sd(d$slope_bias, na.rm = TRUE),
      mean_abs_slope_bias = mean(abs(d$slope_bias), na.rm = TRUE),
      mean_intercept_bias = mean(d$intercept_bias, na.rm = TRUE),
      median_intercept_bias = stats::median(d$intercept_bias, na.rm = TRUE),
      sd_intercept_bias = stats::sd(d$intercept_bias, na.rm = TRUE),
      mean_abs_intercept_bias = mean(abs(d$intercept_bias), na.rm = TRUE),
      mean_sign_acc = mean(d$sign_acc, na.rm = TRUE),
      median_sign_acc = stats::median(d$sign_acc, na.rm = TRUE),
      sd_sign_acc = stats::sd(d$sign_acc, na.rm = TRUE),
      mean_classif_acc = mean(d$classif_acc, na.rm = TRUE),
      mean_elapsed_sec = mean(d$elapsed_sec, na.rm = TRUE),
      success_rate = mean(as.logical(d$ok), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

# Utility: summarize method performance by scenario.
.bmk_suite_method_summary <- function(metrics_long) {
  grp <- split(metrics_long, interaction(metrics_long$scenario, metrics_long$method, drop = TRUE))
  out <- do.call(rbind, lapply(grp, function(d) {
    data.frame(
      scenario = as.character(d$scenario[1]),
      method = as.character(d$method[1]),
      n = nrow(d),
      mean_cor = mean(d$cor, na.rm = TRUE),
      median_cor = stats::median(d$cor, na.rm = TRUE),
      sd_cor = stats::sd(d$cor, na.rm = TRUE),
      mean_rmse = mean(d$rmse, na.rm = TRUE),
      median_rmse = stats::median(d$rmse, na.rm = TRUE),
      sd_rmse = stats::sd(d$rmse, na.rm = TRUE),
      mean_slope_bias = mean(d$slope_bias, na.rm = TRUE),
      median_slope_bias = stats::median(d$slope_bias, na.rm = TRUE),
      sd_slope_bias = stats::sd(d$slope_bias, na.rm = TRUE),
      mean_abs_slope_bias = mean(abs(d$slope_bias), na.rm = TRUE),
      mean_intercept_bias = mean(d$intercept_bias, na.rm = TRUE),
      median_intercept_bias = stats::median(d$intercept_bias, na.rm = TRUE),
      sd_intercept_bias = stats::sd(d$intercept_bias, na.rm = TRUE),
      mean_abs_intercept_bias = mean(abs(d$intercept_bias), na.rm = TRUE),
      mean_sign_acc = mean(d$sign_acc, na.rm = TRUE),
      median_sign_acc = stats::median(d$sign_acc, na.rm = TRUE),
      sd_sign_acc = stats::sd(d$sign_acc, na.rm = TRUE),
      mean_classif_acc = mean(d$classif_acc, na.rm = TRUE),
      mean_elapsed_sec = mean(d$elapsed_sec, na.rm = TRUE),
      success_rate = mean(as.logical(d$ok), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out[order(out$scenario, out$mean_cor, decreasing = TRUE), , drop = FALSE]
}

# Utility: empirical percentile confidence interval.
.bmk_empirical_ci <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) < 1L) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 8))
}

# Utility: flatten nested parameter lists into a single-row metadata record.
.bmk_flatten_named <- function(x, prefix = NULL) {
  out <- list()
  if (is.null(x)) return(out)

  if (is.data.frame(x)) {
    x <- as.list(x)
  }

  if (is.list(x) && !inherits(x, c("POSIXct", "POSIXt", "Date"))) {
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("v", seq_along(x))
    for (i in seq_along(x)) {
      nm <- nms[i]
      key <- if (is.null(prefix) || !nzchar(prefix)) nm else paste(prefix, nm, sep = ".")
      out <- c(out, .bmk_flatten_named(x[[i]], prefix = key))
    }
    return(out)
  }

  key <- .bmk_or(prefix, "value")
  if (length(x) == 1L) {
    out[[key]] <- if (is.factor(x)) as.character(x) else x
  } else {
    out[[key]] <- paste(as.character(x), collapse = ",")
  }
  out
}

# Utility: convert a named list to a one-row data frame.
.bmk_row_df <- function(x) {
  if (!length(x)) return(data.frame())
  out <- lapply(x, function(v) if (is.null(v)) NA else v)
  as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

# Utility: bind data frames with heterogeneous columns by filling missing values.
.bmk_bind_rows <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(NULL)

  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(d) {
    missing <- setdiff(cols, names(d))
    if (length(missing)) {
      for (nm in missing) d[[nm]] <- NA
    }
    d <- d[, cols, drop = FALSE]
    rownames(d) <- NULL
    d
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Utility: attach flattened scenario metadata to a benchmark table.
.bmk_attach_scenario_metadata <- function(df, scenario_design, by = "scenario") {
  if (!is.data.frame(df) || nrow(df) < 1L) return(df)
  if (!is.data.frame(scenario_design) || nrow(scenario_design) < 1L) return(df)
  if (!by %in% names(df) || !by %in% names(scenario_design)) return(df)

  keep <- setdiff(names(scenario_design), by)
  if (!length(keep)) return(df)

  idx <- match(df[[by]], scenario_design[[by]])
  extra <- scenario_design[idx, keep, drop = FALSE]
  rownames(extra) <- NULL
  cbind(df, extra)
}

# Utility: rescale a metric to [0, 1], with optional inversion.
.bmk_rescale01 <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  fin <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(fin)) return(out)

  lo <- min(x[fin], na.rm = TRUE)
  hi <- max(x[fin], na.rm = TRUE)
  if (!is.finite(hi - lo) || (hi - lo) <= 0) {
    out[fin] <- 1
  } else if (isTRUE(higher_better)) {
    out[fin] <- (x[fin] - lo) / (hi - lo)
  } else {
    out[fin] <- (hi - x[fin]) / (hi - lo)
  }
  out
}

# Utility: recommendation profiles for choosing methods under different priorities.
.bmk_recommend_profiles <- function() {
  list(
    balanced = c(
      recovery = 0.35,
      rmse = 0.15,
      sign = 0.10,
      bias = 0.15,
      stability = 0.15,
      speed = 0.05,
      consistency = 0.05
    ),
    accuracy = c(
      recovery = 0.45,
      rmse = 0.20,
      sign = 0.10,
      bias = 0.15,
      stability = 0.05,
      consistency = 0.05
    ),
    stability = c(
      stability = 0.35,
      recovery = 0.25,
      consistency = 0.15,
      rmse = 0.10,
      bias = 0.10,
      speed = 0.05
    ),
    speed = c(
      speed = 0.35,
      recovery = 0.25,
      stability = 0.15,
      rmse = 0.10,
      bias = 0.10,
      sign = 0.05
    ),
    classification = c(
      classification = 0.35,
      recovery = 0.25,
      rmse = 0.10,
      stability = 0.10,
      bias = 0.10,
      speed = 0.10
    )
  )
}

# Utility: explain why a method was recommended within a scenario.
.bmk_recommendation_reason <- function(top, runner_up = NULL) {
  if (is.null(runner_up) || nrow(runner_up) < 1L) {
    return("Only available method in this scenario.")
  }

  reasons <- character(0)
  if (is.finite(top$mean_cor - runner_up$mean_cor) && (top$mean_cor - runner_up$mean_cor) > 0.02) {
    reasons <- c(reasons, "best recovery")
  }
  if (is.finite(runner_up$mean_elapsed_sec) && is.finite(top$mean_elapsed_sec) &&
      top$mean_elapsed_sec < (0.8 * runner_up$mean_elapsed_sec)) {
    reasons <- c(reasons, "materially faster")
  }
  if (is.finite(top$success_rate - runner_up$success_rate) && (top$success_rate - runner_up$success_rate) > 0.05) {
    reasons <- c(reasons, "more stable")
  }

  top_bias <- mean(c(top$mean_abs_slope_bias, top$mean_abs_intercept_bias), na.rm = TRUE)
  runner_bias <- mean(c(runner_up$mean_abs_slope_bias, runner_up$mean_abs_intercept_bias), na.rm = TRUE)
  if (is.finite(top_bias) && is.finite(runner_bias) && top_bias < (0.9 * runner_bias)) {
    reasons <- c(reasons, "lower bias")
  }

  if (!length(reasons)) {
    reasons <- "best overall tradeoff"
  }

  paste(reasons, collapse = ", ")
}

# Utility: scenario-by-profile method recommendations from suite summaries.
.bmk_recommend_methods <- function(scenario_method_summary, profiles = .bmk_recommend_profiles()) {
  if (!is.data.frame(scenario_method_summary) || nrow(scenario_method_summary) < 1L) {
    return(data.frame())
  }

  spl <- split(scenario_method_summary, scenario_method_summary$scenario)
  rows <- unlist(lapply(names(spl), function(sc) {
    d <- spl[[sc]]
    d$score_recovery <- .bmk_rescale01(d$mean_cor, higher_better = TRUE)
    d$score_rmse <- .bmk_rescale01(d$mean_rmse, higher_better = FALSE)
    d$score_sign <- .bmk_rescale01(d$mean_sign_acc, higher_better = TRUE)
    d$score_slope_bias <- .bmk_rescale01(d$mean_abs_slope_bias, higher_better = FALSE)
    d$score_intercept_bias <- .bmk_rescale01(d$mean_abs_intercept_bias, higher_better = FALSE)
    d$score_bias <- rowMeans(d[, c("score_slope_bias", "score_intercept_bias"), drop = FALSE], na.rm = TRUE)
    d$score_stability <- .bmk_rescale01(d$success_rate, higher_better = TRUE)
    d$score_speed <- .bmk_rescale01(d$mean_elapsed_sec, higher_better = FALSE)
    d$score_consistency <- .bmk_rescale01(d$sd_cor, higher_better = FALSE)
    d$score_classification <- .bmk_rescale01(d$mean_classif_acc, higher_better = TRUE)

    lapply(names(profiles), function(profile_name) {
      weights <- profiles[[profile_name]]
      score_cols <- paste0("score_", names(weights))
      weight_vals <- as.numeric(weights)
      score_mat <- as.matrix(d[, score_cols, drop = FALSE])

      overall <- vapply(seq_len(nrow(score_mat)), function(i) {
        ok <- is.finite(score_mat[i, ]) & is.finite(weight_vals)
        if (!any(ok)) return(NA_real_)
        sum(score_mat[i, ok] * weight_vals[ok]) / sum(weight_vals[ok])
      }, numeric(1))

      d$overall_score <- overall
      ord <- order(-d$overall_score, -d$mean_cor, d$mean_rmse, d$mean_elapsed_sec, na.last = TRUE)
      top <- d[ord[1L], , drop = FALSE]
      runner_up <- if (nrow(d) > 1L) d[ord[2L], , drop = FALSE] else NULL
      score_margin <- if (!is.null(runner_up) && nrow(runner_up) > 0L) {
        top$overall_score - runner_up$overall_score
      } else {
        NA_real_
      }

      data.frame(
        scenario = sc,
        profile = profile_name,
        recommended_method = as.character(top$method[1]),
        runner_up_method = if (!is.null(runner_up) && nrow(runner_up) > 0L) as.character(runner_up$method[1]) else NA_character_,
        overall_score = as.numeric(top$overall_score[1]),
        runner_up_score = if (!is.null(runner_up) && nrow(runner_up) > 0L) as.numeric(runner_up$overall_score[1]) else NA_real_,
        score_margin = as.numeric(score_margin[1]),
        mean_cor = as.numeric(top$mean_cor[1]),
        mean_rmse = as.numeric(top$mean_rmse[1]),
        mean_sign_acc = as.numeric(top$mean_sign_acc[1]),
        mean_elapsed_sec = as.numeric(top$mean_elapsed_sec[1]),
        success_rate = as.numeric(top$success_rate[1]),
        reason = .bmk_recommendation_reason(top, runner_up),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Utility: aggregate recommendation winners across profiles.
.bmk_recommendation_summary <- function(recommendations) {
  if (!is.data.frame(recommendations) || nrow(recommendations) < 1L) {
    return(data.frame())
  }

  grp <- split(recommendations, interaction(recommendations$profile, recommendations$recommended_method, drop = TRUE))
  out <- do.call(rbind, lapply(grp, function(d) {
    data.frame(
      profile = as.character(d$profile[1]),
      method = as.character(d$recommended_method[1]),
      n_scenarios = nrow(d),
      mean_score = mean(d$overall_score, na.rm = TRUE),
      mean_margin = mean(d$score_margin, na.rm = TRUE),
      mean_cor = mean(d$mean_cor, na.rm = TRUE),
      mean_rmse = mean(d$mean_rmse, na.rm = TRUE),
      mean_elapsed_sec = mean(d$mean_elapsed_sec, na.rm = TRUE),
      success_rate = mean(d$success_rate, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$profile, out$n_scenarios, out$mean_score, decreasing = TRUE), , drop = FALSE]
}

# Utility: resolve lss() from a loaded development session first.
.bmk_get_lss <- function() {
  if (exists("lss", mode = "function", inherits = TRUE)) {
    return(get("lss", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("fmrilss", quietly = TRUE)) {
    return(get("lss", envir = asNamespace("fmrilss")))
  }
  stop("Could not find lss(). Load or install 'fmrilss' first.", call. = FALSE)
}

# Utility: resolve sbhm_build() from loaded session or namespace.
.bmk_get_sbhm_build <- function() {
  if (exists("sbhm_build", mode = "function", inherits = TRUE)) {
    return(get("sbhm_build", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("fmrilss", quietly = TRUE)) {
    return(get("sbhm_build", envir = asNamespace("fmrilss")))
  }
  stop("Could not find sbhm_build(). Load or install 'fmrilss' first.", call. = FALSE)
}

# Utility: default stglmnet benchmark method library for the internal fmrilss
# backend.
.bmk_default_st_configs <- function(alpha, lambda, overlap_adaptive = TRUE) {
  st_mult <- if (isTRUE(overlap_adaptive)) "multiplicative" else "none"
  list(
    st_pool_compcv_wnever = list(
      overlap_strategy = "none",
      pool_to_mean = TRUE,
      pool_strength = 1.2,
      pool_mean_penalty = 0.2,
      alpha = alpha,
      lambda = lambda,
      whiten = "never",
      cv_args = list(
        type.measure = "composite",
        composite_weights = c(mse = 0.45, correlation = 0.25, reliability = 0.30)
      ),
      standardize = FALSE,
      intercept = FALSE
    ),
    st_pool_compcv = list(
      overlap_strategy = "none",
      pool_to_mean = TRUE,
      pool_strength = 1.2,
      pool_mean_penalty = 0.2,
      alpha = alpha,
      lambda = lambda,
      cv_args = list(
        type.measure = "composite",
        composite_weights = c(mse = 0.45, correlation = 0.25, reliability = 0.30)
      ),
      standardize = FALSE,
      intercept = FALSE
    ),
    st_pool_joint_mult = list(
      overlap_strategy = st_mult,
      overlap_exponent = 2,
      overlap_strength = 1,
      pool_to_mean = TRUE,
      pool_strength = 1.2,
      pool_mean_penalty = 0.2,
      pool_scale_by_overlap = TRUE,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      intercept = FALSE
    ),
    st_pool_joint_mult_1se = list(
      overlap_strategy = st_mult,
      overlap_exponent = 2,
      overlap_strength = 1,
      pool_to_mean = TRUE,
      pool_strength = 1.2,
      pool_mean_penalty = 0.2,
      pool_scale_by_overlap = TRUE,
      alpha = alpha,
      lambda = lambda,
      cv_select = "1se",
      standardize = FALSE,
      intercept = FALSE
    ),
    st_pool_compcv_wB = list(
      overlap_strategy = "none",
      pool_to_mean = TRUE,
      pool_strength = 1.2,
      pool_mean_penalty = 0.2,
      alpha = alpha,
      lambda = lambda,
      cv_args = list(
        type.measure = "composite",
        composite_weights = c(mse = 0.25, correlation = 0.20, reliability = 0.55)
      ),
      standardize = FALSE,
      intercept = FALSE
    ),
    st_pool_mean = list(
      overlap_strategy = "none",
      pool_to_mean = TRUE,
      pool_strength = 1.2,
      pool_mean_penalty = 0.2,
      pool_scale_by_overlap = TRUE,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      intercept = FALSE
    )
  )
}

# Utility: resolve lss_sbhm_design() from loaded session or namespace.
.bmk_get_lss_sbhm_design <- function() {
  if (exists("lss_sbhm_design", mode = "function", inherits = TRUE)) {
    return(get("lss_sbhm_design", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("fmrilss", quietly = TRUE)) {
    return(get("lss_sbhm_design", envir = asNamespace("fmrilss")))
  }
  stop("Could not find lss_sbhm_design(). Load or install 'fmrilss' first.", call. = FALSE)
}

# Utility: resolve hrfals lss_mode_a() from loaded session or namespace.
.bmk_get_hrfals_lss_mode_a <- function() {
  if (exists("lss_mode_a", mode = "function", inherits = TRUE)) {
    return(get("lss_mode_a", mode = "function", inherits = TRUE))
  }
  if (exists("fastlss_shared", mode = "function", inherits = TRUE)) {
    return(get("fastlss_shared", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("hrfals", quietly = TRUE)) {
    ns <- asNamespace("hrfals")
    if (exists("lss_mode_a", mode = "function", envir = ns, inherits = FALSE)) {
      return(get("lss_mode_a", envir = ns))
    }
    if (exists("fastlss_shared", mode = "function", envir = ns, inherits = FALSE)) {
      return(get("fastlss_shared", envir = ns))
    }
  }
  stop("Could not find hrfals lss_mode_a()/fastlss_shared(). Load or install 'hrfals' first.", call. = FALSE)
}

# Utility: resolve hrfals create_cfals_design() from loaded session or namespace.
.bmk_get_hrfals_create_cfals_design <- function() {
  if (exists("create_cfals_design", mode = "function", inherits = TRUE)) {
    return(get("create_cfals_design", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("hrfals", quietly = TRUE)) {
    ns <- asNamespace("hrfals")
    if (exists("create_cfals_design", mode = "function", envir = ns, inherits = FALSE)) {
      return(get("create_cfals_design", envir = ns))
    }
  }
  stop("Could not find hrfals create_cfals_design(). Load or install 'hrfals' first.", call. = FALSE)
}

# Utility: resolve hrfals hrfals_from_design() from loaded session or namespace.
.bmk_get_hrfals_from_design <- function() {
  if (exists("hrfals_from_design", mode = "function", inherits = TRUE)) {
    return(get("hrfals_from_design", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("hrfals", quietly = TRUE)) {
    ns <- asNamespace("hrfals")
    if (exists("hrfals_from_design", mode = "function", envir = ns, inherits = FALSE)) {
      return(get("hrfals_from_design", envir = ns))
    }
  }
  stop("Could not find hrfals hrfals_from_design(). Load or install 'hrfals' first.", call. = FALSE)
}

# Utility: resolve hrfals HRF basis from object or basis name.
.bmk_resolve_hrfals_basis <- function(hrfals_basis) {
  if (inherits(hrfals_basis, "HRF")) return(hrfals_basis)
  if (is.character(hrfals_basis) && length(hrfals_basis) >= 1L) {
    nm <- as.character(hrfals_basis)[1L]
    ns <- asNamespace("fmrihrf")
    cands <- unique(c(
      nm,
      toupper(nm),
      paste0("HRF_", toupper(nm))
    ))
    for (cand in cands) {
      if (exists(cand, envir = ns, mode = "function", inherits = FALSE)) {
        obj <- get(cand, envir = ns)
        if (inherits(obj, "HRF")) return(obj)
      }
    }
    obj <- tryCatch(fmrihrf::make_hrf(nm), error = function(e) NULL)
    if (!is.null(obj) && inherits(obj, "HRF")) return(obj)
    stop("Could not construct hrfals_basis from name: ", nm, call. = FALSE)
  }
  stop("hrfals_basis must be an HRF object or basis name string.", call. = FALSE)
}

# Utility: normalize scenario definitions to a list-of-lists.
.bmk_normalize_scenarios <- function(scenarios) {
  if (is.null(scenarios)) {
    stop("scenarios must be a data.frame or list of scenario lists.", call. = FALSE)
  }

  if (is.data.frame(scenarios)) {
    idx <- seq_len(nrow(scenarios))
    out <- lapply(idx, function(i) {
      vals <- lapply(names(scenarios), function(nm) {
        col <- scenarios[[nm]]
        if (is.list(col)) col[[i]] else col[i]
      })
      stats::setNames(vals, names(scenarios))
    })
    return(out)
  }

  if (!is.list(scenarios) || length(scenarios) < 1L) {
    stop("scenarios must contain at least one scenario.", call. = FALSE)
  }

  if (all(vapply(scenarios, is.list, logical(1)))) {
    return(scenarios)
  }

  list(scenarios)
}

# Preset scenario collections for quick and wide sweeps.
default_benchmark_scenarios <- function(kind = c("quick", "wide")) {
  kind <- match.arg(kind)

  quick <- list(
    baseline_iid = list(name = "baseline_iid"),
    ar1_moderate = list(name = "ar1_moderate", noise_model = "ar1", ar_rho = 0.35),
    isi_exponential = list(name = "isi_exponential", isi_model = "exponential", isi = 3.0),
    hrf_mismatch = list(name = "hrf_mismatch", basis = "spmg1", true_basis = "spmg2")
  )

  if (kind == "quick") return(quick)

  c(
    quick,
    list(
      ar1_prewhitened = list(
        name = "ar1_prewhitened",
        noise_model = "ar1",
        ar_rho = 0.5,
        lss_prewhiten = list(method = "ar", order = 1)
      ),
      ar1_high = list(name = "ar1_high", noise_model = "ar1", ar_rho = 0.6),
      ar1_heterogeneous = list(
        name = "ar1_heterogeneous",
        noise_model = "ar1",
        ar_rho = 0.4,
        ar_rho_sd = 0.18
      ),
      ar2_moderate = list(
        name = "ar2_moderate",
        noise_model = "ar2",
        ar_rho = 0.5,
        ar_phi2 = 0.2
      ),
      ar2_heterogeneous_spiky = list(
        name = "ar2_heterogeneous_spiky",
        noise_model = "ar2",
        ar_rho = 0.45,
        ar_phi2 = 0.3,
        ar_rho_sd = 0.15,
        ar_rho_vox_sd = 0.2,
        noise_spike_prob = 0.01,
        noise_spike_scale = 8,
        noise_spike_shared = TRUE
      ),
      arma11_persistent = list(
        name = "arma11_persistent",
        noise_model = "arma11",
        ar_rho = 0.6,
        arma_theta = -0.45
      ),
      isi_uniform_dense = list(
        name = "isi_uniform_dense",
        isi_model = "uniform",
        isi = 1.0,
        isi_params = list(width = 0.4)
      ),
      isi_lognormal_heavy_tail = list(
        name = "isi_lognormal_heavy_tail",
        isi_model = "lognormal",
        isi = 6.0,
        isi_params = list(sdlog = 0.8)
      ),
      amplitude_lognormal_cv05 = list(
        name = "amplitude_lognormal_cv05",
        trial_amplitude_model = "lognormal",
        trial_amplitude_cv = 0.5
      ),
      amplitude_gamma_cv1 = list(
        name = "amplitude_gamma_cv1",
        trial_amplitude_model = "gamma",
        trial_amplitude_cv = 1.0
      ),
      amplitude_lognormal_cv1p5 = list(
        name = "amplitude_lognormal_cv1p5",
        trial_amplitude_model = "lognormal",
        trial_amplitude_cv = 1.5
      ),
      amplitude_gamma_cv2 = list(
        name = "amplitude_gamma_cv2",
        trial_amplitude_model = "gamma",
        trial_amplitude_cv = 2.0
      ),
      amplitude_lognormal_outlier = list(
        name = "amplitude_lognormal_outlier",
        trial_amplitude_model = "lognormal",
        trial_amplitude_cv = 1.2,
        trial_amplitude_outlier_prob = 0.06,
        trial_amplitude_outlier_scale = 6
      ),
      sparse_laplace_betas = list(
        name = "sparse_laplace_betas",
        beta_dist = "laplace",
        beta_sparsity = 0.5
      ),
      hrf_timing_shift = list(
        name = "hrf_timing_shift",
        true_onset_shift_sd = 1.0
      ),
      hrf_lag_1p5 = list(
        name = "hrf_lag_1p5",
        true_hrf_lag_sec = 1.5
      ),
      hrf_lag_3 = list(
        name = "hrf_lag_3",
        true_hrf_lag_sec = 3.0
      ),
      hrf_widen_1p5 = list(
        name = "hrf_widen_1p5",
        true_hrf_widen_sec = 1.5
      ),
      hrf_shift_and_widen = list(
        name = "hrf_shift_and_widen",
        true_hrf_lag_sec = 2.0,
        true_hrf_widen_sec = 1.5
      ),
      hrf_extreme_lag_4_widen_2 = list(
        name = "hrf_extreme_lag_4_widen_2",
        true_hrf_lag_sec = 4.0,
        true_hrf_widen_sec = 2.0
      ),
      hrf_extreme_basis_shift = list(
        name = "hrf_extreme_basis_shift",
        true_basis = "spmg2",
        true_hrf_lag_sec = 3.0,
        true_hrf_widen_sec = 2.5
      ),
      hrf_voxel_lag_var = list(
        name = "hrf_voxel_lag_var",
        true_hrf_lag_vox_sd = 0.8
      ),
      hrf_voxel_widen_var = list(
        name = "hrf_voxel_widen_var",
        true_hrf_widen_vox_sd = 0.8
      ),
      hrf_trial_lag_var = list(
        name = "hrf_trial_lag_var",
        true_hrf_lag_trial_sd = 0.5
      ),
      hrf_vox_trial_mismatch = list(
        name = "hrf_vox_trial_mismatch",
        true_hrf_lag_sec = 1.5,
        true_hrf_widen_sec = 1.0,
        true_hrf_lag_vox_sd = 0.6,
        true_hrf_widen_vox_sd = 0.6,
        true_hrf_lag_trial_sd = 0.4,
        true_hrf_widen_trial_sd = 0.4
      ),
      snr_high = list(
        name = "snr_high",
        noise_sd = 0.6
      ),
      snr_low = list(
        name = "snr_low",
        noise_sd = 2.4
      ),
      snr_very_low = list(
        name = "snr_very_low",
        noise_sd = 4.0
      ),
      tr_hcp = list(
        name = "tr_hcp",
        tr = 0.72,
        n_tp = 200L
      ),
      tr_standard = list(
        name = "tr_standard",
        tr = 2.0,
        n_tp = 80L
      ),
      combined_hard = list(
        name = "combined_hard",
        noise_model = "arma11",
        ar_rho = 0.55,
        arma_theta = -0.35,
        isi_model = "lognormal",
        isi = 1.1,
        isi_params = list(sdlog = 0.9),
        trial_amplitude_model = "lognormal",
        trial_amplitude_cv = 0.8,
        trial_amplitude_outlier_prob = 0.04,
        trial_amplitude_outlier_scale = 5,
        beta_dist = "laplace",
        beta_sparsity = 0.45,
        true_basis = "spmg2",
        true_onset_shift_sd = 0.8,
        noise_spike_prob = 0.008,
        noise_spike_scale = 7
      )
    )
  )
}

#' Benchmark LSS methods on single-trial beta recovery.
#'
#' This benchmark now supports controlled scenario stressors including:
#' - AR/ARMA noise and innovation tails,
#' - outlier spike contamination,
#' - ISI distributions,
#' - trial amplitude variability,
#' - HRF basis/timing mismatch between data generation and estimation,
#' - voxel/trial HRF heterogeneity in lag and broadening.
benchmark_lss_methods <- function(
  n_reps = 30,
  seed = 20260227L,
  n_tp = 150L,
  n_trials = 40L,
  n_vox = 24L,
  isi = 1.2,
  isi_model = c("spread", "fixed", "uniform", "exponential", "gamma", "lognormal"),
  isi_params = list(),
  tr = 1,
  noise_sd = 1.2,
  noise_model = c("iid", "ar1", "ar2", "arma11"),
  ar_rho = 0.35,
  ar_rho_sd = 0,
  ar_rho_vox_sd = 0,
  ar_phi2 = 0.2,
  arma_theta = -0.3,
  noise_spike_prob = 0,
  noise_spike_scale = 0,
  noise_spike_shared = FALSE,
  noise_innov_dist = c("normal", "t"),
  noise_df = 6,
  spatial_rho = 0,
  beta_sparsity = 0,
  beta_dist = c("normal", "laplace", "t", "lognormal"),
  beta_scale = 1,
  beta_df = 4,
  beta_lognorm_sdlog = 0.6,
  n_conditions = 0L,
  condition_effect = 1,
  trial_amplitude_model = c("none", "lognormal", "gamma"),
  trial_amplitude_cv = 0,
  trial_amplitude_outlier_prob = 0,
  trial_amplitude_outlier_scale = 0,
  basis = "spmg1",
  true_basis = basis,
  onset_jitter_sd = 0,
  true_onset_shift_sd = 0,
  true_hrf_lag_sec = 0,
  true_hrf_widen_sec = 0,
  true_hrf_lag_vox_sd = 0,
  true_hrf_widen_vox_sd = 0,
  true_hrf_lag_trial_sd = 0,
  true_hrf_widen_trial_sd = 0,
  avoid_duplicate_regressors = TRUE,
  duplicate_tol = 1e-10,
  duplicate_max_iter = 12L,
  duplicate_jitter_step = NULL,
  alpha = 0.2,
  lambda = exp(seq(log(0.2), log(0.005), length.out = 25)),
  overlap_adaptive = TRUE,
  st_configs = NULL,
  primary_st = NULL,
  cv_folds = 5L,
  cv_type.measure = c("auto", "mse", "correlation", "reliability", "composite"),
  lss_method = "r_optimized",
  lss_intercept = c("runwise", "global", "none"),
  st_match_lss_nuisance = TRUE,
  lss_prewhiten = NULL,
  include_oasis = FALSE,
  oasis_grid = NULL,
  oasis_prewhiten = lss_prewhiten,
  include_hrfals = FALSE,
  hrfals_modes = c("r", "cpp"),
  hrfals_lambda_ridge = 0,
  hrfals_fit_methods = c("ls_svd_only", "ls_svd_1als", "cf_als"),
  hrfals_basis = "spmg2",
  hrfals_control = list(),
  include_sbhm = FALSE,
  sbhm_rank = 8L,
  sbhm_library_H = NULL,
  sbhm_shifts = c(0.5, 1.0),
  sbhm_match = list(
    shrink = list(tau = 0.15, ref = NULL, snr = NULL),
    topK = 3,
    soft_blend = TRUE,
    blend_margin = 0.08,
    whiten = FALSE,
    sv_floor_rel = 0.05,
    whiten_power = 0.5,
    min_margin = 0.08,
    min_beta_norm = 1e-3,
    orient_ref = TRUE,
    alpha_source = "prepass",
    rank1_min = 0
  ),
  sbhm_oasis = list(
    ridge_mode = "fractional",
    ridge_x = 0.05,
    ridge_b = 0.05
  ),
  sbhm_amplitude = list(
    method = "lss1",
    ridge = list(mode = "fractional", lambda = 0.02),
    ridge_frac = list(x = 0.02, b = 0.02),
    cond_gate = NULL,
    adaptive = list(enable = FALSE, base = 0.02, k0 = 1000, max = 0.08),
    return_se = FALSE
  ),
  sbhm_prewhiten = lss_prewhiten,
  n_runs = 1L,
  resample_events = FALSE,
  fast_mode = FALSE,
  trace = FALSE
) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) {
    stop("Package 'fmridesign' is required for benchmark_lss_methods().", call. = FALSE)
  }
  if (!requireNamespace("fmrihrf", quietly = TRUE)) {
    stop("Package 'fmrihrf' is required for benchmark_lss_methods().", call. = FALSE)
  }
  lss_fn <- .bmk_get_lss()

  if (isTRUE(fast_mode)) {
    if (missing(n_reps) || identical(n_reps, formals()$n_reps)) n_reps <- 2L
    if (missing(n_vox) || identical(n_vox, formals()$n_vox)) n_vox <- 3L
    if (missing(n_trials) || identical(n_trials, formals()$n_trials)) n_trials <- 10L
    if (missing(n_tp) || identical(n_tp, formals()$n_tp)) n_tp <- 60L
  }

  hrfals_lss_mode_a_fn <- if (isTRUE(include_hrfals)) .bmk_get_hrfals_lss_mode_a() else NULL
  hrfals_create_cfals_design_fn <- if (isTRUE(include_hrfals)) .bmk_get_hrfals_create_cfals_design() else NULL
  hrfals_from_design_fn <- if (isTRUE(include_hrfals)) .bmk_get_hrfals_from_design() else NULL
  sbhm_build_fn <- if (isTRUE(include_sbhm)) .bmk_get_sbhm_build() else NULL
  lss_sbhm_design_fn <- if (isTRUE(include_sbhm)) .bmk_get_lss_sbhm_design() else NULL

  noise_model <- match.arg(noise_model)
  noise_innov_dist <- match.arg(noise_innov_dist)
  cv_type.measure <- match.arg(cv_type.measure)
  lss_intercept <- match.arg(lss_intercept)
  st_match_lss_nuisance <- isTRUE(st_match_lss_nuisance)
  if (isTRUE(include_hrfals)) {
    hrfals_modes <- unique(match.arg(hrfals_modes, c("r", "cpp"), several.ok = TRUE))
    hrfals_lambda_ridge <- max(0, as.numeric(hrfals_lambda_ridge))
    hrfals_fit_methods <- unique(match.arg(
      hrfals_fit_methods,
      c("ls_svd_only", "ls_svd_1als", "cf_als"),
      several.ok = TRUE
    ))
    hrfals_basis <- .bmk_resolve_hrfals_basis(hrfals_basis)
    if (fmrihrf::nbasis(hrfals_basis) <= 1L && length(hrfals_fit_methods) > 0L) {
      stop("hrfals_fit_methods require a multi-basis hrfals_basis (nbasis > 1).", call. = FALSE)
    }
  } else {
    hrfals_modes <- character(0)
    hrfals_fit_methods <- character(0)
  }
  isi_model <- match.arg(isi_model)
  beta_dist <- match.arg(beta_dist)
  trial_amplitude_model <- match.arg(trial_amplitude_model)
  ar_rho_vox_sd <- max(0, as.numeric(ar_rho_vox_sd))
  ar_phi2 <- .bmk_clip(as.numeric(ar_phi2), -0.98, 0.98)
  arma_theta <- .bmk_clip(as.numeric(arma_theta), -0.98, 0.98)
  noise_spike_prob <- .bmk_clip(as.numeric(noise_spike_prob), 0, 1)
  noise_spike_scale <- max(0, as.numeric(noise_spike_scale))
  noise_spike_shared <- isTRUE(noise_spike_shared)
  trial_amplitude_outlier_prob <- .bmk_clip(as.numeric(trial_amplitude_outlier_prob), 0, 1)
  trial_amplitude_outlier_scale <- max(0, as.numeric(trial_amplitude_outlier_scale))
  true_hrf_lag_vox_sd <- max(0, as.numeric(true_hrf_lag_vox_sd))
  true_hrf_widen_vox_sd <- max(0, as.numeric(true_hrf_widen_vox_sd))
  true_hrf_lag_trial_sd <- max(0, as.numeric(true_hrf_lag_trial_sd))
  true_hrf_widen_trial_sd <- max(0, as.numeric(true_hrf_widen_trial_sd))
  lss_prewhiten <- .bmk_normalize_prewhiten(lss_prewhiten)
  oasis_prewhiten <- .bmk_normalize_prewhiten(oasis_prewhiten)
  sbhm_prewhiten <- .bmk_normalize_prewhiten(sbhm_prewhiten)

  n_runs <- as.integer(n_runs)
  if (length(n_tp) == 1L) {
    n_tp <- rep(as.integer(n_tp), n_runs)
  } else {
    n_tp <- as.integer(n_tp)
    if (length(n_tp) != n_runs) {
      stop(sprintf("Expected %d run lengths in n_tp, got %d.", n_runs, length(n_tp)), call. = FALSE)
    }
  }
  need_event_model <- isTRUE(include_sbhm) || (isTRUE(include_hrfals) && length(hrfals_fit_methods) > 0L)

  events <- .bmk_make_events(
    n_trials = as.integer(n_trials),
    n_runs = n_runs,
    n_tp = n_tp,
    tr = tr,
    isi = isi,
    isi_model = isi_model,
    isi_params = isi_params,
    jitter_sd = onset_jitter_sd
  )
  if (isTRUE(avoid_duplicate_regressors)) {
    dedup_est <- .bmk_resolve_duplicate_trials(
      events = events,
      n_tp = n_tp,
      tr = tr,
      basis = basis,
      precision = 0.3,
      return_model = need_event_model,
      tol = duplicate_tol,
      max_iter = duplicate_max_iter,
      jitter_step = duplicate_jitter_step
    )
    events <- dedup_est$events
    design_est <- dedup_est$design
    if (isTRUE(trace) && dedup_est$n_fix > 0L) {
      cat(sprintf("resolved duplicate trial regressors in estimate design (%d adjustment rounds)\n", dedup_est$n_fix))
    }
    if (isTRUE(trace) && isTRUE(dedup_est$unresolved)) {
      cat("warning: unresolved duplicate trial regressors remain in estimate design\n")
    }
  } else {
    design_est <- .bmk_trial_design(
      events = events,
      n_tp = n_tp,
      tr = tr,
      basis = basis,
      return_model = need_event_model
    )
  }

  events_true <- .bmk_apply_true_onset_shift(
    events,
    n_tp = n_tp,
    tr = tr,
    shift_sd = true_onset_shift_sd
  )
  design_true <- if (identical(as.character(true_basis), as.character(basis)) && !isTRUE(true_onset_shift_sd > 0)) {
    design_est
  } else if (isTRUE(avoid_duplicate_regressors)) {
    dedup_true <- .bmk_resolve_duplicate_trials(
      events = events_true,
      n_tp = n_tp,
      tr = tr,
      basis = true_basis,
      precision = 0.3,
      return_model = FALSE,
      tol = duplicate_tol,
      max_iter = duplicate_max_iter,
      jitter_step = duplicate_jitter_step
    )
    if (isTRUE(trace) && dedup_true$n_fix > 0L) {
      cat(sprintf("resolved duplicate trial regressors in true design (%d adjustment rounds)\n", dedup_true$n_fix))
    }
    if (isTRUE(trace) && isTRUE(dedup_true$unresolved)) {
      cat("warning: unresolved duplicate trial regressors remain in true design\n")
    }
    dedup_true$design
  } else {
    .bmk_trial_design(events = events_true, n_tp = n_tp, tr = tr, basis = true_basis)
  }

  if (isTRUE(abs(as.numeric(true_hrf_lag_sec)) > 0) || isTRUE(as.numeric(true_hrf_widen_sec) > 0)) {
    design_true$X <- .bmk_apply_true_hrf_mismatch(
      X = design_true$X,
      tr = tr,
      lag_sec = true_hrf_lag_sec,
      widen_sec = true_hrf_widen_sec
    )
  }

  X_trial <- .bmk_scale_design(design_est$X)
  X_true <- .bmk_scale_design(design_true$X)

  if (ncol(X_trial) != ncol(X_true)) {
    stop(
      sprintf("Estimated and true design matrices differ in ncol (%d vs %d).", ncol(X_trial), ncol(X_true)),
      call. = FALSE
    )
  }

  run_id <- as.integer(design_est$blockids)
  n_trial_eff <- ncol(X_trial)
  n_time <- nrow(X_trial)
  Z_lss <- switch(
    lss_intercept,
    none = NULL,
    global = matrix(1, nrow = n_time, ncol = 1L),
    runwise = {
      run_u <- sort(unique(run_id))
      if (length(run_u) < 2L) {
        matrix(1, nrow = n_time, ncol = 1L)
      } else {
        mm <- stats::model.matrix(~ factor(run_id) - 1L)
        colnames(mm) <- paste0("run_", run_u)
        mm
      }
    }
  )
  sbhm_baseline_model <- NULL
  if (isTRUE(include_sbhm) && lss_intercept != "none") {
    bint <- if (lss_intercept == "runwise" && length(unique(run_id)) > 1L) "runwise" else "global"
    sbhm_baseline_model <- fmridesign::baseline_model(
      basis = "constant",
      sframe = design_est$sframe,
      intercept = bint
    )
  }
  sbhm_obj <- NULL
  if (isTRUE(include_sbhm)) {
    sbhm_rank <- max(1L, as.integer(sbhm_rank))
    if (is.null(sbhm_library_H)) {
      tgrid <- fmrihrf::samples(design_est$sframe, global = TRUE)
      sbhm_library_H <- .bmk_sbhm_library(tgrid = tgrid)
    }
    sbhm_obj <- sbhm_build_fn(
      library_H = as.matrix(sbhm_library_H),
      r = sbhm_rank,
      sframe = design_est$sframe,
      normalize = TRUE,
      shifts = sbhm_shifts
    )
  }

  if (is.null(st_configs)) {
    st_configs <- .bmk_default_st_configs(
      alpha = alpha,
      lambda = lambda,
      overlap_adaptive = overlap_adaptive
    )
  }
  if (is.null(names(st_configs)) || any(names(st_configs) == "")) {
    names(st_configs) <- paste0("st_", seq_along(st_configs))
  }
  if (is.null(primary_st)) {
    if ("st_pool_joint_mult_1se" %in% names(st_configs)) {
      primary_st <- "st_pool_joint_mult_1se"
    } else if ("st_pool_joint_mult" %in% names(st_configs)) {
      primary_st <- "st_pool_joint_mult"
    } else {
      primary_st <- names(st_configs)[1L]
    }
  }

  if (isTRUE(include_oasis)) {
    if (is.null(oasis_grid)) {
      oasis_grid <- expand.grid(
        K = c(1L, 2L),
        ridge_mode = c("absolute", "fractional"),
        ridge_x = c(0, 0.005, 0.02),
        ridge_b = c(0, 0.005, 0.02),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    }
  } else {
    oasis_grid <- data.frame(
      K = integer(0),
      ridge_mode = character(0),
      ridge_x = numeric(0),
      ridge_b = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  n_methods <- length(st_configs) +
    1L +
    (if (isTRUE(include_hrfals)) length(hrfals_modes) else 0L) +
    (if (isTRUE(include_hrfals)) length(hrfals_fit_methods) else 0L) +
    (if (isTRUE(include_oasis)) nrow(oasis_grid) else 0L) +
    (if (isTRUE(include_sbhm)) 1L else 0L)
  rows <- vector("list", as.integer(n_reps) * n_methods)
  idx <- 0L
  sbhm_diag_rows <- vector("list", if (isTRUE(include_sbhm)) as.integer(n_reps) else 0L)
  sbhm_diag_idx <- 0L
  true_hrf_diag_rows <- vector("list", as.integer(n_reps))

  for (r in seq_len(as.integer(n_reps))) {
    if (isTRUE(trace)) cat(sprintf("benchmark rep %d/%d\n", r, n_reps))
    set.seed(as.integer(seed) + r - 1L)
    # TODO: implement per-rep event resampling when resample_events = TRUE

    n_cond_eff <- as.integer(n_conditions)
    cond_labels <- NULL
    if (is.finite(n_cond_eff) && n_cond_eff >= 2L) {
      cond_out <- .bmk_sample_beta_conditions(
        n_trial = n_trial_eff,
        n_vox = n_vox,
        n_conditions = n_cond_eff,
        condition_effect = condition_effect,
        beta_dist = beta_dist,
        beta_scale = beta_scale,
        beta_df = beta_df,
        beta_lognorm_sdlog = beta_lognorm_sdlog,
        beta_sparsity = beta_sparsity
      )
      beta_true <- cond_out$beta
      cond_labels <- cond_out$condition_labels
    } else {
      beta_true <- .bmk_sample_beta(
        n_trial = n_trial_eff,
        n_vox = n_vox,
        beta_dist = beta_dist,
        beta_scale = beta_scale,
        beta_df = beta_df,
        beta_lognorm_sdlog = beta_lognorm_sdlog,
        beta_sparsity = beta_sparsity
      )
    }

    amps <- .bmk_trial_amplitudes(
      n_trial = n_trial_eff,
      model = trial_amplitude_model,
      cv = trial_amplitude_cv,
      outlier_prob = trial_amplitude_outlier_prob,
      outlier_scale = trial_amplitude_outlier_scale
    )
    beta_true <- beta_true * matrix(amps, nrow = n_trial_eff, ncol = n_vox)

    E <- .bmk_sim_noise(
      n_time = n_time,
      n_vox = n_vox,
      run_id = run_id,
      noise_sd = noise_sd,
      noise_model = noise_model,
      ar_rho = ar_rho,
      ar_rho_sd = ar_rho_sd,
      ar_rho_vox_sd = ar_rho_vox_sd,
      ar_phi2 = ar_phi2,
      arma_theta = arma_theta,
      noise_spike_prob = noise_spike_prob,
      noise_spike_scale = noise_spike_scale,
      noise_spike_shared = noise_spike_shared,
      noise_innov_dist = noise_innov_dist,
      noise_df = noise_df,
      spatial_rho = spatial_rho
    )

    true_signal <- .bmk_make_true_signal(
      X_true = X_true,
      beta_true = beta_true,
      tr = tr,
      true_hrf_lag_vox_sd = true_hrf_lag_vox_sd,
      true_hrf_widen_vox_sd = true_hrf_widen_vox_sd,
      true_hrf_lag_trial_sd = true_hrf_lag_trial_sd,
      true_hrf_widen_trial_sd = true_hrf_widen_trial_sd
    )
    empirical_snr <- stats::sd(as.numeric(true_signal)) /
      max(stats::sd(as.numeric(E)), .Machine$double.eps)
    hsum <- attr(true_signal, "hrf_heterogeneity_summary")
    if (is.null(hsum)) hsum <- c(mean_abs_lag = NA_real_, max_abs_lag = NA_real_, mean_widen = NA_real_, max_widen = NA_real_)
    true_hrf_diag_rows[[r]] <- data.frame(
      rep = r,
      mean_abs_lag = as.numeric(hsum[["mean_abs_lag"]]),
      max_abs_lag = as.numeric(hsum[["max_abs_lag"]]),
      mean_widen = as.numeric(hsum[["mean_widen"]]),
      max_widen = as.numeric(hsum[["max_widen"]]),
      empirical_snr = empirical_snr,
      stringsAsFactors = FALSE
    )

    Y <- true_signal + E

    for (m in seq_along(st_configs)) {
      method_name <- names(st_configs)[m]
      cfg <- st_configs[[m]]
      cv_args <- .bmk_or(cfg$cv_args, list())
      cfg$cv_args <- NULL

      st_defaults <- list(
        mode = "cv",
        run_id = run_id,
        alpha = alpha,
        lambda = lambda,
        standardize = FALSE,
        intercept = FALSE,
        cv_folds = as.integer(cv_folds),
        cv_type.measure = cv_type.measure,
        cv_fold_scheme = if (n_runs > 1L) "run" else "random"
      )
      st_opts <- utils::modifyList(st_defaults, cfg, keep.null = TRUE)
      if (length(cv_args)) {
        cv_map <- list(
          cv_foldid = cv_args$foldid,
          cv_folds = cv_args$nfolds,
          cv_type.measure = cv_args$type.measure,
          cv_fold_scheme = cv_args$fold_scheme,
          overlap_low_threshold = cv_args$overlap_low_threshold,
          composite_weights = cv_args$composite_weights
        )
        cv_map <- cv_map[!vapply(cv_map, is.null, logical(1))]
        if (length(cv_map)) {
          st_opts <- utils::modifyList(st_opts, cv_map, keep.null = TRUE)
        }
      }

      timing_st <- system.time({
        beta_st <- tryCatch(
          lss_fn(
            Y = Y,
            X = X_trial,
            Z = if (isTRUE(st_match_lss_nuisance)) Z_lss else NULL,
            method = "stglmnet",
            stglmnet = st_opts
          ),
          error = function(e) e
        )
        score <- list(
          cor = NA_real_,
          rmse = NA_real_,
          slope_bias = NA_real_,
          intercept_bias = NA_real_,
          sign_acc = NA_real_,
          classif_acc = NA_real_
        )
        ok <- FALSE
        msg <- NA_character_

        if (!inherits(beta_st, "error")) {
          if (is.list(beta_st) && !is.null(beta_st$beta)) beta_st <- beta_st$beta
          score <- .bmk_score_betas(as.matrix(beta_st), beta_true)
          if (!is.null(cond_labels)) score$classif_acc <- .bmk_classif_accuracy(beta_st, cond_labels)
          ok <- TRUE
        } else {
          msg <- conditionMessage(beta_st)
        }
      })

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        rep = r,
        method = method_name,
        cor = score$cor,
        rmse = score$rmse,
        slope_bias = score$slope_bias,
        intercept_bias = score$intercept_bias,
        sign_acc = score$sign_acc,
        classif_acc = .bmk_or(score$classif_acc, NA_real_),
        elapsed_sec = as.numeric(timing_st[["elapsed"]]),
        ok = ok,
        message = msg,
        stringsAsFactors = FALSE
      )
    }

    lss_args <- list(Y = Y, X = X_trial, Z = Z_lss, method = lss_method)
    if (!is.null(lss_prewhiten)) lss_args$prewhiten <- lss_prewhiten

    timing_lss <- system.time({
      beta_lss <- tryCatch(do.call(lss_fn, lss_args), error = function(e) e)
      score_lss <- list(
        cor = NA_real_,
        rmse = NA_real_,
        slope_bias = NA_real_,
        intercept_bias = NA_real_,
        sign_acc = NA_real_,
        classif_acc = NA_real_
      )
      ok_lss <- FALSE
      msg_lss <- NA_character_

      if (!inherits(beta_lss, "error")) {
        if (is.list(beta_lss) && !is.null(beta_lss$beta)) beta_lss <- beta_lss$beta
        score_lss <- .bmk_score_betas(beta_lss, beta_true)
        if (!is.null(cond_labels)) score_lss$classif_acc <- .bmk_classif_accuracy(beta_lss, cond_labels)
        ok_lss <- TRUE
      } else {
        msg_lss <- conditionMessage(beta_lss)
      }
    })

    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      rep = r,
      method = "fmrilss_lss",
      cor = score_lss$cor,
      rmse = score_lss$rmse,
      slope_bias = score_lss$slope_bias,
      intercept_bias = score_lss$intercept_bias,
      sign_acc = score_lss$sign_acc,
      classif_acc = .bmk_or(score_lss$classif_acc, NA_real_),
      elapsed_sec = as.numeric(timing_lss[["elapsed"]]),
      ok = ok_lss,
      message = msg_lss,
      stringsAsFactors = FALSE
    )

    if (isTRUE(include_hrfals)) {
      A_hrfals <- if (is.null(Z_lss)) matrix(0, nrow = n_time, ncol = 0L) else as.matrix(Z_lss)
      p_hrfals <- numeric(n_time)
      for (hm in hrfals_modes) {
        use_cpp <- identical(hm, "cpp")
        mname <- sprintf("hrfals_fastlss[%s]", hm)
        timing_hrfals <- system.time({
          beta_hrfals <- tryCatch(
            hrfals_lss_mode_a_fn(
              Y = Y,
              A = A_hrfals,
              C = X_trial,
              p_vec = p_hrfals,
              lambda_ridge = hrfals_lambda_ridge,
              use_cpp = use_cpp
            ),
            error = function(e) e
          )
          score_hrfals <- list(
            cor = NA_real_,
            rmse = NA_real_,
            slope_bias = NA_real_,
            intercept_bias = NA_real_,
            sign_acc = NA_real_,
            classif_acc = NA_real_
          )
          ok_hrfals <- FALSE
          msg_hrfals <- NA_character_
          if (!inherits(beta_hrfals, "error")) {
            score_hrfals <- .bmk_score_betas(beta_hrfals, beta_true)
            if (!is.null(cond_labels)) score_hrfals$classif_acc <- .bmk_classif_accuracy(beta_hrfals, cond_labels)
            ok_hrfals <- TRUE
          } else {
            msg_hrfals <- conditionMessage(beta_hrfals)
          }
        })

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          rep = r,
          method = mname,
          cor = score_hrfals$cor,
          rmse = score_hrfals$rmse,
          slope_bias = score_hrfals$slope_bias,
          intercept_bias = score_hrfals$intercept_bias,
          sign_acc = score_hrfals$sign_acc,
          classif_acc = .bmk_or(score_hrfals$classif_acc, NA_real_),
          elapsed_sec = as.numeric(timing_hrfals[["elapsed"]]),
          ok = ok_hrfals,
          message = msg_hrfals,
          stringsAsFactors = FALSE
        )
      }
    }

    if (isTRUE(include_hrfals) && length(hrfals_fit_methods) > 0L) {
      hrfals_design_obj <- tryCatch(
        hrfals_create_cfals_design_fn(
          fmri_data_obj = Y,
          event_model = design_est$event_model,
          hrf_basis = hrfals_basis,
          confound_obj = if (is.null(Z_lss)) NULL else as.matrix(Z_lss),
          design_control = list(standardize_predictors = FALSE, cache_design_blocks = TRUE)
        ),
        error = function(e) e
      )

      for (hm in hrfals_fit_methods) {
        mname <- sprintf("hrfals_fit[%s]", hm)
        timing_hfit <- system.time({
          score_hfit <- list(
            cor = NA_real_,
            rmse = NA_real_,
            slope_bias = NA_real_,
            intercept_bias = NA_real_,
            sign_acc = NA_real_,
            classif_acc = NA_real_
          )
          ok_hfit <- FALSE
          msg_hfit <- NA_character_

          if (!inherits(hrfals_design_obj, "error")) {
            fit_h <- tryCatch(
              hrfals_from_design_fn(
                y = Y,
                design = hrfals_design_obj,
                method = hm,
                control = hrfals_control
              ),
              error = function(e) e
            )
            if (!inherits(fit_h, "error")) {
              beta_hfit <- if (is.list(fit_h) && !is.null(fit_h$beta_amps)) {
                as.matrix(fit_h$beta_amps)
              } else if (is.list(fit_h) && !is.null(fit_h$beta)) {
                as.matrix(fit_h$beta)
              } else {
                NULL
              }
              if (is.null(beta_hfit) || !is.matrix(beta_hfit)) {
                msg_hfit <- "hrfals fit did not return beta_amps matrix."
              } else {
                # Align to trial design ordering when possible.
                if (!is.null(rownames(beta_hfit)) && !is.null(colnames(X_trial))) {
                  idx_align <- match(colnames(X_trial), rownames(beta_hfit))
                  if (all(!is.na(idx_align))) {
                    beta_hfit <- beta_hfit[idx_align, , drop = FALSE]
                  }
                }
                if (nrow(beta_hfit) != n_trial_eff) {
                  msg_hfit <- sprintf(
                    "hrfals beta row mismatch: expected %d rows, got %d.",
                    n_trial_eff, nrow(beta_hfit)
                  )
                } else {
                  score_hfit <- .bmk_score_betas(beta_hfit, beta_true)
                  if (!is.null(cond_labels)) score_hfit$classif_acc <- .bmk_classif_accuracy(beta_hfit, cond_labels)
                  ok_hfit <- TRUE
                }
              }
            } else {
              msg_hfit <- conditionMessage(fit_h)
            }
          } else {
            msg_hfit <- conditionMessage(hrfals_design_obj)
          }
        })

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          rep = r,
          method = mname,
          cor = score_hfit$cor,
          rmse = score_hfit$rmse,
          slope_bias = score_hfit$slope_bias,
          intercept_bias = score_hfit$intercept_bias,
          sign_acc = score_hfit$sign_acc,
          classif_acc = .bmk_or(score_hfit$classif_acc, NA_real_),
          elapsed_sec = as.numeric(timing_hfit[["elapsed"]]),
          ok = ok_hfit,
          message = msg_hfit,
          stringsAsFactors = FALSE
        )
      }
    }

    if (isTRUE(include_oasis)) {
      for (g in seq_len(nrow(oasis_grid))) {
        grid_row <- oasis_grid[g, , drop = FALSE]
        mname <- sprintf(
          "fmrilss_oasis[K=%d|%s|x=%g|b=%g]",
          as.integer(.bmk_or(grid_row$K, 1L)),
          as.character(grid_row$ridge_mode),
          as.numeric(grid_row$ridge_x),
          as.numeric(grid_row$ridge_b)
        )

        lss_oasis_args <- list(
          Y = Y,
          X = X_trial,
          Z = Z_lss,
          method = "oasis",
          oasis = list(
            K = as.integer(.bmk_or(grid_row$K, 1L)),
            infer_K_from_X = FALSE,
            ridge_mode = as.character(grid_row$ridge_mode),
            ridge_x = as.numeric(grid_row$ridge_x),
            ridge_b = as.numeric(grid_row$ridge_b)
          )
        )
        if (!is.null(oasis_prewhiten)) lss_oasis_args$prewhiten <- oasis_prewhiten

        timing_oasis <- system.time({
          beta_oasis <- tryCatch(do.call(lss_fn, lss_oasis_args), error = function(e) e)
          score_oasis <- list(
            cor = NA_real_,
            rmse = NA_real_,
            slope_bias = NA_real_,
            intercept_bias = NA_real_,
            sign_acc = NA_real_,
            classif_acc = NA_real_
          )
          ok_oasis <- FALSE
          msg_oasis <- NA_character_

          if (!inherits(beta_oasis, "error")) {
            if (is.list(beta_oasis) && !is.null(beta_oasis$beta)) beta_oasis <- beta_oasis$beta
            score_oasis <- .bmk_score_betas(beta_oasis, beta_true)
            if (!is.null(cond_labels)) score_oasis$classif_acc <- .bmk_classif_accuracy(beta_oasis, cond_labels)
            ok_oasis <- TRUE
          } else {
            msg_oasis <- conditionMessage(beta_oasis)
          }
        })

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          rep = r,
          method = mname,
          cor = score_oasis$cor,
          rmse = score_oasis$rmse,
          slope_bias = score_oasis$slope_bias,
          intercept_bias = score_oasis$intercept_bias,
          sign_acc = score_oasis$sign_acc,
          classif_acc = .bmk_or(score_oasis$classif_acc, NA_real_),
          elapsed_sec = as.numeric(timing_oasis[["elapsed"]]),
          ok = ok_oasis,
          message = msg_oasis,
          stringsAsFactors = FALSE
        )
      }
    }

    if (isTRUE(include_sbhm)) {
      timing_sbhm <- system.time({
        sbhm_fit <- tryCatch(
          lss_sbhm_design_fn(
            Y = Y,
            sbhm = sbhm_obj,
            event_model = design_est$event_model,
            baseline_model = sbhm_baseline_model,
            prewhiten = sbhm_prewhiten,
            match = sbhm_match,
            oasis = sbhm_oasis,
            amplitude = sbhm_amplitude,
            return = "amplitude",
            validate = FALSE
          ),
          error = function(e) e
        )
      })

      score_sbhm <- list(
        cor = NA_real_,
        rmse = NA_real_,
        slope_bias = NA_real_,
        intercept_bias = NA_real_,
        sign_acc = NA_real_,
        classif_acc = NA_real_
      )
      ok_sbhm <- FALSE
      msg_sbhm <- NA_character_
      sbhm_diag_idx <- sbhm_diag_idx + 1L
      diag_row <- data.frame(
        rep = r,
        ok = FALSE,
        message = NA_character_,
        mean_margin = NA_real_,
        median_margin = NA_real_,
        mean_rho = NA_real_,
        max_rho = NA_real_,
        mean_kappa = NA_real_,
        max_kappa = NA_real_,
        method_global_ls = NA_integer_,
        method_lss1 = NA_integer_,
        method_oasis_voxel = NA_integer_,
        stringsAsFactors = FALSE
      )
      if (!inherits(sbhm_fit, "error")) {
        beta_sbhm <- if (is.list(sbhm_fit) && !is.null(sbhm_fit$amplitude)) {
          sbhm_fit$amplitude
        } else {
          sbhm_fit
        }
        score_sbhm <- .bmk_score_betas(beta_sbhm, beta_true)
        if (!is.null(cond_labels)) score_sbhm$classif_acc <- .bmk_classif_accuracy(beta_sbhm, cond_labels)
        ok_sbhm <- TRUE
        diag_row$ok <- TRUE
        if (is.list(sbhm_fit)) {
          if (!is.null(sbhm_fit$margin)) {
            diag_row$mean_margin <- mean(sbhm_fit$margin, na.rm = TRUE)
            diag_row$median_margin <- stats::median(sbhm_fit$margin, na.rm = TRUE)
          }
          if (!is.null(sbhm_fit$diag) && is.list(sbhm_fit$diag)) {
            if (!is.null(sbhm_fit$diag$rho_max)) {
              diag_row$mean_rho <- mean(sbhm_fit$diag$rho_max, na.rm = TRUE)
              diag_row$max_rho <- max(sbhm_fit$diag$rho_max, na.rm = TRUE)
            }
            if (!is.null(sbhm_fit$diag$kappa)) {
              diag_row$mean_kappa <- mean(sbhm_fit$diag$kappa, na.rm = TRUE)
              diag_row$max_kappa <- max(sbhm_fit$diag$kappa, na.rm = TRUE)
            }
            if (!is.null(sbhm_fit$diag$method_used)) {
              tab <- table(factor(
                as.character(sbhm_fit$diag$method_used),
                levels = c("global_ls", "lss1", "oasis_voxel")
              ))
              diag_row$method_global_ls <- as.integer(tab[["global_ls"]])
              diag_row$method_lss1 <- as.integer(tab[["lss1"]])
              diag_row$method_oasis_voxel <- as.integer(tab[["oasis_voxel"]])
            }
          }
        }
      } else {
        msg_sbhm <- conditionMessage(sbhm_fit)
        diag_row$message <- msg_sbhm
      }
      sbhm_diag_rows[[sbhm_diag_idx]] <- diag_row

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        rep = r,
        method = "fmrilss_sbhm",
        cor = score_sbhm$cor,
        rmse = score_sbhm$rmse,
        slope_bias = score_sbhm$slope_bias,
        intercept_bias = score_sbhm$intercept_bias,
        sign_acc = score_sbhm$sign_acc,
        classif_acc = .bmk_or(score_sbhm$classif_acc, NA_real_),
        elapsed_sec = as.numeric(timing_sbhm[["elapsed"]]),
        ok = ok_sbhm,
        message = msg_sbhm,
        stringsAsFactors = FALSE
      )
    }
  }

  metrics_long <- do.call(rbind, rows[seq_len(idx)])
  method_summary <- .bmk_method_summary(metrics_long)
  sbhm_diag_long <- if (sbhm_diag_idx > 0L) {
    do.call(rbind, sbhm_diag_rows[seq_len(sbhm_diag_idx)])
  } else {
    data.frame()
  }
  true_hrf_diag_long <- do.call(rbind, true_hrf_diag_rows)

  rep_metrics <- matrix(
    NA_real_,
    nrow = as.integer(n_reps),
    ncol = 10L,
    dimnames = list(NULL, c(
      "cor_st", "cor_lss",
      "rmse_st", "rmse_lss",
      "slope_bias_st", "slope_bias_lss",
      "intercept_bias_st", "intercept_bias_lss",
      "sign_acc_st", "sign_acc_lss"
    ))
  )

  st_dat <- metrics_long[
    metrics_long$method == primary_st,
    c("rep", "cor", "rmse", "slope_bias", "intercept_bias", "sign_acc"),
    drop = FALSE
  ]
  lss_dat <- metrics_long[
    metrics_long$method == "fmrilss_lss",
    c("rep", "cor", "rmse", "slope_bias", "intercept_bias", "sign_acc"),
    drop = FALSE
  ]
  if (nrow(st_dat) > 0L && nrow(lss_dat) > 0L) {
    pair <- merge(st_dat, lss_dat, by = "rep", suffixes = c("_st", "_lss"), all = TRUE)
    rep_metrics[pair$rep, ] <- as.matrix(pair[, c(
      "cor_st", "cor_lss",
      "rmse_st", "rmse_lss",
      "slope_bias_st", "slope_bias_lss",
      "intercept_bias_st", "intercept_bias_lss",
      "sign_acc_st", "sign_acc_lss"
    ), drop = FALSE])
  }

  summary <- list(
    mean = colMeans(rep_metrics, na.rm = TRUE),
    median = apply(rep_metrics, 2L, stats::median, na.rm = TRUE),
    sd = apply(rep_metrics, 2L, stats::sd, na.rm = TRUE),
    st_better_corr_rate = mean(rep_metrics[, "cor_st"] > rep_metrics[, "cor_lss"], na.rm = TRUE),
    st_better_rmse_rate = mean(rep_metrics[, "rmse_st"] < rep_metrics[, "rmse_lss"], na.rm = TRUE),
    st_better_sign_acc_rate = mean(rep_metrics[, "sign_acc_st"] > rep_metrics[, "sign_acc_lss"], na.rm = TRUE),
    st_lower_abs_slope_bias_rate = mean(
      abs(rep_metrics[, "slope_bias_st"]) < abs(rep_metrics[, "slope_bias_lss"]),
      na.rm = TRUE
    ),
    st_lower_abs_intercept_bias_rate = mean(
      abs(rep_metrics[, "intercept_bias_st"]) < abs(rep_metrics[, "intercept_bias_lss"]),
      na.rm = TRUE
    ),
    delta_corr_mean = mean(rep_metrics[, "cor_st"] - rep_metrics[, "cor_lss"], na.rm = TRUE),
    delta_rmse_mean = mean(rep_metrics[, "rmse_st"] - rep_metrics[, "rmse_lss"], na.rm = TRUE),
    delta_sign_acc_mean = mean(rep_metrics[, "sign_acc_st"] - rep_metrics[, "sign_acc_lss"], na.rm = TRUE),
    delta_slope_bias_mean = mean(rep_metrics[, "slope_bias_st"] - rep_metrics[, "slope_bias_lss"], na.rm = TRUE),
    delta_abs_slope_bias_mean = mean(
      abs(rep_metrics[, "slope_bias_st"]) - abs(rep_metrics[, "slope_bias_lss"]),
      na.rm = TRUE
    ),
    delta_intercept_bias_mean = mean(
      rep_metrics[, "intercept_bias_st"] - rep_metrics[, "intercept_bias_lss"],
      na.rm = TRUE
    ),
    delta_abs_intercept_bias_mean = mean(
      abs(rep_metrics[, "intercept_bias_st"]) - abs(rep_metrics[, "intercept_bias_lss"]),
      na.rm = TRUE
    )
  )

  out <- list(
    params = list(
      n_reps = n_reps,
      seed = seed,
      n_tp = n_tp,
      n_runs = n_runs,
      n_trials = n_trials,
      n_trial_eff = n_trial_eff,
      n_vox = n_vox,
      isi = isi,
      isi_model = isi_model,
      isi_params = isi_params,
      tr = tr,
      basis = basis,
      true_basis = true_basis,
      onset_jitter_sd = onset_jitter_sd,
      true_onset_shift_sd = true_onset_shift_sd,
      true_hrf_lag_sec = true_hrf_lag_sec,
      true_hrf_widen_sec = true_hrf_widen_sec,
      true_hrf_lag_vox_sd = true_hrf_lag_vox_sd,
      true_hrf_widen_vox_sd = true_hrf_widen_vox_sd,
      true_hrf_lag_trial_sd = true_hrf_lag_trial_sd,
      true_hrf_widen_trial_sd = true_hrf_widen_trial_sd,
      avoid_duplicate_regressors = avoid_duplicate_regressors,
      duplicate_tol = duplicate_tol,
      duplicate_max_iter = duplicate_max_iter,
      duplicate_jitter_step = duplicate_jitter_step,
      noise_sd = noise_sd,
      noise_model = noise_model,
      ar_rho = ar_rho,
      ar_rho_sd = ar_rho_sd,
      ar_rho_vox_sd = ar_rho_vox_sd,
      ar_phi2 = ar_phi2,
      arma_theta = arma_theta,
      noise_spike_prob = noise_spike_prob,
      noise_spike_scale = noise_spike_scale,
      noise_spike_shared = noise_spike_shared,
      noise_innov_dist = noise_innov_dist,
      noise_df = noise_df,
      spatial_rho = spatial_rho,
      beta_sparsity = beta_sparsity,
      beta_dist = beta_dist,
      beta_scale = beta_scale,
      n_conditions = n_conditions,
      condition_effect = condition_effect,
      trial_amplitude_model = trial_amplitude_model,
      trial_amplitude_cv = trial_amplitude_cv,
      trial_amplitude_outlier_prob = trial_amplitude_outlier_prob,
      trial_amplitude_outlier_scale = trial_amplitude_outlier_scale,
      alpha = alpha,
      cv_folds = cv_folds,
      cv_type.measure = cv_type.measure,
      lss_method = lss_method,
      lss_intercept = lss_intercept,
      st_match_lss_nuisance = st_match_lss_nuisance,
      include_sbhm = include_sbhm,
      sbhm_rank = if (isTRUE(include_sbhm)) sbhm_rank else NA_integer_,
      fast_mode = isTRUE(fast_mode),
      resample_events = isTRUE(resample_events)
    ),
    primary_st = primary_st,
    metrics = as.data.frame(rep_metrics),
    summary = summary,
    metrics_long = metrics_long,
    sbhm_diag_long = sbhm_diag_long,
    true_hrf_diag_long = true_hrf_diag_long,
    method_summary = method_summary,
    oasis_grid = oasis_grid
  )

  class(out) <- "lss_benchmark"
  out
}

# Run multiple benchmark scenarios and aggregate outcomes.
benchmark_lss_methods_suite <- function(
  scenarios = default_benchmark_scenarios("wide"),
  base_args = list(),
  base_seed = 20260301L,
  trace = TRUE
) {
  scenario_list <- .bmk_normalize_scenarios(scenarios)
  scenario_keys <- names(scenario_list)
  ns <- length(scenario_list)

  results <- vector("list", ns)
  rows <- vector("list", ns)
  sbhm_diag_rows <- vector("list", ns)
  true_hrf_diag_rows <- vector("list", ns)
  scenario_design_rows <- vector("list", ns)

  for (i in seq_len(ns)) {
    sc <- scenario_list[[i]]

    if (!is.list(sc)) sc <- list()
    key_name <- if (!is.null(scenario_keys) && nzchar(scenario_keys[i])) scenario_keys[i] else NULL
    sc_name <- as.character(.bmk_or(sc$name, .bmk_or(key_name, paste0("scenario_", i))))
    if (!nzchar(sc_name)) sc_name <- paste0("scenario_", i)
    sc$name <- NULL

    args <- utils::modifyList(base_args, sc, keep.null = TRUE)
    if (is.null(args$seed)) args$seed <- as.integer(base_seed + i - 1L)

    if (isTRUE(trace)) cat(sprintf("scenario %d/%d: %s\n", i, ns, sc_name))

    fit <- do.call(benchmark_lss_methods, args)
    fit$scenario <- sc_name
    fit$scenario_args <- sc

    scenario_design_rows[[i]] <- .bmk_row_df(c(
      list(scenario = sc_name),
      .bmk_flatten_named(fit$params)
    ))

    ml <- fit$metrics_long
    ml$scenario <- sc_name
    if (is.data.frame(fit$sbhm_diag_long) && nrow(fit$sbhm_diag_long) > 0L) {
      dl <- fit$sbhm_diag_long
      dl$scenario <- sc_name
      sbhm_diag_rows[[i]] <- dl
    } else {
      sbhm_diag_rows[[i]] <- NULL
    }
    if (is.data.frame(fit$true_hrf_diag_long) && nrow(fit$true_hrf_diag_long) > 0L) {
      th <- fit$true_hrf_diag_long
      th$scenario <- sc_name
      true_hrf_diag_rows[[i]] <- th
    } else {
      true_hrf_diag_rows[[i]] <- NULL
    }

    results[[i]] <- fit
    rows[[i]] <- ml
  }

  metrics_long <- .bmk_bind_rows(rows)
  rownames(metrics_long) <- NULL
  sbhm_diag_long <- .bmk_bind_rows(sbhm_diag_rows)
  if (is.null(sbhm_diag_long)) sbhm_diag_long <- data.frame()
  true_hrf_diag_long <- .bmk_bind_rows(true_hrf_diag_rows)
  if (is.null(true_hrf_diag_long)) true_hrf_diag_long <- data.frame()
  scenario_design <- .bmk_bind_rows(scenario_design_rows)
  if (is.null(scenario_design)) scenario_design <- data.frame()

  scenario_method_summary <- .bmk_suite_method_summary(metrics_long)

  pairwise_summary <- .bmk_bind_rows(lapply(results, function(x) {
    m <- x$metrics
    dcor <- m$cor_st - m$cor_lss
    drmse <- m$rmse_st - m$rmse_lss
    p_cor <- tryCatch(
      stats::wilcox.test(dcor, mu = 0, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    p_rmse <- tryCatch(
      stats::wilcox.test(drmse, mu = 0, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    dsign <- m$sign_acc_st - m$sign_acc_lss
    dslope <- m$slope_bias_st - m$slope_bias_lss
    dabs_slope <- abs(m$slope_bias_st) - abs(m$slope_bias_lss)
    dinter <- m$intercept_bias_st - m$intercept_bias_lss
    dabs_inter <- abs(m$intercept_bias_st) - abs(m$intercept_bias_lss)
    ci_cor <- .bmk_empirical_ci(dcor)
    ci_rmse <- .bmk_empirical_ci(drmse)
    ci_sign <- .bmk_empirical_ci(dsign)
    ci_slope <- .bmk_empirical_ci(dslope)
    ci_abs_slope <- .bmk_empirical_ci(dabs_slope)
    ci_inter <- .bmk_empirical_ci(dinter)
    ci_abs_inter <- .bmk_empirical_ci(dabs_inter)
    data.frame(
      scenario = x$scenario,
      primary_st = x$primary_st,
      mean_cor_st = mean(m$cor_st, na.rm = TRUE),
      mean_cor_lss = mean(m$cor_lss, na.rm = TRUE),
      mean_rmse_st = mean(m$rmse_st, na.rm = TRUE),
      mean_rmse_lss = mean(m$rmse_lss, na.rm = TRUE),
      mean_sign_acc_st = mean(m$sign_acc_st, na.rm = TRUE),
      mean_sign_acc_lss = mean(m$sign_acc_lss, na.rm = TRUE),
      mean_slope_bias_st = mean(m$slope_bias_st, na.rm = TRUE),
      mean_slope_bias_lss = mean(m$slope_bias_lss, na.rm = TRUE),
      mean_intercept_bias_st = mean(m$intercept_bias_st, na.rm = TRUE),
      mean_intercept_bias_lss = mean(m$intercept_bias_lss, na.rm = TRUE),
      delta_cor_mean = mean(dcor, na.rm = TRUE),
      delta_cor_ci_lo = ci_cor[1],
      delta_cor_ci_hi = ci_cor[2],
      p_wilcox_cor = p_cor,
      delta_rmse_mean = mean(drmse, na.rm = TRUE),
      delta_rmse_ci_lo = ci_rmse[1],
      delta_rmse_ci_hi = ci_rmse[2],
      p_wilcox_rmse = p_rmse,
      delta_sign_acc_mean = mean(dsign, na.rm = TRUE),
      delta_sign_acc_ci_lo = ci_sign[1],
      delta_sign_acc_ci_hi = ci_sign[2],
      delta_slope_bias_mean = mean(dslope, na.rm = TRUE),
      delta_slope_bias_ci_lo = ci_slope[1],
      delta_slope_bias_ci_hi = ci_slope[2],
      delta_abs_slope_bias_mean = mean(dabs_slope, na.rm = TRUE),
      delta_abs_slope_bias_ci_lo = ci_abs_slope[1],
      delta_abs_slope_bias_ci_hi = ci_abs_slope[2],
      delta_intercept_bias_mean = mean(dinter, na.rm = TRUE),
      delta_intercept_bias_ci_lo = ci_inter[1],
      delta_intercept_bias_ci_hi = ci_inter[2],
      delta_abs_intercept_bias_mean = mean(dabs_inter, na.rm = TRUE),
      delta_abs_intercept_bias_ci_lo = ci_abs_inter[1],
      delta_abs_intercept_bias_ci_hi = ci_abs_inter[2],
      st_better_corr_rate = mean(m$cor_st > m$cor_lss, na.rm = TRUE),
      st_better_rmse_rate = mean(m$rmse_st < m$rmse_lss, na.rm = TRUE),
      st_better_sign_acc_rate = mean(m$sign_acc_st > m$sign_acc_lss, na.rm = TRUE),
      st_lower_abs_slope_bias_rate = mean(abs(m$slope_bias_st) < abs(m$slope_bias_lss), na.rm = TRUE),
      st_lower_abs_intercept_bias_rate = mean(abs(m$intercept_bias_st) < abs(m$intercept_bias_lss), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  top_by_cor <- .bmk_bind_rows(lapply(split(scenario_method_summary, scenario_method_summary$scenario), function(d) {
    d <- d[order(d$mean_cor, decreasing = TRUE), , drop = FALSE]
    d[1L, c("scenario", "method", "mean_cor", "mean_rmse", "success_rate"), drop = FALSE]
  }))
  rownames(top_by_cor) <- NULL

  metrics_long <- .bmk_attach_scenario_metadata(metrics_long, scenario_design)
  scenario_method_summary <- .bmk_attach_scenario_metadata(scenario_method_summary, scenario_design)
  pairwise_summary <- .bmk_attach_scenario_metadata(pairwise_summary, scenario_design)
  top_by_cor <- .bmk_attach_scenario_metadata(top_by_cor, scenario_design)
  recommendations <- .bmk_recommend_methods(scenario_method_summary)
  recommendations <- .bmk_attach_scenario_metadata(recommendations, scenario_design)
  recommendation_summary <- .bmk_recommendation_summary(recommendations)

  out <- list(
    scenarios = vapply(results, function(x) x$scenario, character(1)),
    scenario_args = lapply(results, function(x) x$scenario_args),
    scenario_design = scenario_design,
    results = results,
    metrics_long = metrics_long,
    sbhm_diag_long = sbhm_diag_long,
    true_hrf_diag_long = true_hrf_diag_long,
    scenario_method_summary = scenario_method_summary,
    pairwise_summary = pairwise_summary,
    top_by_cor = top_by_cor,
    recommendations = recommendations,
    recommendation_summary = recommendation_summary
  )
  class(out) <- "lss_benchmark_suite"
  out
}

# Build a consistent benchmark plotting theme.
bmk_plot_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )
}

# Utility: apply optional row/column facets to a ggplot object.
.bmk_apply_facets <- function(p, facet_rows = NULL, facet_cols = NULL) {
  if (is.null(facet_rows) && is.null(facet_cols)) return(p)
  row_term <- .bmk_or(facet_rows, ".")
  col_term <- .bmk_or(facet_cols, ".")
  p + ggplot2::facet_grid(stats::as.formula(paste(row_term, "~", col_term)))
}

# Plot a benchmark heatmap over parameter space.
bmk_plot_surface_heatmap <- function(
  data,
  x,
  y,
  metric = "mean_cor",
  facet_rows = "method",
  facet_cols = NULL,
  title = NULL,
  subtitle = NULL,
  x_lab = NULL,
  y_lab = NULL,
  fill_lab = NULL,
  low = "#f7fbff",
  high = "#08306b",
  midpoint = NULL,
  base_size = 11
) {
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)
  metric_sym <- rlang::sym(metric)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = !!x_sym,
      y = !!y_sym,
      fill = !!metric_sym
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25)
  p <- .bmk_apply_facets(p, facet_rows = facet_rows, facet_cols = facet_cols)
  p <- if (is.null(midpoint)) {
    p + ggplot2::scale_fill_gradient(low = low, high = high)
  } else {
    p + ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high, midpoint = midpoint)
  }
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = .bmk_or(x_lab, x),
      y = .bmk_or(y_lab, y),
      fill = .bmk_or(fill_lab, metric)
    ) +
    bmk_plot_theme(base_size = base_size)
}

# Plot benchmark curves across a one-dimensional parameter sweep.
bmk_plot_method_curves <- function(
  data,
  x,
  y = "mean_cor",
  color = "method",
  group = "method",
  facet_rows = NULL,
  facet_cols = NULL,
  ymin = NULL,
  ymax = NULL,
  title = NULL,
  subtitle = NULL,
  x_lab = NULL,
  y_lab = NULL,
  base_size = 11
) {
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)
  color_sym <- rlang::sym(color)
  group_sym <- rlang::sym(group)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = !!x_sym,
      y = !!y_sym,
      color = !!color_sym,
      group = !!group_sym
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9, alpha = 0.9) +
    ggplot2::geom_point(size = 2.2)
  if (!is.null(ymin) && !is.null(ymax)) {
    ymin_sym <- rlang::sym(ymin)
    ymax_sym <- rlang::sym(ymax)
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = !!ymin_sym, ymax = !!ymax_sym),
      width = 0.08,
      alpha = 0.45
    )
  }
  p <- .bmk_apply_facets(p, facet_rows = facet_rows, facet_cols = facet_cols)
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = .bmk_or(x_lab, x),
      y = .bmk_or(y_lab, y),
      color = color
    ) +
    bmk_plot_theme(base_size = base_size)
}

# Plot pairwise deltas with empirical confidence intervals.
bmk_plot_pairwise_delta <- function(
  data,
  x = "delta_cor_mean",
  y = "scenario",
  color = "primary_st",
  facet_rows = NULL,
  facet_cols = NULL,
  xmin = "delta_cor_ci_lo",
  xmax = "delta_cor_ci_hi",
  title = NULL,
  subtitle = NULL,
  x_lab = NULL,
  y_lab = NULL,
  base_size = 11
) {
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)
  color_sym <- rlang::sym(color)
  xmin_sym <- rlang::sym(xmin)
  xmax_sym <- rlang::sym(xmax)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = !!x_sym,
      y = !!y_sym,
      color = !!color_sym
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "gray45") +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = !!xmin_sym,
        xend = !!xmax_sym,
        y = !!y_sym,
        yend = !!y_sym
      ),
      alpha = 0.75
    ) +
    ggplot2::geom_point(size = 2.1)
  p <- .bmk_apply_facets(p, facet_rows = facet_rows, facet_cols = facet_cols)
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = .bmk_or(x_lab, x),
      y = .bmk_or(y_lab, y),
      color = color
    ) +
    bmk_plot_theme(base_size = base_size)
}

# Plot recommended methods over parameter space.
bmk_plot_recommendation_tiles <- function(
  data,
  x,
  y,
  fill = "recommended_method",
  label = "recommended_method",
  facet_rows = NULL,
  facet_cols = NULL,
  title = NULL,
  subtitle = NULL,
  x_lab = NULL,
  y_lab = NULL,
  base_size = 11
) {
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)
  fill_sym <- rlang::sym(fill)
  label_sym <- rlang::sym(label)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = !!x_sym,
      y = !!y_sym,
      fill = !!fill_sym,
      label = !!label_sym
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(size = 3.0, lineheight = 0.9)
  p <- .bmk_apply_facets(p, facet_rows = facet_rows, facet_cols = facet_cols)
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = .bmk_or(x_lab, x),
      y = .bmk_or(y_lab, y),
      fill = fill
    ) +
    bmk_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "none")
}

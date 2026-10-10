# GLMdenoise and ridge stages of glmsingle(), processed per HRF group and
# voxel tile.

# Types A and B in one pass over the data.
#
# Type A (ON-OFF): one canonical-HRF regressor for all trials, fitted on the
# stacked runs (one coefficient per voxel shared across runs).
# Type B: single-trial OLS for every library HRF; keep the best HRF per voxel
# (first maximum of R^2, numpy argmax semantics) and its betas.
.glms_fit_types_ab <- function(Ylist, geom, hrf0, library, nuis, tiles, singular) {
  R <- length(Ylist)
  n_vox <- ncol(Ylist[[1]])
  n_hrf <- ncol(library)
  stats_a <- lapply(seq_len(R), function(r) {
    one <- matrix(rowSums(.glms_trial_design(geom$onsets[[r]], geom$n_time[r], hrf0)))
    .glms_design_stats_x(one, nuis[[r]], "k0", solver = FALSE)
  })
  g_a <- sum(vapply(stats_a, function(st) st$G[1, 1], numeric(1)))
  onoff_r2 <- beta_a <- numeric(n_vox)
  stats <- lapply(seq_len(n_hrf), function(h) lapply(seq_len(R), function(r) {
    .glms_design_stats(geom$onsets[[r]], geom$n_time[r], library[, h], nuis[[r]], "k0", singular)
  }))
  fit_r2 <- matrix(NA_real_, n_vox, n_hrf)
  fit_r2_run <- matrix(NA_real_, n_vox, R * n_hrf)  # becomes voxels x runs x HRFs
  beta <- matrix(0, geom$n_trials, n_vox)
  for (vox in tiles) {
    dats <- lapply(seq_len(R), function(r) .glms_data_stats(.glms_cols(Ylist[[r]], vox), nuis[[r]]))
    rhs <- lapply(seq_len(R), function(r) .glms_rhs(stats_a[[r]], dats[[r]]))
    b <- Reduce(`+`, lapply(rhs, function(x) x$b[1, ]))
    bt <- if (g_a > 0) b / g_a else 0 * b
    sse <- Reduce(`+`, lapply(seq_len(R), function(r) {
      .glms_sse(matrix(bt, 1L), rhs[[r]]$c, crossprod(stats_a[[r]]$Ap), dats[[r]]$s)
    }))
    onoff_r2[vox] <- .glms_r2(sse, Reduce(`+`, lapply(dats, `[[`, "s")))
    beta_a[vox] <- bt
    # score every HRF without forming betas ...
    for (h in seq_len(n_hrf)) {
      sc <- .glms_ols_score(stats[[h]], dats)
      fit_r2[vox, h] <- .glms_r2(colSums(sc$sse), colSums(sc$s))
      fit_r2_run[vox, (h - 1L) * R + seq_len(R)] <- t(.glms_r2(sc$sse, sc$s))
    }
    # ... then reconstruct betas only for each voxel's winning HRF
    win <- .glms_argmax_rows(fit_r2[vox, , drop = FALSE])
    for (h in unique(win)) {
      sel <- which(win == h)
      beta[, vox[sel]] <- .glms_fit_ols(stats[[h]], .glms_dats_subset(dats, sel))$beta
    }
  }
  dim(fit_r2_run) <- c(n_vox, R, n_hrf)
  hrf_index <- .glms_argmax_rows(fit_r2)
  hrf_index_run <- matrix(0L, n_vox, R)
  for (r in seq_len(R)) hrf_index_run[, r] <- .glms_argmax_rows(matrix(fit_r2_run[, r, ], n_vox))
  sel <- cbind(seq_len(n_vox), hrf_index)
  r2run <- vapply(seq_len(R), function(r) matrix(fit_r2_run[, r, ], n_vox)[sel], numeric(n_vox))
  list(onoffR2 = onoff_r2, beta_a = beta_a,
       FitHRFR2 = fit_r2, FitHRFR2run = fit_r2_run, HRFindex = hrf_index,
       HRFindexrun = hrf_index_run, R2 = fit_r2[sel],
       R2run = matrix(r2run, n_vox), beta = beta)
}

# GLMdenoise: noise pool, PCs, and the number of PCs by cross-validation.
.glms_denoise <- function(Ylist, geom, library, hrf_index, onoff_r2, meanvol,
                          nuis, poly_degree, extras, max_pcs, pc_stop,
                          noise_pool_brightness, noise_pool_r2, noise_pool_mask, pc_voxel_r2,
                          pc_voxel_mask, pc_cv_all_voxels, nuisance_in_denoise,
                          cv_zero_variance, singular, chunk_size, say) {
  R <- length(Ylist)
  n_vox <- ncol(Ylist[[1]])
  bright <- meanvol > .glms_percentile(meanvol, noise_pool_brightness[1]) * noise_pool_brightness[2]
  noisepool <- bright & !is.na(onoff_r2) & onoff_r2 < noise_pool_r2
  if (!is.null(noise_pool_mask)) noisepool <- noisepool & as.logical(noise_pool_mask)
  pcs <- lapply(seq_len(R), function(r) {
    .glms_noise_pcs(Ylist[[r]], noisepool, nuis[[r]]$Qp, max_pcs, chunk_size)
  })
  out <- list(pcregressors = pcs, noisepool = noisepool, pcnum = 0L,
              xvaltrend = NULL, glmbadness = NULL, pcvoxels = NULL)
  available <- min(vapply(pcs, ncol, integer(1)))
  requested <- if (pc_stop <= 0) -pc_stop else max_pcs
  if (requested > available) {
    warning(sprintf("Limiting the number of noise PCs to %d (noise-pool rank across runs).", available),
            call. = FALSE)
    max_pcs <- available
  }
  if (!available) return(out)
  if (pc_stop <= 0) {
    out$pcnum <- as.integer(min(-pc_stop, max_pcs))
    return(out)
  }
  say("Cross-validating the number of noise regressors")
  mask <- if (is.null(pc_voxel_mask)) rep(TRUE, n_vox) else as.logical(pc_voxel_mask)
  ix <- which(!is.na(onoff_r2) & onoff_r2 > pc_voxel_r2 & mask)
  if (!length(ix)) {
    warning("No voxels passed pc_voxel_r2; using the best 100 voxels.", call. = FALSE)
    ix2 <- which(mask)
    if (!length(ix2)) stop("no voxels are in pc_voxel_mask", call. = FALSE)
    ix <- ix2[order(onoff_r2[ix2], decreasing = TRUE, na.last = TRUE)][seq_len(min(100L, length(ix2)))]
  }
  cv_vox <- if (pc_cv_all_voxels) seq_len(n_vox) else ix
  bad <- matrix(NA_real_, n_vox, max_pcs + 1L)
  bad[cv_vox, ] <- .glms_pc_badness(Ylist, geom, library, hrf_index, cv_vox,
                                    poly_degree, extras, pcs, max_pcs,
                                    nuisance_in_denoise, cv_zero_variance, singular,
                                    chunk_size)
  out$xvaltrend <- -apply(bad[ix, , drop = FALSE], 2L, stats::median)
  out$glmbadness <- bad
  out$pcvoxels <- seq_len(n_vox) %in% ix
  out$pcnum <- .glms_select_pcs(out$xvaltrend, pc_stop)
  out
}

.glms_extras_rule <- function(nuisance_in_denoise) {
  if (identical(nuisance_in_denoise, "always")) function(k) TRUE else function(k) k > 0L
}

# numpy.argmin semantics per row: first minimum; a NaN wins at its position.
.glms_argmin_rows <- function(M) {
  .glms_argmax_rows(-M)
}

.glms_voxel_groups <- function(hrf_index, vox, chunk_size) {
  out <- list()
  for (h in sort(unique(hrf_index[vox]))) {
    vh <- vox[hrf_index[vox] == h]
    for (blk in split(vh, ceiling(seq_along(vh) / chunk_size))) {
      out[[length(out) + 1L]] <- list(h = h, vox = blk)
    }
  }
  out
}

# Cross-validation loss (voxels x (max_pcs + 1)) of OLS fits with 0..max_pcs
# noise PCs, scored against the 0-PC fit, for the voxels in `cv_vox`.
.glms_pc_badness <- function(Ylist, geom, library, hrf_index, cv_vox,
                             poly_degree, extras, pcs, max_pcs,
                             nuisance_in_denoise, cv_zero_variance, singular,
                             chunk_size) {
  R <- length(Ylist)
  ks <- 0:max_pcs
  rule <- .glms_extras_rule(nuisance_in_denoise)
  nuis <- lapply(seq_len(R), function(r) {
    .glms_nuisance_run(geom$n_time[r], poly_degree[r], extras[[r]], pcs[[r]], ks, rule)
  })
  knames <- paste0("k", ks)
  out <- matrix(NA_real_, length(cv_vox), length(ks))
  stats_cache <- list()
  for (grp in .glms_voxel_groups(hrf_index, cv_vox, chunk_size)) {
    key <- as.character(grp$h)
    if (is.null(stats_cache[[key]])) {
      stats_cache[[key]] <- lapply(knames, function(kn) lapply(seq_len(R), function(r) {
        .glms_design_stats(geom$onsets[[r]], geom$n_time[r], library[, grp$h], nuis[[r]], kn, singular)
      }))
    }
    stats_k <- stats_cache[[key]]
    dats <- lapply(seq_len(R), function(r) .glms_data_stats(.glms_cols(Ylist[[r]], grp$vox), nuis[[r]]))
    ref <- .glms_fit_ols(stats_k[[1L]], dats)$beta
    cv <- .glms_cv_compile(geom, ref, cv_zero_variance)
    rows <- match(grp$vox, cv_vox)
    out[rows, 1L] <- .glms_cv_loss_ref(cv, ref[cv$used, , drop = FALSE])
    for (j in seq_along(ks)[-1L]) {
      beta <- .glms_fit_ols(stats_k[[j]], dats)$beta
      out[rows, j] <- .glms_cv_loss(cv, beta[cv$used, , drop = FALSE])
    }
  }
  out
}

# Cross-validation losses (fractions x voxels) of ridge candidates, fused in
# C++: per fraction and run, shrink the eigen coefficients, map them to the
# cross-validated trial rows and accumulate the compiled loss.
.glms_frac_cv_losses <- function(sp, cv, alphas) {
  R <- length(sp$a)
  Vu <- s2 <- a <- vector("list", R)
  offsets <- integer(R)
  off <- 0L; used_off <- 0L
  for (r in seq_len(R)) {
    st <- sp$stats[[r]]
    n <- if (is.null(st$V)) 0L else nrow(st$V)
    keep <- cv$used[cv$used > off & cv$used <= off + n] - off
    offsets[r] <- used_off
    if (length(keep)) {
      Vu[[r]] <- st$V[keep, , drop = FALSE]
      s2[[r]] <- st$s2
      a[[r]] <- sp$a[[r]]
    } else {
      Vu[[r]] <- matrix(0, 0L, 0L); s2[[r]] <- numeric(0); a[[r]] <- matrix(0, 0L, 0L)
    }
    used_off <- used_off + length(keep)
    off <- off + n
  }
  glms_frac_cv_loss(Vu, s2, a, offsets, alphas, cv$mu, cv$isd, cv$session_used,
                    cv$d, cv$M, cv$const, .glms_nt())
}

# Types C (GLMdenoise, unregularised) and D (fractional ridge) for all voxels.
.glms_fit_cd <- function(Ylist, geom, library, hrf_index, poly_degree, extras,
                         pcs, pcnum, nuisance_in_denoise, fracstouse, fracs,
                         ridge, ridge_rescale, cv_zero_variance,
                         ridge_alpha, chunk_size) {
  R <- length(Ylist)
  n_vox <- length(hrf_index)
  N <- geom$n_trials
  rule <- .glms_extras_rule(nuisance_in_denoise)
  nuis <- lapply(seq_len(R), function(r) {
    .glms_nuisance_run(geom$n_time[r], poly_degree[r], extras[[r]],
                       if (pcnum > 0L) pcs[[r]] else NULL, pcnum, rule)
  })
  kname <- paste0("k", pcnum)
  beta_c <- beta_d <- matrix(0, N, n_vox)
  r2_c <- r2_d <- rep(NA_real_, n_vox)
  r2run_c <- r2run_d <- matrix(NA_real_, n_vox, R)
  frac_value <- rep(NA_real_, n_vox)
  cv_select <- ridge && length(fracs) > 1L
  rrbadness <- if (cv_select) matrix(NA_real_, n_vox, length(fracs)) else NULL
  autoscale <- ridge && ridge_rescale && !(length(fracs) == 1L && fracs == 1)
  scaleoffset <- if (autoscale) matrix(NA_real_, n_vox, 2L,
                                       dimnames = list(NULL, c("scale", "offset"))) else NULL
  prepended <- fracs[1L] != 1
  stats_cache <- list()
  for (grp in .glms_voxel_groups(hrf_index, seq_len(n_vox), chunk_size)) {
    key <- as.character(grp$h)
    if (is.null(stats_cache[[key]])) {
      stats_cache[[key]] <- lapply(seq_len(R), function(r) {
        .glms_design_stats(geom$onsets[[r]], geom$n_time[r], library[, grp$h], nuis[[r]], kname,
                           solver = FALSE, spectral = TRUE)
      })
    }
    stats <- stats_cache[[key]]
    v <- grp$vox
    dats <- lapply(seq_len(R), function(r) .glms_data_stats(.glms_cols(Ylist[[r]], v), nuis[[r]], kname))
    sp <- .glms_spectral(stats, dats)
    w1 <- .glms_shrink(sp, rep(0, length(v)))
    b1 <- .glms_ridge_coef(sp, w1)
    q <- .glms_ridge_sse(sp, w1)
    beta_c[, v] <- b1
    r2_c[v] <- .glms_r2(colSums(q$sse), colSums(q$s))
    r2run_c[v, ] <- t(.glms_r2(q$sse, q$s))
    if (!ridge) next

    alphas <- .glms_frac_alphas(sp, fracstouse, ridge_alpha)
    if (cv_select) {
      cv <- .glms_cv_compile(geom, b1, cv_zero_variance)
      bad <- matrix(NA_real_, length(v), length(fracstouse))
      bad[, 1L] <- .glms_cv_loss_ref(cv, b1[cv$used, , drop = FALSE])
      if (length(fracstouse) > 1L) {
        bad[, -1L] <- t(.glms_frac_cv_losses(sp, cv, alphas[-1L, , drop = FALSE]))
      }
      if (prepended) {
        idx <- .glms_argmin_rows(bad[, -1L, drop = FALSE]) + 1L
        rrbadness[v, ] <- bad[, -1L]
      } else {
        idx <- .glms_argmin_rows(bad)
        rrbadness[v, ] <- bad
      }
    } else {
      idx <- rep(length(fracstouse), length(v))
    }
    chosen <- alphas[cbind(idx, seq_along(v))]
    wd <- .glms_shrink(sp, chosen)
    bd <- .glms_ridge_coef(sp, wd)
    q <- .glms_ridge_sse(sp, wd)
    r2_d[v] <- .glms_r2(colSums(q$sse), colSums(q$s))
    r2run_d[v, ] <- t(.glms_r2(q$sse, q$s))
    frac_value[v] <- fracstouse[idx]
    if (autoscale) {
      sc <- .glms_autoscale(bd, b1)
      bd <- sc$fitted
      scaleoffset[v, ] <- sc$h
    }
    beta_d[, v] <- bd
  }
  list(beta_c = beta_c, r2_c = r2_c, r2run_c = r2run_c,
       beta_d = beta_d, r2_d = r2_d, r2run_d = r2run_d,
       frac_value = frac_value, rrbadness = rrbadness,
       scaleoffset = scaleoffset)
}

#' Prewhitened LSS through the weight-matrix kernel
#'
#' Fits the noise model once, then, for every distinct whitening operator
#' (one for global/run pooling, one per parcel or AR bin otherwise), whitens
#' the trial design and confounds with that operator, builds the LSS weight
#' matrix of the whitened design and applies it to the matching whitened
#' voxels. Because each operator gets its own filtered design, voxel- and
#' parcel-specific noise models are supported for a shared design, which the
#' generic prewhitening path cannot do.
#'
#' `pooling = "voxel"` follows the strategy of Nilearn's AR(1) GLM (Worsley et
#' al., 2002): per-voxel autocorrelations are estimated from OLS residuals,
#' voxels are grouped into `voxel_bins` bins of similar autocorrelation, and an
#' AR model is refitted per bin from the pooled residuals of its voxels, which
#' is both cheaper and less noisy than one model per voxel.
#'
#' @param Y Data (n x V).
#' @param X Trial design (n x T).
#' @param Z Fixed regressors (n x F), always present in lss().
#' @param Nuisance Optional nuisance regressors (n x N).
#' @param opts Resolved prewhitening options.
#' @param groups NULL or integer trial group codes.
#' @param method LSS method used to build weights ("r_optimized",
#'   "cpp_optimized" or "cpp").
#' @param ridge Length-2 fractional ridge.
#' @return List with `beta` (T x V) and `whiten_plan`.
#' @keywords internal
#' @noRd
.lss_prewhitened <- function(Y, X, Z, Nuisance, opts, groups = NULL,
                             method = "r_optimized", ridge = c(0, 0)) {
  n <- nrow(Y)
  V <- ncol(Y)
  K <- cbind(Z, Nuisance)

  if (opts$pooling == "parcel" && length(opts$parcels) != V) {
    stop("prewhiten$parcels must have length ncol(Y)", call. = FALSE)
  }
  if (opts$pooling == "run" && length(opts$runs) != n) {
    stop("prewhiten$runs must have length nrow(Y)", call. = FALSE)
  }

  plan <- opts$.whiten_plan
  if (is.null(plan)) {
    X_model <- if (opts$residual_model %in% c("full", "corrected")) {
      X
    } else {
      .aggregate_trials(X, groups)
    }
    noise_design <- cbind(Z, X_model, Nuisance)
    plan <- if (opts$pooling == "voxel") {
      resid <- if (opts$compute_residuals) {
        .noise_residuals(Y, noise_design)
      } else {
        Y
      }
      .fit_binned_voxel_plan(resid, opts)
    } else {
      .fit_noise_plan(Y, noise_design, opts)
    }
  }

  parcels <- if (identical(plan$pooling, "parcel")) {
    as.integer(plan$parcels %||% opts$parcels)
  } else {
    NULL
  }

  # Whiten Y and the joint design [X, K] in one call; parcel plans return one
  # filtered design per parcel.
  design <- cbind(X, K)
  wh <- fmriAR::whiten_apply(plan, X = design, Y = Y, parcels = parcels)
  T_trials <- ncol(X)
  build_w <- function(D) {
    C_res <- .lss_residualize_trials(D[, seq_len(T_trials), drop = FALSE],
                                     D[, -seq_len(T_trials), drop = FALSE])
    if (method == "r_optimized") {
      .lss_weight_matrix(C_res, groups, ridge = ridge)
    } else {
      lss_weight_matrix_cpp(C_res, if (is.null(groups)) integer(0) else groups,
                            1e-12, ridge[1L], ridge[2L])
    }
  }

  if (is.null(parcels)) {
    beta <- crossprod(build_w(wh$X), wh$Y)
  } else {
    beta <- matrix(0, T_trials, V)
    for (key in names(wh$X_by)) {
      cols <- which(parcels == as.integer(key))
      if (!length(cols)) next
      beta[, cols] <- crossprod(build_w(wh$X_by[[key]]), wh$Y[, cols, drop = FALSE])
    }
  }
  list(beta = beta, whiten_plan = plan)
}

#' Fit the noise plan for global, run or parcel pooling
#'
#' For AR models with global or run pooling, fmriAR's estimate depends on the
#' residuals only through their pooled lag products, i.e. through the
#' residual Gram matrix \eqn{E E^\top = M Y Y^\top M}. When there are more
#' voxels than timepoints the plan is therefore fitted from an n-column
#' factor with the same Gram matrix (scaled to the same per-column average),
#' computed from one symmetric rank-k update of Y. This avoids materializing
#' the n x V residual matrix and is exact up to floating point.
#'
#' @param Y Data (n x V).
#' @param design_full Design used to compute OLS residuals.
#' @param opts Resolved prewhitening options.
#' @keywords internal
#' @noRd
.fit_noise_plan <- function(Y, design_full, opts) {
  n <- nrow(Y)
  V <- ncol(Y)
  if (identical(opts$residual_model, "corrected")) {
    # fmriAR corrects the residual autocovariance for projection onto this
    # design; the orthonormal basis spans exactly the projected space,
    # including the intercept added by .noise_basis().
    opts$design <- .noise_basis(design_full, n)
  }
  compressible <- identical(opts$method, "ar") &&
    opts$pooling %in% c("global", "run") && V > 2L * n
  if (compressible) {
    factor <- .pooled_residual_factor(Y, design_full,
                                      residualize = opts$compute_residuals)
    if (!is.null(factor)) return(.fit_whiten_plan(factor, opts))
  }
  resid <- if (opts$compute_residuals) .noise_residuals(Y, design_full) else Y
  .fit_whiten_plan(resid, opts)
}

#' Low-rank factor reproducing the pooled residual lag products
#'
#' Returns F (n x k) with \eqn{F F^\top / k = E E^\top / V}, where E are the
#' OLS residuals of Y (or Y itself when `residualize = FALSE`). Every pooled
#' (column-averaged) lag product, after any per-run centering, is identical
#' for F and E.
#'
#' @keywords internal
#' @noRd
.pooled_residual_factor <- function(Y, design_full, residualize = TRUE) {
  n <- nrow(Y)
  V <- ncol(Y)
  G <- tcrossprod(Y)
  if (residualize) {
    basis <- .noise_basis(design_full, n)
    if (ncol(basis) > 0L) {
      GB <- G %*% basis
      BtGB <- crossprod(basis, GB)
      G <- G - GB %*% t(basis) - basis %*% t(GB) + basis %*% BtGB %*% t(basis)
      G <- (G + t(G)) / 2
    }
  }
  eg <- eigen(G, symmetric = TRUE)
  tol <- n * max(abs(eg$values)) * .Machine$double.eps
  keep <- eg$values > tol
  k <- sum(keep)
  if (k == 0L) return(NULL)
  eg$vectors[, keep, drop = FALSE] *
    rep(sqrt(eg$values[keep] * k / V), each = n)
}

#' Orthonormal basis of the noise-estimation design (with intercept)
#' @keywords internal
#' @noRd
.noise_basis <- function(design_full, n_time) {
  if (is.null(design_full) || ncol(design_full) == 0L) {
    return(matrix(1 / sqrt(n_time), n_time, 1L))
  }
  qr0 <- qr(design_full)
  r1 <- qr.resid(qr0, rep(1, n_time))
  if (sqrt(sum(r1^2)) >= 1e-8) {
    qr0 <- qr(cbind(1, design_full))
  }
  if (qr0$rank == 0L) return(matrix(numeric(), n_time, 0L))
  qr.Q(qr0)[, seq_len(qr0$rank), drop = FALSE]
}

#' OLS residuals used to estimate the noise model
#'
#' Adds an intercept when the constant is not already in the span of the
#' design, then removes the rank-revealed column space with two BLAS-3
#' products (much faster than Householder application via `qr.resid()` for
#' many voxels).
#'
#' @keywords internal
#' @noRd
.noise_residuals <- function(Y, design_full) {
  basis <- .noise_basis(design_full, nrow(Y))
  if (ncol(basis) == 0L) return(Y)
  Y - basis %*% crossprod(basis, Y)
}

#' Fit an fmriAR noise plan from residuals
#' @keywords internal
#' @noRd
.fit_whiten_plan <- function(resid, opts) {
  fmriAR::fit_noise(
    resid = resid,
    runs = opts$runs,
    method = opts$method,
    p = opts$p,
    q = opts$q,
    p_max = opts$p_max,
    exact_first = opts$exact_first,
    pooling = opts$pooling,
    parcels = opts$parcels,
    design = opts$design,
    acvf_correction = opts$acvf_correction,
    correction_max_lag = opts$correction_max_lag
  )
}

#' Voxel-adaptive noise plan from autocorrelation bins
#'
#' Each bin's AR model is estimated from the autocovariance pooled over its
#' member voxels (global pooling restricted to the bin), not from the bin-mean
#' series that fmriAR's parcel pooling uses: bins are groups of voxels with
#' similar autocorrelation, not spatially coherent parcels.
#'
#' @return An fmriAR parcel plan whose parcels are the bins.
#' @keywords internal
#' @noRd
.fit_binned_voxel_plan <- function(resid, opts) {
  bins <- .ar_voxel_bins(resid, opts)
  bin_opts <- opts
  bin_opts$pooling <- "global"
  bin_opts$parcels <- NULL
  ids <- sort(unique(bins))
  plans <- lapply(ids, function(b) {
    .fit_whiten_plan(resid[, bins == b, drop = FALSE], bin_opts)
  })
  key <- as.character(ids)
  phi <- stats::setNames(lapply(plans, function(pl) pl$phi[[1L]]), key)
  theta <- stats::setNames(lapply(plans, function(pl) {
    th <- pl$theta[[1L]]
    if (is.null(th)) numeric(0) else th
  }), key)
  plan <- fmriAR::compat$plan_from_phi(
    phi = phi, theta = theta, runs = opts$runs, parcels = bins,
    pooling = "parcel", exact_first = isTRUE(plans[[1L]]$exact_first),
    method = opts$method
  )
  plan$voxel_bins <- bins
  plan
}

#' Group voxels into bins of similar residual autocorrelation
#'
#' Features are the per-voxel autocorrelations at lags 1..p (p = 2 when the
#' AR order is selected automatically). Bins come from k-means with a
#' deterministic quantile initialization, so results are reproducible
#' without touching the RNG.
#'
#' @return Integer bin labels (length V).
#' @keywords internal
#' @noRd
.ar_voxel_bins <- function(resid, opts) {
  V <- ncol(resid)
  n_bins <- opts$voxel_bins %||% 50L
  lag <- if (is.numeric(opts$p)) as.integer(opts$p) else 2L
  lag <- max(1L, lag + opts$q)
  runs <- opts$runs %||% rep(1L, nrow(resid))
  run_starts <- which(c(TRUE, runs[-1L] != runs[-length(runs)])) - 1L
  feats <- t(voxel_acf_cpp(resid, as.integer(run_starts), lag))

  n_bins <- min(n_bins, V)
  if (n_bins <= 1L) return(rep(1L, V))
  # Quantile initialization along the leading feature direction.
  score <- if (ncol(feats) > 1L) {
    drop(scale(feats, scale = FALSE) %*% stats::prcomp(feats, rank. = 1L)$rotation)
  } else {
    feats[, 1L]
  }
  ord <- order(score)
  init_idx <- ord[pmin(V, floor((seq_len(n_bins) - 0.5) * V / n_bins) + 1L)]
  centers <- unique(feats[init_idx, , drop = FALSE])
  if (nrow(centers) <= 1L) return(rep(1L, V))
  km <- tryCatch(
    suppressWarnings(stats::kmeans(feats, centers = centers, iter.max = 25L)),
    error = function(e) NULL
  )
  if (!is.null(km)) return(as.integer(factor(km$cluster)))
  # Fallback: equal-count bins along the leading direction.
  bins <- integer(V)
  bins[ord] <- ceiling(seq_len(V) * nrow(centers) / V)
  bins
}

#' Low-dimensional trial summary for noise-model residuals
#'
#' Fitting the noise model to residuals of the full trial-wise (LSA) design
#' removes up to one degree of freedom per trial; in rapid designs with many
#' trials this biases the residual autocorrelation strongly downward (often
#' to negative AR coefficients). The aggregate model keeps the confounds and
#' one summed regressor per trial group (or per basis function), which is
#' also what the per-trial LSS models of Nilearn leave in their residuals.
#'
#' @param X Trial design (n x T) or NULL.
#' @param groups NULL or integer codes (length T) for the summed columns.
#' @return n x G matrix (or NULL).
#' @keywords internal
#' @noRd
.aggregate_trials <- function(X, groups = NULL) {
  if (is.null(X) || ncol(X) == 0L) return(X)
  if (is.null(groups)) return(matrix(rowSums(X), ncol = 1L))
  vapply(sort(unique(groups)), function(g) {
    rowSums(X[, groups == g, drop = FALSE])
  }, numeric(nrow(X)))
}

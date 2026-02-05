# File: R/oasis_backend.R
# OASIS backend for fmrilss::lss

#' OASIS backend for fmrilss::lss (internal entry)
#'
#' Add `method="oasis"` to fmrilss::lss(). This path:
#'   - (optionally) builds trial-wise design X via fmrihrf
#'   - residualizes Y (and X downstream) against confounds + Z + other-condition aggregates
#'   - computes all trial betas in one batched pass via the closed-form LSS (exact; or ridge-LSS)
#'   - optionally returns per-trial SEs and design diagnostics
#'
#' @param Y (T x V) numeric matrix
#' @param X (T x N_trials) trial-wise design (if NULL, use oasis$design_spec to build)
#' @param Z (T x K) fixed experimental regressors to be projected out
#' @param Nuisance (T x P) confounds (intercept, motion, drift, aCompCor, ...)
#' @param oasis list of options:
#'    - design_spec: list describing events/HRF to build X via fmrihrf
#'    - K: explicit basis dimension (auto-detected if not provided)
#'    - ridge_mode: "absolute" (default) or "fractional"
#'    - ridge_x, ridge_b: nonnegative ridge on the [a_j, b_j] Gram (default 0 -> exact LSS)
#'    - block_cols: voxel block size (default 4096)
#'    - return_se: logical (default FALSE)
#'    - return_diag: logical (default FALSE)
#'    - whiten: "none" | "ar1" (default "none"); if "ar1", prewhiten Y and design first (DEPRECATED: use prewhiten parameter)
#' @param prewhiten list of prewhitening options using fmriAR (see ?lss for details)
#'
#' @return by default: (N_trials x V) matrix of betas; if `return_se` or `return_diag`, a list
#' @keywords internal
#' @noRd
.lss_oasis <- function(Y, X = NULL, Z = NULL, Nuisance = NULL, oasis = list(), prewhiten = NULL) {
  # Coerce to base matrices early & validate ---------------------------------
  Y        <- .as_base_matrix(Y)
  X        <- .as_base_matrix(X)
  Z        <- .as_base_matrix(Z)
  Nuisance <- .as_base_matrix(Nuisance)

  if (!is.matrix(Y)) stop("Y must be a numeric matrix")
  if (any(!is.finite(Y))) stop("Y contains non-finite values")

  n_time <- nrow(Y)
  V      <- ncol(Y)

  if (!is.null(X)) {
    if (!is.matrix(X)) stop("X must be a matrix")
    if (nrow(X) != n_time) stop("X must have the same number of rows as Y")
    if (any(!is.finite(X))) stop("X contains non-finite values")
  }
  
  if (!is.null(Z)) {
    if (!is.matrix(Z)) stop("Z must be a matrix")
    if (nrow(Z) != n_time) stop("Z must have the same number of rows as Y")
    if (any(!is.finite(Z))) stop("Z contains non-finite values")
  }
  
  if (!is.null(Nuisance)) {
    if (!is.matrix(Nuisance)) stop("Nuisance must be a matrix")
    if (nrow(Nuisance) != n_time) stop("Nuisance must have the same number of rows as Y")
    if (any(!is.finite(Nuisance))) stop("Nuisance contains non-finite values")
  }
  
  # Validate oasis options
  if (!is.list(oasis)) stop("oasis must be a list")
  
  if (!is.null(oasis$ridge_x) && (oasis$ridge_x < 0)) {
    stop("ridge_x must be non-negative")
  }
  if (!is.null(oasis$ridge_b) && (oasis$ridge_b < 0)) {
    stop("ridge_b must be non-negative")
  }
  if (!is.null(oasis$ridge_mode)) {
    oasis$ridge_mode <- match.arg(oasis$ridge_mode, c("absolute", "fractional"))
  }
  # legacy oasis$whiten is dropped; use prewhiten instead

  # Default intercept (aligns with non-OASIS paths)
  if (is.null(Z) && isTRUE(oasis$add_intercept %||% TRUE)) {
    # If using design_spec and there are multiple runs, prefer run-wise intercepts
    if (!is.null(oasis$design_spec)) {
      bl <- tryCatch(fmrihrf::blocklens(oasis$design_spec$sframe), error = function(e) NULL)
      if (!is.null(bl) && length(bl) > 1L) {
        runs <- rep(seq_along(bl), bl)
        Z <- stats::model.matrix(~ 0 + factor(runs))
        colnames(Z) <- paste0("run", seq_along(bl))
      } else {
        Z <- matrix(1, n_time, 1)
        colnames(Z) <- "Intercept"
      }
    } else {
      Z <- matrix(1, n_time, 1)
      colnames(Z) <- "Intercept"
    }
  }

  # 1) Build trial-wise design if needed (from fmrihrf)
  X_other <- NULL
  if (is.null(X)) {
    if (is.null(oasis$design_spec)) {
      stop("Either X or oasis$design_spec must be provided")
    }
    
    # If user supplied an HRF grid, pick best HRF now
    if (!is.null(oasis$design_spec$hrf_grid)) {
      conf_for_grid <- cbind(if (!is.null(Z)) Z, if (!is.null(Nuisance)) Nuisance)
      oasis$design_spec$cond$hrf <- .oasis_pick_hrf_lwu_fast(
        Y, oasis$design_spec, oasis$design_spec$hrf_grid, 
        confounds = conf_for_grid, block_cols = as.integer(oasis$block_cols %||% 4096L)
      )
    }
    
    built <- .oasis_build_X_from_events(oasis$design_spec)
    X       <- built$X_trials
    X_other <- built$X_other
    # Auto-detect K from the built design
    if (is.null(oasis$K) && !is.null(built$K)) {
      oasis$K <- built$K
    }
  }

  # 2) Detect K (basis dimension)
  K <- oasis$K %||% {
    detected_K <- 1L
    if (!is.null(oasis$design_spec$cond$hrf)) {
      detected_K <- tryCatch(fmrihrf::nbasis(oasis$design_spec$cond$hrf), error = function(e) 1L)
    } else if (!is.null(oasis$ntrials) && !is.null(X)) {
      N <- ncol(X); ntr <- as.integer(oasis$ntrials)
      if (N %% ntr != 0L) stop(sprintf("ncol(X)=%d is not divisible by ntrials=%d", N, ntr))
      detected_K <- as.integer(N / ntr)
    } else if (!is.null(X)) {
      N <- ncol(X)
      for (Kcand in c(2L,3L,4L,5L,6L,8L,10L,12L)) {
        if (N %% Kcand != 0L) next
        ntr <- N / Kcand
        trials <- seq_len(min(ntr, 8L))
        ok <- TRUE
        for (i in trials) {
          idx <- ((i-1L)*Kcand + 1L):(i*Kcand)
          B <- X[, idx, drop=FALSE]
          G <- crossprod(B)
          Dn <- 1/sqrt(pmax(diag(G), .Machine$double.eps))
          Cn <- diag(Dn) %*% G %*% diag(Dn)
          if (mean(abs(Cn[upper.tri(Cn)])) < 0.5) { ok <- FALSE; break }
        }
        if (ok) { detected_K <- as.integer(Kcand); break }
      }
    }
    detected_K
  }
  K <- as.integer(K)

  # 3) Nuisance design used for projection: Z + Nuisance + (aggregates of other conditions)
  N_nuis <- cbind(if (!is.null(Z)) Z, if (!is.null(Nuisance)) Nuisance, X_other)

  # 4) Whitening hook (optional)
  if (!is.null(prewhiten) && is.list(prewhiten) && (prewhiten$method %||% "none") != "none") {
    whitened <- .prewhiten_data(Y, X, NULL, N_nuis, prewhiten)
    Y <- whitened$Y_whitened
    X <- whitened$X_whitened
    if (!is.null(whitened$Nuisance_whitened)) N_nuis <- whitened$Nuisance_whitened
  }

  # 4.5) Optional VOXHRF branch: per-voxel HRF via tiny ridge, then LSS with voxel HRFs
  if (isTRUE(oasis$hrf_mode %in% c("voxel_ridge", "voxhrf"))) {
    if (is.null(oasis$design_spec)) {
      stop("hrf_mode='", oasis$hrf_mode, "' requires oasis$design_spec (events + HRF basis)")
    }
    # By this point, if X was NULL, it was built in step (1) and whitened once in step (4).
    # Do not whiten or rebuild here to avoid double-whitening.

    # Fit voxel-wise HRF weights fast, then compute per-trial betas with C++ voxel-HRF engine
    vhrf <- .estimate_voxel_hrf_fast(
      Y = Y,
      X_trials = X,
      design_spec = oasis$design_spec,
      N_nuis = N_nuis,
      K = oasis$K %||% .detect_basis_dimension(X),
      lambda_shape = oasis$lambda_shape %||% 0,
      mu_rough     = oasis$mu_rough %||% 0,
      ref_hrf      = oasis$ref_hrf %||% NULL,
      shrink_global = oasis$shrink_global %||% 0,
      orient_ref   = isTRUE(oasis$orient_ref %||% TRUE)
    )

    events_df <- data.frame(
      onset     = oasis$design_spec$cond$onsets,
      duration  = oasis$design_spec$cond$duration %||% 0,
      condition = "cond"
    )

    # Build TR-level convolved basis list directly from X_trials for exact parity
    built <- .oasis_build_X_from_events(oasis$design_spec)
    X_trials <- built$X_trials
    Kloc <- oasis$K %||% built$K %||% 1L
    idx_by_basis <- lapply(seq_len(Kloc), function(k) seq.int(k, ncol(X_trials), by = Kloc))
    basis_convolved <- lapply(idx_by_basis, function(idx) X_trials[, idx, drop = FALSE])

    # Residualize Y and basis by N_nuis (FWL), keep intercept in Z for design
    n_time <- nrow(Y)
    if (!is.null(N_nuis) && ncol(N_nuis) > 0) {
      qrN <- qr(N_nuis)
      Y_res <- qr.resid(qrN, Y)
      basis_convolved <- lapply(basis_convolved, function(Dk) qr.resid(qrN, Dk))
    } else {
      Y_res <- Y
    }

    # Scale weights by ref_norm so the combined TR design matches canonical OASIS
    coeff_use <- vhrf$coefficients
    if (!is.null(vhrf$ref_norm) && is.finite(vhrf$ref_norm) && vhrf$ref_norm != 0)
      coeff_use <- coeff_use * as.numeric(vhrf$ref_norm)

    beta_mat <- lss_engine_vox_hrf_arma(
      Y = Y_res,
      coeffs = coeff_use,
      basis_convolved = basis_convolved,
      Z = matrix(1, n_time, 1L)
    )
    return(.set_beta_dimnames(beta_mat,
                              trial_names = .default_trial_names(nrow(beta_mat)),
                              voxel_names = colnames(Y)))
  }

  # 5) Branch based on K
  if (K == 1L) {
    # Single-basis path
    pre   <- oasis_precompute_design(X, if (is.null(N_nuis)) matrix(0, nrow(X), 0) else N_nuis)
    mats  <- oasis_AtY_SY_blocked(pre$A, pre$s_all, pre$Q, Y, as.integer(oasis$block_cols %||% 4096L))
    
    # Resolve ridge
    lam   <- .oasis_resolve_ridge(pre, oasis$ridge_x %||% 0, oasis$ridge_b %||% 0,
                                  oasis$ridge_mode %||% "absolute", K = 1L)
    
    # Compute betas
    if (isTRUE(oasis$return_se)) {
      bg <- oasis_betas_gammas(mats$N_Y, mats$S_Y, pre$d, pre$alpha, pre$s,
                               ridge_x = lam$lx, ridge_b = lam$lb)
      B <- bg$beta
    } else {
      B <- oasis_betas_closed_form(mats$N_Y, mats$S_Y, pre$d, pre$alpha, pre$s,
                                   ridge_x = lam$lx, ridge_b = lam$lb)
    }
  } else {
    # Multi-basis path (K > 1)
    pre  <- oasisk_precompute_design(X, if (is.null(N_nuis)) matrix(0, nrow(X), 0) else N_nuis, K)
    mats <- oasisk_products(pre$A, pre$S, pre$Q, Y, as.integer(oasis$block_cols %||% 4096L))
    
    # Resolve ridge
    lam  <- .oasis_resolve_ridge(pre, oasis$ridge_x %||% 0, oasis$ridge_b %||% 0,
                                 oasis$ridge_mode %||% "absolute", K = K)
    
    # Compute betas
    B <- oasisk_betas(pre$D, pre$C, pre$E, mats$N1, mats$SY, 
                     ridge_x = lam$lx, ridge_b = lam$lb)
  }

  # 6) Default: return the bare matrix (back-compat with fmrilss::lss)
  name_rows_cols <- function(Beta) {
    .set_beta_dimnames(Beta,
                       trial_names = .default_trial_names(nrow(Beta)),
                       voxel_names = colnames(Y))
  }
  if (!isTRUE(oasis$return_se) && !isTRUE(oasis$return_diag)) {
    return(name_rows_cols(B))
  }

  # 7) Optional: SEs and diagnostics
  out <- list(beta = name_rows_cols(B))
  if (isTRUE(oasis$return_diag)) {
    if (K == 1L) {
      out$diag <- list(d = pre$d, alpha = pre$alpha, s = pre$s)
    } else {
      out$diag <- list(D = pre$D, C = pre$C, E = pre$E)
    }
  }
  if (isTRUE(oasis$return_se)) {
    # Nuisance rank from the whitened nuisance actually used
    nuis_rank <- if (!is.null(N_nuis) && ncol(N_nuis) > 0L) qr(N_nuis)$rank else 0L
    if (K == 1L) {
      # Ensure RY_norm2 exists
      RY2 <- if (!is.null(mats$RY_norm2)) mats$RY_norm2 else {
        if (!is.null(pre$Q) && nrow(pre$Q) == n_time) {
          Ry <- Y - pre$Q %*% crossprod(pre$Q, Y)
          colSums(Ry^2)
        } else if (!is.null(N_nuis) && ncol(N_nuis) > 0L) {
          Qn <- qr.Q(qr(N_nuis))
          Ry <- Y - Qn %*% crossprod(Qn, Y)
          colSums(Ry^2)
        } else {
          colSums(Y^2)
        }
      }
      dof <- max(1L, n_time - nuis_rank - 2L)
      out$se <- .oasis_se_from_norms(pre$d, pre$alpha, pre$s,
                                     lam$lx, lam$lb,
                                     RY2, bg$beta, bg$gamma, 
                                     mats$N_Y, mats$S_Y, dof)
    } else {
      dof <- max(1L, n_time - nuis_rank - 2L*K)
      RY_norm2 <- if (!is.null(mats$RY_norm2)) mats$RY_norm2 else oasisk_compute_RY_norm2(pre$Q, Y)
      result <- oasisk_betas_se(pre$D, pre$C, pre$E, mats$N1, mats$SY,
                                RY_norm2, ridge_x = lam$lx, ridge_b = lam$lb)
      out$beta <- name_rows_cols(result$beta)
      out$se   <- result$se
    }
  }
  out
}

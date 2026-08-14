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
#'    - K, ntrials, trial_basis_map: required identity contract when a raw
#'      multi-basis X is supplied. The map has columns column, trial, and basis.
#'    - ridge_mode: "fractional" (default) or "absolute"
#'    - ridge_x, ridge_b: nonnegative ridge on the `a_j` / `b_j` Gram
#'      (default 0.05 in fractional mode)
#'    - block_cols: voxel block size (default 4096)
#'    - return_se: logical (default FALSE)
#'    - return_diag: logical (default FALSE)
#' @param prewhiten list of prewhitening options using fmriAR (see ?lss for
#'   details).  The legacy oasis$whiten field is ignored; use this instead.
#'
#' @return by default: (N_trials x V) matrix of betas; if `return_se` or `return_diag`, a list
#' @keywords internal
#' @noRd
.lss_oasis <- function(Y, X = NULL, Z = NULL, Nuisance = NULL, oasis = list(), prewhiten = NULL) {
  X_was_supplied <- !is.null(X)
  X_col_metadata <- attr(X, "col_metadata")
  # Coerce to base matrices early & validate ---------------------------------
  Y        <- .as_base_matrix(Y)
  X        <- .as_base_matrix(X)
  Z        <- .as_base_matrix(Z)
  Nuisance <- .as_base_matrix(Nuisance)

  if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix")
  if (nrow(Y) < 1L || ncol(Y) < 1L) {
    stop("Y must have at least one timepoint and one voxel")
  }
  if (any(!is.finite(Y))) stop("Y contains non-finite values")

  voxel_names <- .validate_or_default_names(
    colnames(Y), ncol(Y), "Voxel_", "Y column names"
  )
  colnames(Y) <- voxel_names

  n_time <- nrow(Y)
  V      <- ncol(Y)

  if (!is.null(X)) {
    if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
    if (nrow(X) < 1L || ncol(X) < 1L) {
      stop("X must have at least one timepoint and one trial regressor")
    }
    if (nrow(X) != n_time) stop("X must have the same number of rows as Y")
    if (any(!is.finite(X))) stop("X contains non-finite values")
    if (!is.null(colnames(X))) {
      .validate_or_default_names(colnames(X), ncol(X), "Trial_", "X column names")
    }
  }
  
  if (!is.null(Z)) {
    if (!is.matrix(Z) || !is.numeric(Z)) stop("Z must be a numeric matrix")
    if (nrow(Z) != n_time) stop("Z must have the same number of rows as Y")
    if (any(!is.finite(Z))) stop("Z contains non-finite values")
  }
  
  if (!is.null(Nuisance)) {
    if (!is.matrix(Nuisance) || !is.numeric(Nuisance)) stop("Nuisance must be a numeric matrix")
    if (nrow(Nuisance) != n_time) stop("Nuisance must have the same number of rows as Y")
    if (any(!is.finite(Nuisance))) stop("Nuisance contains non-finite values")
  }
  
  # Validate oasis options
  .validate_option_names(oasis, .oasis_option_names(), "oasis")
  for (nm in c("return_se", "return_diag", "add_intercept", "infer_K_from_X", "orient_ref")) {
    if (!is.null(oasis[[nm]])) {
      oasis[[nm]] <- .as_scalar_logical(oasis[[nm]], paste0("oasis$", nm))
    }
  }
  oasis$return_se <- oasis$return_se %||% FALSE
  oasis$return_diag <- oasis$return_diag %||% FALSE
  oasis$add_intercept <- oasis$add_intercept %||% TRUE
  if (!is.null(oasis$design_spec) && !is.list(oasis$design_spec)) {
    stop("oasis$design_spec must be a list", call. = FALSE)
  }
  if (!is.null(oasis$hrf_mode)) {
    oasis$hrf_mode <- match.arg(oasis$hrf_mode, c("voxel_ridge", "voxhrf"))
  }
  block_cols <- .as_positive_integer(oasis$block_cols %||% 4096L, "block_cols")
  if (!is.null(oasis$K)) {
    oasis$K <- .as_positive_integer(oasis$K, "K")
  }
  if (!is.null(oasis$ntrials)) {
    oasis$ntrials <- .as_positive_integer(oasis$ntrials, "ntrials")
  }
  if (isTRUE(oasis$infer_K_from_X %||% FALSE)) {
    stop("infer_K_from_X is not a certified path; supply both oasis$K and oasis$ntrials for a raw multi-basis X", call. = FALSE)
  }
  
  if (!is.null(oasis$ridge_x)) {
    oasis$ridge_x <- .as_nonnegative_scalar(oasis$ridge_x, "ridge_x")
  }
  if (!is.null(oasis$ridge_b)) {
    oasis$ridge_b <- .as_nonnegative_scalar(oasis$ridge_b, "ridge_b")
  }
  if (!is.null(oasis$ridge_mode)) {
    oasis$ridge_mode <- match.arg(oasis$ridge_mode, c("absolute", "fractional"))
  }
  for (nm in c("lambda_shape", "mu_rough", "shrink_global")) {
    if (!is.null(oasis[[nm]])) {
      oasis[[nm]] <- .as_nonnegative_scalar(oasis[[nm]], paste0("oasis$", nm))
    }
  }
  if (!is.null(oasis$shrink_global) && oasis$shrink_global > 1) {
    stop("oasis$shrink_global must be between 0 and 1", call. = FALSE)
  }
  if (isTRUE(oasis$return_se %||% FALSE) &&
      !is.null(oasis$design_spec$hrf_grid)) {
    stop(
      "OASIS standard errors are not calibrated after data-adaptive HRF grid selection; return_se must be FALSE",
      call. = FALSE
    )
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
        confounds = conf_for_grid, block_cols = block_cols
      )
    }
    
    built <- .oasis_build_X_from_events(oasis$design_spec)
    built_ntrials <- length(oasis$design_spec$cond$onsets)
    if (!is.null(oasis$K) && oasis$K != built$K) {
      stop(
        sprintf("oasis$K=%d disagrees with design_spec HRF basis dimension %d",
                oasis$K, built$K),
        call. = FALSE
      )
    }
    if (!is.null(oasis$ntrials) && oasis$ntrials != built_ntrials) {
      stop(
        sprintf("oasis$ntrials=%d disagrees with design_spec event count %d",
                oasis$ntrials, built_ntrials),
        call. = FALSE
      )
    }
    X       <- built$X_trials
    X_other <- built$X_other
    if (built$K > 1L) oasis$trial_basis_map <- built$trial_basis_map
    # Auto-detect K from the built design
    if (is.null(oasis$K) && !is.null(built$K)) {
      oasis$K <- built$K
    }
  }

  if (X_was_supplied && is.null(oasis$K) && is.null(oasis$ntrials) &&
      is.null(oasis$trial_basis_map) && !is.null(X_col_metadata)) {
    inferred <- .infer_fmridesign_trial_basis_map(X, X_col_metadata)
    if (!is.null(inferred)) {
      oasis$K <- inferred$K
      oasis$ntrials <- inferred$ntrials
      oasis$trial_basis_map <- inferred$map
    }
  }

  # 2) Detect K (basis dimension)
  K <- oasis$K %||% {
    detected_K <- 1L
    if (!is.null(oasis$design_spec$cond$hrf)) {
      detected_K <- tryCatch(fmrihrf::nbasis(oasis$design_spec$cond$hrf), error = function(e) 1L)
    } else if (!is.null(oasis$ntrials) && !is.null(X)) {
      N <- ncol(X)
      ntr <- .as_positive_integer(oasis$ntrials, "ntrials")
      if (N %% ntr != 0L) stop(sprintf("ncol(X)=%d is not divisible by ntrials=%d", N, ntr))
      detected_K <- as.integer(N / ntr)
    }
    detected_K
  }
  K <- .as_positive_integer(K, "K")

  if (K == 1L && is.null(colnames(X))) {
    colnames(X) <- .default_trial_names(ncol(X))
  }

  if (X_was_supplied && K > 1L) {
    if (is.null(oasis$K) || is.null(oasis$ntrials) || is.null(oasis$trial_basis_map)) {
      stop(
        "Raw multi-basis X requires explicit oasis$K, oasis$ntrials, and ",
        "oasis$trial_basis_map", call. = FALSE
      )
    }
    if (ncol(X) != K * oasis$ntrials) {
      stop(sprintf("ncol(X)=%d must equal oasis$K * oasis$ntrials = %d",
                   ncol(X), K * oasis$ntrials), call. = FALSE)
    }
    mapped <- .validate_oasis_trial_basis_map(
      X, oasis$trial_basis_map, K = K, ntrials = oasis$ntrials
    )
    X <- mapped$X
    oasis$trial_basis_map <- mapped$map
  }

  # 3) Nuisance design used for projection: Z + Nuisance + (aggregates of other conditions)
  N_nuis <- cbind(if (!is.null(Z)) Z, if (!is.null(Nuisance)) Nuisance, X_other)

  # 4) Whitening hook (optional)
  whiten_plan <- NULL
  if (!is.null(prewhiten) && is.list(prewhiten) && (prewhiten$method %||% "none") != "none") {
    whitened <- .prewhiten_data(Y, X, NULL, N_nuis, prewhiten)
    whiten_plan <- whitened$whiten_plan
    Y <- whitened$Y_whitened
    X <- whitened$X_whitened
    if (!is.null(whitened$Nuisance_whitened)) N_nuis <- whitened$Nuisance_whitened
  }

  # Pass a rank-revealing nuisance basis to the C++ kernels. Armadillo's
  # economical QR retains one Q column per input column even when the input is
  # rank deficient, which would otherwise over-project Y and X.
  if (!is.null(N_nuis) && ncol(N_nuis) > 0L) {
    qr_nuis <- qr(N_nuis)
    N_nuis <- if (qr_nuis$rank > 0L) {
      qr.Q(qr_nuis)[, seq_len(qr_nuis$rank), drop = FALSE]
    } else {
      matrix(numeric(), nrow(N_nuis), 0L)
    }
  }

  # 4.5) Optional VOXHRF branch: per-voxel HRF via tiny ridge, then LSS with voxel HRFs
  if (isTRUE(oasis$hrf_mode %in% c("voxel_ridge", "voxhrf"))) {
    if (isTRUE(oasis$return_se %||% FALSE)) {
      stop(
        "Voxel-adaptive HRF uncertainty is not calibrated; return_se must be FALSE",
        call. = FALSE
      )
    }
    if (isTRUE(oasis$return_diag %||% FALSE)) {
      stop(
        "Voxel-adaptive HRF diagnostics are not defined; return_diag must be FALSE",
        call. = FALSE
      )
    }
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
      K = K,
      lambda_shape = oasis$lambda_shape %||% 0,
      mu_rough     = oasis$mu_rough %||% 0,
      ref_hrf      = oasis$ref_hrf %||% NULL,
      shrink_global = oasis$shrink_global %||% 0,
      orient_ref   = isTRUE(oasis$orient_ref %||% TRUE)
    )

    # Split the current design, which is already whitened when prewhitening was
    # requested. Rebuilding from design_spec here would mix whitened responses
    # with unwhitened trial regressors.
    X_trials <- X
    Kloc <- K
    idx_by_basis <- lapply(seq_len(Kloc), function(k) seq.int(k, ncol(X_trials), by = Kloc))
    basis_convolved <- lapply(idx_by_basis, function(idx) X_trials[, idx, drop = FALSE])

    # Residualize Y and basis by N_nuis (FWL), keep intercept in Z for design
    n_time <- nrow(Y)
    if (!is.null(N_nuis) && ncol(N_nuis) > 0) {
      qrN <- qr(N_nuis)
      Y_res <- qr.resid(qrN, Y)
      basis_convolved <- lapply(basis_convolved, function(Dk) qr.resid(qrN, Dk))
      common_fit <- qr.Q(qrN)[, seq_len(qrN$rank), drop = FALSE]
    } else {
      Y_res <- Y
      common_fit <- matrix(numeric(), n_time, 0L)
    }

    # Scale weights by ref_norm so the combined TR design matches canonical OASIS
    coeff_use <- vhrf$coefficients
    if (!is.null(vhrf$ref_norm) && is.finite(vhrf$ref_norm) && vhrf$ref_norm != 0)
      coeff_use <- coeff_use * as.numeric(vhrf$ref_norm)

    beta_mat <- lss_engine_vox_hrf_arma(
      Y = Y_res,
      coeffs = coeff_use,
      basis_convolved = basis_convolved,
      Z = common_fit
    )
    return(.attach_whiten_plan(
      .set_beta_dimnames(
        beta_mat,
        trial_names = .default_trial_names(nrow(beta_mat)),
        voxel_names = colnames(Y)
      ),
      whiten_plan
    ))
  }

  # 5) Branch based on K
  if (K == 1L) {
    # Single-basis path
    pre   <- oasis_precompute_design(X, if (is.null(N_nuis)) matrix(0, nrow(X), 0) else N_nuis)
    mats  <- oasis_AtY_SY_blocked(pre$A, pre$s_all, pre$Q, Y, block_cols)
    
    # Resolve ridge
    lam   <- .oasis_resolve_ridge(
      pre,
      oasis$ridge_x %||% 0.05,
      oasis$ridge_b %||% 0.05,
      oasis$ridge_mode %||% "fractional",
      K = 1L
    )
    single_trial <- length(pre$d) == 1L
    if (single_trial) {
      target_gram <- pre$d[1L] + lam$lx
      if (!is.finite(target_gram) || target_gram <= 0) {
        stop("OASIS target-trial Gram matrix is not positive definite", call. = FALSE)
      }
      B <- mats$N_Y / target_gram
    } else if (lam$lx == 0 && lam$lb == 0) {
      model_rank <- vapply(seq_along(pre$d), function(j) {
        gram <- matrix(c(pre$d[j], pre$alpha[j], pre$alpha[j], pre$s[j]), 2L, 2L)
        qr(gram)$rank
      }, integer(1))
      if (any(model_rank < 2L)) {
        stop("Unpenalized OASIS is undefined for a rank-deficient trial model", call. = FALSE)
      }
    } else {
      penalized_ok <- vapply(seq_along(pre$d), function(j) {
        gram <- matrix(
          c(pre$d[j] + lam$lx, pre$alpha[j],
            pre$alpha[j], pre$s[j] + lam$lb),
          2L, 2L
        )
        !inherits(try(chol(gram), silent = TRUE), "try-error")
      }, logical(1))
      if (any(!penalized_ok)) {
        stop("Penalized OASIS trial Gram matrix is not positive definite", call. = FALSE)
      }
    }

    # Compute betas for ordinary multi-trial models. A one-trial design has
    # no "other trials" regressor, so it was solved above as a one-block model.
    if (!single_trial) {
      if (isTRUE(oasis$return_se)) {
        bg <- oasis_betas_gammas(mats$N_Y, mats$S_Y, pre$d, pre$alpha, pre$s,
                                 ridge_x = lam$lx, ridge_b = lam$lb)
        B <- bg$beta
      } else {
        B <- oasis_betas_closed_form(mats$N_Y, mats$S_Y, pre$d, pre$alpha, pre$s,
                                     ridge_x = lam$lx, ridge_b = lam$lb)
      }
    }
  } else {
    # Multi-basis path (K > 1)
    pre  <- oasisk_precompute_design(X, if (is.null(N_nuis)) matrix(0, nrow(X), 0) else N_nuis, K)
    mats <- oasisk_products(pre$A, pre$S, pre$Q, Y, block_cols)
    
    # Resolve ridge
    lam  <- .oasis_resolve_ridge(
      pre,
      oasis$ridge_x %||% 0.05,
      oasis$ridge_b %||% 0.05,
      oasis$ridge_mode %||% "fractional",
      K = K
    )
    single_trial <- dim(pre$D)[3] == 1L
    if (single_trial) {
      target_gram <- pre$D[, , 1L] + lam$lx * diag(K)
      target_chol <- try(chol(target_gram), silent = TRUE)
      if (inherits(target_chol, "try-error")) {
        stop("OASIS target-trial Gram matrix is not positive definite", call. = FALSE)
      }
      B <- backsolve(
        target_chol,
        forwardsolve(t(target_chol), mats$N1)
      )
    } else if (lam$lx == 0 && lam$lb == 0) {
      model_rank <- vapply(seq_len(dim(pre$D)[3]), function(j) {
        gram <- rbind(cbind(pre$D[, , j], pre$C[, , j]),
                      cbind(t(pre$C[, , j]), pre$E[, , j]))
        qr(gram)$rank
      }, integer(1))
      if (any(model_rank < 2L * K)) {
        stop("Unpenalized OASIS is undefined for a rank-deficient trial model", call. = FALSE)
      }
    } else {
      penalized_ok <- vapply(seq_len(dim(pre$D)[3]), function(j) {
        gram <- rbind(
          cbind(pre$D[, , j] + lam$lx * diag(K), pre$C[, , j]),
          cbind(t(pre$C[, , j]), pre$E[, , j] + lam$lb * diag(K))
        )
        !inherits(try(chol(gram), silent = TRUE), "try-error")
      }, logical(1))
      if (any(!penalized_ok)) {
        stop("Penalized OASIS trial Gram matrix is not positive definite", call. = FALSE)
      }
    }

    if (!single_trial) {
      B <- oasisk_betas(pre$D, pre$C, pre$E, mats$N1, mats$SY,
                       ridge_x = lam$lx, ridge_b = lam$lb)
    }
  }

  # 6) Default: return the bare matrix (back-compat with fmrilss::lss)
  beta_row_names <- if (K == 1L) {
    if (!is.null(colnames(X)) && length(colnames(X)) == nrow(B)) colnames(X) else .default_trial_names(nrow(B))
  } else {
    if (!is.null(colnames(X)) && length(colnames(X)) == nrow(B)) {
      colnames(X)
    } else {
      ntr <- nrow(B) %/% K
      as.vector(t(outer(sprintf("Trial_%d", seq_len(ntr)),
                        sprintf("Basis_%d", seq_len(K)), paste, sep = ":")))
    }
  }
  name_rows_cols <- function(Beta) {
    out <- .set_beta_dimnames(Beta,
                             trial_names = beta_row_names,
                             voxel_names = colnames(Y))
    if (K > 1L && !is.null(oasis$trial_basis_map)) {
      attr(out, "trial_basis_map") <- oasis$trial_basis_map
    }
    out
  }
  if (!isTRUE(oasis$return_se) && !isTRUE(oasis$return_diag)) {
    return(.attach_whiten_plan(name_rows_cols(B), whiten_plan))
  }

  # 7) Optional: SEs and diagnostics
  out <- list(beta = name_rows_cols(B))
  if (isTRUE(oasis$return_diag)) {
    if (K == 1L) {
      out$diag <- list(
        d = stats::setNames(as.numeric(pre$d), beta_row_names),
        alpha = stats::setNames(as.numeric(pre$alpha), beta_row_names),
        s = stats::setNames(as.numeric(pre$s), beta_row_names)
      )
    } else {
      out$diag <- list(D = pre$D, C = pre$C, E = pre$E)
    }
  }
  if (isTRUE(oasis$return_se)) {
    active_whitening <- !is.null(prewhiten) && is.list(prewhiten) &&
      (prewhiten$method %||% "none") != "none"
    if (active_whitening) {
      stop(
        "OASIS standard errors are not calibrated for estimated prewhitening; ",
        "fit without prewhiten or request coefficients only.", call. = FALSE
      )
    }
    if (lam$lx != 0 || lam$lb != 0) {
      stop(
        "Standard errors are currently defined only for unpenalized OASIS. ",
        "Set ridge_x = ridge_b = 0; ridge uncertainty is diagnostic and not returned.",
        call. = FALSE
      )
    }
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
      if (single_trial) {
        dof <- nrow(Y) - nuis_rank - 1L
        if (dof <= 0L) {
          stop("OASIS standard errors require positive residual degrees of freedom", call. = FALSE)
        }
        SSE <- RY2 - 2 * B[1L, ] * mats$N_Y[1L, ] +
          B[1L, ]^2 * pre$d[1L]
        sigma2 <- pmax(SSE / dof, 0)
        out$se <- matrix(sqrt(sigma2 / pre$d[1L]), nrow = 1L)
      } else {
        model_rank <- vapply(seq_along(pre$d), function(j) {
          gram <- matrix(c(pre$d[j], pre$alpha[j], pre$alpha[j], pre$s[j]), 2L, 2L)
          qr(gram)$rank
        }, integer(1))
        if (any(model_rank < 2L)) {
          stop("OASIS standard errors are undefined for a rank-deficient trial model", call. = FALSE)
        }
        dof <- nrow(Y) - nuis_rank - model_rank
        if (any(dof <= 0L)) {
          stop("OASIS standard errors require positive residual degrees of freedom", call. = FALSE)
        }
        out$se <- .oasis_se_from_norms(pre$d, pre$alpha, pre$s,
                                       lam$lx, lam$lb,
                                       RY2, bg$beta, bg$gamma,
                                       mats$N_Y, mats$S_Y, dof)
      }
    } else {
      RY_norm2 <- if (!is.null(mats$RY_norm2)) mats$RY_norm2 else oasisk_compute_RY_norm2(pre$Q, Y)
      if (single_trial) {
        dof <- nrow(Y) - nuis_rank - K
        if (dof <= 0L) {
          stop("OASIS standard errors require positive residual degrees of freedom", call. = FALSE)
        }
        target_inv <- chol2inv(chol(pre$D[, , 1L]))
        SSE <- RY_norm2 - 2 * colSums(B * mats$N1) +
          colSums(B * (pre$D[, , 1L] %*% B))
        sigma2 <- pmax(as.numeric(SSE) / dof, 0)
        out$se <- sqrt(diag(target_inv) %o% sigma2)
      } else {
        model_rank <- vapply(seq_len(dim(pre$D)[3]), function(j) {
          gram <- rbind(cbind(pre$D[, , j], pre$C[, , j]),
                        cbind(t(pre$C[, , j]), pre$E[, , j]))
          qr(gram)$rank
        }, integer(1))
        if (any(model_rank < 2L * K)) {
          stop("OASIS standard errors are undefined for a rank-deficient trial model", call. = FALSE)
        }
        dof <- nrow(Y) - nuis_rank - model_rank
        if (any(dof <= 0L)) {
          stop("OASIS standard errors require positive residual degrees of freedom", call. = FALSE)
        }
        result <- oasisk_betas_se(pre$D, pre$C, pre$E, mats$N1, mats$SY,
                                  RY_norm2, dof,
                                  ridge_x = lam$lx, ridge_b = lam$lb)
        out$se <- name_rows_cols(result$se)
      }
    }
    out$se <- name_rows_cols(out$se)
  }
  .attach_whiten_plan(out, whiten_plan)
}

# Infer a safe multi-basis identity map from an unmodified fmridesign matrix.
# fmridesign currently emits basis-major columns, while fmrilss returns the
# canonical trial-major/basis-within-trial layout.
.infer_fmridesign_trial_basis_map <- function(X, metadata) {
  metadata <- as.data.frame(metadata)
  if (nrow(metadata) != ncol(X) || !"term_tag" %in% names(metadata)) {
    stop("fmridesign col_metadata is inconsistent with X", call. = FALSE)
  }
  term_tags <- unique(as.character(metadata$term_tag))
  term_tags <- term_tags[!is.na(term_tags) & nzchar(term_tags)]
  if (length(term_tags) != 1L) {
    stop(
      "A fmridesign X with multiple event terms cannot be inferred safely; ",
      "use lss_design() or supply an explicit trial-only X and ",
      "oasis$trial_basis_map.",
      call. = FALSE
    )
  }
  column_names <- colnames(X)
  if (is.null(column_names) || anyNA(column_names) || any(!nzchar(column_names)) ||
      anyDuplicated(column_names)) {
    stop("fmridesign X must have complete unique column names", call. = FALSE)
  }
  has_basis_suffix <- grepl("_b[0-9]+$", column_names)
  if (!any(has_basis_suffix)) return(NULL)
  if (!all(has_basis_suffix)) {
    stop("fmridesign X has inconsistent basis suffixes", call. = FALSE)
  }
  basis <- as.integer(sub("^.*_b([0-9]+)$", "\\1", column_names))
  trial_name <- sub("_b[0-9]+$", "", column_names)
  trial_levels <- unique(trial_name)
  K <- max(basis)
  if (!identical(sort(unique(basis)), seq_len(K))) {
    stop("fmridesign X basis suffixes must cover a contiguous sequence", call. = FALSE)
  }
  trial <- match(trial_name, trial_levels)
  key <- paste(trial, basis, sep = "\r")
  if (anyDuplicated(key) || length(key) != length(trial_levels) * K) {
    stop("fmridesign X does not form a complete trial-by-basis mapping", call. = FALSE)
  }
  output_name <- paste(
    trial_name, paste0("basis_", basis), sep = ":"
  )
  list(
    K = K,
    ntrials = length(trial_levels),
    map = data.frame(
      column = column_names,
      trial = trial,
      basis = basis,
      output_name = output_name,
      stringsAsFactors = FALSE
    )
  )
}

.validate_oasis_trial_basis_map <- function(X, map, K, ntrials) {
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(!nzchar(colnames(X))) ||
      anyDuplicated(colnames(X))) {
    stop("Raw multi-basis X must have unique, non-empty column names", call. = FALSE)
  }
  if (!is.data.frame(map) ||
      !all(c("column", "trial", "basis") %in% names(map)) ||
      nrow(map) != ncol(X)) {
    stop(
      "oasis$trial_basis_map must have one row per X column and fields ",
      "column, trial, and basis", call. = FALSE
    )
  }
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  map$column <- as.character(map$column)
  if (anyNA(map$column) || any(!nzchar(map$column)) || anyDuplicated(map$column) ||
      !setequal(map$column, colnames(X))) {
    stop("trial_basis_map$column must identify every X column exactly once", call. = FALSE)
  }
  map <- map[match(colnames(X), map$column), , drop = FALSE]
  trial_numeric <- suppressWarnings(as.numeric(as.character(map$trial)))
  basis_numeric <- suppressWarnings(as.numeric(as.character(map$basis)))
  if (anyNA(trial_numeric) || anyNA(basis_numeric) ||
      any(!is.finite(trial_numeric)) || any(!is.finite(basis_numeric)) ||
      any(trial_numeric != round(trial_numeric)) ||
      any(basis_numeric != round(basis_numeric))) {
    stop("trial_basis_map trial/basis indices must be exact integers", call. = FALSE)
  }
  trial <- as.integer(trial_numeric)
  basis <- as.integer(basis_numeric)
  if (
      !identical(sort(unique(trial)), seq_len(ntrials)) ||
      !identical(sort(unique(basis)), seq_len(K))) {
    stop("trial_basis_map trial/basis indices must cover 1:ntrials and 1:K", call. = FALSE)
  }
  key <- paste(trial, basis, sep = "\r")
  expected <- expand.grid(
    basis = seq_len(K), trial = seq_len(ntrials),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )[, c("trial", "basis")]
  expected_key <- paste(expected$trial, expected$basis, sep = "\r")
  order_idx <- match(expected_key, key)
  if (anyNA(order_idx) || anyDuplicated(key)) {
    stop("trial_basis_map must contain each trial/basis pair exactly once", call. = FALSE)
  }
  map$trial <- trial
  map$basis <- basis
  output_name <- if ("output_name" %in% names(map)) as.character(map$output_name) else {
    paste0("Trial_", trial, ":Basis_", basis)
  }
  if (anyNA(output_name) || any(!nzchar(output_name)) || anyDuplicated(output_name)) {
    stop("trial_basis_map output names must be unique and non-empty", call. = FALSE)
  }
  map$output_name <- output_name
  X_out <- X[, order_idx, drop = FALSE]
  colnames(X_out) <- output_name[order_idx]
  map_out <- map[order_idx, , drop = FALSE]
  rownames(map_out) <- NULL
  list(
    X = X_out,
    map = map_out
  )
}

# File: R/oasis_glue.R
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
#' @param prewhiten list of prewhitening options using fmriAR (see \code{?lss}
#'   and \code{\link{prewhiten_options}} for details).  The legacy
#'   \code{oasis$whiten} field is ignored; use this parameter instead.
#'
#' @return by default: (N_trials x V) matrix of betas; if `return_se` or `return_diag`, a list
#' @keywords internal
.lss_oasis <- function(Y, X = NULL, Z = NULL, Nuisance = NULL, oasis = list(), prewhiten = NULL) {
  # Coerce to base matrices early & validate ---------------------------------
  to_mat <- function(M) {
    if (is.null(M)) return(NULL)
    if (inherits(M, "Matrix")) return(as.matrix(M))
    if (is.data.frame(M))     return(as.matrix(M))
    M
  }
  Y        <- to_mat(Y)
  X        <- to_mat(X)
  Z        <- to_mat(Z)
  Nuisance <- to_mat(Nuisance)

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
    if (is.null(rownames(beta_mat))) rownames(beta_mat) <- sprintf("Trial_%d", seq_len(nrow(beta_mat)))
    if (!is.null(colnames(Y)) && is.null(colnames(beta_mat))) colnames(beta_mat) <- colnames(Y)
    return(beta_mat)
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
    if (is.null(rownames(Beta))) rownames(Beta) <- sprintf("Trial_%d", seq_len(nrow(Beta)))
    if (!is.null(colnames(Y)) && is.null(colnames(Beta))) colnames(Beta) <- colnames(Y)
    Beta
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

# --- Helper: Resolve fractional vs absolute ridge ---
.oasis_resolve_ridge <- function(pre, ridge_x, ridge_b, ridge_mode = "absolute", K = 1L) {
  ridge_mode <- match.arg(ridge_mode, c("absolute","fractional"))
  
  if (ridge_mode == "absolute") {
    return(list(lx = as.numeric(ridge_x), lb = as.numeric(ridge_b)))
  }
  
  # Fractional mode: scale by mean design energy
  if (K == 1L) {
    # Scale by mean design energy (units: a'a and b'b)
    mx <- mean(pre$d)
    mb <- mean(pre$s)
    return(list(lx = as.numeric(ridge_x) * mx,
                lb = as.numeric(ridge_b) * mb))
  } else {
    # Scale by mean per-basis energy
    DD <- pre$D
    EE <- pre$E
    N  <- dim(DD)[3]
    mx <- mean(vapply(seq_len(N), function(j) mean(diag(DD[,,j])), numeric(1)))
    mb <- mean(vapply(seq_len(N), function(j) mean(diag(EE[,,j])), numeric(1)))
    return(list(lx = as.numeric(ridge_x) * mx,
                lb = as.numeric(ridge_b) * mb))
  }
}

# --- Helper: per-trial SEs from 2x2 normal equations and SSE ---
.oasis_se_from_norms <- function(d, alpha, s, ridge_x, ridge_b,
                                 RY_norm2, B, G, N_Y, S_Y, dof) {
  # For each trial j and voxel v:
  #   SSE_jv = ||RY||^2 - 2*(β*n1 + γ*n2) + β^2*d + γ^2*e + 2βγ*c
  #   sigma2_jv = SSE_jv / dof
  #   Var(beta) = sigma2_jv * (G^{-1})_{11}, where G = [[d+lambda_x, c],[c, e+lambda_b]]
  N <- nrow(B)
  V <- ncol(B)
  se <- matrix(NA_real_, N, V)

  for (j in seq_len(N)) {
    dj <- d[j] + ridge_x
    ej <- s[j] + ridge_b
    cj <- alpha[j]
    n1 <- N_Y[j, ]
    n2 <- S_Y - n1
    beta  <- B[j, ]
    gamma <- G[j, ]

    SSE <- RY_norm2 - 2*(beta * n1 + gamma * n2) + 
           (beta^2) * d[j] + (gamma^2) * s[j] + 2*beta*gamma*cj
    sigma2 <- pmax(SSE / dof, 0)
    
    # (G^{-1})_{11} = (e) / (d*e - c^2) with ridge included
    denom <- pmax(dj * ej - cj * cj, .Machine$double.eps)
    g11   <- ej / denom
    se[j, ] <- sqrt(sigma2 * g11)
  }
  se
}

# --- Helper: simple AR(1) whitening (optional) ---
## legacy simple AR(1) whitener removed; use prewhiten instead

# --- Helper: build X (trial-wise) and X_other (aggregates for other conditions) via fmrihrf ---
#' @keywords internal
.oasis_build_X_from_events <- function(spec) {
  if (is.null(spec)) stop("design_spec must be provided when X is NULL.")
  
  sframe <- spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)   # global seconds grid

  # Minimal safety for multi-run: warn if onsets look run‑relative
  # Heuristic: in multi-run, if max(onsets) is no larger than a single run's
  # duration in seconds, users likely supplied run-relative onsets.
  bl <- tryCatch(fmrihrf::blocklens(sframe), error = function(e) NULL)
  if (!is.null(bl) && length(bl) > 1L && !is.null(spec$cond$onsets)) {
    TR <- as.numeric(median(diff(times)))
    run_sec <- max(bl) * TR
    onv <- tryCatch(as.numeric(unlist(spec$cond$onsets)), error = function(e) NA_real_)
    max_onset <- suppressWarnings(max(onv, na.rm = TRUE))
    if (is.finite(max_onset) && max_onset <= run_sec + 1e-6) {
      warning(
        paste0(
          "Onsets appear run-relative but design_spec expects global time.\n",
          "For multi-run designs, prefer lss_design() with fmridesign::event_model(),\n",
          "or convert onsets to global seconds (offset each run)."
        ), call. = FALSE
      )
    }
  }
  
  # Get HRF object
  hrf_obj <- spec$cond$hrf %||% fmrihrf::make_hrf("spmg1")
  
  # Detect K from HRF
  K <- tryCatch(fmrihrf::nbasis(hrf_obj), error = function(e) 1L)
  
  # Build trial-wise design
  fac    <- factor(seq_along(spec$cond$onsets))
  rset   <- fmrihrf::regressor_set(onsets   = spec$cond$onsets,
                                   fac      = fac,
                                   hrf      = hrf_obj,
                                   duration = spec$cond$duration %||% 0,
                                   amplitude= spec$cond$amplitude %||% 1,
                                   span     = spec$cond$span %||% 40,
                                   summate  = TRUE)
  X_trials <- fmrihrf::evaluate(rset, grid = times, 
                                precision = spec$precision %||% 0.1,
                                method = spec$method %||% "conv")
  X_trials <- if (inherits(X_trials, "Matrix")) as.matrix(X_trials) else X_trials

  # Aggregates for "other" conditions (one col each) -> nuisances
  X_other <- NULL
  if (length(spec$others)) {
    X_other <- do.call(cbind, lapply(spec$others, function(oc) {
      rr <- fmrihrf::regressor(onsets   = oc$onsets,
                               hrf      = oc$hrf %||% hrf_obj,
                               duration = oc$duration %||% 0,
                               amplitude= oc$amplitude %||% 1,
                               span     = oc$span %||% 40,
                               summate  = TRUE)
      x_eval <- fmrihrf::evaluate(rr, times,
                                  precision = spec$precision %||% 0.1,
                                  method = spec$method %||% "conv")
      if (inherits(x_eval, "Matrix")) x_eval <- as.matrix(x_eval)
      # For multi-basis HRFs, aggregate columns to a single condition column
      if (is.matrix(x_eval)) rowSums(x_eval) else as.numeric(x_eval)
    }))
    if (is.null(dim(X_other))) X_other <- matrix(X_other, ncol = 1)
  }

  return(list(X_trials = X_trials, X_other = X_other, K = K))
}

# --- Helper: LWU HRF grid selection via matched-filter score ---
#' @keywords internal
.oasis_pick_hrf_lwu_fast <- function(Y, design_spec, hrf_grid, confounds = NULL, 
                                      block_cols = 4096L) {
  # Build aggregate regressor for all trials with each HRF candidate
  sframe <- design_spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)
  onsets <- design_spec$cond$onsets
  
  # Residualize Y against confounds once
  if (!is.null(confounds) && ncol(confounds) > 0) {
    Q <- qr.Q(qr(confounds))
    Y <- Y - Q %*% crossprod(Q, Y)
  }
  
  scores <- numeric(length(hrf_grid))
  
  for (i in seq_along(hrf_grid)) {
    # Build aggregate regressor with this HRF
    r <- fmrihrf::regressor(
      onsets = onsets,
      hrf = hrf_grid[[i]],
      duration = design_spec$cond$duration %||% 0,
      amplitude = design_spec$cond$amplitude %||% 1,
      span = design_spec$cond$span %||% 40,
      summate = TRUE
    )
    x_eval <- fmrihrf::evaluate(r, times, 
                               precision = design_spec$precision %||% 0.1,
                               method = design_spec$method %||% "conv")
    
    # Handle multi-basis HRFs - sum across basis functions for aggregate
    if (is.matrix(x_eval)) {
      x <- rowSums(x_eval)
    } else {
      x <- as.numeric(x_eval)
    }
    
    # Matched-filter score: sum of correlations
    x_norm <- matrix(x / sqrt(sum(x^2)), ncol = 1)
    
    # Process in blocks for memory efficiency
    score <- 0
    for (v_start in seq(1, ncol(Y), by = block_cols)) {
      v_end <- min(v_start + block_cols - 1, ncol(Y))
      Y_block <- Y[, v_start:v_end, drop = FALSE]
      
      # Correlation with each voxel
      cors <- t(x_norm) %*% Y_block / 
              sqrt(colSums(Y_block^2))
      score <- score + sum(abs(cors), na.rm = TRUE)
    }
    scores[i] <- score
  }
  
  # Return HRF with highest score
  best_idx <- which.max(scores)
  hrf_grid[[best_idx]]
}

# --- Helper: detect basis dimension from design matrix ---
.detect_basis_dimension <- function(X) {
  N <- ncol(X)
  if (N <= 1) return(1L)
  
  # Check for common basis dimensions
  for (K in c(2, 3, 4, 5, 6, 8, 10, 12)) {
    if (N %% K == 0) {
      # Check if columns group into K-sized blocks with high within-block correlation
      n_trials <- N / K
      block_cors <- numeric(n_trials)
      
      for (i in seq_len(n_trials)) {
        block_start <- (i - 1) * K + 1
        block_end <- i * K
        block <- X[, block_start:block_end, drop = FALSE]
        
        if (K > 1) {
          cors <- cor(block)
          diag(cors) <- NA
          block_cors[i] <- mean(abs(cors), na.rm = TRUE)
        }
      }
      
      # If blocks show high internal correlation, likely found K
      if (mean(block_cors) > 0.5) {
        return(as.integer(K))
      }
    }
  }
  
  # Default to 1 if no clear pattern
  return(1L)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Internal: Fast per-voxel HRF estimation via aggregated per-basis regressors ---
#' @keywords internal
.estimate_voxel_hrf_fast <- function(Y, X_trials, design_spec, N_nuis = NULL, K = NULL,
                                     lambda_shape = 0, mu_rough = 0,
                                     ref_hrf = NULL, shrink_global = 0, orient_ref = TRUE) {
  # Coerce
  to_mat <- function(M) if (is.null(M)) NULL else (if (inherits(M, "Matrix")) as.matrix(M) else as.matrix(M))
  Y <- to_mat(Y)
  X_trials <- to_mat(X_trials)
  N_nuis <- to_mat(N_nuis)

  # Prefer HRF nbasis for K; fall back to detection when needed
  K_hrf <- tryCatch(fmrihrf::nbasis(design_spec$cond$hrf), error = function(e) NULL)
  if (is.null(K)) {
    K <- if (!is.null(K_hrf)) as.integer(K_hrf) else .detect_basis_dimension(X_trials)
  } else {
    K <- as.integer(K)
  }
  if ((ncol(X_trials) %% K) != 0L) {
    K <- .detect_basis_dimension(X_trials)
  }
  stopifnot(ncol(X_trials) %% K == 0L)

  # 1) Build aggregated per-basis regressors A [T x K] by summing K-blocked columns
  idx_by_basis <- lapply(seq_len(K), function(k) seq.int(k, ncol(X_trials), by = K))
  A <- do.call(cbind, lapply(idx_by_basis, function(idx) rowSums(X_trials[, idx, drop = FALSE])))

  # 2) Residualize Y and A by confounds (Z+Nuisance+others already combined into N_nuis)
  if (!is.null(N_nuis) && ncol(N_nuis) > 0) {
    qrN <- qr(N_nuis)
    Y <- qr.resid(qrN, Y)
    A <- qr.resid(qrN, A)
  }

  # 3) Shape-prior ridge in basis space
  sframe <- design_spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)

  # Evaluate basis kernels on TR grid using a single onset at 0
  span <- design_spec$cond$span %||% (max(times) - min(times)) %||% 40
  r <- fmrihrf::regressor(onsets = 0, hrf = design_spec$cond$hrf, duration = 0, span = span)
  B_time <- fmrihrf::evaluate(r, grid = times, precision = design_spec$precision %||% 0.1,
                              method = design_spec$method %||% "conv")
  if (inherits(B_time, "Matrix")) B_time <- as.matrix(B_time)
  if (is.vector(B_time)) B_time <- cbind(B_time)
  # If evaluated basis dimension disagrees with K, rebuild A accordingly
  K_eval <- ncol(B_time)
  if (K_eval != K) {
    K <- as.integer(K_eval)
    stopifnot(ncol(X_trials) %% K == 0L)
    idx_by_basis <- lapply(seq_len(K), function(k) seq.int(k, ncol(X_trials), by = K))
    A <- do.call(cbind, lapply(idx_by_basis, function(idx) rowSums(X_trials[, idx, drop = FALSE])))
  }

  # Reference HRF projection -> W0 (K-vector)
  if (is.null(ref_hrf)) {
    # Default: use first basis kernel as reference shape
    ref_hrf <- B_time[, 1]
  } else {
    # Coerce to a vector of same length as time grid if needed
    ref_hrf <- as.numeric(ref_hrf)
    if (length(ref_hrf) != nrow(B_time)) {
      # If provided on a different grid, evaluate HRF object if possible; otherwise fallback to first kernel
      ref_hrf <- tryCatch({
        rr <- fmrihrf::regressor(onsets = 0, hrf = ref_hrf, duration = 0, span = span)
        as.numeric(fmrihrf::evaluate(rr, grid = times, precision = design_spec$precision %||% 0.1,
                                     method = design_spec$method %||% "conv"))
      }, error = function(e) B_time[, 1])
    }
  }
  BtB <- crossprod(B_time)                          # K x K
  # Small ridge on BtB to avoid singularity in projection
  W0 <- drop(solve(BtB + 1e-8 * diag(ncol(B_time)), crossprod(B_time, ref_hrf)))

  # Roughness penalty R in basis space: R = (D2 B)^T (D2 B)
  if (mu_rough > 0) {
    D2 <- diff(diag(nrow(B_time)), differences = 2)
    R  <- crossprod(D2 %*% B_time)
  } else {
    R <- matrix(0, ncol(B_time), ncol(B_time))
  }

  # 4) Solve tiny KxK ridge system batched over voxels
  AtA <- crossprod(A)                              # K x K
  AtY <- crossprod(A, Y)                           # K x V
  P   <- AtA + lambda_shape * diag(K) + mu_rough * R
  rhs <- AtY + lambda_shape * W0 %*% matrix(1, 1, ncol(Y))
  W   <- solve(P, rhs)
  W   <- matrix(W, nrow = K, ncol = ncol(Y))        # ensure matrix even when K=1

  # 5) Normalize each voxel's HRF shape to unit energy on the time grid
  # energy_v = sqrt( diag(W^T (B^T B) W) ) computed without forming H = B W
  M <- BtB
  MW <- M %*% W                                   # K x V
  s2 <- sqrt(pmax(colSums(W * MW), .Machine$double.eps))
  W  <- sweep(W, 2L, s2, "/")
  W  <- matrix(W, nrow = K, ncol = ncol(Y))

  ref_norm <- sqrt(sum((B_time %*% W0)^2))

  # 6) Optional global shrinkage for stability
  if (shrink_global > 0) {
    Wbar <- rowMeans(W)
    W <- (1 - shrink_global) * W + shrink_global * matrix(Wbar, nrow = K, ncol = ncol(W))
    W <- matrix(W, nrow = K, ncol = ncol(W))
  }

  # 7) Optional orientation: align with reference HRF so betas are positive-interpretable
  if (isTRUE(orient_ref)) {
    cvec <- as.numeric(crossprod(B_time, ref_hrf))   # length K
    if (length(cvec) != nrow(W)) {
      stop("Orientation reference length mismatch with basis dimension")
    }
    dots <- as.numeric(crossprod(cvec, W))          # length V
    flip <- which(dots < 0)
    if (length(flip)) W[, flip] <- -W[, flip, drop = FALSE]
  }

  structure(list(
    coefficients = W,                          # K x V
    basis       = design_spec$cond$hrf,        # fmrihrf basis object
    conditions  = "cond",
    ref_norm    = as.numeric(ref_norm)
  ), class = "VoxelHRF")
}

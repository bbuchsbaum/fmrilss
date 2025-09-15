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
#'    - whiten: "none" | "ar1" (default "none"); if "ar1", prewhiten Y and design first
#'
#' @return by default: (N_trials x V) matrix of betas; if `return_se` or `return_diag`, a list
#' @keywords internal
.lss_oasis <- function(Y, X = NULL, Z = NULL, Nuisance = NULL, oasis = list()) {
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

  if (!is.matrix(Y)) stop("Y must be a matrix")
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
  if (!is.null(oasis$whiten)) {
    oasis$whiten <- match.arg(tolower(oasis$whiten), c("none", "ar1"))
  }

  # Default intercept (aligns with non-OASIS paths)
  if (is.null(Z) && isTRUE(oasis$add_intercept %||% TRUE)) {
    Z <- matrix(1, n_time, 1)
    colnames(Z) <- "Intercept"
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
    if (!is.null(oasis$design_spec$cond$hrf)) {
      return(tryCatch(fmrihrf::nbasis(oasis$design_spec$cond$hrf), error = function(e) 1L))
    }
    if (!is.null(oasis$ntrials) && !is.null(X)) {
      N <- ncol(X); ntr <- as.integer(oasis$ntrials)
      if (N %% ntr != 0L) stop(sprintf("ncol(X)=%d is not divisible by ntrials=%d", N, ntr))
      return(as.integer(N / ntr))
    }
    if (!is.null(X)) {
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
        if (ok) return(as.integer(Kcand))
      }
    }
    1L
  }
  K <- as.integer(K)

  # 3) Nuisance design used for projection: Z + Nuisance + (aggregates of other conditions)
  N_nuis <- cbind(if (!is.null(Z)) Z, if (!is.null(Nuisance)) Nuisance, X_other)

  # 4) Whitening hook (optional)
  if (isTRUE(tolower(oasis$whiten %||% "none") == "ar1")) {
    w <- .oasis_ar1_whitener(Y, X_nuis = N_nuis)
    Y <- w$Wy
    if (!is.null(N_nuis) && ncol(N_nuis) > 0) N_nuis <- w$W_apply(N_nuis)
    if (!is.null(X))      X      <- w$W_apply(X)
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
.oasis_ar1_whitener <- function(Y, X_nuis = NULL) {
  # Estimate AR(1) parameter and apply variance-preserving whitening
  n_time <- nrow(Y)

  # Optionally deflate nuisance before rho estimation
  Y_for_rho <- if (!is.null(X_nuis) && is.matrix(X_nuis) && ncol(X_nuis) > 0L) {
    Qn <- qr.Q(qr(X_nuis))
    Y - Qn %*% crossprod(Qn, Y)
  } else Y

  y1  <- Y_for_rho[-n_time, , drop=FALSE]
  y2  <- Y_for_rho[-1,     , drop=FALSE]
  rho <- mean(colSums(y1 * y2) / pmax(colSums(y1 * y1), .Machine$double.eps))
  rho <- min(max(rho, -0.98), 0.98)

  Wy        <- Y
  Wy[1, ]   <- sqrt(1 - rho^2) * Y[1, ]
  Wy[-1, ]  <- Y[-1, ] - rho * Y[-n_time, ]

  W_apply <- function(X) {
    if (is.null(X)) return(NULL)
    if (!is.matrix(X)) X <- as.matrix(X)
    WX        <- X
    WX[1, ]   <- sqrt(1 - rho^2) * X[1, ]
    WX[-1, ]  <- X[-1, ] - rho * X[-n_time, ]
    WX
  }
  list(W_apply = W_apply, Wy = Wy, rho = rho)
}

# --- Helper: build X (trial-wise) and X_other (aggregates for other conditions) via fmrihrf ---
#' @keywords internal
.oasis_build_X_from_events <- function(spec) {
  if (is.null(spec)) stop("design_spec must be provided when X is NULL.")
  
  sframe <- spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)   # global seconds grid
  
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
      as.numeric(fmrihrf::evaluate(rr, times, 
                                   precision = spec$precision %||% 0.1,
                                   method = spec$method %||% "conv"))
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

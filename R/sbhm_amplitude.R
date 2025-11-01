# Internal helpers for SBHM amplitude estimation (single-shape GLM / LSS1 / OASIS-K1)

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Residualize columns of M against Z (FWL)
#' @keywords internal
.sbhm_resid <- function(M, Z) {
  if (is.null(Z)) return(M)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (ncol(Z) == 0L) return(M)
  qrZ <- qr(Z)
  qr.resid(qrZ, M)
}

#' Build per-trial basis regressors (each T×r) from SBHM basis + design_spec
#' @keywords internal
.sbhm_build_trial_regs <- function(sbhm, design_spec) {
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  os <- design_spec$cond
  ntrials <- length(os$onsets)
  regs <- vector("list", ntrials)
  for (t in seq_len(ntrials)) {
    rr_t <- fmrihrf::regressor(onsets = os$onsets[t], hrf = hrf_B,
                               duration = os$duration %||% 0,
                               span     = os$span     %||% sbhm$span,
                               summate  = FALSE)
    Xt <- fmrihrf::evaluate(rr_t, grid = sbhm$tgrid,
                            precision = design_spec$precision %||% 0.1,
                            method    = design_spec$method    %||% "conv")
    if (inherits(Xt, "Matrix")) Xt <- as.matrix(Xt)
    if (!is.matrix(Xt)) Xt <- cbind(Xt)
    regs[[t]] <- Xt
  }
  regs
}

#' Optionally prewhiten Y, regs (stacked), intercept and nuisance
#' @keywords internal
.sbhm_prewhiten <- function(Y, regs, Zint, Nuisance, prewhiten) {
  if (is.null(prewhiten) || (is.list(prewhiten) && (prewhiten$method %||% "none") == "none")) {
    return(list(Yw = Y, Zw = Zint, Nw = Nuisance, regs_w = regs, applied = FALSE))
  }
  X_stack <- do.call(cbind, regs)
  pw <- .prewhiten_data(Y, X = X_stack, Z = Zint, Nuisance = Nuisance, prewhiten = prewhiten)
  Xw <- pw$X_whitened
  r <- ncol(regs[[1]]); ntrials <- length(regs)
  regs_w <- vector("list", ntrials)
  for (j in seq_len(ntrials)) {
    cols <- ((j - 1L) * r + 1L):(j * r)
    regs_w[[j]] <- Xw[, cols, drop = FALSE]
  }
  list(Yw = pw$Y_whitened, Zw = pw$Z_whitened, Nw = pw$Nuisance_whitened, regs_w = regs_w, applied = TRUE)
}

#' Resolve ridge value from spec (absolute or fractional)
#' @keywords internal
.sbhm_resolve_ridge <- function(G, ridge) {
  if (is.null(ridge)) return(0)
  if (is.numeric(ridge)) return(as.numeric(ridge))
  mode <- ridge$mode %||% "fractional"
  lam  <- as.numeric(ridge$lambda %||% 0.02)
  if (mode == "absolute") return(lam)
  lam * mean(diag(G))
}

#' Stable linear solve for multiple RHS with tiny ridge
#' @keywords internal
.sbhm_solve <- function(G, B, ridge = 0) {
  G <- as.matrix(G)
  if (!is.matrix(B)) B <- as.matrix(B)
  if (ridge > 0) {
    diag(G) <- diag(G) + ridge
  }
  out <- tryCatch({
    R <- chol(G)
    backsolve(R, forwardsolve(t(R), B))
  }, error = function(e) {
    qr.solve(G, B)
  })
  out
}

#' Single-shape GLM amplitudes given matched coordinates (global LS)
#' @keywords internal
sbhm_amplitude_ls <- function(Y, sbhm, design_spec, alpha_hat,
                              Nuisance = NULL,
                              ridge = list(mode = "fractional", lambda = 0.02),
                              prewhiten = NULL,
                              return_se = FALSE) {
  stopifnot(is.matrix(Y), is.matrix(alpha_hat))
  Tlen <- nrow(Y); V <- ncol(Y)
  r <- nrow(alpha_hat)
  if (ncol(alpha_hat) != V) stop("alpha_hat must be r×V to match Y")

  regs <- .sbhm_build_trial_regs(sbhm, design_spec)
  ntrials <- length(regs)

  # Build nuisance design: intercept + provided + others (aggregates)
  Zint <- matrix(1, Tlen, 1)
  N_mat <- cbind(Zint, if (!is.null(Nuisance)) Nuisance)
  if (is.null(N_mat)) N_mat <- matrix(0, Tlen, 0)

  # Optional prewhitening
  if (!is.null(prewhiten)) {
    pw <- .sbhm_prewhiten(Y, regs, Zint, Nuisance, prewhiten)
    Y <- pw$Yw; regs <- pw$regs_w; Zint <- pw$Zw; N_mat <- cbind(pw$Zw, pw$Nw)
  }
  # Residualize Y once
  Y_res <- .sbhm_resid(Y, N_mat)

  amps <- matrix(0, ntrials, V)
  se_mat <- if (isTRUE(return_se)) matrix(NA_real_, ntrials, V) else NULL
  for (v in seq_len(V)) {
    # Assemble X_v = [regs[[t]] %*% alpha_hat[, v]]_t  (T×ntrials)
    Xcols <- lapply(regs, function(Xt) as.numeric(Xt %*% alpha_hat[, v]))
    Xv <- do.call(cbind, Xcols)
    Xv <- .sbhm_resid(Xv, N_mat)

    G <- crossprod(Xv)
    # Support per-voxel absolute ridge vector
    lam <- if (is.numeric(ridge) && length(ridge) == V) ridge[v] else .sbhm_resolve_ridge(G, ridge)
    B <- crossprod(Xv, Y_res[, v, drop = FALSE])
    beta_v <- as.numeric(.sbhm_solve(G, B, ridge = lam))
    amps[, v] <- beta_v
    if (isTRUE(return_se)) {
      # OLS-style SEs (approximate if ridge>0): se_j = sqrt(sigma2 * (G+lamI)^{-1}_{jj})
      invG <- tryCatch({
        R <- chol(G + diag(lam, nrow(G)))
        inv <- chol2inv(R)
        inv
      }, error = function(e) tryCatch(solve(G + diag(lam, nrow(G))), error = function(e2) diag(NA_real_, nrow(G))))
      resid <- as.numeric(Y_res[, v] - Xv %*% beta_v)
      dof <- max(1L, nrow(Xv) - ncol(Xv))
      sigma2 <- sum(resid^2) / dof
      se_mat[, v] <- sqrt(pmax(0, sigma2) * pmax(0, diag(invG)))
    }
  }
  if (isTRUE(return_se)) list(beta = amps, se = se_mat) else amps
}

#' Single-shape LSS amplitudes (2x2 per trial) given matched coordinates
#' Uses one regressor for the trial of interest (A_j α̂) and one for all other trials
#' (Σ_{t≠j} A_t α̂), with optional fractional ridge on the 2x2 Gram for stability.
#' @keywords internal
sbhm_amplitude_lss1 <- function(Y, sbhm, design_spec, alpha_hat,
                                Nuisance = NULL,
                                ridge_frac = list(x = 0.02, b = 0.02),
                                prewhiten = NULL,
                                return_se = FALSE) {
  stopifnot(is.matrix(Y), is.matrix(alpha_hat))
  Tlen <- nrow(Y); V <- ncol(Y)
  r <- nrow(alpha_hat)
  if (ncol(alpha_hat) != V) stop("alpha_hat must be r×V to match Y")

  # Build basis HRF and per-trial regressors (T×r each)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  os <- design_spec$cond
  ntrials <- length(os$onsets)
  regs <- vector("list", ntrials)
  for (t in seq_len(ntrials)) {
    rr_t <- fmrihrf::regressor(onsets = os$onsets[t], hrf = hrf_B,
                               duration = os$duration %||% 0,
                               span     = os$span     %||% sbhm$span,
                               summate  = FALSE)
    Xt <- fmrihrf::evaluate(rr_t, grid = sbhm$tgrid,
                            precision = design_spec$precision %||% 0.1,
                            method    = design_spec$method    %||% "conv")
    if (inherits(Xt, "Matrix")) Xt <- as.matrix(Xt)
    if (!is.matrix(Xt)) Xt <- cbind(Xt)
    regs[[t]] <- Xt  # T×r
  }

  # Nuisance matrix (intercept + provided + others)
  spec2 <- design_spec
  spec2$cond <- spec2$cond %||% list()
  spec2$cond$hrf <- hrf_B
  built <- .oasis_build_X_from_events(spec2)
  X_other <- built$X_other
  Zint <- matrix(1, Tlen, 1)
  N_mat <- cbind(Zint, if (!is.null(Nuisance)) Nuisance,
                 if (!is.null(X_other)) X_other)
  if (is.null(N_mat)) N_mat <- matrix(0, Tlen, 0)

  # Optional prewhitening
  if (!is.null(prewhiten)) {
    X_stack <- do.call(cbind, regs)
    if (!is.null(X_other)) X_stack <- cbind(X_stack, X_other)
    pw <- .prewhiten_data(Y, X = X_stack, Z = Zint, Nuisance = Nuisance, prewhiten = prewhiten)
    Y <- pw$Y_whitened
    Zint <- pw$Z_whitened
    Nuisance <- pw$Nuisance_whitened
    Xw <- pw$X_whitened
    regs <- lapply(seq_len(ntrials), function(j) {
      rc <- ncol(regs[[1]])
      cols <- ((j - 1L) * rc + 1L):(j * rc)
      Xw[, cols, drop = FALSE]
    })
    if (!is.null(X_other)) {
      X_other <- Xw[, (ntrials * ncol(regs[[1]]) + 1L):ncol(Xw), drop = FALSE]
    }
    N_mat <- cbind(Zint, if (!is.null(Nuisance)) Nuisance,
                   if (!is.null(X_other)) X_other)
  }
  Y_res <- .sbhm_resid(Y, N_mat)

  amps <- matrix(0, ntrials, V)
  se_mat <- if (isTRUE(return_se)) matrix(NA_real_, ntrials, V) else NULL
  for (v in seq_len(V)) {
    # Build trial and others columns per voxel
    x_cols <- lapply(regs, function(Xt) as.numeric(Xt %*% alpha_hat[, v]))
    X_sum  <- Reduce(`+`, x_cols)
    # Residualize columns once
    x_cols_res <- lapply(x_cols, function(x) .sbhm_resid(matrix(x, Tlen, 1), N_mat))
    X_sum_res  <- .sbhm_resid(matrix(X_sum, Tlen, 1), N_mat)
    yv <- Y_res[, v, drop = FALSE]

    for (j in seq_len(ntrials)) {
      x1 <- x_cols_res[[j]]                   # T×1 (trial j)
      x2 <- X_sum_res - x1                    # T×1 (others)
      # 2x2 Gram with fractional ridge
      d  <- sum(x1 * x1)
      e  <- sum(x2 * x2)
      c  <- sum(x1 * x2)
      n1 <- sum(x1 * yv)
      n2 <- sum(x2 * yv)
      lam_x <- (ridge_frac$x %||% 0.01) * d
      lam_b <- (ridge_frac$b %||% 0.01) * e
      D <- (d + lam_x) * (e + lam_b) - c * c
      if (D <= 0) {
        amps[j, v] <- 0
      } else {
        beta_j  <- ((e + lam_b) * n1 - c * n2) / D
        gamma_j <- (((d + lam_x) * n2) - c * n1) / D
        amps[j, v] <- as.numeric(beta_j)
        if (isTRUE(return_se)) {
          # SSE and SE via (G_aug)^{-1}
          # SSE = ||y||^2 - 2(β n1 + γ n2) + β^2 d + γ^2 e + 2βγc
          y_norm2 <- sum(yv * yv)
          sse <- y_norm2 - 2 * (beta_j * n1 + gamma_j * n2) +
                 (beta_j^2) * d + (gamma_j^2) * e + 2 * beta_j * gamma_j * c
          dof <- max(1L, nrow(x1) - 2L)
          sigma2 <- sse / dof
          # (G_aug)^{-1}_{11} = (e + lam_b)/D
          g11 <- (e + lam_b) / D
          se_mat[j, v] <- sqrt(pmax(0, sigma2) * pmax(0, g11))
        }
      }
    }
  }
  if (isTRUE(return_se)) list(beta = amps, se = se_mat) else amps
}

#' OASIS K=1 amplitudes per voxel with matched HRF columns.
#' Returns list(beta, se)
#' @keywords internal
sbhm_amplitude_oasis_k1 <- function(Y, sbhm, design_spec, alpha_hat,
                                    Nuisance = NULL,
                                    ridge_frac = list(x = 0.02, b = 0.02),
                                    prewhiten = NULL, return_se = TRUE) {
  stopifnot(is.matrix(Y), is.matrix(alpha_hat))
  V <- ncol(Y)
  regs <- .sbhm_build_trial_regs(sbhm, design_spec)
  ntrials <- length(regs)
  BETA <- matrix(NA_real_, ntrials, V)
  SE   <- if (isTRUE(return_se)) matrix(NA_real_, ntrials, V) else NULL
  for (v in seq_len(V)) {
    Xv <- do.call(cbind, lapply(regs, function(Xt) as.numeric(Xt %*% alpha_hat[, v])))
    res <- .lss_oasis(
      Y = Y[, v, drop = FALSE], X = Xv, Z = NULL, Nuisance = Nuisance,
      oasis = list(K = 1L,
                   ridge_mode = "fractional",
                   ridge_x = ridge_frac$x %||% 0,
                   ridge_b = ridge_frac$b %||% 0,
                   return_se = isTRUE(return_se),
                   add_intercept = TRUE),
      prewhiten = prewhiten
    )
    if (isTRUE(return_se)) {
      BETA[, v] <- res$beta
      SE[, v]   <- res$se
    } else {
      BETA[, v] <- res
    }
  }
  list(beta = BETA, se = SE)
}

# --- Helpers for gating/ridge tuning (internal) ------------------------------

#' Auto-set cond_gate threshold from observed diagnostics
#' @keywords internal
.sbhm_auto_gate_threshold <- function(diag_values, metric = c("kappa","rho"), q = 0.99, hard_min = NULL) {
  metric <- match.arg(metric)
  thr <- as.numeric(stats::quantile(diag_values, probs = q, na.rm = TRUE))
  if (!is.null(hard_min)) thr <- max(thr, hard_min)
  thr
}

#' Adaptive fractional ridge per voxel from conditioning
#' @keywords internal
.sbhm_adaptive_ridge_gls <- function(kappa_vec, base = 0.02, k0 = 1000, max_lam = 0.08) {
  # Simple linear scale: base * (kappa/k0), clipped
  lam <- base * (kappa_vec / k0)
  pmin(pmax(lam, base), max_lam)
}

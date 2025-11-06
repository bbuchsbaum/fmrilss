#' End-to-End LSS with Shared-Basis HRF Matching (SBHM)
#'
#' Orchestrates the SBHM pipeline: (1) prepass aggregate fit in the learned
#' shared basis, (2) cosine matching to a library of HRFs represented in the
#' same basis, (3) Least Squares Separate (OASIS) with the SBHM basis to obtain
#' trial-wise r-dimensional coefficients, and (4) projection of those
#' coefficients onto the matched coordinates to produce scalar amplitudes.
#'
#' @param Y Numeric matrix T×V of fMRI time series.
#' @param sbhm SBHM object from `sbhm_build()`.
#' @param design_spec List for design construction (same as `oasis$design_spec`):
#'   `list(sframe=..., cond=list(onsets=..., duration=0, span=...), others=list(...))`.
#'   The HRF in `cond$hrf` is ignored and replaced with the SBHM basis HRF.
#' @param Nuisance Optional T×P nuisance regressors.
#' @param prewhiten Optional prewhitening options (see `?lss`).
#' @param prepass Optional list forwarded to `sbhm_prepass()` (e.g., ridge, data_fac).
#' @param match Optional list forwarded to `sbhm_match()` (e.g., shrink, topK, whiten, orient_ref).
#'   Additional fields handled here:
#'   - `soft_blend` logical (default FALSE): when `topK > 1`, blend the top-K
#'     library coordinates per voxel using softmax weights returned by `sbhm_match()`.
#'     If `blend_margin` is provided, blending is only applied to voxels with
#'     `margin < blend_margin`; others use the hard top-1 assignment.
#'   - `blend_margin` optional numeric threshold on the matching margin for
#'     conditional blending.
#' @param oasis Optional list forwarded to `lss(..., method="oasis")`. `K` is set to
#'   `ncol(sbhm$B)` if not provided, and `design_spec` is injected automatically.
#' @param amplitude List controlling the scalar amplitude stage. Fields:
#'   - `method`: one of "global_ls" (default), "lss1", "oasis_voxel".
#'   - `ridge`: for `global_ls`, either numeric (absolute) or list(mode, lambda).
#'   - `ridge_frac`: for `lss1`/`oasis_voxel`, list(x, b) fractional ridge.
#'   - `cond_gate`: optional auto-fallback rule, e.g., list(metric="rho", thr=0.999, fallback="lss1").
#' @param return One of `"amplitude"`, `"coefficients"`, or `"both"` (default `"amplitude"`).
#'
#' @return A list with components:
#'   - `amplitude` ntrials×V matrix (when requested)
#'   - `coeffs_r` r×ntrials×V array of trial-wise coefficients (when requested)
#'   - `matched_idx` length-V integer indices into the library
#'   - `margin` length-V confidence margins (top1 - top2 cosine)
#'   - `alpha_coords` r×V matched coordinates per voxel
#'   - `diag` list with `r`, `ntrials`, and `times`
#'
#' @examples
#' \dontrun{
#'   library(fmrihrf)
#'   set.seed(3)
#'   Tlen <- 180; V <- 4
#'   sframe <- sampling_frame(blocklens = Tlen, TR = 1)
#'   # Simple HRF library matrix
#'   H <- cbind(exp(-seq(0, 30, length.out = Tlen)/5),
#'              exp(-seq(0, 30, length.out = Tlen)/7))
#'   sbhm <- sbhm_build(library_H = H, r = 4, sframe = sframe, normalize = TRUE)
#'   # Events
#'   onsets <- seq(8, 140, by = 12)
#'   design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
#'   # Simulate Y
#'   hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
#'   rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
#'   Xr <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
#'   alpha_true <- rnorm(ncol(sbhm$B))
#'   Y <- matrix(rnorm(Tlen*V, sd = .6), Tlen, V)
#'   Y[,1] <- Y[,1] + Xr %*% alpha_true
#'   out <- lss_sbhm(Y, sbhm, design_spec)
#'   names(out)
#' }
#'
#' @export
lss_sbhm <- function(Y, sbhm, design_spec,
                     Nuisance = NULL,
                     prewhiten = NULL,
                     prepass = list(),
                     match = list(shrink = list(tau = 0, ref = NULL, snr = NULL),
                                  topK = 1, whiten = TRUE, orient_ref = TRUE),
                     oasis = list(),
                     amplitude = list(method = "global_ls",
                                      ridge = list(mode = "fractional", lambda = 0.02),
                                      ridge_frac = list(x = 0.02, b = 0.02),
                                      cond_gate = NULL,
                                      adaptive = list(enable = FALSE, base = 0.02, k0 = 1000, max = 0.08),
                                      return_se = FALSE),
                     return = c("amplitude", "coefficients", "both")) {

  return <- match.arg(return)
  stopifnot(is.matrix(Y), is.list(sbhm), is.list(design_spec))
  r <- ncol(sbhm$B)

  # 1) Prepass aggregate fit in shared basis ---------------------------------
  pre <- sbhm_prepass(
    Y = Y,
    sbhm = sbhm,
    design_spec = design_spec,
    Nuisance = Nuisance,
    prewhiten = prewhiten,
    ridge = prepass$ridge %||% NULL,
    data_fac = prepass$data_fac %||% NULL
  )

  # 2) Matching voxel shapes to library in coefficient space ------------------
  m <- sbhm_match(
    beta_bar   = pre$beta_bar,
    S          = sbhm$S,
    A          = sbhm$A,
    shrink     = match$shrink %||% list(tau = 0, ref = sbhm$ref$alpha_ref, snr = NULL),
    topK       = match$topK %||% 1,
    whiten     = if (is.null(match$whiten)) TRUE else isTRUE(match$whiten),
    orient_ref = if (is.null(match$orient_ref)) TRUE else isTRUE(match$orient_ref)
  )

  # Hard assignment by default; optional soft blending across top-K
  alpha_hat <- m$alpha_hat
  if ((match$topK %||% 1) > 1 && isTRUE(match$soft_blend %||% FALSE)) {
    V <- ncol(Y)
    r <- ncol(sbhm$B)
    alpha_soft <- matrix(0, nrow = r, ncol = V)
    blend_thr <- match$blend_margin %||% NULL
    for (v in seq_len(V)) {
      # If a margin threshold is given and margin is high, keep hard assignment
      # Supports +/-Inf thresholds: margin >= -Inf is always TRUE (keep hard),
      # margin >= +Inf is always FALSE (always blend)
      if (!is.null(blend_thr) && (m$margin[v] >= blend_thr)) {
        alpha_soft[, v] <- m$alpha_hat[, v]
      } else {
        idx <- m$topK_idx[, v]
        w   <- m$weights[, v]
        alpha_soft[, v] <- as.numeric(sbhm$A[, idx, drop = FALSE] %*% w)
      }
    }
    alpha_hat <- alpha_soft
  }

  # 3) OASIS LSS with the SBHM basis (K=r) -----------------------------------
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  spec  <- design_spec
  spec$cond <- spec$cond %||% list()
  spec$cond$hrf <- hrf_B
  os   <- oasis
  os$design_spec <- spec
  os$K <- os$K %||% r

  BetaMat <- lss(
    Y = Y,
    X = NULL,
    Z = NULL,
    Nuisance = Nuisance,
    method = "oasis",
    oasis = os,
    prewhiten = prewhiten
  )
  # BetaMat dims: (ntrials*K) × V with K=r
  if (!is.matrix(BetaMat)) stop("Unexpected return from lss(method='oasis')")
  V <- ncol(BetaMat)
  if ((nrow(BetaMat) %% r) != 0L) stop("Row count not divisible by r in OASIS betas")
  ntrials <- as.integer(nrow(BetaMat) / r)
  beta_rt <- array(BetaMat, dim = c(r, ntrials, V))

  # 4) Scalar amplitudes with method + optional auto-fallback -----------------
  amp_method <- match.arg(amplitude$method %||% "global_ls",
                          c("global_ls","lss1","oasis_voxel"))
  ridge_gls  <- amplitude$ridge %||% list(mode = "fractional", lambda = 0.02)
  ridge_frac <- amplitude$ridge_frac %||% list(x = 0.02, b = 0.02)

  # Conditioning diagnostics in (optionally) whitened space
  regs_diag <- .sbhm_build_trial_regs(sbhm, design_spec)
  Zint_diag <- matrix(1, nrow(Y), 1)
  if (!is.null(prewhiten)) {
    pw <- .sbhm_prewhiten(Y, regs_diag, Zint_diag, Nuisance, prewhiten)
    Y_diag <- pw$Yw; regs_diag <- pw$regs_w; Zint_diag <- pw$Zw; Nuis_diag <- pw$Nw
  } else {
    Y_diag <- Y; Nuis_diag <- Nuisance
  }
  N_mat_diag <- cbind(Zint_diag, if (!is.null(Nuis_diag)) Nuis_diag)
  V <- ncol(Y)
  diag_rho <- rep(NA_real_, V); diag_kappa <- rep(NA_real_, V)
  .cond_metrics <- function(regs, alpha_v, Nmat) {
    # Build trial-by-trial design for this voxel’s matched shape
    Xv <- do.call(cbind, lapply(regs, function(Xt) as.numeric(Xt %*% alpha_v)))
    Xv <- .sbhm_resid(Xv, Nmat)

    # Conditioning via singular values (scale-invariant)
    xs <- sweep(Xv, 2L, sqrt(colSums(Xv^2)) + 1e-12, "/")
    sv <- svd(xs, nu = 0, nv = 0)$d
    kappa <- if (length(sv) > 1L) max(sv) / pmax(min(sv), 1e-8) else 1

    # Robust overlap metric: maximum pairwise absolute cosine between columns
    # This tends to be higher (and more stable) under heavy overlap than
    # correlating with the sum of the other trials.
    if (ncol(Xv) <= 1L) {
      rho_max <- 0
    } else {
      G <- crossprod(Xv)                                 # trial×trial inner products
      d <- sqrt(pmax(1e-12, diag(G)))
      # Normalized correlation/cosine matrix
      C <- G / outer(d, d, "*")
      diag(C) <- 0
      rho_max <- max(abs(C[upper.tri(C)]), na.rm = TRUE)
      if (!is.finite(rho_max)) rho_max <- 0
    }

    c(rho = rho_max, kappa = kappa)
  }
  for (v in seq_len(V)) {
    met <- .cond_metrics(regs_diag, alpha_hat[, v], N_mat_diag)
    diag_rho[v]   <- met[["rho"]]
    diag_kappa[v] <- met[["kappa"]]
  }

  # Resolve amplitude policy defaults
  amp_method <- match.arg((amplitude$method %||% "global_ls"), c("global_ls","lss1","oasis_voxel"))
  ridge_gls  <- amplitude$ridge %||% list(mode = "fractional", lambda = 0.02)
  ridge_frac <- amplitude$ridge_frac %||% list(x = 0.02, b = 0.02)
  return_se  <- isTRUE(amplitude$return_se %||% FALSE)

  method_used <- rep(amp_method, V)
  if (!is.null(amplitude$cond_gate) && amp_method == "global_ls") {
    gate <- amplitude$cond_gate
    metric <- gate$metric %||% "rho"
    thr    <- as.numeric(gate$thr %||% 0.999)
    fallback <- match.arg(gate$fallback %||% "lss1", c("lss1","oasis_voxel"))
    trig <- if (metric == "kappa") (diag_kappa > (gate$kappa_thr %||% 1e3)) else (diag_rho > thr)
    method_used[which(trig)] <- fallback
  }

  amps <- matrix(NA_real_, ntrials, V)
  se_gls <- se_lss1 <- NULL
  if (any(method_used == "global_ls")) {
    sel <- method_used == "global_ls"
    # Adaptive ridge per voxel if requested
    ridge_use <- ridge_gls
    if (isTRUE(amplitude$adaptive$enable %||% FALSE)) {
      lam_vec <- .sbhm_adaptive_ridge_gls(diag_kappa[sel],
                                          base = amplitude$adaptive$base %||% 0.02,
                                          k0   = amplitude$adaptive$k0   %||% 1000,
                                          max_lam = amplitude$adaptive$max %||% 0.08)
      ridge_use <- lam_vec
    }
    res_gls <- sbhm_amplitude_ls(Y, sbhm, design_spec, alpha_hat[, sel, drop = FALSE],
                                 Nuisance = Nuisance, ridge = ridge_use, prewhiten = prewhiten,
                                 return_se = return_se)
    if (is.list(res_gls)) { amps[, sel] <- res_gls$beta; se_gls <- res_gls$se } else { amps[, sel] <- res_gls }
  }
  if (any(method_used == "lss1")) {
    sel <- method_used == "lss1"
    res_lss1 <- sbhm_amplitude_lss1(Y, sbhm, design_spec, alpha_hat[, sel, drop = FALSE],
                                    Nuisance = Nuisance, ridge_frac = ridge_frac,
                                    prewhiten = prewhiten, return_se = return_se)
    if (is.list(res_lss1)) { amps[, sel] <- res_lss1$beta; se_lss1 <- res_lss1$se } else { amps[, sel] <- res_lss1 }
  }
  oasis_se <- NULL
  if (any(method_used == "oasis_voxel")) {
    sel <- method_used == "oasis_voxel"
    res_o <- sbhm_amplitude_oasis_k1(Y, sbhm, design_spec, alpha_hat[, sel, drop = FALSE],
                                     Nuisance = Nuisance, ridge_frac = ridge_frac,
                                     prewhiten = prewhiten, return_se = return_se)
    amps[, sel] <- res_o$beta
    if (return_se) {
      # Build full SE matrix (NA for non-oasis voxels)
      se_full <- matrix(NA_real_, ntrials, V)
      se_full[, sel] <- res_o$se
      oasis_se <- se_full
    }
  }

  out <- list(
    matched_idx = m$idx,
    margin      = m$margin,
    alpha_coords= alpha_hat,
    diag        = list(r = r, ntrials = ntrials, times = sbhm$tgrid,
                       method_used = method_used, rho_max = diag_rho, kappa = diag_kappa,
                       ridge_gls = if (is.list(ridge_gls)) ridge_gls else as.list(ridge_gls),
                       ridge_frac = ridge_frac,
                       se = if (return_se) {
                         se_all <- matrix(NA_real_, ntrials, V)
                         if (!is.null(se_gls))   se_all[, method_used == "global_ls"]   <- se_gls[,  method_used == "global_ls",   drop = FALSE]
                         if (!is.null(se_lss1))  se_all[, method_used == "lss1"]        <- se_lss1[, method_used == "lss1",        drop = FALSE]
                         if (!is.null(oasis_se)) se_all[, method_used == "oasis_voxel"] <- oasis_se[, method_used == "oasis_voxel", drop = FALSE]
                         se_all
                       } else NULL)
  )
  if ((match$topK %||% 1) > 1) {
    out$topK_idx <- m$topK_idx
    out$weights  <- m$weights
    out$alpha_mode <- if (isTRUE(match$soft_blend %||% FALSE)) "soft" else "hard"
  }
  if (return %in% c("amplitude", "both")) out$amplitude <- amps
  if (return %in% c("coefficients", "both")) out$coeffs_r <- beta_rt

  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

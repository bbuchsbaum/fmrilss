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

  # 4) Projection to scalar amplitudes ---------------------------------------
  amps <- sbhm_project(beta_rt, alpha_hat)

  out <- list(
    matched_idx = m$idx,
    margin      = m$margin,
    alpha_coords= alpha_hat,
    diag        = list(r = r, ntrials = ntrials, times = sbhm$tgrid)
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

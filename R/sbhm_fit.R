#' SBHM Prepass: Aggregate Fit in Shared Basis
#'
#' Compute per-voxel coefficients in the shared SBHM basis by fitting a single
#' aggregate GLM with one regressor per basis column (trials summed), optionally
#' residualizing by nuisances and prewhitening. This produces `beta_bar` (r×V)
#' that you can feed to `sbhm_match()`.
#'
#' @param Y Numeric matrix T×V of fMRI time series.
#' @param sbhm SBHM object from `sbhm_build()` (must contain B, S, A, tgrid, span).
#' @param design_spec List describing events (same shape as `oasis$design_spec`).
#'   Must contain `sframe` and `cond` with `onsets` (and optional `duration`,
#'   `amplitude`, `span`). `cond$hrf` is ignored and replaced with `sbhm_hrf`.
#'   Optional `others` (list of other conditions) will be aggregated as nuisances.
#' @param Nuisance Optional T×P nuisance matrix (motion, drift, etc.).
#' @param prewhiten Optional fmriAR prewhitening options (see `?lss`). If provided,
#'   Y and design are prewhitened together.
#' @param ridge Optional list for targeted ridge shrinkage in the prepass solve:
#'   - `mode`: "fractional" (default) or "absolute". Fractional scales by mean(diag(G)).
#'   - `lambda`: nonnegative scalar (default 0.01 in fractional mode).
#'   - `alpha_ref`: r-vector to shrink towards (default `sbhm$ref$alpha_ref`).
#' @param data_fac Optional list for external factorization: `scores` (T×q),
#'   `loadings` (q×V). If provided, computes X'Y via (X'Scores) × Loadings. In this
#'   PR2 version, prewhitening is not applied when `data_fac` is used.
#'
#' @return List with:
#'   - `beta_bar` r×V aggregate coefficients
#'   - `A_agg`   T×r aggregated per-basis design (after any residualization/whitening)
#'   - `G`       r×r crossprod of A_agg
#'   - `diag`    list with K=r, ntrials, times, used_prewhiten
#'
#' @examples
#' \dontrun{
#'   library(fmrihrf)
#'   set.seed(1)
#'   Tlen <- 120; V <- 5; r <- 4
#'   sframe <- sampling_frame(blocklens = Tlen, TR = 1)
#'   H <- cbind(exp(-seq(0, 30, length.out = Tlen)/4),
#'              exp(-seq(0, 30, length.out = Tlen)/6))
#'   sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
#'   onsets <- seq(5, 95, by = 10)
#'   design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
#'   hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
#'   rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
#'   X <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
#'   betas_true <- matrix(rnorm(r), r)
#'   Y <- matrix(rnorm(Tlen*V, sd = 0.5), Tlen, V)
#'   Y[,1] <- Y[,1] + X %*% betas_true
#'   pre <- sbhm_prepass(Y, sbhm, design_spec)
#'   str(pre)
#' }
#'
#'
#' @details
#' Notes:
#' - Aggregated per-basis regressors can be highly collinear, making G = A' A
#'   ill-conditioned. A small ridge is recommended for stability. The default
#'   uses fractional mode with `lambda = 0.01` (scaled by mean(diag(G))).
#' - When `data_fac` is provided (factorized data path), prewhitening is skipped
#'   in this version; both dense and factorized paths perform nuisance
#'   residualization consistently.
#'
#' @export
sbhm_prepass <- function(Y, sbhm, design_spec,
                         Nuisance = NULL,
                         prewhiten = NULL,
                         ridge = list(mode = "fractional", lambda = 0.01, alpha_ref = NULL),
                         data_fac = NULL) {
  stopifnot(is.matrix(Y), is.list(sbhm), is.list(design_spec))

  # Validate required SBHM fields early for clearer errors
  required_fields <- c("B", "S", "A", "tgrid", "span", "ref")
  if (!all(required_fields %in% names(sbhm))) {
    missing <- setdiff(required_fields, names(sbhm))
    stop("sbhm object missing required fields: ", paste(missing, collapse = ", "))
  }
  r <- ncol(sbhm$B)

  # 1) Build aggregated per-basis design (A: T×r) directly using fmrihrf
  #    Avoids constructing full trial-wise matrix X_trials (memory/time saver)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  spec <- design_spec
  spec$cond <- spec$cond %||% list()
  spec$cond$hrf <- hrf_B

  onsets <- spec$cond$onsets
  if (is.null(onsets)) stop("design_spec$cond$onsets must be provided")
  # Aggregated regressor set with K=r columns (one per basis), trials summed
  rr_agg <- fmrihrf::regressor(
    onsets   = onsets,
    hrf      = hrf_B,
    duration = spec$cond$duration %||% 0,
    amplitude= spec$cond$amplitude %||% 1,
    span     = spec$cond$span %||% 40,
    summate  = TRUE
  )
  A <- fmrihrf::evaluate(rr_agg, grid = sbhm$tgrid,
                         precision = spec$precision %||% 0.1,
                         method = spec$method %||% "conv")
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  if (!is.matrix(A)) A <- cbind(A)

  # Aggregates for "other" conditions (one col each) -> nuisances
  X_other <- NULL
  if (length(spec$others)) {
    X_other <- do.call(cbind, lapply(spec$others, function(oc) {
      rr <- fmrihrf::regressor(onsets   = oc$onsets,
                               hrf      = oc$hrf %||% hrf_B,
                               duration = oc$duration %||% 0,
                               amplitude= oc$amplitude %||% 1,
                               span     = oc$span %||% 40,
                               summate  = TRUE)
      x_eval <- fmrihrf::evaluate(rr, sbhm$tgrid,
                                  precision = spec$precision %||% 0.1,
                                  method = spec$method %||% "conv")
      if (inherits(x_eval, "Matrix")) x_eval <- as.matrix(x_eval)
      # For multi-basis HRFs, aggregate columns to a single condition column
      if (is.matrix(x_eval)) rowSums(x_eval) else as.numeric(x_eval)
    }))
    if (is.null(dim(X_other))) X_other <- matrix(X_other, ncol = 1)
  }
  K <- r
  ntrials <- as.integer(length(onsets))

  # 3) Nuisance combination (intercept + provided + others)
  Zint <- matrix(1, nrow(Y), 1)
  N_nuis <- cbind(Zint,
                  if (!is.null(Nuisance)) Nuisance,
                  if (!is.null(X_other))  X_other)

  used_prewhiten <- FALSE

  # 4) Optional prewhitening (dense Y only in PR2)
  if (!is.null(prewhiten) && is.null(data_fac)) {
    whitened <- .prewhiten_data(Y, A, NULL, N_nuis, prewhiten)
    Yw <- whitened$Y_whitened
    Aw <- whitened$X_whitened
    Nw <- whitened$Nuisance_whitened
    used_prewhiten <- TRUE
  } else {
    if (!is.null(prewhiten) && !is.null(data_fac)) {
      warning("prewhitening requested but data_fac provided: skipping prewhitening for factorized path")
    }
    Yw <- Y
    Aw <- A
    Nw <- N_nuis
  }

  # 5) Residualize by nuisances (FRWL) for dense path; for factorized path, apply QR to A and Scores
  if (is.null(data_fac)) {
    if (!is.null(Nw) && ncol(Nw) > 0) {
      qrN <- qr(Nw)
      Y_res <- qr.resid(qrN, Yw)
      A_res <- qr.resid(qrN, Aw)
    } else {
      Y_res <- Yw
      A_res <- Aw
    }
    G <- crossprod(A_res)          # r×r
    F <- crossprod(A_res, Y_res)   # r×V
  } else {
    # Factorized path (no prewhitening in PR2): residualize A and Scores by N_nuis
    Scores <- as.matrix(data_fac$scores)
    Load   <- as.matrix(data_fac$loadings)
    if (!is.null(N_nuis) && ncol(N_nuis) > 0) {
      qrN <- qr(N_nuis)
      A_res <- qr.resid(qrN, A)
      S_res <- qr.resid(qrN, Scores)
    } else {
      A_res <- A
      S_res <- Scores
    }
    G <- crossprod(A_res)
    XS <- crossprod(A_res, S_res)  # r×q
    F <- XS %*% Load               # r×V
  }

  # 6) Targeted ridge (optional) towards alpha_ref
  ridge$mode   <- ridge$mode %||% "fractional"
  ridge$lambda <- ridge$lambda %||% 0.01
  alpha_ref <- ridge$alpha_ref %||% sbhm$ref$alpha_ref
  lam <- as.numeric(ridge$lambda)
  if (tolower(ridge$mode) == "fractional") {
    lam <- lam * mean(diag(G))
  }
  if (lam > 0) {
    beta_bar <- qr.solve(G + diag(lam, r), F + alpha_ref %*% matrix(lam, nrow = 1, ncol = ncol(F)))
  } else {
    beta_bar <- qr.solve(G, F)
  }

  list(
    beta_bar = beta_bar,
    A_agg    = A_res,
    G        = G,
    diag     = list(K = K, ntrials = ntrials, times = sbhm$tgrid, used_prewhiten = used_prewhiten)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

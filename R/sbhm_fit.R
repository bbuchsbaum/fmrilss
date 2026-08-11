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
#'   - `alpha_ref`: r-vector to shrink towards (default zero vector).
#' @param data_fac Optional list for external factorization: `scores` (T×q),
#'   `loadings` (q×V). If provided, computes X'Y via (X'Scores) × Loadings. In this
#'   The original full `Y` is still required by `lss_sbhm()` for later stages.
#'   Active estimated prewhitening cannot be combined with `data_fac`.
#'
#' @return List with:
#'   - `beta_bar` r×V aggregate coefficients
#'   - `A_agg`   T×r aggregated per-basis design (after any residualization/whitening)
#'   - `G`       r×r crossprod of A_agg
#'   - `diag`    list with K=r, ntrials, times, used_prewhiten
#'
#' @examples
#' \donttest{
#'   library(fmrihrf)
#'   set.seed(1)
#'   Tlen <- 120; V <- 5; r <- 4
#'   sframe <- sampling_frame(blocklens = Tlen, TR = 1)
#'   H <- cbind(exp(-seq(0, 30, length.out = Tlen)/4),
#'              exp(-seq(0, 30, length.out = Tlen)/6))
#'   sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
#'   r <- ncol(sbhm$B)
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
#' - When `data_fac` is provided, scores must be T×q and loadings q×V. Both
#'   dense and factorized paths perform nuisance residualization consistently;
#'   active prewhitening fails explicitly because the supplied factorization is
#'   not a factorization of the estimated whitened data.
#' @noRd
#'
# Build a run-safe, trial-major SBHM design once and share it across the
# prepass, matching diagnostics, and scalar-amplitude stages.
.sbhm_events_from_condition <- function(condition, sframe, label = "target") {
  if (!is.list(condition) || is.null(condition$onsets)) {
    stop(sprintf("%s condition must provide onsets", label), call. = FALSE)
  }
  onsets <- suppressWarnings(as.numeric(condition$onsets))
  n <- length(onsets)
  if (n < 1L || any(!is.finite(onsets))) {
    stop(sprintf("%s onsets must contain finite values", label), call. = FALSE)
  }
  recycle_event_field <- function(x, default, field) {
    x <- x %||% default
    if (!is.numeric(x) || !(length(x) %in% c(1L, n)) || any(!is.finite(x))) {
      stop(sprintf("%s %s must have length 1 or %d and be finite", label, field, n),
           call. = FALSE)
    }
    rep_len(as.numeric(x), n)
  }
  duration <- recycle_event_field(condition$duration, 0, "duration")
  if (any(duration < 0)) {
    stop(sprintf("%s duration values must be nonnegative", label), call. = FALSE)
  }
  amplitude <- recycle_event_field(condition$amplitude, 1, "amplitude")
  if (identical(label, "target") && any(amplitude == 0)) {
    stop("target amplitude values must be non-zero", call. = FALSE)
  }
  events <- data.frame(
    onset = onsets,
    duration = duration,
    condition = rep.int(label, n),
    amplitude = amplitude,
    stringsAsFactors = FALSE
  )
  if (!is.null(condition$run)) events$run <- condition$run
  .voxhrf_validate_events(events, sframe)
  events
}

.sbhm_build_design <- function(sbhm, design_spec) {
  if (!is.list(design_spec) || !inherits(design_spec$sframe, "sampling_frame")) {
    stop("design_spec must contain an fmrihrf sampling_frame", call. = FALSE)
  }
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  target_events <- .sbhm_events_from_condition(
    design_spec$cond, design_spec$sframe, "target"
  )
  built <- .voxhrf_trial_basis(
    target_events, hrf_B, design_spec$sframe,
    precision = design_spec$precision %||% 0.1,
    method = design_spec$method %||% "conv",
    span = design_spec$cond$span %||% sbhm$span
  )
  r <- ncol(sbhm$B)
  if (built$K != r) stop("SBHM design basis dimension is inconsistent", call. = FALSE)
  ntrials <- nrow(target_events)
  X_trials <- built$X
  regs <- lapply(seq_len(ntrials), function(j) {
    cols <- ((j - 1L) * r + 1L):(j * r)
    X_trials[, cols, drop = FALSE]
  })
  A_agg <- do.call(cbind, lapply(seq_len(r), function(k) {
    rowSums(X_trials[, seq.int(k, ncol(X_trials), by = r), drop = FALSE])
  }))
  colnames(A_agg) <- paste0("basis_", seq_len(r))

  X_other <- NULL
  if (length(design_spec$others)) {
    other_blocks <- lapply(seq_along(design_spec$others), function(i) {
      oc <- design_spec$others[[i]]
      events <- .sbhm_events_from_condition(
        oc, design_spec$sframe, paste0("other_", i)
      )
      other_hrf <- oc$hrf %||% hrf_B
      ob <- .voxhrf_trial_basis(
        events, other_hrf, design_spec$sframe,
        precision = design_spec$precision %||% 0.1,
        method = design_spec$method %||% "conv",
        span = oc$span %||% sbhm$span
      )
      do.call(cbind, lapply(seq_len(ob$K), function(k) {
        rowSums(ob$X[, seq.int(k, ncol(ob$X), by = ob$K), drop = FALSE])
      }))
    })
    X_other <- do.call(cbind, other_blocks)
    if (is.null(dim(X_other))) X_other <- matrix(X_other, ncol = 1L)
    colnames(X_other) <- make.unique(paste0("other_", seq_len(ncol(X_other))))
  }

  trial_names <- design_spec$cond$trial_names %||% paste0("trial_", seq_len(ntrials))
  trial_names <- as.character(trial_names)
  if (length(trial_names) != ntrials || anyNA(trial_names) ||
      any(!nzchar(trial_names)) || anyDuplicated(trial_names)) {
    stop("target trial names must be complete and unique", call. = FALSE)
  }
  output_names <- if (r == 1L) trial_names else {
    as.vector(t(outer(trial_names, paste0("basis_", seq_len(r)), paste, sep = ":")))
  }
  colnames(X_trials) <- output_names
  map <- data.frame(
    output_row = seq_along(output_names),
    column = output_names,
    trial = rep(seq_len(ntrials), each = r),
    trial_name = rep(trial_names, each = r),
    basis = rep(seq_len(r), times = ntrials),
    run = rep(target_events$run %||% rep.int(1L, ntrials), each = r),
    name = output_names,
    output_name = output_names,
    stringsAsFactors = FALSE
  )
  blocks <- fmrihrf::blockids(design_spec$sframe)
  run_levels <- sort(unique(blocks))
  intercepts <- vapply(
    run_levels, function(run) as.numeric(blocks == run), numeric(length(blocks))
  )
  if (is.null(dim(intercepts))) intercepts <- matrix(intercepts, ncol = 1L)
  colnames(intercepts) <- paste0("run_", run_levels, "_intercept")

  list(
    X_trials = X_trials,
    regs = regs,
    A_agg = A_agg,
    X_other = X_other,
    map = map,
    trial_names = trial_names,
    events = target_events,
    intercepts = intercepts,
    K = r,
    ntrials = ntrials
  )
}

# Prewhitening filters must restart at every sampling-frame run boundary.
# Supply the canonical segmentation when callers omit it, and reject an
# explicit segmentation that disagrees with the event/design frame.
.sbhm_run_safe_prewhiten <- function(prewhiten, sframe) {
  if (is.null(prewhiten)) {
    return(prewhiten)
  }
  if (!is.list(prewhiten)) {
    stop("prewhiten must be a list of prewhitening options", call. = FALSE)
  }
  trusted_plan <- inherits(prewhiten, "fmrilss_internal_prewhiten")
  prewhiten <- .resolve_prewhiten_options(prewhiten, internal = trusted_plan)
  if (identical(prewhiten$method, "none")) return(prewhiten)
  blocks <- fmrihrf::blockids(sframe)
  if (is.null(prewhiten$runs)) {
    prewhiten$runs <- blocks
    return(prewhiten)
  }
  runs <- prewhiten$runs
  if (length(runs) != length(blocks) || anyNA(runs)) {
    stop("prewhiten$runs must contain one complete run label per scan", call. = FALSE)
  }
  run_codes <- match(runs, unique(runs))
  block_codes <- match(blocks, unique(blocks))
  if (!identical(as.integer(run_codes), as.integer(block_codes))) {
    stop("prewhiten$runs must match the sampling-frame run boundaries", call. = FALSE)
  }
  prewhiten
}

#' SBHM Prepass: Aggregate Fit in a Shared Basis
#'
#' Fit the trial-aggregated SBHM basis to every voxel after projecting
#' run-specific intercepts, supplied nuisance columns, and complete modeled
#' other-condition spans. The same run-safe design builder is used by the full
#' SBHM pipeline.
#'
#' @param Y Finite numeric T by V response matrix.
#' @param sbhm Object returned by [sbhm_build()].
#' @param design_spec Event specification with `sframe` and `cond`; multi-run
#'   inputs use run-relative onsets and an explicit `cond$run` vector.
#' @param Nuisance Optional T by P nuisance matrix.
#' @param prewhiten Optional prewhitening options passed to the package
#'   prewhitening layer. Run labels are inferred from `design_spec$sframe` when
#'   omitted; supplied labels must match those sampling-frame boundaries.
#' @param ridge List with `mode`, nonnegative `lambda`, and optional
#'   rank-length `alpha_ref`.
#' @param data_fac Optional exact or approximate factorization with `scores`
#'   T by q and `loadings` q by V. Active prewhitening is unsupported with this
#'   shortcut, and the original full `Y` remains required by [lss_sbhm()]. If
#'   factor axes are named, scores and loadings must provide the same complete
#'   unique names and are aligned before multiplication. When `Y` is named,
#'   loadings must provide the same complete unique voxel names and are reordered
#'   to `Y`.
#'
#' @return A list containing named `beta_bar` (rank by voxel), the residualized
#'   aggregate design `A_agg`, its Gram matrix `G`, and design diagnostics and
#'   identity maps.
#' @export
sbhm_prepass <- function(Y, sbhm, design_spec,
                         Nuisance = NULL,
                         prewhiten = NULL,
                         ridge = list(mode = "fractional", lambda = 0.01, alpha_ref = NULL),
                         data_fac = NULL) {
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 1L || ncol(Y) < 1L ||
      any(!is.finite(Y))) {
    stop("Y must be a non-empty finite numeric matrix", call. = FALSE)
  }
  if (!is.null(colnames(Y)) &&
      (anyNA(colnames(Y)) || any(!nzchar(colnames(Y))) || anyDuplicated(colnames(Y)))) {
    stop("Y voxel names must be complete and unique", call. = FALSE)
  }
  if (!is.list(sbhm) || !is.list(design_spec)) {
    stop("sbhm and design_spec must be lists", call. = FALSE)
  }

  # Validate required SBHM fields early for clearer errors
  required_fields <- c("B", "S", "A", "tgrid", "span", "ref")
  if (!all(required_fields %in% names(sbhm))) {
    missing <- setdiff(required_fields, names(sbhm))
    stop("sbhm object missing required fields: ", paste(missing, collapse = ", "))
  }
  r <- ncol(sbhm$B)

  built <- .sbhm_build_design(sbhm, design_spec)
  prewhiten <- .sbhm_run_safe_prewhiten(prewhiten, design_spec$sframe)
  A <- built$A_agg
  X_other <- built$X_other
  K <- r
  ntrials <- built$ntrials
  if (nrow(A) != nrow(Y)) stop("SBHM design rows must match Y", call. = FALSE)

  # 3) Nuisance combination (intercept + provided + others)
  N_nuis <- cbind(built$intercepts,
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
    if (!is.null(prewhiten) &&
        !identical(prewhiten$method %||% "none", "none")) {
      stop("factorized SBHM prepass does not support active prewhitening", call. = FALSE)
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
    if (!is.numeric(Scores) || !is.numeric(Load) || any(!is.finite(Scores)) ||
        any(!is.finite(Load)) || nrow(Scores) != nrow(Y) ||
        nrow(Load) != ncol(Scores) || ncol(Load) != ncol(Y)) {
      stop(
        "data_fac requires finite scores (T x q) and loadings (q x V) matching Y",
        call. = FALSE
      )
    }
    score_names <- colnames(Scores)
    factor_names <- rownames(Load)
    if (xor(is.null(score_names), is.null(factor_names))) {
      stop("factor names must be supplied on both scores and loadings or neither",
           call. = FALSE)
    }
    if (!is.null(score_names)) {
      if (anyNA(score_names) || anyNA(factor_names) || any(!nzchar(score_names)) ||
          any(!nzchar(factor_names)) || anyDuplicated(score_names) ||
          anyDuplicated(factor_names) || !setequal(score_names, factor_names)) {
        stop("scores and loadings factor names must be complete, unique, and identical",
             call. = FALSE)
      }
      Load <- Load[match(score_names, factor_names), , drop = FALSE]
    }
    y_names <- colnames(Y)
    load_names <- colnames(Load)
    if (xor(is.null(y_names), is.null(load_names))) {
      stop("named Y requires loadings with the same complete voxel names",
           call. = FALSE)
    }
    if (!is.null(y_names)) {
      if (anyNA(load_names) || any(!nzchar(load_names)) ||
          anyDuplicated(load_names) || !setequal(y_names, load_names)) {
        stop("loadings voxel names must be complete, unique, and match Y",
             call. = FALSE)
      }
      Load <- Load[, match(y_names, load_names), drop = FALSE]
    }
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
  ridge <- ridge %||% list(mode = "fractional", lambda = 0.01, alpha_ref = NULL)
  if (!is.list(ridge)) stop("ridge must be a list", call. = FALSE)
  ridge$mode   <- match.arg(ridge$mode %||% "fractional", c("fractional", "absolute"))
  ridge$lambda <- ridge$lambda %||% 0.01
  alpha_ref <- ridge$alpha_ref %||% rep(0, r)
  alpha_ref <- as.numeric(alpha_ref)
  if (length(alpha_ref) == 1L) alpha_ref <- rep(alpha_ref, r)
  if (length(alpha_ref) != r) stop("ridge$alpha_ref must have length r")
  lam <- as.numeric(ridge$lambda)
  if (length(lam) != 1L || !is.finite(lam) || lam < 0) {
    stop("ridge$lambda must be one nonnegative finite value", call. = FALSE)
  }
  if (tolower(ridge$mode) == "fractional") {
    lam <- lam * mean(diag(G))
  }
  if (lam > 0) {
    beta_bar <- qr.solve(G + diag(lam, r), F + alpha_ref %*% matrix(lam, nrow = 1, ncol = ncol(F)))
  } else {
    beta_bar <- qr.solve(G, F)
  }
  basis_names <- rownames(sbhm$A) %||% paste0("basis_", seq_len(r))
  voxel_names <- colnames(Y) %||% paste0("voxel_", seq_len(ncol(Y)))
  dimnames(beta_bar) <- list(basis_names, voxel_names)
  dimnames(G) <- list(basis_names, basis_names)

  list(
    beta_bar = beta_bar,
    A_agg    = A_res,
    G        = G,
    diag     = list(
      K = K, ntrials = ntrials, times = sbhm$tgrid,
      used_prewhiten = used_prewhiten,
      trial_basis_map = built$map,
      trial_names = built$trial_names
    )
  )
}

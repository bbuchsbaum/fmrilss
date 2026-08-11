#' VoxelHRF object
#'
#' Simple list-based S3 class returned by \code{estimate_voxel_hrf} containing
#' voxel-wise HRF basis coefficients and related metadata.
#' @name VoxelHRF
#' @section Stored fields:
#' `coefficients` contains one normalized HRF-shape column per voxel;
#' `amplitude_scale` records the removed positive-peak scale; `basis` stores the
#' HRF basis; `conditions` records observed labels while `condition_pooling`
#' states that all events estimate one pooled shape; `sframe` preserves physical scan
#' timing; and `normalization` plus `coefficient_units` state the coefficient
#' scale. Voxel names, when supplied on `Y`, identify coefficient columns and
#' amplitude-scale entries.
#' @return No value itself. This topic documents the structure returned by
#'   `estimate_voxel_hrf()`.
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' Y <- matrix(rnorm(100), 50, 2)
#' events <- data.frame(onset = c(5, 25), duration = 1, condition = "A")
#' basis <- fmrihrf::HRF_SPMG1
#' sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
#' est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
#' class(est)
#' }
NULL

#' LSSBeta object
#'
#' Simple list-based S3 class returned by \code{lss_with_hrf} containing
#' trial-wise beta estimates. `as.matrix()` materializes C++-backed estimates
#' as a dense matrix while preserving voxel names, sampling-frame,
#' normalization, and unit metadata.
#' @name LSSBeta
#' @section Stored fields:
#' `betas` is the file-backed trial-by-voxel matrix; `dimnames` preserves trial
#' and voxel identity; `sframe`, `normalization`, `units`, and
#' `event_amplitude` and `event_duration` define timing and coefficient
#' interpretation; and
#' `engine_requested`, `engine_used`, and `chunk_size` record execution. Use
#' `as.matrix()` for the documented dense extraction path; these metadata are
#' retained as matrix attributes.
#' @return No value itself. This topic documents the object returned by
#'   `lss_with_hrf(..., engine = "C++")`.
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' Y <- matrix(rnorm(100), 50, 2)
#' events <- data.frame(onset = c(5, 25), duration = 1, condition = "A")
#' basis <- fmrihrf::HRF_SPMG1
#' sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
#' est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
#' fit <- lss_with_hrf(Y, events, est, engine = "C++", verbose = FALSE)
#' class(fit)
#' }
NULL

#' Estimate Voxel-wise HRF Basis Coefficients
#'
#' Fits a common-amplitude GLM to estimate one pooled HRF shape for every
#' voxel. The HRF basis and response are residualized against the complete
#' fixed and nuisance span before fitting. Estimation fails when that span
#' makes the HRF basis unidentifiable or rank deficient.
#'
#' @param Y Numeric matrix of BOLD data (time x voxels).
#' @param events Data frame with \code{onset}, \code{duration} and
#'   \code{condition} columns. For a multi-run sampling frame, a \code{run}
#'   column is required and onsets are physical times relative to that run.
#' @param basis HRF object from the \code{fmrihrf} package.
#' @param nuisance_regs Optional finite numeric matrix of nuisance regressors.
#' @param sframe Explicit \code{fmrihrf} sampling frame defining scan times and
#'   TR. Required; onset and duration values use its physical-time units.
#' @param fixed_regs Optional finite numeric matrix of fixed/common regressors. An
#'   intercept is added when it is not already in their span.
#'
#' @return A \link{VoxelHRF} object containing at least:
#'   \item{coefficients}{Matrix of positive-peak-normalized HRF shape weights
#'     with one row per basis function and one column per voxel.}
#'   \item{amplitude_scale}{The signed scale removed from each raw pooled-fit
#'     coefficient column.}
#'   \item{basis}{The HRF basis object used.}
#'   \item{conditions}{Observed event labels. Labels are metadata: all events
#'     are pooled into one shape per voxel.}
#'   \item{condition_pooling}{The literal string "all-events".}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(rnorm(100), 50, 2)
#' events <- data.frame(onset = c(5, 25), duration = 1,
#'                      condition = "A")
#' basis <- fmrihrf::HRF_SPMG1
#' sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
#' times <- fmrihrf::samples(sframe, global = TRUE)
#' rset <- fmrihrf::regressor_set(onsets = events$onset,
#'                                fac = factor(rep("all events", nrow(events))),
#'                                hrf = basis, duration = events$duration,
#'                                span = 30)
#' X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
#' coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
#' Y <- X %*% coef + Y * 0.1
#' est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
#' str(est)
#' }
#' @export
estimate_voxel_hrf <- function(Y, events, basis, nuisance_regs = NULL,
                               sframe = NULL, fixed_regs = NULL) {
  # Input validation
  .voxhrf_validate_response(Y)

  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition")
  }

  if (!inherits(basis, "HRF")) {
    stop("basis must be an 'HRF' object from fmrihrf")
  }

  if (!is.null(nuisance_regs)) {
    if (!is.matrix(nuisance_regs) || !is.numeric(nuisance_regs)) {
      stop("nuisance_regs must be a numeric matrix")
    }
    if (nrow(nuisance_regs) != nrow(Y)) {
      stop("nuisance_regs must have same number of rows as Y")
    }
    if (any(!is.finite(nuisance_regs))) {
      stop("nuisance_regs must contain only finite values", call. = FALSE)
    }
  }

  if (!is.null(fixed_regs) &&
      (!is.matrix(fixed_regs) || !is.numeric(fixed_regs) || nrow(fixed_regs) != nrow(Y))) {
    stop("fixed_regs must be a numeric matrix with same number of rows as Y", call. = FALSE)
  }
  if (!is.null(fixed_regs) && any(!is.finite(fixed_regs))) {
    stop("fixed_regs must contain only finite values", call. = FALSE)
  }
  if (is.null(sframe) || !inherits(sframe, "sampling_frame")) {
    stop("sframe must be an explicit fmrihrf sampling_frame", call. = FALSE)
  }
  times <- fmrihrf::samples(sframe, global = TRUE)
  if (length(times) != nrow(Y)) {
    stop("sframe must contain exactly nrow(Y) scan times", call. = FALSE)
  }

  # Construct each trial within its own run, then aggregate trial columns by
  # basis. This preserves the pooled-shape estimand without convolving across
  # acquisition boundaries.
  trial_design <- .voxhrf_trial_basis(events, basis, sframe)
  X_basis <- do.call(cbind, lapply(trial_design$basis_convolved, rowSums))

  fixed_use <- .voxhrf_fixed_design(fixed_regs, nrow(Y), sframe)
  common_design <- cbind(fixed_use, nuisance_regs)
  common_span <- .voxhrf_orthonormal_span(common_design)
  residualize <- function(value) {
    if (!ncol(common_span)) return(value)
    value - common_span %*% crossprod(common_span, value)
  }
  X_residual <- residualize(X_basis)
  Y_residual <- residualize(Y)

  original_norm <- vapply(seq_len(ncol(X_basis)), function(column) {
    .voxhrf_stable_norm(X_basis[, column])
  }, numeric(1))
  residual_norm <- vapply(seq_len(ncol(X_residual)), function(column) {
    .voxhrf_stable_norm(X_residual[, column])
  }, numeric(1))
  if (any(original_norm == 0) ||
      any(residual_norm / original_norm <= 1e-8)) {
    stop(
      "voxel HRF basis is not identifiable after projecting the common design",
      call. = FALSE
    )
  }
  residual_unit <- sweep(X_residual, 2L, residual_norm, "/")
  singular_values <- svd(residual_unit, nu = 0L, nv = 0L)$d
  if (length(singular_values) < ncol(X_residual) ||
      min(singular_values) / max(singular_values) <= 1e-8) {
    stop(
      "voxel HRF basis is rank-deficient after projecting the common design",
      call. = FALSE
    )
  }

  coef_basis_unit <- qr.coef(qr(residual_unit), Y_residual)
  coef_basis_raw <- sweep(coef_basis_unit, 1L, residual_norm, "/")
  if (!is.matrix(coef_basis_raw)) {
    coef_basis_raw <- matrix(coef_basis_raw, nrow = ncol(X_basis))
  }
  normalized <- .normalize_voxel_hrf_coefficients(
    coef_basis_raw, basis,
    span = if (!is.null(attr(basis, "span"))) attr(basis, "span") else 30,
    precision = 0.05
  )
  colnames(normalized$coefficients) <- colnames(Y)
  names(normalized$amplitude_scale) <- colnames(Y)

  result <- list(
    coefficients = normalized$coefficients,
    amplitude_scale = normalized$amplitude_scale,
    basis = basis,
    conditions = unique(as.character(events$condition)),
    sframe = sframe,
    condition_pooling = "all-events",
    normalization = "positive-peak",
    coefficient_units = "unit-peak HRF shape weights"
  )
  class(result) <- "VoxelHRF"
  result
}

.voxhrf_stable_norm <- function(value) {
  scale <- max(abs(value))
  if (!is.finite(scale) || scale == 0) return(0)
  scale * sqrt(sum((value / scale)^2))
}

.voxhrf_orthonormal_span <- function(design, tolerance = 1e-10) {
  design <- as.matrix(design)
  if (!ncol(design)) return(matrix(numeric(), nrow(design), 0L))
  normalized <- matrix(0, nrow(design), ncol(design))
  keep <- logical(ncol(design))
  for (column in seq_len(ncol(design))) {
    column_norm <- .voxhrf_stable_norm(design[, column])
    if (column_norm > 0) {
      normalized[, column] <- design[, column] / column_norm
      keep[column] <- TRUE
    }
  }
  normalized <- normalized[, keep, drop = FALSE]
  if (!ncol(normalized)) return(matrix(numeric(), nrow(design), 0L))
  decomposition <- qr(normalized, tol = tolerance, LAPACK = FALSE)
  if (!decomposition$rank) return(matrix(numeric(), nrow(design), 0L))
  qr.Q(decomposition, complete = FALSE)[
    , seq_len(decomposition$rank), drop = FALSE
  ]
}

.voxhrf_fixed_design <- function(fixed_regs, n_time, sframe = NULL) {
  fixed <- if (is.null(fixed_regs)) matrix(numeric(), n_time, 0L) else as.matrix(fixed_regs)
  block <- if (is.null(sframe)) rep.int(1L, n_time) else fmrihrf::blockids(sframe)
  baseline <- if (length(unique(block)) == 1L) {
    matrix(1, nrow = n_time, ncol = 1L)
  } else {
    stats::model.matrix(~ 0 + factor(block))
  }
  colnames(baseline) <- if (ncol(baseline) == 1L) "Intercept" else {
    paste0("run_", seq_len(ncol(baseline)))
  }
  for (j in seq_len(ncol(baseline))) {
    candidate <- baseline[, j]
    fixed_span <- .voxhrf_orthonormal_span(fixed)
    residual <- if (ncol(fixed_span)) {
      candidate - fixed_span %*% crossprod(fixed_span, candidate)
    } else {
      candidate
    }
    in_span <- .voxhrf_stable_norm(residual) < 1e-8
    if (!in_span) fixed <- cbind(fixed, baseline[, j, drop = FALSE])
  }
  if (is.null(colnames(fixed))) colnames(fixed) <- paste0("fixed_", seq_len(ncol(fixed)))
  unnamed <- !nzchar(colnames(fixed))
  colnames(fixed)[unnamed] <- paste0("fixed_", which(unnamed))
  fixed
}

.normalize_voxel_hrf_coefficients <- function(coefficients, basis, span, precision = 0.05) {
  grid <- seq(0, span, by = precision)
  impulse <- fmrihrf::regressor(onsets = 0, hrf = basis, duration = 0, span = span)
  H <- fmrihrf::evaluate(impulse, grid, precision = precision, method = "conv")
  if (inherits(H, "Matrix")) H <- as.matrix(H)
  if (!is.matrix(H)) H <- matrix(H, ncol = 1L)
  if (ncol(H) != nrow(coefficients)) {
    stop("basis dimension does not match coefficient rows", call. = FALSE)
  }
  H <- sweep(H, 2L, H[1L, ], "-")
  reference <- H[, 1L]
  waveforms <- H %*% coefficients
  orientation <- sign(as.numeric(crossprod(reference, waveforms)))
  orientation[!is.finite(orientation) | orientation == 0] <- 1
  oriented <- sweep(waveforms, 2L, orientation, "*")
  peak <- apply(oriented, 2L, max)
  if (any(!is.finite(peak) | peak <= sqrt(.Machine$double.eps))) {
    stop("voxel HRF has no identifiable positive peak", call. = FALSE)
  }
  list(
    coefficients = sweep(sweep(coefficients, 2L, orientation, "*"), 2L, peak, "/"),
    amplitude_scale = orientation * peak
  )
}

#' Perform LSS using Voxel-wise HRFs
#'
#' Computes trial-wise beta estimates using voxel-specific HRFs.
#'
#' @param Y Numeric matrix of BOLD data (time x voxels).
#' @param events Data frame with \code{onset}, \code{duration} and
#'   \code{condition} columns. For a multi-run sampling frame, a \code{run}
#'   column is required and onsets are physical times relative to that run.
#' @param hrf_estimates A \link{VoxelHRF} object returned by \code{estimate_voxel_hrf}.
#' @param nuisance_regs Optional finite numeric matrix of nuisance regressors.
#' @param engine Computational engine: "R" for the pure-R implementation
#'   (default), or "C++" to request the compiled backend. The returned metadata
#'   records the backend actually used after fallback.
#' @param chunk_size Number of voxels to process per batch (C++ engine only).
#' @param verbose Logical; emit progress messages.
#' @param backing_dir Directory for bigmemory backing files. If NULL, a
#'   temporary directory is used (C++ engine only).
#' @param sframe Explicit sampling frame. Defaults to the frame stored in
#'   \code{hrf_estimates}; no unit-TR fallback is used.
#' @param fixed_regs Optional finite numeric matrix of fixed/common regressors. An
#'   intercept is added when absent.
#'
#' @return An object of class \link{LSSBeta} for the C++ request, or a numeric
#'   matrix (n_trials x n_vox) for the R engine. Matrix outputs carry
#'   sampling-frame, normalization, coefficient-unit, event-amplitude,
#'   event-duration, requested-engine, realized-engine, and chunk-size metadata.
#'   With a positive-peak-normalized shape, a zero-duration unit-amplitude event
#'   has a peak-response-amplitude coefficient. Otherwise the result is a
#'   coefficient on the supplied duration- and amplitude-coded event design.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(rnorm(100), 50, 2)
#' events <- data.frame(onset = c(5, 25), duration = 1,
#'                      condition = "A")
#' basis <- fmrihrf::HRF_SPMG1
#' sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
#' times <- fmrihrf::samples(sframe, global = TRUE)
#' rset <- fmrihrf::regressor_set(onsets = events$onset,
#'                                fac = factor(rep("all events", nrow(events))),
#'                                hrf = basis, duration = events$duration,
#'                                span = 30)
#' X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
#' coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
#' Y <- X %*% coef + Y * 0.1
#' est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
#' betas <- lss_with_hrf(Y, events, est, verbose = FALSE, engine = "R")
#' dim(betas)
#' }
#' @importFrom bigmemory filebacked.big.matrix
#' @export
lss_with_hrf <- function(Y, events, hrf_estimates, nuisance_regs = NULL,
                         engine = "R", chunk_size = 5000, verbose = TRUE,
                         backing_dir = NULL, sframe = NULL, fixed_regs = NULL) {
  engine <- match.arg(engine, c("R", "C++"))

  .voxhrf_validate_response(Y)

  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition")
  }

  coefficients <- .voxhrf_align_coefficients(Y, hrf_estimates)

  if (!is.null(nuisance_regs)) {
    if (!is.matrix(nuisance_regs) || !is.numeric(nuisance_regs) ||
        nrow(nuisance_regs) != nrow(Y)) {
      stop("nuisance_regs must be a numeric matrix with same number of rows as Y")
    }
    if (any(!is.finite(nuisance_regs))) {
      stop("nuisance_regs must contain only finite values", call. = FALSE)
    }
  }

  if (!is.null(fixed_regs) &&
      (!is.matrix(fixed_regs) || !is.numeric(fixed_regs) || nrow(fixed_regs) != nrow(Y))) {
    stop("fixed_regs must be a numeric matrix with same number of rows as Y", call. = FALSE)
  }
  if (!is.null(fixed_regs) && any(!is.finite(fixed_regs))) {
    stop("fixed_regs must contain only finite values", call. = FALSE)
  }

  chunk_size <- .as_positive_integer(chunk_size, "chunk_size")

  sframe <- sframe %||% hrf_estimates$sframe
  if (is.null(sframe) || !inherits(sframe, "sampling_frame")) {
    stop("sframe must be supplied explicitly or stored in hrf_estimates", call. = FALSE)
  }
  times <- fmrihrf::samples(sframe, global = TRUE)
  if (length(times) != nrow(Y)) {
    stop("sframe must contain exactly nrow(Y) scan times", call. = FALSE)
  }

  built <- .voxhrf_trial_basis(events, hrf_estimates$basis, sframe)
  event_amplitude <- if ("amplitude" %in% names(events)) {
    as.numeric(events$amplitude)
  } else {
    rep.int(1, nrow(events))
  }
  if (any(event_amplitude == 0)) {
    stop(
      "events amplitude values must be nonzero for trial-wise voxel-HRF LSS",
      call. = FALSE
    )
  }
  event_duration <- as.numeric(events$duration)
  fixed_use <- .voxhrf_fixed_design(fixed_regs, nrow(Y), sframe)
  normalization <- hrf_estimates$normalization %||% "unspecified"
  units <- if (identical(normalization, "positive-peak")) {
    if (all(event_amplitude == 1) && all(event_duration == 0)) {
      "peak-response amplitude"
    } else {
      "coefficient on supplied event design with unit-peak HRF shape"
    }
  } else {
    "coefficient on supplied event design relative to supplied HRF shape"
  }

  solve_voxels <- function(voxel_index, method, show_progress = FALSE) {
    lss_with_hrf_pure_r(
      Y = Y[, voxel_index, drop = FALSE],
      onset_idx = seq_len(nrow(events)),
      durations = rep.int(0L, nrow(events)),
      hrf_basis_kernels = matrix(0, 1L, nrow(coefficients)),
      coefficients = coefficients[, voxel_index, drop = FALSE],
      Z = fixed_use,
      Nuisance = nuisance_regs,
      verbose = show_progress,
      method = method,
      basis_convolved = built$basis_convolved
    )
  }

  decorate_dense <- function(beta, engine_used, chunk_size_used) {
    rownames(beta) <- paste0("trial_", seq_len(nrow(events)))
    colnames(beta) <- colnames(Y)
    attr(beta, "sframe") <- sframe
    attr(beta, "normalization") <- normalization
    attr(beta, "units") <- units
    attr(beta, "event_amplitude") <- event_amplitude
    attr(beta, "event_duration") <- event_duration
    attr(beta, "engine_requested") <- engine
    attr(beta, "engine_used") <- engine_used
    attr(beta, "chunk_size") <- chunk_size_used
    beta
  }

  if (engine == "R") {
    betas_dense <- solve_voxels(seq_len(ncol(Y)), "r", verbose)
    return(decorate_dense(
      betas_dense, attr(betas_dense, "engine_used") %||% "r", ncol(Y)
    ))
  }

  if (is.null(backing_dir)) {
    backing_dir <- tempdir()
  }
  bfile <- tempfile("betas", tmpdir = backing_dir, fileext = ".bin")
  descfile <- tempfile("betas", tmpdir = backing_dir, fileext = ".desc")
  betas <- bigmemory::filebacked.big.matrix(
    nrow = nrow(events), ncol = ncol(Y), type = "double",
    backingfile = basename(bfile), descriptorfile = basename(descfile),
    backingpath = dirname(bfile)
  )

  chunks <- split(seq_len(ncol(Y)), ceiling(seq_len(ncol(Y)) / chunk_size))
  engines_used <- character(length(chunks))
  for (chunk_id in seq_along(chunks)) {
    voxel_index <- chunks[[chunk_id]]
    chunk_beta <- solve_voxels(voxel_index, "cpp_arma", FALSE)
    betas[, voxel_index] <- unname(chunk_beta)
    engines_used[chunk_id] <- attr(chunk_beta, "engine_used") %||% "r"
    if (verbose) {
      message("Processed voxel chunk ", chunk_id, " / ", length(chunks))
    }
  }
  engines_used <- unique(engines_used)
  engine_used <- if (length(engines_used) == 1L) engines_used else {
    paste(sort(engines_used), collapse = "+")
  }

  result <- list(
    betas = betas,
    dimnames = list(paste0("trial_", seq_len(nrow(events))), colnames(Y)),
    sframe = sframe,
    normalization = normalization,
    units = units,
    event_amplitude = event_amplitude,
    event_duration = event_duration,
    engine_requested = engine,
    engine_used = engine_used,
    chunk_size = chunk_size
  )
  class(result) <- "LSSBeta"
  result
}

.voxhrf_trial_basis <- function(events, basis, sframe, precision = 0.1,
                                method = "conv", span = NULL) {
  times <- fmrihrf::samples(sframe, global = FALSE)
  blocks <- fmrihrf::blockids(sframe)
  event_run <- .voxhrf_validate_events(events, sframe)
  span <- span %||% if (!is.null(attr(basis, "span"))) attr(basis, "span") else 30
  span <- as.numeric(span)
  if (length(span) != 1L || !is.finite(span) || span <= 0) {
    stop("HRF span must be one positive finite value", call. = FALSE)
  }
  K <- tryCatch(fmrihrf::nbasis(basis), error = function(e) 1L)
  X <- matrix(0, nrow = length(times), ncol = nrow(events) * K)
  for (j in seq_len(nrow(events))) {
    run_rows <- which(blocks == event_run[j])
    rr <- fmrihrf::regressor(
      onsets = events$onset[j], hrf = basis,
      duration = events$duration[j],
      amplitude = if ("amplitude" %in% names(events)) events$amplitude[j] else 1,
      span = span, summate = FALSE
    )
    xj <- fmrihrf::evaluate(
      rr, grid = times[run_rows], precision = precision, method = method
    )
    if (inherits(xj, "Matrix")) xj <- as.matrix(xj)
    if (!is.matrix(xj)) xj <- matrix(xj, ncol = 1L)
    if (ncol(xj) != K) {
      stop("fmrihrf trial design does not match the HRF basis dimension", call. = FALSE)
    }
    cols <- ((j - 1L) * K + 1L):(j * K)
    X[run_rows, cols] <- xj
  }
  colnames(X) <- as.vector(t(outer(
    paste0("trial_", seq_len(nrow(events))), paste0("basis_", seq_len(K)),
    paste, sep = ":"
  )))
  list(
    X = X,
    basis_convolved = lapply(seq_len(K), function(k) {
      X[, seq.int(k, ncol(X), by = K), drop = FALSE]
    }),
    K = K
  )
}

#' @export
as.matrix.LSSBeta <- function(x, ...) {
  out <- as.matrix(x$betas[])
  dimnames(out) <- x$dimnames
  attr(out, "sframe") <- x$sframe
  attr(out, "normalization") <- x$normalization
  attr(out, "units") <- x$units
  attr(out, "event_amplitude") <- x$event_amplitude
  attr(out, "event_duration") <- x$event_duration
  attr(out, "engine_requested") <- x$engine_requested
  attr(out, "engine_used") <- x$engine_used
  attr(out, "chunk_size") <- x$chunk_size
  out
}

# Internal helper: compute VOXHRF betas using an explicit sampling frame.
# Returns a dense numeric matrix (n_trials x n_vox)
.voxhrf_betas_cpp_arma <- function(Y, onsets, durations, hrf_estimates,
                                    nuisance_regs = NULL, sframe = NULL,
                                    fixed_regs = NULL, runs = NULL) {
  coefficients <- .voxhrf_align_coefficients(Y, hrf_estimates)
  sframe <- sframe %||% hrf_estimates$sframe
  if (is.null(sframe) || !inherits(sframe, "sampling_frame")) {
    stop("sframe must be supplied explicitly or stored in hrf_estimates", call. = FALSE)
  }
  events <- data.frame(onset = onsets, duration = durations, condition = "condition")
  if (!is.null(runs)) events$run <- runs
  built <- .voxhrf_trial_basis(events, hrf_estimates$basis, sframe)
  lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = seq_along(onsets),
    durations = rep.int(0L, length(onsets)),
    hrf_basis_kernels = matrix(0, 1L, nrow(coefficients)),
    coefficients = coefficients,
    Z = .voxhrf_fixed_design(fixed_regs, nrow(Y), sframe),
    Nuisance = nuisance_regs,
    method = "cpp_arma",
    verbose = FALSE,
    basis_convolved = built$basis_convolved
  )
}

.voxhrf_align_coefficients <- function(Y, hrf_estimates) {
  if (!inherits(hrf_estimates, "VoxelHRF") ||
      !is.matrix(hrf_estimates$coefficients) ||
      !is.numeric(hrf_estimates$coefficients)) {
    stop("hrf_estimates must be a VoxelHRF object with a numeric coefficient matrix")
  }
  coefficients <- hrf_estimates$coefficients
  if (any(!is.finite(coefficients))) {
    stop("hrf_estimates coefficients must be finite", call. = FALSE)
  }
  if (!inherits(hrf_estimates$basis, "HRF")) {
    stop("hrf_estimates must contain an fmrihrf HRF basis", call. = FALSE)
  }
  expected_k <- tryCatch(
    fmrihrf::nbasis(hrf_estimates$basis),
    error = function(e) NA_integer_
  )
  if (!is.finite(expected_k) || nrow(coefficients) != expected_k) {
    stop("hrf_estimates coefficient rows must match the HRF basis dimension", call. = FALSE)
  }
  if (ncol(coefficients) != ncol(Y)) {
    stop("hrf_estimates must have one coefficient column per Y voxel", call. = FALSE)
  }
  y_names <- colnames(Y)
  coefficient_names <- colnames(coefficients)
  if (xor(is.null(y_names), is.null(coefficient_names))) {
    stop("Y and hrf_estimates must either both have voxel names or both be unnamed", call. = FALSE)
  }
  if (!is.null(y_names)) {
    if (anyNA(y_names) || any(!nzchar(y_names)) || anyDuplicated(y_names) ||
        anyNA(coefficient_names) || any(!nzchar(coefficient_names)) ||
        anyDuplicated(coefficient_names) || !setequal(y_names, coefficient_names)) {
      stop("Y and hrf_estimates must have the same unique voxel names", call. = FALSE)
    }
    coefficients <- coefficients[, match(y_names, coefficient_names), drop = FALSE]
  }
  coefficients
}

.voxhrf_validate_response <- function(Y) {
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix", call. = FALSE)
  }
  if (nrow(Y) < 1L || ncol(Y) < 1L || any(!is.finite(Y))) {
    stop("Y must be a non-empty matrix of finite values", call. = FALSE)
  }
  voxel_names <- colnames(Y)
  if (!is.null(voxel_names) &&
      (anyNA(voxel_names) || any(!nzchar(voxel_names)) || anyDuplicated(voxel_names))) {
    stop("Y voxel names must be complete and unique", call. = FALSE)
  }
  invisible(TRUE)
}

.voxhrf_validate_events <- function(events, sframe) {
  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition", call. = FALSE)
  }
  if (nrow(events) < 1L) {
    stop("events must contain at least one row", call. = FALSE)
  }
  if (!is.numeric(events$onset) || !is.numeric(events$duration) ||
      any(!is.finite(events$onset)) || any(!is.finite(events$duration))) {
    stop("events onset and duration values must be finite numeric values", call. = FALSE)
  }
  if (any(events$duration < 0)) {
    stop("events duration values must be nonnegative", call. = FALSE)
  }
  condition <- as.character(events$condition)
  if (anyNA(condition) || any(!nzchar(condition))) {
    stop("events condition values must be complete and non-empty", call. = FALSE)
  }
  if ("amplitude" %in% names(events) &&
      (!is.numeric(events$amplitude) || any(!is.finite(events$amplitude)))) {
    stop("events amplitude values must be finite numeric values", call. = FALSE)
  }

  blocks <- fmrihrf::blockids(sframe)
  n_blocks <- length(unique(blocks))
  if (n_blocks > 1L && !"run" %in% names(events)) {
    stop("multi-run voxel-HRF events require an explicit run column", call. = FALSE)
  }
  event_run_numeric <- if ("run" %in% names(events)) {
    suppressWarnings(as.numeric(as.character(events$run)))
  } else {
    rep.int(1L, nrow(events))
  }
  if (anyNA(event_run_numeric) || any(!is.finite(event_run_numeric)) ||
      any(event_run_numeric != round(event_run_numeric))) {
    stop("events$run must contain exact integer run indices", call. = FALSE)
  }
  event_run <- as.integer(event_run_numeric)
  if (any(!event_run %in% seq_len(n_blocks))) {
    stop("events$run must contain sampling-frame run indices", call. = FALSE)
  }

  local_times <- fmrihrf::samples(sframe, global = FALSE)
  for (run in seq_len(n_blocks)) {
    event_rows <- which(event_run == run)
    if (!length(event_rows)) next
    scan_rows <- which(blocks == run)
    run_times <- local_times[scan_rows]
    half_step <- if (length(run_times) > 1L) {
      stats::median(diff(run_times)) / 2
    } else {
      0
    }
    lower <- min(run_times) - half_step
    upper <- max(run_times) + half_step
    tolerance <- 100 * .Machine$double.eps * max(1, abs(lower), abs(upper))
    outside <- events$onset[event_rows] < lower - tolerance |
      events$onset[event_rows] > upper + tolerance
    if (any(outside)) {
      stop("events$onset must lie within its sampling-frame run", call. = FALSE)
    }
  }
  event_run
}

#' Rank-1 GLM: joint voxel-wise HRF and trial amplitude estimation
#'
#' Fits the rank-1 GLM of Pedregosa et al. (2015): within each voxel one HRF,
#' expressed in an HRF basis, is shared by all trials, and each trial has its
#' own amplitude. With `model = "separate"` (R1-GLMS, the default and the best
#' performer in that paper) every trial is fitted in its own least-squares
#' separate (LSS) model, so the amplitudes are LSS estimates under the
#' voxel's estimated HRF. With `model = "joint"` (R1-GLM) all trials enter one
#' least-squares-all model.
#'
#' @details
#' For trial \eqn{i} with basis-convolved design block \eqn{X_i} (n x K) and
#' the sum of all blocks \eqn{T_X}, the separate model minimizes, per voxel,
#' \deqn{\sum_i \| y - \beta_i X_i h - r_i (T_X - X_i) h \|^2}
#' over the HRF coefficients \eqn{h}, trial amplitudes \eqn{\beta} and
#' other-trial amplitudes \eqn{r}; the joint model minimizes
#' \eqn{\| y - \sum_i \beta_i X_i h \|^2}. Fixed regressors (run intercepts are
#' added when absent) and nuisance regressors are projected from the design
#' and data first, which gives each trial-wise model its own confound
#' coefficients, as in [lss()]. (The original formulation shares one set of
#' confound coefficients across the trial-wise models.)
#'
#' The fit uses exact alternating least squares rather than the quasi-Newton
#' solver of the paper: given \eqn{h}, the amplitudes are ordinary LSS (or
#' LSA) estimates; given the amplitudes, \eqn{h} solves a K-dimensional least
#' squares problem. Each step minimizes its block exactly, so the objective
#' never increases. All quantities either step needs are K x K Gram blocks of
#' the residualized design and the product \eqn{X^\top Y}, so after that one
#' product the separate model costs \eqn{O(T K^2)} per voxel and iteration,
#' independent of the number of scans. The joint model costs
#' \eqn{O(T^2 K^2 + T^3)} per voxel and iteration and is intended for regions
#' of interest.
#'
#' The problem is bilinear: \eqn{h} and the amplitudes are determined only up
#' to a common scale and sign. After fitting, each voxel's HRF is oriented to
#' correlate positively with `ref_hrf` and scaled to a unit positive peak, and
#' the amplitudes absorb that scale, so they are in peak-response units for
#' unit-amplitude, zero-duration events, like [lss_with_hrf()]. The returned
#' `hrf` is a [VoxelHRF] object and can be passed to [lss_with_hrf()] to
#' estimate trial amplitudes in new data (for example a held-out run) with the
#' learned HRFs.
#'
#' The objective is not jointly convex. `init = "aggregate"` starts each voxel
#' from the pooled-shape fit of [estimate_voxel_hrf()] (one common amplitude
#' for all trials); `init = "reference"` starts every voxel from the
#' projection of `ref_hrf` onto the basis, which is more stable for voxels
#' whose mean response is near zero. A K x V matrix of starting coefficients
#' is also accepted.
#'
#' @param Y Numeric matrix of BOLD data (time x voxels).
#' @param events Data frame with `onset`, `duration` and `condition` columns
#'   (and `run` for multi-run sampling frames, with run-relative onsets), as
#'   for [estimate_voxel_hrf()]. Each row is one trial. Conditions are labels
#'   only: all trials share each voxel's HRF.
#' @param basis An `fmrihrf` HRF basis, e.g. `fmrihrf::HRF_SPMG3` or an FIR
#'   basis from `fmrihrf::hrf_fir_generator()`/`fmrihrf::HRF_FIR`.
#' @param sframe An `fmrihrf` sampling frame with `nrow(Y)` scans.
#' @param nuisance_regs Optional numeric matrix of nuisance regressors.
#' @param fixed_regs Optional numeric matrix of fixed regressors (run
#'   intercepts are added when not already spanned).
#' @param model `"separate"` (R1-GLMS) or `"joint"` (R1-GLM).
#' @param trial_groups Optional vector with one condition label per trial
#'   (row of `events`), e.g. `events$condition`, for `model = "separate"`.
#'   Each trial-wise model then has one "other trials" regressor per group
#'   (LSS-N, as in [lss()]). Because the HRF is learned mostly through these
#'   group regressors, supply `trial_groups` whenever conditions may differ in
#'   mean response, in particular in sign. See Details.
#' @param init `"aggregate"`, `"reference"`, or a K x V numeric matrix (or
#'   length-K vector) of starting HRF coefficients.
#' @param ref_hrf HRF used for the `"reference"` initialization and to orient
#'   the sign of each estimated HRF. Defaults to `fmrihrf::HRF_SPMG1`.
#' @param prewhiten Optional prewhitening options (see [prewhiten_options()]).
#'   Global and run pooling are supported; the noise model is estimated with
#'   one aggregate regressor per basis function.
#' @param solver `"als"` (default): exact alternating least squares, run in
#'   parallel over voxels with OpenMP. `"lbfgs"`: joint quasi-Newton
#'   optimization of HRF and amplitudes with L-BFGS-B, as in Pedregosa et al.
#'   (2015), on the same objective and starting point (serial; separate model
#'   only). See Details.
#' @param max_iter Maximum number of alternating iterations (for `"lbfgs"`,
#'   at least 1000 quasi-Newton iterations are allowed).
#' @param tol Relative objective change used as the convergence criterion
#'   (for `"lbfgs"`, passed as `factr = tol / .Machine$double.eps`).
#'
#' @return A list with
#' \describe{
#'   \item{beta}{Trial amplitudes (trials x voxels), in peak-response units.}
#'   \item{other}{For `model = "separate"`, the amplitude of the pooled
#'     other-trials regressor in each trial's model (trials x voxels);
#'     `NULL` for the joint model.}
#'   \item{hrf}{A [VoxelHRF] object with the unit-peak HRF coefficients
#'     (K x voxels), the removed `amplitude_scale`, `basis` and `sframe`.}
#'   \item{objective}{Final residual sum of squares per voxel (summed over
#'     the trial-wise models for `"separate"`).}
#'   \item{iterations, converged}{Per-voxel iteration counts and convergence
#'     flags.}
#'   \item{degenerate}{Voxels whose estimated HRF has no positive peak after
#'     orientation (scaled by its largest absolute value instead).}
#'   \item{model}{The fitted model.}
#' }
#' When prewhitening is applied the fitted plan is attached as the
#' `whiten_plan` attribute.
#'
#' @references
#' Pedregosa, F., Eickenberg, M., Ciuciu, P., Thirion, B., & Gramfort, A.
#' (2015). Data-driven HRF estimation for encoding and decoding models.
#' NeuroImage, 104, 209-220.
#'
#' @seealso [estimate_voxel_hrf()], [lss_with_hrf()], [lss_sbhm()]
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' set.seed(1)
#' sframe <- fmrihrf::sampling_frame(blocklens = 200, TR = 1)
#' events <- data.frame(onset = seq(10, 180, by = 9), duration = 0,
#'                      condition = "A")
#' basis <- fmrihrf::HRF_SPMG3
#' X <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)$X
#' h <- c(1, 0.4, -0.2)
#' beta <- matrix(rnorm(nrow(events) * 5, 1, 0.3), nrow(events), 5)
#' Y <- X %*% kronecker(beta, h) + matrix(rnorm(200 * 5, sd = 0.5), 200, 5)
#' fit <- lss_rank1(Y, events, basis, sframe)
#' dim(fit$beta)
#' fit$hrf$coefficients[, 1]
#' }
#' @export
lss_rank1 <- function(Y, events, basis, sframe, nuisance_regs = NULL,
                      fixed_regs = NULL, model = c("separate", "joint"),
                      trial_groups = NULL,
                      init = "aggregate", ref_hrf = NULL, prewhiten = NULL,
                      solver = c("als", "lbfgs"), max_iter = 100L, tol = 1e-7) {
  model <- match.arg(model)
  solver <- match.arg(solver)
  if (solver == "lbfgs" && model != "separate") {
    stop("solver = 'lbfgs' is available for model = 'separate'", call. = FALSE)
  }
  groups <- .lss_group_codes(trial_groups, nrow(events))
  group_labels <- if (is.null(groups)) NULL else {
    as.character(unique(if (is.factor(trial_groups)) as.character(trial_groups) else trial_groups))
  }
  if (!is.null(groups) && model != "separate") {
    stop("trial_groups applies to model = 'separate'", call. = FALSE)
  }
  if (!is.null(groups) && solver == "lbfgs") {
    stop("solver = 'lbfgs' supports one pooled other-trials regressor; ",
         "use solver = 'als' with trial_groups", call. = FALSE)
  }
  .voxhrf_validate_response(Y)
  if (!inherits(basis, "HRF")) {
    stop("basis must be an 'HRF' object from fmrihrf", call. = FALSE)
  }
  if (is.null(sframe) || !inherits(sframe, "sampling_frame")) {
    stop("sframe must be an explicit fmrihrf sampling_frame", call. = FALSE)
  }
  if (length(fmrihrf::samples(sframe, global = TRUE)) != nrow(Y)) {
    stop("sframe must contain exactly nrow(Y) scan times", call. = FALSE)
  }
  for (nm in c("nuisance_regs", "fixed_regs")) {
    M <- get(nm)
    if (!is.null(M) && (!is.matrix(M) || !is.numeric(M) || nrow(M) != nrow(Y) ||
                        any(!is.finite(M)))) {
      stop(nm, " must be a finite numeric matrix with nrow(Y) rows", call. = FALSE)
    }
  }
  max_iter <- .as_positive_integer(max_iter, "max_iter")
  tol <- .as_nonnegative_scalar(tol, "tol")
  ref_hrf <- ref_hrf %||% fmrihrf::HRF_SPMG1
  if (!inherits(ref_hrf, "HRF")) stop("ref_hrf must be an fmrihrf HRF", call. = FALSE)

  built <- .voxhrf_trial_basis(events, basis, sframe)
  X <- built$X
  K <- built$K
  n_trials <- nrow(events)
  fixed <- .voxhrf_fixed_design(fixed_regs, nrow(Y), sframe)

  # Optional prewhitening: one aggregate regressor per basis function feeds
  # the noise model, then data and every design are filtered alike.
  whiten_plan <- NULL
  if (!is.null(prewhiten) && !identical(prewhiten$method %||% "ar", "none")) {
    opts <- .resolve_prewhiten_options(prewhiten, internal = FALSE)
    if (opts$pooling %in% c("voxel", "parcel")) {
      stop("lss_rank1() supports prewhiten pooling = 'global' or 'run'",
           call. = FALSE)
    }
    X_noise <- do.call(cbind, lapply(built$basis_convolved, rowSums))
    wh <- .prewhiten_data(Y, X, fixed, nuisance_regs, opts, X_noise = X_noise)
    whiten_plan <- wh$whiten_plan
    Y <- wh$Y_whitened
    X <- wh$X_whitened
    fixed <- wh$Z_whitened
    if (!is.null(nuisance_regs)) nuisance_regs <- wh$Nuisance_whitened
  }

  # Residualize the trial design against the common span. Y itself is never
  # residualized: U = X_res'Y = X_res'Y_res.
  common <- .voxhrf_orthonormal_span(cbind(fixed, nuisance_regs))
  resid <- function(M) {
    if (!ncol(common)) return(M)
    M - common %*% crossprod(common, M)
  }
  Xr <- resid(X)
  blocks <- lapply(seq_len(n_trials), function(i) {
    Xr[, ((i - 1L) * K + 1L):(i * K), drop = FALSE]
  })
  TX <- Reduce(`+`, blocks)
  GTT <- crossprod(TX)
  ev <- eigen(GTT, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10 * max(ev)) {
    stop("HRF basis is not identifiable after projecting the common design",
         call. = FALSE)
  }

  U <- crossprod(Xr, Y)
  yy <- colSums(Y^2) - if (ncol(common)) colSums(crossprod(common, Y)^2) else 0

  # Basis waveforms on a fine grid, for initialization and normalization
  span <- if (!is.null(attr(basis, "span"))) attr(basis, "span") else 30
  wf <- .rank1_waveforms(basis, ref_hrf, span)
  H0 <- .rank1_init(init, K, ncol(Y), U, GTT, n_trials, wf)

  fit <- if (model == "separate") {
    G <- if (is.null(groups)) 1L else max(groups)
    codes <- if (is.null(groups)) rep(1L, n_trials) else groups
    A <- lapply(seq_len(G), function(g) Reduce(`+`, blocks[codes == g]))
    A_all <- do.call(cbind, A)
    GA <- crossprod(A_all)
    Gii <- array(0, c(K, K, n_trials))
    S <- array(0, c(K, K, n_trials * G))
    for (i in seq_len(n_trials)) {
      Gii[, , i] <- crossprod(blocks[[i]])
      XA <- crossprod(blocks[[i]], A_all)
      for (g in seq_len(G)) {
        S[, , (i - 1L) * G + g] <- XA[, ((g - 1L) * K + 1L):(g * K), drop = FALSE]
      }
    }
    zero_codes <- as.integer(codes - 1L)
    if (solver == "als") {
      r1glms_fit_cpp(U, Gii, S, GA, zero_codes, yy, H0, max_iter, tol)
    } else {
      # Same starting point as ALS: h0 and the LSS amplitudes given h0.
      start <- r1glms_fit_cpp(U, Gii, S, GA, zero_codes, yy, H0, 0L, tol)
      r1glms_lbfgs_cpp(U, Gii, S, GA, yy, start$h, start$beta, start$other,
                       maxit = max(max_iter, 1000L), factr = tol / .Machine$double.eps)
    }
  } else {
    r1glm_fit_cpp(U, crossprod(Xr), yy, H0, max_iter, tol)
  }

  # Resolve the scale/sign ambiguity: unit positive peak, oriented to ref_hrf.
  waves <- wf$H %*% fit$h
  orientation <- sign(drop(crossprod(wf$ref, waves)))
  orientation[!is.finite(orientation) | orientation == 0] <- 1
  peak <- apply(sweep(waves, 2L, orientation, "*"), 2L, max)
  amp <- apply(abs(waves), 2L, max)
  degenerate <- !is.finite(peak) | peak <= 1e-8 * pmax(amp, .Machine$double.xmin)
  peak[degenerate] <- amp[degenerate]
  peak[!is.finite(peak) | peak == 0] <- 1
  scale <- orientation * peak

  voxel_names <- colnames(Y)
  trial_names <- paste0("trial_", seq_len(n_trials))
  coefficients <- sweep(fit$h, 2L, scale, "/")
  dimnames(coefficients) <- list(paste0("basis_", seq_len(K)), voxel_names)
  beta <- sweep(fit$beta, 2L, scale, "*")
  dimnames(beta) <- list(trial_names, voxel_names)
  other <- if (model == "separate") {
    out <- sweep(fit$other, 2L, scale, "*")
    if (is.null(groups)) {
      dimnames(out) <- list(trial_names, voxel_names)
    } else {
      # rows are trial-major (trial i, group g) -> trials x groups x voxels
      G <- max(groups)
      out <- aperm(array(out, c(G, n_trials, ncol(out))), c(2L, 1L, 3L))
      dimnames(out) <- list(trial_names, group_labels, voxel_names)
    }
    out
  } else {
    NULL
  }
  names(scale) <- voxel_names

  hrf <- structure(list(
    coefficients = coefficients,
    amplitude_scale = scale,
    basis = basis,
    conditions = unique(as.character(events$condition)),
    sframe = sframe,
    condition_pooling = "all-events",
    normalization = "positive-peak",
    coefficient_units = "unit-peak HRF shape weights",
    estimator = paste0("rank1_", model)
  ), class = "VoxelHRF")

  result <- list(
    beta = beta,
    other = other,
    hrf = hrf,
    objective = stats::setNames(drop(fit$objective), voxel_names),
    iterations = stats::setNames(as.integer(fit$iterations), voxel_names),
    converged = stats::setNames(as.logical(fit$converged), voxel_names),
    degenerate = stats::setNames(degenerate, voxel_names),
    model = model
  )
  .attach_whiten_plan(result, whiten_plan)
}

#' Basis and reference waveforms on a fine grid
#' @keywords internal
#' @noRd
.rank1_waveforms <- function(basis, ref_hrf, span, precision = 0.05) {
  grid <- seq(0, span, by = precision)
  eval_hrf <- function(h) {
    impulse <- fmrihrf::regressor(onsets = 0, hrf = h, duration = 0, span = span)
    out <- fmrihrf::evaluate(impulse, grid, precision = precision, method = "conv")
    if (inherits(out, "Matrix")) out <- as.matrix(out)
    if (!is.matrix(out)) out <- matrix(out, ncol = 1L)
    sweep(out, 2L, out[1L, ], "-")
  }
  H <- eval_hrf(basis)
  ref <- eval_hrf(ref_hrf)[, 1L]
  list(H = H, ref = ref, ref_coef = qr.coef(qr(H), ref))
}

#' Starting HRF coefficients (K x V)
#' @keywords internal
#' @noRd
.rank1_init <- function(init, K, V, U, GTT, n_trials, wf) {
  if (is.numeric(init)) {
    H0 <- if (is.matrix(init)) init else matrix(init, K, V)
    if (!identical(dim(H0), c(K, V)) || any(!is.finite(H0))) {
      stop("numeric init must be a finite K x V matrix or length-K vector",
           call. = FALSE)
    }
    return(H0)
  }
  init <- match.arg(init, c("aggregate", "reference"))
  ref <- wf$ref_coef
  ref[!is.finite(ref)] <- 0
  if (init == "reference" || all(ref == 0)) {
    if (all(ref == 0)) ref <- c(1, rep(0, K - 1L))
    return(matrix(ref, K, V))
  }
  # Pooled-shape fit: T_X h = y, i.e. GTT h = sum_i U_i
  TY <- Reduce(`+`, lapply(seq_len(n_trials), function(i) {
    U[((i - 1L) * K + 1L):(i * K), , drop = FALSE]
  }))
  H0 <- solve(GTT, TY)
  bad <- colSums(H0^2) == 0 | !is.finite(colSums(H0))
  H0[, bad] <- ref
  H0
}

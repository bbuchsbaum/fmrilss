#' SBHM Pipeline with fmridesign Models
#'
#' Run the SBHM end-to-end pipeline using fmridesign's `event_model` and
#' optional `baseline_model`, mirroring the convenience of `lss_design()` but
#' producing SBHM shape coefficients and (optionally) scalar event-design trial
#' coefficients in the component named `amplitude`.
#'
#' This function wraps `lss_sbhm()` by locating the unique trialwise event term,
#' rebuilding it with the SBHM basis, and retaining run identities and
#' heterogeneous durations. Non-trial event terms and `baseline_model` columns
#' enter the fixed nuisance span. Multi-run regressors are constructed within
#' runs so HRF tails cannot cross run boundaries.
#'
#' @param Y Numeric matrix T×V of fMRI time series (timepoints × voxels).
#' @param sbhm SBHM object as returned by `sbhm_build()`.
#' @param event_model An `event_model` from `fmridesign::event_model()` defining
#'   the trial structure (typically created with `trialwise()`).
#' @param baseline_model Optional `baseline_model` from
#'   `fmridesign::baseline_model()`. Its drift, block, and nuisance terms are
#'   projected out as confounds.
#' @param prewhiten Optional prewhitening options (see `?lss`).
#' @param prepass Optional list forwarded to `sbhm_prepass()`.
#' @param match Optional list forwarded to `sbhm_match()`.
#' @param oasis Optional SBHM OASIS override list. Supported fields are
#'   `ridge_mode`, `ridge_x`, `ridge_b`, and `block_cols`; false-valued
#'   `return_se` and `return_diag` are accepted, while true values fail because
#'   SBHM uncertainty/design diagnostics are not exposed. The basis rank,
#'   trial count, intercept policy, and trial/basis map are fixed internally
#'   from `sbhm` and the event model.
#' @param amplitude Amplitude options (see `?lss_sbhm`).
#' @param return One of `"amplitude"`, `"coefficients"`, or `"both"`.
#' @param validate Logical; when TRUE, performs basic checks (sampling frame
#'   compatibility, temporal alignment) analogous to `lss_design()`.
#' @param ... Reserved for future use.
#'
#' @return Same return contract as `lss_sbhm()`.
#'
#' @examplesIf requireNamespace("fmridesign", quietly = TRUE)
#' \donttest{
#'   library(fmridesign)
#'   sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
#'   trials <- data.frame(onset = c(5, 20, 35), run = 1)
#'   emod <- event_model(onset ~ trialwise(basis = "spmg1"), data = trials,
#'                       block = ~run, sampling_frame = sframe)
#'   times <- fmrihrf::samples(sframe, global = TRUE)
#'   H <- cbind(exp(-times / 5), exp(-times / 7))
#'   sbhm <- sbhm_build(library_H = H, r = 2, sframe = sframe,
#'                      normalize = TRUE, baseline = NULL)
#'   Y <- matrix(rnorm(60 * 4), 60, 4)
#'   out <- lss_sbhm_design(Y, sbhm, emod)
#' }
#'
#' @export
lss_sbhm_design <- function(Y, sbhm, event_model, baseline_model = NULL,
                            prewhiten = NULL,
                            prepass = list(),
                            match = list(),
                            oasis = list(),
                            amplitude = list(),
                            return = c("amplitude", "coefficients", "both"),
                            validate = TRUE,
                            ...) {

  return <- match.arg(return)
  validate <- .as_scalar_logical(validate, "validate")
  dots <- list(...)
  if (length(dots)) {
    stop("... is reserved; unsupported argument(s): ",
         paste(names(dots) %||% rep("<unnamed>", length(dots)), collapse = ", "),
         call. = FALSE)
  }

  # Basic validation ----------------------------------------------------------
  if (!requireNamespace("fmridesign", quietly = TRUE)) {
    stop("Package 'fmridesign' is required. Install it to use lss_sbhm_design().",
         call. = FALSE)
  }
  if (missing(event_model) || is.null(event_model) || !inherits(event_model, "event_model")) {
    stop("event_model must be an 'event_model' from fmridesign::event_model()", call. = FALSE)
  }
  if (!is.null(baseline_model) && !inherits(baseline_model, "baseline_model")) {
    stop("baseline_model must be a 'baseline_model' from fmridesign::baseline_model()", call. = FALSE)
  }

  sframe <- event_model$sampling_frame
  if (!is.null(baseline_model)) {
    if (!identical(baseline_model$sampling_frame, sframe)) {
      stop("event_model and baseline_model must use the same sampling_frame", call. = FALSE)
    }
  }

  # Temporal alignment (optional)
  if (validate) {
    expected_scans <- sum(fmrihrf::blocklens(sframe))
    if (nrow(Y) != expected_scans) {
      stop(sprintf("Y has %d rows but sampling_frame expects %d scans",
                   nrow(Y), expected_scans), call. = FALSE)
    }
  }

  # Identify the unique trialwise target term and preserve all other event
  # terms as fixed regressors.
  event_dm <- fmridesign::design_matrix(event_model)
  prepared <- .prepare_fmridesign_lss(event_dm, event_model)
  spec <- .sbhm_design_spec_from_event_model(event_model, sbhm, prepared)

  # Map baseline_model -> Nuisance (project out)
  Nuisance <- prepared$fixed
  if (!is.null(baseline_model)) {
    bm_term_mats <- fmridesign::term_matrices(baseline_model)
    # Project all baseline pieces as nuisances for SBHM (Z vs Nuis identical in OASIS projection)
    nuis_list <- list()
    if ("drift" %in% names(bm_term_mats)) nuis_list$drift <- bm_term_mats$drift
    if ("block" %in% names(bm_term_mats)) nuis_list$block <- bm_term_mats$block
    if ("nuisance" %in% names(bm_term_mats)) nuis_list$nuis <- bm_term_mats$nuisance
    if (length(nuis_list) > 0) {
      Nuisance <- cbind(Nuisance, as.matrix(do.call(cbind, nuis_list)))
    }
  }

  # Call the existing SBHM pipeline with design_spec -------------------------
  out <- lss_sbhm(
    Y = Y, sbhm = sbhm, design_spec = spec,
    Nuisance = Nuisance,
    prewhiten = prewhiten,
    prepass = prepass,
    match = match,
    oasis = oasis,
    amplitude = amplitude,
    return = return
  )

  # Attach provenance
  attr(out, "event_model") <- event_model
  attr(out, "baseline_model") <- baseline_model
  attr(out, "sampling_frame") <- sframe
  attr(out, "method") <- "lss_sbhm_design"
  attr(out, "trial_basis_map") <- out$trial_basis_map
  out
}

# Internal: Convert event_model -> design_spec using SBHM basis ---------------
#' @keywords internal
.sbhm_design_spec_from_event_model <- function(event_model, sbhm, prepared = NULL) {
  sframe <- event_model$sampling_frame
  terms <- event_model$terms
  if (length(terms) < 1) stop("event_model has no terms to build a design from")
  trial_tags <- names(Filter(function(term) {
    ids <- term$condition_ids
    onsets <- term$onsets
    length(ids) == length(onsets) && length(ids) > 0L && !anyDuplicated(ids)
  }, terms))
  if (length(trial_tags) != 1L) {
    stop("lss_sbhm_design() requires exactly one trialwise event term", call. = FALSE)
  }
  t1 <- terms[[trial_tags]]

  # Onsets/durations and run ids
  onsets_run <- as.numeric(t1$onsets)
  durs_run   <- if (!is.null(t1$durations)) as.numeric(t1$durations) else rep(0, length(onsets_run))
  runs       <- as.integer(t1$blockids)
  if (anyNA(onsets_run) || anyNA(runs)) stop("Could not extract onsets/blockids from event_model")
  if (length(durs_run) != length(onsets_run) || any(!is.finite(durs_run)) || any(durs_run < 0)) {
    stop("Could not extract valid trial durations from event_model", call. = FALSE)
  }
  if (is.null(prepared)) {
    prepared <- .prepare_fmridesign_lss(
      fmridesign::design_matrix(event_model), event_model
    )
  }
  trial_names <- unique(prepared$map$trial_name)
  if (length(trial_names) != length(onsets_run)) {
    stop("fmridesign trial identities disagree with event timing metadata", call. = FALSE)
  }

  list(
    sframe = sframe,
    cond = list(
      onsets = onsets_run,
      duration = durs_run,
      amplitude = rep.int(1, length(onsets_run)),
      run = runs,
      trial_names = trial_names,
      span = if (is.null(sbhm$span)) 40 else sbhm$span
    ),
    precision = if (is.null(event_model$model_spec$precision)) 0.1 else event_model$model_spec$precision,
    method = "conv"
  )
}

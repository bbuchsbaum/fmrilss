#' SBHM Pipeline with fmridesign Models
#'
#' Run the SBHM end-to-end pipeline using fmridesign's `event_model` and
#' optional `baseline_model`, mirroring the convenience of `lss_design()` but
#' producing SBHM coefficients and (optionally) scalar amplitudes.
#'
#' This function wraps `lss_sbhm()` by converting the `event_model` into an
#' OASIS `design_spec` that uses the SBHM basis HRF, and by mapping
#' `baseline_model` terms to nuisance regressors for projection.
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
#' @param oasis Optional list forwarded to `lss(..., method = "oasis")`. `K`
#'   defaults to `ncol(sbhm$B)`.
#' @param amplitude Amplitude options (see `?lss_sbhm`).
#' @param return One of `"amplitude"`, `"coefficients"`, or `"both"`.
#' @param validate Logical; when TRUE, performs basic checks (sampling frame
#'   compatibility, temporal alignment) analogous to `lss_design()`.
#' @param ... Reserved for future use.
#'
#' @return Same return contract as `lss_sbhm()`.
#'
#' @examples
#' \dontrun{
#'   library(fmridesign)
#'   sframe <- fmrihrf::sampling_frame(blocklens = c(150, 150), TR = 2)
#'   trials <- data.frame(onset = c(10,30,50, 10,30,50), run = rep(1:2, each=3))
#'   emod <- event_model(onset ~ trialwise(basis = "spmg1"), data = trials,
#'                       block = ~run, sampling_frame = sframe)
#'   # sbhm <- sbhm_build(...)
#'   Y <- matrix(rnorm(300*100), 300, 100)
#'   out <- lss_sbhm_design(Y, sbhm, emod)
#' }
#'
#' @export
lss_sbhm_design <- function(Y, sbhm, event_model, baseline_model = NULL,
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
                            return = c("amplitude", "coefficients", "both"),
                            validate = TRUE,
                            ...) {

  return <- match.arg(return)

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

  # Build design_spec (SBHM basis) from event_model --------------------------
  spec <- .sbhm_design_spec_from_event_model(event_model, sbhm)

  # Map baseline_model -> Nuisance (project out)
  Nuisance <- NULL
  if (!is.null(baseline_model)) {
    bm_term_mats <- fmridesign::term_matrices(baseline_model)
    # Project all baseline pieces as nuisances for SBHM (Z vs Nuis identical in OASIS projection)
    nuis_list <- list()
    if ("drift" %in% names(bm_term_mats)) nuis_list$drift <- bm_term_mats$drift
    if ("block" %in% names(bm_term_mats)) nuis_list$block <- bm_term_mats$block
    if ("nuisance" %in% names(bm_term_mats)) nuis_list$nuis <- bm_term_mats$nuisance
    if (length(nuis_list) > 0) Nuisance <- as.matrix(do.call(cbind, nuis_list))
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
  out
}

# Internal: Convert event_model -> design_spec using SBHM basis ---------------
#' @keywords internal
.sbhm_design_spec_from_event_model <- function(event_model, sbhm) {
  sframe <- event_model$sampling_frame
  # Extract from first term (assumes trialwise design)
  terms <- event_model$terms
  if (length(terms) < 1) stop("event_model has no terms to build a design from")
  t1 <- terms[[1]]

  # Onsets/durations and run ids
  onsets_run <- as.numeric(t1$onsets)
  durs_run   <- if (!is.null(t1$durations)) as.numeric(t1$durations) else rep(0, length(onsets_run))
  runs       <- as.integer(t1$blockids)
  if (anyNA(onsets_run) || anyNA(runs)) stop("Could not extract onsets/blockids from event_model")

  # Convert run-relative to global seconds
  bl <- fmrihrf::blocklens(sframe)
  TRvec <- tryCatch(as.numeric(sframe$TR), error = function(e) NULL)
  if (is.null(TRvec) || length(TRvec) == 0L) {
    times <- fmrihrf::samples(sframe, global = TRUE)
    TR <- as.numeric(median(diff(times)))
    TRvec <- rep(TR, length(bl))
  } else if (length(TRvec) == 1L) {
    TRvec <- rep(TRvec, length(bl))
  }
  run_starts <- c(0, cumsum(head(bl * TRvec, -1)))
  onsets_global <- onsets_run + run_starts[pmax(1, runs)]

  # SBHM basis HRF
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

  # Condense durations to a scalar (SBHM builds per-trial regressors internally)
  dur_scalar <- 0
  if (length(durs_run) > 0L) {
    unq <- unique(as.numeric(durs_run))
    if (length(unq) == 1L) dur_scalar <- unq
  }

  list(
    sframe = sframe,
    cond = list(
      onsets = onsets_global,
      duration = dur_scalar,
      hrf = hrf_B,
      span = if (is.null(sbhm$span)) 40 else sbhm$span
    ),
    precision = if (is.null(event_model$model_spec$precision)) 0.1 else event_model$model_spec$precision,
    method = "conv"
  )
}

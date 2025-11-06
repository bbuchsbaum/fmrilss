#' LSS Analysis with fmridesign Objects
#'
#' Perform Least Squares Separate (LSS) analysis using event_model and
#' baseline_model objects from the fmridesign package. This provides a
#' streamlined interface for complex designs with multi-condition, parametric
#' modulators, and structured nuisance handling.
#'
#' @param Y Numeric matrix of fMRI data (timepoints × voxels).
#' @param event_model An event_model object from \code{fmridesign::event_model()}.
#'   This defines the trial-wise or condition-wise task design. For LSS, typically
#'   created with \code{trialwise()} to generate one regressor per trial.
#' @param baseline_model Optional baseline_model object from
#'   \code{fmridesign::baseline_model()}. Defines drift correction, block intercepts,
#'   and nuisance regressors. If NULL, basic baseline intercepts are
#'   auto-injected: per-run intercepts derived from `blockids` (or the
#'   sampling frame) are used to ensure proper baseline modeling.
#' @param method LSS method to use. Currently only "oasis" is supported for
#'   event_model integration.
#' @param oasis List of OASIS-specific options (see \code{?lss} for details).
#'   Note: \code{design_spec} is not used when providing event_model.
#' @param prewhiten Optional prewhitening specification (see \code{?lss}).
#' @param blockids Optional block/run identifiers for event_model. If NULL,
#'   extracted from event_model$blockids.
#' @param validate Logical. If TRUE (default), performs validation checks on
#'   design compatibility, collinearity, and temporal alignment.
#' @param ... Additional arguments passed to the underlying LSS method.
#'
#' @return Matrix of trial-wise beta estimates (trials × voxels), or
#'   (trials × basis_functions) × voxels for multi-basis HRFs.
#'
#' @details
#' \strong{Design Specification:}
#'
#' The \code{event_model} should typically use \code{trialwise()} for LSS:
#' \preformatted{
#'   emod <- event_model(onset ~ trialwise(basis = "spmg1"),
#'                       data = events,
#'                       block = ~run,
#'                       sampling_frame = sframe)
#' }
#'
#' For factorial designs (e.g., estimating condition-level betas separately):
#' \preformatted{
#'   emod <- event_model(onset ~ hrf(condition),
#'                       data = events,
#'                       block = ~run,
#'                       sampling_frame = sframe)
#' }
#'
#' \strong{Baseline Model:}
#'
#' If provided, baseline_model components are mapped as follows:
#' \itemize{
#'   \item \code{drift} and \code{block} terms → Z parameter (fixed effects)
#'   \item \code{nuisance} term → Nuisance parameter (confounds)
#' }
#'
#' \strong{Multi-Run Handling:}
#'
#' Both event_model and baseline_model must use the same \code{sampling_frame}.
#' Run structure is automatically respected. Event onsets should be run-relative
#' (resetting to 0 each run) as per fmridesign convention - conversion to global
#' time is handled automatically.
#'
#' \strong{Validation:}
#'
#' When \code{validate = TRUE}, the function checks:
#' \itemize{
#'   \item Temporal alignment: nrow(Y) matches total scans in sampling_frame
#'   \item Collinearity: Design matrix condition number < 30
#'   \item Compatibility: event_model and baseline_model use same sampling_frame
#' }
#'
#' @examples
#' \dontrun{
#' library(fmridesign)
#' library(fmrihrf)
#'
#' # 1. Define temporal structure
#' sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)
#'
#' # 2. Create trial data (run-relative onsets)
#' trials <- data.frame(
#'   onset = c(10, 30, 50, 70, 90, 110,   # Run 1
#'             10, 30, 50, 70, 90, 110),   # Run 2
#'   run = rep(1:2, each = 6)
#' )
#'
#' # 3. Build event model with trialwise
#' emod <- event_model(
#'   onset ~ trialwise(basis = "spmg1"),
#'   data = trials,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#'
#' # 4. Build baseline model
#' motion <- list(
#'   matrix(rnorm(150 * 6), 150, 6),  # Run 1: 6 motion params
#'   matrix(rnorm(150 * 6), 150, 6)   # Run 2: 6 motion params
#' )
#' bmodel <- baseline_model(
#'   basis = "bs",
#'   degree = 5,
#'   sframe = sframe,
#'   nuisance_list = motion
#' )
#'
#' # 5. Run LSS
#' Y <- matrix(rnorm(300 * 1000), 300, 1000)  # 300 scans, 1000 voxels
#' beta <- lss_design(Y, emod, bmodel, method = "oasis")
#'
#' # Output: 12 × 1000 (12 trials, 1000 voxels)
#' dim(beta)
#' }
#'
#' @seealso
#' \code{\link{lss}} for the traditional matrix-based interface,
#' \code{fmridesign::event_model} for event model creation,
#' \code{fmridesign::baseline_model} for baseline model creation
#'
#' @export
lss_design <- function(Y,
                       event_model,
                       baseline_model = NULL,
                       method = "oasis",
                       oasis = list(),
                       prewhiten = NULL,
                       blockids = NULL,
                       validate = TRUE,
                       ...) {

  # ---- Input Validation ----

  # Check fmridesign availability
  if (!requireNamespace("fmridesign", quietly = TRUE)) {
    stop("Package 'fmridesign' is required for lss_design().\n",
         "Install with: remotes::install_github('bbuchsbaum/fmridesign')",
         call. = FALSE)
  }

  # Check event_model
  if (missing(event_model) || is.null(event_model)) {
    stop("event_model is required. Create with fmridesign::event_model()",
         call. = FALSE)
  }
  if (!inherits(event_model, "event_model")) {
    stop("event_model must be an 'event_model' object from fmridesign",
         call. = FALSE)
  }

  # Check baseline_model if provided
  if (!is.null(baseline_model) && !inherits(baseline_model, "baseline_model")) {
    stop("baseline_model must be a 'baseline_model' object from fmridesign",
         call. = FALSE)
  }

  # Currently only OASIS method supported
  method <- match.arg(method, "oasis")

  # ---- Extract Components ----

  # Get sampling frame
  sframe <- event_model$sampling_frame

  # Validate baseline_model uses same sampling frame
  if (!is.null(baseline_model)) {
    if (!identical(baseline_model$sampling_frame, sframe)) {
      stop("event_model and baseline_model must use the same sampling_frame",
           call. = FALSE)
    }
  }

  # Extract blockids
  if (is.null(blockids)) {
    blockids <- event_model$blockids
  }

  # ---- Build Design Matrices ----

  # Extract X from event_model
  # Use design_matrix() on the whole event_model (returns pre-computed matrix)
  X <- as.matrix(fmridesign::design_matrix(event_model))

  # Extract Z and Nuisance from baseline_model
  Z <- NULL
  Nuisance <- NULL

  if (!is.null(baseline_model)) {
    # Get term matrices using term_matrices() for cleaner extraction
    # term_matrices() returns a named list where keys are guaranteed to be:
    # "drift" (always), "block" (conditional), "nuisance" (conditional)
    bm_term_mats <- fmridesign::term_matrices(baseline_model)

    # Combine drift and block terms into Z
    # These are fixed effect regressors (scanner drift and run intercepts)
    z_terms <- list()
    if ("drift" %in% names(bm_term_mats)) {
      z_terms$drift <- bm_term_mats$drift
    }
    if ("block" %in% names(bm_term_mats)) {
      z_terms$block <- bm_term_mats$block
    }

    if (length(z_terms) > 0) {
      Z <- as.matrix(do.call(cbind, z_terms))
    }

    # Extract nuisance term (motion, physiology, etc.)
    # These are confounds to be projected out
    if ("nuisance" %in% names(bm_term_mats)) {
      Nuisance <- as.matrix(bm_term_mats$nuisance)
    }
  } else {
    # No baseline model provided: inject per-run intercepts into Z
    bl <- fmrihrf::blocklens(sframe)
    expected_scans <- sum(bl)
    runs <- if (!is.null(blockids) && length(blockids) == expected_scans) {
      as.integer(blockids)
    } else {
      rep(seq_along(bl), bl)
    }
    # If only a single run level, just add a single intercept column
    n_levels <- length(unique(runs))
    if (n_levels <= 1L) {
      Z <- matrix(1, nrow = length(runs), ncol = 1)
      colnames(Z) <- "Intercept_run1"
    } else {
      Z <- stats::model.matrix(~ 0 + factor(runs))
      colnames(Z) <- paste0("Intercept_run", seq_len(ncol(Z)))
    }
  }

  # ---- Validation ----

  if (validate) {
    # Basic design compatibility checks
    .validate_design_models(event_model, baseline_model)
    # Check temporal alignment
    expected_scans <- sum(fmrihrf::blocklens(sframe))
    if (nrow(Y) != expected_scans) {
      stop(sprintf("Y has %d rows but sampling_frame expects %d scans",
                   nrow(Y), expected_scans),
           call. = FALSE)
    }

    if (nrow(X) != expected_scans) {
      stop(sprintf("Design matrix has %d rows but sampling_frame expects %d scans",
                   nrow(X), expected_scans),
           call. = FALSE)
    }

    # Check collinearity (combined design)
    X_full <- X
    if (!is.null(Z)) X_full <- cbind(X_full, Z)
    if (!is.null(Nuisance)) X_full <- cbind(X_full, Nuisance)

    cond_num <- tryCatch(kappa(X_full), error = function(e) NA_real_)
    if (is.finite(cond_num) && cond_num > 30) {
      warning(sprintf(
        "High collinearity detected (condition number = %.1f). Consider ridge via oasis$ridge_*",
        cond_num
      ), call. = FALSE)
    }

    # Detect multi-basis and extract K from event_model/X
    K <- 1L
    ntr <- tryCatch(length(event_model$blockids), error = function(e) 0L)
    if (!is.null(ntr) && ntr > 0L) {
      N <- ncol(X)
      if (!is.null(N) && N %% ntr == 0L) {
        Kcand <- as.integer(N / ntr)
        if (Kcand > 1L && Kcand <= 12L) K <- Kcand
      }
    }
    if (isTRUE(K == 1L)) {
      K <- tryCatch(.detect_basis_dimension(X), error = function(e) 1L)
    }
    if (!is.null(K) && K > 1 && is.null(oasis$K)) {
      message(sprintf("Detected multi-basis HRF with K = %d", K))
      oasis$K <- K
    }
  }

  # ---- Call lss() ----

  result <- lss(
    Y = Y,
    X = X,
    Z = Z,
    Nuisance = Nuisance,
    method = method,
    oasis = oasis,
    prewhiten = prewhiten,
    ...
  )

  # ---- Attach Metadata ----
  if (isTRUE(validate)) {
    attr(result, "event_model") <- event_model
    attr(result, "baseline_model") <- baseline_model
    attr(result, "sampling_frame") <- sframe
    attr(result, "method") <- "lss_design"
  }

  return(result)
}

 


#' Validate design compatibility
#' @keywords internal
#' @noRd
.validate_design_models <- function(event_model, baseline_model = NULL) {
  issues <- character(0)

  # Check sampling frame consistency
  if (!is.null(baseline_model)) {
    if (!identical(event_model$sampling_frame, baseline_model$sampling_frame)) {
      issues <- c(issues, "event_model and baseline_model have different sampling_frames")
    }
  }

  # Check for empty design
  dm <- fmridesign::design_matrix(event_model)
  if (ncol(dm) == 0) {
    issues <- c(issues, "event_model produces empty design matrix")
  }

  if (length(issues) > 0) {
    stop("Design validation failed:\n  ", paste(issues, collapse = "\n  "),
         call. = FALSE)
  }

  invisible(TRUE)
}

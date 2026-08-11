#' LSS Analysis with fmridesign Objects
#'
#' Perform Least Squares Separate (LSS) analysis using event_model and
#' baseline_model objects from the fmridesign package. This provides a
#' streamlined interface for complex designs with multi-condition, parametric
#' modulators, and structured nuisance handling.
#'
#' @param Y Numeric matrix of fMRI data (timepoints × voxels).
#' @param event_model An event_model object from \code{fmridesign::event_model()}.
#'   It must contain exactly one \code{trialwise()} target term. Any additional
#'   condition-level event terms are included as common fixed regressors.
#' @param baseline_model Optional baseline_model object from
#'   \code{fmridesign::baseline_model()}. Defines drift correction, block intercepts,
#'   and nuisance regressors. If NULL, basic baseline intercepts are
#'   auto-injected: per-run intercepts derived from `blockids` (or the
#'   sampling frame) are used to ensure proper baseline modeling.
#' @param method LSS method to use. Currently only "oasis" is supported for
#'   event_model integration.
#' @param oasis List of OASIS-specific options: ridge regularization
#'   (\code{ridge_x}, \code{ridge_b}, \code{ridge_mode}), standard errors
#'   (\code{return_se}), etc.  See \code{\link{oasis_options}} and the
#'   Details section of \code{\link{lss}} for the full list.
#'   Note: \code{design_spec} must not be supplied when providing event_model
#'   and is rejected at the adapter boundary;
#'   \code{oasis$whiten} is deprecated — use \code{prewhiten} instead.
#' @param prewhiten Optional prewhitening specification as a list (or
#'   \code{NULL} for no whitening).  Controls temporal autocorrelation
#'   correction via the \pkg{fmriAR} package.  Key fields: \code{method}
#'   (\code{"ar"}, \code{"arma"}, \code{"none"}), \code{p} (AR order or
#'   \code{"auto"}), \code{pooling} (\code{"global"}, \code{"voxel"},
#'   \code{"run"}, \code{"parcel"}), and \code{runs}/\code{parcels}.
#'   See \code{\link{prewhiten_options}} and \code{\link{lss}} for full
#'   details and examples.
#' @param blockids Optional complete exact-integer block/run identifiers, one
#'   per scan. If NULL, run intercepts are derived from the sampling frame.
#' @param validate Logical. If TRUE (default), performs validation checks on
#'   design compatibility, collinearity, and temporal alignment.
#' @param ... Additional arguments passed to the underlying LSS method.
#'
#' @return Normally a trial-by-voxel beta matrix, or a
#'   (trial × basis)-by-voxel matrix for multi-basis HRFs. When OASIS
#'   `return_diag = TRUE` or `return_se = TRUE`, returns
#'   `list(beta, diag?, se?)` with the matrix/list shapes documented in [lss()].
#'   Adapter metadata are attached to the returned object; multi-basis beta and
#'   SE matrices retain the canonical `trial_basis_map`.
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
#' Non-trial event terms are supported as common fixed regressors, for example
#' a parametric modulator alongside the unique trialwise target term:
#' \preformatted{
#'   emod <- event_model(onset ~ trialwise(basis = "spmg1") + hrf(RT),
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
#' Event onsets should be run-relative (resetting to 0 each run) as per
#' fmridesign convention; conversion to global time and run boundaries in the
#' design are handled automatically. A joint call still uses one OASIS
#' other-trial aggregate across all runs, so its coefficients are cross-run
#' pooled. Fit each run separately when run-local LSS coefficients are required.
#'
#' \strong{Prewhitening:}
#'
#' Use the \code{prewhiten} parameter (not the \code{oasis} list) for temporal
#' autocorrelation correction. For active multi-run whitening, this adapter
#' infers the scan-level run segmentation from the sampling frame. An explicit
#' \code{prewhiten$runs} vector is allowed only when it encodes those same
#' boundaries. The event model's \code{blockids} usually has one value per event
#' and is not the required scan-level vector.
#' See \code{\link{lss}} and \code{\link{prewhiten_options}} for full details.
#'
#' \strong{Validation:}
#'
#' When \code{validate = TRUE}, the function checks:
#' \itemize{
#'   \item Temporal alignment: nrow(Y) matches total scans in sampling_frame
#'   \item Collinearity: emits a warning when the full assembled design has a
#'     condition number above 30 and all effective ridge penalties are zero.
#'     The condition number is scale-dependent and is not returned.
#'   \item Compatibility: event_model and baseline_model use same sampling_frame
#' }
#'
#' @examplesIf requireNamespace("fmridesign", quietly = TRUE)
#' \donttest{
#' library(fmridesign)
#' library(fmrihrf)
#'
#' sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)
#'
#' trials <- data.frame(
#'   onset = c(10, 30, 50, 70, 90, 110,
#'             10, 30, 50, 70, 90, 110),
#'   run = rep(1:2, each = 6)
#' )
#'
#' emod <- event_model(
#'   onset ~ trialwise(basis = "spmg1"),
#'   data = trials,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#'
#' motion <- list(
#'   matrix(rnorm(150 * 6), 150, 6),
#'   matrix(rnorm(150 * 6), 150, 6)
#' )
#' bmodel <- baseline_model(
#'   basis = "bs",
#'   degree = 5,
#'   sframe = sframe,
#'   nuisance_list = motion
#' )
#'
#' Y <- matrix(rnorm(300 * 1000), 300, 1000)
#' beta <- lss_design(Y, emod, bmodel, method = "oasis")
#'
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

  validate <- .as_scalar_logical(validate, "validate")
  if (!is.list(oasis)) stop("oasis must be a list", call. = FALSE)
  if (!is.null(oasis$design_spec)) {
    stop("oasis$design_spec must not be supplied to lss_design(); event_model defines the design",
         call. = FALSE)
  }
  if (!is.null(oasis$K)) oasis$K <- .as_positive_integer(oasis$K, "oasis$K")

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
  if (!is.null(blockids)) {
    blockids <- .as_integer_ids(blockids, "blockids")
    frame_lengths <- fmrihrf::blocklens(sframe)
    expected_scans <- sum(frame_lengths)
    if (length(blockids) != expected_scans) {
      stop("blockids must contain one run identifier per scan", call. = FALSE)
    }
    supplied_codes <- match(blockids, unique(blockids))
    frame_codes <- rep(seq_along(frame_lengths), frame_lengths)
    if (!identical(as.integer(supplied_codes), as.integer(frame_codes))) {
      stop("blockids must match the sampling-frame run boundaries", call. = FALSE)
    }
  } else {
    blockids <- event_model$blockids
  }

  if (!is.null(prewhiten)) {
    frame_runs <- rep(
      seq_along(fmrihrf::blocklens(sframe)), fmrihrf::blocklens(sframe)
    )
    if (identical(prewhiten$pooling %||% "global", "run") &&
        is.null(prewhiten$runs)) {
      prewhiten$runs <- frame_runs
    }
    prewhiten <- .resolve_prewhiten_options(prewhiten, internal = FALSE)
    if (prewhiten$method != "none") {
      if (is.null(prewhiten$runs)) {
        prewhiten$runs <- frame_runs
      } else {
        supplied_codes <- match(prewhiten$runs, unique(prewhiten$runs))
        if (length(supplied_codes) != length(frame_runs) ||
            !identical(as.integer(supplied_codes), as.integer(frame_runs))) {
          stop("prewhiten$runs must match the sampling-frame run boundaries",
               call. = FALSE)
        }
      }
    }
  }

  # ---- Build Design Matrices ----

  # Extract the event design with its column metadata. LSS targets only the
  # trialwise term; other event terms are fixed regressors.
  event_dm <- fmridesign::design_matrix(event_model)
  prepared <- .prepare_fmridesign_lss(event_dm, event_model)
  X <- prepared$X
  event_fixed <- prepared$fixed
  if (prepared$K > 1L) {
    if (!is.null(oasis$K) && oasis$K != prepared$K) {
      stop("oasis$K disagrees with the fmridesign basis metadata", call. = FALSE)
    }
    oasis$K <- prepared$K
    oasis$ntrials <- prepared$ntrials
    oasis$trial_basis_map <- prepared$map[, c("column", "trial", "basis", "output_name")]
  }

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

  if (!is.null(event_fixed)) {
    Z <- cbind(Z, event_fixed)
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
    ridge_values <- suppressWarnings(as.numeric(c(
      oasis$ridge_x %||% 0.05,
      oasis$ridge_b %||% 0.05,
      oasis$ridge %||% numeric()
    )))
    ridge_already_set <- any(is.finite(ridge_values) & ridge_values > 0)
    if (is.finite(cond_num) && cond_num > 30 && !ridge_already_set) {
      warning(sprintf(
        "High collinearity detected (condition number = %.1f). Consider ridge via oasis$ridge_*",
        cond_num
      ), call. = FALSE)
    }

    if (prepared$K > 1L) {
      message(sprintf("Detected multi-basis HRF with K = %d", prepared$K))
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
  attr(result, "event_model") <- event_model
  attr(result, "baseline_model") <- baseline_model
  attr(result, "sampling_frame") <- sframe
  attr(result, "method") <- "lss_design"
  attr(result, "trial_basis_map") <- prepared$map

  return(result)
}

# Convert a fmridesign event design into the explicit LSS target contract.
# Trialwise columns are reordered to trial-major/basis-within-trial order;
# all other event terms are fixed regressors.
.prepare_fmridesign_lss <- function(event_dm, event_model) {
  X_all <- as.matrix(event_dm)
  metadata <- attr(event_dm, "col_metadata")
  if (is.null(metadata)) {
    stop("fmridesign design matrix lacks col_metadata; cannot identify LSS targets", call. = FALSE)
  }
  metadata <- as.data.frame(metadata)
  if (!"term_tag" %in% names(metadata) || nrow(metadata) != ncol(X_all)) {
    stop("fmridesign col_metadata is inconsistent with the design matrix", call. = FALSE)
  }

  trial_tags <- names(Filter(function(term) {
    ids <- term$condition_ids
    onsets <- term$onsets
    length(ids) == length(onsets) && length(ids) > 0L && !anyDuplicated(ids)
  }, event_model$terms))
  if (length(trial_tags) != 1L) {
    stop("lss_design() requires exactly one trialwise event term", call. = FALSE)
  }

  target_idx <- which(metadata$term_tag == trial_tags)
  if (!length(target_idx)) {
    stop("No columns were mapped to the trialwise event term", call. = FALSE)
  }
  fixed_idx <- setdiff(seq_len(ncol(X_all)), target_idx)
  target_names <- colnames(X_all)[target_idx]
  if (is.null(target_names) || any(!nzchar(target_names))) {
    stop("Trialwise design columns must have stable names", call. = FALSE)
  }

  has_basis_suffix <- grepl("_b[0-9]+$", target_names)
  if (any(has_basis_suffix) && !all(has_basis_suffix)) {
    stop("Inconsistent basis suffixes in fmridesign trial columns", call. = FALSE)
  }
  trial_term <- event_model$terms[[trial_tags]]
  expected_K <- tryCatch(
    fmrihrf::nbasis(attr(trial_term, "hrfspec")$hrf),
    error = function(e) 1L
  )
  if (expected_K > 1L && !all(has_basis_suffix)) {
    stop(
      "Multi-basis fmridesign columns lack basis identity in their names; ",
      "cannot construct a safe trial/basis map", call. = FALSE
    )
  }
  if (all(has_basis_suffix)) {
    basis <- as.integer(sub("^.*_b([0-9]+)$", "\\1", target_names))
    trial <- sub("_b[0-9]+$", "", target_names)
  } else {
    basis <- rep.int(1L, length(target_names))
    trial <- target_names
  }
  condition_levels <- as.character(trial_term$condition_levels)
  expected_trial_names <- paste0(trial_tags, "_", condition_levels)
  if (length(condition_levels)) {
    if (!setequal(trial, expected_trial_names)) {
      stop(
        "fmridesign trial names disagree with the event-term condition metadata",
        call. = FALSE
      )
    }
    trial_levels <- expected_trial_names
  } else {
    trial_levels <- sort(unique(trial))
  }
  K <- max(basis)
  if (K != expected_K) {
    stop("fmridesign basis names disagree with the event-term HRF metadata", call. = FALSE)
  }
  expected <- data.frame(
    trial = rep(trial_levels, each = K),
    basis = rep(seq_len(K), times = length(trial_levels)),
    stringsAsFactors = FALSE
  )
  key <- paste(trial, basis, sep = "\r")
  expected_key <- paste(expected$trial, expected$basis, sep = "\r")
  order_idx <- match(expected_key, key)
  if (anyNA(order_idx) || anyDuplicated(key) || length(order_idx) != length(target_idx)) {
    stop("Trial/basis columns do not form a complete rectangular mapping", call. = FALSE)
  }

  X <- X_all[, target_idx[order_idx], drop = FALSE]
  row_names <- if (K == 1L) expected$trial else {
    paste(expected$trial, sprintf("basis_%d", expected$basis), sep = ":")
  }
  colnames(X) <- row_names
  fixed <- if (length(fixed_idx)) X_all[, fixed_idx, drop = FALSE] else NULL
  map <- data.frame(
    output_row = seq_along(row_names),
    column = row_names,
    trial = rep(seq_along(trial_levels), each = K),
    trial_name = expected$trial,
    basis = expected$basis,
    source_column = target_idx[order_idx],
    name = row_names,
    output_name = row_names,
    stringsAsFactors = FALSE
  )
  list(X = X, fixed = fixed, K = K, ntrials = length(trial_levels), map = map)
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

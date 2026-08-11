# Internal utilities used across backends.
# Keep this file early in collation order (alphabetical) so helpers are available.

#' Coerce supported inputs to a base matrix
#'
#' Accepts base matrices, data.frames, and Matrix objects.
#' Returns NULL for NULL input.
#'
#' @keywords internal
#' @noRd
.as_base_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "Matrix")) return(as.matrix(x))
  if (is.data.frame(x)) return(as.matrix(x))
  x
}

#' Validate and coerce a positive scalar integer
#'
#' @keywords internal
#' @noRd
.as_positive_integer <- function(x, name) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= 1 && x == floor(x) && x <= .Machine$integer.max
  if (!valid) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  as.integer(x)
}

#' Validate and coerce a non-negative scalar integer
#' @keywords internal
#' @noRd
.as_nonnegative_integer <- function(x, name) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= 0 && x == floor(x) && x <= .Machine$integer.max
  if (!valid) {
    stop(name, " must be a non-negative integer", call. = FALSE)
  }
  as.integer(x)
}

#' Validate a scalar logical without coercion
#' @keywords internal
#' @noRd
.as_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

#' Reject unknown named list fields
#' @keywords internal
#' @noRd
.validate_option_names <- function(x, allowed, name) {
  if (!is.list(x)) stop(name, " must be a list", call. = FALSE)
  nms <- names(x)
  if (length(x) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))) {
    stop(name, " must be a fully named list", call. = FALSE)
  }
  unknown <- setdiff(nms, allowed)
  if (length(unknown)) {
    stop(name, " contains unknown option(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  invisible(x)
}

#' Certified OASIS option names
#' @keywords internal
#' @noRd
.oasis_option_names <- function() {
  c(
    "design_spec", "K", "ntrials", "trial_basis_map", "ridge_mode",
    "ridge_x", "ridge_b", "block_cols", "return_se", "return_diag",
    "add_intercept", "hrf_mode", "infer_K_from_X", "lambda_shape",
    "mu_rough", "ref_hrf", "shrink_global", "orient_ref", "whiten"
  )
}

#' Certified prewhitening option names
#' @keywords internal
#' @noRd
.prewhiten_option_names <- function(internal = FALSE) {
  out <- c(
    "method", "p", "q", "p_max", "pooling", "runs", "parcels",
    "exact_first", "compute_residuals"
  )
  if (internal) c(out, ".whiten_plan") else out
}

#' Validate supplied identities or generate canonical names
#' @keywords internal
#' @noRd
.validate_or_default_names <- function(nms, n, prefix, name) {
  if (is.null(nms)) return(sprintf("%s%d", prefix, seq_len(n)))
  valid <- length(nms) == n && !anyNA(nms) && all(nzchar(nms)) &&
    !anyDuplicated(nms)
  if (!valid) {
    stop(name, " must be complete and unique when supplied", call. = FALSE)
  }
  nms
}

#' Validate an integer identity vector without lossy coercion
#' @keywords internal
#' @noRd
.as_integer_ids <- function(x, name) {
  valid <- is.numeric(x) && length(x) > 0L && !anyNA(x) &&
    all(is.finite(x)) && all(x == floor(x)) &&
    all(abs(x) <= .Machine$integer.max)
  if (!valid) {
    stop(name, " must contain complete exact integer identifiers", call. = FALSE)
  }
  as.integer(x)
}

#' Validate and coerce a non-negative numeric scalar
#'
#' @keywords internal
#' @noRd
.as_nonnegative_scalar <- function(x, name) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 0
  if (!valid) {
    stop(name, " must be a non-negative scalar", call. = FALSE)
  }
  as.numeric(x)
}

#' Null-coalescing helper
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Default trial names: Trial_1, Trial_2, ...
#' @keywords internal
#' @noRd
.default_trial_names <- function(n) {
  sprintf("Trial_%d", seq_len(n))
}

#' Trial names from matrix columns, falling back to Trial_*
#' @keywords internal
#' @noRd
.trial_names_from_cols <- function(X, n = NULL) {
  if (is.null(n)) n <- ncol(X)
  if (!is.null(X) && !is.null(colnames(X))) return(colnames(X))
  .default_trial_names(n)
}

#' Set dimnames for a beta matrix consistently
#' @keywords internal
#' @noRd
.set_beta_dimnames <- function(beta, trial_names = NULL, voxel_names = NULL) {
  if (is.null(trial_names)) trial_names <- .default_trial_names(nrow(beta))
  if (!is.null(trial_names) && is.null(rownames(beta))) rownames(beta) <- trial_names
  if (!is.null(voxel_names) && is.null(colnames(beta))) colnames(beta) <- voxel_names
  beta
}

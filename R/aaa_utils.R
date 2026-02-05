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

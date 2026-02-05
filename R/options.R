#' Option constructors for nested interfaces
#'
#' These helpers create validated option lists for `lss()` and friends.
#'
#' @name fmrilss_options
NULL

#' Construct OASIS options
#'
#' Convenience constructor for the `oasis=` list accepted by `lss(method="oasis")`.
#' Unknown fields are allowed via `...` for forward compatibility.
#'
#' @param design_spec Optional design spec list used to build `X` via `fmrihrf`.
#' @param K Optional basis dimension override.
#' @param ridge_mode `"absolute"` (default) or `"fractional"`.
#' @param ridge_x,ridge_b Non-negative ridge penalties.
#' @param block_cols Voxel block size for blocked products.
#' @param return_se Logical; return standard errors.
#' @param return_diag Logical; return diagnostics.
#' @param add_intercept Logical; add intercept when `Z` is NULL.
#' @param hrf_mode Optional mode (e.g. `"voxhrf"`); advanced use.
#' @param ... Additional options.
#'
#' @return A list with class `"fmrilss_oasis_options"`.
#' @export
oasis_options <- function(
  design_spec = NULL,
  K = NULL,
  ridge_mode = c("absolute", "fractional"),
  ridge_x = 0,
  ridge_b = 0,
  block_cols = 4096L,
  return_se = FALSE,
  return_diag = FALSE,
  add_intercept = TRUE,
  hrf_mode = NULL,
  ...
) {
  ridge_mode <- match.arg(ridge_mode)
  if (!is.null(K)) K <- as.integer(K)
  if (!is.numeric(ridge_x) || length(ridge_x) != 1L || ridge_x < 0) stop("ridge_x must be a non-negative scalar")
  if (!is.numeric(ridge_b) || length(ridge_b) != 1L || ridge_b < 0) stop("ridge_b must be a non-negative scalar")
  block_cols <- as.integer(block_cols)
  if (!is.finite(block_cols) || block_cols < 1L) stop("block_cols must be a positive integer")

  opts <- list(
    design_spec = design_spec,
    K = K,
    ridge_mode = ridge_mode,
    ridge_x = as.numeric(ridge_x),
    ridge_b = as.numeric(ridge_b),
    block_cols = block_cols,
    return_se = isTRUE(return_se),
    return_diag = isTRUE(return_diag),
    add_intercept = isTRUE(add_intercept),
    hrf_mode = hrf_mode
  )

  extra <- list(...)
  if (length(extra)) opts <- utils::modifyList(opts, extra)
  class(opts) <- c("fmrilss_oasis_options", "list")
  opts
}

#' Construct prewhitening options
#'
#' Convenience constructor for the `prewhiten=` list accepted by `lss()`.
#'
#' @param method `"none"`, `"ar"`, or `"arma"`.
#' @param p AR order or `"auto"`.
#' @param q MA order for ARMA.
#' @param p_max Maximum AR order when `p="auto"`.
#' @param pooling `"global"`, `"voxel"`, `"run"`, or `"parcel"`.
#' @param runs Optional run identifiers.
#' @param parcels Optional parcel ids.
#' @param exact_first `"ar1"` or `"none"`.
#' @param compute_residuals Logical.
#'
#' @return A list with class `"fmrilss_prewhiten_options"`.
#' @export
prewhiten_options <- function(
  method = c("none", "ar", "arma"),
  p = "auto",
  q = 0L,
  p_max = 6L,
  pooling = c("global", "voxel", "run", "parcel"),
  runs = NULL,
  parcels = NULL,
  exact_first = c("ar1", "none"),
  compute_residuals = TRUE
) {
  method <- match.arg(method)
  pooling <- match.arg(pooling)
  exact_first <- match.arg(exact_first)

  opts <- list(
    method = method,
    p = p,
    q = as.integer(q),
    p_max = as.integer(p_max),
    pooling = pooling,
    runs = runs,
    parcels = parcels,
    exact_first = exact_first,
    compute_residuals = isTRUE(compute_residuals)
  )
  class(opts) <- c("fmrilss_prewhiten_options", "list")
  opts
}

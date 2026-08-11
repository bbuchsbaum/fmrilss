#' Option constructors for nested interfaces
#'
#' These helpers create validated option lists for `lss()` and friends.
#'
#' @name fmrilss_options
#' @return No value itself. This topic groups the documented constructors
#'   `stglmnet_options()`, `oasis_options()`, and `prewhiten_options()`.
#'
#' @examples
#' stglmnet_options(mode = "fixed", lambda = 0.1)
#' oasis_options(ridge_x = 0.1, ridge_b = 0.1)
#' prewhiten_options(method = "ar", p = 1)
NULL

#' Construct stglmnet backend options
#'
#' Convenience constructor for the `stglmnet=` list accepted by
#' `lss(method = "stglmnet")`.
#' Advanced fields accepted through `...` are validated; unknown fields are
#' rejected so misspelled scientific controls cannot be silently ignored.
#'
#' @param mode `"cv"` (default) selects lambda by internal cross-validation,
#'   while `"fixed"` uses the supplied lambda sequence or the smallest fitted
#'   lambda when no scalar is provided.
#' @param alpha Elastic-net mixing parameter passed to `glmnet`.
#' @param lambda Optional lambda sequence (or scalar in fixed mode).
#' @param overlap_strategy Trial-overlap penalty mapping. One of `"none"`,
#'   `"multiplicative"`, `"additive"`, `"hybrid"`, or `"threshold"`.
#' @param pool_to_mean Logical; reparameterize trial effects into a pooled mean
#'   plus orthogonal contrasts.
#' @param pool_strength Penalty multiplier applied to pooled contrasts.
#' @param pool_mean_penalty Penalty applied to the pooled mean coefficient.
#' @param whiten One of `"inherit"` (default), `"auto"`, `"never"`, or
#'   `"always"`. `"inherit"` uses the top-level `prewhiten=` argument only.
#' @param cv_folds Number of folds used when `mode = "cv"`.
#' @param cv_type.measure Cross-validation objective.
#' @param cv_select Lambda selection rule in CV mode. `"optimal"` uses the
#'   best-scoring lambda, `"1se"` applies the one-standard-error rule.
#' @param return_fit Logical; when `TRUE`, `lss(method="stglmnet")` returns a
#'   list containing `beta`, fit metadata, and the selected lambda.
#' @param ... Certified advanced backend fields such as graph-pooling,
#'   overlap-strength, fold, and whitening controls. Unknown names are rejected.
#'
#' @return A list with class `"fmrilss_stglmnet_options"`.
#' @examples
#' stglmnet_options(mode = "fixed", lambda = 0.05, alpha = 0.5)
#' @export
stglmnet_options <- function(
  mode = c("cv", "fixed"),
  alpha = 0.2,
  lambda = NULL,
  overlap_strategy = c("none", "multiplicative", "additive", "hybrid", "threshold"),
  pool_to_mean = FALSE,
  pool_strength = 1,
  pool_mean_penalty = 0,
  whiten = c("inherit", "auto", "never", "always"),
  cv_folds = 5L,
  cv_type.measure = c("auto", "mse", "correlation", "reliability", "composite"),
  cv_select = c("optimal", "1se"),
  return_fit = FALSE,
  ...
) {
  mode <- match.arg(mode)
  overlap_strategy <- match.arg(overlap_strategy)
  whiten <- match.arg(whiten)
  cv_type.measure <- match.arg(cv_type.measure)
  cv_select <- match.arg(cv_select)

  pool_to_mean <- .as_scalar_logical(pool_to_mean, "pool_to_mean")
  return_fit <- .as_scalar_logical(return_fit, "return_fit")
  cv_folds <- .as_positive_integer(cv_folds, "cv_folds")

  opts <- list(
    mode = mode,
    alpha = as.numeric(alpha),
    lambda = lambda,
    overlap_strategy = overlap_strategy,
    pool_to_mean = pool_to_mean,
    pool_strength = as.numeric(pool_strength),
    pool_mean_penalty = as.numeric(pool_mean_penalty),
    whiten = whiten,
    cv_folds = cv_folds,
    cv_type.measure = cv_type.measure,
    cv_select = cv_select,
    return_fit = return_fit
  )

  extra <- list(...)
  if (length(extra)) opts <- utils::modifyList(opts, extra)
  opts <- .stg_resolve_options(opts)
  class(opts) <- c("fmrilss_stglmnet_options", "list")
  opts
}

#' Construct OASIS options
#'
#' Convenience constructor for the `oasis=` list accepted by `lss(method="oasis")`.
#' Unknown fields are rejected so misspelled scientific controls cannot be
#' silently ignored.
#'
#' @param design_spec Optional design spec list used to build `X` via `fmrihrf`.
#' @param K Optional basis dimension override.
#' @param ntrials Number of trials for a raw multi-basis design.
#' @param trial_basis_map For a raw multi-basis design, a data frame with one
#'   row per `X` column and fields `column`, `trial`, and `basis`.
#' @param ridge_mode `"fractional"` (default) or `"absolute"`.
#' @param ridge_x,ridge_b Non-negative ridge penalties (defaults 0.05 each).
#' @param block_cols Voxel block size for blocked products.
#' @param return_se Logical; return standard errors.
#' @param return_diag Logical; return diagnostics.
#' @param add_intercept Logical; add intercept when `Z` is NULL.
#' @param hrf_mode Optional mode (e.g. `"voxhrf"`); advanced use.
#' @param ... Certified advanced options: `infer_K_from_X`, `lambda_shape`,
#'   `mu_rough`, `ref_hrf`, `shrink_global`, `orient_ref`, and the deprecated
#'   `whiten` compatibility field. Unknown names are rejected.
#'
#' @return A list with class `"fmrilss_oasis_options"`.
#' @examples
#' oasis_options(ridge_mode = "fractional", ridge_x = 0.1, ridge_b = 0.1)
#' @export
oasis_options <- function(
  design_spec = NULL,
  K = NULL,
  ntrials = NULL,
  trial_basis_map = NULL,
  ridge_mode = c("fractional", "absolute"),
  ridge_x = 0.05,
  ridge_b = 0.05,
  block_cols = 4096L,
  return_se = FALSE,
  return_diag = FALSE,
  add_intercept = TRUE,
  hrf_mode = NULL,
  ...
) {
  ridge_mode <- match.arg(ridge_mode)
  if (!is.null(K)) K <- .as_positive_integer(K, "K")
  if (!is.null(ntrials)) ntrials <- .as_positive_integer(ntrials, "ntrials")
  ridge_x <- .as_nonnegative_scalar(ridge_x, "ridge_x")
  ridge_b <- .as_nonnegative_scalar(ridge_b, "ridge_b")
  block_cols <- .as_positive_integer(block_cols, "block_cols")
  return_se <- .as_scalar_logical(return_se, "return_se")
  return_diag <- .as_scalar_logical(return_diag, "return_diag")
  add_intercept <- .as_scalar_logical(add_intercept, "add_intercept")

  if (!is.null(hrf_mode)) {
    hrf_mode <- match.arg(hrf_mode, c("voxel_ridge", "voxhrf"))
  }

  extra <- list(...)
  advanced <- c(
    "infer_K_from_X", "lambda_shape", "mu_rough", "ref_hrf",
    "shrink_global", "orient_ref", "whiten"
  )
  .validate_option_names(extra, advanced, "oasis options supplied through ...")
  if (!is.null(extra$infer_K_from_X)) {
    extra$infer_K_from_X <- .as_scalar_logical(extra$infer_K_from_X, "infer_K_from_X")
  }
  if (!is.null(extra$orient_ref)) {
    extra$orient_ref <- .as_scalar_logical(extra$orient_ref, "orient_ref")
  }
  for (nm in intersect(c("lambda_shape", "mu_rough", "shrink_global"), names(extra))) {
    extra[[nm]] <- .as_nonnegative_scalar(extra[[nm]], nm)
  }
  if (!is.null(extra$shrink_global) && extra$shrink_global > 1) {
    stop("shrink_global must be between 0 and 1", call. = FALSE)
  }

  opts <- list(
    design_spec = design_spec,
    K = K,
    ntrials = ntrials,
    trial_basis_map = trial_basis_map,
    ridge_mode = ridge_mode,
    ridge_x = ridge_x,
    ridge_b = ridge_b,
    block_cols = block_cols,
    return_se = return_se,
    return_diag = return_diag,
    add_intercept = add_intercept,
    hrf_mode = hrf_mode
  )

  if (length(extra)) opts <- utils::modifyList(opts, extra)
  if (isTRUE(opts$return_se) && (opts$ridge_x != 0 || opts$ridge_b != 0)) {
    stop(
      "return_se requires ridge_x = ridge_b = 0; ridge uncertainty is diagnostic-only",
      call. = FALSE
    )
  }
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
#' @param runs Optional complete exact-integer run identifiers. Required for
#'   `pooling = "run"`; execution requires one value per timepoint.
#' @param parcels Optional complete exact-integer parcel identifiers. Required
#'   for `pooling = "parcel"`; execution requires one value per voxel.
#' @param exact_first `"ar1"` or `"none"`.
#' @param compute_residuals Logical.
#'
#' @return A list with class `"fmrilss_prewhiten_options"`.
#' @examples
#' prewhiten_options(method = "ar", p = 1, pooling = "run",
#'                   runs = rep(1:2, each = 50))
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

  if (identical(p, "auto")) {
    p_out <- p
  } else {
    p_out <- .as_positive_integer(p, "p")
  }
  q <- .as_nonnegative_integer(q, "q")
  p_max <- .as_positive_integer(p_max, "p_max")
  compute_residuals <- .as_scalar_logical(compute_residuals, "compute_residuals")
  if (method == "arma" && q == 0L) {
    stop("q must be a positive integer when method = 'arma'", call. = FALSE)
  }
  if (method == "ar" && q != 0L) {
    stop("q must be 0 when method = 'ar'", call. = FALSE)
  }
  if (pooling == "run" && (is.null(runs) || !length(runs) || anyNA(runs))) {
    stop("runs must be supplied and complete when pooling = 'run'", call. = FALSE)
  }
  if (pooling == "parcel" && (is.null(parcels) || !length(parcels) || anyNA(parcels))) {
    stop("parcels must be supplied and complete when pooling = 'parcel'", call. = FALSE)
  }
  if (!is.null(runs)) runs <- .as_integer_ids(runs, "runs")
  if (!is.null(parcels)) parcels <- .as_integer_ids(parcels, "parcels")

  opts <- list(
    method = method,
    p = p_out,
    q = q,
    p_max = p_max,
    pooling = pooling,
    runs = runs,
    parcels = parcels,
    exact_first = exact_first,
    compute_residuals = compute_residuals
  )
  class(opts) <- c("fmrilss_prewhiten_options", "list")
  opts
}

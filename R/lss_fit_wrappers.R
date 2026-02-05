#' Convenience wrappers for modern `lss()` usage
#'
#' These functions provide a modern-signature entry point for methods that
#' historically required a block-design (`bdes`) object.
#'
#' @name lss_fit_wrappers
NULL

#' Naive LSS with modern signature
#'
#' Equivalent to calling `lss(..., method = "naive")`.
#'
#' @param Y Numeric matrix (timepoints x voxels).
#' @param X Trial design matrix (timepoints x trials).
#' @param Z Optional experimental regressors.
#' @param Nuisance Optional nuisance regressors to project out.
#' @param prewhiten Optional prewhitening options list (see `prewhiten_options()`).
#'
#' @return A numeric matrix (trials x voxels) of beta estimates.
#' @export
lss_naive_fit <- function(Y, X, Z = NULL, Nuisance = NULL, prewhiten = NULL) {
  lss(Y = Y, X = X, Z = Z, Nuisance = Nuisance, method = "naive", prewhiten = prewhiten)
}

#' Optimized LSS with modern signature
#'
#' This is a convenience wrapper around `lss()` that selects one of the optimized
#' implementations.
#'
#' @param Y Numeric matrix (timepoints x voxels).
#' @param X Trial design matrix (timepoints x trials).
#' @param Z Optional experimental regressors.
#' @param Nuisance Optional nuisance regressors to project out.
#' @param engine `"cpp"` (default) for `method="cpp_optimized"` or `"r"` for
#'   `method="r_optimized"`.
#' @param block_size Block size used by the C++ optimized path.
#' @param prewhiten Optional prewhitening options list (see `prewhiten_options()`).
#'
#' @return A numeric matrix (trials x voxels) of beta estimates.
#' @export
lss_optimized_fit <- function(Y, X, Z = NULL, Nuisance = NULL,
                              engine = c("cpp", "r"),
                              block_size = 96,
                              prewhiten = NULL) {
  engine <- match.arg(engine)
  method <- if (engine == "cpp") "cpp_optimized" else "r_optimized"
  lss(
    Y = Y, X = X, Z = Z, Nuisance = Nuisance,
    method = method,
    block_size = block_size,
    prewhiten = prewhiten
  )
}


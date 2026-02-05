#' Mixed Model Solver
#'
#' Solves mixed models with random effects using REML or ML estimation.
#' This function provides a unified interface to mixed model estimation,
#' similar to the lss/lsa functions in this package.
#'
#' @param Y Response vector or matrix. If a matrix, each column is treated as 
#'   a separate response variable.
#' @param X Design matrix for fixed effects. If NULL, defaults to intercept only.
#' @param Z Design matrix for random effects. If NULL, defaults to identity matrix.
#' @param K Kinship matrix for random effects. If NULL, defaults to identity matrix.
#' @param Nuisance An alias for X, provided for consistency with lss/lsa interface.
#'   If both X and Nuisance are provided, X takes precedence.
#' @param method Character string specifying the estimation method:
#'   \itemize{
#'     \item "REML" - Restricted Maximum Likelihood (default)
#'     \item "ML" - Maximum Likelihood
#'   }
#' @param bounds Numeric vector of length 2 specifying bounds for variance 
#'   component optimization. Defaults to c(1e-9, 1e9).
#' @param SE Logical, whether to compute and return standard errors. Defaults to FALSE.
#' @param return_Hinv Logical, whether to return the inverse of the H matrix. 
#'   Defaults to FALSE.
#'
#' @return A list containing:
#'   \item{Vu}{Estimated variance component for random effects.}
#'   \item{Ve}{Estimated variance component for residuals.}
#'   \item{beta}{Estimated fixed effects coefficients.}
#'   \item{u}{Estimated random effects coefficients.}
#'   \item{LL}{Log-likelihood of the model.}
#'   \item{beta.SE}{Standard errors of fixed effects coefficients (if SE = TRUE).}
#'   \item{u.SE}{Standard errors of random effects coefficients (if SE = TRUE).}
#'   \item{Hinv}{Inverse of H matrix (if return_Hinv = TRUE).}
#'
#' @details
#' This function fits the mixed model: Y = X*beta + Z*u + error, where
#' u ~ N(0, Vu*K) and error ~ N(0, Ve*I). The variance components Vu and Ve
#' are estimated using REML or ML.
#'
#' @examples
#' \dontrun{
#' # Example with random data
#' set.seed(123)
#' n <- 100
#' Y <- rnorm(n)
#' Z <- matrix(rnorm(n * 5), n, 5)
#' K <- diag(5)
#' X <- matrix(1, n, 1)
#' 
#' # Fit mixed model
#' result <- mixed_solve(Y, X, Z, K)
#' }
#' @export
mixed_solve <- function(Y, X = NULL, Z = NULL, K = NULL, Nuisance = NULL,
                        method = c("REML", "ML"),
                        bounds = c(1e-9, 1e9),
                        SE = FALSE,
                        return_Hinv = FALSE) {
  
  method <- match.arg(method)
  
  # Input validation
  if (!is.vector(Y) && !is.matrix(Y)) {
    stop("Y must be a vector or matrix")
  }
  
  # Handle Nuisance parameter (alias for X)
  if (is.null(X) && !is.null(Nuisance)) {
    X <- Nuisance
  }

  # Guard against saturated / ill-posed fixed-effect designs early, to avoid
  # cryptic C++ errors and provide stable error contracts.
  n_obs <- if (is.matrix(Y)) nrow(Y) else length(Y)
  if (!is.null(X)) {
    X <- as.matrix(X)
    if (nrow(X) != n_obs) {
      stop("X must have the same number of rows as Y", call. = FALSE)
    }
    n_filtered <- if (is.matrix(Y)) sum(stats::complete.cases(Y)) else sum(!is.na(Y))
    if (n_filtered <= ncol(X)) {
      stop("Need more non-NA observations than columns in X", call. = FALSE)
    }
  }
  
  # Use the C++ implementation (currently the only available method)
  result <- .mixed_solve_cpp(Y, X, Z, K, method, bounds, SE, return_Hinv)
  
  return(result)
}

#' Mixed Model Solver using C++
#'
#' C++ implementation of the mixed model solver. This function is typically
#' called through the main `mixed_solve` function rather than directly.
#'
#' @param Y Response vector.
#' @param X Design matrix for fixed effects (default: intercept only).
#' @param Z Design matrix for random effects (default: identity matrix).
#' @param K Kinship matrix (default: identity matrix).
#' @param method Optimization method, either "REML" or "ML".
#' @param bounds Bounds for the optimizer.
#' @param SE Logical, whether to return standard errors.
#' @param return_Hinv Logical, whether to return the inverse of H.
#' @return A list with mixed model results.
#' @keywords internal
.mixed_solve_cpp <- function(Y, X = NULL, Z = NULL, K = NULL,
                             method = "REML",
                             bounds = c(1e-9, 1e9),
                             SE = FALSE,
                             return_Hinv = FALSE) {

  result <- mixed_solve_internal(Y, Z, K, X, method, bounds, SE, return_Hinv)
  return(result)
}

#' @rdname mixed_solve
#' @export
mixed_solve_cpp <- function(Y, X = NULL, Z = NULL, K = NULL, Nuisance = NULL,
                            method = c("REML", "ML"),
                            bounds = c(1e-9, 1e9),
                            SE = FALSE,
                            return_Hinv = FALSE) {
  # Backward compatibility function - calls the main mixed_solve
  mixed_solve(Y, X, Z, K, Nuisance, method, bounds, SE, return_Hinv)
}

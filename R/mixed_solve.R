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
#' set.seed(123)
#' n <- 100
#' Y <- rnorm(n)
#' Z <- matrix(rnorm(n * 5), n, 5)
#' K <- diag(5)
#' X <- matrix(1, n, 1)
#'
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
  if ((!is.vector(Y) && !is.matrix(Y)) || !is.numeric(Y)) {
    stop("Y must be a numeric vector or matrix")
  }
  if (is.matrix(Y) && ncol(Y) < 1L) {
    stop("Y must contain at least one response column", call. = FALSE)
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
    n_filtered <- if (is.matrix(Y)) NULL else sum(is.finite(Y))
    if (!is.null(n_filtered) && n_filtered <= ncol(X)) {
      stop("Need more non-NA observations than columns in X", call. = FALSE)
    }
  }
  if (!is.null(Z) && (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != n_obs)) {
    stop("Z must be a numeric matrix with the same number of rows as Y", call. = FALSE)
  }

  if (is.matrix(Y)) {
    fits <- lapply(seq_len(ncol(Y)), function(j) {
      mixed_solve(
        Y = Y[, j], X = X, Z = Z, K = K, method = method,
        bounds = bounds, SE = SE, return_Hinv = return_Hinv
      )
    })
    response_names <- colnames(Y) %||% paste0("Response_", seq_len(ncol(Y)))
    combine_rows <- function(field) {
      out <- t(vapply(fits, `[[`, numeric(length(fits[[1]][[field]])), field))
      rownames(out) <- response_names
      out
    }
    result <- list(
      Vu = stats::setNames(vapply(fits, `[[`, numeric(1), "Vu"), response_names),
      Ve = stats::setNames(vapply(fits, `[[`, numeric(1), "Ve"), response_names),
      beta = combine_rows("beta"),
      u = combine_rows("u"),
      LL = stats::setNames(vapply(fits, `[[`, numeric(1), "LL"), response_names)
    )
    if (SE) {
      result$beta.SE <- combine_rows("beta.SE")
      result$u.SE <- combine_rows("u.SE")
    }
    if (return_Hinv) {
      result$Hinv <- stats::setNames(lapply(fits, `[[`, "Hinv"), response_names)
    }
    return(result)
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

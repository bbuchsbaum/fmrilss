% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixed_solve.R
\name{mixed_solve}
\alias{mixed_solve}
\alias{mixed_solve_cpp}
\title{Mixed Model Solver}
\usage{
mixed_solve(
  Y,
  X = NULL,
  Z = NULL,
  K = NULL,
  Nuisance = NULL,
  method = c("REML", "ML"),
  bounds = c(1e-09, 1e+09),
  SE = FALSE,
  return_Hinv = FALSE
)

mixed_solve_cpp(
  Y,
  X = NULL,
  Z = NULL,
  K = NULL,
  Nuisance = NULL,
  method = c("REML", "ML"),
  bounds = c(1e-09, 1e+09),
  SE = FALSE,
  return_Hinv = FALSE
)
}
\arguments{
\item{Y}{Response vector or matrix. If a matrix, each column is treated as 
a separate response variable.}

\item{X}{Design matrix for fixed effects. If NULL, defaults to intercept only.}

\item{Z}{Design matrix for random effects. If NULL, defaults to identity matrix.}

\item{K}{Kinship matrix for random effects. If NULL, defaults to identity matrix.}

\item{Nuisance}{An alias for X, provided for consistency with lss/lsa interface.
If both X and Nuisance are provided, X takes precedence.}

\item{method}{Character string specifying the estimation method:
\itemize{
  \item "REML" - Restricted Maximum Likelihood (default)
  \item "ML" - Maximum Likelihood
}}

\item{bounds}{Numeric vector of length 2 specifying bounds for variance 
component optimization. Defaults to c(1e-9, 1e9).}

\item{SE}{Logical, whether to compute and return standard errors. Defaults to FALSE.}

\item{return_Hinv}{Logical, whether to return the inverse of the H matrix. 
Defaults to FALSE.}
}
\value{
A list containing:
  \item{Vu}{Estimated variance component for random effects.}
  \item{Ve}{Estimated variance component for residuals.}
  \item{beta}{Estimated fixed effects coefficients.}
  \item{u}{Estimated random effects coefficients.}
  \item{LL}{Log-likelihood of the model.}
  \item{beta.SE}{Standard errors of fixed effects coefficients (if SE = TRUE).}
  \item{u.SE}{Standard errors of random effects coefficients (if SE = TRUE).}
  \item{Hinv}{Inverse of H matrix (if return_Hinv = TRUE).}
}
\description{
Solves mixed models with random effects using REML or ML estimation.
This function provides a unified interface to mixed model estimation,
similar to the lss/lsa functions in this package.
}
\details{
This function fits the mixed model: Y = X*beta + Z*u + error, where
u ~ N(0, Vu*K) and error ~ N(0, Ve*I). The variance components Vu and Ve
are estimated using REML or ML.
}
\examples{
\dontrun{
# Example with random data
set.seed(123)
n <- 100
Y <- rnorm(n)
Z <- matrix(rnorm(n * 5), n, 5)
K <- diag(5)
X <- matrix(1, n, 1)

# Fit mixed model
result <- mixed_solve(Y, X, Z, K)
}
}

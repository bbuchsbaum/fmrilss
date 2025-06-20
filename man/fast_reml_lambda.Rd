% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fast_reml_lambda}
\alias{fast_reml_lambda}
\title{Fast analytical REML estimation for single variance component}
\arguments{
\item{omega}{Transformed response vector Q'y}

\item{theta}{Transformed eigenvalues}

\item{tol}{Convergence tolerance for Newton iterations}

\item{max_iter}{Maximum Newton iterations}
}
\value{
Estimated variance ratio λ
}
\description{
For a single variance component model, the REML estimate of λ = σe²/σu²
has a closed-form solution that can be computed efficiently.
}

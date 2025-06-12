#' fmrilss: Least Squares Separate (LSS) Analysis for fMRI Data
#'
#' This package implements efficient least squares separate (LSS) analysis for
#' functional magnetic resonance imaging (fMRI) data. LSS is used to estimate
#' trial-by-trial activation patterns in event-related fMRI designs.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{lss}}: Main function for performing LSS analysis
#'   \item \code{\link{lss_naive}}: Naive LSS implementation for reference
#'   \item \code{\link{project_confounds}}: R implementation for projecting out confounds
#'   \item \code{\link{project_confounds_cpp}}: Fast C++ confound projection
#'   \item \code{\link{lss_beta_cpp}}: Vectorized C++ LSS beta computation
#'   \item \code{\link{get_data_matrix}}: Helper function for data extraction
#' }
#'
#' @section Features:
#' \itemize{
#'   \item Optimized C++ implementation using vectorized matrix algebra
#'   \item Memory-efficient projection without forming Q matrices
#'   \item Cholesky decomposition for numerical stability
#'   \item Fallback R implementation with QR decomposition
#'   \item Support for various design matrix configurations
#'   \item Robust numerical handling for edge cases
#'   \item OpenMP support for multi-core processing
#' }
#'
#' @docType package
#' @name fmrilss-package
#' @aliases fmrilss
#' @useDynLib fmrilss, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases median sd var
NULL

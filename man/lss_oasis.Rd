\name{lss_oasis}
\alias{lss_oasis}
\title{OASIS single-trial solver inside \code{fmrilss::lss}}
\description{
Adds \code{method="oasis"} to \code{fmrilss::lss()}, providing an exact, fast closed-form
Least-Squares Separate (LSS) estimator for single-trial betas. Supports both single-basis
(K=1) and multi-basis (K>1) HRFs. Optionally applies fractional ridge on the LSS Gram for 
stability, computes standard errors, and constructs trial-wise designs via \pkg{fmrihrf}, 
including HRF grid searches.
}
\usage{
lss(Y, X = NULL, Z = NULL, Nuisance = NULL,
    method = "oasis",
    oasis = list(), ...)
}
\arguments{
  \item{Y}{T x V matrix of time series data.}
  \item{X}{T x N matrix of trial-wise regressors for a single condition. If \code{NULL},
           provide \code{oasis$design_spec} to build via \pkg{fmrihrf}.}
  \item{Z}{T x K matrix of fixed regressors (e.g., other effects of interest) to be projected out.}
  \item{Nuisance}{T x P matrix of confounds (intercept, motion, drifts, etc.).}
  \item{method}{Must be \code{"oasis"} to use this backend.}
  \item{oasis}{List of options:
    \describe{
      \item{\code{design_spec}}{List with fields:
        \code{sframe} (from fmrihrf::sampling_frame),
        \code{cond} (list with onsets, duration, amplitude, hrf, span),
        \code{others} (list of other conditions),
        \code{method} ("fft"|"conv")}
      \item{\code{K}}{Explicit basis dimension (auto-detected from HRF if not provided)}
      \item{\code{ridge_mode}}{"absolute" (default) or "fractional"}
      \item{\code{ridge_x}, \code{ridge_b}}{Ridge strengths on the LSS Gram (default 0)}
      \item{\code{block_cols}}{Voxel block size for batched products (default 4096)}
      \item{\code{return_se}}{Logical; if TRUE, also return per-trial standard errors}
      \item{\code{return_diag}}{Logical; if TRUE, also return design diagnostics}
      \item{\code{whiten}}{Character; "none" (default) or "ar1"}
    }}
}
\value{
By default, a \eqn{N \times V} matrix of single-trial betas (or \eqn{(N \times K) \times V} 
for multi-basis HRFs).
If \code{return_se} or \code{return_diag} is TRUE, a list with elements:
\code{$beta}, \code{$se}, and/or \code{$diag}.
}
\details{
OASIS computes exact LSS betas without per-trial GLMs. After projecting out nuisances once,
it uses a closed-form two-regressor solution for each trial (K=1), or a small block solve
for multi-basis HRFs (K>1), batched across voxels.

\subsection{Ridge Regularization}{
Setting \code{ridge_x} or \code{ridge_b} > 0 applies ridge to stabilize highly collinear
designs. When \code{ridge_mode="fractional"}, ridge values are scaled by mean design energy:
\itemize{
  \item For K=1: \eqn{\lambda_x = \rho_x \times \text{mean}(||a_j||^2)}
  \item For K>1: \eqn{\lambda_x = \rho_x \times \text{mean}(\text{tr}(D_j)/K)}
}
}

\subsection{Multi-basis HRF Support}{
Automatically detects K from \code{fmrihrf::nbasis(hrf)} and switches between:
\itemize{
  \item K=1: Fast 2x2 closed form (SPMG1, single LWU, Gamma)
  \item K>1: Block solver (SPMG3, FIR, B-spline, Fourier)
}
}
}
\examples{
\dontrun{
# Single-basis HRF (SPMG1)
library(fmrihrf)
sframe <- sampling_frame(blocklens = c(150,150), TR = 1.0)
beta <- lss(Y, X = NULL, method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(onsets = c(10, 22, 35, 47, 65, 88),
                  hrf = make_hrf("spmg1"))
    ),
    ridge_mode = "fractional",
    ridge_x = 0.05, ridge_b = 0.05
  ))

# Multi-basis HRF (SPMG3, K=3)
beta <- lss(Y, X = NULL, method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(onsets = c(10, 22, 35, 47, 65, 88),
                  hrf = HRF_SPMG3)
    )
  ))

# FIR basis (K=12)
beta <- lss(Y, X = NULL, method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(onsets = c(10, 22, 35, 47, 65, 88),
                  hrf = hrf_fir_generator(nbasis = 12, span = 24))
    )
  ))
}
}
\seealso{
\code{\link[fmrihrf]{regressor_set}}, \code{\link[fmrihrf]{evaluate}}, 
\code{\link[fmrihrf]{sampling_frame}}, \code{\link{lss}}
}
\author{The \pkg{fmrilss} authors}
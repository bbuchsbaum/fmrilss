\name{lss_naive}
\alias{lss_naive}
\title{Naive Least Squares Separate (LSS) Analysis}
\usage{
lss_naive(Y = NULL, bdes, dset = NULL)
}
\arguments{
\item{Y}{A numeric matrix where rows are timepoints and columns are voxels/features.
If NULL, the function will attempt to extract data from \code{dset}.}

\item{bdes}{A list containing design matrices with components:
\itemize{
\item \code{dmat_base}: Base design matrix (e.g., intercept, drift terms)
\item \code{dmat_fixed}: Fixed effects design matrix (optional)
\item \code{dmat_ran}: Random/trial design matrix for LSS analysis
\item \code{fixed_ind}: Indices for fixed effects (optional)
}}

\item{dset}{Optional dataset object. If provided and Y is NULL, data will be
extracted using \code{get_data_matrix}.}
}
\value{
A numeric matrix with dimensions (n_events x n_voxels) containing
the LSS beta estimates for each trial and voxel.
}
\description{
Performs LSS analysis using the naive approach where each trial model is fit
separately. This is the conceptually simplest implementation but less efficient
than the optimized \code{\link{lss}} function.
}
\details{
This function implements the naive LSS approach where for each trial, a separate
GLM is fitted that includes:
\itemize{
\item All base regressors (intercept, drift, etc.)
\item All fixed effects regressors (if any)
\item Only the current trial's regressor from the trial design matrix
}

While less efficient than the optimized \code{\link{lss}} function, this
implementation is conceptually simpler and can serve as a reference or for
validation purposes.
}
\examples{
\dontrun{
# Using same setup as lss() examples
beta_estimates_naive <- lss_naive(Y = Y, bdes = bdes)

# Compare with optimized version
beta_estimates_fast <- lss(Y = Y, bdes = bdes)
max(abs(beta_estimates_naive - beta_estimates_fast))
}
}
\seealso{
\code{\link{lss}} for the optimized implementation
} 
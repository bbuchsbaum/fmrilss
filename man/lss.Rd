\name{lss}
\alias{lss}
\title{Least Squares Separate (LSS) Analysis}
\usage{
lss(Y = NULL, bdes, dset = NULL, use_cpp = TRUE)
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

\item{use_cpp}{Logical. If TRUE (default), uses C++ implementation for speed.
If FALSE, uses pure R implementation.}
}
\value{
A numeric matrix with dimensions (n_events x n_voxels) containing
the LSS beta estimates for each trial and voxel.
}
\description{
Performs least squares separate (LSS) analysis for fMRI data to estimate
trial-by-trial activation patterns in event-related designs.
}
\details{
LSS (Least Squares Separate) is a method for estimating single-trial activation
patterns in event-related fMRI designs. For each trial, it estimates the activation
while treating all other trials as confounds.

The function supports both R and C++ implementations. The C++ version is
recommended for large datasets due to computational efficiency.
}
\examples{
\dontrun{
# Create example design matrices
n_timepoints <- 100
n_trials <- 20
n_voxels <- 1000

# Base design (intercept + linear trend)
dmat_base <- cbind(1, 1:n_timepoints)

# Trial design matrix (one column per trial)
dmat_ran <- matrix(0, n_timepoints, n_trials)
for(i in 1:n_trials) {
  trial_onset <- sample(1:(n_timepoints-10), 1)
  dmat_ran[trial_onset:(trial_onset+5), i] <- 1
}

# Create design list
bdes <- list(
  dmat_base = dmat_base,
  dmat_ran = dmat_ran,
  fixed_ind = NULL
)

# Simulate data
Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

# Run LSS analysis
beta_estimates <- lss(Y = Y, bdes = bdes)
}
}
\references{
Mumford, J.A., et al. (2012). Deconvolving BOLD activation in event-related
designs for multivoxel pattern classification analyses. NeuroImage, 59(3), 2636-2643.
} 
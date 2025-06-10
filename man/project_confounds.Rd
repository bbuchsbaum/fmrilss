\name{project_confounds}
\alias{project_confounds}
\title{Project Out Confound Variables}
\usage{
project_confounds(X)
}
\arguments{
\item{X}{Confound design matrix (n x p) where n is number of timepoints
and p is number of confound regressors}
}
\value{
Projection matrix Q (n x n) that projects out the column space of X
}
\description{
Computes the orthogonal projection matrix Q = I - X(X'X)^(-1)X' that projects
out the space spanned by confound regressors X. This is useful for advanced
users who want to cache and reuse projection matrices.
}
\details{
This function uses QR decomposition for numerical stability instead of
computing the Moore-Penrose pseudoinverse directly. The resulting matrix
Q can be applied to data to remove the influence of confound regressors.
}
\examples{
\dontrun{
# Create confound matrix (intercept + linear trend)
n <- 100
X_confounds <- cbind(1, 1:n)

# Get projection matrix
Q <- project_confounds(X_confounds)

# Apply to data to remove confounds
Y_clean <- Q \%*\% Y_raw
}
} 
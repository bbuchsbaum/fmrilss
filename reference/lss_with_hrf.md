# Perform LSS using Voxel-wise HRFs

Computes trial-wise beta estimates using voxel-specific HRFs.

## Usage

``` r
lss_with_hrf(
  Y,
  events,
  hrf_estimates,
  nuisance_regs = NULL,
  engine = "R",
  chunk_size = 5000,
  verbose = TRUE,
  backing_dir = NULL
)
```

## Arguments

- Y:

  Numeric matrix of BOLD data (time x voxels).

- events:

  Data frame with `onset`, `duration` and `condition` columns.

- hrf_estimates:

  A
  [VoxelHRF](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
  object returned by `estimate_voxel_hrf`.

- nuisance_regs:

  Optional numeric matrix of nuisance regressors.

- engine:

  Computational engine: "R" for pure R implementation (default), "C++"
  for optimized C++ (experimental).

- chunk_size:

  Number of voxels to process per batch (C++ engine only).

- verbose:

  Logical; display progress bar.

- backing_dir:

  Directory for bigmemory backing files. If NULL, a temporary directory
  is used (C++ engine only).

## Value

An object of class
[LSSBeta](https://bbuchsbaum.github.io/fmrilss/reference/LSSBeta.md) for
C++ engine, or a numeric matrix (n_trials x n_vox) for R engine.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1,
                     condition = "A")
basis <- fmrihrf::hrf_gamma()
sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
times <- fmrihrf::samples(sframe, global = TRUE)
rset <- fmrihrf::regressor_set(onsets = events$onset,
                               fac = factor(1:nrow(events)),
                               hrf = basis, duration = events$duration,
                               span = 30)
X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
Y <- X %*% coef + Y * 0.1
est <- estimate_voxel_hrf(Y, events, basis)
betas <- lss_with_hrf(Y, events, est, verbose = FALSE, engine = "R")
dim(betas)
} # }
```

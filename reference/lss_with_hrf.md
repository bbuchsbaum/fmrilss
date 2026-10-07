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
  backing_dir = NULL,
  sframe = NULL,
  fixed_regs = NULL
)
```

## Arguments

- Y:

  Numeric matrix of BOLD data (time x voxels).

- events:

  Data frame with `onset`, `duration` and `condition` columns. For a
  multi-run sampling frame, a `run` column is required and onsets are
  physical times relative to that run.

- hrf_estimates:

  A
  [VoxelHRF](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
  object returned by `estimate_voxel_hrf`.

- nuisance_regs:

  Optional finite numeric matrix of nuisance regressors.

- engine:

  Computational engine: "R" for the pure-R implementation (default), or
  "C++" to request the compiled backend. The returned metadata records
  the backend actually used after fallback.

- chunk_size:

  Number of voxels to process per batch (C++ engine only).

- verbose:

  Logical; emit progress messages.

- backing_dir:

  Directory for bigmemory backing files. If NULL, a temporary directory
  is used (C++ engine only).

- sframe:

  Explicit sampling frame. Defaults to the frame stored in
  `hrf_estimates`; no unit-TR fallback is used.

- fixed_regs:

  Optional finite numeric matrix of fixed/common regressors. An
  intercept is added when absent.

## Value

An object of class
[LSSBeta](https://bbuchsbaum.github.io/fmrilss/reference/LSSBeta.md) for
the C++ request, or a numeric matrix (n_trials x n_vox) for the R
engine. Matrix outputs carry sampling-frame, normalization,
coefficient-unit, event-amplitude, event-duration, requested-engine,
realized-engine, and chunk-size metadata. With a
positive-peak-normalized shape, a zero-duration unit-amplitude event has
a peak-response-amplitude coefficient. Otherwise the result is a
coefficient on the supplied duration- and amplitude-coded event design.
The `degenerate` metadata repeats the per-voxel flags of
`hrf_estimates`: flagged voxels have no positive-peak shape, so their
coefficients are not in peak-response units.

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1,
                     condition = "A")
basis <- fmrihrf::HRF_SPMG1
sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
times <- fmrihrf::samples(sframe, global = TRUE)
rset <- fmrihrf::regressor_set(onsets = events$onset,
                               fac = factor(rep("all events", nrow(events))),
                               hrf = basis, duration = events$duration,
                               span = 30)
X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
Y <- X %*% coef + Y * 0.1
est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
betas <- lss_with_hrf(Y, events, est, verbose = FALSE, engine = "R")
dim(betas)
#> [1] 2 2
# }
```

# Fused Single-Pass LSS Solver (C++)

Computes Least Squares-Separate (LSS) beta estimates by residualizing
the trial design against the confounds, forming the n x T LSS weight
matrix once, and applying it to the data with a single matrix product.
The data matrix is never residualized: each weight vector already lies
in the confound residual space.

## Usage

``` r
lss_fused_optim_cpp(
  X,
  Y,
  C,
  block_size = 96L,
  groups = NULL,
  use_omp = TRUE,
  ridge_x = 0,
  ridge_b = 0
)
```

## Arguments

- X:

  The confound regressor matrix (n x k).

- Y:

  The data matrix (n x V).

- C:

  The trial-wise design matrix (n x T).

- block_size:

  The number of voxels per OpenMP block when `use_omp = TRUE`.

- groups:

  Optional 1-based integer trial group codes (LSS-N); NULL for a single
  pooled "other trials" regressor.

- use_omp:

  Logical; distribute voxel blocks across OpenMP threads. Useful with a
  single-threaded BLAS. With a multithreaded BLAS a single matrix
  product is faster.

- ridge_x, ridge_b:

  Fractional ridge penalties on the trial and other-trial coefficients.

## Value

A T x V matrix of LSS beta estimates.

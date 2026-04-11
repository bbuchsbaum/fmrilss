# Fused Single-Pass LSS Solver (C++)

This function computes Least Squares-Separate (LSS) beta estimates using
a memory-efficient, single-pass algorithm. It fuses the projection and
estimation steps, processing voxels in parallel blocks to maximize cache
efficiency.

## Usage

``` r
lss_fused_optim_cpp(X, Y, C, block_size = 96L)
```

## Arguments

- X:

  The nuisance regressor matrix (confounds).

- Y:

  The data matrix (e.g., fMRI data).

- C:

  The trial-wise design matrix.

- block_size:

  The number of voxels to process in each parallel block.

## Value

A matrix of LSS beta estimates.

# OASIS backend for fmrilss::lss (internal entry)

Add `method="oasis"` to fmrilss::lss(). This path:

- (optionally) builds trial-wise design X via fmrihrf

- residualizes Y (and X downstream) against confounds + Z +
  other-condition aggregates

- computes all trial betas in one batched pass via the closed-form LSS
  (exact; or ridge-LSS)

- optionally returns per-trial SEs and design diagnostics

## Usage

``` r
.lss_oasis(
  Y,
  X = NULL,
  Z = NULL,
  Nuisance = NULL,
  oasis = list(),
  prewhiten = NULL
)
```

## Arguments

- Y:

  (T x V) numeric matrix

- X:

  (T x N_trials) trial-wise design (if NULL, use oasis\$design_spec to
  build)

- Z:

  (T x K) fixed experimental regressors to be projected out

- Nuisance:

  (T x P) confounds (intercept, motion, drift, aCompCor, ...)

- oasis:

  list of options:

  - design_spec: list describing events/HRF to build X via fmrihrf

  - K: explicit basis dimension (recommended when X is supplied
    directly)

  - infer_K_from_X: logical; if TRUE and K is missing, infer basis
    dimension heuristically from X (default FALSE for safety)

  - ridge_mode: "fractional" (default) or "absolute"

  - ridge_x, ridge_b: nonnegative ridge on the `a_j` / `b_j` Gram
    (default 0.05 in fractional mode)

  - block_cols: voxel block size (default 4096)

  - return_se: logical (default FALSE)

  - return_diag: logical (default FALSE)

- prewhiten:

  list of prewhitening options using fmriAR (see
  [`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) and
  [`prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
  for details). The legacy `oasis$whiten` field is ignored; use this
  parameter instead.

## Value

by default: (N_trials x V) matrix of betas; if `return_se` or
`return_diag`, a list

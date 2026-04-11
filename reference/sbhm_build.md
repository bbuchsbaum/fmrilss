# Build a Shared-Basis HRF Library (SBHM)

Learn a low-rank shared time basis from a parameterized HRF library
using \`fmrihrf::hrf_library()\`. The library is evaluated on the TR
grid, optionally baseline-removed and L2-normalized, and decomposed via
SVD into \`B = U_r\` (shared basis), singular values \`S\`, and library
coordinates \`A = diag(S)

## Usage

``` r
sbhm_build(
  library_spec = NULL,
  library_H = NULL,
  r = 6,
  sframe = NULL,
  tgrid = NULL,
  span = 32,
  normalize = TRUE,
  baseline = c(0, 0.5),
  shifts = NULL,
  ref = c("mean", "spmg1")
)
```

## Arguments

- library_spec:

  Either NULL (when \`library_H\` is provided) or a list with: -
  \`fun\`: a function compatible with \`fmrihrf::hrf_library(fun, pgrid,
  ...)\` that returns an \`fmrihrf\` HRF object when called with
  parameters. - \`pgrid\`: a data.frame of parameter combinations (see
  examples). - \`span\`: numeric, HRF span in seconds (default
  \`span\`). - \`precision\`: numeric, evaluation precision (default 0.1
  sec). - \`method\`: evaluation method for \`fmrihrf::evaluate()\`
  (default "conv"). - \`extras\`: optional list of additional arguments
  passed to \`hrf_library\`.

- library_H:

  Optional precomputed TxK matrix of candidate HRFs, already aligned to
  the TR grid \`tgrid\` (or \`sframe\`). Mutually exclusive with
  \`library_spec\`.

- r:

  Target rank for the shared basis (default 6). Clipped to \`min(T,
  K)\`.

- sframe:

  Optional \`fmrihrf::sampling_frame\`, used to derive the global time
  grid when \`tgrid\` is not provided.

- tgrid:

  Optional numeric vector of global times (in seconds). If provided,
  takes precedence over \`sframe\`.

- span:

  HRF span in seconds (default 32). Used for reference HRF when needed.

- normalize:

  Logical, L2-normalize library columns (default TRUE).

- baseline:

  Numeric length-2 vector specifying a time window (in seconds) used for
  baseline removal (column-wise mean subtraction within this window).
  Set to NULL to skip.

- shifts:

  Optional numeric vector of time shifts (in seconds). When provided,
  shifted copies of the library are added by linear interpolation on
  \`tgrid\`.

- ref:

  Reference for coefficient-space shrinkage and orientation. One of
  \`"mean"\` (default) or \`"spmg1"\`. If \`"spmg1"\`, the SPMG1 HRF is
  projected onto the learned basis to form \`alpha_ref\`.

## Value

A list with components: - \`B\` (Txr): shared orthonormal time basis -
\`S\` (length r): singular values - \`A\` (rxK): coordinates of library
HRFs in the shared basis - \`tgrid\`: the global time grid used
(seconds) - \`span\`: span used for reference HRF - \`ref\`: list with
\`alpha_ref\` (length r) and \`name\` - \`meta\`: list with \`r\`,
\`K\`, \`normalize\`, \`baseline\`

## Examples

``` r
if (FALSE) { # \dontrun{
  library(fmrihrf)
  param_grid <- expand.grid(shape = c(6, 8, 10), rate = c(0.9, 1.0, 1.1))
  gamma_fun  <- function(shape, rate) fmrihrf::as_hrf(
    fmrihrf::hrf_gamma, params = list(shape = shape, rate = rate)
  )

  sframe <- fmrihrf::sampling_frame(blocklens = 200, TR = 1)
  sbhm <- sbhm_build(
    library_spec = list(fun = gamma_fun, pgrid = param_grid, span = 32),
    r = 6, sframe = sframe, baseline = c(0, 0.5)
  )

  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
} # }
```

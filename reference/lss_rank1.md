# Rank-1 GLM: joint voxel-wise HRF and trial amplitude estimation

Fits the rank-1 GLM of Pedregosa et al. (2015): within each voxel one
HRF, expressed in an HRF basis, is shared by all trials, and each trial
has its own amplitude. With `model = "separate"` (R1-GLMS, the default
and the best performer in that paper) every trial is fitted in its own
least-squares separate (LSS) model, so the amplitudes are LSS estimates
under the voxel's estimated HRF. With `model = "joint"` (R1-GLM) all
trials enter one least-squares-all model.

## Usage

``` r
lss_rank1(
  Y,
  events,
  basis,
  sframe,
  nuisance_regs = NULL,
  fixed_regs = NULL,
  model = c("separate", "joint"),
  trial_groups = NULL,
  init = "aggregate",
  ref_hrf = NULL,
  prewhiten = NULL,
  solver = c("als", "lbfgs"),
  max_iter = 100L,
  tol = 1e-07
)
```

## Arguments

- Y:

  Numeric matrix of BOLD data (time x voxels).

- events:

  Data frame with `onset`, `duration` and `condition` columns (and `run`
  for multi-run sampling frames, with run-relative onsets), as for
  [`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md).
  Each row is one trial. Conditions are labels only: all trials share
  each voxel's HRF.

- basis:

  An `fmrihrf` HRF basis, e.g.
  [`fmrihrf::HRF_SPMG3`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF_objects.html)
  or an FIR basis from
  [`fmrihrf::hrf_fir_generator()`](https://bbuchsbaum.github.io/fmrihrf/reference/hrf_fir_generator.html)/[`fmrihrf::HRF_FIR`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF_objects.html).

- sframe:

  An `fmrihrf` sampling frame with `nrow(Y)` scans.

- nuisance_regs:

  Optional numeric matrix of nuisance regressors.

- fixed_regs:

  Optional numeric matrix of fixed regressors (run intercepts are added
  when not already spanned).

- model:

  `"separate"` (R1-GLMS) or `"joint"` (R1-GLM).

- trial_groups:

  Optional vector with one condition label per trial (row of `events`),
  e.g. `events$condition`, for `model = "separate"`. Each trial-wise
  model then has one "other trials" regressor per group (LSS-N, as in
  [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)).
  Because the HRF is learned mostly through these group regressors,
  supply `trial_groups` whenever conditions may differ in mean response,
  in particular in sign. See Details.

- init:

  `"aggregate"`, `"reference"`, or a K x V numeric matrix (or length-K
  vector) of starting HRF coefficients.

- ref_hrf:

  HRF used for the `"reference"` initialization and to orient the sign
  of each estimated HRF. Defaults to
  [`fmrihrf::HRF_SPMG1`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF_objects.html).

- prewhiten:

  Optional prewhitening options (see
  [`prewhiten_options()`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)).
  Global and run pooling are supported; the noise model is estimated
  with one aggregate regressor per basis function.

- solver:

  `"als"` (default): exact alternating least squares, run in parallel
  over voxels with OpenMP. `"lbfgs"`: joint quasi-Newton optimization of
  HRF and amplitudes with L-BFGS-B, as in Pedregosa et al. (2015), on
  the same objective and starting point (serial; separate model only).
  See Details.

- max_iter:

  Maximum number of alternating iterations (for `"lbfgs"`, at least 1000
  quasi-Newton iterations are allowed).

- tol:

  Relative objective change used as the convergence criterion (for
  `"lbfgs"`, passed as `factr = tol / .Machine$double.eps`).

## Value

A list with

- beta:

  Trial amplitudes (trials x voxels), in peak-response units.

- other:

  For `model = "separate"`, the amplitude of the pooled other-trials
  regressor in each trial's model (trials x voxels); `NULL` for the
  joint model.

- hrf:

  A
  [VoxelHRF](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
  object with the unit-peak HRF coefficients (K x voxels), the removed
  `amplitude_scale`, `basis` and `sframe`.

- objective:

  Final residual sum of squares per voxel (summed over the trial-wise
  models for `"separate"`).

- iterations, converged:

  Per-voxel iteration counts and convergence flags.

- degenerate:

  Voxels whose estimated HRF has no positive peak of at least 5% of its
  largest absolute deflection after orientation (scaled by that
  deflection instead, so their amplitudes are not in peak-response
  units). Also stored in `hrf`. The policy is shared with
  [`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md):
  such voxels are flagged and reported in one warning rather than
  failing the fit.

- model:

  The fitted model.

When prewhitening is applied the fitted plan is attached as the
`whiten_plan` attribute.

## Details

For trial \\i\\ with basis-convolved design block \\X_i\\ (n x K) and
the sum of all blocks \\T_X\\, the separate model minimizes, per voxel,
\$\$\sum_i \\ y - \beta_i X_i h - r_i (T_X - X_i) h \\^2\$\$ over the
HRF coefficients \\h\\, trial amplitudes \\\beta\\ and other-trial
amplitudes \\r\\; the joint model minimizes \\\\ y - \sum_i \beta_i X_i
h \\^2\\. Fixed regressors (run intercepts are added when absent) and
nuisance regressors are projected from the design and data first, which
gives each trial-wise model its own confound coefficients, as in
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md). (The
original formulation shares one set of confound coefficients across the
trial-wise models.)

The fit uses exact alternating least squares rather than the
quasi-Newton solver of the paper: given \\h\\, the amplitudes are
ordinary LSS (or LSA) estimates; given the amplitudes, \\h\\ solves a
K-dimensional least squares problem. Each step minimizes its block
exactly, so the objective never increases. All quantities either step
needs are K x K Gram blocks of the residualized design and the product
\\X^\top Y\\, so after that one product the separate model costs \\O(T
K^2)\\ per voxel and iteration, independent of the number of scans. The
joint model costs \\O(T^2 K^2 + T^3)\\ per voxel and iteration and is
intended for regions of interest.

The problem is bilinear: \\h\\ and the amplitudes are determined only up
to a common scale and sign. After fitting, each voxel's HRF is oriented
to correlate positively with `ref_hrf` and scaled to a unit positive
peak, and the amplitudes absorb that scale, so they are in peak-response
units for unit-amplitude, zero-duration events, like
[`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md).
The returned `hrf` is a
[VoxelHRF](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
object and can be passed to
[`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md)
to estimate trial amplitudes in new data (for example a held-out run)
with the learned HRFs.

The objective is not jointly convex. `init = "aggregate"` starts each
voxel from the pooled-shape fit of
[`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md)
(one common amplitude for all trials); `init = "reference"` starts every
voxel from the projection of `ref_hrf` onto the basis, which is more
stable for voxels whose mean response is near zero. A K x V matrix of
starting coefficients is also accepted.

## References

Pedregosa, F., Eickenberg, M., Ciuciu, P., Thirion, B., & Gramfort, A.
(2015). Data-driven HRF estimation for encoding and decoding models.
NeuroImage, 104, 209-220.

## See also

[`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md),
[`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md),
[`lss_sbhm()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md)

## Examples

``` r
# \donttest{
set.seed(1)
sframe <- fmrihrf::sampling_frame(blocklens = 200, TR = 1)
events <- data.frame(onset = seq(10, 180, by = 9), duration = 0,
                     condition = "A")
basis <- fmrihrf::HRF_SPMG3
X <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)$X
h <- c(1, 0.4, -0.2)
beta <- matrix(rnorm(nrow(events) * 5, 1, 0.3), nrow(events), 5)
Y <- X %*% kronecker(beta, h) + matrix(rnorm(200 * 5, sd = 0.5), 200, 5)
fit <- lss_rank1(Y, events, basis, sframe)
dim(fit$beta)
#> [1] 19  5
fit$hrf$coefficients[, 1]
#>   basis_1   basis_2   basis_3 
#>  5.217199  4.814784 -1.088945 
# }
```

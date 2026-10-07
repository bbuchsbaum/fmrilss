# Rank-1 GLM: Learning One HRF per Voxel While Estimating Trials

The haemodynamic response differs between voxels, and a canonical HRF
that is wrong for a voxel biases every trial estimate in it. Estimating
a free HRF per trial is hopelessly noisy. The *rank-1 GLM* of Pedregosa
et al. (2015) sits in between: within a voxel, every trial shares one
HRF, expressed in a small basis, and each trial has its own amplitude.
[`lss_rank1()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_rank1.md)
implements that model, including its least-squares-separate variant,
with a solver designed for whole-brain data.

This article explains the model, how fmrilss fits it and why that solver
was chosen over the quasi-Newton method of the original paper, and which
variant to use for which design. Read
[`vignette("voxel-wise-hrf")`](https://bbuchsbaum.github.io/fmrilss/articles/voxel-wise-hrf.md)
first for HRF bases and the `VoxelHRF` object.

## The model

Let $`y`$ be one voxel’s time series of $`n`$ scans, and let the design
have $`T`$ trials. For trial $`i`$, let $`X_i`$ ($`n \times K`$) be the
trial’s onsets convolved with each of $`K`$ basis functions, $`\beta_i`$
the trial’s amplitude, $`h \in \mathbb{R}^K`$ the voxel’s HRF
coefficients and $`T_X = \sum_j X_j`$ the summed blocks of all trials.
[`lss_rank1()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_rank1.md)
fits, per voxel, either

- `model = "joint"` (R1-GLM): one least-squares-all (LSA) model,
  ``` math
   \min_{h, \beta} \; \Big\| y - \sum_i \beta_i X_i h \Big\|^2, 
  ```
- `model = "separate"` (R1-GLMS, the default): one
  least-squares-separate (LSS) model per trial, all sharing $`h`$, where
  $`r_i`$ is the amplitude of the “other trials” regressor in trial
  $`i`$’s model,
  ``` math
   \min_{h, \beta, r} \; \sum_i \big\| y - \beta_i X_i h - r_i (T_X - X_i) h \big\|^2 . 
  ```

With `trial_groups`, the single “other trials” term of the separate
model becomes one term per group,
$`\sum_g r_{ig} (A_g - [g = g_i] X_i) h`$, where $`A_g`$ is the summed
blocks of group $`g`$, $`g_i`$ is trial $`i`$’s group, and $`[g = g_i]`$
is 1 for that group and 0 otherwise. With $`G`$ groups, this is the
LSS-N design of
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md).

Fixed regressors (`fixed_regs`, plus run intercepts when they are
absent) and nuisance regressors (`nuisance_regs`) are projected out of
the data and design first, which gives each trial-wise model its own
confound coefficients as in
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md). The
fitted HRF is oriented to correlate positively with a reference
(`ref_hrf`, the canonical HRF by default) and scaled to a unit positive
peak; the amplitudes absorb that scale.

## How `lss_rank1()` fits the model, and why not L-BFGS-B

The objectives are bilinear: linear in $`h`$ for fixed amplitudes, and
linear in the amplitudes for fixed $`h`$. Pedregosa et al. minimize them
jointly with L-BFGS-B, using a box constraint and a penalty to control
the scale ambiguity $`h \to c\,h,\ \beta \to \beta / c`$: multiplying
the HRF by any nonzero $`c`$ and dividing the amplitudes by $`c`$ leaves
the fit unchanged.
[`lss_rank1()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_rank1.md)
instead uses **exact alternating least squares (ALS)**:

1.  *Amplitude step.* For fixed $`h`$, the trial regressors $`X_i h`$
    are fixed, so the amplitudes are ordinary LSS (or LSS-N, or LSA)
    estimates: one $`2 \times 2`$ (or $`(G+1) \times (G+1)`$) solve per
    trial.
2.  *HRF step.* For fixed amplitudes the objective is a least-squares
    problem in $`K`$ unknowns: one $`K \times K`$ solve.
3.  Rescale $`h`$ to unit norm (the objective is unchanged) and repeat
    until the relative change in the objective is below `tol`.

Each step minimizes the objective exactly over its block, so the
objective never increases, there is no step size or line search, and the
scale ambiguity is removed by step 3 instead of by constraints.

Both steps need only small Gram blocks. Writing
$`G_{ii} = X_i^\top X_i`$, $`S_i = X_i^\top T_X`$,
$`G_{TT} = T_X^\top T_X`$ and $`U_i = X_i^\top y`$, the amplitude step
for trial $`i`$ uses
$`h^\top G_{ii} h,\ h^\top S_i h,\ h^\top G_{TT} h,\ h^\top U_i`$, and
the HRF step solves
``` math
\Big[\sum_i a_i^2 G_{ii} + a_i r_i (S_i + S_i^\top) + r_i^2 G_{TT}\Big] h =
\sum_i a_i U_i + r_i\, T_X^\top y,
\qquad a_i = \beta_i - r_i .
```
So after one matrix product $`X^\top Y`$ for all voxels, an iteration of
the separate model costs $`O(T K^2)`$ per voxel, independent of the
number of scans $`n`$, and voxels are processed in parallel with OpenMP.
The joint model’s amplitude step is a $`T \times T`$ solve,
$`O(T^2 K^2 + T^3)`$ per voxel, so it suits regions of interest.

### Benchmark: ALS against L-BFGS-B

To decide between the two solvers on evidence, fmrilss also ships the
paper’s approach as `solver = "lbfgs"`: L-BFGS-B over $`(h, \beta, r)`$
with the analytic gradient of the *same* objective, computed from the
same Gram blocks, and started from the same point.
`bench/rank1/run_rank1_bench.R` in the package sources compares them
single-threaded on 1,000 voxels whose HRFs peak between 4 and 7.5 s, for
a moderate design (inter-trial interval 3–8 s) and a rapid, noisy one
(1.5–4 s), two 30-second bases (SPMG3 and a 15-bin FIR), and both
initializations (`init = "aggregate"` and `init = "reference"`,
described under the cautions below):

| Basis | Tolerance | ALS seconds | L-BFGS-B seconds | ALS iterations | L-BFGS-B evaluations | Worst objective gap, ALS | Worst objective gap, L-BFGS-B |
|:---|:---|:---|:---|:---|:---|:---|:---|
| SPMG3 | default | 0.14-0.19 | 0.48-8.6 | 2-4 | 9-164 | 9e-9 | 1e-4 |
| SPMG3 | tight | 0.12-0.19 | 8.2-37 | 3-6 | 275-750 | 9e-12 | 2e-7 |
| FIR (15) | default | 0.24-0.58 | 2.0-30 | 3-4 | 9-195 | 1e-8 | 4e-3 |
| FIR (15) | tight | 0.29-0.45 | 24-87 | 4-6 | 304-588 | 8e-12 | 2e-3 |

Ranges over the two designs and two initializations (1,000 voxels, one
thread). Seconds are for all 1,000 voxels; iteration and evaluation
counts are per-voxel medians. Default tolerance is tol = 1e-7 for ALS
(the lss_rank1() default) and factr = 1e7 for L-BFGS-B (the optim()
default); tight is tol = 1e-10 and factr = 10. The objective gap is
relative to the best objective any solver reached for that voxel.
{.table}

Run to tight tolerances, ALS reaches the lowest objective either solver
found in every voxel (worst relative gap about $`10^{-11}`$). L-BFGS-B
reaches the same solution in almost every voxel, but in a few it stops
short (worst gap about $`2 \times 10^{-7}`$ with SPMG3 and
$`2 \times 10^{-3}`$ with FIR). ALS gets there in a handful of
iterations: 3–8 times faster than L-BFGS-B at its default tolerance when
started from `init = "aggregate"`, and roughly 25–200 times faster when
started from `init = "reference"` or run to a tight tolerance. At the
default tolerance, ALS is also closer to the optimum. Two properties
explain this. The amplitude step solves for all $`2T`$ amplitudes
exactly, which removes most of the coupling a gradient method must
discover iteration by iteration. And the scale ambiguity leaves the
joint objective with a flat, badly conditioned direction that slows
quasi-Newton methods but does not affect ALS, which renormalizes $`h`$.
ALS also parallelizes over voxels, whereas R’s L-BFGS-B is not
thread-safe. `solver = "lbfgs"` remains available, single-threaded, for
comparison.

The live check below fits 60 simulated voxels with both solvers at tight
tolerances and reports the largest objective and amplitude differences
between them:

``` r

simulate_voxels <- function(n_voxels, iti = c(3, 8), amplitude = c("positive", "opposite"),
                            noise_sd = 0.8, n = 300, seed = 1) {
  amplitude <- match.arg(amplitude)
  set.seed(seed)
  sframe <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  onsets <- cumsum(c(12, runif(400, iti[1], iti[2])))
  onsets <- onsets[onsets < n - 30]
  events <- data.frame(
    onset = onsets, duration = 0,
    condition = rep(c("A", "B"), length.out = length(onsets))
  )
  # each voxel gets a double-gamma HRF peaking between 4 and 7.5 s
  peak <- runif(n_voxels, 4, 7.5)
  grid <- seq(0, 30, by = 0.1)
  hrf_of <- function(p) {
    a <- p / 0.9 + 1
    h <- dgamma(grid, a, scale = 0.9) - dgamma(grid, a + 10, scale = 0.9) / 6
    h / max(h)
  }
  stick <- function(o) {
    u <- numeric(n * 10)
    u[round(o * 10) + 1] <- 1
    u
  }
  sticks <- vapply(onsets, stick, numeric(n * 10))
  amp_mean <- if (amplitude == "positive") 1 else ifelse(events$condition == "A", 1, -1)
  betas <- matrix(rnorm(length(onsets) * n_voxels, amp_mean, 0.5),
                  length(onsets), n_voxels)
  Y <- vapply(seq_len(n_voxels), function(v) {
    drive <- drop(sticks %*% betas[, v])
    bold <- stats::filter(c(numeric(length(grid) - 1), drive), hrf_of(peak[v]),
                          sides = 1)[-seq_len(length(grid) - 1)]
    bold[seq(1, n * 10, by = 10)]
  }, numeric(n))
  noise <- apply(matrix(rnorm(n * n_voxels), n, n_voxels), 2, function(e) {
    as.numeric(stats::filter(e, 0.3, "recursive"))
  })
  list(Y = Y + noise_sd * noise, events = events, sframe = sframe,
       betas = betas, peak = peak)
}

sim <- simulate_voxels(60)
# The simulated HRFs last 30 s; give the basis the same support (fmrihrf's
# built-in HRF_SPMG3 declares 24 s and would truncate the undershoot).
basis <- fmrihrf::gen_hrf(fmrihrf::HRF_SPMG3, span = 30)
time_als <- system.time(
  als <- lss_rank1(sim$Y, sim$events, basis, sim$sframe, tol = 1e-10)
)[["elapsed"]]
time_lbfgs <- system.time(
  lbfgs <- lss_rank1(sim$Y, sim$events, basis, sim$sframe, solver = "lbfgs",
                     tol = 10 * .Machine$double.eps)
)[["elapsed"]]
c(
  max_relative_objective_difference =
    max(abs(als$objective - lbfgs$objective) / als$objective),
  max_beta_difference = max(abs(als$beta - lbfgs$beta)),
  als_seconds = time_als,
  lbfgs_seconds = time_lbfgs
)
#> max_relative_objective_difference               max_beta_difference 
#>                      1.832845e-09                      9.227377e-04 
#>                       als_seconds                     lbfgs_seconds 
#>                      7.400000e-02                      5.340000e-01
```

## Which variant to use

The same benchmark compares the estimators by the mean per-voxel
correlation between estimated and true trial amplitudes, across four
amplitude regimes. The oracle rows know each voxel’s true HRF. Numbers
in parentheses below are correlations from this table.

| Estimator | Positive amplitudes | Positive, rapid and noisy | Zero-mean amplitudes | Opposite-signed conditions |
|:---|---:|---:|---:|---:|
| Oracle LSS (true HRF) | 0.633 | 0.390 | 0.807 | 0.649 |
| Oracle LSS-N (true HRF) | 0.632 | 0.375 | 0.807 | 0.880 |
| LSS, canonical HRF | 0.559 | 0.341 | 0.729 | 0.505 |
| LSS-N, canonical HRF | 0.557 | 0.329 | 0.728 | 0.691 |
| estimate_voxel_hrf() + lss_with_hrf() | 0.620 | 0.368 | 0.453 | 0.147 |
| lss_rank1(), separate | 0.620 | 0.368 | 0.599 | 0.270 |
| lss_rank1(), separate + trial_groups | 0.616 | 0.339 | 0.592 | 0.876 |
| lss_rank1(), joint | 0.571 | 0.117 | 0.767 | 0.640 |

Trial-amplitude recovery (2,000 voxels, 30-second SPMG3 basis). {.table}

- **Trials that share a response sign** (the common activation setting):
  the separate model recovers much of the gap between the canonical HRF
  and the oracle (about 80% of it in the moderate design and 55% in the
  rapid, noisy one), in about 0.1–0.2 s per thousand voxels on one
  thread (first table).
- **Conditions with different mean responses, especially opposite
  signs:** supply `trial_groups`. In the separate model the HRF is
  learned mainly through the “other trials” regressors. With one pooled
  regressor, opposite responses cancel, the HRF is learned poorly, and
  the amplitude correlation is only 0.27. With one regressor per
  condition, each condition informs the shared HRF with its own
  amplitude (0.88, matching the LSS-N oracle).
- **Responses that vary around zero within a condition:** the pooled
  regressors carry little information about the HRF, and the joint
  model, which learns from every trial’s amplitude, is better (0.77).
  Prefer it when the design is not rapid and the region is small enough
  for its cost.
- **Rapid designs:** the joint model inherits LSA’s variance and
  degrades sharply (0.12); use the separate model.
- **Basis size:** in these simulations the 3-function SPMG3 basis
  outperformed a 15-bin FIR basis for amplitude recovery. More flexible
  bases pay in variance, as Pedregosa et al. also found for decoding.

Two further cautions. First, the separate model inherits LSS’s
approximation: with overlapping trials its pooled “other trials”
regressor cannot represent neighbouring trials’ different amplitudes, so
it is biased even without noise. Second, the objectives are not convex.
`init = "aggregate"` (the pooled shape of
[`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md))
and `init = "reference"` usually reach the same solution, but voxels
with little signal can differ.

## Learning HRFs on one run and applying them to another

The `hrf` element is a `VoxelHRF` object, so the learned HRFs can
estimate trials in new data with
[`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md).
Fitting the HRF on training data and applying it to held-out data avoids
circularity in encoding and decoding analyses.

``` r

train <- simulate_voxels(40, seed = 2)
fit <- lss_rank1(train$Y, train$events, basis, train$sframe)
str(fit$hrf$coefficients)
#>  num [1:3, 1:40] 5.88 -16.99 -7.1 6.36 -5.43 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:3] "basis_1" "basis_2" "basis_3"
#>   ..$ : NULL

# Re-estimating trials in the training data with the learned HRFs reproduces
# the separate model's amplitudes exactly.
refit <- lss_with_hrf(train$Y, train$events, fit$hrf, verbose = FALSE)
max(abs(unclass(refit)[seq_len(nrow(train$events)), ] - fit$beta))
#> [1] 1.64313e-14
```

For held-out data, call `lss_with_hrf(Y_test, events_test, fit$hrf)`
with the test run’s data and events. By default it uses the sampling
frame stored in the `VoxelHRF` object, which is the training run’s; if
the test run differs in length or TR, pass its own frame through
`sframe`
([`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md)
stops when the frame does not match `nrow(Y_test)`).

## Reading the result

``` r

names(fit)
#> [1] "beta"       "other"      "hrf"        "objective"  "iterations"
#> [6] "converged"  "degenerate" "model"
summary(fit$iterations)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    2.00    2.00    3.00    2.55    3.00    3.00
all(fit$converged)
#> [1] TRUE
```

- `beta`: trial amplitudes in peak-response units (trials x voxels).
- `other`: the “other trials” amplitudes of each trial’s model (trials x
  voxels, or trials x groups x voxels with `trial_groups`; separate
  model only).
- `hrf`: the unit-peak HRF coefficients as a `VoxelHRF` object.
- `objective`, `iterations`, `converged`: per-voxel fit diagnostics.
- `degenerate`: voxels whose estimated HRF has no positive peak of at
  least 5% of its largest absolute deflection after orientation,
  typically voxels without signal. They are scaled by their largest
  absolute value instead, so their amplitudes are not in peak-response
  units, and a warning reports how many there are. This is the same
  policy and normalization code as
  [`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md);
  the flags are also stored in `fit$hrf$degenerate` and passed on by
  [`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md).

`prewhiten` accepts global or run pooling; the noise model is estimated
from one aggregate regressor per basis function (see
[`vignette("prewhitening")`](https://bbuchsbaum.github.io/fmrilss/articles/prewhitening.md)).

## Reference

Pedregosa, F., Eickenberg, M., Ciuciu, P., Thirion, B., & Gramfort, A.
(2015). Data-driven HRF estimation for encoding and decoding models.
*NeuroImage*, 104, 209–220.

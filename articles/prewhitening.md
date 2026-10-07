# Prewhitening: Where the Noise Model Comes From

fMRI noise is temporally autocorrelated. Prewhitening estimates that
autocorrelation, builds a filter that removes it, and applies the same
filter to the data and to every design matrix, so that ordinary least
squares on the filtered system is generalized least squares on the
original one. When the noise model is right, trial estimates become more
precise.

The filter is only as good as the autocorrelation estimate, and that
estimate has to come from *residuals of some fitted model*: the noise
itself is never observed. This article explains why the choice of that
model matters a great deal for single-trial analyses, what the three
choices offered by `prewhiten$residual_model` do, and how to pick one.
Read it after
[`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md).

## Why the residual model matters

Let $`M = I - D(D^\top D)^{-1}D^\top`$ be the residual-forming matrix of
a fitted design $`D`$ with $`q`$ columns and $`n`$ scans. The residuals
are $`r = M y`$. If the noise were white with variance $`\sigma^2`$, the
expected lag-$`k`$ autocovariance of the residuals would be

``` math
\mathbb{E}\Big[\tfrac{1}{n}\textstyle\sum_t r_t r_{t+k}\Big] =
\tfrac{\sigma^2}{n}\operatorname{tr}(M L_k),
```

where $`L_k`$ is the lag-$`k`$ shift matrix (ones on the $`k`$-th
superdiagonal, zeros elsewhere). White noise has zero autocovariance at
every lag $`k \ge 1`$, but this expectation is not zero. Fitting $`D`$
absorbs part of the noise into the fitted values, and because
HRF-convolved trial regressors are smooth, the part they absorb is the
slowly varying, positively autocorrelated part. What is left looks
rougher than the noise really is, so the autocorrelation estimate is
biased **downward**. The size of the bias grows with $`q / n`$.

That is harmless for a block design with a handful of regressors. It is
not harmless for single-trial models, where the natural “full” design,
the least-squares-all (LSA) model with one column per trial, can have a
column for every two or three scans.

The opposite mistake is to fit too little. If the residual model leaves
real signal in the residuals, for example trial-to-trial amplitude
variability, that smooth signal is mistaken for autocorrelated noise and
the estimate is biased **upward**.

## The three residual models

| `residual_model` | Trial part of the fitted design | Bias of the autocorrelation estimate | Cost | Pooling |
|----|----|----|----|----|
| `"aggregate"` (default) | One summed regressor per `trial_groups` level (or per basis function) | Upward when trial-to-trial variability is large relative to noise | Fast | All (voxel and parcel: weight-matrix methods only) |
| `"full"` | Every trial regressor (LSA) | Downward, severe when trials are a sizeable fraction of scans | Fast | All (voxel and parcel: weight-matrix methods only) |
| `"corrected"` | Every trial regressor, then fmriAR’s residual-bias correction | Approximately unbiased | Seconds; grows with the scan count, independent of voxels | `"global"`, `"run"` |

In every case the confounds (`Z`, `Nuisance`, and an intercept if none
is present) are part of the fitted design.

- **`"aggregate"`** is the model the per-trial GLMs of Nilearn and
  NiBetaSeries effectively estimate their noise from: the
  condition-level response is removed, trial-by-trial deviations from it
  are not.
- **`"full"`** was the only behaviour in fmrilss 0.2.0 and earlier. Use
  it to reproduce earlier results.
- **`"corrected"`** fits the full model and then inverts the linear map
  $`\gamma \mapsto \operatorname{tr}(M L_k M \Gamma)`$ from true to
  expected residual autocovariance
  ([`fmriAR::acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.html)),
  where $`\gamma`$ is the noise autocovariance sequence and $`\Gamma`$
  the noise covariance matrix built from it. Neither the over-fitting of
  the full model nor unmodelled trial variability then distorts the
  estimate. It requires `method = "ar"`. fmrilss builds the correction
  design for you. Supplying `prewhiten$design` yourself gives the same
  correction; `residual_model` then defaults to `"full"`, and any other
  value is an error.

## A simulation

The simulation below draws AR(1) noise with $`\phi = 0.4`$ and unit
variance, and trial amplitudes with mean 1 and either modest (SD 0.5) or
large (SD 2) trial-to-trial variability, in a rapid design (inter-trial
interval 1.5–4 s, about one trial per three scans) and a slow design
(8–12 s).

``` r

canonical_hrf <- function(t) {
  h <- dgamma(t, 6) - dgamma(t, 16) / 6
  h / max(h)
}

simulate_run <- function(iti, trial_sd, n = 300, n_voxels = 200, phi = 0.4,
                         seed = 1) {
  set.seed(seed)
  onsets <- cumsum(c(20, runif(1000, iti[1], iti[2])))
  onsets <- onsets[onsets < n - 25]
  X <- vapply(onsets, function(onset) {
    lag <- seq_len(n) - 1 - onset
    ifelse(lag >= 0, canonical_hrf(pmax(lag, 0)), 0)
  }, numeric(n))
  betas <- matrix(rnorm(ncol(X) * n_voxels, 1, trial_sd), ncol(X), n_voxels)
  innovations <- matrix(rnorm(n * n_voxels, sd = sqrt(1 - phi^2)), n, n_voxels)
  noise <- apply(innovations, 2, function(e) {
    as.numeric(stats::filter(e, phi, method = "recursive"))
  })
  list(Y = X %*% betas + noise, X = X, betas = betas)
}

compare_models <- function(label, run) {
  rmse <- function(estimate) sqrt(mean((estimate - run$betas)^2))
  fit <- function(model) {
    lss(run$Y, run$X,
        prewhiten = list(method = "ar", p = 1, residual_model = model))
  }
  fits <- lapply(c(aggregate = "aggregate", full = "full",
                   corrected = "corrected"), fit)
  phi_hat <- vapply(fits, function(f) attr(f, "whiten_plan")$phi[[1]], 1)
  data.frame(
    design = label,
    trials = ncol(run$X),
    phi_full = phi_hat[["full"]],
    phi_aggregate = phi_hat[["aggregate"]],
    phi_corrected = phi_hat[["corrected"]],
    rmse_ols = rmse(lss(run$Y, run$X)),
    rmse_full = rmse(fits$full),
    rmse_aggregate = rmse(fits$aggregate),
    rmse_corrected = rmse(fits$corrected)
  )
}

scenarios <- rbind(
  compare_models("rapid, trial SD 0.5", simulate_run(c(1.5, 4), 0.5)),
  compare_models("slow, trial SD 0.5", simulate_run(c(8, 12), 0.5)),
  compare_models("rapid, trial SD 2", simulate_run(c(1.5, 4), 2)),
  compare_models("slow, trial SD 2", simulate_run(c(8, 12), 2))
)
```

| Design | Trials | phi: full | phi: aggregate | phi: corrected | RMSE: OLS | RMSE: full | RMSE: aggregate | RMSE: corrected |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| rapid, trial SD 0.5 | 93 | -0.12 | 0.52 | 0.38 | 0.89 | 0.89 | 0.84 | 0.86 |
| slow, trial SD 0.5 | 26 | 0.27 | 0.44 | 0.40 | 0.74 | 0.73 | 0.73 | 0.73 |
| rapid, trial SD 2 | 93 | -0.12 | 0.85 | 0.38 | 2.15 | 2.17 | 1.59 | 1.99 |
| slow, trial SD 2 | 26 | 0.27 | 0.72 | 0.40 | 0.78 | 0.78 | 0.86 | 0.77 |

AR(1) estimates (truth 0.4) and trial-beta RMSE for 300 scans and 200
voxels. {.table}

Reading the table:

- **`"full"` is biased whenever the trial count is not small relative to
  the scan count.** In the rapid design (one trial column per three
  scans) it estimates a *negative* $`\phi`$ for positively
  autocorrelated noise, and the resulting filter gives betas no better
  (slightly worse) than no prewhitening at all (the OLS column). Even in
  the slow design, with about one column per eleven scans, it
  underestimates $`\phi`$ by about a third.
- **`"aggregate"` works well when trial-to-trial variability is modest
  relative to the noise**, which is the common situation for
  single-trial fMRI: its estimate is close in the slow design and
  somewhat high in the rapid one, and it gives the lowest or
  equal-lowest RMSE in both. When the variability is large, the
  deviations it leaves in the residuals inflate $`\phi`$ substantially;
  in the slow design that over-whitening costs accuracy. (In the rapid,
  high-variability design the over-whitened LSS betas happen to have
  lower error, because LSS already mixes neighbouring trials’ deviations
  into each estimate. Treat that as a coincidence of this design, not a
  reason to over-whiten.)
- **`"corrected"` recovers $`\phi`$ in every scenario.** Its RMSE is the
  lowest or close to it in every row except the rapid, high-variability
  one, where the over-whitened `"aggregate"` betas described above have
  lower error.

## Choosing a residual model

1.  Start with the default, `"aggregate"`. It is fast, works with every
    pooling mode, including voxel-adaptive pooling, and matches what
    Nilearn does.

2.  If you can afford it, compare its estimate with `"corrected"` once
    on your data, for example on a subset of voxels or one run:

    ``` r

    est <- function(model) {
      fit <- lss(Y[, voxel_subset], X, Z, Nuisance,
                 prewhiten = list(method = "ar", p = 1, residual_model = model))
      attr(fit, "whiten_plan")$phi
    }
    est("aggregate")
    est("corrected")
    ```

    Close agreement means the default is fine. A clearly larger
    `"aggregate"` estimate means trial-to-trial variability is large
    relative to the noise; prefer `"corrected"` (with global or run
    pooling) for that dataset.

3.  Avoid `"full"` except to reproduce earlier results.

The correction’s cost does not depend on the voxel count but grows
quickly with the number of scans: in our measurements it took about 0.8
s for 400 scans and 150 trials, and about 8 s for 1,200 scans and 300
trials. It is computed once per fit. The correction is designed for
high-pass-filtered designs: without high-pass filtering, the lag budget
it needs (`correction_max_lag`) can become impractically large. See
`correction_max_lag` in
[`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) and
[`fmriAR::fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.html).

## How the aggregate columns are chosen

The aggregate model sums trial regressors so that the condition-level
response is removed from the residuals:

| Caller | Aggregate trial columns |
|----|----|
| [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) with `trial_groups` | One summed regressor per group |
| [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) without `trial_groups` | One regressor summing all trials |
| `lss(method = "oasis")`, multi-basis | One summed regressor per basis function |
| [`lss_sbhm()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md), SBHM amplitude fits | One summed regressor per shared basis function (plus other conditions) |
| [`lss_rank1()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_rank1.md) | One summed regressor per basis function |

Passing `trial_groups` therefore also improves the noise model whenever
conditions differ in their mean response.

## Pooling and the residual model

`residual_model` decides which residuals are used; `pooling` decides how
the autocorrelation is shared across voxels:

- `"global"` estimates one model from all voxels, and `"run"` one model
  per run from all voxels. All three residual models are available.
- `"voxel"`, available only with the weight-matrix methods
  `r_optimized`, `cpp_optimized` and `cpp`, groups voxels into
  `voxel_bins` bins of similar residual autocorrelation, refits one AR
  model per bin, and estimates each bin’s betas from a design filtered
  with that bin’s model. `"parcel"` does the same for user-supplied
  parcels. These use `"aggregate"` or `"full"` residuals; fmriAR’s
  correction is not available for them.

Every prewhitened result stores the fitted noise model as its
`whiten_plan` attribute. Here the default `"aggregate"` model is fitted
to the rapid design with trial SD 0.5, so `plan$phi` matches the first
row of the table:

``` r

run <- simulate_run(c(1.5, 4), 0.5)
fit <- lss(run$Y, run$X, prewhiten = prewhiten_options(method = "ar", p = 1))
plan <- attr(fit, "whiten_plan")
plan$phi
#> [[1]]
#> [1] 0.5160224
```

## Where to go next

- [`?prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
  and [`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)
  list every prewhitening field.
- [`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
  covers prewhitening with OASIS and multi-basis designs.
- [`vignette("lss_with_fmridesign")`](https://bbuchsbaum.github.io/fmrilss/articles/lss_with_fmridesign.md)
  covers multi-run designs and run-aware pooling.

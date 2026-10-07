# GLMsingle in fmrilss: library HRFs, GLMdenoise and fractional ridge

## What GLMsingle estimates

GLMsingle (Prince et al., 2022) estimates one response amplitude per
trial and voxel. It adds three data-driven steps to a single-trial GLM:

1.  **A library HRF per voxel.** Each voxel uses whichever of 20
    canonical HRF shapes best explains its time series (type B).
2.  **GLMdenoise.** Principal components of the residual time series in
    a “noise pool” of task-unresponsive voxels become nuisance
    regressors. The number of components is chosen by cross-validation
    (type C).
3.  **Fractional ridge regression.** Each voxel gets its own amount of
    shrinkage, chosen by cross-validation, and the result is rescaled to
    the unregularised betas (type D).

Both cross-validations compare single-trial betas for the **same
condition in different runs**, so the design must repeat conditions
across runs.

[`glmsingle()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle.md)
computes the same estimator as the reference GLMsingle implementation,
but organises the work so that each run is solved separately, data are
touched once per stage, and every candidate model is scored from small
summary matrices.

## A simulated experiment

Four runs of 110 TRs. Each run shows trials from 8 conditions; every
condition repeats across runs. Voxels differ in HRF shape, and all
voxels share slow structured noise, which GLMdenoise can learn from
voxels that do not respond to the task.

``` r

library(fmrilss)
set.seed(2024)
tr <- 1; stimdur <- 3
n_runs <- 4; n_time <- 110; n_vox <- 120; n_cond <- 8
lib <- glmsingle_hrf_library(stimdur, tr)
true_hrf <- sample(ncol(lib), n_vox, replace = TRUE)
responsive <- seq_len(n_vox) <= 80
cond_mean <- matrix(rnorm(n_cond * n_vox, 1.5, 1), n_cond)
load <- matrix(rnorm(3 * n_vox, 0, 1.5), 3)

design <- Y <- true_beta <- vector("list", n_runs)
for (r in seq_len(n_runs)) {
  onsets <- seq(3, n_time - 15, by = 4) + sample(0:1, 1)
  cond <- sample(rep_len(seq_len(n_cond), length(onsets)))
  D <- matrix(0, n_time, n_cond)
  D[cbind(onsets + 1, cond)] <- 1
  amp <- cond_mean[cond, ] + matrix(rnorm(length(onsets) * n_vox, 0, 0.5), length(onsets))
  amp[, !responsive] <- 0
  signal <- vapply(seq_len(n_vox), function(v) {
    s <- numeric(n_time); s[onsets + 1] <- amp[, v]
    stats::convolve(s, rev(lib[, true_hrf[v]]), type = "open")[seq_len(n_time)]
  }, numeric(n_time))
  latent <- apply(matrix(rnorm(3 * n_time), n_time), 2, cumsum) * 0.4
  noise <- matrix(rnorm(n_time * n_vox), n_time) + latent %*% load
  Y[[r]] <- 1000 * (1 + 0.01 * (signal + noise))
  design[[r]] <- D
  true_beta[[r]] <- amp
}
true_beta <- do.call(rbind, true_beta)  # trials x voxels, chronological
```

`design` is in GLMsingle’s format: one time x condition 0/1 matrix per
run, with a 1 at each trial onset. A data frame with columns `run`,
`onset` (in seconds) and `condition` works too, as does an fmridesign
event model via
[`glmsingle_design()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle_design.md).

## Fitting

``` r

fit <- glmsingle(Y, design, tr = tr, stimdur = stimdur, verbose = FALSE)
fit
#> <glmsingle_fit>
#>   4 runs, 96 trials, 8 conditions, 120 voxels
#>   models: A (ON-OFF), B (HRF library), C (GLMdenoise), D (fractional ridge) 
#>   noise PCs selected: 3
#>   median ridge fraction: 0.70
#>   elapsed: 0.4 s
```

Each model type is a list with GLMsingle’s field names. Betas are trials
x voxels in chronological trial order and, by default, in percent signal
change.

``` r

dim(coef(fit))            # type D
#> [1]  96 120
fit$typed$pcnum           # number of GLMdenoise components
#> [1] 3
table(fit$typed$FRACvalue)
#> 
#> 0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8 
#>   22    5    2    1    2    2    8    4    2    1    1    2    4    6    7    9 
#> 0.85  0.9 0.95    1 
#>   15   13   12    2
```

## Did the data-driven steps help?

The truth is known here, so each model type can be scored by how well
its betas correlate with the true trial amplitudes in responsive voxels.
Betas are in percent signal change and the simulation scales signal by
1% of the baseline, so the scale matches.

``` r

score <- function(B) {
  mean(vapply(which(responsive), function(v) stats::cor(B[, v], true_beta[, v]), numeric(1)))
}
round(c(B = score(coef(fit, "b")), C = score(coef(fit, "c")), D = score(coef(fit, "d"))), 3)
#>     B     C     D 
#> 0.334 0.601 0.660
```

HRF selection can also be checked directly. Neighbouring library HRFs
are very similar, so a near miss is still a good shape:

``` r

mean(abs(fit$typeb$HRFindex[responsive] - true_hrf[responsive]) <= 2)
#> [1] 0.7875
```

## Choices that differ from the reference implementation

The defaults reproduce GLMsingle (Python, commit `1ab54a6`) except where
GLMsingle is internally inconsistent. Each difference is an argument
whose other value reproduces GLMsingle:

| Argument | Default | GLMsingle behaviour |
|----|----|----|
| `extras_in_denoise` | `"always"`: user nuisance regressors (e.g. motion) are in every GLMdenoise and ridge fit | `"with_pcs"`: dropped when zero PCs are used, including from the cross-validation reference |
| `zero_sd_cv` | `"zero"`: voxels whose reference betas have zero variance carry no cross-validation weight | `"python"`: treated inconsistently across candidates by an in-place division helper |
| `singular` | `"error"`, as GLMsingle | `"pinv"` gives a minimum-norm solution |
| `frac_alpha` | `"fracridge"`: GLMsingle’s grid interpolation of the ridge penalty | `"exact"` solves the fraction equation exactly (experimental) |

The automatic ON-OFF R^2 threshold uses a deterministic Gaussian-mixture
fit with a small variance floor (as in GLMsingle’s MATLAB code). Pass
`brain_r2` and `pc_r2_cutoff` to fix the thresholds yourself.

`full_glmbadness = TRUE` fills the PC cross-validation diagnostic for
every voxel. By default only the voxels that decide the number of PCs
are evaluated, which gives identical estimates.

## How agreement with GLMsingle is checked

The package tests compare
[`glmsingle()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle.md)
with pinned Python GLMsingle on 11 simulated scenarios: default
settings, extra regressors, two sessions with grouped folds, unequal run
lengths, unrepeated conditions, rapid events, near-collinear regressors,
an all-zero voxel, a single ridge fraction, and no HRF library. Every
HRF choice, number of PCs and ridge fraction matches. Betas and R^2
agree to about 1e-6, the precision of GLMsingle’s single-precision
arithmetic. For rapid designs the agreement is limited by GLMsingle’s
float32 normal equations (relative error grows with the square of the
design’s condition number). Separate tests check every stage in double
precision against literal re-implementations of GLMsingle’s stacked
computations.

## Reference

Prince, J. S., Charest, I., Kurzawski, J. W., Pyles, J. A., Tarr, M. J.,
& Kay, K. N. (2022). Improving the accuracy of single-trial fMRI
response estimates using GLMsingle. *eLife*, 11, e77599.

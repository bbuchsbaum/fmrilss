# Voxel-wise HRF Modeling with fmrilss

## Why HRF Variability Matters

Vascular properties, neurovascular coupling, and acquisition protocols
all influence the hemodynamic response. Relying on a single canonical
HRF shape can bias trial-wise estimates, especially across brain regions
or clinical populations.

This vignette shows how to estimate voxel-specific HRFs and incorporate
them into LSS. You should be familiar with the core LSS workflow
([`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md))
and `fmrihrf` basics.

**Alternative approach:** For a library-constrained method that is
usually faster and more stable, see
[`vignette("sbhm")`](https://bbuchsbaum.github.io/fmrilss/articles/sbhm.md).

## Setup

``` r
library(fmrilss)
library(fmrihrf)
set.seed(123)
```

## Simulate Data with Variable HRFs

### Experiment parameters

You need a sampling frame, jittered onsets, and a small set of voxels
whose HRFs differ.

``` r
n_time <- 200; n_vox <- 5; TR <- 1.0
sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = TR)
grid   <- fmrihrf::samples(sframe, global = TRUE)
```

``` r
isi    <- runif(200, min = 3, max = 9)
onsets <- cumsum(c(10, isi))
onsets <- onsets[onsets < (n_time - 20)]
n_trials <- length(onsets)
```

### Voxel-specific HRFs

Each voxel gets a slightly shifted and scaled version of the canonical
HRF. This mimics spatial variation in vascular properties.

We create each voxel’s HRF by wrapping the canonical shape with a peak
shift and width scaling.

``` r
voxel_hrfs <- lapply(1:n_vox, function(v) {
  peak_shift  <- (v - 3) * 0.5        # -1 to +1 s
  width_scale <- 1 + (v - 3) * 0.1    #  0.8 to 1.2
  fmrihrf::HRF(
    fun   = function(t) HRF_SPMG1(t - peak_shift) * width_scale,
    name  = paste0("voxel_", v),
    span  = attr(HRF_SPMG1, "span"),
    nbasis = 1L
  )
})
```

### Generate signal and noise

For each voxel, you convolve the trial onsets with that voxel’s HRF,
scale by true betas, then add AR(1) noise.

``` r
true_betas <- matrix(rnorm(n_trials * n_vox, mean = 1, sd = 0.3),
                     nrow = n_trials, ncol = n_vox)
```

``` r
Y <- matrix(0, n_time, n_vox)
for (v in 1:n_vox) {
  rset <- fmrihrf::regressor_set(
    onsets = onsets, fac = factor(seq_len(n_trials)),
    hrf = voxel_hrfs[[v]], duration = 0, span = 30, summate = FALSE
  )
  Xv <- fmrihrf::evaluate(rset, grid = grid, precision = 0.1, method = "conv")
  Y[, v] <- as.matrix(Xv) %*% true_betas[, v]
}
```

``` r
noise_sd <- 0.5; ar_coef <- 0.3
for (v in 1:n_vox) {
  e <- rnorm(n_time, sd = noise_sd)
  noise <- as.numeric(stats::filter(e, filter = ar_coef, method = "recursive"))
  Y[, v] <- Y[, v] + noise
}
colnames(Y) <- paste0("V", 1:n_vox)
```

### Visualise the design

``` r
rset_vis <- fmrihrf::regressor_set(
  onsets = onsets, fac = factor(seq_len(n_trials)),
  hrf = HRF_SPMG1, duration = 0, span = 30, summate = FALSE)
X_vis <- as.matrix(fmrihrf::evaluate(rset_vis, grid = grid, precision = 0.1, method = "conv"))
image(seq_len(nrow(X_vis)), seq_len(ncol(X_vis)), X_vis,
      col = hcl.colors(64, "BluGrn"), xlab = "Time (TR)", ylab = "Trial",
      main = "Trial-wise design (jittered ISIs)")
```

![Trial-wise design
heatmap.](voxel-wise-hrf_files/figure-html/design-heatmap-1.png)

## Standard LSS with Canonical HRF

Standard LSS assumes every voxel shares the same canonical HRF. When
that assumption is wrong, you get biased betas.

``` r
rset_can <- fmrihrf::regressor_set(
  onsets = onsets, fac = factor(seq_len(n_trials)),
  hrf = HRF_SPMG1, duration = 0, span = 30, summate = FALSE)
X_can <- as.matrix(fmrihrf::evaluate(rset_can, grid = grid, precision = 0.1, method = "conv"))
standard_betas <- lss(Y, X_can, method = "r_optimized")
```

## Estimate Voxel-Specific HRFs

### Multi-basis GLM

The SPMG3 basis set includes the canonical HRF plus its temporal and
dispersion derivatives. Fitting a GLM with this basis set lets you
estimate how each voxel’s HRF deviates from canonical.

``` r
rset_mb <- fmrihrf::regressor_set(
  onsets = onsets, fac = factor(rep(1, n_trials)),
  hrf = HRF_SPMG3, duration = 0, span = 30, summate = TRUE)
X_mb <- as.matrix(fmrihrf::evaluate(rset_mb, grid = grid, precision = 0.1, method = "conv"))
```

``` r
hrf_weights <- sapply(1:n_vox, function(v) coef(lm(Y[, v] ~ X_mb - 1)))
cat("Basis weights (3 x", n_vox, "voxels):\n")
#> Basis weights (3 x 5 voxels):
print(round(hrf_weights, 2))
#>       [,1] [,2]  [,3]  [,4]  [,5]
#> X_mb1 0.74 0.85  1.04  1.16  1.17
#> X_mb2 0.46 0.49  0.32 -0.46 -0.81
#> X_mb3 0.14 0.23 -0.31  0.24  0.30
```

The first row captures the canonical amplitude; rows 2–3 capture latency
and width shifts.

## Apply Voxel-Specific HRFs in LSS

With basis weights in hand, you can build a voxel-specific design by
weighting the three basis columns for each trial. Notice that this is a
per-voxel loop: each voxel gets its own tailored design matrix.

``` r
voxel_betas <- matrix(NA, n_trials, n_vox)
```

For each voxel, build a tailored design matrix by weighting the three
basis columns with that voxel’s estimated HRF coefficients, then run
LSS.

``` r
for (v in 1:n_vox) {
  Xv <- X_can * 0
  for (tr in 1:n_trials) {
    rset_tr <- fmrihrf::regressor_set(onsets = onsets[tr], fac = factor(1),
      hrf = HRF_SPMG3, duration = 0, span = 30, summate = FALSE)
    cols <- as.matrix(fmrihrf::evaluate(rset_tr, grid = grid, precision = 0.1, method = "conv"))
    Xv[, tr] <- cols %*% hrf_weights[, v]
  }
  voxel_betas[, v] <- lss(Y[, v, drop = FALSE], Xv, method = "r_optimized")
}
```

For production analyses,
[`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md)
wraps this loop with optional C++ acceleration. See
[`?lss_with_hrf`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md).

## OASIS Alternative

OASIS handles multi-basis HRFs and LSS in a single call. You pass
`HRF_SPMG3` and it solves for all basis components simultaneously, with
optional ridge regularization for stability.

``` r
oasis_betas <- lss(
  Y, X = NULL, method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(onsets = onsets, hrf = HRF_SPMG3, span = 30)
    ),
    ridge_mode = "fractional", ridge_x = 0.01, ridge_b = 0.01
  )
)
oasis_canonical <- oasis_betas[seq(1, nrow(oasis_betas), by = 3), ]
```

See
[`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
for ridge tuning and standard-error computation.

## Compare Methods

### Compute accuracy metrics

``` r
cors  <- c(Standard  = cor(as.vector(standard_betas), as.vector(true_betas)),
           VoxelHRF  = cor(as.vector(voxel_betas),    as.vector(true_betas)),
           OASIS     = cor(as.vector(oasis_canonical), as.vector(true_betas)))
rmses <- c(Standard  = sqrt(mean((standard_betas - true_betas)^2)),
           VoxelHRF  = sqrt(mean((voxel_betas    - true_betas)^2)),
           OASIS     = sqrt(mean((oasis_canonical - true_betas)^2)))
comparison_summary <- data.frame(
  Correlation = round(cors, 3),
  RMSE = round(rmses, 3)
)
comparison_summary
#>          Correlation  RMSE
#> Standard       0.733 0.279
#> VoxelHRF       0.802 0.214
#> OASIS          0.697 0.322
stopifnot(all(is.finite(cors)), all(is.finite(rmses)))
```

These metrics give you a direct accuracy check on the simulated data
rather than relying on a visual impression alone.

### Scatter plots

Points closer to the diagonal mean better recovery.

``` r
rng <- range(true_betas)
cls <- rep(1:n_vox, each = n_trials)
```

Each panel shows estimated versus true betas; colour distinguishes
voxels.

``` r
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (i in seq_along(cors)) {
  est <- list(standard_betas, voxel_betas, oasis_canonical)[[i]]
  plot(true_betas, est, pch = 19, col = cls, xlim = rng, ylim = rng,
       xlab = "True", ylab = "Estimated",
       main = paste0(names(cors)[i], " (r=", round(cors[i], 2), ")"))
  abline(0, 1, lty = 2, col = "gray")
}
```

![Scatter plots comparing true vs estimated betas for each
method.](voxel-wise-hrf_files/figure-html/compare-plots-1.png)

## Next Steps

**When to use voxel-wise HRFs.** You should consider this approach when
regions differ in vascular architecture (motor vs. visual cortex), when
studying populations with altered neurovascular coupling (aging,
clinical), or when high-resolution acquisitions expose laminar-level
variation.

**Choosing a method.** Standard LSS with a canonical HRF is sufficient
for homogeneous responses. Voxel-wise HRF LSS improves accuracy when HRF
heterogeneity is expected. OASIS is preferred for rapid-event designs or
when you want a single-step solve with HRF flexibility and ridge
control.

**Further reading:**

- [`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md)
  – LSS basics
- [`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
  – OASIS solver with ridge, multi-basis HRFs, and standard errors
- [`vignette("sbhm")`](https://bbuchsbaum.github.io/fmrilss/articles/sbhm.md)
  – Shared-Basis HRF Matching for efficient voxel-specific HRFs
- [`?estimate_voxel_hrf`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md)
  and
  [`?lss_with_hrf`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md)
  – production-ready voxel-wise HRF workflow

&nbsp;

    #> R version 4.5.3 (2026-03-11)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] fmrihrf_0.3.0 fmrilss_0.1.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] Matrix_1.7-4        gtable_0.3.6        jsonlite_2.0.0     
    #>  [4] dplyr_1.2.1         compiler_4.5.3      Rcpp_1.1.1         
    #>  [7] tidyselect_1.2.1    assertthat_0.2.1    jquerylib_0.1.4    
    #> [10] splines_4.5.3       systemfonts_1.3.2   scales_1.4.0       
    #> [13] textshaping_1.0.5   uuid_1.2-2          yaml_2.3.12        
    #> [16] fastmap_1.2.0       lattice_0.22-9      ggplot2_4.0.2      
    #> [19] R6_2.6.1            generics_0.1.4      knitr_1.51         
    #> [22] htmlwidgets_1.6.4   tibble_3.3.1        desc_1.4.3         
    #> [25] bslib_0.10.0        pillar_1.11.1       RColorBrewer_1.1-3 
    #> [28] rlang_1.2.0         cachem_1.1.0        xfun_0.57          
    #> [31] fs_2.0.1            sass_0.4.10         S7_0.2.1           
    #> [34] otel_0.2.0          memoise_2.0.1       cli_3.6.6          
    #> [37] pkgdown_2.2.0       magrittr_2.0.5      digest_0.6.39      
    #> [40] grid_4.5.3          bigmemory.sri_0.1.8 bigmemory_4.6.4    
    #> [43] lifecycle_1.0.5     vctrs_0.7.2         evaluate_1.0.5     
    #> [46] glue_1.8.0          numDeriv_2016.8-1.1 fmriAR_0.3.1       
    #> [49] farver_2.1.2        ragg_1.5.2          purrr_1.2.2        
    #> [52] rmarkdown_2.31      albersdown_1.0.0    tools_4.5.3        
    #> [55] pkgconfig_2.0.3     htmltools_0.5.9

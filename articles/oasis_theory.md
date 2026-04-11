# OASIS Theory: Algebra and Implementation Details

## Motivation

Optimized Analytic Single-pass Inverse Solution (OASIS) extends Least
Squares Separate (LSS) estimation through algebraic reformulation that
enables single-pass computation of all trial estimates. For a
single-basis HRF (K = 1) without ridge, OASIS reduces exactly to the
closed-form LSS solution; the per‑trial 2x2 normal equations are the
same. The value of OASIS is in batching those solves efficiently and
generalizing the same algebra to multi‑basis HRFs (2Kx2K) with optional
ridge and diagnostics. This document provides the mathematical
foundation and implementation details.

Read this after
[`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
if you want the linear algebra behind the user-facing API. The code
chunks below are there to make the main scaling claims executable rather
than purely narrative.

### Prerequisites

This vignette assumes familiarity with: - QR decomposition and
orthogonal projection matrices - Ridge regression and regularization -
Matrix calculus and linear algebra - The standard LSS formulation

### Computational Intuition

Standard LSS requires N separate GLM fits for N trials, each
involving: 1. Matrix assembly: O(T²) operations 2. QR decomposition:
O(T³) operations 3. Back-substitution: O(T²) operations

Total complexity: O(NT³) for N trials

OASIS recognizes that these N models share substantial structure. By
factoring out common computations, OASIS reduces complexity to: 1.
Single QR decomposition: O(T³) 2. Shared projections: O(NT²) 3.
Per-trial solutions: O(N)

Total complexity: O(T³ + NT²), a significant reduction when N is large.

### Visual Comparison of Computational Scaling

``` r
# Demonstrate computational scaling
N_trials <- c(10, 50, 100, 200, 500, 1000)
T_points <- 200  # Fixed number of timepoints

# Simplified complexity models (arbitrary units)
classical_ops <- N_trials * T_points^3 / 1e6  # O(NT³)
oasis_ops <- (T_points^3 + N_trials * T_points^2) / 1e6  # O(T³ + NT²)

# Create comparison plot
plot(N_trials, classical_ops, type='l', col='red', lwd=2,
     xlab='Number of Trials', ylab='Computational Operations (millions)',
     main='Computational Complexity: Classical LSS vs OASIS',
     ylim=c(0, max(classical_ops)))
lines(N_trials, oasis_ops, col='blue', lwd=2)
legend('topleft', c('Classical LSS', 'OASIS'),
       col=c('red', 'blue'), lwd=2, bty='n')

# Add shaded region showing computational savings
polygon(c(N_trials, rev(N_trials)),
        c(classical_ops, rev(oasis_ops)),
        col=rgb(0.2, 0.8, 0.2, 0.3), border=NA)
text(500, mean(c(classical_ops[4], oasis_ops[4])),
     'Computational\nSavings', col='darkgreen')
```

![Line plot showing classical LSS scaling linearly with trials while
OASIS remains nearly
flat.](oasis_theory_files/figure-html/complexity-comparison-1.png)

Computational complexity: Classical LSS vs OASIS

``` r
complexity_summary <- data.frame(
  Trials = N_trials,
  Classical = classical_ops,
  OASIS = oasis_ops,
  SpeedupRatio = classical_ops / oasis_ops
)
complexity_summary
#>   Trials Classical OASIS SpeedupRatio
#> 1     10        80   8.4      9.52381
#> 2     50       400  10.0     40.00000
#> 3    100       800  12.0     66.66667
#> 4    200      1600  16.0    100.00000
#> 5    500      4000  28.0    142.85714
#> 6   1000      8000  48.0    166.66667
```

Code references point to `R/oasis_glue.R` and `src/oasis_core.cpp`
implementations.

Notation used throughout:

- $Y \in {\mathbb{R}}^{T \times V}$: voxel data ($T$ time points, $V$
  voxels)
- $X = \left\lbrack x_{1},\ldots,x_{N} \right\rbrack \in {\mathbb{R}}^{T \times N}$:
  trial regressors for one condition
- $Z \in {\mathbb{R}}^{T \times K_{z}}$: nuisance/experimental
  regressors shared across trials
- $R = I - QQ^{T}$: orthogonal projector removing nuisance effects ($Q$
  comes from QR factorisation of
  $\left\lbrack Z,\text{others} \right\rbrack$)
- Inner products are denoted $\langle a,b\rangle = a^{T}b$

We first treat the single-basis case (one regressor per trial) before
generalizing to multi-basis HRFs.

## Classical LSS Recap

Classical LSS fits, for each trial $j$, a GLM with design
$\left\lbrack x_{j},b_{j},Z \right\rbrack$, where
$b_{j} = \sum_{i \neq j}x_{i}$. Solving each model independently costs
$\mathcal{O}(N)$ QR factorizations. Algebraically, the trial-specific
beta can be expressed as

$${\widehat{\beta}}_{j} = \frac{\langle Rx_{j},RY\rangle - \frac{\langle Rx_{j},Rb_{j}\rangle}{\parallel Rb_{j} \parallel^{2}}\langle Rb_{j},RY\rangle}{\parallel Rx_{j} \parallel^{2} - \frac{\langle Rx_{j},Rb_{j}\rangle^{2}}{\parallel Rb_{j} \parallel^{2}}}.$$

OASIS extracts and reuses the common computational components
(projections, norms, cross-products) across all trials, computing each
only once.

## Single-Basis OASIS Algebra

After residualising against nuisance regressors we define:

- $a_{j} = Rx_{j}$
- $s = \sum_{j = 1}^{N}a_{j}$
- $d_{j} = \parallel a_{j} \parallel^{2}$
- $\alpha_{j} = \langle a_{j},s - a_{j}\rangle$
- $s_{j} = \parallel s - a_{j} \parallel^{2}$

Let $n_{jv} = \langle a_{j},RY_{\cdot v}\rangle$ and
$m_{v} = \langle s,RY_{\cdot v}\rangle$. The pair
$\left( \beta_{j},\gamma_{j} \right)$ solving the 2x2 system for trial
$j$ and voxel $v$ is obtained from

$$G_{j}\begin{bmatrix}
\beta_{jv} \\
\gamma_{jv}
\end{bmatrix} = \begin{bmatrix}
n_{jv} \\
{m_{v} - n_{jv}}
\end{bmatrix},\quad G_{j} = \begin{bmatrix}
{d_{j} + \lambda_{x}} & \alpha_{j} \\
\alpha_{j} & {s_{j} + \lambda_{b}}
\end{bmatrix},$$

with ridge penalties $\lambda_{x},\lambda_{b} \geq 0$. The inverse of
$G_{j}$ is analytic, so  
$$\beta_{jv} = \frac{\left( s_{j} + \lambda_{b} \right)n_{jv} - \alpha_{j}\left( m_{v} - n_{jv} \right)}{\left( d_{j} + \lambda_{x} \right)\left( s_{j} + \lambda_{b} \right) - \alpha_{j}^{2}}.$$

This is exactly what `oasis_betas_closed_form()` implements (C++ file
`src/oasis_core.cpp`). The precomputation step
`oasis_precompute_design()` produces $a_{j},s,d_{j},\alpha_{j},s_{j}$
once, while `oasis_AtY_SY_blocked()` streams through voxels to obtain
$n_{jv}$ and $m_{v}$.

### Fractional Ridge

`oasis$ridge_mode = "fractional"` sets
$\lambda_{x} = \eta_{x} \cdot \bar{d}$ and
$\lambda_{b} = \eta_{b} \cdot \bar{s}$, where $\bar{d}$ and $\bar{s}$
are means of $d_{j}$ and $s_{j}$. The helper `.oasis_resolve_ridge()`
implements this scaling. Absolute ridge uses the supplied values
directly.

### Standard Errors

Given $G_{j}^{- 1}$ and residual norm $\parallel RY \parallel^{2}$, the
variance of $\beta_{jv}$ is

$$\operatorname{Var}\left( {\widehat{\beta}}_{jv} \right) = \sigma_{jv}^{2}\left( G_{j}^{- 1} \right)_{11},\quad\sigma_{jv}^{2} = \frac{\text{SSE}_{jv}}{\text{dof}},$$

with

$$\text{SSE}_{jv} = \parallel RY_{\cdot v} \parallel^{2} - 2\left( \beta_{jv}n_{jv} + \gamma_{jv}\left( m_{v} - n_{jv} \right) \right) + d_{j}\beta_{jv}^{2} + s_{j}\gamma_{jv}^{2} + 2\alpha_{j}\beta_{jv}\gamma_{jv}.$$

`.oasis_se_from_norms()` realises this computation, reusing $n_{jv}$,
$m_{v}$ and the cached design scalars.

## Multi-Basis Extension

When the HRF contributes $K > 1$ basis functions, each trial has columns
$A_{j} \in {\mathbb{R}}^{T \times K}$. Define

- $S = \sum_{j}A_{j}$
- $D_{j} = A_{j}^{T}A_{j}$
- $C_{j} = A_{j}^{T}\left( S - A_{j} \right)$
- $E_{j} = \left( S - A_{j} \right)^{T}\left( S - A_{j} \right)$

Per voxel we need $N1 = A^{T}RY$ (stacked $N$ blocks of size $K$) and
$SY = S^{T}RY$. The block system is

$$\begin{bmatrix}
{D_{j} + \lambda_{x}I} & C_{j} \\
C_{j}^{T} & {E_{j} + \lambda_{b}I}
\end{bmatrix}\begin{bmatrix}
B_{jv} \\
\Gamma_{jv}
\end{bmatrix} = \begin{bmatrix}
{N1_{jv}} \\
{SY_{v} - N1_{jv}}
\end{bmatrix},$$

where $B_{jv} \in {\mathbb{R}}^{K}$. `oasisk_betas()` solves this 2Kx2K
system via Cholesky factorisation. Ridge again adds $\lambda_{x}I$ and
$\lambda_{b}I$ to the block diagonals. Compared to the single-basis
path, only the shapes of the cached matrices differ; the solve is still
analytic per trial/voxel block.

The companion `oasisk_betas_se()` extends the SSE/variance calculation
to the multi-basis case, using the same building blocks.

## HRF-Aware Design Construction

OASIS can construct $X$ on the fly from event specifications.
`.oasis_build_X_from_events()` uses
[`fmrihrf::regressor_set()`](https://bbuchsbaum.github.io/fmrihrf/reference/regressor_set.html)
to generate trial-wise columns (and optional other-condition aggregates)
given:

- `cond$onsets`: per-trial onset times
- `cond$hrf`: HRF object (canonical, FIR, multi-basis, user-defined)
- `cond$span`, `precision`, `method`: convolution controls

This design is then residualised against nuisance regressors and fed
into the algebra above. Because the HRF definition enters directly,
switching HRFs or running grid searches automatically regenerates a
matching design. When you provide an explicit `X`, OASIS skips this step
and assumes you have already encoded the HRF in the matrix.

## AR(1) Whitening

`oasis$whiten = "ar1"` estimates a common AR(1) coefficient from
residualised data. `.oasis_ar1_whitener()` computes $\rho$ and applies
Toeplitz-safe differencing:

$${\widetilde{y}}_{t} = \begin{cases}
{\sqrt{1 - \rho^{2}}y_{1}} & {t = 1,} \\
{y_{t} - \rho y_{t - 1}} & {t > 1.}
\end{cases}$$

The same transformation is applied to $X$ and nuisance regressors before
the standard OASIS algebra runs. Whitening preserves the single-pass
benefits because the transformed data are treated exactly like the
original inputs.

## Diagnostics Output

When `oasis$return_diag = TRUE`, OASIS returns the precomputed design
scalars:

- Single-basis: $d_{j},\alpha_{j},s_{j}$ (from
  `oasis_precompute_design()`)
- Multi-basis: $D_{j},C_{j},E_{j}$ (from `oasisk_precompute_design()`)

These matrices are useful for checking trial collinearity, energy, and
the effect of ridge scaling.

## Algorithm Summary

Putting everything together, the single-basis solver proceeds as
follows:

1.  Residualise $Y$ and $X$ against nuisance regressors, optionally with
    whitening.
2.  Compute $a_{j},s,d_{j},\alpha_{j},s_{j}$
    (`oasis_precompute_design`).
3.  Stream through voxels in blocks, forming $N_{Y} = A^{T}RY$ and
    $S_{Y} = s^{T}RY$ (`oasis_AtY_SY_blocked`).
4.  Apply ridge scaling (absolute or fractional) to obtain
    $\lambda_{x},\lambda_{b}$.
5.  For each trial, evaluate the closed-form $\beta_{jv}$ (and
    $\gamma_{jv}$ if SEs requested).
6.  Optionally compute SEs and diagnostics.

The multi-basis path swaps steps 2–5 for their block equivalents. In
both cases, the cost is dominated by the single projection of $Y$ and
the matrix–vector multiplies in step 3, giving $\mathcal{O}(TV)$
complexity with a small trial-dependent overhead.

## Complexity and Memory

- Projection / whitening: $\mathcal{O}\left( TVK_{z} \right)$
  arithmetic, $\mathcal{O}\left( TK_{z} \right)$ memory for confounds
- Precomputation: $\mathcal{O}(TN)$
- Products (blocked): $\mathcal{O}(TV)$ with block size tuning
- Closed-form solves: $\mathcal{O}(NV)$ with negligible constants (2x2
  or 2Kx2K systems)

Compared to classical LSS ($N$ separate regressions), OASIS shaves off
repeated projections and linear solves, yielding substantial speedups
when $N$ or $V$ is large.

## Next Steps

- [`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
  — practical OASIS usage with ridge, multi-basis HRFs, and standard
  errors
- [`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md)
  — foundational LSS concepts
- [`vignette("sbhm")`](https://bbuchsbaum.github.io/fmrilss/articles/sbhm.md)
  — library-constrained voxel-specific HRFs

## References

- Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012).
  Deconvolving BOLD activation in event-related designs for multivoxel
  pattern classification analyses. *NeuroImage*, 59(3), 2636–2643.
- fmrilss source files `R/oasis_glue.R` and `src/oasis_core.cpp` (for
  implementation alignment).

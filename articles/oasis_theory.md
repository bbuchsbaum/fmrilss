# OASIS Theory: What Is Reused, Solved, and Returned

OASIS computes least squares separate (LSS) estimates by reusing the
calculations shared across trial models, then solving the small
trial-specific systems in batches. With zero ridge penalties, it gives
the ordinary LSS estimator. Nonzero penalties give a regularized LSS
estimator.

This article defines what OASIS estimates, derives the linear systems it
solves for one and for several HRF basis functions per trial, and maps
each mathematical object to the current implementation. Read it after
[`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md)
and
[`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md);
those articles introduce the practical meaning of LSS and the
user-facing OASIS API.

## Scope and notation

The dimensions below are used throughout.

| Symbol | Meaning                                        |
|:------:|:-----------------------------------------------|
|   T    | time points                                    |
|   V    | response columns or voxels                     |
|   N    | trials in the condition being estimated        |
|   K    | HRF basis columns per trial                    |
|  K_z   | columns supplied in the common nuisance design |
|   p    | rank of the common nuisance design; p \<= K_z  |
|   B    | voxel block size used for data products        |

Symbols and dimensions used in the OASIS derivation. {.table}

Let $`Y \in \mathbb{R}^{T \times V}`$ be the data and
$`C \in \mathbb{R}^{T \times K_z}`$ the common design: intercepts,
drifts, nuisance variables, and aggregates for conditions not being
estimated. Only the column space of $`C`$ matters. If
$`Q \in \mathbb{R}^{T \times p}`$ is an orthonormal basis for that
space, then the residualizing projection

``` math
R = I_T - QQ^\mathsf{T}
```

removes from any vector its component in that space. The implementation
rank-reduces $`C`$ before passing it to the compiled kernels, so
redundant nuisance columns do not add arbitrary projection directions.

For trial $`j`$, let $`A_j \in \mathbb{R}^{T \times K}`$ contain its HRF
basis columns and define

``` math
B_j = \sum_{i \ne j} A_i.
```

so $`B_j`$ is the sum of the basis columns of all other trials. The
Frisch–Waugh–Lovell theorem says that the coefficients on $`A_j`$ and
$`B_j`$ are unchanged if $`C`$ is dropped from the model and $`Y`$,
$`A_j`$, and $`B_j`$ are residualized by $`R`$ instead. After this
reduction, the trial-wise model is

``` math
RY = RA_j\,\boldsymbol{\beta}_j + RB_j\,\boldsymbol{\gamma}_j + R\varepsilon_j,
```

where
$`\boldsymbol{\beta}_j,\boldsymbol{\gamma}_j \in \mathbb{R}^{K \times V}`$.
OASIS returns $`\boldsymbol{\beta}_j`$. For $`K=1`$, each trial
contributes one scalar coefficient per voxel. For $`K>1`$, each trial
contributes $`K`$ rows of basis coefficients. Converting those
coefficients to a single response amplitude requires a separate
definition and normalization.

## One-basis system

For $`K=1`$, $`A_j`$ is a single column $`x_j`$. Write $`a_j=Rx_j`$,
$`b_j=\sum_{i\ne j}a_i`$, and $`s=\sum_i a_i`$. The design scalars,
which depend on the design but not on $`Y`$ and are therefore computed
once for all voxels, are

``` math
d_j=a_j^\mathsf{T}a_j,\qquad
c_j=a_j^\mathsf{T}b_j,\qquad
e_j=b_j^\mathsf{T}b_j.
```

The C++ source, and the diagnostics returned with `return_diag = TRUE`,
name these quantities `d`, `alpha`, and `s`, respectively; note that the
code’s `s` is $`e_j`$, not the summed column $`s`$. For all voxels at
once, define $`n_{1j}=a_j^\mathsf{T}RY`$ and
$`n_{2j}=b_j^\mathsf{T}RY=s^\mathsf{T}RY-n_{1j}`$. With ridge penalties
$`\lambda_x`$ on the trial coefficient and $`\lambda_b`$ on the
other-trials coefficient (both zero for ordinary LSS; see [Ridge changes
the estimator](#ridge-changes-the-estimator)), the penalized normal
equations are

``` math
\begin{bmatrix}
d_j+\lambda_x & c_j \\
c_j & e_j+\lambda_b
\end{bmatrix}
\begin{bmatrix}\beta_j\\\gamma_j\end{bmatrix}
=
\begin{bmatrix}n_{1j}\\n_{2j}\end{bmatrix}.
```

Here each right-hand-side row has length $`V`$. Inverting the
$`2\times2`$ matrix gives

``` math
\widehat\beta_j =
\frac{(e_j+\lambda_b)n_{1j}-c_jn_{2j}}
{(d_j+\lambda_x)(e_j+\lambda_b)-c_j^2}.
```

When both penalties are zero and the two-column trial model has full
rank, this is exactly the ordinary LSS coefficient. A rank-deficient
unpenalized model is not identifiable, and the public OASIS path,
`lss(..., method = "oasis")`, stops with an error rather than silently
adding a small diagonal jitter. This $`2\times2`$ form assumes $`N>1`$;
with one trial there is no “other trials” column, and OASIS solves the
one-column target model instead.

## Multi-basis system

For $`K>1`$, define the residualized blocks $`\widetilde A_j=RA_j`$,
$`\widetilde S=\sum_j\widetilde A_j`$, and
$`\widetilde B_j=\widetilde S-\widetilde A_j`$. The cached $`K\times K`$
blocks are

``` math
D_j=\widetilde A_j^\mathsf{T}\widetilde A_j,\qquad
C_j=\widetilde A_j^\mathsf{T}\widetilde B_j,\qquad
E_j=\widetilde B_j^\mathsf{T}\widetilde B_j.
```

With $`N_{1j}=\widetilde A_j^\mathsf{T}RY`$ and
$`N_{2j}=\widetilde S^\mathsf{T}RY-N_{1j}`$, OASIS solves

``` math
\begin{bmatrix}
D_j+\lambda_xI_K & C_j \\
C_j^\mathsf{T} & E_j+\lambda_bI_K
\end{bmatrix}
\begin{bmatrix}\boldsymbol{\beta}_j\\\boldsymbol{\gamma}_j\end{bmatrix}
=
\begin{bmatrix}N_{1j}\\N_{2j}\end{bmatrix}.
```

The implementation factorizes one $`2K\times2K`$ Gram matrix per trial
and uses all $`V`$ voxel columns as right-hand sides. The unpenalized
public path requires every such Gram matrix to have rank $`2K`$. For
$`N=1`$, the corresponding target-only system is $`K\times K`$.

## An executable equality check

We check the implementation against independently assembled trial-wise
GLMs. The fixed example below compares every trial, basis coefficient,
and voxel for both $`K=1`$ and $`K=3`$.

Set both penalties to zero for ordinary LSS; OASIS uses fractional ridge
by default.

``` r

unpenalized <- oasis_options(
  ridge_mode = "absolute",
  ridge_x = 0,
  ridge_b = 0,
  return_se = TRUE
)
fit1 <- lss(Y1, X1, Z = Z, method = "oasis", oasis = unpenalized)
dim(fit1$beta)
#> [1] 6 5
```

| Basis_dimension | Maximum_beta_error | Maximum_SE_error |
|----------------:|-------------------:|-----------------:|
|               1 |                  0 |                0 |
|               3 |                  0 |                0 |

Maximum discrepancy from independently assembled trial-wise GLMs.
{.table}

The vignette build stops if either maximum discrepancy reaches
$`10^{-10}`$. These checks use a correctly specified fixed design with
common coefficients for the summed trial signal and independent Gaussian
errors. They verify numerical agreement with the corresponding GLMs in
these examples. Uncertainty for ridge, estimated whitening, and
data-adaptive HRF selection requires additional methods.

## Ridge changes the estimator

For each trial and voxel, the ridge estimate minimizes

``` math
\lVert RY-RA_j\boldsymbol{\beta}_j-RB_j\boldsymbol{\gamma}_j\rVert_F^2
+\lambda_x\lVert\boldsymbol{\beta}_j\rVert_F^2
+\lambda_b\lVert\boldsymbol{\gamma}_j\rVert_F^2.
```

With `ridge_mode = "absolute"`, `ridge_x` and `ridge_b` are the
penalties in that expression. With `ridge_mode = "fractional"`, the
implementation multiplies them by mean residualized design energies,
where a column’s energy is its squared norm after residualization:

- for $`K=1`$, the means of $`d_j`$ and $`e_j`$;
- for $`K>1`$, the mean diagonal entries of $`D_j`$ and $`E_j`$,
  averaged over trials.

Because of this scaling, multiplying the whole trial design by a common
constant leaves the intended relative penalty unchanged. Fractional
scaling does not make arbitrary, separate rescalings of individual
multi-basis columns equivalent. Ridge can stabilize a poorly conditioned
system, but its coefficients are penalized estimates rather than
ordinary LSS coefficients.

The package default is fractional ridge with `ridge_x = ridge_b = 0.05`.
Request zero penalties explicitly when exact unpenalized LSS is the
target. The public path checks that each penalized multi-trial Gram
matrix is positive definite before solving it. When $`N=1`$, there is no
identifiable “other trials” coefficient, so OASIS solves only the target
block; `ridge_b` is then irrelevant.

## What the standard errors mean

The model-based standard errors below apply under the following
conditions:

- `ridge_x = ridge_b = 0`;
- for $`N>1`$, the trial-specific $`2K`$-column model is full rank; for
  $`N=1`$, the target $`K`$-column model is full rank;
- prewhitening is not estimated in the same call;
- the chosen HRF design is treated as fixed;
- residual degrees of freedom are positive;
- conditional errors have spherical covariance,
  $`\operatorname{Var}(\varepsilon_{\cdot v}\mid A_j,B_j,C)=\sigma_{jv}^2I_T`$
  for each voxel $`v`$.

The API checks the rank, penalty, and fitting-mode restrictions. The
error covariance assumption must be assessed for the data; the function
cannot verify it from the design alone.

For trial $`j`$, let $`G_j`$ be its identifiable unpenalized model Gram
matrix: $`2K\times2K`$ when $`N>1`$ and $`K\times K`$ when $`N=1`$. The
implementation computes

``` math
\widehat\sigma^2_{jv}=\frac{\operatorname{SSE}_{jv}}
{T-p-\operatorname{rank}(G_j)}
```

and takes the square roots of the relevant diagonal entries of
$`\widehat\sigma^2_{jv}G_j^{-1}`$. Thus the reported values are the same
conditional, homoskedastic-model standard errors as in the corresponding
full GLM under uncorrelated, constant-variance temporal errors.
Gaussianity is an additional requirement for exact finite-sample t
inference. The values are not ridge standard errors, not
heteroskedasticity- or autocorrelation-robust errors, and not calibrated
for uncertainty introduced by estimating a whitening model or selecting
an HRF from the same data.
[`oasis_options()`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
rejects ridge together with `return_se = TRUE`; the backend also rejects
estimated prewhitening and voxel-adaptive HRF modes for this request.

## Prewhitening is set by `prewhiten`, not `oasis`

Temporal whitening is handled by the top-level `prewhiten` argument and
the shared `fmriAR` integration. The same estimated linear
transformation is applied to $`Y`$, the trial design, and the common
design before the OASIS algebra runs.

``` r

fit_ar1 <- lss(
  Y,
  X,
  Z = Z,
  method = "oasis",
  oasis = oasis_options(),
  prewhiten = prewhiten_options(method = "ar", p = 1)
)
```

For multiple runs, supply `pooling = "run"` and the run labels. The
legacy `oasis$whiten` field is deprecated and ignored with a warning;
use `prewhiten` instead. Coefficients remain available after estimated
whitening, but `return_se = TRUE` fails closed: the call stops with an
error instead of returning standard errors, because uncertainty for
feasible GLS (GLS with an estimated noise covariance) is not calibrated
here.

## Design construction and trial/basis identity

There are two supported routes into the same solver:

1.  `X` supplies an already convolved trial design.
2.  `oasis$design_spec` asks the backend to construct one from event
    onsets with `fmrihrf`.

The event route (2) preserves event duration and amplitude and adds
aggregate regressors for other conditions to the common nuisance span.
In multi-run low-level design specifications, onsets are global seconds
rather than within-run times;
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
with a
[`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
is the safer run-aware interface.

A raw multi-basis `X` cannot reveal trial identity from dimensions
alone. It therefore requires explicit `K`, `ntrials`, and a
`trial_basis_map` identifying every uniquely named column. The output
contains $`NK`$ rows in canonical trial-by-basis order. Supplying an
incorrect grouping would change the model being fitted, so ambiguous
trial/basis mappings are rejected.

If `design_spec$hrf_grid` is supplied, the implementation uses the
observed data to choose a candidate HRF from the grid before fitting.
That is a model-selection step. Treat the resulting coefficients as
conditional on the selected design. Because ordinary post-selection
standard errors are not established by the fixed-design formula above,
this route rejects `return_se = TRUE`.

The vignette also checks output and diagnostic dimensions, preservation
of trial/basis identity after column permutations, and the expected
response to rescaling. It checks that unsupported standard-error
requests produce errors.

## Computational cost and memory

For multiple trials, each full LSS design has $`p+2K`$ effective
columns: $`p`$ common columns, $`K`$ target columns, and $`K`$
aggregate-other columns. The table below lists the leading costs for a
dense design with $`K_z\le T`$ and nuisance rank $`p`$.

| Stage                         |   K = 1    |     general K      |
|:------------------------------|:----------:|:------------------:|
| Rank-reveal common design     | O(T K_z^2) |     O(T K_z^2)     |
| Residualize all trial columns |  O(T p N)  |     O(T p N K)     |
| Build trial Gram terms        |   O(T N)   |     O(T N K^2)     |
| Residualize all voxel data    |  O(T p V)  |      O(T p V)      |
| Trial-data cross-products     |  O(T N V)  |     O(T N K V)     |
| Solve all trial systems       |   O(N V)   | O(N K^3 + N K^2 V) |

Leading arithmetic terms in the current dense OASIS implementation.
{.table}

The $`O(TNKV)`$ cross-product is usually the largest OASIS term when
both the trial and voxel counts are large. `block_cols = B` limits the
temporary voxel block to $`O(TB)`$, but it does not remove the
$`NK\times V`$ products or the $`NK\times V`$ output. Other persistent
storage includes $`O(Tp)`$ for the nuisance basis, $`O(TNK)`$ for the
residualized trial design, and $`O(NK^2)`$ for multi-basis Gram blocks.

A direct implementation that refits each trial model separately performs
$`N`$ factorizations of $`T\times(p+2K)`$ designs and repeatedly forms
trial-specific products with $`Y`$. OASIS reuses the nuisance
projection, aggregate design, Gram terms, and batched products. The time
saved depends on $`T,N,V,K,p`$, BLAS, memory bandwidth, block size, and
the implementation used for comparison; the operation counts do not
imply a fixed speedup.

## Implementation map and return values

These names are internal implementation landmarks, not additional public
APIs.

| Responsibility | Location |
|:---|:---|
| Validate options, assemble common span, dispatch | R/oasis_backend.R: .lss_oasis |
| Build event-based trial and other-condition designs | R/oasis_design.R: .oasis_build_X_from_events |
| Resolve absolute or fractional ridge | R/oasis_ridge_se.R: .oasis_resolve_ridge |
| One-basis design cache and blocked products | src/oasis_core.cpp: oasis_precompute_design, oasis_AtY_SY_blocked |
| One-basis coefficient solve | src/oasis_core.cpp: oasis_betas_closed_form |
| Multi-basis design cache and blocked products | src/oasis_core.cpp: oasisk_precompute_design, oasisk_products |
| Multi-basis coefficient and SE solves | src/oasis_core.cpp: oasisk_betas, oasisk_betas_se |

Internal implementation landmarks for each algebraic stage. {.table}

By default, the public result is an $`NK\times V`$ coefficient matrix
(an $`N\times V`$ matrix when $`K=1`$). Setting `return_se = TRUE` or
`return_diag = TRUE` changes the result to `list(beta, diag?, se?)`;
`beta` and `se` have the same $`NK\times V`$ shape. For $`K=1`$, `d`,
`alpha`, and `s` are length-$`N`$ diagnostic vectors. For $`K>1`$, `D`,
`C`, and `E` are $`K\times K\times N`$ arrays. These are unpenalized
cross-products on the residualized design scale, after any requested
whitening. They can expose low energy or collinearity, but they are not
a condition-number report and do not by themselves validate a chosen
ridge.

## Boundaries of the result

- Exact equality with classical LSS requires zero ridge, a fixed design,
  and full-rank trial models.
- Multi-basis output is a vector of basis coefficients per trial and
  voxel; any scalar amplitude needs an explicit, identified
  normalization rule.
- Estimated whitening changes the fitted space and adds uncertainty not
  covered by the conditional standard-error formula.
- HRF grid selection and voxel-adaptive HRFs are data-adaptive
  procedures; fixed-design uncertainty does not automatically survive
  selection.
- Blocking controls temporary memory. It cannot make storage smaller
  than the requested $`NK\times V`$ result.

## Next steps

- Continue to
  [`vignette("lss_with_fmridesign")`](https://bbuchsbaum.github.io/fmrilss/articles/lss_with_fmridesign.md)
  for run-aware event-table construction and trial/basis identity.
- Return to
  [`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
  for fitting and diagnostics, or
  [`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md)
  for the basic LSS workflow.
- Later articles cover explicitly normalized voxel-specific HRF shapes
  in
  [`vignette("voxel-wise-hrf")`](https://bbuchsbaum.github.io/fmrilss/articles/voxel-wise-hrf.md)
  and library-constrained shapes and trial coefficients in
  [`vignette("sbhm")`](https://bbuchsbaum.github.io/fmrilss/articles/sbhm.md).

## Reference

Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012).
Deconvolving BOLD activation in event-related designs for multivoxel
pattern classification analyses. *NeuroImage*, 59(3), 2636–2643.
<https://doi.org/10.1016/j.neuroimage.2011.08.076>

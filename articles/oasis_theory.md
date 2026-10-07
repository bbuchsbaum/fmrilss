# OASIS Theory: What Is Reused, Solved, and Returned

OASIS is an algebraic implementation of least-squares-separate (LSS). It
does not define a different unpenalized estimator. It identifies the
work that all trial-wise LSS models share, computes that work once, and
solves the small trial-specific systems in batches. With ridge
penalties, it instead computes a penalized LSS estimator.

This article establishes the estimand, derives the one- and multi-basis
systems, and maps each mathematical object to the current
implementation. Read it after
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
space, then

``` math
R = I_T - QQ^\mathsf{T}
```

removes it. The implementation rank-reduces $`C`$ before passing it to
the compiled kernels, so redundant nuisance columns do not add arbitrary
projection directions.

For trial $`j`$, let $`A_j \in \mathbb{R}^{T \times K}`$ contain its HRF
basis columns and define

``` math
B_j = \sum_{i \ne j} A_i.
```

After applying the Frisch–Waugh–Lovell reduction, the trial-wise model
is

``` math
RY = RA_j\,\boldsymbol{\beta}_j + RB_j\,\boldsymbol{\gamma}_j + R\varepsilon_j,
```

where
$`\boldsymbol{\beta}_j,\boldsymbol{\gamma}_j \in \mathbb{R}^{K \times V}`$.
OASIS returns $`\boldsymbol{\beta}_j`$. For $`K=1`$, each trial
contributes one scalar coefficient per voxel. For $`K>1`$, the $`K`$
returned rows are basis coefficients; they are not automatically a
single response-amplitude estimate.

## One-basis system

Write $`a_j=Rx_j`$, $`b_j=\sum_{i\ne j}a_i`$, and $`s=\sum_i a_i`$. The
reusable design scalars are

``` math
d_j=a_j^\mathsf{T}a_j,\qquad
c_j=a_j^\mathsf{T}b_j,\qquad
e_j=b_j^\mathsf{T}b_j.
```

The source names these quantities `d`, `alpha`, and `s`, respectively.
For all voxels at once, define $`n_{1j}=a_j^\mathsf{T}RY`$ and
$`n_{2j}=b_j^\mathsf{T}RY=s^\mathsf{T}RY-n_{1j}`$. The penalized normal
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
unpenalized model is not identifiable, and the public OASIS path fails
rather than silently adding jitter. This $`2\times2`$ form assumes
$`N>1`$; with one trial there is no “other trials” column, and OASIS
solves the one-column target model instead.

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

The strongest check of the algebra is not agreement between two OASIS
helpers; it is agreement with independently assembled trial-wise GLMs.
The following fixed example checks every trial, basis coefficient, and
voxel for both $`K=1`$ and $`K=3`$.

The exact public call for ordinary, unpenalized LSS is short; the zeros
are important because the OASIS default is fractionally penalized.

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

These checks use a correctly specified fixed design with common
coefficients for the summed trial signal and independent Gaussian
errors. They establish implementation equality with the corresponding
GLMs. They do not establish uncertainty for ridge, estimated whitening,
or data-adaptive HRF selection.

## Ridge changes the estimator

For each trial and voxel, ridge minimizes

``` math
\lVert RY-RA_j\boldsymbol{\beta}_j-RB_j\boldsymbol{\gamma}_j\rVert_F^2
+\lambda_x\lVert\boldsymbol{\beta}_j\rVert_F^2
+\lambda_b\lVert\boldsymbol{\gamma}_j\rVert_F^2.
```

With `ridge_mode = "absolute"`, `ridge_x` and `ridge_b` are the
penalties in that expression. With `ridge_mode = "fractional"`, the
implementation multiplies them by mean residualized design energies:

- for $`K=1`$, the means of $`d_j`$ and $`e_j`$;
- for $`K>1`$, the mean diagonal entries of $`D_j`$ and $`E_j`$,
  averaged over trials.

Fractional scaling preserves the intended relative penalty under a
common rescaling of the whole trial design. It does not make arbitrary,
separate rescalings of multi-basis columns equivalent. Ridge can
stabilize a poorly conditioned system, but its coefficients are
penalized estimates rather than ordinary LSS coefficients.

The package default is fractional ridge with `ridge_x = ridge_b = 0.05`.
Request zero penalties explicitly when exact unpenalized LSS is the
target. The public path checks that each penalized multi-trial Gram
matrix is positive definite before solving it. When $`N=1`$, there is no
identifiable “other trials” coefficient, so OASIS solves only the target
block; `ridge_b` is then irrelevant.

## What the standard errors mean

OASIS returns model-based standard errors only when all of the following
hold:

- `ridge_x = ridge_b = 0`;
- for $`N>1`$, the trial-specific $`2K`$-column model is full rank; for
  $`N=1`$, the target $`K`$-column model is full rank;
- prewhitening is not estimated in the same call;
- the chosen HRF design is treated as fixed;
- conditional errors have spherical covariance,
  $`\operatorname{Var}(\varepsilon_{\cdot v}\mid A_j,B_j,C)=\sigma_{jv}^2I_T`$
  for each voxel $`v`$.

For trial $`j`$, let $`G_j`$ be its identifiable unpenalized model Gram
matrix: $`2K\times2K`$ when $`N>1`$ and $`K\times K`$ when $`N=1`$. The
implementation computes

``` math
\widehat\sigma^2_{jv}=\frac{\operatorname{SSE}_{jv}}
{T-p-\operatorname{rank}(G_j)}
```

and takes the relevant diagonal of $`\widehat\sigma^2_{jv}G_j^{-1}`$.
Thus the reported values are the same conditional, homoskedastic-model
standard errors as in the corresponding full GLM under uncorrelated,
constant-variance temporal errors. Gaussianity is an additional
requirement for exact finite-sample t inference. The values are not
ridge standard errors, not heteroskedasticity- or autocorrelation-robust
errors, and not calibrated for uncertainty introduced by estimating a
whitening model or selecting an HRF from the same data.
[`oasis_options()`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
rejects ridge together with `return_se = TRUE`; the backend also rejects
estimated prewhitening and voxel-adaptive HRF modes for this request.

## Prewhitening belongs outside `oasis=`

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
legacy `oasis$whiten` field is deprecated, ignored, and emits a note; it
must not be used to describe the current algorithm. Coefficients remain
available after estimated whitening, but `return_se = TRUE` fails closed
because feasible-GLS uncertainty is not calibrated here.

## Design construction and identity

There are two supported routes into the same solver:

1.  `X` supplies an already convolved trial design.
2.  `oasis$design_spec` asks the backend to construct one with
    `fmrihrf`.

The event route preserves event duration and amplitude and adds
aggregates for other conditions to the common nuisance span. In
multi-run low-level design specifications, onsets are global seconds;
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
with a
[`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
is the safer run-aware interface.

A raw multi-basis `X` cannot reveal trial identity from dimensions
alone. It therefore requires explicit `K`, `ntrials`, and a
`trial_basis_map` identifying every uniquely named column. The output
contains $`NK`$ rows in canonical trial-by-basis order. Supplying an
incorrect grouping would change the estimand, so ambiguous inference is
rejected.

If `design_spec$hrf_grid` is supplied, the implementation chooses a
candidate using the observed data before fitting. That is a
model-selection step. Treat the resulting coefficients as conditional on
the selected design. Because ordinary post-selection standard errors are
not established by the fixed-design formula above, `return_se = TRUE`
fails closed on this route.

Hidden executable contracts also bind output shapes, diagnostic shapes,
permutation-safe multi-basis identity, event-built multi-basis
dimensions, scale equivariance, and the principal fail-closed inference
boundaries.

## Work and memory inventory

The former shorthand “one $`T\times T`$ QR instead of $`N`$ such QRs” is
not the relevant comparison. An LSS design has $`p+2K`$ effective
columns, not $`T`$ columns. The table below exposes every leading term
instead. It assumes $`K_z\le T`$, a dense design, and a nuisance rank
$`p`$.

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
projection, aggregate design, Gram terms, and batched products. That
structural reuse is exact; a universal wall-clock speedup is not.
Runtime depends on $`T,N,V,K,p`$, BLAS, memory bandwidth, block size,
and the competing implementation, so this article makes no synthetic
timing claim.

## Implementation map

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
- For a backward reference,
  [`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md)
  contains the practical fitting and diagnostics workflow, and
  [`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md)
  contains the foundational LSS introduction.
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

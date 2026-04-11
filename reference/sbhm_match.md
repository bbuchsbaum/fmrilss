# Match Voxels to Library HRFs in Shared Basis (SBHM)

Given per-voxel aggregate coefficients `beta_bar` in the shared basis
`B`, and library coordinates `A`, perform shape-only matching in a
whitened, L2-normalized coefficient space (cosine similarity).
Optionally apply a simple shrinkage of `beta_bar` towards a reference
coordinate before matching, and an orientation fix.

## Usage

``` r
sbhm_match(
  beta_bar,
  S,
  A,
  shrink = list(tau = 0, ref = NULL, snr = NULL),
  topK = 1,
  whiten = TRUE,
  sv_floor_rel = 1e-06,
  whiten_power = 1,
  orient_ref = TRUE
)
```

## Arguments

- beta_bar:

  Numeric matrix (r×V): per-voxel coefficients from a prepass GLM in the
  SBHM basis. Columns correspond to voxels.

- S:

  Numeric vector (length r): singular values of the library SVD.

- A:

  Numeric matrix (r×K): coordinates of library HRFs in SBHM basis.

- shrink:

  List with optional shrinkage options:

  - `tau` numeric \>=0: global strength (default 0, i.e., no shrinkage).

  - `ref` numeric length-r vector (alpha_ref) or NULL. If NULL, uses the
    mean of A columns. Shrinkage is: beta_bar \<- (1-lambda) beta_bar +
    lambda ref.

  - `snr` optional numeric length-V estimates; if provided, per-voxel
    lambda_v = tau/(snr_v + tau). Otherwise lambda = tau.

- topK:

  Integer, return top-K scores/weights if \>1 (default 1).

- whiten:

  Logical, divide coefficients by S before normalization (default TRUE).

- sv_floor_rel:

  Relative singular-value floor used when `whiten=TRUE` (default
  `1e-6`).

- whiten_power:

  Numeric in `[0, 1]` controlling whitening strength when `whiten=TRUE`.
  Uses division by `S^whiten_power` (`1` = full whitening, `0.5` =
  partial whitening, `0` = no whitening). Default `1`.

- orient_ref:

  Logical, flip beta_bar columns when their dot with `ref` is negative
  before matching (default TRUE).

## Value

A list with:

- `idx` length-V integer indices of best-matching library HRF (1..K)

- `margin` length-V numeric: score(top1) - score(top2)

- `alpha_hat` r×V matrix: the selected library coordinates (unwhitened,
  unnormalized)

- `scores` optional K×V cosine score matrix (returned when topK \> 1)

- `weights` optional top-K weights per voxel (when topK \> 1)

## Examples

``` r
if (FALSE) { # \dontrun{
  set.seed(42)
  r <- 4; K <- 12; V <- 3
  A <- matrix(rnorm(r*K), r, K)
  S <- seq(1, 0.2, length.out = r)
  alpha2 <- A[,2]
  beta_bar <- cbind(alpha2 + rnorm(r, sd = 0.1),
                    A[,7] + rnorm(r, sd = 0.1),
                    A[,10] + rnorm(r, sd = 0.1))
  m <- sbhm_match(beta_bar, S, A)
  m$idx
} # }
```

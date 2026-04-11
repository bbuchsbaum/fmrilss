# Prewhiten fMRI data using AR/ARMA models

Internal function providing a unified interface for prewhitening fMRI
data using the fmriAR package. Supports various AR/ARMA models with
flexible parameter estimation strategies.

## Usage

``` r
.prewhiten_data(Y, X = NULL, Z = NULL, Nuisance = NULL, prewhiten = list())
```

## Arguments

- Y:

  Numeric matrix (timepoints x voxels) of fMRI data

- X:

  Optional design matrix for trials (timepoints x trials)

- Z:

  Optional experimental design matrix (timepoints x regressors)

- Nuisance:

  Optional nuisance regressors (timepoints x nuisance)

- prewhiten:

  List of prewhitening options:

  method

  :   Character: "ar" (default), "arma", or "none"

  p

  :   Integer or "auto": AR order (default "auto")

  q

  :   Integer: MA order for ARMA (default 0)

  p_max

  :   Integer: Maximum AR order when p="auto" (default 6)

  pooling

  :   Character: "global" (default), "voxel", "run", or "parcel"

  runs

  :   Integer vector: Run identifiers for run-aware estimation

  parcels

  :   Integer vector: Parcel memberships for parcel-based pooling

  exact_first

  :   Character: "ar1" or "none" for exact AR(1) scaling

  compute_residuals

  :   Logical: Whether to compute residuals first (default TRUE)

## Value

List containing:

- Y_whitened:

  Whitened data matrix

- X_whitened:

  Whitened trial design matrix (if provided)

- Z_whitened:

  Whitened experimental design (if provided)

- Nuisance_whitened:

  Whitened nuisance regressors (if provided)

- whiten_plan:

  fmriAR plan object for diagnostics

- applied:

  Logical: Whether whitening was applied

# PRD: ITEM Utility Layer for fmrilss

Status: Draft  
Target release: next minor (`0.2.0`)

## Problem

`fmrilss` has strong LS-S/LS-A/OASIS trial-estimation workflows but no native ITEM preparation/fitting utilities for covariance-aware trial decoding. Users currently have to re-implement `X_t`, `U`, and fold-wise ITEM inversion externally.

## Goal

Add a compact ITEM helper layer that lets users move from trial-wise design and betas to ITEM-ready crossvalidated decoding inputs and scores.

## Non-goals

- Searchlight orchestration.
- Group-level inference/statistics maps.
- Replacing existing `rMVPA` decoding infrastructure.

## Users

- Methods developers prototyping covariance-aware decoding.
- Analysts using `fmridesign` + `fmrilss` who need ITEM-compatible outputs.

## Scope (MVP)

1. `item_build_design(...)`
- Build and validate trial-wise design metadata.
- Output includes `X_t`, `T_target`, `run_id`, `trial_info`.

2. `item_compute_u(X_t, V = NULL, v_type = c("cov", "precision"), ridge = 0, ...)`
- Compute `U = (X_t' V^{-1} X_t)^{-1}` using stable solver paths.
- `V = NULL` implies identity.
- `V` accepts dense matrix, sparse `Matrix`, or run-block list.
- Dense/sparse `V` must be `n_time x n_time`.

3. `item_fit(Gamma_train, T_train, U_train, ridge = 0)`
- Fit ITEM weights with GLS:
  `W_hat = (Gamma' U^{-1} Gamma)^{-1} Gamma' U^{-1} T_target`.

4. `item_predict(Gamma_test, W_hat)`
- Predict `T_hat = Gamma_test %*% W_hat`.

5. `item_cv(Gamma, T_target, U, run_id, mode = c("classification", "regression"), metric = NULL)`
- Deterministic LOSO CV (`sort(unique(run_id))`).
- Return per-fold and aggregate metrics.
- Classification tie handling fixed (`ties.method = "first"` + fixed class order).

6. `item_from_lsa(...)`
- Convenience wrapper: run LS-A (`lsa()` path), compute `U`, return ITEM bundle.

## Object contract

S3 class: `item_bundle`

Required fields:
- `Gamma`
- `X_t`
- `T_target`
- `U` and/or `U_by_run`
- `run_id`
- `meta`
- `diagnostics`

Optional fields:
- `C_transform` (for `X = X_t %*% C_transform`)
- `trial_id`
- `trial_hash`
- `trial_info`

## Functional requirements

- Strict dimension checks with actionable errors.
- Deterministic output for fixed seed/folds.
- Fold slicing helper (`item_slice_fold(bundle, test_run)`) used by `item_cv()`.
- Alignment guards across `Gamma`, `T_target`, `U`, `run_id`.
- Optional `trial_id` + hash checks.
- Support classification and regression targets via `T_target`.

## Non-functional requirements

- Numerical stability under collinearity.
- Accept base matrices and `Matrix` inputs.
- Reasonable performance for approximately 100 to 1000 trials.

## Validation plan

- Algebra tests: dimensions, symmetry/PSD expectations, fold submatrix correctness.
- Numerical stress tests: collinearity fallback path still returns finite outputs.
- Null simulations: classification near chance, regression near zero correlation.
- Signal simulations: performance improves with SNR.
- `R CMD check` clean at completion.

## Acceptance criteria

- New APIs exported and documented.
- Vignette section for ITEM preparation path.
- `item_cv()` runs end-to-end on simulated data and reproduces expected null behavior.
- No breaking changes to existing `lss()`/`lsa()`/`lss_design()` workflows.

## Risks

- Instability in `U` solves for highly collinear designs.
- Memory cost of dense `U` for large trial counts.
- Trial-order mismatches across `Gamma` and `T_target`.

## Mitigations

- Ridge + solver fallback (`chol -> svd -> pinv`) with warnings and diagnostics.
- Memory mode supporting run-block `U_by_run` lists.
- Explicit alignment assertions and optional hash-based trial checks.

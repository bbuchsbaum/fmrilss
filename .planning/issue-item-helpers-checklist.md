# ITEM Helpers: Implementation Checklist

## API and contracts

- [ ] Confirm final API signatures:
  - [ ] `item_build_design()`
  - [ ] `item_compute_u()`
  - [ ] `item_fit()`
  - [ ] `item_predict()`
  - [ ] `item_slice_fold()`
  - [ ] `item_cv()`
  - [ ] `item_from_lsa()`
- [ ] Confirm naming convention: `C_transform` vs `T_target`.
- [ ] Confirm `item_bundle` field contract and minimal required fields.

## Core implementation

- [ ] Add strict dimension/alignment checks with actionable errors.
- [ ] Support `V = NULL`, dense/sparse `V`, and run-block `V` list in `item_compute_u()`.
- [ ] Implement stable solver fallbacks (`chol`, `svd`, `pinv`) with diagnostics.
- [ ] Ensure no explicit inverse formulas in the implementation (use solves).
- [ ] Support both full `U` and run-block `U_by_run` memory mode.
- [ ] Add deterministic LOSO fold ordering (`sort(unique(run_id))`).
- [ ] Add deterministic classification tie handling.
- [ ] Add alignment guards with optional `trial_id` + hash checks.

## Testing

- [ ] Algebra tests:
  - [ ] Dimensions and symmetry
  - [ ] PSD expectations
  - [ ] Fold slicing correctness
- [ ] Numerical stress tests:
  - [ ] Collinearity fallback path activated
  - [ ] Finite outputs returned
- [ ] Simulation tests:
  - [ ] Null classification near chance
  - [ ] Null regression near zero correlation
  - [ ] Signal performance increases with SNR

## Documentation

- [ ] Roxygen docs for all new exported functions.
- [ ] `NAMESPACE` exports for all ITEM APIs.
- [ ] Vignette section: "Preparing ITEM inputs with fmrilss".
- [ ] Add NEWS entry for ITEM helpers.

## Release gates

- [ ] Run `devtools::test(filter = "item")` clean.
- [ ] Run full test suite clean.
- [ ] `R CMD check` clean.
- [ ] No regressions in existing `lss()` / `lsa()` / `lss_design()` workflows.

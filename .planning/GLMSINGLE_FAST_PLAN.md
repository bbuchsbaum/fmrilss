# Fast GLMsingle in fmrilss: due diligence and implementation plan

Status: planning. Reference implementation: GLMsingle Python, commit `1ab54a6`
(`cvnlab/GLMsingle`), with `fracridge` 3.0.

Goal: a `glmsingle_fast()` in fmrilss that reproduces GLMsingle's estimates
(types A–D) within floating-point tolerance and runs substantially faster.
Changes that alter the estimator (continuous HRFs, SURE, cross-relation
initialisation) are a separate, gated research track.

---

## 1. Bottom line

The proposal's core is sound. I verified each source-code claim it makes against
the pinned commit, and checked the two central exactness claims numerically
against GLMsingle's own `fracridge` and `calcbadness` code
(`.planning/glmsingle_checks/verify_core_identities.py`).

| Check | Result |
|---|---|
| Run-blocked fractional ridge vs `fracridge` on the stacked GLMsingle design | coef rel. err 8e-13, alpha rel. err 6e-13 |
| Compiled CV loss (`Σ d_i (z_i − m_i)^2 + c`) vs `calcbadness`, LORO and grouped folds, 2 sessions | rel. err 3e-16 |
| `fracridge` grid-interpolated alpha vs exact root of the proposal's fraction equation | median 1.2%, max 4% apart; achieved fraction off by ≤ 0.006 |

Two corrections are needed before building:

1. **`fracridge` does not solve the fraction equation.** It interpolates
   `log(1+α)` on a fixed grid (`10^k`, step 0.2 decades, spanning
   `1e-2·s_min²` to `1e4·s_max²`). The proposal's "solve the monotone scalar
   equation" therefore gives a slightly different estimator. To stay compatible
   we must reproduce the grid, using the global `s_min` and `s_max` over all run
   blocks, and the interpolation. The run-blocked version with the reproduced
   grid matched to 1e-12. An exact root solve can be offered as a non-default
   option.
2. **GLMsingle computes in float32.** It casts data and designs to float32, so
   comparing against the Python reference can only reach about 1e-4 relative
   error. Argmax/argmin decisions (HRF index, number of PCs, fraction index)
   will flip on near-ties. We therefore need two comparison tiers (§4).

The largest speedups come from not touching time series repeatedly, more than
from factorisation. The proposal underplays this. Details are in §3.

---

## 2. Claim-by-claim assessment

| § | Claim | Verdict | Notes |
|---|---|---|---|
| 1 | Each run's design carries columns for all trials | **Verified** | `designSINGLE[run]` is `T_r × numtrials` (glmsingle.py:583-612) |
| 1 | Fitter called once per fraction | **Verified** | Type C/D loop calls `glm_estimatemodel` per frac (glmsingle.py:~1530). Each call redoes the SVD, the dense projection and the full prediction |
| 1 | Dense `T×T` projectors | **Verified** | `make_projection_matrix` → `combinedmatrix @ data` and `polymatrix @ data` on every call |
| 2 | Block-diagonal design; pooled fraction preserves estimator | **Verified, with correction** | Correct, but the fraction→α map must replicate the `fracridge` grid (§1). Also replicate `tol=1e-10` zeroing of singular values per block |
| 2 | Cost goes from `O((Rn)^3)` to `O(Rn^3)` | Correct but minor | At `Rn≈750` the global SVD costs well under 1 s. The real waste is the per-call `T²V` projections and the time-series prediction (§3) |
| 3 | Low-rank nuisance application `Y − Q(QᵀY)` | **Verified, high value** | Removes the `T²V` term |
| 3 | Sherman–Morrison PC-prefix updates | Correct (I re-derived it) but **not needed** | Run blocks are about 60×60, so a fresh Cholesky per PC count costs microseconds. The win is reusing `QᵀY` prefixes. Re-solving is simpler and avoids the cancellation risk the proposal flags |
| 3 | Fit uses poly+extra projection; R² uses poly-only data and the raw design | **Verified** | glm_estimatemodel.py:528-562, 849-874. `R2 = 100·(1 − Σ‖P_poly(y − X_raw β)‖² / Σ‖P_poly y‖²)`, with no mean subtraction. It is computable from `n×n` and `n×V` statistics |
| 4 | Compile pairwise repeat CV into fixed targets and weights | **Verified exactly** | Normalisation is fixed from `results[0]` per session (calcbadness.py). Holds for the GLMdenoise stage too, where the reference is the 0-PC fit |
| 4 | Reconstruct only trials with `d_i > 0` during selection | Correct | Saves a little; the norm still uses all trials |
| 5 | Banded Gram plus Woodbury nuisance correction | Correct but **low priority** | The HRF is about 25–30 TRs long and the ISI is about 3 TRs, so the half-bandwidth is around 10 of about 60 trials. Dense BLAS wins. Defer |
| 5 | Delay-invariant autocorrelation | Correct but **not applicable** | The GLMsingle library is not a pure delay family |
| 6 | HRF basis compression with a projector-distance certificate | Correct bound, **probably useless in practice** | R² differences between library HRFs are often below `δ_h‖y‖²`. With the fast path, the exact 20-HRF sweep is cheap anyway (§3). Defer |
| 7 | VarPro profile objective and envelope gradient | Correct (I re-derived the gradient) | It changes the estimator, so it belongs in the research track. The HRF-scale normalisation point is important |
| 8 | Polyphase cross-relation shape initialiser | Mathematically correct (SIMO cross-relation) | Applies to NSD-like regular grids only. Polynomial projection breaks commutativity. Research only |
| 9 | SURE in amplitude space; evidence criterion | The formula is correct for fixed λ, and the fixed-fraction Jacobian caveat is right | It changes the selection target, so it is research only |
| 10 | Reject rank-one SVD of an unrestricted basis fit | Agree | |
| 11 | 43.5× kernel timing | Not reproducible here (no script) | Treat it as a hypothesis. It is a kernel benchmark, not end-to-end |

### Things the proposal missed

* **Matched-filter `XᵀY`.** Each un-projected trial column has only about
  `L ≈ 25` nonzeros. So `X_hᵀY` is a sparse cross-correlation costing
  `O(n·L·V)` per HRF per run, compared with `O(n·T·V)` dense. Then
  `XᵀMY = XᵀY − (XᵀQ)(QᵀY)`, which is exact. This makes the 20-HRF library
  sweep nearly free relative to reading the data.
* **No time-series predictions at all.** Every R², R2run, and
  fit/residual quantity GLMsingle reports is a quadratic form in
  `{β, XᵀP y, XᵀP X, ‖P y‖²}`. `glm_predictresponses` is never needed.
* **The GLMdenoise CV only needs the PC-selection voxels.** `glmbadness` is
  computed for every voxel, but `xvaltrend` uses only `ix`: voxels with
  `onoffR2 > pcR2cutoff` that also pass the mask, or a top-100 fallback.
  Restricting to `ix` is exact for every estimate. The only cost is the
  full-volume `glmbadness` diagnostic, which we can compute on request.
* **Betas only for the winning HRF.** The type-B stage stores betas for all
  `nh` HRFs and then discards all but one. The fast path scores all HRFs, then
  solves once.
* **All fractions share one grid pass.** In `fracridge`, the per-voxel cost for
  F fractions is one `newlen` evaluation plus one interpolation. So the whole
  path costs about the same as one fraction.

Rough cost model per run: GLMsingle does about 1 + 20 + 11·(HRF groups) + 2·20
`glm_estimatemodel` passes, each touching `T²V` (projection) plus `nTV`
(prediction). The fast path touches the data once for `QᵀY` and `‖P y‖²`, plus
20 sparse `n·L·V` correlations. Everything after that is `n`-space. For
`T≈220, n≈63, L≈25`, a ≥10× end-to-end speedup is plausible. It still needs
to be measured (task 4.3).

---

## 3. Where the work lives in fmrilss

New files, kept separate from the LSS code paths:

```
R/glmsingle.R              # glmsingle_fast() user API, options, outputs
R/glmsingle_ref.R          # slow, literal R port (oracle; internal)
R/glmsingle_design.R       # run geometry: onsets, trial↔run↔condition maps, HRF library
R/glmsingle_nuisance.R     # polynomial basis, extra regressors, PCs, orthonormal Q per run
R/glmsingle_cv.R           # compiled CV targets/weights, session normalisation
R/glmsingle_select.R       # findtailthreshold, select_noise_regressors, autoscale
src/glmsingle_kernels.cpp  # matched-filter XᵀY, run-block spectral ridge path, quadratic-form R²
inst/extdata/glmsingle_hrflibrary.tsv   # vendored (BSD-3, attribute Kendrick Kay)
tools/glmsingle_ref/       # pinned Python harness that writes golden fixtures
tests/testthat/test-glmsingle-*.R
```

Reusable pieces already in fmrilss: the `oasis` ridge/SE plumbing in
`R/oasis_backend.R` (pattern only), and the SBHM SVD-library code in
`R/sbhm_build.R` for the deferred §6 work. Prewhitening is out of scope for
parity because GLMsingle does not prewhiten.

---

## 4. Definition of "equally accurate"

* **Tier A, port fidelity:** `glmsingle_ref` (R, float64) vs Python GLMsingle
  (float32). Continuous outputs must agree to about 1e-4 relative. Discrete
  choices (HRFindex, pcnum, FRACindex) must agree everywhere except voxels
  whose winning margin is below a float32-scaled tie threshold, and those
  voxels are reported.
* **Tier B, optimisation exactness:** `glmsingle_fast` vs `glmsingle_ref`, both
  float64. Continuous outputs must agree to ≤ 1e-9 relative, and discrete
  choices must agree exactly except at exact ties. This is the contract that
  proves the speedups do not change the estimator.
* **Tier C, statistical accuracy:** only needed for research-track changes.
  Use a ground-truth simulation as in proposal §12.

---

## 5. Granular plan

Each task lists its deliverable and when it counts as done. Phases 0–4 are the
deliverable. Phase 5 is optional research.

### Phase 0: environment and reference harness

- **0.1 R toolchain.** The container has no R. Install R ≥ 4.3 plus
  Rcpp/RcppArmadillo/testthat/devtools, and add a SessionStart hook so cloud
  sessions can build and test. *Done when* `devtools::test()` passes on main.
- **0.2 Pin the Python reference.** `tools/glmsingle_ref/requirements.txt`
  pins GLMsingle@1ab54a6, fracridge 3.0, numpy and scipy. Run with
  `jit=False` or with numba pinned, and record which in the fixture metadata.
  *Done when* a script runs GLMsingle end to end on a synthetic dataset.
- **0.3 Fixture generator.** `make_fixtures.py` writes inputs and every
  intermediate: onoffR2, meanvol, FitHRFR2(+run), HRFindex, pcregressors,
  glmbadness, xvaltrend, pcnum, rrbadness, FRACvalue, scaleoffset, and
  betasmd A–D. Use little-endian `.bin` files plus a JSON manifest, read by a
  small R helper (no new dependencies). Cover 3–4 small scenarios: default; extra
  regressors; 2 sessions + custom xvalscheme; unequal run lengths. Keep each
  under 2 MB. *Done when* the fixtures are committed and regenerate
  deterministically.
- **0.4 Ground-truth simulator** `sim_glmsingle_data()` (tests/helper).
  Following proposal §12: condition means, a within-condition behavioural
  covariate, idiosyncratic trial variation, voxel-specific library HRFs,
  structured noise from a noise pool, and drift. *Done when* it is
  seed-deterministic and documented.
- **0.5 Baseline timing.** Time Python GLMsingle on (a) the simulator at
  realistic scale (12 runs × ~220 TRs × 50k voxels) and (b) the GLMsingle
  example NSD subset, if it can be downloaded. Record per-stage wall time.
  *Done when* the numbers are in `bench/glmsingle_baseline.md`.

### Phase 1: literal R port (`glmsingle_ref`, the oracle)

This phase is deliberately slow and stays line-by-line traceable to Python.
It is internal and not exported.

- **1.1 Design bookkeeping:** designSINGLE, stimorder, validcolumns, stimix,
  condcounts/condinruns, and the repeat checks that switch off
  wantglmdenoise/wantfracridge.
- **1.2 Polynomial basis** (`constructpolynomialmatrix`, default
  `maxpolydeg`) and the projection matrices.
- **1.3 HRF library:** vendor the TSV, port `getcanonicalhrflibrary`
  (stimdur convolution, pchip resampling to TR, peak normalisation), and
  `getcanonicalhrf`. Verify pchip against scipy on fixtures; R's
  `signal::pchip` is not a dependency, so write a small pchip.
- **1.4 `glm_estimatemodel` 'assume' path:** fit, R², R2run, meanvol,
  percent-BOLD scaling, and `olsmatrix2` zero-column handling.
- **1.5 Port `fracridge`** with the exact grid and interpolation (§1).
- **1.6 Stages:** FIR diagnostic (only what downstream needs), type A
  ON-OFF, `findtailthreshold` (GMM: port carefully or match the threshold on
  fixtures), noise pool, PC derivation (SVD of `T×T` per run, std
  scaling), type B library sweep, `calcbadness`, `select_noise_regressors`
  (pcstop), type C/D with the frac=1 prepend rule, and autoscale.
- **1.7 Tier A tests** against the 0.3 fixtures. *Done when* all scenarios pass
  the Tier A tolerance, with the tie-voxel list reported and empty or tiny.

Out of v1 scope (error clearly if requested): `hrfmodel='optimize'`,
`wantlss`, bootstrap/xval resampling modes, figure output, hdf5 output.

### Phase 2: fast kernels (each tested against the oracle at Tier B)

- **2.1 Run geometry object:** per run, onset indices, trial→column map,
  condition IDs, session, `T_r`. It is immutable and built once.
- **2.2 Nuisance basis per run:** orthonormal `Q_poly`, then
  `Q_full = [Q_poly, Q_extra, Q_pc(1..k)]` via Gram–Schmidt with
  re-orthogonalisation, so PC prefixes are nested. Unit test:
  `I − QQᵀ` equals `make_projection_matrix([poly, extra, pcs[:k]])` to 1e-12,
  including rank-deficient extras.
- **2.3 Data pass (C++, the only O(T·V) step):** per run and voxel tile,
  compute `QᵀY` (all PC columns), `‖y‖²`, and the matched-filter `X_hᵀY` for
  all `nh` library HRFs via sparse correlation. Unit test against the dense
  products at 1e-12. Optionally use OpenMP over voxel tiles (reuse
  `src/Makevars`).
- **2.4 Design statistics per (run, HRF, nuisance prefix):** `XᵀX`, `XᵀQ`,
  hence `G = XᵀX − (XᵀQ)(XᵀQ)ᵀ`, and `XᵀMY = XᵀY − (XᵀQ)(QᵀY)`. These are
  independent of voxels and cheap. Also keep the poly-only versions for R².
- **2.5 OLS solve and R² from quadratic forms:** `β = G⁺ b` (matching the
  `olsmatrix2` zero-column rule), with SSE and SST from sufficient statistics,
  per run and pooled. Test against oracle R², R2run, and betas.
- **2.6 Run-blocked spectral ridge path:** eigendecompose each `G_r`, pool the
  spectra, build the global `fracridge` grid, compute every fraction's α per
  voxel in one pass, and reconstruct coefficients only for the requested trial
  subset. Test against the oracle `fracridge` (expected ~1e-12, as already
  shown in Python).
- **2.7 Compiled CV:** from `(xvalscheme, validcolumns, stimix, session)`
  build the sparse `W` once, giving `d`, `used` and the condition/fold sums
  needed for `m`. Do not materialise `W` for large designs; use
  condition×run totals. Normalisation comes from the reference fit. Test
  against the ported `calcbadness` on random inputs, including zero-SD voxels
  (`zerodiv → 0`), conditions absent from training, and multi-run test folds.
- **2.8 Autoscale and percent-BOLD** as vectorised per-voxel 2×2 OLS. Match the
  `h[0] < 0 → [1, 0]` rule.

### Phase 3: assemble `glmsingle_fast()`

- **3.1 API:** `glmsingle_fast(design, data, stimdur, tr, opt = list())`.
  `data` is a list of `T_r × V` matrices (masking happens upstream).
  Option names mirror GLMsingle (`wantlibrary`, `wantglmdenoise`,
  `wantfracridge`, `fracs`, `n_pcs`, `pcstop`, `xvalscheme`,
  `sessionindicator`, `extra_regressors`, `maxpolydeg`, `brainthresh`,
  `brainR2`, `pcR2cutoff(mask)`, `wantpercentbold`, `wantautoscale`,
  `chunklen`). Validate them with the existing `.validate_option_names`
  helpers. The return value is a list with `typea..typed` matching GLMsingle
  field names.
- **3.2 Stage A + noise pool + PCs:** reuse the 2.3 data pass. PC derivation
  needs the poly-projected noise-pool time series, a bounded `T×|pool|`
  object, so it stays dense.
- **3.3 Stage B:** score all HRFs from the stored statistics, take the argmax,
  then compute betas once for the winner.
- **3.4 GLMdenoise selection:** for `ix` voxels only, by default, compute
  `n_pcs+1` OLS fits from prefix statistics, then compiled CV, then
  `xvaltrend` and `pcnum`. `opt$want_full_glmbadness = TRUE` computes it for
  all voxels to give an exact diagnostic.
- **3.5 Stages C/D:** C is one OLS fit at `pcnum`. D is 2.6 for all fractions,
  then 2.7 on the used trials, then the per-voxel argmin, then full
  reconstruction for the chosen fraction, then autoscale.
- **3.6 Voxel tiling and memory:** process in `chunklen` tiles, and group
  voxels by HRF index within a tile, as GLMsingle does. Peak memory should be
  O(tile) except for outputs.
- **3.7 Tier B suite:** `test-glmsingle-fast-vs-ref.R` covers all scenarios
  and stage-by-stage intermediates, not just final betas.

### Phase 4: validation, benchmarking, docs

- **4.1 Tier A end-to-end:** run `glmsingle_fast` against the Python fixtures.
  Run the example NSD subset locally (not in CI) and report beta correlation,
  max relative difference, and discrete-choice agreement maps.
- **4.2 Simulation accuracy report:** same estimator, so the metrics should be
  identical to Python up to ties. Record total and within-condition beta
  error, behavioural-effect recovery, and HRF selection accuracy. This becomes
  the baseline for Phase 5.
- **4.3 End-to-end benchmark** (`bench/glmsingle_fast.R`): per-stage wall time
  vs 0.5, at 1 thread and at N threads, plus peak RSS. *Target:* ≥ 5× end to
  end at 1 thread. *Stretch:* ≥ 10×.
- **4.4 Docs:** roxygen, a `glmsingle.Rmd` vignette (usage, parity evidence,
  timing table), NEWS, pkgdown reference entry, and LICENSE attribution for
  the vendored HRF library and ported algorithms (BSD-3 notice in
  `inst/COPYRIGHTS`).
- **4.5 `R CMD check`** clean, with test runtime under 60 s for the fast suite.

### Phase 5: research track (gated on Phases 0–4)

Each item ships as an opt-in option and needs a Tier C report on the 0.4
simulator, with held-out runs, showing it is no worse on within-condition and
behavioural-effect error.

- **5.1 HRF selection uncertainty:** store the R² margin and profile
  curvature per voxel, and optionally shrink toward neighbourhood HRFs where
  the profile is flat (§7).
- **5.2 Continuous HRF refinement** via the VarPro objective and gradient on
  the library-spanning basis, with fixed-peak normalisation (§7).
- **5.3 Cross-relation initialiser** for regular-grid designs (§8), used only
  as a candidate generator.
- **5.4 SURE/evidence fraction selection** as a fast mode for designs without
  repeats (§9), including the fixed-fraction Jacobian term.
- **5.5 Compressed library screening with a certificate** (§6), only if
  profiling shows the library sweep matters after 2.3.

---

## 6. Risks

| Risk | Mitigation |
|---|---|
| Discrete-choice flips at Tier A from float32 | Tie-margin reporting; Tier B is the real exactness gate |
| `findtailthreshold` (GMM fit) is hard to match bit-for-bit | Port it, but allow passing `brainR2`/`pcR2cutoff`. Fixtures record the Python threshold so later stages can be tested in isolation |
| Rank-deficient nuisance (custom extras collinear with polys/PCs) | Pivoted QR with a tolerance matching `make_projection_matrix`'s `inv` behaviour; dedicated test |
| Speedup smaller than hoped because data I/O dominates | Benchmark per stage in 0.5 first. If I/O dominates, the gain is still the removed `~70×` passes over the data |
| Python GLMsingle drifts upstream | Pin the commit and record it in fixture metadata |

## 7. Open decisions for the maintainer

1. Should the deliverable live in fmrilss (this plan) or in a separate package?
   fmrilss fits: shared Rcpp toolchain and HRF tooling.
2. Should parity be with Python GLMsingle (this plan) or MATLAB? They differ in
   small ways.
3. Is the `ix`-restricted GLMdenoise CV acceptable as the default? It is exact
   for all estimates and drops only the full-volume `glmbadness` diagnostic.

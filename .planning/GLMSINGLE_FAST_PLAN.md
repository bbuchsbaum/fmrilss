# Fast GLMsingle in fmrilss: due diligence and implementation plan

**Status:** plan v3.1. It incorporates external review 1 (disposition log in
§8) and the maintainer's policy decisions (§9).
**Reference implementation:** GLMsingle Python, `cvnlab/GLMsingle` at
`1ab54a6` (2025-11-09, current HEAD), with `fracridge` 3.0.
**Evidence:** `.planning/glmsingle_checks/`
- `verify_core_identities.py` — our original checks.
- `review1/` — the reviewer's reproducibility bundle, unmodified and
  hash-verified (`python3 verify_manifest.py`). It is re-run here on
  2026-10-07; logs are in `review1/rerun_2026-10-07/`.

## 0. Goal and decisions

Build `glmsingle()` in fmrilss. It must reproduce GLMsingle's single-trial
estimates (types A–D) under an explicit numerical contract (§4) and run
substantially faster. Estimator changes are a separate, gated research track
(Phase 8).

| Decision | Resolution |
|---|---|
| Where it lives | **fmrilss.** It shares the Rcpp/Armadillo/OpenMP toolchain, the fmrihrf/fmridesign integration, and the test and bench infrastructure |
| Parity target | **Python GLMsingle at `1ab54a6`.** It can be scripted without a licence and includes Python's latest fix (`19e6617`). MATLAB-only later commits are audited in §2.4 |
| PC-count CV on the selection subset only | **Yes, by default.** Estimates are unchanged. The full-volume `glmbadness` diagnostic is available on request. Both a lean-output and a matched-output mode are benchmarked (task 7.3) |
| Upstream quirks (§2.3) | **Default = the most sensible behaviour, not bug-for-bug parity.** Where upstream is internally inconsistent or the two ports disagree (S3, S6), fmrilss uses the defensible choice by default. `policy = "upstream"` reproduces pinned Python exactly, for comparisons and for the parity tests (§3.4, §9) |

## 1. Bottom line

The run-blocked architecture is sound. Both independent reviews reproduce the
core identities: pooled-spectrum `fracridge` matches to about 1e-12, and the
compiled repeat CV matches to about 3e-16. **"On these fixtures" is the
operative qualification.** What changed in v3:

1. **A conditioning-aware numerical contract** replaces the blanket 1e-4 / 1e-9
   tolerances. The reviewer's stress tests (re-run here) show that
   float32-vs-float64 disagreement reaches **2.96%** at κ(X)≈686 (ISI 1.2 s),
   and that Python's float32 path itself becomes meaningless at κ(X)≈2.5e4.
   See §4.
2. **Explicit source-compatibility policies.** v2 contained several wrong
   assumptions: a QR rank cutoff, `G⁺` for `olsmatrix2`, a universal nuisance
   basis, and a plain 2×2 autoscale. All are corrected in §2.3 and §3.4.
3. **Execution is explicitly staged.** v2 implied PC projections could happen
   in one initial data pass, but the PCs depend on stage A. See §3.5.
4. **PC-prefix optimisation is back on the shortlist.** Applying cached solves
   to thousands of right-hand sides is the real cost. Re-run here, the
   alternatives were 4.0–8.6× faster than repeated cached solves.
5. **HRF support comes from the sampled kernels.** It is 24–66 TRs, not 25.
   The "nearly free" claim for the library sweep is withdrawn.
6. **Delivery is reordered into vertical slices** (one working path first),
   and every release gets a statistical-accuracy gate, not only the research
   track.

## 2. Due diligence

### 2.1 Original proposal: claim-by-claim (updated)

| Proposal § | Claim | Verdict |
|---|---|---|
| 1 | Per-run design carries all-trial columns; fitter called per fraction; dense `T×T` projectors | **Verified** in source |
| 2 | Block-diagonal design; pooled fraction preserves the estimator | **Verified**, provided the `fracridge` grid and interpolation are reproduced with *global* spectral bounds. Calibrating a fraction separately in each run is a different estimator |
| 2 | `O((Rn)^3) → O(Rn^3)` | Correct. The kernel benchmark (now reproducible, see §2.2) shows 43.5–52× **for that kernel only** |
| 3 | Low-rank nuisance application | **Verified, high value** |
| 3 | Sherman–Morrison PC-prefix updates | **Revised:** a profiling-shortlist candidate. Rank-one β updates and batched prefix maps beat cached solves 4–8.6× on 1k–8k voxels (rerun). Cancellation guard required |
| 3 | Fit vs reporting projections differ | **Verified.** Score-only SSE with custom nuisances matches explicit predictions to 6e-16 (review1 probe) |
| 4 | Compiled repeat CV | **Verified.** Keep the constant `C` when reproducing reported `glmbadness`/`rrbadness` |
| 5 | Banded Gram + Woodbury | Correct; **defer** (dense wins at n≈60) |
| 6 | Compressed library with certificate | Correct bound; **defer** |
| 7–9 | VarPro, cross-relation initialiser, SURE/evidence | Correct mathematics; they change the estimator → **Phase 8** |
| 10 | Reject rank-one SVD shortcut | Agree |

### 2.2 Performance findings (corrected)

* **Sparse raw products.** Compute
  `X_hᵀ M_Q Y = X_hᵀY − (X_hᵀQ)(QᵀY)`, using the compact support of the *raw*
  trial columns. The residualised design is dense. Support comes from the
  sampled kernel: for stimdur 3 s it is 24–26 samples at TR 2.0, 36–39 at TR
  1.333, 48–52 at TR 1.0, and 61–66 at TR 0.8 (review1). The cost is about
  `H·n·L` MACs per voxel per run (≈ 50k for H=20, n=63, L=40). That is cheap
  relative to the dense passes, but not free. Implement it as a **filter bank**:
  per onset, contract an `L × V_tile` data window with the `H × L` bank,
  batching windows. Benchmark it against dense GEMM.
* **R² from sufficient statistics.** With `P` = poly-only projection, keep
  `s = yᵀPy`, `c = XᵀPy` and `G_P = XᵀPX`. Then
  `SSE = s − 2βᵀc + βᵀG_Pβ` and `R² = 100(1 − SSE/s)`. Pooled R² uses pooled
  sums of squares. This is scoped to the reported R²/R2run. Residual time
  courses, if ever requested, need more.
* **Score every HRF, reconstruct the winner.** Scoring needs the forward solve
  (`z = L⁻¹b`). Only the winner needs the back-solve and stored β. The
  custom-nuisance score correction is required (review1 formula).
* **Selection-subset PC CV.** The subset is every voxel with
  `onoffR2 > pcR2cutoff` that passes the mask, or a top-100 fallback when that
  set is empty. It can be thousands of voxels. It is distinct from the noise
  pool.
* **Fractions are not free.** Sharing the decomposition and grid removes
  repeated factorisation. Per-fraction shrinkage, candidate reconstruction and
  CV remain. Reconstruct only CV-used trial rows during selection.
* **Noise covariance.** PC derivation accumulates `T×T` covariance over the
  noise pool, costing `O(T²·V_pool)`. Use tiled symmetric rank-k updates
  (`syrk`) and never materialise the full normalised pool.
* **The 43.5× kernel benchmark** (`review1/original/check_algebra.py`)
  reproduces: 52× in our rerun. Its scope is fixed-penalty Gram, solve and
  reconstruction. If that kernel is half the runtime, a 45× gain yields only
  ≈ 2× end to end. The **≥ 5× end-to-end target** is an acceptance target, not
  a prediction.

### 2.3 Source-compatibility findings

Each finding is verified against `1ab54a6`, gets a fixture (task 0.3), and has
a policy in §3.4.

| ID | Behaviour in pinned source | Consequence for us | Found by |
|---|---|---|---|
| S1 | `make_projection_matrix` → `olsmatrix(mode=0)`: drop exactly-zero columns, unit-normalise, then `np.linalg.pinv` of the normalised Gram (numpy default `rcond`) | A conventional QR cutoff gives a different nuisance subspace. For `[q0, q1, q1+1e-9·q2]`, QR keeps rank 3 while the reference removes rank 2 (projector difference ≈ 1.0) | review1 |
| S2 | `olsmatrix2`: drop exactly-zero columns, then `np.linalg.solve`. Raises `LinAlgError` on singular Gram | `G⁺` (v2 wording) is a behaviour change | review1 |
| S3 | User `extra_regressors` are appended **only when the PC count > 0**, in both PC-CV and C/D. **MATLAB differs:** `GLMestimatesingletrial.m` (lines ~1167 and ~1355) *never* adds user extras in PC-CV or C/D; only the PCs are added. Both ports use extras in FIR/A/B | (a) If `pcnum = 0`, final C/D fits omit the user's extras. (b) The 0-PC fit is also the **held-out reference for every PC count** in `calcbadness`, so with extras supplied the entire PC-selection curve is scored against an extras-free reference. (c) Whenever extras are supplied, the two ports give different C/D estimates. A universal `[poly, extras, PCs₁..k]` basis matches neither port. **fmrilss default: extras in every fit** (§9) | review1 (a); us (b, c) |
| S4 | Autoscale uses `olsmatrix` (normalised-Gram pinv) on float32 `[β_f, 1]` | A plain 2×2 solve fails for constant candidate β. The reference returns a defined `(scale, offset)` | review1 |
| S5 | Type C/D R² is stored from the selected-fraction fit **before** autoscale and percent-BOLD scaling | R² must not be recomputed from the final betas | review1 |
| S6 | `zerodiv(..., wantcaution=0)` aliases `tmp = y` and sets zero divisors to 1 **in place**. MATLAB passes by value, so it has no such mutation; this is likely a Python port artefact | In `calcbadness`, zero-SD voxels give `z = 0` for `results[0]` but `(x − μ)/1` for later candidates. "Zero all candidates" is wrong (CV scores `[0, 2]` vs `[0, 0]`) | review1 |
| S7 | Type A: one ON–OFF column per run, fitted **stacked**, so one coefficient is shared across runs | Needs pooled scalar statistics, not per-run fits | review1 |
| S8 | `glmsingle.py` imports `select_noise_regressors` from `utils/`. A **different** function of the same name in `ols/make_poly_matrix.py` loops over `range(1, n)` and can never return 0 | The port must use the `utils/` version. Add a test that `pcnum = 0` is reachable | us |
| S9 | `findtailthreshold`: unseeded sklearn GMM (`n_init = 3`, `reg_covar = 0`) | Python results vary between runs. Fixtures record the thresholds used, and tests inject them | us |
| S10 | Mixed precision: float32 design/data products, then some float64 promotion, then a cast back | A float64 port cannot match bitwise; see §4 | review1 |

### 2.4 MATLAB-only commits after the Python fix

| Commit | Change | Plan |
|---|---|---|
| `a07f152` | Polynomial construction avoids the `T×T` projector | Memory only; no action |
| `91e5b7e` | GMM `RegularizationValue = eps` | Adopt it in our deterministic GMM as the one intentional deviation (§9, item 3) |

## 3. Architecture

### 3.1 Principles

1. **Layers with one-way dependencies.** Stages are pure functions on plain
   lists/matrices, with no I/O, global RNG or hidden state. The orchestrator is
   thin, and every stage boundary can inject Python intermediates.
2. **Geometry, voxel statistics and selection are separate objects with
   declared lifetimes** (§3.6).
3. **fmrilss conventions:** time × voxel `Y` with run identity from a
   `sampling_frame` or run ids, trials × voxel betas, an `*_options()`
   constructor, and the existing validators. GLMsingle names are kept for
   options and output fields.
4. **Numerical policy is a named module** (§3.4), not a set of scattered
   choices. Every rank decision, singular case and degenerate divisor goes
   through it.
5. **Hot kernels are in C++.** Orchestration and selection logic stay in R.
   No new hard dependencies.

### 3.2 Module layout

```
R/glmsingle.R            glmsingle(), glmsingle_design(): public API, staged orchestration
R/glmsingle_options.R    glmsingle_options() (defaults from glmsingle/defaults.py)
R/glmsingle_policy.R     numerical policies (§3.4): rank rules, singular handling, zerodiv, autoscale operator
R/glmsingle_runs.R       run geometry: onset TRs, trial↔run↔condition↔session maps, folds
R/glmsingle_hrf.R        HRF library: vendored TSV, stimdur conv, pchip→TR, peak-normalise, exact support
R/glmsingle_nuisance.R   per-run bases (policy-aware), PC derivation via tiled covariance
R/glmsingle_stats.R      wrappers for data-pass kernels; sufficient-statistic R²
R/glmsingle_ridge.R      run-blocked OLS/fracridge path, pooled grid + interpolation
R/glmsingle_select.R     compiled CV, tail threshold, select_noise_regressors, autoscale
R/glmsingle_methods.R    "glmsingle_fit": print/summary/coef/as.matrix
src/glmsingle_kernels.cpp  filter-bank XᵀY, QᵀY, syrk covariance, block spectral path, prefix maps
src/blocked_products.h     shared residualise-and-multiply helper (also OASIS)
inst/extdata/glmsingle_hrflibrary.tsv   vendored, BSD-3 (inst/COPYRIGHTS)
tools/glmsingle_ref/       pinned Python env + fixture writer (Rbuildignored)
tests/testthat/helper-glmsingle.R        fixture reader, simulator, compact float64 spec references
```

Dependency order: `options → policy → runs → hrf → nuisance → stats → ridge →
select → glmsingle()`. Budget: about 1,300 lines of R and about 500 of C++,
excluding tests. Split any file approaching 400 lines.

### 3.3 Reuse map

| Existing | Use | Action |
|---|---|---|
| `.voxhrf_orthonormal_span()` (R/voxel_hrf.R) | Orthonormal bases | Promote to `.orthonormal_span(X, policy = c("qr", "glmsingle"))` in aaa_utils. `"qr"` keeps current behaviour for the existing callers; `"glmsingle"` implements S1 |
| `oasis_AtY_SY_blocked`, `oasisk_compute_RY_norm2` | Blocked residualise-and-multiply | Extract into `blocked_products.h`. OASIS migrates separately, gated by its tests |
| `.as_*` validators, `.validate_option_names`, `oasis_options()` pattern | Options | Reuse / follow |
| `.set_beta_dimnames`, `.default_trial_names` | Output naming | Reuse |
| `lss_design()` + `.validate_design_models()` | `glmsingle_design()` fmridesign front-end | Mirror it; error if onsets are off the TR grid |
| `fmrihrf::sampling_frame/blocklens/blockids`; `fmrihrf::evaluate` | Run structure; non-parity user HRF libraries | Reuse |
| Makevars OpenMP; `bench/` conventions | Threading; benchmark script | Reuse |
| `.item_safe_solve()` | — | Not used: its chol→svd→pinv fallbacks contradict S2/S4 |
| `generate_rapid_design()` | — | Not used: continuous onsets, global `set.seed` |
| fmriAR prewhitening | — | Out of scope. The nuisance layer leaves a hook |
| `lss()` dispatch | — | Separate entry point (conditions, runs, sessions, four model types) |

### 3.4 Numerical policy (`glmsingle_policy.R`)

There are two policy sets, selected with `glmsingle_options(policy = )`:

* **`"fmrilss"` (default):** pinned behaviour everywhere, except where pinned
  behaviour is inconsistent or the ports disagree (S3, S6, GMM floor). There
  the defensible choice is used.
* **`"upstream"`:** reproduces pinned Python `1ab54a6` exactly, including S3,
  S6 and the unregularised GMM. It is used for parity testing and for users
  who need to match existing GLMsingle outputs.

The policy set used is recorded in `fit$options`. Individual policies can be
overridden for diagnostics.

| Policy | `"fmrilss"` (default) | `"upstream"` |
|---|---|---|
| Nuisance rank (S1) | Drop exact-zero columns; normalise; SVD of the normalised design with cutoff equivalent to `pinv(rcond = 1e-15)` on the Gram, i.e. singular values of X below `sqrt(1e-15)·s_max` are dropped | Same |
| Extras in PC-CV and C/D (S3) | **Always included:** nuisance basis `[poly, extras, PCs₁..k]` for every k, including k = 0, the CV reference, and final C/D | Python: extras only when k > 0 |
| Trial OLS singularity (S2) | Drop exact-zero columns, then Cholesky; on failure, error with a diagnostic naming the colliding trials. `"pinv"` opt-in | Same |
| Tail-threshold GMM | Deterministic EM, 3 fixed restarts, eps variance floor (MATLAB `91e5b7e`) | Thresholds injected from fixtures, since Python's GMM is unseeded |
| Autoscale (S4) | Normalised-Gram pinv on `[β_f, 1]`, `h[0] < 0 → (1, 0)` | Same |
| Divisor (S6) | **Zero for all candidates:** zero-SD voxels contribute nothing to CV, as in MATLAB | Emulate Python's in-place `zerodiv` mutation |
| Precision | float64 throughout; per-run mean-centring of `Y` before products (exact under polynomial projection, and it reduces cancellation in `XᵀY − (XᵀQ)(QᵀY)`) | — |
| Gram construction | Form the residualised design `A_r = X_r − Q(QᵀX_r)` explicitly per run (voxel-independent, cheap), with `G = AᵀA`. This avoids the `XᵀX − (XᵀQ)(XᵀQ)ᵀ` cancellation. Use a per-block QR/SVD path when `κ(A_r)` exceeds a threshold | — |
| Conditioning diagnostics | Record `κ(A_r)` per run/HRF, decision margins, and the policy used, in `fit$diagnostics` | — |

### 3.5 Staged execution

The stages run in dependency order. Each is tiled over voxels.

1. **Compile static geometry.** Run/trial maps, exact sampled HRFs and their
   support, polynomial bases, extras, residualised designs `A_r` per HRF, and
   design caches.
2. **Stage A, plus B where convenient, per voxel tile.** Mean volume, pooled
   ON–OFF statistics (S7), filter-bank `X_hᵀY`, poly-only sufficient
   statistics, HRF scores, winners and requested diagnostics.
3. **Noise pool and PCs.** Tail thresholds, then the pool, then tiled `syrk`
   accumulation of the per-run `T×T` covariance over detrended, normalised
   pool voxels, then eigendecomposition, then std scaling.
4. **PC-count selection on the selection subset only.** `QᵀY` prefixes, then
   prefix fits (cached solves first; prefix maps if profiling warrants), then
   compiled CV, the median trend and `pcstop`. This yields a global `pcnum`.
5. **Final C/D per tile.** Recompute or reuse winning-HRF products within the
   memory budget. Then the fraction path, CV on used rows, argmin,
   reconstruction of winners, R² recorded pre-scale (S5), autoscale, and
   percent-BOLD.

### 3.6 Memory budget

| Object | Lifetime | Size driver |
|---|---|---|
| Input `Y` (R matrix) | Whole call | `T·V·8` B. A user-owned input, not counted against our budget |
| Design cache (`A_r`, Gram factors per HRF/run) | Whole call | `H·R·(T_r·n + n²)`; small |
| Voxel statistics (`XᵀY` all HRFs, `QᵀY`, norms) | Per tile only | `H·N·V_tile·8`. **All 20 HRFs × 750 trials × 50k voxels would be 6 GB, so these are never held globally** |
| Noise covariance | Stage 3 | `R·T²·8`; small |
| Outputs (betas A–D, R², indices) | Whole call | `N·V·8` per requested type. Types not requested are not kept |

`chunklen` (the tile size) defaults to a value derived from a
`memory_limit_gb` option. Peak memory is reported in `fit$timing`.

## 4. Numerical contract

Three questions are kept separate: **source behaviour**, **mathematical
equivalence** and **statistical accuracy**.

**References.** These are mutually checking:
- **(R1)** Pinned Python with recorded intermediates. This is the source of
  truth for behaviour. Unpatched, it is the reference for `policy = "upstream"`.
  With the reviewed `fmrilss_policy.patch` applied, it is the reference for
  the default (task 0.3).
- **(R2)** Compact float64 spec references in `helper-glmsingle.R`: explicit
  stacked designs, literal `calcbadness` loops, per-fraction `fracridge`, dense
  projectors. These are test-only.
- **(R3)** Independent numerical checks: direct QR/SVD at a supplied penalty,
  residual orthogonality, and rank diagnostics.

R2 agreeing with our fast path does not prove R1 parity. Both R2 and R3 are
written from the source and the policy table, and are reviewed against R1
fixtures.

**Tolerances are conditioning-aware.** For a quantity computed from run block
`r`:

* **vs R2/R3 (float64):** `tol = c₆₄ · κ(A_r) · ε₆₄`, floored at 1e-12.
* **vs R1 (float32 reference):** `tol = c₃₂ · κ(A_r) · ε₃₂`. Blocks with
  `κ(A_r) · ε₃₂ > 0.1` are classed **reference-unstable**. They are reported,
  not asserted. The reviewer's κ≈2.5e4 float32 case (251% disagreement)
  falls in this class.
* `c₆₄` and `c₃₂` are fixed once from the fixture suite. They are recorded in
  the test helper and never tuned per test.

**Decisions are certified by margin.** For each argmax/argmin (HRF index,
fraction index, pcnum), compute the winning margin and an error bound from the
tolerance above. A decision is **certified** if the margin exceeds the bound.
Uncertified decisions must be reported, and are excluded from the equality
assertion but not from the report.

**Stage-conditional vs end-to-end parity.**
* *Stage-conditional:* each stage runs on R1's recorded inputs for that stage
  (thresholds, noise pool, PCs, pcnum, HRF index), so local errors cannot
  cascade.
* *End-to-end adaptive:* the complete automatic pipeline is compared to R1. Any
  divergence must trace to an uncertified or reference-unstable decision. The
  harness then re-runs downstream with R1's decision injected to confirm this.
  Global decisions are compared as wholes:
  - threshold values;
  - noise-pool set difference;
  - PC subspaces via principal angles, which is sign- and rotation-invariant
    when eigenvalues are near-degenerate;
  - `pcnum`.

**Statistical accuracy (Tier C) for every release.** On the ground-truth
simulator, report:
- trial-beta MSE and bias;
- preservation of within-condition variation, including recovery of the
  behavioural covariate;
- HRF selection accuracy;
- repeat reliability, which on its own can reward over-shrinkage.

For v1 the requirement is "identical to Python up to certified ties". For
Phase 8 the requirement is "no worse on held-out runs".

## 5. Plan

Each task is about one reviewable commit with a done-criterion. Phases 2–5 are
**vertical slices**: each ends with a working, tested path.

### Phase 0: freeze behaviour

- **0.1 R toolchain** in the container, plus a SessionStart hook. *Done when*
  `devtools::test()` passes on the branch.
- **0.2 Pinned Python env** (`tools/glmsingle_ref/requirements.txt`, matching
  review1's pins), `jit=False`. *Done when* GLMsingle runs end to end.
- **0.3 Fixture writer.** It records inputs, every stage intermediate,
  thresholds used, and κ per block. Scenarios:
  (a) defaults;
  (b) extras + `pcnum > 0`;
  (c) **extras + `pcnum = 0`** (S3);
  (d) 2 sessions + grouped `xvalscheme`;
  (e) unequal run lengths;
  (f) unrepeated conditions;
  (g) strongly overlapping events (ISI ≈ 1.2 s);
  (h) near-collinear extras (S1);
  (i) a constant-response voxel (S6) and a constant-candidate autoscale (S4);
  (j) a zero-variance/all-zero voxel.
  Fixtures are under 2 MB each. Each scenario is written twice:
  - by **pinned Python**, which is the reference for `policy = "upstream"`;
  - by **pinned Python plus `tools/glmsingle_ref/fmrilss_policy.patch`**, which
    is the reference for the default policy. The patch is minimal and reviewed:
    extras are included at every k, and `zerodiv` copies its divisor.

  Scenarios without extras or zero-SD voxels give identical outputs under both,
  and that identity is asserted. *Done when* they regenerate deterministically
  with thresholds pinned.
- **0.4 Ground-truth simulator** (test helper, local RNG). *Done when* it is
  documented and deterministic.
- **0.5 Baseline** per-stage timing and peak RSS of Python at 12 runs × 220
  TRs × 50k voxels, and on the GLMsingle example data if available.
  *Done when* the numbers are in `bench/results/glmsingle_baseline.md`.
- **0.6 Policy and contract freeze.** Review §3.4 and §4 against the fixtures
  and set `c₆₄` and `c₃₂`. *Done when* the policy table is signed off.

### Phase 1: groundwork

- **1.1** `.orthonormal_span(policy = )` promotion. Existing tests stay green;
  an S1 fixture test is added.
- **1.2** Extract `blocked_products.h`. The OASIS tests stay green.
- **1.3** `glmsingle_options()` and `glmsingle_policy.R`, with unit tests for
  each policy against the S1–S6 fixtures.

### Phase 2: slice 1, fixed HRF and fixed nuisance → type D

The HRF index and `pcnum` are injected from fixtures.

- **2.1** Run geometry (designSINGLE trial order, `validcolumns`, `stimix`,
  folds, repeat checks).
- **2.2** Nuisance bases for a given k under the S3 policy. The residualised
  designs `A_r`.
- **2.3** Run-blocked OLS and the `fracridge` path: pooled spectra, global
  grid, interpolation, `tol` zeroing, and an opt-in `exact_alpha`.
- **2.4** Compiled CV, including S6 emulation and the constant `C`.
- **2.5** Autoscale (S4), percent-BOLD, and R² recorded pre-scale (S5).
- *Gate:* stage-conditional parity vs R1 and R2 on all scenarios. Tier C
  metrics match Python.

### Phase 3: slice 2, HRF selection → type B

- **3.1** HRF library port with exact support (vs R1 at 1e-6 on float32
  values).
- **3.2** Filter-bank `X_hᵀY` kernel (vs dense at the R2 tolerance).
- **3.3** Sufficient-statistic R²/R2run and score-only selection (forward
  solve), with winner reconstruction.
- *Gate:* type B parity. HRF-index decisions are certified or reported.

### Phase 4: slice 3, adaptive denoising → types A and C, full automatic D

- **4.1** Type A with pooled scalar statistics (S7). Deterministic GMM tail
  threshold with an eps floor.
- **4.2** Noise pool and PCs via tiled `syrk` covariance. Compare subspaces by
  principal angles.
- **4.3** Selection-subset PC CV with cached solves, `utils`
  `select_noise_regressors` (S8), and `want_full_glmbadness`. *Test:* the
  restricted and full paths give identical `xvaltrend` and `pcnum`.
- *Gate:* end-to-end adaptive parity per §4. Every divergence is traced to an
  uncertified or reference-unstable decision.

### Phase 5: assembly

- **5.1** `glmsingle()` staged orchestrator (§3.5) with the memory budget (§3.6)
  and timing and diagnostics capture.
- **5.2** `glmsingle_design()` fmridesign front-end.
- **5.3** `glmsingle_fit` methods.

### Phase 6: optimise measured bottlenecks only

Profile first (task 7.3 tooling), then choose from:
- PC-prefix maps or rank-one updates (with a cancellation guard and fallback);
- filter-bank batching vs GEMM;
- OpenMP over tiles;
- CV-row-only reconstruction.

Each change must keep every Phase 2–4 gate green.

### Phase 7: validation, benchmarking, docs

- **7.1 Parity report:** all scenarios plus example data. Certified-decision
  rates, reference-unstable blocks, and divergence traces.
- **7.2 Tier C report** for v1 (the baseline for Phase 8). It includes
  `"fmrilss"` vs `"upstream"` on a simulator with motion-like extras that
  correlate with the noise. The default must be no worse than upstream on beta
  MSE, within-condition variation and behavioural-effect recovery. If it is
  worse, revisit the S3 default.
- **7.3 Matched benchmark** (`bench/run_glmsingle_benchmark.R`): same outputs
  and diagnostics, same thread count, same input precision, I/O excluded on
  both sides. Per-stage time and peak RSS at 1 and N threads, in lean and
  matched-output modes. *Target:* ≥ 5× end to end at 1 thread.
- **7.4** Optional OASIS migration to `blocked_products.h`.
- **7.5 Docs:** vignette (usage, a "Differences from GLMsingle" section listing every `"fmrilss"` vs `"upstream"` policy, parity evidence, timing), pkgdown,
  NEWS, `inst/COPYRIGHTS` (GLMsingle BSD-3, fracridge BSD-2).
- **7.6** `R CMD check` clean; the GLMsingle test suite runs in under 60 s.

Out of v1 scope (each errors clearly): `hrfmodel = "optimize"`, `wantlss`,
resampling modes, figures, hdf5.

### Phase 8: research track (opt-in; gated on a Tier C held-out report)

- 8.1 HRF uncertainty and shrinkage.
- 8.2 VarPro continuous HRF.
- 8.3 Cross-relation initialiser.
- 8.4 SURE/evidence fraction selection.
- 8.5 Compressed-library screening.
- 8.6 `exact_alpha` as a candidate path. It needs Tier C evidence, not just
  closer fraction attainment.

## 6. Risks

| Risk | Mitigation |
|---|---|
| Python float32 unstable on ill-conditioned designs | Reference-unstable class; report rather than assert; float64 spec references (R2/R3) are the algebraic gate |
| Global decisions cascade (threshold → pool → PCs → pcnum) | Stage-conditional parity plus traced end-to-end divergences |
| Shared misreading of the source in R2 and the fast path | R1 fixtures for every policy edge (0.3 c, h–j); policy table reviewed in 0.6 |
| Speedup below target because unoptimised stages dominate | Per-stage baseline (0.5) and profiling-driven Phase 6 |
| Memory blow-up from retained statistics | Per-tile lifetimes (§3.6); `memory_limit_gb` |
| Upstream changes | Pinned commit in the fixture manifest; re-audit on releases. Default-policy deviations are documented in the vignette |

## 7. Evidence index

| Claim | Evidence |
|---|---|
| Pooled `fracridge` ≡ stacked | `verify_core_identities.py`; `review1/.../run_reviewer_checks.py` |
| Interpolated vs exact α (≤ 4–5%) | Same, plus `review1/.../check_precision.py` |
| Precision stress (2.96% at κ≈686) | `review1/recorded_results/precision_results.json`; rerun log |
| S1, S2, S4, score-only SSE, HRF support | `review1/tmp/glmsingle/source_edge_probe.py` |
| S6 zerodiv mutation | `review1/rerun_2026-10-07/run_reviewer_checks.txt` |
| S3b, S8 | Source lines cited in §2.3 (`glmsingle.py` PC-CV loop; `ols/make_poly_matrix.py`) |
| PC-prefix 4–8.6× | `review1/tmp/glmsingle/due_diligence/pc_prefix_benchmark.py`; rerun log |
| 43.5–52× kernel | `review1/original/check_algebra.py`; rerun log |

## 8. Review 1: disposition log

| # | Reviewer point | Disposition |
|---|---|---|
| 1 | Core identities reproduce, but only on these fixtures | **Accepted.** Contract is conditioning-aware (§4) |
| 2 | Interpolation correction valid; calibrate α globally across runs | **Accepted** (already in v2); exact α stays opt-in research (8.6) |
| 3 | Blanket float32/float64 tolerances invalid | **Accepted.** §4 rewritten. Reproduced: 2.96% at κ≈686, and float32 meaningless at κ≈2.5e4 |
| 4 | A float64 oracle cannot prove Python parity; use three references; margin-based decisions; global decisions | **Accepted** (§4 R1–R3, certified decisions, stage-conditional vs end-to-end) |
| 5A | Nuisance projector rank policy | **Accepted** (S1, policy) |
| 5B | `olsmatrix2` uses solve, not `G⁺` | **Accepted** (S2, policy) |
| 5C | Extras omitted at zero PCs | **Accepted and extended:** the 0-PC fit is also the CV reference for every k (S3b). Maintainer decision needed (§9) |
| 5D | Autoscale rank handling; R² before autoscale; zerodiv mutation; type A pooled | **Accepted** (S4–S7) |
| 6 | Gram cancellation; need a stable fallback | **Accepted.** Explicit residualised design plus conditional QR/SVD (§3.4) |
| 7 | Sparse products: use actual support; not free; filter bank | **Accepted** (§2.2) |
| 8 | Score-only still needs the forward solve | **Accepted** |
| 9 | Selection subset may be thousands; benchmark lean vs matched output | **Accepted** (7.3) |
| 10 | Do not dismiss PC-prefix optimisation | **Accepted.** Rerun gave 4.0–8.6×. Phase 6 shortlist; cached solves first |
| 11 | Fractions do not cost the same as one | **Accepted** (v2 wording withdrawn) |
| 12 | 43.5× scope; ≥ 5× is a target, not a prediction | **Accepted** |
| 13 | Data flow must be staged; noise covariance `O(T²V_pool)` via `syrk` | **Accepted** (§3.5) |
| 14 | Concrete memory budget | **Accepted** (§3.6) |
| 15 | Shorten the full R port; small oracle; vertical slices | **Accepted.** v2 had already replaced the full port with test-only spec references; v3 adopts the slice ordering |
| 16 | Done-criteria: conditioning and adaptive coverage, PC sign ambiguity, decision margins, matched benchmarks, Tier C for v1 | **Accepted** (0.3, §4, 4.2, 7.3, 7.2) |

**Note:** the reviewer audited v1 (`40cb4cf`). Points 15 and part of 9 were
already addressed in v2 (`1dcaabf`).

## 9. Decisions

| # | Item | Decision |
|---|---|---|
| 1 | **S3 extras policy** | **Always include user extras**, in every PC-CV fit (including the k = 0 reference) and in final C/D. Rationale: (i) user nuisances such as motion are part of the noise model, and dropping them when no PCs are chosen leaves known confounds in the betas; (ii) it keeps the PC models nested (`[poly, extras] ⊂ [poly, extras, PC₁] ⊂ …`), so the CV curve compares like with like and the k = 0 reference is the same model family; (iii) it removes Python's discontinuity between k = 0 and k = 1; (iv) it is consistent with stages A and B, which already use the extras in both ports. MATLAB's "never in C/D" would let the GLMdenoise PCs silently replace user-specified nuisances, which is less defensible. `policy = "upstream"` reproduces Python. The choice is verified by Tier C (7.2) |
| 2 | **S2 singular trial Gram** | Error, as upstream does, with a diagnostic; `"pinv"` opt-in |
| 3 | **GMM eps floor** (MATLAB `91e5b7e`) | Adopt, with deterministic restarts |
| 4 | **S6 zero-SD voxels in CV** | Zero for all candidates (MATLAB semantics). Python's in-place divisor mutation is an aliasing artefact that scores constant voxels inconsistently across candidates |
| 5 | **Upstream issues** | None will be filed. Deviations are documented in the vignette's "Differences from GLMsingle" section |

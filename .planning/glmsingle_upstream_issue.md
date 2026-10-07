# Draft issue for cvnlab/GLMsingle

**Target:** https://github.com/cvnlab/GLMsingle/issues (to be filed by the
maintainer). All references are to commit `1ab54a65edd3ea41a6133d4b4ecb78a9c7296684`.

---

**Title:** Python and MATLAB handle `extra_regressors` differently in the GLMdenoise and fracridge stages (types C/D)

Hi, and thanks for GLMsingle. While checking the Python port line by line
against the MATLAB code, we found a difference in how user-supplied extra
regressors are used. We'd like to know which behaviour is intended.

**MATLAB** (`matlab/GLMestimatesingletrial.m`): in the PC cross-validation
loop (around line 1168) and the type C/D fits (around line 1356),
`optA.extraregressors` starts as empty cells, and **only the PCs** are
appended. User `extraregressors` are passed to the FIR, ON-OFF and FitHRF
stages (around lines 753, 849, 931 and 990), but not to GLMdenoise or
fracridge.

**Python** (`glmsingle/glmsingle.py`): in the PC cross-validation loop
(line 1346) and the type C/D fits (line 1541), user `extra_regressors` are
concatenated with the PCs, but **only inside `if n_pc > 0` /
`if pcnum > 0`**. As a result:

1. When `pcnum == 0` is selected, the final type C/D fits include neither the
   PCs nor the user's extra regressors. When `pcnum > 0`, they include both.
2. In the PC cross-validation, `results0[0]` (the 0-PC fit) is fitted without
   the extra regressors, but every `n_pc > 0` fit includes them.
   `calcbadness` uses `results[0]` as the held-out reference for every
   candidate (line 1369), so with extra regressors supplied, the whole
   PC-selection curve is scored against a reference that omits them.

So whenever `extra_regressors` is supplied, the two ports give different type
C/D estimates. Python's behaviour also changes discontinuously between 0 and
1 PCs.

**Question:** what is the intended behaviour for types C/D?
(a) Extra regressors are deliberately left out of GLMdenoise/fracridge, as in
MATLAB.
(b) They should be included at every PC count, including 0.
(c) Python's current behaviour is intended.

Two smaller Python/MATLAB differences we noticed along the way:

- `glmsingle/utils/zerodiv.py` (lines 48 and 56) sets `tmp = y` and then
  `tmp[bad] = 1`, which changes the caller's divisor in place. In
  `calcbadness`, the per-session `sd` is therefore changed after the first
  candidate. For voxels with zero SD, `results[0]` gets z = 0, but later
  candidates get `(x − mean)/1`. MATLAB passes by value, so it zeroes all of
  them.
- Two functions are named `select_noise_regressors`.
  `glmsingle/ols/make_poly_matrix.py:54` loops over `range(1, n)`, so it can
  never return 0 PCs. `glmsingle/utils/select_noise_regressors.py:4` (the
  one `glmsingle.py` imports) can. The unused copy may be worth removing to
  avoid confusion.

We can provide minimal reproducing scripts if that would help.

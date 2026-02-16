Vignette Review Checklist

Scope
- Files: `vignettes/getting_started.Rmd`, `vignettes/oasis_method.Rmd` (incl. OASIS‑VOXHRF section), `vignettes/oasis_theory.Rmd`, `vignettes/voxel-wise-hrf.Rmd`
- Build system: `knitr`/`rmarkdown` (ensure `VignetteBuilder: knitr` in `DESCRIPTION`)

Success Criteria
- Quality: clear goals, runnable end-to-end, readable figures, no dead ends.
- Style: consistent headers, chunk options, code formatting, terminology.
- Accuracy: outputs match code; claims align with package functions.
- Hygiene: clean builds, no spelling errors, no broken links/citations.
- Reproducibility: clean-session renders; seeds set; deterministic outputs; reasonable runtime.

Automated Checks
- Build vignettes
  - R: `devtools::check(build_vignettes = TRUE)` or per-file `rmarkdown::render()`.
- Spelling
  - R: `spelling::spell_check_files(list.files("vignettes", "*.Rmd", full.names = TRUE))`
  - Custom words in `inst/WORDLIST`.
- Lint/format code in Rmd
  - R: `styler::style_dir("vignettes", filetype = c("Rmd"))`
  - R: `lintr::lint_dir("vignettes")`
- Links/URLs
  - Render to HTML then check links; package: `urlchecker::url_check()`.

Quality Review
- Structure: intro (purpose, prerequisites, data, runtime); clear sections; closing with summary and `sessionInfo()`.
- Figures/tables: set fig width/height/dpi; captions; labeled axes/units; consistent theme.
- Copy: active voice; consistent tense; unify terminology (HRF, OASIS, LSS); remove TODOs.

Style Standards
- Headers: Title Case for H1/H2; sentence case elsewhere.
- Chunk defaults: `knitr::opts_chunk$set(cache = TRUE, dpi = 150, fig.width = 6, fig.height = 4, warning = FALSE, message = FALSE)`.
- Code: tidyverse style via `styler`; consistent pipe usage; meaningful names.
- References: consistent citation style; inline code/backticks; use `here::here()` if paths.

Accuracy & Reproducibility
- Fresh session renders: `callr::r(function() rmarkdown::render("vignettes/<file>.Rmd"))`.
- Seeds set at top; verify stable outputs across two renders.
- Verify numbers/plots against toy data or tests; ensure alignment with `R/lss.R`, `R/oasis_glue.R` (VOXHRF branch), and `R/oasis_hrf_recovery.R`.
- Document data assumptions (sampling rate, TR, units) and match code.

Performance & Practicality
- Runtime target: < 2–3 minutes per vignette; cache heavy chunks.
- Memory: limit large objects; remove intermediates; downsample examples if needed.
- Failure modes: graceful messages when optional deps are missing.

CI & Tooling
- CI: build vignettes in GitHub Actions (r-lib/actions, `build-vignettes: true`).
- Pre-commit hooks: styler, lintr, spelling, EoF fixer, trailing whitespace.
- PR checklist: clean build, spellcheck clear, links ok, runtime within budget.

## Test environments

- Local: macOS Sonoma 14.3, R 4.5.1, aarch64-apple-darwin20
- GitHub Actions: R-devel, R-release, and R-oldrel on Linux; R-release on
  macOS and Windows

## R CMD check results

Local `R CMD check --as-cran` results:

- 0 errors
- 1 warning
- 1 note

The warning is emitted from R's installed `R_ext/Boolean.h` when Homebrew clang
20 encounters an unknown diagnostic group (`-Wfixed-enum-extension`). No
package-owned source warning is reported.

The note reports that the locally installed HTML Tidy is too old to validate
the generated HTML manual. The PDF manual and all package documentation pass.

## Optional dependency

`fmridesign` is used only in optional design-aware workflows. It is listed in
`Suggests` and is available from the repository declared in
`Additional_repositories`. Calls, tests, and vignette evaluation are guarded
with `requireNamespace()`. The installed-package test suite and the
`lss_with_fmridesign` vignette have also been run with `fmridesign` absent.

## Submission

This is the first CRAN submission of fmrilss.

# Contributing

## Workflow

- Open an issue or draft pull request before larger changes.
- Keep changes narrowly scoped and explain user-facing impact in the PR.
- Run package checks for the paths you touched before pushing.

## Development Checks

- Regenerate documentation after roxygen changes.
- Run targeted tests when changing R or C++ code.
- Run `R CMD check` or the equivalent CI workflow before requesting
  review.

## Style

- Prefer small, well-named functions over adding complexity to existing
  ones.
- Keep examples minimal and deterministic.
- Do not commit generated build artifacts or local check directories.

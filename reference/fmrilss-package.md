# fmrilss: Least Squares Separate (LSS) Analysis for fMRI Data

This package implements efficient least squares separate (LSS) analysis
for functional magnetic resonance imaging (fMRI) data. LSS is used to
estimate trial-by-trial activation patterns in event-related fMRI
designs.

## Main functions

- [`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md): Main
  function for performing LSS analysis

- [`lss_naive`](https://bbuchsbaum.github.io/fmrilss/reference/lss_naive.md):
  Naive LSS implementation for reference

- [`project_confounds`](https://bbuchsbaum.github.io/fmrilss/reference/project_confounds.md):
  R implementation for projecting out confounds

- [`project_confounds_cpp`](https://bbuchsbaum.github.io/fmrilss/reference/project_confounds_cpp.md):
  Fast C++ confound projection

- [`lss_beta_cpp`](https://bbuchsbaum.github.io/fmrilss/reference/lss_beta_cpp.md):
  Vectorized C++ LSS beta computation

- [`get_data_matrix`](https://bbuchsbaum.github.io/fmrilss/reference/get_data_matrix.md):
  Helper function for data extraction

## Features

- Optimized C++ implementation using vectorized matrix algebra

- Memory-efficient projection without forming Q matrices

- Cholesky decomposition for numerical stability

- Fallback R implementation with QR decomposition

- Support for various design matrix configurations

- Robust numerical handling for edge cases

- OpenMP support for multi-core processing

## See also

Useful links:

- <https://bbuchsbaum.github.io/fmrilss/>

- <https://github.com/bbuchsbaum/fmrilss>

- Report bugs at <https://github.com/bbuchsbaum/fmrilss/issues>

## Author

**Maintainer**: Brad Buchsbaum <brad.buchsbaum@gmail.com>

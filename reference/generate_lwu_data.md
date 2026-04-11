# Generate Synthetic fMRI Data with LWU HRF

Creates synthetic fMRI time series using specified LWU HRF parameters

## Usage

``` r
generate_lwu_data(
  onsets,
  tau = 6,
  sigma = 2.5,
  rho = 0.35,
  TR = 1,
  total_time = 300,
  n_voxels = 10,
  amplitudes = NULL,
  noise_sd = 0.2,
  seed = NULL
)
```

## Arguments

- onsets:

  Vector of event onset times in seconds

- tau:

  LWU tau parameter (time-to-peak)

- sigma:

  LWU sigma parameter (width)

- rho:

  LWU rho parameter (undershoot amplitude)

- TR:

  Repetition time in seconds

- total_time:

  Total scan time in seconds

- n_voxels:

  Number of voxels to simulate

- amplitudes:

  Event amplitudes (scalar or vector)

- noise_sd:

  Standard deviation of noise

- seed:

  Random seed

## Value

List with Y (data matrix), true_hrf, true_betas, and design info

## Examples

``` r
if (FALSE) { # \dontrun{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
dim(sim$Y)
} # }
```

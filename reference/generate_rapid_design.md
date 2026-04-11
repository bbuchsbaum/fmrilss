# OASIS HRF Recovery Testing Functions

Functions to test OASIS's ability to recover HRF parameters from rapid
event-related designs with overlapping HRFs. Generate Rapid
Event-Related Design

## Usage

``` r
generate_rapid_design(
  n_events = 25,
  total_time = 300,
  min_isi = 2,
  max_isi = 4,
  seed = NULL
)
```

## Arguments

- n_events:

  Number of events to generate

- total_time:

  Total time in seconds

- min_isi:

  Minimum inter-stimulus interval in seconds

- max_isi:

  Maximum inter-stimulus interval in seconds

- seed:

  Random seed for reproducibility

## Value

Numeric vector of event onset times

## Details

Creates a rapid event-related design with specified ISI range

## Examples

``` r
generate_rapid_design(n_events = 5, total_time = 40, seed = 1)
#> [1]  5.000000  7.531017 10.275265 13.420972 17.237387
```

# Convert old-style whitening options to new format

Internal helper to maintain backward compatibility with the old
oasis\$whiten = "ar1" syntax by converting to new prewhiten format.

## Usage

``` r
.convert_legacy_whiten(oasis_opts)
```

## Arguments

- oasis_opts:

  List of OASIS options potentially containing whiten field

## Value

List of prewhiten options or NULL if no whitening

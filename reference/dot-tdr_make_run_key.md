# Generate a unique run key for cache directory isolation

Creates a key from timestamp + a random hex string to avoid collisions
across concurrent runs and across datasets.

## Usage

``` r
.tdr_make_run_key()
```

## Value

Character scalar, e.g. "run_20260302_143012_a7f3c1b2"

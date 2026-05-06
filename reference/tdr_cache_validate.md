# Validate that all on-disk cache files are intact

Walks the `@cellmap` sub-slots for path strings and checks file
existence. This function is called automatically (in quiet mode) when
entering
[`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md),
[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md),
and
[`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)
so that broken caches are caught early.

## Usage

``` r
tdr_cache_validate(.tdr.obj, .verbose = TRUE)
```

## Arguments

- .tdr.obj:

  A tinydenseR object.

- .verbose:

  Logical; if `TRUE`, print a summary. Default `TRUE`.

## Value

The (unmodified) `.tdr.obj`, invisibly. Stops with an error if any
cached file is missing.

## Details

Because the cache is ephemeral (stored under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) and cleaned up at
session end), heavyweight checksum verification is no longer performed.
Only file existence and schema version are checked.

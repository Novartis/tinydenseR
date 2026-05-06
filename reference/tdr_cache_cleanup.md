# Remove all cached files for a tinydenseR object

Deletes the on-disk cache directory, deregisters it from the session-end
cleanup finalizer, and clears the path strings from `@cellmap`. After
cleanup, the cellmap entries will be `NULL` and must be regenerated via
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

## Usage

``` r
tdr_cache_cleanup(.tdr.obj)
```

## Arguments

- .tdr.obj:

  A tinydenseR object.

## Value

Updated `.tdr.obj` with cache removed.

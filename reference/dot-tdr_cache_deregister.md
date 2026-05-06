# Remove a single cache directory from the cleanup registry

Called by
[`tdr_cache_cleanup()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_cleanup.md)
after manually deleting a cache so that the finalizer does not attempt a
redundant [`unlink()`](https://rdrr.io/r/base/unlink.html).

## Usage

``` r
.tdr_cache_deregister(cache_dir)
```

## Arguments

- cache_dir:

  Character path to deregister.

## Value

Invisible `NULL`.

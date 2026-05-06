# Register a cache directory for automatic removal at session end

Adds `cache_dir` to a package-private registry and, on first call,
attaches a `reg.finalizer(..., onexit = TRUE)` to a package-private
environment so that all registered directories are deleted when the R
session terminates.

## Usage

``` r
.tdr_cache_register_cleanup(cache_dir)
```

## Arguments

- cache_dir:

  Character path of the cache directory to register.

## Value

Invisible `NULL`.

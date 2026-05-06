# Return the root temporary cache directory for this R session

All
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
cache runs are stored under `<tempdir()>/tinydenseR_cache/`. The
directory is created lazily on first call.

## Usage

``` r
.tdr_cache_root()
```

## Value

Character scalar — the session-level cache root path.

# Print a human-readable summary of the on-disk cache state

Reports whether on-disk caching is active, the cache directory, schema
version, number of cached slots and files, and total size on disk.
Useful for interactive inspection and debugging.

## Usage

``` r
tdr_cache_info(.tdr.obj)
```

## Arguments

- .tdr.obj:

  A tinydenseR object.

## Value

A list (invisible) with components `active`, `root`, `schema_v`, `slots`
(character vector of cached slot names), `n_samples`, `n_files`, and
`total_bytes`.

## See also

[`tdr_cache_validate`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_validate.md),
[`tdr_cache_cleanup`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_cleanup.md),
[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running get.map() with on-disk caching
tdr_cache_info(lm.cells)
# On-disk cache: ACTIVE
# Directory:     /tmp/RtmpXXXXXX/tinydenseR_cache/run_20260302_143012_a7f3c1b2
# Schema:        v3
# Slots:         clustering, celltyping, nearest.lm, fuzzy.graphs
# Samples:       6
# Files:         24  (4 slots x 6 samples)
# Total size:    142.3 MB
} # }
```

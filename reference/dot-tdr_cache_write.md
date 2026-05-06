# Write one cached slot to disk (crash-safe via atomic rename)

Writes `object` to an RDS file inside `cache_dir/slot_name/` using an
atomic rename pattern (write to a `.rds.tmp` file, then
[`file.rename()`](https://rdrr.io/r/base/files.html) to the final
`.rds`). When the filelock package is installed, a shared file lock is
held during the rename to prevent concurrent writes from clobbering each
other.

## Usage

``` r
.tdr_cache_write(object, cache_dir, slot_name, sample_name, compress = FALSE)
```

## Arguments

- object:

  R object to serialize.

- cache_dir:

  Root cache directory for the current run.

- slot_name:

  Character – one of `"clustering"`, `"celltyping"`, `"nearest.lm"`,
  `"fuzzy.graphs"`.

- sample_name:

  Character – sample identifier (used as file stem).

- compress:

  Passed to `saveRDS`. Default `FALSE` for speed.

## Value

An attributed path string: a character scalar with `attr(, "schema_v")`
and `attr(, "bytes")` set.

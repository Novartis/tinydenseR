# Read a cached slot from disk

Accepts either a new-style attributed path string (schema v3+) or a
legacy metadata list (schema v2) with a `$path` element.

## Usage

``` r
.tdr_cache_read(meta_record)
```

## Arguments

- meta_record:

  An attributed path string (character scalar with `attr(, "schema_v")`
  and `attr(, "bytes")`), OR a legacy metadata list with `$path` and
  `$schema_v`.

## Value

The deserialized R object.

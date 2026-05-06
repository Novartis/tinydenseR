# Quiet internal cache validation for DE entry points

Runs the same checks as
[`tdr_cache_validate()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_validate.md)
but silently; it only stops on errors and never emits messages. Designed
to be injected at the top of
[`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md),
[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md),
etc.

## Usage

``` r
.tdr_cache_validate_quiet(.tdr.obj)
```

## Arguments

- .tdr.obj:

  A tinydenseR object.

## Value

`NULL` invisibly.

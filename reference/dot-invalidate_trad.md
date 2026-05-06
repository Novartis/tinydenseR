# Invalidate stale traditional analysis fits

Sets `@results$lm[[model]]$trad[[.annot.type]]` to `NULL` for every
stored model and emits a warning.

## Usage

``` r
.invalidate_trad(.tdr.obj, .annot.type = "celltyping")
```

## Arguments

- .tdr.obj:

  A TDRObj.

- .annot.type:

  Character: `"celltyping"` (default) or `"clustering"`.

## Value

The modified `.tdr.obj`.

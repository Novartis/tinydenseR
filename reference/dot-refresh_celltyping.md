# Refresh all celltyping-dependent slots

Inspects the TDRObj and updates every populated celltyping-dependent
slot. This is called at the end of
[`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
so that late-bound celltyping (after
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md))
produces a fully consistent object.

## Usage

``` r
.refresh_celltyping(.tdr.obj, .verbose = FALSE)
```

## Arguments

- .tdr.obj:

  A TDRObj with `@landmark.annot$celltyping$ids` already set to the new
  values.

- .verbose:

  Logical; print progress messages.

## Value

The modified `.tdr.obj`.

## Details

Slots that have not yet been computed (e.g.
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
has not run) are simply skipped — the object remains valid at whatever
pipeline stage it is in.

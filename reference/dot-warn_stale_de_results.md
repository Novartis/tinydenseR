# Warn about potentially stale pseudobulk / marker DE results

Checks whether any
[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
results (design or marker mode) were computed with `.id.from` matching
`.annot.type` and emits a warning if so.

## Usage

``` r
.warn_stale_de_results(.tdr.obj, .annot.type = "celltyping")
```

## Arguments

- .tdr.obj:

  A TDRObj.

- .annot.type:

  Character: `"celltyping"` (default) or `"clustering"`.

## Value

The (unmodified) `.tdr.obj`, invisibly.

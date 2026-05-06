# Recompute composition matrices from cell-level annotation IDs

Aggregates cell-level assignments into samples x populations count and
percentage matrices.

## Usage

``` r
.recompute_composition(.tdr.obj, .annot.type = "celltyping")
```

## Arguments

- .tdr.obj:

  A TDRObj with freshly relabeled annotation IDs.

- .annot.type:

  Character: `"celltyping"` (default) or `"clustering"`.

## Value

The modified `.tdr.obj`.

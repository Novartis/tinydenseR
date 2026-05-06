# List stored celltyping solutions

Returns the names of all stored celltyping solutions (excluding the
active `$ids` slot).

## Usage

``` r
list_celltyping_solutions(x, ...)

# S3 method for class 'TDRObj'
list_celltyping_solutions(x, ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, or SingleCellExperiment object.

- ...:

  Additional arguments passed to methods.

## Value

A character vector of solution names.

## See also

[`set_active_celltyping`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md),
[`import_cell_annotations`](https://opensource.nibr.com/tinydenseR/reference/import_cell_annotations.md)

## Examples

``` r
if (FALSE) { # \dontrun{
list_celltyping_solutions(tdr.obj)
# [1] "cell_type_l1" "cell_type_l2" "azimuth_pred"
} # }
```

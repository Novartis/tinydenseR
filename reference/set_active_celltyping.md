# Set active celltyping solution

Switches the active celltyping to a previously stored solution and
refreshes all celltyping-dependent downstream slots.

## Usage

``` r
# S3 method for class 'Seurat'
set_active_celltyping(x, ...)

# S3 method for class 'SingleCellExperiment'
set_active_celltyping(x, ...)

set_active_celltyping(x, ...)

# S3 method for class 'TDRObj'
set_active_celltyping(x, .column.name, .verbose = TRUE, ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object.

- ...:

  Additional arguments passed to methods.

- .column.name:

  Character: name of the stored celltyping solution to activate.

- .verbose:

  Logical (default TRUE).

## Value

Updated object with the selected celltyping as active `$ids`.

# Store a TDRObj inside a container object

Store a TDRObj inside a container object

## Usage

``` r
SetTDR(x, tdr, ...)

# Default S3 method
SetTDR(x, tdr, ...)

# S3 method for class 'Seurat'
SetTDR(x, tdr, ...)

# S3 method for class 'SingleCellExperiment'
SetTDR(x, tdr, ...)

# S3 method for class 'SummarizedExperiment'
SetTDR(x, tdr, ...)
```

## Arguments

- x:

  The container (Seurat, SCE, or TDRObj).

- tdr:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  to store.

- ...:

  Additional arguments (currently unused).

## Value

The container with the TDRObj stored.

## Methods (by class)

- `SetTDR(default)`: Default: if x is a TDRObj, returns tdr directly

- `SetTDR(Seurat)`: Store TDRObj in Seurat Misc slot

- `SetTDR(SingleCellExperiment)`: Store TDRObj in SCE metadata

- `SetTDR(SummarizedExperiment)`: Store TDRObj in SummarizedExperiment
  metadata

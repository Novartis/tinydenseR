# Extract a TDRObj from a container object

S3 generic that retrieves a
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
from the object in which
[`RunTDR`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)
stored it.

## Usage

``` r
GetTDR(x, ...)

# Default S3 method
GetTDR(x, ...)

# S3 method for class 'Seurat'
GetTDR(x, ...)

# S3 method for class 'SingleCellExperiment'
GetTDR(x, ...)

# S3 method for class 'SummarizedExperiment'
GetTDR(x, ...)
```

## Arguments

- x:

  An object that may contain a TDRObj.

- ...:

  Additional arguments (currently unused).

## Value

A
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

## Methods (by class)

- `GetTDR(default)`: Default method – returns a TDRObj as-is, errors
  otherwise

- `GetTDR(Seurat)`: Extract TDRObj from a Seurat object's Misc slot

- `GetTDR(SingleCellExperiment)`: Extract TDRObj from a
  SingleCellExperiment's metadata

- `GetTDR(SummarizedExperiment)`: Extract TDRObj from a
  SummarizedExperiment's metadata

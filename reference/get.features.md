# Graph-based feature discovery for landmarks

Identifies the most characteristic features (genes or markers) for each
landmark by determining which features contribute most strongly to the
landmark's defining principal component. This creates landmark-specific
feature signatures useful for interpretation.

## Usage

``` r
get.features(x, ...)

# S3 method for class 'TDRObj'
get.features(x, ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object.

- ...:

  Additional arguments passed to methods.

## Value

Updated `.tdr.obj` with `$interact.plot$lm.features` containing:

- `res`: List of data frames (one per landmark) with top 10 feature
  importances

- `html`: HTML-formatted tables for interactive plotting

Feature importance is the signed loading from the landmark's top
principal component.

## Details

**Algorithm:**

1.  For each landmark, identify the PC with highest absolute value (its
    "defining" PC)

2.  Extract that PC's rotation vector (feature loadings)

3.  Multiply by sign of landmark's PC score to preserve directionality

4.  Rank features by signed loading; keep top 10

Positive loadings indicate features that increase along the PC
direction, negative indicate features that decrease. The signed values
reflect how the landmark is characterized relative to the population
average.

Results are stored for use in `plotPCA` and `plotUMAP`, which displays
feature signatures when hovering over landmarks when
`.hover.stats = "marker"`.

## See also

[`plotPCA`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md),
[`plotUMAP`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md),
[`get.landmarks`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After landmark selection and graph construction
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph() |>
  get.features()

# View feature signature for first landmark
head(lm.cells$interact.plot$lm.features$res[[1]])

# Use in interactive plotting
plotPCA(lm.cells, .hover.stats = "marker")
} # }
```

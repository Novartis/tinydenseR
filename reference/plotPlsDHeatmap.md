# Plot plsD Expression Heatmap

Creates an expression heatmap with landmarks as columns and top-loaded
features as rows. Landmarks are ordered by plsD scores (or custom
ordering), and features are ranked by absolute regression loading. For
RNA data, expression is properly normalized, log-transformed, and
centered.

## Usage

``` r
plotPlsDHeatmap(x, ...)

# S3 method for class 'TDRObj'
plotPlsDHeatmap(
  x,
  .coef.col,
  .plsD.dim = 1,
  .model.name = "default",
  .n.features = 50,
  .order.by = "dens.contrast",
  .add.annot = NULL,
  .order.decreasing = FALSE,
  .viridis.options.annot = c("cividis", "rocket", "inferno", "mako", "magma"),
  .annot.panel.width = 4,
  .annot.panel.height = 0.15,
  .panel.width = 4,
  .panel.height = 3,
  .feature.font.size = 7,
  .show.landmark.labels = FALSE,
  .label.substr.rm = "",
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object after
  [`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md).

- ...:

  Additional arguments passed to methods.

- .coef.col:

  Character: coefficient name matching a slot in `.tdr.obj$plsD`.

- .plsD.dim:

  Integer or integer vector: plsD component(s) to use for: (1) ranking
  features by absolute loading, and (2) ordering landmarks (if
  `.order.by = NULL`). If vector (e.g., `c(1, 2)`), the first dimension
  is used for feature ranking, and landmarks are sorted by all
  dimensions sequentially. Default 1.

- .model.name:

  Character: name of the fitted model (default "default").

- .n.features:

  Integer: maximum number of features to display. For RNA, features are
  ranked by absolute loading and top `.n.features` are shown (half from
  positive loadings, half from negative). For cytometry, all markers are
  shown. Default 50.

- .order.by:

  Controls landmark ordering. One of:

  `"dens.contrast"`

  :   (default) Order by the raw (uncentered) density contrast
      coefficient from `.coef.col`, placing genuinely unaffected
      landmarks (raw Y \\\approx\\ 0) in the middle of the heatmap
      rather than at an arbitrary position determined by the mean.

  `"plsD.dim"`

  :   Order by the selected `.plsD.dim` score(s).

  A numeric matrix with column names

  :   Landmarks are sorted by columns sequentially (first column
      primary, etc.) and each column is displayed as an annotation strip
      with the column name as legend title.

- .add.annot:

  Numeric matrix with column names: displayed as annotation strips.
  These are added as additional annotation strips after the density
  contrast and before plsD score strips.

- .order.decreasing:

  Logical: sort landmarks in decreasing order? Default FALSE.

- .viridis.options.annot:

  Character vector: viridis color options for annotation strips. Cycled
  if fewer options than strips. Default
  `c("cividis", "rocket", "inferno", "mako", "magma")`.

- .annot.panel.width:

  Numeric: width of annotation strips in inches. Default 4.

- .annot.panel.height:

  Numeric: height of each annotation strip in inches. Default 0.15.

- .panel.width:

  Numeric: width of expression heatmap panel in inches. Default 4.

- .panel.height:

  Numeric: height of expression heatmap panel in inches. Default 3.

- .feature.font.size:

  Numeric: font size for feature labels. Default 7.

- .show.landmark.labels:

  Logical: show landmark IDs on x-axis? Default FALSE.

- .label.substr.rm:

  Character substring to remove from density contrast label (default
  "").

## Value

A ggplot2 object (gtable composition with annotation strips).

## Details

**Feature selection and ordering**: Features are *selected* by absolute
regression loading for the first element of `.plsD.dim`. For RNA, the
top `.n.features/2` positive and top `.n.features/2` negative-loaded
features are selected. Selected features are then *ordered* by signed
loading (highest positive at top, most negative at bottom).

**Loadings**: Uses OLS regression loadings (Xc ~ score_k): each loading
is the slope of a per-gene regression on the component score, capturing
marker-component association. Positive loadings = genes upregulated in
high-score (high-Y) landmarks; negative loadings = genes upregulated in
low-score landmarks. In datasets with structural score balancing (see
`get.plsD` Details), large-magnitude loadings at either end of the
ranking may partially reflect the geometric mean-zero constraint rather
than genuine differential expression. Use the raw Y annotation strip and
the `plotPlsD` scatter to assess whether extreme-loading features
correspond to landmarks with genuine density contrast signal.

Caution: depending on the direction of the contrast and the component,
large-magnitude features at either the positive or negative end of the
loading ranking may reflect structural score balancing (geometric
mean-zero constraint) rather than genuine biology. Always
cross-reference with the raw Y annotation strip: a feature with a large
negative loading that is highly expressed in landmarks with near-zero
raw Y is likely characterizing a structural counterweight region.

**Expression normalization** (RNA only):

1.  Computes row sums (per-landmark library size) from full
    `raw.landmarks`

2.  Subsets to top features

3.  Applies size factor normalization

4.  Log2-transforms: `log2(x + 1)`

5.  Centers each feature (row) to mean 0

**Annotation strips**:

1.  Density contrast: raw (uncentered) log2 fold change from the linear
    model. Zero = no density change. Positive = enriched in the contrast
    direction; negative = depleted.

2.  plsD scores for all dimensions in `.plsD.dim`

3.  Custom annotations (`.order.by` matrix columns, `.add.annot`
    columns)

## See also

[`get.plsD`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)
for computing plsD,
[`plotPlsD`](https://opensource.nibr.com/tinydenseR/reference/plotPlsD.md)
for score visualization

## Examples

``` r
if (FALSE) { # \dontrun{
# After running plsD
lm.obj <- get.plsD(lm.obj, .coef.col = "Infection")

# Basic heatmap using plsD1
plotPlsDHeatmap(lm.obj, .coef.col = "Infection", .plsD.dim = 1)

# Order by plsD dimension scores
plotPlsDHeatmap(lm.obj, .coef.col = "Infection", .plsD.dim = 1,
                 .order.by = "plsD.dim")
} # }
```

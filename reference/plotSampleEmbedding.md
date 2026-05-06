# Plot Sample Embedding from get.embedding

Visualizes sample-level embeddings computed by
[`get.embedding`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md).
Supports three embedding types: supervised partial-effect PCA (pePC),
unsupervised PCA on landmark densities (pca), and diffusion map
trajectory (traj).

## Usage

``` r
plotSampleEmbedding(x, ...)

# S3 method for class 'TDRObj'
plotSampleEmbedding(
  x,
  .embedding = "pePC",
  .sup.embed.slot = NULL,
  .color.by = NULL,
  .x.by = NULL,
  .pc.x = 1,
  .pc.y = 2,
  .cat.feature.color = Color.Palette[1, 1:5],
  .point.size = 1,
  .panel.size = 2,
  .midpoint = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  and
  [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md).
  All embeddings are stored in `.tdr.obj$map$embedding`: unsupervised
  (`pca`, `traj`) and supervised (`pePC$<name>`).

- ...:

  Additional arguments passed to methods.

- .embedding:

  Character string specifying which embedding type to plot. One of:

  "pePC"

  :   (default) Supervised partial-effect PCA. Requires
      `.sup.embed.slot`.

  "pca"

  :   Unsupervised PCA on log-transformed landmark densities. Uses
      `.tdr.obj$map$embedding$pca`.

  "traj"

  :   Diffusion map trajectory embedding. Uses
      `.tdr.obj$map$embedding$traj`.

- .sup.embed.slot:

  Character string specifying which supervised embedding to plot. Only
  used when `.embedding = "pePC"`. Must match a name in
  `.tdr.obj$map$embedding$pePC` (contrast name for FWL method, term name
  for nested model). Ignored (with warning) when `.embedding` is "pca"
  or "traj".

- .color.by:

  Character string specifying which metadata column to use for coloring
  points. Defaults to `.sup.embed.slot` if it exists in metadata (for
  pePC), otherwise uses the first metadata column.

- .x.by:

  Character string specifying which metadata column to use for the
  x-axis in single-PC embeddings. Defaults to `NULL`, which uses
  `.sup.embed.slot` if it exists in metadata. Required when
  `.sup.embed.slot` is not a metadata column (e.g., FWL contrast names).

- .pc.x:

  Integer specifying which PC/DC to plot on x-axis (default 1). Only
  used when embedding has more than 1 dimension.

- .pc.y:

  Integer specifying which PC/DC to plot on y-axis (default 2). Only
  used when embedding has more than 1 dimension.

- .cat.feature.color:

  Character vector of colors for categorical features. Defaults to first
  5 colors from `Color.Palette` row 1.

- .point.size:

  Numeric point size (default 1).

- .panel.size:

  Numeric panel size in inches (default 2).

- .midpoint:

  Numeric midpoint for continuous color scale. Defaults to median of
  `.color.by` column.

## Value

A `ggplot` object showing the sample embedding colored by the metadata
variable specified by `.color.by`.

## Details

This function plots the sample coordinates from embeddings computed by
[`get.embedding`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md):

- pePC (supervised):

  Partial-effect PCA isolating variation attributable to a specific
  contrast or set of covariates. Requires running
  [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md)
  with `.contrast.of.interest` or `.red.model`. Stored in
  `.tdr.obj$map$embedding$pePC$<name>`.

- pca (unsupervised):

  Standard PCA on log-transformed landmark densities. Good for
  exploratory visualization and QC. Stored in
  `.tdr.obj$map$embedding$pca`.

- traj (unsupervised):

  Diffusion map trajectory capturing continuous variation. Useful for
  developmental or temporal trajectories. Stored in
  `.tdr.obj$map$embedding$traj`.

For pePC embeddings with rank \> 1 (e.g., from nested model comparison
with multiple covariates), you can plot different PC combinations using
`.pc.x` and `.pc.y`.

When pePC embedding has only 1 PC (e.g., from a single contrast), the
function creates either a scatter plot (for continuous x-axis variable)
or a boxplot (for categorical x-axis variable), with the x-axis
determined by `.x.by` or `.sup.embed.slot`.

## See also

[`get.embedding`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md)
for computing embeddings,
[`plotSamplePCA`](https://opensource.nibr.com/tinydenseR/reference/plotSamplePCA.md)
for legacy PCA visualization

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute unsupervised embeddings
lm.obj <- get.embedding(.tdr.obj = lm.obj)

# Plot unsupervised PCA
plotSampleEmbedding(lm.obj, .embedding = "pca", .color.by = "Condition")

# Plot diffusion map trajectory
plotSampleEmbedding(lm.obj, .embedding = "traj", .color.by = "Timepoint")

# Compute supervised embedding
lm.obj <- get.embedding(lm.obj, .contrast.of.interest = "TrtVsCtrl")

# Plot supervised pePC embedding
plotSampleEmbedding(lm.obj, .embedding = "pePC",
                    .sup.embed.slot = "TrtVsCtrl", .color.by = "Group")
} # }
```

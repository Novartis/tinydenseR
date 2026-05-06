# Differential Density Testing

Performs landmark-based differential density testing using limma on size
factor-normalized and log-transformed fuzzy densities. Tests which
landmarks (cell states) change in density between conditions. More
sensitive than traditional cluster-level testing because landmarks can
capture within-cluster heterogeneity. Uses PCA-weighted q-values that
leverage the correlation structure among landmarks to improve
statistical power.

## Usage

``` r
get.lm(x, ...)

# S3 method for class 'TDRObj'
get.lm(
  x,
  .design,
  .contrasts = NULL,
  .block = NULL,
  .model.name = "default",
  .force.recalc = FALSE,
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

- ...:

  Additional arguments passed to methods.

- .design:

  Design matrix specifying the experimental design. Rows correspond to
  samples (matching `.tdr.obj$cells`), columns to coefficients. Create
  with [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

- .contrasts:

  Optional contrast matrix for specific comparisons. Each column defines
  one contrast. Create with
  [`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html).
  If NULL, tests all coefficients in `.design`.

- .block:

  Optional character: column name in `.tdr.obj$metadata` for blocking
  factor (e.g., "Donor", "Batch"). Accounts for within-block correlation
  using `duplicateCorrelation`.

- .model.name:

  Character string naming this model fit (default "default"). Results
  are stored in `.tdr.obj$map$lm[[.model.name]]`. Use different names to
  store multiple model fits (e.g., full vs reduced models for nested
  model comparisons).

- .force.recalc:

  Logical: if TRUE, overwrite existing results in the specified slot
  (default FALSE). If FALSE and slot already exists, an error is thrown.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

The modified `.tdr.obj` with results stored in
`.tdr.obj$map$lm[[.model.name]]`:

- `fit`:

  limma MArrayLM object after eBayes with moderated statistics,
  containing:

- `fit$coefficients`:

  Log fold changes (landmarks x coefficients), nested in fit

- `fit$p.value`:

  P-values (landmarks x coefficients), nested in fit

- `fit$density.weighted.bh.fdr`:

  Density-weighted BH FDR adjustment (landmarks x coefficients)

- `fit$pca.weighted.q`:

  PCA-weighted q-values (landmarks x coefficients)

- `trad`:

  List with traditional cluster-level DA results

- `trad$clustering$fit`:

  limma fit for cluster percentages with `$adj.p` slot

- `trad$celltyping$fit`:

  limma fit for celltype percentages with `$adj.p` slot (if celltyping
  exists)

## Details

The landmark-based DA testing workflow:

1.  Computes log2(fuzzy density + 0.5) for each landmark across samples

2.  Fits linear model: `lmFit(y ~ design)`

3.  If blocking specified: estimates within-block correlation via
    `duplicateCorrelation`

4.  Applies contrasts if provided

5.  Performs empirical Bayes moderation: `eBayes()`

6.  Computes density-weighted FDR and PCA-weighted q-values

7.  Performs traditional cluster/celltype-level DA for comparison

**Advantages over cluster-level testing:**

- Detects subtle shifts within clusters

- No arbitrary clustering thresholds

- Continuous representation of cell state space

- Powered by limma's variance shrinkage

**Design matrix considerations:**

- Include intercept for standard comparisons

- No-intercept models with continuous covariates assume zero-centering
  (warning issued)

- Must be full rank (checked automatically)

**PCA-weighted q-value procedure:**

Standard FDR methods assume independent tests, but landmarks are
correlated in the cell state space. To account for this:

1.  **Residualize PCs**: Regress out design matrix from landmark PCA
    embeddings to obtain covariates independent of experimental factors

2.  **Estimate pi0**: Use
    [`swfdr::lm_pi0`](https://rdrr.io/pkg/swfdr/man/lm_pi0.html) to
    estimate the proportion of true nulls conditional on PC coordinates.
    Landmarks in similar cell states have correlated p-values; PCA
    captures this spatial structure

3.  **Conservative safeguard**: If global pi0 \< 0.6, apply pi0 floor of
    0.6 to prevent anti-conservative estimates in high-signal regimes

4.  **Compute q-values**: Use
    [`swfdr::lm_qvalue`](https://rdrr.io/pkg/swfdr/man/lm_qvalue.html)
    to calculate spatially-weighted FDR that properly accounts for
    landmark correlation structure

This approach is more powerful than standard BH adjustment when
landmarks exhibit correlated expression, while maintaining proper FDR
control. Falls back to standard q-value or BH for small test sets
(\<1000 landmarks).

## See also

[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
(required predecessor),
[`plotBeeswarm`](https://opensource.nibr.com/tinydenseR/reference/plotBeeswarm.md)
for visualization,
[`plotTradStats`](https://opensource.nibr.com/tinydenseR/reference/plotTradStats.md)
for comparison with cluster-level tests

## Examples

``` r
if (FALSE) { # \dontrun{
# After mapping
lm.obj <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |> get.graph() |> get.map()

# Simple two-group comparison (results stored in lm.obj$map$lm$default)
design <- model.matrix(~ Condition, data = .meta)
lm.obj <- get.lm(lm.obj, .design = design)

# Visualize results
plotBeeswarm(lm.obj, .coefs = "ConditionB")

# Complex design with contrasts
design <- model.matrix(~ 0 + Group, data = .meta)
contrasts <- limma::makeContrasts(
  TvsC = GroupTreatment - GroupControl,
  T1vsT2 = GroupTreatment1 - GroupTreatment2,
  levels = design
)
lm.obj <- get.lm(lm.obj, .design = design, .contrasts = contrasts,
                    .model.name = "full", .force.recalc = TRUE)

# With blocking for paired samples
design <- model.matrix(~ Timepoint, data = .meta)
lm.obj <- get.lm(lm.obj, .design = design, .block = "Subject",
                    .model.name = "timepoint")

# Nested model comparison: fit reduced model
red.design <- model.matrix(~ 1, data = .meta)  # intercept only
lm.obj <- get.lm(lm.obj, .design = red.design, .model.name = "reduced")

# Access results by model name
lm.obj$map$lm$full$fit$coefficients
lm.obj$map$lm$reduced$fit$coefficients
} # }
```

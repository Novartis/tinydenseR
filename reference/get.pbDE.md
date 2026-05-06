# Pseudobulk Differential Expression Analysis

Performs pseudobulk differential expression (DE) analysis for
genes/markers with two modes.

## Usage

``` r
get.pbDE(x, ...)

# S3 method for class 'TDRObj'
get.pbDE(
  x,
  .source = NULL,
  .mode = NULL,
  .design = NULL,
  .contrasts = NULL,
  .block = NULL,
  .geneset.ls = NULL,
  .id = NULL,
  .id.idx = NULL,
  .id2 = "..all.other.landmarks..",
  .id2.idx = NULL,
  .id.from = NULL,
  .model.name = "default",
  .result.name = NULL,
  .population.name = NULL,
  .comparison.name = NULL,
  .force.recalc = FALSE,
  .verbose = TRUE,
  .label.confidence = 0.5,
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

- .source:

  The raw data object for non-file backends. `NULL` (default) for the
  files backend; otherwise a Seurat, SingleCellExperiment, or anndataR
  AnnData object. Used by `.get_sample_matrix()` to retrieve per-sample
  expression matrices.

- .mode:

  Character: analysis mode. One of `"design"` or `"marker"`. If `NULL`
  (default), auto-detected from arguments:

  - `.design` provided \\\Rightarrow\\ `"design"`

  - `.id`/`.id.idx` provided without `.design` \\\Rightarrow\\
    `"marker"`

- .design:

  Design matrix specifying experimental design (design mode only). Rows
  = samples, columns = coefficients.

- .contrasts:

  Optional contrast matrix for specific comparisons (design mode only).
  Create with
  [`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html).
  If NULL, tests all `.design` coefficients.

- .block:

  Optional character: column name in `.tdr.obj$metadata` for blocking
  factor (design mode only, e.g., "Donor"). Accounts for within-block
  correlation.

- .geneset.ls:

  Optional named list of character vectors defining gene sets for GSVA
  enrichment analysis. Only for RNA data. Example:
  `list("Tcell" = c("CD3D", "CD3E"), "Bcell" = c("CD19", "MS4A1"))`.

- .id:

  Optional character vector of cluster/celltype IDs. In design mode,
  restricts analysis to cells matching these IDs. In marker mode,
  defines group 1 (test group).

- .id.idx:

  Optional integer vector specifying landmark indices. In design mode,
  restricts analysis to cells confidently assigned to these landmarks.
  In marker mode, defines group 1 landmark indices.

- .id2:

  Character vector of cluster/celltype IDs for group 2 (marker mode
  only). Default `"..all.other.landmarks.."` compares group 1 to all
  other cells. Can specify specific IDs for pairwise comparisons.

- .id2.idx:

  Optional integer vector specifying landmark indices for group 2
  (marker mode only). When provided, takes priority over `.id2`.

- .id.from:

  Character: `"clustering"` or `"celltyping"`. Source of IDs in `.id`
  and `.id2`. Default `NULL` (resolved to `"clustering"` when needed).

- .model.name:

  Character string naming this model fit (default `"default"`).

- .result.name:

  Character string naming this result. In design mode defaults to
  `"all"`; in marker mode auto-generated from `.id` and `.id2`. Used as
  the storage key: `$pbDE[[.model.name]][[.result.name]]` (design) or
  `$markerDE[[.model.name]][[.result.name]]` (marker).

- .population.name:

  `NULL`. Deprecated alias for `.result.name` (design mode).

- .comparison.name:

  `NULL`. Deprecated alias for `.result.name` (marker mode).

- .force.recalc:

  Logical: if TRUE, overwrite existing results in the specified slot
  (default FALSE).

- .verbose:

  Logical: print progress messages? Default TRUE.

- .label.confidence:

  Numeric scalar in `[0,1]` controlling the minimum posterior confidence
  required to assign a cell to a set of target landmarks (used when
  `.id.idx` or `.id2.idx` is provided). Default 0.5.

## Value

The modified `.tdr.obj` with results stored depending on mode:

### Design mode

Results in `.tdr.obj$pbDE[[.model.name]][[.result.name]]`:

- `coefficients`:

  Log fold change matrix (features x coefficients)

- `p.value`:

  P-values (features x coefficients)

- `adj.p`:

  FDR-adjusted p-values (features x coefficients)

- `smpl.outlier`:

  Logical vector indicating outlier samples

- `id.idx`:

  Per-sample list of cell indices used

- `n.pseudo`:

  Integer vector of pseudobulk cell counts per sample

- `geneset`:

  (RNA + `.geneset.ls`) GSVA results

### Marker mode

Results in `.tdr.obj$markerDE[[.model.name]][[.result.name]]`:

- `coefficients`:

  Log fold changes (features x coefficients)

- `p.value`:

  P-values (features x coefficients)

- `adj.p`:

  FDR-adjusted p-values (features x coefficients)

- `smpl.outlier.1`:

  Logical: samples excluded from group 1

- `smpl.outlier.2`:

  Logical: samples excluded from group 2

- `id1.idx`:

  Per-sample cell indices for group 1

- `id2.idx`:

  Per-sample cell indices for group 2

- `n.pseudo1`:

  Pseudobulk cell counts per sample for group 1

- `n.pseudo2`:

  Pseudobulk cell counts per sample for group 2

- `geneset`:

  (RNA + `.geneset.ls`) GSVA results

## Details

**Design mode** (`.mode = "design"`): Aggregates cells into pseudobulk
samples using fuzzy landmark membership, then uses limma-voom (RNA) or
limma (cytometry) to test for DE across experimental conditions. Uses a
user-supplied design matrix with samples as replicates.

**Marker mode** (`.mode = "marker"`): Identifies marker genes/proteins
distinguishing one cell population from another (or all others) via a
within-sample paired comparison. Unlike design mode which tests
experimental conditions, marker mode compares cell populations to find
defining features.

### Design mode

Tests for DE across experimental conditions using a user-supplied design
matrix:

1.  **Cell selection**: If `.id` specified, select matching cells. If
    `.id.idx` specified, use fuzzy confidence thresholding. Otherwise
    use all cells.

2.  **Pseudobulk aggregation**: Fuzzy-weighted expression per sample.

3.  **Outlier removal**: (RNA only) Exclude samples with \<10\\

4.  **Normalization**: TMM + voom (RNA) or as-is (cytometry).

5.  **Linear modeling**:
    [`limma::lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) with
    optional blocking.

6.  **Empirical Bayes**: `limma::eBayes(robust = TRUE)`.

7.  **GSVA** (optional, RNA only): Gene set variation analysis.

### Marker mode

Compares cell populations via a within-sample paired design
(`~ .ids + .pairs`):

1.  **Cell selection**: Extract cells for group 1 (`.id`) and group 2
    (`.id2`).

2.  **Pseudobulk aggregation**: Independent aggregation per group per
    sample.

3.  **Outlier removal**: (RNA only) Cross-group outlier flagging.

4.  **Paired comparison**:
    [`limma::lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) with
    sample-as-pair blocking.

Positive logFC = higher in group 1, negative = higher in group 2.

## See also

[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
(required),
[`plotPbDE`](https://opensource.nibr.com/tinydenseR/reference/plotPbDE.md),
[`plotMarkerDE`](https://opensource.nibr.com/tinydenseR/reference/plotMarkerDE.md),
[`get.plsD`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After mapping
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |> get.graph() |> get.map()

# --- Design mode (default when .design is provided) ---
design <- model.matrix(~ Condition, data = .meta)
lm.cells <- get.pbDE(lm.cells, .design = design)
plotPbDE(lm.cells, .coefs = "ConditionB")

# DE within specific cell type
lm.cells <- get.pbDE(lm.cells, .design = design,
                     .id = c("1", "2", "3"),
                     .id.from = "clustering",
                     .result.name = "tcells")

# --- Marker mode (auto-detected when .id is provided without .design) ---
lm.cells <- get.pbDE(lm.cells, .mode = "marker",
                     .id = "cluster.3",
                     .result.name = "cluster3_markers")
plotMarkerDE(lm.cells, .comparison.name = "cluster3_markers")

# Pairwise comparison
lm.cells <- get.pbDE(lm.cells, .mode = "marker",
                     .id = c("cluster.1", "cluster.3"),
                     .id2 = c("cluster.2", "cluster.4"),
                     .result.name = "cd4_vs_cd8")

# Access results
lm.cells$pbDE$default$all$adj.p
lm.cells$markerDE$default$cluster3_markers$coefficients
} # }
```

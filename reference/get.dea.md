# Pseudobulk Differential Expression Analysis (Deprecated)

**\[deprecated\]**

`get.dea()` has been renamed to
[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
for clarity. This function is provided for backward compatibility and
will be removed in a future version.

## Usage

``` r
get.dea(
  .tdr.obj,
  .design,
  .contrasts = NULL,
  .block = NULL,
  .geneset.ls = NULL,
  .id.idx = NULL,
  .id = NULL,
  .id.from = NULL,
  .verbose = TRUE,
  .label.confidence = 0.5
)
```

## Arguments

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

- .id.idx:

  Optional integer vector specifying landmark indices. In design mode,
  restricts analysis to cells confidently assigned to these landmarks.
  In marker mode, defines group 1 landmark indices.

- .id:

  Optional character vector of cluster/celltype IDs. In design mode,
  restricts analysis to cells matching these IDs. In marker mode,
  defines group 1 (test group).

- .id.from:

  Character: `"clustering"` or `"celltyping"`. Source of IDs in `.id`
  and `.id2`. Default `NULL` (resolved to `"clustering"` when needed).

- .verbose:

  Logical: print progress messages? Default TRUE.

- .label.confidence:

  Numeric scalar in `[0,1]` controlling the minimum posterior confidence
  required to assign a cell to a set of target landmarks (used when
  `.id.idx` or `.id2.idx` is provided). Default 0.5.

## Value

A list containing DE analysis results (legacy format). Use
[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
instead which stores results in `.tdr.obj$pbDE` and returns the modified
object.

## See also

[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)

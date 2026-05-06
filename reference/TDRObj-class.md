# TDRObj: S4 class for tinydenseR analysis objects

The core data structure for tinydenseR landmark-based analysis. Contains
expression data references, metadata, dimensionality reduction results,
graph structure, and differential expression results.

## Usage

``` r
# S4 method for class 'TDRObj'
x$name

# S4 method for class 'TDRObj'
x$name <- value

# S4 method for class 'TDRObj'
names(x)

# S4 method for class 'TDRObj'
show(object)
```

## Arguments

- x:

  A TDRObj object.

- name:

  A character string naming the slot to access.

- value:

  The value to assign to the slot.

- object:

  A TDRObj object (used in show method).

## Slots

- `cells`:

  list. Named list of per-sample file paths to expression matrices.

- `metadata`:

  data.frame. Sample-level metadata.

- `config`:

  list. Run parameters: key, sampling, assay.type, markers, n.threads.

- `integration`:

  list. Trained projection models and batch variables (harmony.var,
  harmony.obj, symphony.obj, umap.model).

- `assay`:

  list. Landmark expression layers (L x features matrices): raw (raw
  counts), expr (normalized/log expression), scaled (Z-scored).

- `landmark.embed`:

  list. Landmark-space coordinate matrices; each entry has a \$coord
  slot. Contains pca, le, and umap sub-lists.

- `landmark.annot`:

  list. Per-landmark categorical annotations (factor, length L).
  Contains clustering and celltyping sub-lists, each with an \$ids
  factor.

- `graphs`:

  list. Landmark-landmark connectivity matrices: adj.matrix, snn,
  fgraph.

- `density`:

  list. Fuzzy density analytics populated by
  [`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).
  Contains five sub-elements:

  raw

  :   L × N matrix of pre-normalization fuzzy density sums.

  norm

  :   L × N matrix after size-factor normalization:
      `norm = t(t(raw) / size.factors)`.

  log.norm

  :   L × N matrix: `log2(norm + 0.5)`.

  size.factors

  :   Named numeric(N): \\n_j / \bar{n}\\, guaranteed to have mean 1.

  composition

  :   List of clustering/celltyping cell count and percentage matrices
      (samples × clusters).

  Access via
  [`get.density`](https://opensource.nibr.com/tinydenseR/reference/get.density.md)
  or the `$` accessor.

- `sample.embed`:

  list. Sample-level embeddings (N x k matrices), each with \$coord.
  Contains pca, traj, and pepc sub-lists.

- `cellmap`:

  list. Per-cell, per-sample data in unified structure: clustering\$ids,
  celltyping\$ids (named per-sample lists with optional named
  solutions), nearest.lm, fuzzy.graphs. Each sample entry is either an
  in-memory R object or an attributed path string for on-disk cache.

- `results`:

  list. All statistical outputs: lm, pb, marker, spec, nmf, pls,
  features.

## See also

[`as.SummarizedExperiment.TDRObj`](https://opensource.nibr.com/tinydenseR/reference/as.SummarizedExperiment.TDRObj.md)
for converting to SummarizedExperiment;
[`GetTDR`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md),
[`SetTDR`](https://opensource.nibr.com/tinydenseR/reference/SetTDR.md)
for container extraction/storage.

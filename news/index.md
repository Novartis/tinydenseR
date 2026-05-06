# Changelog

## tinydenseR 0.0.3.0001

### Framework overview

tinydenseR is a clustering-independent framework for sample-level
modeling of single-cell data. It constructs a landmark-by-sample fuzzy
density matrix from UMAP-derived cell–landmark connection strengths and
supports four downstream analysis modes within a single analytical
scaffold:

1.  **Differential cell-state density analysis** — landmark-level linear
    modeling of density shifts across conditions
2.  **pePC (partial-effect Principal Component projections)** —
    supervised quantitative sample embedding that isolates variation
    attributable to a specified effect
3.  **plsD (graph-diffused, density contrast-aligned partial least
    squares Decomposition)** — multivariate feature interpretation
    identifying expression programs aligned with density contrasts
4.  **Pseudobulk differential expression** — manifold-informed,
    connection- strength-weighted pseudobulk aggregation with design and
    marker modes

### S4 TDRObj class

The data container is a formal **S4 class** with **12 slots**, providing
a structured representation for landmark-based single-cell analysis
objects. A backward-compatible `$` accessor preserves legacy syntax
(`obj$pbDE`, `obj$markerDE`, etc.) while all internal code uses
`@`-access.

### Multi-backend support

All analysis and plot functions dispatch via S3 generics with methods
for multiple input formats:

| Backend | Format |
|----|----|
| Seurat | v4, v5 (on-disk or in-memory) |
| SingleCellExperiment | On-disk or in-memory |
| H5AnnData | On-disk `.h5ad` files (via anndataR / rhdf5 + BPCells) |
| BPCells | IterableMatrix (on-disk) |
| dgCMatrix | Sparse matrix (in-memory) |
| DelayedMatrix | Converted to BPCells for on-disk efficiency |
| flowSet / cytoset | Flow, mass, and spectral cytometry (flowWorkspace) |

- **[`RunTDR()`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)**:
  Unified entry point that sets up a `TDRObj` from any supported
  container.
- **[`GetTDR()`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md)
  /
  [`SetTDR()`](https://opensource.nibr.com/tinydenseR/reference/SetTDR.md)**:
  Extract / embed the `TDRObj` within a container.
- S3 methods follow a three-step pattern: extract → compute on TDRObj →
  re-embed, so all functions work identically on any backend.

### Manifold-informed pseudobulk DE

[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
supports **two analysis modes**:

| Mode | Purpose |
|----|----|
| **design** | Pseudobulk DE across experimental conditions |
| **marker** | Marker identification comparing cell populations within samples |

Pseudobulk aggregation uses connection-strength-weighted sums over
landmark neighborhoods rather than hard cluster or cell-type membership,
yielding a manifold-informed pseudobulk representation while retaining
samples as the unit of inference.

### plsD: graph-diffused partial least squares decomposition

[`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)
identifies graph-smoothed expression programs maximally aligned with a
fitted density contrast. Key properties:

- Sparse-implicit NIPALS PLS1 implementation (no dense intermediates)
- Concordance-filtered and raw loadings
- Post-hoc sign alignment with density contrast
- Designed for interpretation and hypothesis generation, not formal
  inference

### celltyping: multi-solution support

[`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
accepts **two input modes**:

| Mode  | Input type               | Semantics                            |
|-------|--------------------------|--------------------------------------|
| **A** | Named `list`             | Cluster-to-label mapping             |
| **B** | Named `character` vector | Per-cell labels for the full dataset |

Multiple named celltyping solutions can coexist. Switch active solutions
with
[`set_active_celltyping()`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md)
and list available solutions with
[`list_celltyping_solutions()`](https://opensource.nibr.com/tinydenseR/reference/list_celltyping_solutions.md).

### Reclustering and multi-resolution workflows

[`recluster()`](https://opensource.nibr.com/tinydenseR/reference/recluster.md)
re-runs Leiden community detection with new parameters and refreshes all
clustering-dependent downstream slots. Multiple clustering solutions can
be stored and switched via
[`set_active_clustering()`](https://opensource.nibr.com/tinydenseR/reference/set_active_clustering.md).

### New features

- **[`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md)**:
  Sample-level embeddings (PCA, diffusion-map trajectory, and pePC) from
  the landmark-by-sample density matrix.
- **[`get.cellmap()`](https://opensource.nibr.com/tinydenseR/reference/get.cellmap.md)**:
  Accessor for per-cell mapping data by sample and slot.
- **`sim_trajectory_tdr`**: Bundled simulated scRNA-seq dataset with
  condition-dependent differential abundance for examples.
- **[`simulate_DA_data()`](https://opensource.nibr.com/tinydenseR/reference/simulate_DA_data.md)
  /
  [`simulate_DE_data()`](https://opensource.nibr.com/tinydenseR/reference/simulate_DE_data.md)**:
  Generate synthetic cytometry datasets for benchmarking.

### Vignettes

- **Clinical scRNA-seq analysis**: End-to-end workflow with on-disk
  data.
- **Clinical cytometry analysis**: Flow cytometry with longitudinal
  design.
- **Multi-backend workflows**: Demonstrates all supported input formats.
- **Package architecture**: Object structure and analysis pipeline
  overview.

------------------------------------------------------------------------

## tinydenseR 0.0.2.0001

### On-disk caching for `get.map()`

[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
supports on-disk caching (`.cache.on.disk = TRUE`, default), serializing
large per-sample data to disk to reduce peak RAM usage for large
datasets. All downstream functions read from cache transparently.

Cache management:
[`tdr_cache_info()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_info.md),
[`tdr_cache_validate()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_validate.md),
[`tdr_cache_cleanup()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_cleanup.md).

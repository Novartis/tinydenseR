# tinydenseR 0.0.3.0001

## Removed: `.cl.ct.to.ign` argument from `get.map()`

The `.cl.ct.to.ign` argument has been removed from `get.map()`. This argument
previously allowed excluding a single cluster or cell type from density and
compositional statistics while still mapping those cells. This created internal
inconsistencies between the mapped data and the reported statistics.

**Migration:** Analysts must now perform QC and cell-type/cluster exclusion
**upstream**, before calling `get.map()`. For example, remove unwanted
populations (e.g., erythrocytes in PBMC data) by subsetting cells before
running the pipeline.

Passing `.cl.ct.to.ign` to `get.map()` now produces an informative error
directing users to the upstream-QC workflow.

## Removed: `specDE` and `nmfDE` decomposition methods

`get.specDE()`, `get.nmfDE()`, and their associated plotting functions
(`plotSpecDE()`, `plotSpecDEHeatmap()`, `plotNmfDE()`, `plotNmfDEHeatmap()`)
have been removed from the package. **`get.plsD()`** (partial least squares
decomposition) is the recommended and sole supported decomposition method
in this family. It is substantially more robust and easier to interpret.

The `$specDE` and `$nmfDE` accessor shims on `TDRObj` are retained for
backward compatibility with previously serialized objects.

## S4 TDRObj class (breaking internal change)

The data container has been re-implemented as a formal **S4 class** and
restructured from 17 slots to **12 slots**.  A backward-compatible `$` shim
preserves the old accessor syntax (`obj$pbDE`, `obj$markerDE`, …) while all
internal code now uses `@`-access.

Key slot mapping changes:

| Old path | New slot |
|----------|----------|
| `$graph$clustering` | `@landmark.annot$clustering` |
| `$graph$celltyping` | `@landmark.annot$celltyping` |
| `$graph$nearest.landmarks` | `@cellmap$nearest.lm` |
| `$graph$fuzzy.graph` | `@cellmap$fuzzy.graph` |

## Multi-backend S3 dispatch

All analysis and plot functions now dispatch via S3 generics with methods for
**Seurat**, **SingleCellExperiment**, **HDF5AnnData** (anndataR), **cytoset**
(flowWorkspace), **dgCMatrix**, **DelayedMatrix**, and **IterableMatrix**
(BPCells).

- **`RunTDR()`**: Unified entry point that sets up a `TDRObj` from any
  supported container and runs the full pipeline.
- **`GetTDR()` / `SetTDR()`**: Extract / embed the `TDRObj` within a
  container (e.g., `Seurat@misc$TDRObj`).
- S3 methods follow a three-step pattern: extract → compute on TDRObj →
  re-embed, so all functions (`get.landmarks()`, `get.graph()`, `get.map()`,
  `get.pbDE()`, `get.plsD()`, all `plot*()`
  functions, etc.) work identically on any backend.
- **`get.meta()`** dispatchers added for all backends.
- **`get.features()`**: Renamed from `get.lm.features.stats()`.

### BPCells on-disk backend

`RunTDR.HDF5AnnData`, `RunTDR.DelayedMatrix`, and
`RunTDR.SingleCellExperiment` now use **BPCells** for zero-copy on-disk
matrix access, replacing the previous rhdf5-based subsetting.

### Cytoset support

`RunTDR.cytoset` added for flowWorkspace cytometry containers, with
index-based cell access.

## Unified `get.pbDE()` with design and marker modes

`get.pbDE()` now supports **two analysis modes** via the `.mode` parameter:

| Mode | Purpose | Triggered by |
|------|---------|-------------|
| **design** | Pseudobulk DE across experimental conditions | `.design` provided |
| **marker** | Marker identification comparing cell populations | `.id`/`.id.idx` provided without `.design` |

### Key changes

- **Auto-detection**: `.mode` is inferred from arguments when not specified.
  Providing `.design` selects design mode; providing `.id`/`.id.idx` without
  `.design` selects marker mode.
- **New parameters**: `.mode`, `.id2`, `.id2.idx`, `.result.name`.
- **Deprecated aliases**: `.population.name` and `.comparison.name` are
  accepted with a deprecation warning; use `.result.name` instead.
- **`get.markerDE()` soft-deprecated**: Now a thin wrapper dispatching to
  `get.pbDE(.mode = "marker")` with parameter mapping (`.id1` → `.id`,
  `.id1.idx` → `.id.idx`, `.comparison.name` → `.result.name`).
  Will be removed in a future release.
- **`get.dea()` and `get.marker()`**: Updated deprecated wrappers to call
  `get.pbDE()` directly.
- **Internal helpers**: Extracted `.tdr_resolve_cell_idx()`,
  `.tdr_pseudobulk_aggregate_rna()`, and `.tdr_pseudobulk_aggregate_cyto()`
  to reduce duplication between modes.

### Storage (unchanged)

- Design mode → `$pbDE[[.model.name]][[.result.name]]`
- Marker mode → `$markerDE[[.model.name]][[.result.name]]`

### Migration

```r
# Before (still works, with deprecation warning):
get.markerDE(obj, .id1 = "cluster.3", .comparison.name = "c3_markers")

# After:
get.pbDE(obj, .mode = "marker", .id = "cluster.3", .result.name = "c3_markers")
```

## `celltyping()` Mode B: per-cell labels

`celltyping()` now accepts **two mutually exclusive input modes**,
dispatched automatically based on the type of `.celltyping.map`:

| Mode | Input type | Semantics |
|------|-----------|-----------|
| **A** (existing) | Named `list` | Cluster-to-label mapping; all clusters must be covered |
| **B** (new) | Named `character` vector | Per-cell labels for the full dataset; landmarks are matched automatically |

```r
# cell_labels: named character vector, names = original cell IDs,
# values = cell-type labels, covering ALL cells across all samples.
obj <- celltyping(obj, cell_labels)
```

A new `$landmark.annot$celltyping$mode` slot records `"cluster_map"` (A)
or `"cell_labels"` (B).  Late celltyping after `get.map()` is also now
supported.

## Unified cyto landmark naming

Cytometry landmark rownames now use the same `paste0(sample, "_", cell_id)`
convention as RNA, preventing cross-sample name collisions and enabling
Mode B cell matching.

## Rename: `plsDE` → `plsD`

All references to `plsDE` (partial least squares differential expression)
have been renamed to **`plsD`** (partial least squares decomposition) to
better reflect the method's purpose.  Affected: `get.plsD()`,
`plotPlsD()`, `plotPlsDHeatmap()`, `$plsD` accessor, vignettes, tests.

### Other `plsD` improvements

- Sparse-implicit NIPALS PLS1 — eliminates ~3.2 GB dense intermediates.
- Raw and concordance-filtered loadings now exposed in results.
- Removed `y.weighted.loadings` (sign-direction issue).
- Fixed model name in `plotPlsDHeatmap` annotation.

## Cache infrastructure v2

Ephemeral temp-directory storage with session-end cleanup, replacing the
previous persistent cache.  Internal schema field updated to `v2`.

## Label transfer refactor

Label transfer logic extracted from `stats.model.R` into a dedicated
`label.transfer.R` module.  Default confidence thresholds adjusted.
Seurat v5 cell access and landmark extraction improved.

## New features

- **`get.cellmap()`**: Accessor for per-cell mapping data by sample and slot.
- **`.label.substr.rm`** parameter for customizing DE heatmap coefficient
  labels in `plotPbDE()`, `plotMarkerDE()`, and related heatmap functions.
- **`sim_trajectory_tdr`**: Bundled simulated scRNA-seq dataset with
  condition-dependent differential abundance for vignettes and examples.
- **COVID-19 PBMC analysis vignette**: End-to-end clinical scRNA-seq workflow.
- **Multi-backend workflows vignette** and **package architecture vignette**.

## Bug fixes

- `.warn_stale_de_results()` now correctly iterates both nesting levels
  (model → comparison) for marker DE results.
- `get.markerDE()` `.id2.idx` parameter now correctly takes priority over
  `.id2` when both are supplied.
- `plotPlsDHeatmap` uses correct model name in annotation.
- README: `update = FALSE` added to `BiocManager::install()` to avoid
  interactive prompts.

---

# tinydenseR 0.0.2.0001

## On-disk caching for `get.map()` (new feature)

`get.map()` now defaults to **on-disk caching** (`.cache.on.disk = TRUE`),
serializing four large per-sample slots to disk instead of holding them
in memory:

- `clustering$ids`
- `celltyping$ids`
- `nearest.landmarks`
- `fuzzy.graph`

This dramatically reduces peak RAM usage for large single-cell datasets.
All downstream functions (`goi.summary()`, `get.lm()`, `get.pbDE()`,
`get.plsD()`) read from the cache
transparently.

### Cache management functions

| Function | Purpose |
|---|---|
| `tdr_cache_info()` | Print a human-readable summary of cache state (directory, size, slots, samples) |
| `tdr_cache_validate()` | Verify that all cached files exist (and optionally check MD5 checksums with `.verify.checksum = TRUE`) |
| `tdr_cache_cleanup()` | Delete the cache directory and strip metadata from the object |

### Hardening

- **MD5 checksums** are stored per file at write time and can be verified
  with `tdr_cache_validate(.verify.checksum = TRUE)` or
  `.tdr_cache_read(meta, verify = TRUE)`.
- **Schema versioning** — cache metadata carries a `schema_v` field; reads
  error immediately when a version mismatch is detected.
- **Crash-safe writes** — files are written to a `.rds.tmp` temp file and
  atomically renamed; orphaned temp files are cleaned up automatically.
- **File-level locking** — when the `filelock` package is installed
  (added to `Suggests`), a lock is acquired during the atomic rename to
  prevent concurrent-write clobbering.
- **Automatic validation** — `get.lm()`, `get.pbDE()`, and
  `.validate.DE.inputs()` (used by `get.plsD()`) silently validate the cache on entry.

### Disabling caching

Pass `.cache.on.disk = FALSE` to `get.map()` to keep all data in memory
(original behaviour).

# Identify landmarks via leverage score sampling

Selects representative "landmark" cells from the full dataset using a
two-pass leverage score sampling strategy. Landmarks capture the
representative cells while being computationally tractable for
downstream graph construction and mapping.

## Usage

``` r
get.landmarks(x, ...)

# S3 method for class 'TDRObj'
get.landmarks(
  x,
  .source = NULL,
  .verbose = TRUE,
  .seed = 123,
  .nHVG = 5000,
  .nPC = 30,
  .exc.vdj.mito.ribo.genes.from.hvg = TRUE,
  .force.in = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  initialized with `setup.tdr.obj`.

- ...:

  Additional arguments passed to methods.

- .source:

  The raw data object for non-file backends. `NULL` (default) for the
  files backend; otherwise a Seurat, SingleCellExperiment, or anndataR
  AnnData object. Used by `.get_sample_matrix()` to retrieve per-sample
  expression matrices.

- .verbose:

  Logical, print progress messages. Default TRUE.

- .seed:

  Integer for reproducibility. Default 123.

- .nHVG:

  Integer, number of highly variable genes to select for RNA data.
  Default 5000. Ignored for cytometry. Higher values capture more
  variation but increase computation time.

- .nPC:

  Integer, number of principal components for dimensionality reduction.
  Default 30. Must be less than the number of cells in smallest sample.

- .exc.vdj.mito.ribo.genes.from.hvg:

  Logical, whether to exclude V(D)J variable-region genes
  (`TR[ABDG][VDJ]`, `IG[KHL][VDJ]`), mitochondrial genes (`MT-`), and
  ribosomal protein genes (`RPS/RPL/RPLP/RPSA`) from HVG selection (RNA
  only). Default TRUE. Constant-region genes (e.g. TRAC, IGHG, IGKC) are
  intentionally retained as they carry cell-identity signal. Recommended
  to avoid technical/biological noise dominating variation.

- .force.in:

  Character vector of gene names to force into the feature set
  regardless of variance (RNA only). Useful for known markers. Default
  NULL.

## Value

Updated `.tdr.obj` with populated fields:

- `$raw.landmarks`:

  Raw counts matrix for landmarks (landmarks × features)

- `$landmarks`:

  Processed landmark expression on selected features (landmarks ×
  features):

  - RNA: PCA-reconstructed denoised expression (log2-scale after library
    size normalization)

  - Cytometry: Original marker values on selected markers

- `$scaled.landmarks`:

  Z-scored landmark expression (landmarks × features, for
  visualization/heatmaps)

- `$pca`:

  List containing PCA results:

  - `$embed` - PC coordinates for landmarks (landmarks × PCs)

  - `$rotation` - Feature loadings (features × PCs)

  - `$center` - Feature means (length = \# features)

  - `$scale` - Feature standard deviations (length = \# features)

  - `$sdev` - Standard deviations of PCs (length = \# PCs)

  - `$HVG` - Selected feature names (character vector)

  If Harmony used: `$embed` and `$rotation` are
  Harmony-corrected/approximated

- `$integration$harmony.obj`:

  Symphony reference object (if `.harmony.var` specified), used for
  batch-corrected mapping of query cells

## Details

**Two-pass landmark selection algorithm:**

**Pass 1 - Initial sampling:**

1.  For RNA: normalize, log-transform, select top HVGs per sample

2.  For cytometry: use specified markers

3.  Compute sample-specific PCA

4.  Calculate leverage scores (sum of squared PC loadings per cell)

5.  Sample landmarks proportionally to leverage scores

**Pass 2 - Refinement:**

1.  Pool landmarks from all samples

2.  Compute dataset-wide PCA on pooled landmarks

3.  Project ALL cells onto this shared PC space

4.  Recalculate leverage scores using shared PCA

5.  Resample landmarks with improved scores (final set)

This two-pass approach ensures landmarks are representative of global
(not just sample-specific) variation patterns. Leverage score sampling
prioritizes cells in high-variance regions while maintaining diversity.

**Optional Harmony integration:** If `.harmony.var` was specified in
`setup.tdr.obj`, performs batch correction on landmark PC/SVD
embeddings. This creates a Symphony reference object for mapping query
cells in a batch-corrected space. Supported for both RNA and cytometry
assay types. For cytometry, Harmony corrects batch effects in the full
SVD embedding of the marker matrix (one dimension per marker). Cytometry
data should be pre-transformed (e.g., arcsinh, logicle) before entering
the tinydenseR pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{
# Typical workflow (from README)
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks(.nHVG = 500, .nPC = 3)

# RNA with more PCs and custom HVGs
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks(.nPC = 50, .nHVG = 3000)

# Force specific markers into feature set
lm.cells <- get.landmarks(lm.cells, 
                          .force.in = c("CD3D", "CD4", "CD8A"))
} # }
```

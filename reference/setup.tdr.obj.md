# Initialize tinydenseR object for landmark-based analysis

Creates and validates the main tinydenseR object structure that will
hold expression data, metadata, and analysis results throughout the
workflow. This is the required first step before landmark
identification.

`setup.lm.obj()` is deprecated; use `setup.tdr.obj()` instead.

## Usage

``` r
setup.tdr.obj(
  .cells,
  .meta,
  .markers = NULL,
  .harmony.var = NULL,
  .assay.type = "cyto",
  .celltype.vec = NULL,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 }
)

setup.lm.obj(...)
```

## Arguments

- .cells:

  A named list of file paths (character strings) pointing to RDS files,
  each containing one expression matrix per sample. For RNA: sparse
  matrix (dgCMatrix) with genes as rows and cells as columns. For
  cytometry: matrix with cells as rows and markers as columns. List
  names become sample identifiers.

- .meta:

  A data.frame with sample-level metadata. Rownames must match names in
  `.cells` exactly. Used for batch correction and downstream modeling.

- .markers:

  Character vector of marker names to use for landmark identification.
  Only applicable for `.assay.type = "cyto"`. If NULL, uses all markers
  from the first sample. Ignored for RNA data (uses HVG selection
  instead).

- .harmony.var:

  Character vector of column names from `.meta` to use for Harmony batch
  correction. Supported for both `"RNA"` and `"cyto"` assay types. For
  cytometry, Harmony operates on the SVD embedding of the (centered,
  scaled) marker matrix. If NULL, no batch correction is performed.

- .assay.type:

  Character string: "cyto" for cytometry (default) or "RNA" for
  scRNA-seq. Determines normalization strategy and feature selection
  approach.

- .celltype.vec:

  Optional named character vector mapping cell IDs to cell type labels.

- .verbose:

  Logical, whether to print progress messages. Default TRUE.

- .seed:

  Integer for random seed to ensure reproducibility. Default 123.

- .prop.landmarks:

  Numeric between 0 and 1 specifying proportion of cells to sample as
  landmarks. Default 0.1 (10%). Total landmarks capped at about 5000
  regardless.

- .n.threads:

  Integer for parallel processing. Default automatically detects maximum
  available threads (using BLAS settings on HPC, or `detectCores()`
  locally).

- ...:

  Arguments passed to `setup.tdr.obj()`.

## Value

A list object with initialized structure containing:

- `$cells`:

  Input file paths

- `$metadata`:

  Sample metadata with added cell count columns

- `$config`:

  Configuration parameters:

  `$key`

  :   Named vector mapping future landmarks to samples

  `$sampling`

  :   Sampling parameters including `n.cells`, `n.perSample`

  `$assay.type`

  :   Assay type ("cyto" or "RNA")

  `$markers`

  :   Marker list (cytometry only)

  `$n.threads`

  :   Number of threads for parallel processing

- `$integration`:

  Integration/batch correction results:

  `$harmony.var`

  :   Batch variables (if provided)

  `$harmony.obj`

  :   Symphony reference object (populated by `get.landmarks`)

- Empty slots:

  landmarks, scaled.landmarks, raw.landmarks, pca, graph, map, etc. -
  populated by downstream functions

## Details

This function:

1.  Validates input data structure and compatibility

2.  Calculates the number of landmarks to sample per sample (max 5000
    total)

3.  Creates a "key" vector mapping landmarks to samples

4.  Initializes empty slots for downstream analyses (PCA, graph,
    mapping)

5.  Performs quality checks (warns if sample sizes vary \>10-fold)

The landmark sampling strategy aims for proportional representation
across samples while capping total landmarks at 5000 for computational
efficiency. Large samples are capped to prevent domination and ensure
adequate representation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Prepare data (from README workflow)
.meta <- get.meta(.obj = sim_trajectory, .sample.var = "Sample")
.cells <- get.cells(.exprs = sim_trajectory, 
                    .meta = .meta,
                    .sample.var = "Sample")

# Basic setup for RNA data
tdr.obj <- setup.tdr.obj(
  .cells = .cells,
  .meta = .meta,
  .assay.type = "RNA",
  .prop.landmarks = 0.15
)

# Cytometry workflow with marker selection
tdr.obj <- setup.tdr.obj(
  .cells = .cells,
  .meta = .meta,
  .markers = c("CD3", "CD4", "CD8", "CD19"),
  .assay.type = "cyto"
)

# RNA workflow with batch correction
tdr.obj <- setup.tdr.obj(
  .cells = .cells,
  .meta = .meta,
  .harmony.var = "batch",
  .assay.type = "RNA"
)
} # }
```

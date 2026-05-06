# Create .cells Object from SingleCellExperiment Object

Converts a SingleCellExperiment (SCE) object into tinydenseR's .cells
format. SCE is the standard Bioconductor container for single-cell data.
Each sample's cells are extracted and saved as temporary RDS files.

## Usage

``` r
get.cells.SCE(
  .sce.obj,
  .meta,
  .sample.var,
  .assay = "counts",
  .min.cells.per.sample = 10,
  .compress = FALSE,
  .verbose = TRUE
)
```

## Arguments

- .sce.obj:

  SingleCellExperiment object containing expression data in assays.

- .meta:

  Data frame with sample-level metadata. Rownames must be sample IDs
  matching values in `.sce.obj@colData[[.sample.var]]`. Required for
  filtering valid samples.

- .sample.var:

  Character: column name in `colData(.sce.obj)` identifying which sample
  each cell belongs to (e.g., "sample_id", "Sample").

- .assay:

  Character: assay name to extract. Default "counts". Can also be
  "logcounts", "normcounts", or any assay in `assayNames(.sce.obj)`.

- .min.cells.per.sample:

  Integer: minimum cells required per sample. Default 10.

- .compress:

  Logical: compress RDS files? Default FALSE for faster I/O.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Named list of file paths to temporary RDS files, one per sample. Only
includes samples meeting minimum cell threshold and present in `.meta`.

## Details

SingleCellExperiment is the Bioconductor standard for single-cell data,
used by many analysis packages. This function extracts cells
sample-by-sample and converts to sparse dgCMatrix format for
memory-efficient storage.

The function verifies that `colnames(.sce.obj)` exist, as these are
required for proper cell tracking.

## See also

[`get.cells`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
for automatic format detection,
[`get.cells.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat.md)
for Seurat objects,
[`setup.tdr.obj`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)
for next workflow step

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
trajectory_data <- fetch_trajectory_data()
sim_trajectory.meta <- trajectory_data$meta |>
  dplyr::select(Condition, Replicate, Sample) |>
  unique()
rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
sim_trajectory <- trajectory_data$SCE

# Convert SCE directly to .cells format
cells <- get.cells.SCE(.sce.obj = sim_trajectory,
                       .meta = sim_trajectory.meta,
                       .sample.var = "Sample")

# Use in tinydenseR workflow
lm.obj <- setup.tdr.obj(.cells = cells,
                       .meta = sim_trajectory.meta,
                       .assay.type = "RNA")
} # }
```

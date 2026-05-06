# Create .cells Object from Seurat v4 Object

Converts a Seurat v4 object into tinydenseR's .cells format. Unlike v5,
Seurat v4 stores all data in a single slot per assay. Each sample's
cells are extracted and saved as temporary RDS files.

## Usage

``` r
get.cells.Seurat(
  .seurat.obj,
  .meta,
  .sample.var,
  .assay = "RNA",
  .slot = "counts",
  .min.cells.per.sample = 10,
  .compress = FALSE,
  .verbose = TRUE
)
```

## Arguments

- .seurat.obj:

  Seurat v4 object containing expression data. Must use Assay class (not
  Assay5).

- .meta:

  Data frame with sample-level metadata. Rownames must be sample IDs
  matching values in `.seurat.obj@meta.data[[.sample.var]]`. Required
  for filtering valid samples.

- .sample.var:

  Character: column name in `.seurat.obj@meta.data` identifying which
  sample each cell belongs to (e.g., "sample_id", "Sample",
  "orig.ident").

- .assay:

  Character: assay name to extract. Default "RNA".

- .slot:

  Character: data slot to extract. Default "counts". Can also be "data"
  (normalized) or "scale.data" (scaled).

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

Seurat v4 uses a simpler assay structure than v5, with all data stored
in a single matrix per slot. This function extracts cells
sample-by-sample and converts to sparse dgCMatrix format for
memory-efficient storage.

Compatible with Seurat versions \< 5.0.0 or when using legacy Assay
class in v5.

## See also

[`get.cells`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
for automatic format detection,
[`get.cells.Seurat5`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat5.md)
for Seurat v5 objects,
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

# Create Seurat v4 object
counts.A_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
                 trajectory_data$SCE$Sample == "A_R1"]
counts.B_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
                 trajectory_data$SCE$Sample == "B_R1"]
combined.counts <- cbind(counts.A_R1, counts.B_R1)

cell.meta <- data.frame(
  sample_id = c(rep("A_R1", ncol(counts.A_R1)),
                rep("B_R1", ncol(counts.B_R1))),
  row.names = colnames(combined.counts)
)

seurat.obj <- CreateSeuratObject(counts = combined.counts,
                                 meta.data = cell.meta)

# Convert to .cells format
cells <- get.cells.Seurat(.seurat.obj = seurat.obj,
                          .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
                          .sample.var = "sample_id")

# Use in tinydenseR workflow
lm.obj <- setup.tdr.obj(.cells = cells,
                       .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
                       .assay.type = "RNA")
} # }
```

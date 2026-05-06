# Create .cells Object from Seurat v5 Object

Converts a Seurat v5 object with layered data into tinydenseR's .cells
format. Can handle Seurat's v5 assay structure where data is split
across multiple layers if all cells in a single sample are contained
within a single layer (for example, samples 1:3 in layer `x` and sample
4:6 in layer `y`). If cells within a sample are split across layers, an
error is thrown.

## Usage

``` r
get.cells.Seurat5(
  .seurat.obj,
  .meta,
  .sample.var,
  .assay = "RNA",
  .layer.pattern = "counts",
  .layer.name.source = NULL,
  .min.cells.per.sample = 10,
  .compress = FALSE,
  .verbose = TRUE,
  .dest.path = NULL
)
```

## Arguments

- .seurat.obj:

  Seurat v5 object containing layered RNA data. Must use Assay5 class.

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

- .layer.pattern:

  Character: pattern to match layer names. Default "counts" matches
  layers containing count data (e.g., "counts.Sample1",
  "counts.Sample2"). Uses `fixed = TRUE` matching.

- .layer.name.source:

  Character: optional metadata column mapping each cell to a layer
  label. If provided, only cells where this column matches the current
  layer label are included for that sample. This allows for more
  flexible layer naming schemes. Default NULL (no additional filtering).

- .min.cells.per.sample:

  Integer: minimum cells required per sample. Samples with fewer cells
  are excluded. Default 10.

- .compress:

  Logical: compress RDS files? Default FALSE for faster I/O.

- .verbose:

  Logical: print progress messages? Default TRUE.

- .dest.path:

  Character: optional directory path to save `.cells` RDS files. If
  NULL, files are saved as temporary files.

## Value

Named list of file paths to temporary RDS files, one per sample. Only
includes samples meeting minimum cell threshold and present in `.meta`.

## Details

Seurat v5 uses a layered assay structure where data for different
samples may be stored in separate layers. This function:

1.  Identifies all layers matching `.layer.pattern`

2.  Finds common genes across all layers

3.  Validates each sample maps to one layer (when multiple layers are
    present)

4.  Extracts and combines sample cells from matching layer data into one
    sparse matrix per sample

5.  Saves as temporary RDS files for memory-efficient processing

Only samples appearing in both the Seurat object and `.meta` rownames
are processed.

## See also

[`get.cells`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
for automatic format detection,
[`get.cells.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat.md)
for Seurat v4 objects,
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

# Create Seurat v5 object from SCE data
# Extract counts for two samples
counts.A_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
                 trajectory_data$SCE$Sample == "A_R1"]
counts.B_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
                 trajectory_data$SCE$Sample == "B_R1"]
combined.counts <- cbind(counts.A_R1, counts.B_R1)

# Build cell metadata
cell.meta <- data.frame(
  sample_id = c(rep("A_R1", ncol(counts.A_R1)),
                rep("B_R1", ncol(counts.B_R1))),
  row.names = colnames(combined.counts)
)

# Create Seurat v5 object
seurat.obj <- CreateSeuratObject(counts = combined.counts,
                                 meta.data = cell.meta)
seurat.obj[["RNA"]] <- as(seurat.obj[["RNA"]], "Assay5")
seurat.obj <- JoinLayers(seurat.obj)

# Convert to .cells format
cells <- get.cells.Seurat5(.seurat.obj = seurat.obj,
                           .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
                           .sample.var = "sample_id")

# Use in tinydenseR workflow
lm.obj <- setup.tdr.obj(.cells = cells,
                       .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
                       .assay.type = "RNA")
} # }
```

# Create .cells Object from Count Matrices

Converts a list of count matrices into tinydenseR's standardized .cells
format. Each matrix represents one sample and is saved as a temporary
RDS file for memory-efficient processing. Use this when you already have
count data loaded in R memory.

## Usage

``` r
get.cells.list.mat(.count.mat.list, .meta, .compress = FALSE, .verbose = TRUE)
```

## Arguments

- .count.mat.list:

  Named list of count matrices, one per sample. Names must match
  `rownames(.meta)`. Each matrix can be dense, sparse (dgCMatrix), or
  similar matrix-like object with features as rows and cells as columns.

- .meta:

  Data frame with sample-level metadata. Rownames must match
  `names(.count.mat.list)`. Should include experimental variables like
  Condition, Replicate, etc.

- .compress:

  Logical: compress RDS files? Default FALSE for faster I/O. Set TRUE to
  save disk space for large datasets.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Named list of file paths to temporary RDS files, one per sample.
Structure suitable for `setup.tdr.obj(.cells = ...)`.

## Details

This function creates temporary RDS files containing each sample's count
matrix. These files persist for the R session and allow tinydenseR to
process data without loading all samples into memory simultaneously.

All matrices in `.count.mat.list` should have:

- Same feature names (rownames) across samples

- Cell barcodes as column names

- Raw or normalized counts (specify via
  `setup.tdr.obj(.assay.type = ...)`)

## See also

[`get.cells`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
for automatic format detection,
[`get.cells.SCE`](https://opensource.nibr.com/tinydenseR/reference/get.cells.SCE.md)
for SingleCellExperiment input,
[`setup.tdr.obj`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)
for next step in workflow

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example trajectory data
trajectory_data <- fetch_trajectory_data()
sim_trajectory.meta <- trajectory_data$meta
sim_trajectory <- trajectory_data$SCE
  
# Prepare sample-level metadata
sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
  unique()
rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample

# Extract count matrices for two samples
count.matrices <- list(
  A_R1 = SingleCellExperiment::counts(sim_trajectory)[, 
           sim_trajectory$Sample == "A_R1"],
  B_R1 = SingleCellExperiment::counts(sim_trajectory)[, 
           sim_trajectory$Sample == "B_R1"]
)

# Create .cells object
cells <- get.cells.list.mat(.count.mat.list = count.matrices,
                            .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ])

# Use in tinydenseR workflow
lm.obj <- setup.tdr.obj(.cells = cells,
                       .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
                       .assay.type = "RNA")
} # }
```

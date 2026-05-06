# Simulate differential expression (DE) flow cytometry data

Generates simulated flow cytometry data with differential expression of
Marker2 in a target cell type across treatment groups, batch effects,
and multiple shift magnitudes. Writes per-sample FCS files and returns
sample- and cell-level metadata.

## Usage

``` r
simulate_DE_data(
  groups = c("Baseline", "Activation"),
  batches = c("Batch1", "Batch2"),
  sd_shifts = c(0.5, 1, 2),
  samples_per_group = 6,
  mean_cells = 50000,
  sd_cells = 500,
  seed = 42,
  output_dir = file.path(tempdir(), "sim_DE_fcs")
)
```

## Arguments

- groups:

  Character vector of treatment group names.

- batches:

  Character vector of batch names.

- sd_shifts:

  Numeric vector of standard deviation shifts applied to Marker2 in the
  non-Baseline group for target cells.

- samples_per_group:

  Integer. Number of samples per group per shift.

- mean_cells:

  Integer. Mean number of cells per sample.

- sd_cells:

  Numeric. Standard deviation of cell count.

- seed:

  Integer. Random seed for reproducibility.

- output_dir:

  Character. Directory path where FCS files are written.

## Value

A list with two elements:

- sample_meta:

  A `data.frame` with one row per sample and columns `Sample`,
  `Treatment`, `Batch`, `SD_Shift`, `fcs_path`.

- cell_meta:

  A `data.frame` with one row per cell and columns `Sample`,
  `Treatment`, `Batch`, `SD_Shift`, `CellType`.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_DE_data(
  samples_per_group = 2,
  mean_cells = 1000,
  sd_cells = 100,
  output_dir = file.path(tempdir(), "sim_DE_fcs")
)
head(sim$sample_meta)
head(sim$cell_meta)
} # }
```

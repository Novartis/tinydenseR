# Simulate differential abundance (DA) flow cytometry data

Generates simulated flow cytometry data with differential abundance of a
target cell type across treatment groups, batch effects, and multiple
abundance settings. Writes per-sample FCS files and returns sample- and
cell-level metadata.

## Usage

``` r
simulate_DA_data(
  groups = c("Baseline", "Depletion"),
  batches = c("Batch1", "Batch2"),
  settings = c(0.005, 0.05, 0.5),
  samples_per_group = 6,
  mean_cells = 50000,
  sd_cells = 500,
  seed = 42,
  output_dir = file.path(tempdir(), "sim_DA_fcs")
)
```

## Arguments

- groups:

  Character vector of treatment group names.

- batches:

  Character vector of batch names.

- settings:

  Numeric vector of proportions for the Baseline group. The non-Baseline
  group uses half the proportion.

- samples_per_group:

  Integer. Number of samples per group per setting.

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
  `Treatment`, `Batch`, `Setting`, `fcs_path`.

- cell_meta:

  A `data.frame` with one row per cell and columns `Sample`,
  `Treatment`, `Batch`, `Setting`, `CellType`.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_DA_data(
  samples_per_group = 2,
  mean_cells = 1000,
  sd_cells = 100,
  output_dir = file.path(tempdir(), "sim_DA_fcs")
)
head(sim$sample_meta)
head(sim$cell_meta)
} # }
```

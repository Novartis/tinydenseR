# Fetch trajectory simulation dataset from miloR

Downloads and extracts the sim_trajectory dataset from the miloR package
GitHub repository. This dataset contains simulated single-cell
trajectory data with two conditions (A and B) across three replicates
each.

## Usage

``` r
fetch_trajectory_data()
```

## Source

Data originally from: Dann, E., Henderson, N.C., Teichmann, S.A. et al.
Differential abundance testing on single-cell data using k-nearest
neighbor graphs. Nat Biotechnol (2021).
https://doi.org/10.1038/s41587-021-01033-z

miloR GitHub repository: https://github.com/MarioniLab/miloR Direct data
link:
https://github.com/MarioniLab/miloR/blob/bdecaebeb5595545f9b9c8f2defb321519d98a70/data/sim_trajectory.RData

## Value

A list containing:

- meta:

  Sample metadata with Condition, Replicate, and Sample columns

- SCE:

  SingleCellExperiment object containing the count data

## Note

The downloaded data is subject to GPL v3 license terms from miloR
package

## Examples

``` r
if (FALSE) { # \dontrun{
# Fetch the trajectory data
sim_trajectory_data <- fetch_trajectory_data()

# Extract components
sim_trajectory.meta <- sim_trajectory_data$meta
sim_trajectory <- sim_trajectory_data$SCE
} # }
```

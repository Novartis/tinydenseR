# Detect if running on High-Performance Computing cluster (internal)

Checks environment variables to determine if code is executing on an HPC
system with a job scheduler (SLURM, LSF, PBS, SGE, OAR). Used to
optimize thread detection strategy.

## Usage

``` r
is.hpc()
```

## Value

Logical: TRUE if HPC scheduler detected, FALSE otherwise

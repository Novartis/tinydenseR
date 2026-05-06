# Shared nearest neighbors via fast Jaccard index calculation

Computes a shared nearest neighbor (SNN) graph by calculating Jaccard
similarity between cells based on their k-nearest neighbors. The Jaccard
index J(A,B) = \|A∩B\| / \|A∪B\| measures the overlap of neighbor sets.
This implementation uses sparse matrix operations for efficiency.

## Usage

``` r
fast.jaccard.r(.adj.matrix, .prune = 1/15)
```

## Arguments

- .adj.matrix:

  A non-symmetric sparse adjacency matrix from `get.adj.matrix`, where
  entry (j,i)=1 indicates j is among i's nearest neighbors (neighbors in
  rows, cells in columns).

- .prune:

  A numeric threshold (default 1/15 ≈ 0.067) for pruning weak edges.
  Jaccard values below this are set to zero for sparsity.

## Value

A symmetric sparse matrix of Jaccard indices representing the SNN graph.
Higher values indicate greater neighbor overlap. Entry (i,j) gives the
Jaccard similarity between cell i and cell j based on their shared
k-nearest neighbors. Matrix is pruned to remove weak connections and
forced symmetric.

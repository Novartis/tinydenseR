# Sparse matrix representation of nearest neighbors

Converts a nearest neighbor index matrix into a sparse adjacency matrix.
The resulting matrix is non-symmetric (i can be neighbor of j without j
being neighbor of i). Note that the matrix orientation places neighbors
in rows and source cells in columns.

## Usage

``` r
get.adj.matrix(.nn.idx)
```

## Arguments

- .nn.idx:

  A matrix where each row represents a cell and each column represents a
  nearest neighbor index. Element (i,k) is the index of the k-th nearest
  neighbor of cell i.

## Value

A non-symmetric sparse matrix (dgCMatrix) of dimensions n×n where n is
the number of cells. Entry (j,i) = 1 if j is among i's nearest
neighbors, 0 otherwise. That is, row j contains 1s in columns
corresponding to cells that have j as a neighbor.

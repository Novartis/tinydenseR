# Leiden clustering with straggler absorption

Performs community detection on a similarity graph using the Leiden
algorithm, followed by reassignment of small "straggler" clusters to
better-connected neighbors. This improves cluster quality by preventing
tiny, isolated groups.

## Usage

``` r
leiden.cluster(
  .tdr.obj = NULL,
  .sim.matrix,
  .resolution.parameter,
  .small.size = floor((nrow(x = .sim.matrix)/200)),
  .verbose = TRUE,
  .seed = 123
)
```

## Arguments

- .tdr.obj:

  Optional tinydenseR object from
  [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md)
  and
  [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md).
  If provided, uses Laplacian Eigenmap embedding for initialization (if
  computed successfully), otherwise falls back to PCA embedding. If
  NULL, computes 3D embedding from `.sim.matrix` via truncated SVD.

- .sim.matrix:

  Square similarity matrix (typically Jaccard similarity of k-NN
  graphs). Higher values indicate stronger connections between nodes.
  Self-loops (diagonal) are automatically removed.

- .resolution.parameter:

  Numeric controlling cluster granularity in Leiden CPM algorithm.
  Higher values (e.g., 0.0005-0.005) yield more/smaller clusters. Lower
  values yield fewer/larger clusters. Optimal value depends on data
  scale.

- .small.size:

  Integer threshold for straggler clusters. Clusters with fewer cells
  are absorbed into neighbors. Default: 0.5% of total nodes
  (floor(n/200)).

- .verbose:

  Logical, print progress and cluster counts. Default TRUE.

- .seed:

  Integer for reproducibility (affects k-means and Leiden). Default 123.

## Value

Factor vector of cluster assignments with levels ordered by appearance.
Labels are formatted as "cluster.01", "cluster.02", etc.

## Details

**Algorithm steps:**

1.  Initialize with k-means (up to 25 centers) on embedding to avoid
    singletons

2.  Run Leiden community detection with CPM objective function

3.  Identify "straggler" clusters smaller than `.small.size`

4.  For each straggler, compute mean connectivity to all keeper clusters

5.  Reassign straggler cells to the most connected keeper cluster

6.  Relabel clusters sequentially

**Initialization strategy:** Uses Laplacian Eigenmap embedding (if
computed and available in `.tdr.obj$graph$LE$embed`), otherwise PCA
embedding (up to first 3 components), or computes 3D embedding from
similarity matrix via truncated SVD (`irlba`) when `.tdr.obj` is NULL.
K-means initialization helps Leiden start from reasonable communities
rather than singleton nodes, improving stability.

**Straggler absorption:** Small clusters (\< 0.5% of total by default)
often represent noise or transitional states. Rather than keeping them
as separate clusters, they're merged with the most similar neighboring
cluster based on mean edge weights. For sparse similarity matrices, the
mean is computed by normalizing the sum of stored values by the full
submatrix size to properly account for implicit zeros.

## Examples

``` r
if (FALSE) { # \dontrun{
# After building similarity matrix
clusters <- leiden.cluster(
  .tdr.obj = lm.obj,
  .sim.matrix = jaccard_sim,
  .resolution.parameter = 0.001,
  .small.size = 10  # Merge clusters < 10 cells
)
table(clusters)
} # }
```

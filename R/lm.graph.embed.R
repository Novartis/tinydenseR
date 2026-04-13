#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# Licensed under the MIT License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# https://www.mit.edu/~amini/LICENSE.md
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#####

#' Leiden clustering with straggler absorption
#'
#' Performs community detection on a similarity graph using the Leiden algorithm,
#' followed by reassignment of small "straggler" clusters to better-connected
#' neighbors. This improves cluster quality by preventing tiny, isolated groups.
#'
#' @details
#' **Algorithm steps:**
#' \enumerate{
#'   \item Initialize with k-means (up to 25 centers) on embedding to avoid singletons
#'   \item Run Leiden community detection with CPM objective function
#'   \item Identify "straggler" clusters smaller than \code{.small.size}
#'   \item For each straggler, compute mean connectivity to all keeper clusters
#'   \item Reassign straggler cells to the most connected keeper cluster
#'   \item Relabel clusters sequentially
#' }
#' 
#' **Initialization strategy:**
#' Uses Laplacian Eigenmap embedding (if computed and available in \code{.tdr.obj$graph$LE$embed}),
#' otherwise PCA embedding (up to first 3 components), or computes 3D embedding from
#' similarity matrix via truncated SVD (\code{irlba}) when \code{.tdr.obj} is NULL.
#' K-means initialization helps Leiden start from reasonable communities rather than
#' singleton nodes, improving stability.
#' 
#' **Straggler absorption:**
#' Small clusters (< 0.5% of total by default) often represent noise or transitional
#' states. Rather than keeping them as separate clusters, they're merged with the
#' most similar neighboring cluster based on mean edge weights. For sparse similarity
#' matrices, the mean is computed by normalizing the sum of stored values by the full
#' submatrix size to properly account for implicit zeros.
#'
#' @param .tdr.obj Optional tinydenseR object from \code{get.landmarks()} and
#'   \code{get.graph()}. If provided, uses Laplacian Eigenmap embedding for
#'   initialization (if computed successfully), otherwise falls back to PCA embedding.
#'   If NULL, computes 3D embedding from \code{.sim.matrix} via truncated SVD.
#' @param .sim.matrix Square similarity matrix (typically Jaccard similarity of
#'   k-NN graphs). Higher values indicate stronger connections between nodes.
#'   Self-loops (diagonal) are automatically removed.
#' @param .resolution.parameter Numeric controlling cluster granularity in Leiden
#'   CPM algorithm. Higher values (e.g., 0.0005-0.005) yield more/smaller clusters.
#'   Lower values yield fewer/larger clusters. Optimal value depends on data scale.
#' @param .small.size Integer threshold for straggler clusters. Clusters with fewer
#'   cells are absorbed into neighbors. Default: 0.5% of total nodes (floor(n/200)).
#' @param .verbose Logical, print progress and cluster counts. Default TRUE.
#' @param .seed Integer for reproducibility (affects k-means and Leiden). Default 123.
#' 
#' @return Factor vector of cluster assignments with levels ordered by appearance.
#'   Labels are formatted as "cluster.01", "cluster.02", etc.
#' 
#' @examples
#' \dontrun{
#' # After building similarity matrix
#' clusters <- leiden.cluster(
#'   .tdr.obj = lm.obj,
#'   .sim.matrix = jaccard_sim,
#'   .resolution.parameter = 0.001,
#'   .small.size = 10  # Merge clusters < 10 cells
#' )
#' table(clusters)
#' }
leiden.cluster <-
  function(.tdr.obj = NULL,
           .sim.matrix,
           .resolution.parameter,
           .small.size = (nrow(x = .sim.matrix) / 200) |> floor(),
           .verbose = TRUE,
           .seed = 123) {
    
    # Determine initialization embedding: LE > PCA > computed from similarity
    .init.embed <- 
      if(!is.null(x = .tdr.obj)){
        if(is.null(.tdr.obj@landmark.embed$le$coord) || isTRUE(is.na(x = .tdr.obj@landmark.embed$le$coord[1,1]))) {
          # Fallback to PCA if LE not available
          .tdr.obj@landmark.embed$pca$coord[,1:min(3, ncol(x = .tdr.obj@landmark.embed$pca$coord)), drop = FALSE]
        } else {
          .tdr.obj@landmark.embed$le$coord
        }
      } else {
        # No lm.obj: compute 3D embedding from similarity matrix
        irlba::irlba(A = .sim.matrix,
                     nv = 3,
                     nu = 3)$u
      }  
    
    set.seed(seed = .seed)
    
    if(dim(x = .sim.matrix) |>
       (\(x)
        !identical(x = x[1], 
                   y = x[2])
       )()){
      stop("Similarity matrix must be square.\n",
           "Dimensions: ", nrow(.sim.matrix), " x ", ncol(.sim.matrix))
    }
    
    # Remove self-loops (diagonal)
    Matrix::diag(x = .sim.matrix) <-
      0
    
    # Convert to igraph object
    g <-
      igraph::graph_from_adjacency_matrix(
        adjmatrix = .sim.matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE,
        add.colnames = NULL,
        add.rownames = NA)
    
    # Initialize with k-means to avoid starting Leiden from singletons
    # Up to 25 clusters or half the data size, whichever is smaller
    kres <-
      stats::kmeans(x = .init.embed,
                    centers = min(25,
                                  (nrow(x = .init.embed) / 2) |>
                                    ceiling()),
                    nstart = 10,
                    iter.max = 100)
    
    # Handle k-means convergence failure (ifault=4: no convergence in iter.max)
    # Switch to MacQueen algorithm which is more robust but slower
    # See: https://stackoverflow.com/a/30055776
    if(kres$ifault == 4){
      kres <-
        stats::kmeans(x = .init.embed,
                      centers = kres$centers,
                      nstart = 10,
                      algorithm = "MacQueen")
    }
    
    m0 <-
      as.integer(x = kres$cluster)
    
    # Run Leiden clustering with CPM objective (Constant Potts Model)
    # beta=0.01: small randomness for reproducibility
    # n_iterations=10: refinement passes
    ids <-
      igraph::cluster_leiden(graph = g,
                             objective_function = "CPM",
                             weights = igraph::E(graph = g)$weight,
                             resolution = .resolution.parameter,
                             beta = 0.01,
                             initial_membership = m0,
                             n_iterations = 10,
                             vertex_weights = NULL) |>
      igraph::membership() |>
      as.integer() |>
      # Format cluster IDs with zero-padding
      (\(x)
       formatC(x = x,
               width = nchar(x = x) |>
                 max(),
               format = "d",
               flag = "0")
      )() |>
      (\(x)
       paste0("cluster.",
              x)
      )()
    
    if(isTRUE(x = .verbose)){
      message(paste(
        length(x = unique(x = ids)),
        "clusters identified."
      ))}
    
    # Straggler absorption: reassign small clusters to better-connected neighbors
    if (isTRUE(.verbose)) {
      message("checking for stragglers.")
    }
    
    # Identify clusters smaller than threshold
    tbl <- table(ids)
    stragglers <- names(tbl)[tbl < .small.size]
    keepers    <- setdiff(names(tbl), stragglers)
    
    if (length(stragglers) > 0 && length(keepers) > 0) {
      
      # Ensure similarity matrix has row/col names for indexing
      if (is.null(rownames(.sim.matrix))) {
        # If you have real cell IDs, use those instead of seq_len
        rownames(.sim.matrix) <- colnames(.sim.matrix) <- seq_len(nrow(.sim.matrix))
      }
      
      # For each straggler cluster, find the keeper with highest mean connectivity
      # For each straggler cluster, find the keeper with highest mean connectivity
      for (i in stragglers) {
        i.cells <- which(ids == i)
        # Build connectivity vector: mean edge weight to each keeper cluster
        connectivity <- numeric(length(keepers))
        names(connectivity) <- keepers
        
        for (j in keepers) {
          j.cells <- which(ids == j)
          # Extract submatrix of straggler-to-keeper edges
          sub.sim.matrix <- .sim.matrix[i.cells, j.cells, drop = FALSE]
          
          if (inherits(sub.sim.matrix, "sparseMatrix")) {
            # For sparse: mean of non-zero entries scaled by matrix size
            connectivity[j] <- sum(sub.sim.matrix@x) / (nrow(sub.sim.matrix) * ncol(sub.sim.matrix))
          } else {
            connectivity[j] <- mean(sub.sim.matrix)
          }
        }
        
        # Assign to most connected keeper (break ties randomly)
        m  <- max(connectivity, na.rm = TRUE)
        mi <- which(connectivity == m)
        closest.cluster <- names(connectivity)[ if (length(mi) == 1) mi else sample(mi, 1) ]
        
        ids[i.cells] <- closest.cluster
      }
      
      if (isTRUE(.verbose)) {
        message(length(stragglers), " stragglers absorbed; ",
                length(unique(ids)), " final clusters.")
      }
    }
    
    # Relabel clusters sequentially starting from 1 with formatted names
    uniq <- 
      unique(x = ids)
    map  <- 
      seq_along(along.with = uniq) |>
      stats::setNames(nm = uniq)
    lab  <- 
      paste0("cluster.", "%0", nchar(x = length(x = uniq)), "d") |>
      sprintf(map[ids])
    ids  <- 
      factor(x = lab, 
             levels = unique(x = lab))
    
    if(isTRUE(x = .verbose)){
      table(ids) |>
        print()
    }
    
    return(ids)
  }


#' Sparse matrix representation of nearest neighbors
#' 
#' Converts a nearest neighbor index matrix into a sparse adjacency matrix. The resulting matrix 
#' is non-symmetric (i can be neighbor of j without j being neighbor of i). Note that the matrix 
#' orientation places neighbors in rows and source cells in columns.
#' 
#' @param .nn.idx A matrix where each row represents a cell and each column represents a nearest 
#'   neighbor index. Element (i,k) is the index of the k-th nearest neighbor of cell i.
#'   
#' @return A non-symmetric sparse matrix (dgCMatrix) of dimensions n×n where n is the number of 
#'   cells. Entry (j,i) = 1 if j is among i's nearest neighbors, 0 otherwise. That is, row j 
#'   contains 1s in columns corresponding to cells that have j as a neighbor.
#'
#' @keywords internal
#'
get.adj.matrix <-
  function(.nn.idx){
    # Convert NN index matrix to sparse: each row spreads to its NN columns
    Matrix::sparseMatrix(i = as.vector(x = .nn.idx),  # Target cells (NN indices)
                         j = nrow(x = .nn.idx) |>      # Source cells (repeated k times)
                           seq_len() |>
                           rep(times = ncol(x = .nn.idx)),
                         x = 1,                         # Binary adjacency
                         dims = nrow(x = .nn.idx) |>
                           rep(times = 2),              # Square matrix
                         repr = "C",                    # Compressed column format
                         dimnames = list(rownames(x = .nn.idx),
                                         rownames(x = .nn.idx)))
  }

#' Shared nearest neighbors via fast Jaccard index calculation
#'
#' Computes a shared nearest neighbor (SNN) graph by calculating Jaccard similarity between 
#' cells based on their k-nearest neighbors. The Jaccard index J(A,B) = |A∩B| / |A∪B| measures 
#' the overlap of neighbor sets. This implementation uses sparse matrix operations for efficiency.
#'
#' @param .adj.matrix A non-symmetric sparse adjacency matrix from \code{get.adj.matrix}, where 
#'   entry (j,i)=1 indicates j is among i's nearest neighbors (neighbors in rows, cells in columns).
#' @param .prune A numeric threshold (default 1/15 ≈ 0.067) for pruning weak edges. Jaccard 
#'   values below this are set to zero for sparsity.
#'   
#' @return A symmetric sparse matrix of Jaccard indices representing the SNN graph. Higher values 
#'   indicate greater neighbor overlap. Entry (i,j) gives the Jaccard similarity between cell i 
#'   and cell j based on their shared k-nearest neighbors. Matrix is pruned to remove weak 
#'   connections and forced symmetric.
#'   
#' @keywords internal
#'
fast.jaccard.r <-
  function(
    .adj.matrix,
    .prune = 1/15
  ){
    
    # Crossprod gives intersection counts: adj^T × adj
    adj.cp <-
      Matrix::crossprod(x = .adj.matrix)
    
    # Convert intersection to Jaccard: J = intersect / union
    # Union = 2k - intersect (both have k neighbors, minus overlap)
    adj.cp@x <-
      1 / (((sum(x = .adj.matrix[,1]) * 2) / adj.cp@x) - 1)
    
    # Prune weak edges, symmetrize, convert to general format
    adj.cp <-
      Matrix::drop0(x = adj.cp,
                    tol = .prune) |>
      Matrix::forceSymmetric() |>
      methods::as(Class = "generalMatrix")
    
    return(adj.cp)
    
  }

#' Graph embedding of landmarks
#'
#' Builds k-NN graph and trains UMAP model for landmark cells. This creates the graph structure 
#' and transformation model required for generating fuzzy density matrices in \code{get.map}. 
#' Clustering is performed as a secondary step for visualization and interpretation, but is not 
#' the primary objective.
#'
#' @param .k Integer number of nearest neighbors for graph construction (default 20). Higher values 
#'   create more connected graphs and smoother UMAP embeddings. Used for k-NN graph construction 
#'   and passed to \code{uwot::umap} as \code{n_neighbors}.
#' @param .scale Logical indicating whether to scale features before UMAP. Defaults to FALSE
#'   when Harmony batch correction is active (any assay type) or for RNA assays (PCA already
#'   scaled). Defaults to TRUE for cytometry without Harmony (raw marker input to UMAP).
#' @param .verbose Logical for progress messages (default TRUE).
#' @param .seed Integer seed for reproducibility of UMAP and clustering (default 123).
#' @param .cl.method Character specifying clustering method: "snn" (shared nearest neighbors via 
#'   Jaccard similarity) or "fgraph" (fuzzy graph from UMAP). Default "snn".
#' @param .cl.resolution.parameter Numeric controlling cluster granularity (default 0.8). Higher 
#'   values produce more fine-grained clusters.
#' @param .small.size Integer threshold for straggler clusters (default 0.5% of landmarks). Clusters 
#'   smaller than this are absorbed into nearest large cluster.
#'   
#' @return Updated \code{.tdr.obj} with \code{$graph} component containing:
#'   \itemize{
#'     \item \code{uwot}: UMAP model with multiple components:
#'       \itemize{
#'         \item \code{$embedding}: 2D UMAP coordinates for visualization
#'         \item \code{$nn$euclidean$idx}: k-NN indices matrix (landmarks × k neighbors)
#'         \item \code{$nn$euclidean$dist}: k-NN distance matrix
#'         \item \code{$fgraph}: Fuzzy graph (landmark-landmark probabilistic edges) for optional clustering
#'         \item \code{$model}: Trained UMAP model for projecting query cells in \code{get.map}
#'       }
#'     \item \code{adj.matrix}: Sparse k-NN adjacency matrix (non-symmetric).
#'     \item \code{snn}: Shared nearest neighbor graph via Jaccard similarity.
#'     \item \code{LE}: Laplacian Eigenmap components computed from symmetrized k-NN graph, including 
#'       \code{$W.sym} (symmetrized adjacency), \code{$L} (Laplacian), \code{$Disqrt} (degree matrix), 
#'       \code{$vectors}, \code{$values}, \code{$converged}, \code{$nontriv} (non-trivial eigenvalue indices), 
#'       \code{$elbow} (selected dimensionality), \code{$keep} (retained eigenvectors), and \code{$embed} 
#'       (final normalized embedding). Used for clustering initialization only.
#'     \item \code{clustering}: List containing \code{$ids} (factor of cluster assignments), \code{$median.exprs} 
#'       (matrix of mean expression per cluster), and \code{$pheatmap} (heatmap object).
#'   }
#'   
#' @details
#' The function orchestrates graph-based analysis in the following order:
#' 
#' \strong{1. k-NN Graph Construction:}
#' Builds k-nearest neighbor graph from PCA/Harmony embeddings and converts it to multiple 
#' representations for different downstream applications:
#' \itemize{
#'   \item \code{adj.matrix}: Non-symmetric sparse adjacency matrix (used to compute SNN graph)
#'   \item \code{snn}: Shared nearest neighbor graph via Jaccard similarity (used for clustering)
#'   \item \code{W.sym}: Symmetrized binary adjacency (used for Laplacian Eigenmap spectral analysis)
#' }
#' 
#' \strong{2. UMAP Model Training:}
#' Trains UMAP transformation on landmark cells with \code{ret_model=TRUE} to enable projection of 
#' query cells later. The UMAP model includes the fuzzy graph (landmark-landmark probabilistic edges) which 
#' is used for optional clustering and, after projection in \code{get.map}, generates cell-landmark 
#' edge weights essential for computing fuzzy density matrices. The UMAP model is the primary output 
#' of \code{get.graph}.
#' 
#' \strong{3. Laplacian Eigenmap (LE):}
#' Computes spectral embedding by solving the generalized eigenvalue problem of the graph Laplacian. 
#' Dimensionality is automatically selected via elbow detection on eigenvalue spectrum. The LE 
#' embedding is currently used only for warm-start initialization in clustering.
#' 
#' \strong{4. Clustering (Secondary):}
#' Applies Leiden algorithm to identify communities for visualization and interpretation. Uses SNN 
#' graph (Jaccard similarity) or UMAP fuzzy graph depending on \code{.cl.method}. Small "straggler" 
#' clusters are absorbed into neighbors. While useful for exploration, clustering is not required 
#' for downstream statistical analysis, unless traditional analysis is requested later in \code{get.lm}.
#' 
#' @seealso \code{\link{setup.tdr.obj}}, \code{\link{get.landmarks}}, \code{\link{lm.cluster}}
#' 
#' @examples
#' \dontrun{
#' # Typical workflow after landmark selection
#' lm.cells <- setup.tdr.obj(.cells = .cells, 
#'                           .meta = .meta,
#'                           .assay.type = "RNA") |>
#'   get.landmarks(.nHVG = 500, .nPC = 3) |>
#'   get.graph(.k = 10)
#' 
#' # Higher resolution for finer clusters
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph(.cl.resolution.parameter = 200)
#' 
#' # Use fuzzy graph for clustering instead of SNN
#' lm.cells <- get.graph(lm.cells, .cl.method = "fgraph")
#' }
#' 
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object.
#' @param ... Additional arguments passed to methods.
#' @export
get.graph <- function(x, ...) UseMethod("get.graph")

#' @rdname get.graph
#' @export
get.graph.TDRObj <-
  function(x,
           .k = 20,
           .scale = if(!is.null(x@integration$harmony.obj) || x@config$assay.type == "RNA") FALSE else TRUE,
           .verbose = TRUE,
           .seed = 123,
           .cl.method = "snn",
           .cl.resolution.parameter = 0.8,
           .small.size = floor(x = nrow(x = x@assay$expr) / 200),
           ...){
    .tdr.obj <- x
    
    # Validate clustering method
    .cl.method <-
      match.arg(arg = .cl.method,
                choices = c("fgraph","snn"))
    
    # Compute UMAP embedding with k-NN graph
    set.seed(seed = .seed)
    .tdr.obj@integration$umap.model <-
      uwot::umap(X = if(!is.null(x = .tdr.obj@integration$harmony.obj) || .tdr.obj@config$assay.type == "RNA") .tdr.obj@landmark.embed$pca$coord else .tdr.obj@assay$expr,
                 n_neighbors = .k,
                 n_components = 2,           # 2D for visualization
                 n_epochs = 500,             # Training iterations
                 scale = .scale,
                 pca = NULL,                 # Already have PCA/LE
                 verbose = isTRUE(x = .verbose),
                 ret_model = TRUE,           # Keep model for query projection
                 batch = TRUE,
                 seed = .seed,
                 n_threads = .tdr.obj@config$n.threads,
                 fast_sgd = FALSE,
                 n_sgd_threads = .tdr.obj@config$n.threads,
                 ret_extra = c("fgraph",     # Fuzzy graph (for optional clustering)
                               "nn"))         # Nearest neighbors
    
    # Extract embedding and fuzzy graph from the uwot result into their new slots
    .tdr.obj@landmark.embed$umap <- list(coord = .tdr.obj@integration$umap.model$embedding)
    .tdr.obj@graphs$fgraph <- .tdr.obj@integration$umap.model$fgraph
    
    colnames(x = .tdr.obj@landmark.embed$umap$coord) <-
      c("umap.1",
        "umap.2")
    
    # Convert k-NN indices to sparse adjacency matrix
    .tdr.obj@graphs$adj.matrix <-
      get.adj.matrix(.nn.idx = .tdr.obj@integration$umap.model$nn$euclidean$idx)
    
    # Compute SNN graph via Jaccard similarity of neighbor overlaps
    .tdr.obj@graphs$snn <-
      fast.jaccard.r(.adj.matrix = .tdr.obj@graphs$adj.matrix,
                     .prune = 1/15) |>
      (\(x)
       `dimnames<-`(x = x,
                    value = list(rownames(x = .tdr.obj@assay$expr),
                                 rownames(x = .tdr.obj@assay$expr)))
      )()
    
    if(isTRUE(x = .verbose)){
      message("getting Laplacian Eigenmap")
    }
    
    # Symmetrize adjacency: W[i,j] = 1 if i is NN of j OR j is NN of i
    # See uwot implementation: https://github.com/jlmelville/uwot/blob/f9e576e97d9df44d48be2cc559412282838dc4a5/R/init.R
    .tdr.obj@landmark.embed$le$W.sym <-
      (((.tdr.obj@graphs$adj.matrix > 0) |
          (Matrix::t(x = .tdr.obj@graphs$adj.matrix) > 0)) * 1) |>
      Matrix::Matrix(sparse = TRUE)
    
    # Target dimensions: match PCA dimensionality
    target_k <-
      ncol(x = .tdr.obj@landmark.embed$pca$coord)
    nv <-
      min(target_k, 
          nrow(x = .tdr.obj@landmark.embed$pca$coord) - 1L) # Eigenspectrum limited to n-1
    stopifnot(nv >= 2L)
    
    # Form modified graph Laplacian: L = D^(-1/2) (D - W) D^(-1/2)
    .tdr.obj@landmark.embed$le <- 
      c(.tdr.obj@landmark.embed$le,
        form_modified_laplacian(A = .tdr.obj@landmark.embed$le$W.sym,
                                ret_d = TRUE))
    
    # Compute truncated SVD of Laplacian for nv components
    .tdr.obj@landmark.embed$le <-
      c(.tdr.obj@landmark.embed$le,
        irlba_spectral_tsvd(L = .tdr.obj@landmark.embed$le$L, 
                            n = nv))
    
    # Verify convergence and valid eigenvectors
    ok <-
      all(c("vectors", "values", "converged") %in% 
            names(x = .tdr.obj@landmark.embed$le)) &&
      is.matrix(x = .tdr.obj@landmark.embed$le$vectors) &&
      (ncol(x = .tdr.obj@landmark.embed$le$vectors) >= nv) &&
      isTRUE(x = .tdr.obj@landmark.embed$le$converged)
    
    if(!ok) {
      
      message("Laplacian Eigenmap failed to converge.")
      
      .tdr.obj@landmark.embed$le$coord <- 
        matrix(data = NA_real_,
               nrow = 1,
               ncol = 1)
      
    } else {
      
      # Filter out trivial eigenvalues (near-zero, correspond to disconnected components)
      tol <- 1e-6
      .tdr.obj@landmark.embed$le$nontriv <-
        which(x = .tdr.obj@landmark.embed$le$values > tol)
      
      if(length(x = .tdr.obj@landmark.embed$le$nontriv) == 0L) {
        
        message("No non-trivial Laplacian eigenvalues found (graph may be empty or singular).")
        
        .tdr.obj@landmark.embed$le$coord <- 
          matrix(data = NA_real_,
                 nrow = 1,
                 ncol = 1)
        
      } else {
        
        # Determine embedding dimensionality via elbow detection
        if(length(x = .tdr.obj@landmark.embed$le$nontriv) < 4L){
          
          k <- length(x = .tdr.obj@landmark.embed$le$nontriv)  # Too few for elbow, use all
          
        } else {
          
          .tdr.obj@landmark.embed$le$elbow <-
            elbow.sec.deriv(x = .tdr.obj@landmark.embed$le$values[.tdr.obj@landmark.embed$le$nontriv],
                            smooth = TRUE,
                            df = NULL,
                            sort.order = "asc")$index
          
          # Use elbow dimensions, bounded by target and available eigenvectors
          k <- 
            min(.tdr.obj@landmark.embed$le$elbow, 
                target_k, 
                length(x = .tdr.obj@landmark.embed$le$nontriv)) |>
            max(2L)
          
        }
        
        # Extract first k non-trivial eigenvectors
        .tdr.obj@landmark.embed$le$keep <-
          .tdr.obj@landmark.embed$le$nontriv[seq_len(length.out = k)]
        
        # Normalize eigenvectors: multiply by D^(-1/2) and L2-normalize columns
        .tdr.obj@landmark.embed$le$coord <-
          (.tdr.obj@landmark.embed$le$Disqrt * .tdr.obj@landmark.embed$le$vectors[,.tdr.obj@landmark.embed$le$keep, drop = FALSE]) |> 
          (\(x)
           sweep(x = x, 
                 MARGIN = 2, 
                 STATS = Matrix::colSums(x * x) |>
                   sqrt(), 
                 FUN = `/`)
          )() |>
          (\(x)
           `dimnames<-`(x = x,
                        value = list(rownames(x = .tdr.obj@landmark.embed$pca$coord),
                                     paste0("LE", 1:ncol(x = x))))
          )()
      }
    }
    
    if(isTRUE(x = .verbose)){
      message("clustering")
    }
    
    # Apply Leiden clustering with straggler absorption
    .tdr.obj <-
      lm.cluster(.tdr.obj,
                 .cl.method = .cl.method,
                 .cl.resolution.parameter = .cl.resolution.parameter,
                 .seed = .seed,
                 .verbose = .verbose,
                 .small.size = .small.size)
    
    return(.tdr.obj)
    
  }

#' Mapping cells to landmarks
#'
#' Projects all cells onto the landmark graph to compute fuzzy graph edge weights between cells and landmarks.
#' In addition, transfers cluster/cell type labels from landmarks to all cells.
#' 
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object with \code{$graph} component populated by \code{get.graph}.
#' @param .source The raw data object for non-file backends. \code{NULL} (default) for 
#'   the files backend; otherwise a Seurat, SingleCellExperiment, or anndataR AnnData object. 
#'   Used by \code{.get_sample_matrix()} to retrieve per-sample expression matrices.
#' @param .ref.obj Optional Symphony reference object for cell type annotation. Must have 
#'   \code{Z_corr} field (harmony-corrected embeddings) and metadata with cell type labels. 
#'   Only compatible with RNA assays. Replaces any existing \code{$graph$celltyping}.
#' @param .celltype.col.name Column name in \code{.ref.obj$meta_data} containing cell type labels 
#'   (default "cell_type"). Only relevant when \code{.ref.obj} is provided.
#' @param .verbose Logical for progress messages (default TRUE).
#' @param .seed Integer seed for reproducibility (default 123).
#' @param .label.confidence Numeric scalar in \code{[0,1]} controlling the minimum posterior confidence required
#'   to assign a cell to a landmark‑derived cluster/celltype label.
#' @param .cache.on.disk Logical (default TRUE). When \code{TRUE}, four large per-sample
#'   slots (\code{clustering$ids}, \code{celltyping$ids}, \code{nearest.lm},
#'   \code{fuzzy.graphs}) are serialized to disk as uncompressed RDS files and stored
#'   in \code{@cellmap} as attributed path strings.  Downstream accessors
#'   (e.g.\ in \code{get.pbDE}, \code{goi.summary}) read them back lazily on a
#'   per-sample basis.  Set to \code{FALSE} to keep everything in memory.
#'
#'   Cache files are stored under the system temporary directory
#'   (\code{tempdir()}) and are automatically removed when the R session
#'   ends via a registered finalizer.  This means the cache is
#'   \strong{ephemeral} and never persists across R sessions.  There are
#'   no implications for reproducibility since the cache only stores
#'   intermediate results that are recomputed deterministically.
#'   
#' @return Updated \code{.tdr.obj} with \code{$map} component containing:
#'   \itemize{
#'     \item \code{fdens}: Matrix of fuzzy graph densities (landmarks × samples). Each entry is 
#'       the sum of cell-landmark edge weights for that sample, normalized by sample size.
#'     \item \code{Y}: Matrix of log2-transformed densities (landmarks × samples): 
#'       \code{log2(fdens + 0.5)}. Used by \code{get.lm()} for linear modeling and 
#'       \code{get.embedding()} for unsupervised embeddings.
#'     \item \code{clustering$ids}: List of named character vectors (one per sample) with cluster 
#'       assignments for all cells.
#'     \item \code{clustering$cell.count}: Matrix (samples × clusters) of cell counts per cluster 
#'       per sample. Used for "traditional" compositional statistics.
#'     \item \code{clustering$cell.perc}: Matrix (samples × clusters) of percentage of cells per 
#'       cluster per sample. Used for "traditional" compositional statistics.
#'     \item \code{celltyping$ids}: List of named character vectors (one per sample) with cell type 
#'       assignments (only if celltyping available or \code{.ref.obj} provided).
#'     \item \code{celltyping$cell.count}: Matrix (samples × cell types) of cell counts per cell type 
#'       per sample. Used for "traditional" compositional statistics.
#'     \item \code{celltyping$cell.perc}: Matrix (samples × cell types) of percentage of cells per 
#'       cell type per sample. Used for "traditional" compositional statistics.
#'     \item \code{nearest.lm}: List of matrices (one per sample) with nearest landmark indices 
#'       for all cells from UMAP transform.
#'   }
#'
#'   When \code{.cache.on.disk = TRUE}, the four cell-level slots above are stored as
#'   attributed path strings (with \code{schema_v} and \code{bytes} attributes) in
#'   \code{@cellmap} rather than in-memory objects.  The cache root is stored in
#'   \code{@config$.cache.root}.  Use \code{tdr_cache_cleanup()} to remove cached files.
#'   If \code{.ref.obj} provided, also updates \code{$graph$celltyping} with list containing \code{$ids} 
#'   (factor of cell type assignments for landmarks), \code{$median.exprs} (matrix of mean expression 
#'   per cell type), and \code{$pheatmap} (heatmap object).
#'   
#' @details
#' \strong{Workflow Overview:}
#' 
#' For each sample, the function:
#' 1. Loads expression data and normalizes (size factors for RNA, marker subset for cytometry)
#' 2. Projects to PCA/Harmony space (matching landmark processing)
#' 3. Uses landmark UMAP model to compute fuzzy graph (cell-landmark edge weights) and find nearest landmarks
#' 4. Assigns clusters/cell types by confidence-thresholded voting
#' 5. Aggregates fuzzy graph edge weights into landmark densities per sample
#'
#' \strong{Label transfer confidence model:}
#' - Without `.ref.obj`: label confidence is the normalized fuzzy-mass ratio,
#'   \eqn{\mathrm{conf}(c,\ell)=\sum_{m\in\ell}w_{c,m}/\sum_m w_{c,m}},
#'   where \eqn{w_{c,m}} are UMAP-derived cell-landmark connection strengths.
#' - With `.ref.obj`: label confidence is kNN voting frequency in reference space,
#'   \eqn{\mathrm{conf}(c,\ell)=N_{c,\ell}/k} with \eqn{k=10} nearest neighbors.
#' - In both modes, a label is accepted only if
#'   \eqn{\mathrm{conf}(c,\ell) \ge {.label.confidence}}; otherwise
#'   the cell is labeled `"..low.confidence.."`.
#' 
#' \strong{Reference-Based Cell Typing:}
#' 
#' When \code{.ref.obj} is provided:
#' - Expression is mapped to reference via Symphony
#' - Cell types assigned by kNN voting (k = 10) in reference embedding
#' - Landmark cell types updated and used for visualization/statistics
#' - Replaces any existing manual cell type annotations
#' 
#' \strong{Fuzzy Graph Densities:}
#' 
#' The \code{fdens} matrix quantifies how strongly each landmark is connected to cells in each 
#' sample. High values indicate the landmark's neighborhood is enriched in that sample. This 
#' forms the basis for differential density testing in \code{get.lm}.
#' 
#' @seealso \code{\link{get.graph}}, \code{\link{get.lm}}, \code{\link{celltyping}}
#' 
#' @examples
#' \dontrun{
#' # Complete workflow with mapping
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500) |>
#'   get.graph() |>
#'   get.map()
#' 
#' # Use Symphony reference for cell typing (RNA data only)
#' ref <- readRDS("pbmc_reference.rds")
#' lm.cells <- get.map(lm.cells, 
#'                     .ref.obj = ref, 
#'                     .celltype.col.name = "cell_type")
#' }
#' 
#' @param ... Additional arguments passed to methods.
#' @export
get.map <- function(x, ...) UseMethod("get.map")

#' @rdname get.map
#' @export
get.map.TDRObj <-
  function(x,
           .source = NULL,
           .ref.obj = NULL,
           .celltype.col.name = "cell_type",
           .verbose = TRUE,
           .seed = 123,
           .label.confidence = 0.5,
           .cache.on.disk = TRUE,
           ...){
    .tdr.obj <- x

    .tdr_validate_label_confidence(.label.confidence)
    
    # R CMD check appeasement for non-standard evaluation in dplyr and collapse
    ref.idx <- N <- cell.pop <- id <- value <- ri <- i <- j <- landmark <- cell <- label <- x <- confidence <-
      NULL
    
    # Reject removed argument .cl.ct.to.ign if passed via ...
    if (".cl.ct.to.ign" %in% names(list(...))) {
      stop(
        "The '.cl.ct.to.ign' argument has been removed from get.map().\n",
        "Cell type / cluster exclusion must now be performed upstream ",
        "(e.g., by subsetting cells before running the pipeline).\n",
        "See ?get.map for the updated workflow."
      )
    }
    
    # Persist label confidence for downstream refresh (e.g. late celltyping)
    .tdr.obj@config$label.confidence <- .label.confidence
    
    # Validate and setup Symphony reference if provided
    if(!is.null(x = .ref.obj)){
      
      if(.tdr.obj@config$assay.type != "RNA"){
        stop("Cell typing with reference object is only supported for assay type RNA.\nCytometry data requires manual annotation via celltyping().")
      }
      
      if(inherits(x = .ref.obj,
                  what = "list")){
        
        if(is.null(x = .ref.obj$Z_corr)){
          stop(".ref.obj must be a Symphony object with Z_corr field populated.\nEnsure reference was built with symphony::buildReference().")
          
        }
        
      } else {
        
        stop(".ref.obj must be a Symphony object (list) with Z_corr field.\nProvided: ", class(.ref.obj)[1])
        
      }
      
      if(!(.celltype.col.name %in% colnames(x = .ref.obj$meta_data))){
        
        stop("Column name '", .celltype.col.name,
             "' not found in .ref.obj$meta_data.\nAvailable columns: ",
             paste(colnames(.ref.obj$meta_data), collapse = ", "))
        
      }
      
      warning("Using Symphony reference object for cell typing. Removing previous celltyping from .tdr.obj.")
      
      if(isTRUE(x = .verbose)){
        message("building Annoy index for Symphony reference object")
      }
      
      .tdr.obj@integration$symphony.obj <-
        .ref.obj
      
      # Build k-NN index on reference embeddings for fast query lookup
      set.seed(seed = .seed)
      .tdr.obj@integration$symphony.obj$ref.knn.idx <-
        Matrix::t(x = .tdr.obj@integration$symphony.obj$Z_corr) |>
        annoy_build(metric = "euclidean", 
                    n_trees = 50,
                    verbose = .verbose)
      
      .tdr.obj@integration$symphony.obj$celltype.col.name <-
        .celltype.col.name
      
      .tdr.obj@landmark.annot$celltyping <- NULL
      
    }
    
    # Determine cell type levels for fdens filtering (all populations kept)
    .ct.to.keep <-
      if(is.null(x = .ref.obj)){
        levels(x = .tdr.obj@landmark.annot$celltyping$ids)
      } else {
        unique(x = .ref.obj$meta_data[[.celltype.col.name]])
      }
    
    set.seed(seed = .seed)
    
    # ── Set up on-disk cache directory ──
    backend <- .tdr.obj@config$backend
    if (is.null(backend)) backend <- "files"

    .cache.meta <- NULL
    if (isTRUE(x = .cache.on.disk)) {
      .run_key <- .tdr_make_run_key()
      .run_cache_dir <- file.path(.tdr_cache_root(), .run_key)
      dir.create(.run_cache_dir, recursive = TRUE, showWarnings = FALSE)
      .tdr_cache_sweep_orphans(.run_cache_dir)
      .tdr_cache_register_cleanup(.run_cache_dir)
      
      if (isTRUE(x = .verbose)) {
        message("-> On-disk caching enabled: ", .run_cache_dir)
      }
    }
    
    # Process each sample: load data, normalize, project to UMAP, assign labels
    if(isTRUE(x = .verbose)){
      message("-> Mapping ", length(.tdr.obj@cells), " samples to landmark space...")
      .map_start <- Sys.time()
    }
    
    res <-
      seq_along(along.with = .tdr.obj@cells) |>
      stats::setNames(nm = names(x = .tdr.obj@cells)) |> #_[1:2] |>
      lapply(FUN = function(cells.idx){
        
        set.seed(seed = .seed)
        
        mat.exprs <-
          .get_sample_matrix(.source, .tdr.obj, cells.idx)
        
        if(.tdr.obj@config$assay.type == "RNA"){
          
          mat.exprs <-
            methods::as(object = mat.exprs, 
              Class = "dgCMatrix") |>
            Matrix::t() |>
            (\(x)
             # size factor normalization to landmark library-size scale
             x / (Matrix::rowSums(x = x) /
                    mean(x = Matrix::rowSums(x = .tdr.obj@assay$raw)))
            )()
          
          mat.exprs@x <-
            log2(x = mat.exprs@x + 1)
          
          temp.cell.names <-
            paste0(names(x = .tdr.obj@cells)[cells.idx],
                   "_",
                   rownames(x = mat.exprs))
          
        } else {
          
          mat.exprs <-
            mat.exprs[,.tdr.obj@config$markers]
          
          temp.cell.names <-
            paste0(names(x = .tdr.obj@cells)[cells.idx],
                   "_",
                   rownames(x = mat.exprs))
          
        }
        
        cells.of.interest <-
          !(temp.cell.names %in% rownames(x = .tdr.obj@assay$expr))
        
        if(!is.null(x = .tdr.obj@integration$harmony.obj)){
          
          # Map query
          mat.exprs <-
            symphony::mapQuery(exp_query = Matrix::t(x = mat.exprs[,.tdr.obj@landmark.embed$pca$HVG]),
                               metadata_query = matrix(data = 1,
                                                       nrow = nrow(x = mat.exprs),
                                                       ncol = 2),
                               ref_obj = .tdr.obj@integration$harmony.obj,
                               vars = NULL,
                               verbose = .verbose,
                               do_normalize = FALSE,
                               do_umap = FALSE)$Z |>
            Matrix::t() |> 
            (\(x)
             `rownames<-`(x = x,
                          value = rownames(x = mat.exprs))
            )()
          
        } else if(.tdr.obj@config$assay.type == "RNA"){
          
          mat.exprs <-
            (((Matrix::t(x = mat.exprs[,.tdr.obj@landmark.embed$pca$HVG]) - .tdr.obj@landmark.embed$pca$center) /
                .tdr.obj@landmark.embed$pca$scale) |>
               Matrix::t()) %*%
            .tdr.obj@landmark.embed$pca$rotation
          
        }
        
        if(inherits(x = mat.exprs,
                    what = "dgeMatrix")){
          mat.exprs <-
            matrix(data = mat.exprs@x,
                   nrow = nrow(x = mat.exprs),
                   ncol = ncol(x = mat.exprs),
                   dimnames = dimnames(x = mat.exprs))
        }
        
        res2 <-
          uwot::umap_transform(
            X = mat.exprs,
            model = .tdr.obj@integration$umap.model,
            init_weighted = TRUE,
            search_k = NULL,
            tmpdir = tempdir(),
            n_epochs = 0,
            n_threads = .tdr.obj@config$n.threads,
            n_sgd_threads = .tdr.obj@config$n.threads,
            grain_size = 1,
            verbose = isTRUE(x = .verbose),
            init = "average",
            batch = TRUE,
            learning_rate = NULL,
            opt_args = NULL,
            epoch_callback = NULL,
            ret_extra = c("fgraph","nn"),
            seed = .seed
          )
        
        if(isTRUE(x = .verbose)){
          message("assigning cells to clusters")
        }
        
        res2$cell.clustering <-
          .tdr_transfer_labels(
            .method = "fuzzy",
            .label.confidence = .label.confidence,
            .n.cells = nrow(x = mat.exprs),
            .cell.names = rownames(x = mat.exprs),
            .fgraph = res2$fgraph,
            .landmark.labels = .tdr.obj@landmark.annot$clustering$ids
          )
        
        if(is.null(x = .tdr.obj@integration$symphony.obj)){
          
          if(!is.null(x = .tdr.obj@landmark.annot$celltyping$ids)){
            
            if(isTRUE(x = .verbose)){
              
              message("assigning cells to celltypes")
              
            }
            
            res2$cell.celltyping <-
              .tdr_transfer_labels(
                .method = "fuzzy",
                .label.confidence = .label.confidence,
                .n.cells = nrow(x = mat.exprs),
                .cell.names = rownames(x = mat.exprs),
                .fgraph = res2$fgraph,
                .landmark.labels = .tdr.obj@landmark.annot$celltyping$ids
              )
            
          }
          
        } else {
          
          if(isTRUE(x = .verbose)){
            message("assigning cells to celltypes using symphony reference obj.")
          }
          
          raw.exprs <-
            .get_sample_matrix(.source, .tdr.obj, cells.idx) |>
            methods::as(Class = "dgCMatrix")
          
          # Map query
          query <-
            symphony::mapQuery(exp_query = raw.exprs,
                               metadata_query = matrix(data = 1,
                                                       nrow = ncol(x = raw.exprs),
                                                       ncol = 2),
                               ref_obj = .tdr.obj@integration$symphony.obj,
                               vars = NULL,
                               verbose = FALSE,
                               do_normalize = TRUE,
                               do_umap = FALSE)$Z |>
            Matrix::t()
          
          nn <-
            annoy_search(
              X = query,
              k = 10,
              ann = .tdr.obj@integration$symphony.obj$ref.knn.idx,
              #search_k = search_k,
              #prep_data = TRUE,
              #tmpdir = tmpdir,
              n_threads = .tdr.obj@config$n.threads,
              #grain_size = grain_size,
              verbose = FALSE
            )
          
          #res2$cell.celltyping <-
          #  as.character(x = .tdr.obj@integration$symphony.obj$meta_data[[
          #    .tdr.obj@integration$symphony.obj$celltype.col.name
          #  ]])[nn$idx[,1,drop = TRUE]] |>
          #  stats::setNames(nm = colnames(x = raw.exprs))  
          
          #res2$lm.celltyping <-
          #  as.character(x = .tdr.obj@integration$symphony.obj$meta_data[[
          #    .tdr.obj@integration$symphony.obj$celltype.col.name
          #  ]])[nn$idx[!cells.of.interest,1,drop = TRUE]] |>
          #  stats::setNames(nm = colnames(x = raw.exprs)[!cells.of.interest])  
          
          res2$cell.celltyping <-
            .tdr_transfer_labels(
              .method = "knn_vote",
              .label.confidence = .label.confidence,
              .n.cells = ncol(x = raw.exprs),
              .cell.names = colnames(x = raw.exprs),
              .nn.idx = nn$idx,
              .ref.labels = .tdr.obj@integration$symphony.obj$meta_data[[
                .tdr.obj@integration$symphony.obj$celltype.col.name
              ]]
            )
          
          res2$lm.celltyping <-
            res2$cell.celltyping[!cells.of.interest]
          
        }
        
        
        # Prefix cell names with sample name for both RNA and cyto
        names(x = res2$cell.clustering) <-
          paste0(names(x = .tdr.obj@cells)[cells.idx],
                 "_",
                 names(x = res2$cell.clustering))
        
        if(!is.null(x = res2$cell.celltyping)){
          names(x = res2$cell.celltyping) <-
            paste0(names(x = .tdr.obj@cells)[cells.idx],
                   "_",
                   names(x = res2$cell.celltyping))
        }
        
        if(!is.null(x = res2$lm.celltyping)){
          names(x = res2$lm.celltyping) <-
            paste0(names(x = .tdr.obj@cells)[cells.idx],
                   "_",
                   names(x = res2$lm.celltyping))
        }
        
        # ── Compute this sample's fdens column (before eviction) ──
        smpl.name <- names(x = .tdr.obj@cells)[cells.idx]
        
        smpl.fdens <-
          if(ncol(x = res2$fgraph) == nrow(x = .tdr.obj@assay$expr)){
            Matrix::colSums(x = res2$fgraph[
              if(!is.null(x = .ct.to.keep)){
                (res2$cell.celltyping %in% .ct.to.keep)
              } else {
                1:nrow(x = res2$fgraph)
              }, , drop = FALSE])
          } else {
            Matrix::rowSums(x = res2$fgraph[,
              if(!is.null(x = .ct.to.keep)){
                (res2$cell.celltyping %in% .ct.to.keep)
              } else {
                1:nrow(x = res2$fgraph)
              }, drop = FALSE])
          }
        
        smpl.n.cells <- nrow(x = res2$embedding)
        
        # ── Compute streaming cell.count for this sample ──
        smpl.cl.count <-
          table(res2$cell.clustering) |>
          (\(x) stats::setNames(as.vector(x), names(x)))()
        
        smpl.ct.count <-
          if(!is.null(x = res2$cell.celltyping)){
            table(res2$cell.celltyping) |>
              (\(x) stats::setNames(as.vector(x), names(x)))()
          } else {
            NULL
          }
        
        # ── Cache large slots to disk or keep in-memory ──
        smpl.cache.records <- NULL
        if (isTRUE(x = .cache.on.disk)) {
          
          smpl.cache.records <- list(
            clustering = .tdr_cache_write(
              object = res2$cell.clustering,
              cache_dir = .run_cache_dir,
              slot_name = "clustering",
              sample_name = smpl.name),
            nearest.lm = .tdr_cache_write(
              object = res2$nn$euclidean$idx,
              cache_dir = .run_cache_dir,
              slot_name = "nearest.lm",
              sample_name = smpl.name),
            fuzzy.graphs = .tdr_cache_write(
              object = res2$fgraph,
              cache_dir = .run_cache_dir,
              slot_name = "fuzzy.graphs",
              sample_name = smpl.name)
          )
          
          if(!is.null(x = res2$cell.celltyping)){
            smpl.cache.records$celltyping <- .tdr_cache_write(
              object = res2$cell.celltyping,
              cache_dir = .run_cache_dir,
              slot_name = "celltyping",
              sample_name = smpl.name)
          }
          
        }
        
        if(isTRUE(x = .verbose)){
          .show_progress(current = cells.idx, 
                         total = length(.tdr.obj@cells),
                         item_label = names(x = .tdr.obj@cells)[cells.idx],
                         start_time = .map_start)
        }
        
        # Return lightweight result (large objects are on disk if caching enabled)
        list(
          fdens.col      = smpl.fdens,
          n.cells        = smpl.n.cells,
          cl.count       = smpl.cl.count,
          ct.count       = smpl.ct.count,
          cell.clustering  = res2$cell.clustering,
          cell.celltyping  = res2$cell.celltyping,
          lm.celltyping    = res2$lm.celltyping,
          fgraph           = res2$fgraph,
          nn.idx           = res2$nn$euclidean$idx,
          cache.records    = smpl.cache.records
        )
        
      })
    
    if(!is.null(x = .tdr.obj@integration$symphony.obj)){
      
      .tdr.obj@landmark.annot$celltyping <-
        vector(mode = "list",
               length = 0)
      
      .tdr.obj@landmark.annot$celltyping$ids <-
        seq_along(along.with = res) |>
        stats::setNames(nm = names(x = res)) |>
        lapply(FUN = function(smpl.idx){
          res[[smpl.idx]]$lm.celltyping
        }) |>
        unlist(use.names = FALSE)
      
      names(x = .tdr.obj@landmark.annot$celltyping$ids) <-
        seq_along(along.with = res) |>
        stats::setNames(nm = names(x = res)) |>
        lapply(FUN = function(smpl.idx){
          names(x = res[[smpl.idx]]$lm.celltyping)
        }) |>
        unlist(use.names = FALSE)
      
      .tdr.obj@landmark.annot$celltyping$ids <-
        .tdr.obj@landmark.annot$celltyping$ids[
          rownames(x = .tdr.obj@assay$expr)
        ] |>
        as.factor()
      
      # Store as named solution
      .tdr.obj@landmark.annot$celltyping[[.celltype.col.name]] <-
        .tdr.obj@landmark.annot$celltyping$ids
      
    }
    
    # ── Assemble fdens from pre-computed per-sample columns ──
    fdens <-
      lapply(X = res, FUN = `[[`, "fdens.col") |>
      do.call(what = cbind) |>
      (\(x)
       `rownames<-`(x = x,
                    value = rownames(x = .tdr.obj@assay$expr))
      )() |>
      Matrix::t() |>
      (\(x) {
        n.cells <- vapply(res, `[[`, numeric(1), "n.cells")
        x / (n.cells / mean(n.cells))
      })() |>
      Matrix::t()
    
    # Compute log2-transformed landmark densities for statistical modeling
    Y <-
      log2(x = fdens + 0.5)
    
    # ── Assemble @cellmap: cache-aware ──
    if (isTRUE(x = .cache.on.disk)) {
      
      # Populate @cellmap with path strings from per-sample records
      clustering_ids <- stats::setNames(
        vector("list", length(res)), names(res))
      nearest_lm     <- clustering_ids
      fuzzy_graphs   <- clustering_ids
      celltyping_ids <- clustering_ids
      has_celltyping <- FALSE
      
      for (sn in names(res)) {
        recs <- res[[sn]]$cache.records
        clustering_ids[[sn]] <- recs$clustering
        nearest_lm[[sn]]    <- recs$nearest.lm
        fuzzy_graphs[[sn]]  <- recs$fuzzy.graphs
        if (!is.null(recs$celltyping)) {
          celltyping_ids[[sn]] <- recs$celltyping
          has_celltyping <- TRUE
        }
      }
      
      .tdr.obj@density$fdens <- fdens
      .tdr.obj@density$Y <- Y
      .tdr.obj@cellmap$clustering$ids  <- clustering_ids
      .tdr.obj@cellmap$nearest.lm      <- nearest_lm
      .tdr.obj@cellmap$fuzzy.graphs    <- fuzzy_graphs
      if (has_celltyping) {
        .tdr.obj@cellmap$celltyping$ids <- celltyping_ids
      }
      .tdr.obj@config$.cache.root <- .run_cache_dir
      
    } else {
      
      # In-memory path (original behaviour)
      cell.fg <-
        lapply(X = res, FUN = `[[`, "fgraph")
      
      cell.nlmn <-
        lapply(X = res, FUN = `[[`, "nn.idx")
      
      cell.clustering <-
        lapply(X = res, FUN = `[[`, "cell.clustering")
      
      if(!is.null(x = .tdr.obj@landmark.annot$celltyping$ids)){
        cell.celltyping <-
          lapply(X = res, FUN = `[[`, "cell.celltyping")
      } else {
        cell.celltyping <- NULL
      }
      
      .tdr.obj@density$fdens <- fdens
      .tdr.obj@density$Y <- Y
      .tdr.obj@cellmap$clustering$ids  <- cell.clustering
      .tdr.obj@cellmap$celltyping$ids  <- cell.celltyping
      .tdr.obj@cellmap$nearest.lm      <- cell.nlmn
      .tdr.obj@cellmap$fuzzy.graphs    <- cell.fg
    }
    
    # ── Compute cell.count / cell.perc from streaming counts ──
    .tdr.obj@density$composition$clustering$cell.count <-
      lapply(X = res, FUN = `[[`, "cl.count") |>
      stats::setNames(nm = names(res)) |>
      dplyr::bind_rows(.id = "sample") |>
      as.data.frame() |>
      (\(x)
       `rownames<-`(x = as.matrix(x = x[,colnames(x = x) != "sample"]),
                    value = x$sample)
      )() |>
      (\(x)
       x[,colnames(x = x) |>
           sort()]
      )()
    
    .tdr.obj@density$composition$clustering$cell.count[
      is.na(x = .tdr.obj@density$composition$clustering$cell.count)
    ] <- 0
    
    .tdr.obj@density$composition$clustering$cell.perc <-
      (.tdr.obj@density$composition$clustering$cell.count * 100) /
      Matrix::rowSums(x = .tdr.obj@density$composition$clustering$cell.count)
    
    if(!is.null(x = .tdr.obj@landmark.annot$celltyping$ids)){
      
      .tdr.obj@density$composition$celltyping$cell.count <-
        lapply(X = res, FUN = `[[`, "ct.count") |>
        stats::setNames(nm = names(res)) |>
        dplyr::bind_rows(.id = "sample") |>
        as.data.frame() |>
        (\(x)
         `rownames<-`(x = as.matrix(x = x[,colnames(x = x) != "sample"]),
                      value = x$sample)
        )()  |>
        (\(x)
         x[,colnames(x = x) |>
             sort()]
        )()
      
      .tdr.obj@density$composition$celltyping$cell.count[
        is.na(x = .tdr.obj@density$composition$celltyping$cell.count)
      ] <- 0
      
      .tdr.obj@density$composition$celltyping$cell.perc <-
        (.tdr.obj@density$composition$celltyping$cell.count * 100) /
        Matrix::rowSums(x = .tdr.obj@density$composition$celltyping$cell.count)
      
    }
    
    .tdr.obj@density$ignored <- NULL
    
    return(.tdr.obj)
    
  }

#' Leiden clustering of landmarks
#'
#' Applies Leiden community detection to partition landmarks into clusters based on either 
#' the SNN graph or UMAP fuzzy graph. Clusters are used for visualization and as a coarser 
#' grouping for statistical testing. This is a wrapper around \code{leiden.cluster}.
#'
#' @param .cl.method Character specifying clustering method: "snn" (shared nearest neighbors) 
#'   or "fgraph" (fuzzy graph from UMAP). Default "snn".
#' @param .cl.resolution.parameter Numeric resolution for Leiden clustering (default 0.8). 
#'   Internally scaled by 1e-3 for CPM objective. Higher values 
#'   yield finer/more clusters.
#' @param .seed Integer seed for reproducibility (default 123).
#' @param .verbose Logical for progress messages (default TRUE).
#' @param .small.size Integer threshold for straggler absorption (default 3). Clusters smaller 
#'   than this are merged into most connected neighbors.
#'   
#' @return Updated \code{.tdr.obj} with \code{$graph$clustering} containing:
#'   \itemize{
#'     \item \code{ids}: Factor of cluster assignments for each landmark
#'     \item \code{median.exprs}: Matrix of mean expression per cluster. For RNA: top features 
#'       from PCA rotation (6 genes per PC: 3 highest, 3 lowest loadings). For cytometry: all markers.
#'     \item \code{pheatmap}: Heatmap object showing cluster profiles
#'   }
#'   
#' @details
#' For RNA data, the heatmap displays top contributors to each principal component (6 genes 
#' per PC: 3 most positive, 3 most negative loadings). For cytometry, all markers are shown.
#' 
#' The fuzzy graph method ("fgraph") uses probabilistic UMAP edges (pruned with tolerance 1/20 
#' to remove weak connections) and may produce more visually coherent clusters in UMAP space. 
#' The SNN method uses Jaccard similarity of k-NN overlaps and often gives more stable results 
#' for downstream statistics.
#' 
#' @seealso \code{\link{leiden.cluster}}, \code{\link{get.graph}}
#' 
#' @examples
#' \dontrun{
#' # After graph construction
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#' 
#' # Extract clustering (called internally by get.graph)
#' lm.cells <- lm.cluster(lm.cells)
#' 
#' # Higher resolution clustering
#' lm.cells <- lm.cluster(lm.cells, .cl.resolution.parameter = 200)
#' 
#' # Use fuzzy graph instead of SNN
#' lm.cells <- lm.cluster(lm.cells, .cl.method = "fgraph")
#' }
#' 
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object.
#' @param ... Additional arguments passed to methods.
#' @export
lm.cluster <- function(x, ...) UseMethod("lm.cluster")

#' @rdname lm.cluster
#' @export
lm.cluster.TDRObj <-
  function(x,
           .cl.method = "snn",
           .cl.resolution.parameter = 0.8,
           .seed = 123,
           .verbose = TRUE,
           .small.size = 3,
           .column.name = NULL,
           ...){
    .tdr.obj <- x
    
    # R CMD check appeasement
    cell.pop <- NULL
    
    # Default column name from resolution parameter
    if (is.null(.column.name)) {
      .column.name <- paste0("leiden.res.", .cl.resolution.parameter)
    }
    if (.column.name == "ids") {
      stop("'.column.name' cannot be 'ids' (reserved for the active solution).")
    }
    
    set.seed(seed = .seed)
    
    if(is.null(.tdr.obj@graphs$adj.matrix)){
      stop("Graph component missing. Run get.graph() before clustering.")
    }
    
    # Validate clustering method
    .cl.method <-
      match.arg(arg = .cl.method,
                choices = c("fgraph","snn"))
    
    # Run Leiden clustering on selected similarity matrix
    .tdr.obj@landmark.annot$clustering$ids <-
      leiden.cluster(
        .tdr.obj = .tdr.obj,
        .sim.matrix = if(.cl.method == "fgraph"){
          Matrix::drop0(x = .tdr.obj@graphs$fgraph,   # Prune weak fuzzy edges
                        tol = 1/20)
        } else {
          .tdr.obj@graphs$snn                              # Use SNN graph
        },
        .resolution.parameter = .cl.resolution.parameter * 1e-3,  # Scale for CPM
        .small.size = .small.size,
        .verbose = .verbose)
    
    # Store as named solution
    .tdr.obj@landmark.annot$clustering[[.column.name]] <-
      .tdr.obj@landmark.annot$clustering$ids
    
    return(.tdr.obj)
  }

#' Refresh all clustering-dependent downstream slots
#'
#' Mirrors \code{.refresh_celltyping} but for clustering annotations.
#' Called after reclustering or switching the active clustering solution.
#'
#' @param .tdr.obj A TDRObj with \code{@landmark.annot$clustering$ids}
#'   already set to the new values.
#' @param .verbose Logical; print progress messages.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.refresh_clustering <- function(.tdr.obj, .verbose = FALSE) {
  
  # --- Cell-level IDs & composition (only if get.map() has run) ---
  if (!is.null(.tdr.obj@density$fdens)) {
    
    if (isTRUE(.verbose)) {
      message("-> get.map() results detected; refreshing cell-level cluster IDs and composition...")
    }
    
    .tdr.obj <- .relabel_cellmap(.tdr.obj, .annot.type = "clustering", .verbose = .verbose)
    .tdr.obj <- .recompute_composition(.tdr.obj, .annot.type = "clustering")
    
    if (isTRUE(.verbose)) {
      message("-> Cell-level cluster IDs and composition matrices refreshed.")
    }
  }
  
  # --- trad clustering fit (only if get.lm() was run) ---
  if (!is.null(.tdr.obj@results$lm)) {
    .tdr.obj <- .invalidate_trad(.tdr.obj, .annot.type = "clustering")
  }
  
  # --- Warn about potentially stale pbDE / markerDE results ---
  .warn_stale_de_results(.tdr.obj, .annot.type = "clustering")
  
  # --- Warn if celltyping exists (may now be inconsistent) ---
  if (!is.null(.tdr.obj@landmark.annot$celltyping$ids)) {
    warning("Celltyping annotations exist and may be inconsistent with the new clustering. ",
            "Consider re-running celltyping().", call. = FALSE)
  }
  
  .tdr.obj
}

#' Recluster landmarks
#'
#' Re-runs Leiden community detection with new parameters and refreshes all
#' clustering-dependent downstream slots (cell-level IDs, composition matrices,
#' traditional fits). This is the recommended way to explore different
#' resolutions after the initial \code{get.graph()} call.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object with \code{get.graph()} already run.
#' @param .cl.resolution.parameter Numeric resolution for Leiden CPM (default 0.8).
#'   Internally scaled by 1e-3. Higher = finer clusters.
#' @param .cl.method Character: \code{"snn"} or \code{"fgraph"} (default \code{"snn"}).
#' @param .small.size Integer threshold for straggler absorption (default 3).
#' @param .column.name Character name for storing the solution (default
#'   \code{paste0("leiden.res.", .cl.resolution.parameter)}). Cannot be \code{"ids"}.
#' @param .seed Integer for reproducibility (default 123).
#' @param .verbose Logical (default TRUE).
#' @param ... Additional arguments passed to methods.
#' @return Updated object with new active clustering and refreshed downstream slots.
#' @export
recluster <- function(x, ...) UseMethod("recluster")

#' @rdname recluster
#' @export
recluster.TDRObj <-
  function(x,
           .cl.resolution.parameter = 0.8,
           .cl.method = "snn",
           .small.size = 3,
           .column.name = NULL,
           .seed = 123,
           .verbose = TRUE,
           ...) {
    
    .tdr.obj <- lm.cluster.TDRObj(
      x,
      .cl.method = .cl.method,
      .cl.resolution.parameter = .cl.resolution.parameter,
      .seed = .seed,
      .verbose = .verbose,
      .small.size = .small.size,
      .column.name = .column.name
    )
    
    .tdr.obj <- .refresh_clustering(.tdr.obj, .verbose = .verbose)
    
    return(.tdr.obj)
  }

#' Set active clustering solution
#'
#' Switches the active clustering to a previously stored solution and
#' refreshes all clustering-dependent downstream slots.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object.
#' @param .column.name Character: name of the stored clustering solution to activate.
#' @param .verbose Logical (default TRUE).
#' @param ... Additional arguments passed to methods.
#' @return Updated object with the selected clustering as active \code{$ids}.
#' @export
set_active_clustering <- function(x, ...) UseMethod("set_active_clustering")

#' @rdname set_active_clustering
#' @export
set_active_clustering.TDRObj <-
  function(x,
           .column.name,
           .verbose = TRUE,
           ...) {
    .tdr.obj <- x
    
    if (.column.name == "ids") {
      stop("'.column.name' cannot be 'ids'. Specify a named solution.")
    }
    
    stored <- .tdr.obj@landmark.annot$clustering[[.column.name]]
    if (is.null(stored)) {
      avail <- setdiff(names(.tdr.obj@landmark.annot$clustering), "ids")
      stop("Clustering solution '", .column.name, "' not found.\n",
           "Available solutions: ", paste(avail, collapse = ", "))
    }
    
    .tdr.obj@landmark.annot$clustering$ids <- stored
    .tdr.obj <- .refresh_clustering(.tdr.obj, .verbose = .verbose)
    
    return(.tdr.obj)
  }

#' Graph-based feature discovery for landmarks
#'
#' Identifies the most characteristic features (genes or markers) for each landmark by 
#' determining which features contribute most strongly to the landmark's defining principal 
#' component. This creates landmark-specific feature signatures useful for interpretation.
#'
#' @return Updated \code{.tdr.obj} with \code{$interact.plot$lm.features} containing:
#'   \itemize{
#'     \item \code{res}: List of data frames (one per landmark) with top 10 feature importances
#'     \item \code{html}: HTML-formatted tables for interactive plotting
#'   }
#'   Feature importance is the signed loading from the landmark's top principal component.
#'   
#' @details
#' \strong{Algorithm:}
#' 
#' 1. For each landmark, identify the PC with highest absolute value (its "defining" PC)
#' 2. Extract that PC's rotation vector (feature loadings)
#' 3. Multiply by sign of landmark's PC score to preserve directionality
#' 4. Rank features by signed loading; keep top 10
#' 
#' Positive loadings indicate features that increase along the PC direction, negative indicate 
#' features that decrease. The signed values reflect how the landmark is characterized relative 
#' to the population average.
#' 
#' Results are stored for use in \code{plotPCA} and \code{plotUMAP}, which displays feature signatures 
#' when hovering over landmarks when \code{.hover.stats = "marker"}.
#' 
#' @seealso \code{\link{plotPCA}}, \code{\link{plotUMAP}}, \code{\link{get.landmarks}}
#' 
#' @examples
#' \dontrun{
#' # After landmark selection and graph construction
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.features()
#' 
#' # View feature signature for first landmark
#' head(lm.cells$interact.plot$lm.features$res[[1]])
#' 
#' # Use in interactive plotting
#' plotPCA(lm.cells, .hover.stats = "marker")
#' }
#' 
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object.
#' @param ... Additional arguments passed to methods.
#' @export
get.features<- function(x, ...) UseMethod("get.features")

#' @rdname get.features
#' @method get.features TDRObj
#' @export
get.features.TDRObj <-
  function(
    x,
    ...
  ){
    .tdr.obj <- x
    
    # Identify the top (dominant) PC for each landmark
    top.comp <-
      .tdr.obj@landmark.embed$pca$coord |>
      abs() |>
      as.matrix() |>
      (\(x)
       matrixStats::rowRanks(x = x,
                             ties.method = "max") ==   # Max ensures ties get highest rank
         ncol(x = x)                                   # Compare to last column (highest rank)
      )() |>
      which(arr.ind = TRUE) |>
      (\(x)
       x[order(x = x[,"row"]),"col"]                   # Extract column indices
      )()
    
    # Get sign of each landmark's score in its top PC
    top.comp.sign <-
      diag(x = .tdr.obj@landmark.embed$pca$coord[,top.comp]) |>
      sign()
    
    # Signed feature loadings: rotation × sign for each landmark's top PC
    coefs <-
      ((Matrix::t(x = .tdr.obj@landmark.embed$pca$rotation[,top.comp]) * top.comp.sign) |>
         Matrix::t())
    
    # Find top 10 features for each landmark (by signed loading)
    hits <-
      coefs |>
      (\(x)
       matrixStats::colRanks(x = x,
                             ties.method = "max",
                             preserveShape = TRUE) >
         (nrow(x = x) - 10)                            # Top 10 = rank > n-10
      )() |>
      which(arr.ind = TRUE) |>
      as.data.frame()
    
    # Build result list: one data frame per landmark
    res <-
      nrow(x = .tdr.obj@assay$expr) |>
      seq_len() |>
      stats::setNames(nm = rownames(x = .tdr.obj@assay$expr)) |>
      lapply(FUN = function(lm.idx){
        
        coefs[hits$row[hits$col == lm.idx],            # Features for this landmark
              lm.idx] |>
          sort(decreasing = TRUE) |>                   # Order by importance
          as.data.frame() |>
          (\(x)
           `colnames<-`(x = x,
                        value = "feature importance")
          )()
        
      })
    
    # Generate HTML tables for interactive display
    html <-
      lapply(X = res,
             FUN = knitr::kable,
             format = "html")
    
    .tdr.obj@results$features$lm.features <-
      list(res = res,
           html = html)
    
    return(.tdr.obj)
  }

# Internal functions imported from uwot package
# These are not exported by uwot but used for Laplacian Eigenmap computation

#' @keywords internal
annoy_build <-
  getFromNamespace(x = "annoy_build",
                   ns = "uwot")

#' @keywords internal
annoy_search <-
  getFromNamespace(x = "annoy_search",
                   ns = "uwot")

#' @keywords internal
annoy_nn <-
  getFromNamespace(x = "annoy_nn",
                   ns = "uwot")

#' @keywords internal
form_modified_laplacian <-
  getFromNamespace(x = "form_modified_laplacian",
                   ns = "uwot")

#' @keywords internal
irlba_spectral_tsvd <-
  getFromNamespace(x = "irlba_spectral_tsvd",
                   ns = "uwot")


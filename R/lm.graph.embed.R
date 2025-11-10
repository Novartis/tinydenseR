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
#' Uses Laplacian Eigenmap embedding (if available in \code{.lm.obj}), otherwise
#' PCA embedding (first 3 components), or computes embedding from similarity matrix.
#' K-means initialization helps Leiden start from reasonable communities rather than
#' singleton nodes, improving stability.
#' 
#' **Straggler absorption:**
#' Small clusters (< 0.5% of total by default) often represent noise or transitional
#' states. Rather than keeping them as separate clusters, they're merged with the
#' most similar neighboring cluster based on mean edge weights.
#'
#' @param .lm.obj Optional tinydenseR object from \code{get.landmarks()} and
#'   \code{get.graph()}. If provided, uses Laplacian Eigenmap embedding for
#'   initialization. If NULL, computes 3D embedding from \code{.sim.matrix}.
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
#'   .lm.obj = lm.obj,
#'   .sim.matrix = jaccard_sim,
#'   .resolution.parameter = 0.001,
#'   .small.size = 10  # Merge clusters < 10 cells
#' )
#' table(clusters)
#' }
leiden.cluster <-
  function(.lm.obj = NULL,
           .sim.matrix,
           .resolution.parameter,
           .small.size = floor(nrow(.sim.matrix) / 200),
           .verbose = TRUE,
           .seed = 123) {
    
    # Determine initialization embedding: LE > PCA > computed from similarity
    .init.embed <- 
      if(!is.null(x = .lm.obj)){
        if(is.na(x = .lm.obj$graph$LE$embed[1,1])) {
          # Fallback to PCA if LE not available
          .lm.obj$pca$embed[,1:min(3, ncol(x = .lm.obj$pca$embed)), drop = FALSE]
        } else {
          .lm.obj$graph$LE$embed
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
      kres$cluster |>
      as.integer()
    
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
#' Converts a nearest neighbor index matrix into a sparse adjacency matrix. Each entry (i,j) 
#' indicates that cell j is a nearest neighbor of cell i. The resulting matrix is non-symmetric 
#' (i can be neighbor of j without j being neighbor of i).
#' 
#' @param .nn.idx A matrix where each row represents a cell and each column represents a nearest 
#'   neighbor index. Element (i,k) is the index of the k-th nearest neighbor of cell i.
#'   
#' @return A non-symmetric sparse matrix (dgCMatrix) of dimensions n×n where n is the number of 
#'   cells. Entry (i,j) = 1 if j is among i's nearest neighbors, 0 otherwise.
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
#'   entry (i,j)=1 indicates j is a nearest neighbor of i.
#' @param .prune A numeric threshold (default 1/15 ≈ 0.067) for pruning weak edges. Jaccard 
#'   values below this are set to zero for sparsity.
#'   
#' @return A symmetric sparse matrix of Jaccard indices representing the SNN graph. Higher values 
#'   indicate greater neighbor overlap. Matrix is pruned and forced symmetric.
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
#' Performs the core graph analysis workflow: UMAP embedding, Laplacian Eigenmap (LE) computation, 
#' and clustering of landmark cells. This creates the low-dimensional representation used for 
#' visualization and cluster identification.
#'
#' @param .lm.obj A tinydenseR object initialized with \code{setup.lm.obj} and processed with 
#'   \code{get.landmarks}. Must contain \code{$lm} (landmark matrix) and \code{$pca} components.
#' @param .k Integer number of nearest neighbors for graph construction (default 20). Higher values 
#'   create more connected graphs with broader cluster definitions.
#' @param .scale Logical indicating whether to scale features before UMAP. Defaults to FALSE for 
#'   RNA assays (PCA already scaled) and TRUE for cytometry (feature scales vary).
#' @param .verbose Logical for progress messages (default TRUE).
#' @param .seed Integer seed for reproducibility of UMAP and clustering (default 123).
#' @param .cl.method Character specifying clustering method: "snn" (shared nearest neighbors via 
#'   Jaccard similarity) or "fgraph" (fuzzy graph from UMAP). Default "snn".
#' @param .cl.resolution.parameter Numeric controlling cluster granularity (default 0.8). Higher 
#'   values produce more fine-grained clusters.
#' @param .small.size Integer threshold for straggler clusters (default 0.5% of landmarks). Clusters 
#'   smaller than this are absorbed into nearest large cluster.
#'   
#' @return Updated \code{.lm.obj} with \code{$graph} component containing:
#'   \itemize{
#'     \item \code{uwot}: UMAP model with embedding, nearest neighbors, fuzzy graph
#'     \item \code{LE}: Laplacian Eigenmap embedding and eigenspectrum
#'     \item \code{snn}: Shared nearest neighbor graph
#'     \item \code{clustering}: Cluster assignments from \code{lm.cluster}
#'   }
#'   
#' @details
#' The function orchestrates several graph-based analyses:
#' 
#' \strong{1. UMAP Embedding:}
#' Computes 2D embedding for visualization using harmony-corrected PCA (if available) or raw data.
#' Returns the UMAP model for later projection of query cells.
#' 
#' \strong{2. Graph Construction:}
#' Builds k-NN graph and computes SNN graph via Jaccard similarity of neighbor overlaps.
#' 
#' \strong{3. Laplacian Eigenmap:}
#' Computes spectral embedding by solving the generalized eigenvalue problem of the graph Laplacian.
#' Automatically selects dimensionality via elbow detection on eigenvalue spectrum. Provides 
#' alternative to PCA that respects local graph structure.
#' 
#' \strong{4. Clustering:}
#' Applies Leiden algorithm to identify communities, with straggler absorption for robustness.
#' 
#' @seealso \code{\link{setup.lm.obj}}, \code{\link{get.landmarks}}, \code{\link{lm.cluster}}
#' 
#' @examples
#' \dontrun{
#' # Typical workflow after landmark selection
#' lm.cells <- setup.lm.obj(.cells = .cells, 
#'                           .meta = .meta,
#'                           .assay.type = "RNA") |>
#'   get.landmarks(.nHVG = 500, .nPC = 3) |>
#'   get.graph(.k = 10)
#' 
#' # Higher resolution for finer clusters
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph(.cl.resolution.parameter = 200)
#' 
#' # Use fuzzy graph for clustering instead of SNN
#' lm.cells <- get.graph(lm.cells, .cl.method = "fgraph")
#' }
#' 
#' @export
get.graph <-
  function(.lm.obj,
           .k = 20,
           .scale = if(.lm.obj$assay.type == "RNA") FALSE else TRUE,
           .verbose = TRUE,
           .seed = 123,
           .cl.method = "snn",
           .cl.resolution.parameter = 0.8,
           .small.size = floor(x = nrow(x = .lm.obj$lm) / 200)){
    
    # Validate clustering method
    .cl.method <-
      match.arg(arg = .cl.method,
                choices = c("fgraph","snn"))
    
    # Compute UMAP embedding with k-NN graph
    set.seed(seed = .seed)
    .lm.obj$graph$uwot <-
      uwot::umap(X = if(!is.null(x = .lm.obj$harmony.obj) || .lm.obj$assay.type == "RNA") .lm.obj$pca$embed else .lm.obj$lm,
                 n_neighbors = .k,
                 n_components = 2,           # 2D for visualization
                 n_epochs = 500,             # Training iterations
                 scale = .scale,
                 pca = NULL,                 # Already have PCA/LE
                 verbose = isTRUE(x = .verbose),
                 ret_model = TRUE,           # Keep model for query projection
                 batch = TRUE,
                 seed = .seed,
                 n_threads = .lm.obj$n.threads,
                 fast_sgd = FALSE,
                 n_sgd_threads = .lm.obj$n.threads,
                 ret_extra = c("fgraph",     # Fuzzy graph (for optional clustering)
                               "nn"))         # Nearest neighbors
    
    colnames(x = .lm.obj$graph$uwot$embedding) <-
      c("umap.1",
        "umap.2")
    
    # Convert k-NN indices to sparse adjacency matrix
    .lm.obj$graph$adj.matrix <-
      get.adj.matrix(.nn.idx = .lm.obj$graph$uwot$nn$euclidean$idx)
    
    # Compute SNN graph via Jaccard similarity of neighbor overlaps
    .lm.obj$graph$snn <-
      fast.jaccard.r(.adj.matrix = .lm.obj$graph$adj.matrix,
                     .prune = 1/15) |>
      (\(x)
       `dimnames<-`(x = x,
                    value = list(rownames(x = .lm.obj$lm),
                                 rownames(x = .lm.obj$lm)))
      )()
    
    if(isTRUE(x = .verbose)){
      message("getting Laplacian Eigenmap")
    }
    
    # Symmetrize adjacency: W[i,j] = 1 if i is NN of j OR j is NN of i
    # See uwot implementation: https://github.com/jlmelville/uwot/blob/f9e576e97d9df44d48be2cc559412282838dc4a5/R/init.R
    .lm.obj$graph$LE$W.sym <-
      (((.lm.obj$graph$adj.matrix > 0) |
          (Matrix::t(x = .lm.obj$graph$adj.matrix) > 0)) * 1) |>
      Matrix::Matrix(sparse = TRUE)
    
    # Target dimensions: match PCA dimensionality
    target_k <-
      ncol(x = .lm.obj$pca$embed)
    nv <-
      min(target_k, 
          nrow(x = .lm.obj$pca$embed) - 1L) # Eigenspectrum limited to n-1
    stopifnot(nv >= 2L)
    
    # Form modified graph Laplacian: L = D^(-1/2) (D - W) D^(-1/2)
    .lm.obj$graph$LE <- 
      c(.lm.obj$graph$LE,
        form_modified_laplacian(A = .lm.obj$graph$LE$W.sym,
                                ret_d = TRUE))
    
    # Compute truncated SVD of Laplacian for nv components
    .lm.obj$graph$LE <-
      c(.lm.obj$graph$LE,
        irlba_spectral_tsvd(L = .lm.obj$graph$LE$L, 
                            n = nv))
    
    # Verify convergence and valid eigenvectors
    ok <-
      all(c("vectors", "values", "converged") %in% 
            names(x = .lm.obj$graph$LE)) &&
      is.matrix(x = .lm.obj$graph$LE$vectors) &&
      (ncol(x = .lm.obj$graph$LE$vectors) >= nv) &&
      isTRUE(x = .lm.obj$graph$LE$converged)
    
    if(!ok) {
      
      message("Laplacian Eigenmap failed to converge.")
      
      .lm.obj$graph$LE$embed <- 
        matrix(data = NA_real_,
               nrow = 1,
               ncol = 1)
      
    } else {
      
      # Filter out trivial eigenvalues (near-zero, correspond to disconnected components)
      tol <- 1e-6
      .lm.obj$graph$LE$nontriv <-
        which(x = .lm.obj$graph$LE$values > tol)
      
      if(length(x = .lm.obj$graph$LE$nontriv) == 0L) {
        
        message("No non-trivial Laplacian eigenvalues found (graph may be empty or singular).")
        
        .lm.obj$graph$LE$embed <- 
          matrix(data = NA_real_,
                 nrow = 1,
                 ncol = 1)
        
      } else {
        
        # Determine embedding dimensionality via elbow detection
        if(length(x = .lm.obj$graph$LE$nontriv) < 4L){
          
          k <- length(x = .lm.obj$graph$LE$nontriv)  # Too few for elbow, use all
          
        } else {
          
          .lm.obj$graph$LE$elbow <-
            elbow.sec.deriv(x = .lm.obj$graph$LE$values[.lm.obj$graph$LE$nontriv],
                            smooth = TRUE,
                            df = NULL,
                            sort.order = "asc")$index
          
          # Use elbow dimensions, bounded by target and available eigenvectors
          k <- 
            min(.lm.obj$graph$LE$elbow, 
                target_k, 
                length(x = .lm.obj$graph$LE$nontriv)) |>
            max(2L)
          
        }
        
        # Extract first k non-trivial eigenvectors
        .lm.obj$graph$LE$keep <-
          .lm.obj$graph$LE$nontriv[seq_len(length.out = k)]
        
        # Normalize eigenvectors: multiply by D^(-1/2) and L2-normalize columns
        .lm.obj$graph$LE$embed <-
          (.lm.obj$graph$LE$Disqrt * .lm.obj$graph$LE$vectors[,.lm.obj$graph$LE$keep, drop = FALSE]) |> 
          (\(x)
           sweep(x = x, 
                 MARGIN = 2, 
                 STATS = Matrix::colSums(x * x) |>
                   sqrt(), 
                 FUN = `/`)
          )() |>
          (\(x)
           `dimnames<-`(x = x,
                        value = list(rownames(x = .lm.obj$pca$embed),
                                     paste0("LE", 1:ncol(x = x))))
          )()
      }
    }
    
    if(isTRUE(x = .verbose)){
      message("clustering")
    }
    
    # Apply Leiden clustering with straggler absorption
    .lm.obj <-
      lm.cluster(.lm.obj = .lm.obj,
                 .cl.method = .cl.method,
                 .cl.resolution.parameter = .cl.resolution.parameter,
                 .seed = .seed,
                 .verbose = .verbose,
                 .small.size = .small.size)
    
    return(.lm.obj)
    
  }

#' Mapping cells to landmarks
#'
#' Projects all cells onto the landmark graph to compute probabilities of cell-landmark neighborhood.
#' In addition, transfers cluster/cell type labels from landmarks to all cells.
#' 
#'
#' @param .lm.obj A tinydenseR object with \code{$graph} component populated by \code{get.graph}.
#' @param .ref.obj Optional Symphony reference object for cell type annotation. Must have 
#'   \code{Z_corr} field (harmony-corrected embeddings) and metadata with cell type labels. 
#'   Only compatible with RNA assays. Replaces any existing \code{$graph$celltyping}.
#' @param .integrate.vars Character vector of batch variables for integration (metadata columns). 
#'   Currently not used in function body.
#' @param .celltype.col.name Column name in \code{.ref.obj$meta_data} containing cell type labels 
#'   (default "cell_type"). Only relevant when \code{.ref.obj} is provided.
#' @param .cl.ct.to.ign Optional cluster or cell type name to exclude from downstream statistics. 
#'   Use for biologically irrelevant populations (e.g., erythrocytes in PBMC, which likely 
#'   represent technical artifacts). Cells in this population are mapped but excluded from density 
#'   calculations. Only one value allowed; merge multiple populations first via \code{celltyping}.
#' @param .irrel.par Character vector of irrelevant parameters to exclude. Currently not used.
#' @param .verbose Logical for progress messages (default TRUE).
#' @param .seed Integer seed for reproducibility (default 123).
#'   
#' @return Updated \code{.lm.obj} with \code{$map} component containing:
#'   \itemize{
#'     \item \code{fdens}: Matrix of fuzzy graph densities (landmarks × samples). Each entry is 
#'       the sum of cell-landmark edge weights for that sample, excluding ignored populations.
#'     \item \code{cell.clustering}: Named vector of cluster assignments for all cells.
#'     \item \code{cell.celltyping}: Named vector of cell type assignments (if celltyping available).
#'     \item \code{nn.idx}: Nearest landmark indices for all cells.
#'   }
#'   If \code{.ref.obj} provided, also updates \code{$graph$celltyping} with reference-based labels.
#'   
#' @details
#' \strong{Workflow Overview:}
#' 
#' For each sample, the function:
#' 1. Loads expression data and normalizes (size factors for RNA, marker subset for cytometry)
#' 2. Projects to PCA/Harmony space (matching landmark processing)
#' 3. Transforms to UMAP coordinates using the landmark UMAP model
#' 4. Assigns clusters/cell types by nearest landmark
#' 5. Computes fuzzy graph densities (cell-landmark edge weights)
#' 
#' \strong{Reference-Based Cell Typing:}
#' 
#' When \code{.ref.obj} is provided:
#' - Expression is mapped to reference via Symphony
#' - Cell types assigned by nearest neighbor in reference embedding
#' - Landmark cell types updated and used for visualization/statistics
#' - Replaces any existing manual cell type annotations
#' 
#' \strong{Fuzzy Graph Densities:}
#' 
#' The \code{fdens} matrix quantifies how strongly each landmark is connected to cells in each 
#' sample. High values indicate the landmark's neighborhood is enriched in that sample. This 
#' forms the basis for differential abundance testing in \code{get.stats}.
#' 
#' @seealso \code{\link{get.graph}}, \code{\link{get.stats}}, \code{\link{celltyping}}
#' 
#' @examples
#' \dontrun{
#' # Complete workflow with mapping
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500) |>
#'   get.graph() |>
#'   get.map()
#' 
#' # Exclude erythrocytes from PBMC statistics
#' lm.cells <- get.map(lm.cells, .cl.ct.to.ign = "Erythrocytes")
#' 
#' # Use Symphony reference for cell typing (RNA data only)
#' ref <- readRDS("pbmc_reference.rds")
#' lm.cells <- get.map(lm.cells, 
#'                     .ref.obj = ref, 
#'                     .celltype.col.name = "cell_type")
#' }
#' 
#' @export
get.map <-
  function(.lm.obj,
           .ref.obj = NULL,
           .integrate.vars = NULL,
           .celltype.col.name = "cell_type",
           .cl.ct.to.ign = NULL,
           .irrel.par = NULL,
           .verbose = TRUE,
           .seed = 123){
    
    # R CMD check appeasement for non-standard evaluation in dplyr
    cell.pop <- id <- value <- ri <- NULL
    
    # Validate and setup Symphony reference if provided
    if(!is.null(x = .ref.obj)){
      
      if(.lm.obj$assay.type != "RNA"){
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
      
      warning("Using Symphony reference object for cell typing. Removing previous celltyping from .lm.obj.")
      
      if(isTRUE(x = .verbose)){
        message("building Annoy index for Symphony reference object")
      }
      
      .lm.obj$symphony.obj <-
        .ref.obj
      
      # Build k-NN index on reference embeddings for fast query lookup
      set.seed(seed = .seed)
      .lm.obj$symphony.obj$ref.knn.idx <-
        Matrix::t(x = .lm.obj$symphony.obj$Z_corr) |>
        annoy_build(metric = "euclidean", 
                    n_trees = 50,
                    verbose = .verbose)
      
      .lm.obj$symphony.obj$celltype.col.name <-
        .celltype.col.name
      
      .lm.obj$graph$celltyping <- NULL
      
    }
    
    # Validate and identify populations to exclude from statistics
    if(!is.null(x = .cl.ct.to.ign)){
      
      if(length(x = .cl.ct.to.ign) > 1){
        stop("Please provide only one cluster or cell type to ignore.\nTo exclude multiple populations, merge them first using celltyping().")
      }
      
      # Check if specified population exists in clustering, celltyping, or reference
      if(!(.cl.ct.to.ign %in% levels(x = .lm.obj$graph$clustering$ids)) &
         !((.cl.ct.to.ign %in% levels(x = .lm.obj$graph$celltyping$ids)) |
           ((.cl.ct.to.ign %in% .ref.obj$meta_data[[.celltype.col.name]])))){
        stop("The cluster or cell type '", .cl.ct.to.ign, "' was not found.\nCheck clustering IDs, celltyping IDs, or reference metadata.")
      }
      
      # Build keep lists excluding the ignored population
      .cl.to.keep <-
        levels(x = .lm.obj$graph$clustering$ids)[
          !(levels(x = .lm.obj$graph$clustering$ids) %in%
              .cl.ct.to.ign)
        ]
      
      .ct.to.keep <-
        if(is.null(x = .ref.obj)){
          levels(x = .lm.obj$graph$celltyping$ids)[
            !(levels(x = .lm.obj$graph$celltyping$ids) %in%
                .cl.ct.to.ign)
          ]
        } else {
          unique(x = .ref.obj$meta_data[[.celltype.col.name]])[
            !(unique(x = .ref.obj$meta_data[[.celltype.col.name]]) %in%
                .cl.ct.to.ign)
          ]
        }
      
    } else {
      
      # Keep all populations
      .cl.to.keep <-
        levels(x = .lm.obj$graph$clustering$ids)
      
      .ct.to.keep <-
        if(is.null(x = .ref.obj)){
          levels(x = .lm.obj$graph$celltyping$ids)
        } else {
          unique(x = .ref.obj$meta_data[[.celltype.col.name]])
        }
      
    }
    
    set.seed(seed = .seed)
    
    # Process each sample: load data, normalize, project to UMAP, assign labels
    res <-
      seq_along(along.with = .lm.obj$cells) |>
      stats::setNames(nm = names(x = .lm.obj$cells)) |> #_[1:2] |>
      lapply(FUN = function(cells.idx){
        
        set.seed(seed = .seed)
        
        mat.exprs <-
          readRDS(file = .lm.obj$cells[[cells.idx]])
        
        if(.lm.obj$assay.type == "RNA"){
          
          mat.exprs <-
            Matrix::t(x = mat.exprs) |>
            (\(x)
             # size factor normalization, taking into consideration size factor of landmarks
             (x / (Matrix::rowSums(x = x) /
                     mean(x = Matrix::rowSums(x = x)))) * 
               (mean(x = Matrix::rowSums(x = x)) / mean(x = Matrix::rowSums(x = .lm.obj$raw.lm)))
            )()
          
          mat.exprs@x <-
            log2(x = mat.exprs@x + 1)
          
          temp.cell.names <-
            paste0(names(x = .lm.obj$cells)[cells.idx],
                   "_",
                   rownames(x = mat.exprs))
          
        } else {
          
          mat.exprs <-
            mat.exprs[,.lm.obj$markers]
          
          temp.cell.names <-
            rownames(x = mat.exprs)
          
        }
        
        cells.of.interest <-
          !(temp.cell.names %in% rownames(x = .lm.obj$lm))
        
        if(!is.null(x = .lm.obj$harmony.obj)){
          
          # Map query
          mat.exprs <-
            symphony::mapQuery(exp_query = Matrix::t(x = mat.exprs[,.lm.obj$pca$HVG]),
                               metadata_query = matrix(data = 1,
                                                       nrow = nrow(x = mat.exprs),
                                                       ncol = 2),
                               ref_obj = .lm.obj$harmony.obj,
                               vars = NULL,
                               verbose = .verbose,
                               do_normalize = FALSE,
                               do_umap = FALSE)$Z |>
            Matrix::t() |> 
            (\(x)
             `rownames<-`(x = x,
                          value = rownames(x = mat.exprs))
            )()
          
        } else if(.lm.obj$assay.type == "RNA"){
          
          mat.exprs <-
            (((Matrix::t(x = mat.exprs[,.lm.obj$pca$HVG]) - .lm.obj$pca$center) /
                .lm.obj$pca$scale) |>
               Matrix::t()) %*%
            .lm.obj$pca$rotation
          
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
            model = .lm.obj$graph$uwot,
            init_weighted = TRUE,
            search_k = NULL,
            tmpdir = tempdir(),
            n_epochs = 0,
            n_threads = .lm.obj$n.threads,
            n_sgd_threads = .lm.obj$n.threads,
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
          as.character(x = .lm.obj$graph$clustering$ids[res2$nn$euclidean$idx[,1]]) |>
          stats::setNames(nm = rownames(x = mat.exprs))
        
        if(is.null(x = .lm.obj$symphony.obj)){
          
          if(!is.null(x = .lm.obj$graph$celltyping)){
            
            if(isTRUE(x = .verbose)){
              
              message("assigning cells to celltypes")
              
            }
            
            res2$cell.celltyping <-
              as.character(x = .lm.obj$graph$celltyping$ids[res2$nn$euclidean$idx[,1]]) |>
              stats::setNames(nm = rownames(x = mat.exprs))  
            
          }
          
        } else {
          
          if(isTRUE(x = .verbose)){
            message("assigning cells to celltypes using symphony reference obj.")
          }
          
          raw.exprs <-
            readRDS(file = .lm.obj$cells[[cells.idx]])
          
          # Map query
          query <-
            symphony::mapQuery(exp_query = raw.exprs,
                               metadata_query = matrix(data = 1,
                                                       nrow = ncol(x = raw.exprs),
                                                       ncol = 2),
                               ref_obj = .lm.obj$symphony.obj,
                               vars = NULL,
                               verbose = FALSE,
                               do_normalize = TRUE,
                               do_umap = FALSE)$Z |>
            Matrix::t()
          
          nn <-
            annoy_search(
              X = query,
              k = 1,
              ann = .lm.obj$symphony.obj$ref.knn.idx,
              #search_k = search_k,
              #prep_data = TRUE,
              #tmpdir = tmpdir,
              n_threads = .lm.obj$n.threads,
              #grain_size = grain_size,
              verbose = FALSE
            )
          
          res2$cell.celltyping <-
            as.character(x = .lm.obj$symphony.obj$meta_data[[
              .lm.obj$symphony.obj$celltype.col.name
            ]])[nn$idx[,1,drop = TRUE]] |>
            stats::setNames(nm = colnames(x = raw.exprs))  
          
          res2$lm.celltyping <-
            as.character(x = .lm.obj$symphony.obj$meta_data[[
              .lm.obj$symphony.obj$celltype.col.name
            ]])[nn$idx[!cells.of.interest,1,drop = TRUE]] |>
            stats::setNames(nm = colnames(x = raw.exprs)[!cells.of.interest])  
          
        }
        
        
        if(.lm.obj$assay.type == "RNA"){
          
          names(x = res2$cell.clustering) <-
            paste0(names(x = .lm.obj$cells)[cells.idx],
                   "_",
                   names(x = res2$cell.clustering))
          
          if(!is.null(x = res2$cell.celltyping)){
            names(x = res2$cell.celltyping) <-
              paste0(names(x = .lm.obj$cells)[cells.idx],
                     "_",
                     names(x = res2$cell.celltyping))
          }
          
          if(!is.null(x = res2$lm.celltyping)){
            names(x = res2$lm.celltyping) <-
              paste0(names(x = .lm.obj$cells)[cells.idx],
                     "_",
                     names(x = res2$lm.celltyping))
          }
          
        }
        
        if(isTRUE(x = .verbose)){
          message(paste0("mapping progress: ",
                         round(x = cells.idx * 100 / length(x = .lm.obj$cells),
                               digits = 2),
                         "%"))
        }
        
        return(res2)
        
      })
    
    if(!is.null(x = .lm.obj$symphony.obj)){
      
      .lm.obj$graph$celltyping <-
        vector(mode = "list",
               length = 0)
      
      .lm.obj$graph$celltyping$ids <-
        seq_along(along.with = res) |>
        stats::setNames(nm = names(x = res)) |>
        lapply(FUN = function(smpl.idx){
          res[[smpl.idx]]$lm.celltyping
        }) |>
        unlist(use.names = FALSE)
      
      names(x = .lm.obj$graph$celltyping$ids) <-
        seq_along(along.with = res) |>
        stats::setNames(nm = names(x = res)) |>
        lapply(FUN = function(smpl.idx){
          names(x = res[[smpl.idx]]$lm.celltyping)
        }) |>
        unlist(use.names = FALSE)
      
      .lm.obj$graph$celltyping$ids <-
        .lm.obj$graph$celltyping$ids[
          rownames(x = .lm.obj$lm)
        ] |>
        as.factor()
      
      top <-
        apply(X = .lm.obj$pca$rotation,
              MARGIN = 2,
              FUN = function(PC.rot){
                
                order(PC.rot,
                      decreasing = TRUE) |>
                  (\(x)
                   c(utils::head(x = x, 
                                 n = 3),
                     utils::tail(x = x,
                                 n = 3))
                  )()
                
              }) |>
        as.vector() |> 
        (\(x)
         rownames(x = .lm.obj$pca$rotation)[x]
        )() |>
        unique()
      
      .lm.obj$graph$celltyping$mean.exprs <-
        (if(.lm.obj$assay.type == "RNA") .lm.obj$scaled.lm[,top] else .lm.obj$lm) |>
        dplyr::as_tibble() |>
        cbind(cell.pop = as.character(x = .lm.obj$graph$celltyping$ids)) |>
        dplyr::group_by(cell.pop) |>
        dplyr::summarize_all(.funs = mean) |>
        as.data.frame() |>
        (\(x)
         `rownames<-`(x = x[,colnames(x = x) != "cell.pop"],
                      value = x$cell.pop)
        )() |>
        as.matrix()
      
      .lm.obj$graph$celltyping$pheatmap <-
        pheatmap::pheatmap(mat = .lm.obj$graph$celltyping$mean.exprs,
                           color = grDevices::colorRampPalette(
                             unname(obj =
                                      Color.Palette[1,c(1,6,2)]))(100),
                           kmeans_k = NA,
                           breaks = NA,
                           border_color = NA,
                           scale = "none",
                           angle_col = 90,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           cellwidth = 20,
                           cellheight = 20,
                           treeheight_row = 20,
                           treeheight_col = 20,
                           silent = TRUE)
      
    }
    
    fdens <-
      seq_along(along.with = res) |>
      stats::setNames(nm = names(x = res)) |>
      lapply(FUN = function(smpl.idx){
        
        if(ncol(x = res[[smpl.idx]]$fgraph) == nrow(x = .lm.obj$lm)){ # NOT IDEAL! This if/else usage was only introduced to avoid bug in uwot reported in https://github.com/jlmelville/uwot/issues/129
          Matrix::colSums(x = res[[smpl.idx]]$fgraph[
            if(length(x = .cl.to.keep) !=
               (levels(x = .lm.obj$graph$clustering$ids) |>
                length())){
              (res[[smpl.idx]]$cell.clustering %in% .cl.to.keep)
            } else {
              if(!is.null(x = .ct.to.keep)){
                (res[[smpl.idx]]$cell.celltyping %in% .ct.to.keep)
              } else {
                1:nrow(x = res[[smpl.idx]]$fgraph)
              }
            },])
        } else {
          Matrix::rowSums(x = res[[smpl.idx]]$fgraph[,
                                                     if(length(x = .cl.to.keep) !=
                                                        (levels(x = .lm.obj$graph$clustering$ids) |>
                                                         length())){
                                                       (res[[smpl.idx]]$cell.clustering %in% .cl.to.keep)
                                                     } else {
                                                       if(!is.null(x = .ct.to.keep)){
                                                         (res[[smpl.idx]]$cell.celltyping %in% .ct.to.keep)
                                                       } else {
                                                         1:nrow(x = res[[smpl.idx]]$fgraph)
                                                       }
                                                     }])
        }
        
      }) |>
      do.call(what = cbind) |>
      (\(x)
       `rownames<-`(x = x,
                    value = rownames(x = .lm.obj$lm))
      )() |>
      Matrix::t() |>
      (\(x)
       seq_along(along.with = res) |>
         stats::setNames(nm = names(x = res)) |>
         lapply(FUN = function(smpl.idx){
           nrow(x = res[[smpl.idx]]$embedding)
         }) |>
         unlist() |>
         (\(n.cells)
          x / (n.cells / mean(x = n.cells))
         )()
      )() |>
      Matrix::t()
    
    cell.nlmn <-
      lapply(X = res,
             FUN = function(smpl){
               smpl$nn$euclidean$idx
             })
    
    cell.clustering <-
      lapply(X = res,
             FUN = function(smpl){
               smpl$cell.clustering
             })
    
    if(!is.null(x = .lm.obj$graph$celltyping)){
      
      cell.celltyping <-
        lapply(X = res,
               FUN = function(smpl){
                 smpl$cell.celltyping
               })
      
    } else {
      cell.celltyping <- NULL
    }
    
    .lm.obj$map <-
      list(fdens = fdens,
           clustering = list(ids = cell.clustering),
           celltyping = list(ids = cell.celltyping),
           nearest.landmarks = cell.nlmn)
    
    .lm.obj$map$clustering$cell.count <-
      seq_along(along.with = .lm.obj$map$clustering$ids) |>
      stats::setNames(nm = names(x = .lm.obj$map$clustering$ids)) |>
      lapply(FUN = function(smpl.idx){
        
        table(.lm.obj$map$clustering$ids[[smpl.idx]]) |>
          (\(x)
           stats::setNames(object = as.vector(x = x),
                           nm = names(x = x))
          )()
      }) |>
      dplyr::bind_rows(.id = "sample") |>
      as.data.frame() |>
      (\(x)
       `rownames<-`(x = as.matrix(x = x[,colnames(x = x) != "sample"]),
                    value = x$sample)
      )() |>
      (\(x)
       x[,colnames(x = x) |>
           sort()]
      )() |> 
      (\(x)
       x[,!(colnames(x = x) %in% .cl.ct.to.ign)]
      )()
    
    .lm.obj$map$clustering$cell.count[
      is.na(x = .lm.obj$map$clustering$cell.count)
    ] <- 0
    
    .lm.obj$map$clustering$cell.perc <-
      (.lm.obj$map$clustering$cell.count * 100) /
      Matrix::rowSums(x = .lm.obj$map$clustering$cell.count)
    
    if(!is.null(x = .lm.obj$graph$celltyping)){
      
      .lm.obj$map$celltyping$cell.count <-
        seq_along(along.with = .lm.obj$map$celltyping$ids) |>
        stats::setNames(nm = names(x = .lm.obj$map$celltyping$ids)) |>
        lapply(X = ,
               FUN = function(smpl.idx){
                 
                 table(.lm.obj$map$celltyping$ids[[smpl.idx]]) |>
                   (\(x)
                    stats::setNames(object = as.vector(x = x),
                                    nm = names(x = x))
                   )()
                 
               }) |>
        dplyr::bind_rows(.id = "sample") |>
        as.data.frame() |>
        (\(x)
         `rownames<-`(x = as.matrix(x = x[,colnames(x = x) != "sample"]),
                      value = x$sample)
        )()  |>
        (\(x)
         x[,colnames(x = x) |>
             sort()]
        )() |> 
        (\(x)
         x[,!(colnames(x = x) %in% .cl.ct.to.ign)]
        )()
      
      .lm.obj$map$celltyping$cell.count[
        is.na(x = .lm.obj$map$celltyping$cell.count)
      ] <- 0
      
      .lm.obj$map$celltyping$cell.perc <-
        (.lm.obj$map$celltyping$cell.count * 100) /
        Matrix::rowSums(x = .lm.obj$map$celltyping$cell.count)
      
    }
    
    .lm.obj$map$cl.ct.to.ign <-
      .cl.ct.to.ign
    
    return(.lm.obj)
    
  }

#' Leiden clustering of landmarks
#'
#' Applies Leiden community detection to partition landmarks into clusters based on either 
#' the SNN graph or UMAP fuzzy graph. Clusters are used for visualization and as a coarser 
#' grouping for statistical testing. This is a wrapper around \code{leiden.cluster}.
#'
#' @param .lm.obj A tinydenseR object with \code{$graph} component from \code{get.graph}.
#' @param .cl.method Character specifying clustering method: "snn" (shared nearest neighbors) 
#'   or "fgraph" (fuzzy graph from UMAP). Default "snn".
#' @param .cl.resolution.parameter Numeric resolution for Leiden clustering (default 0.8). 
#'   Internally scaled by 1e-3 for CPM objective. Higher values yield finer clusters.
#' @param .seed Integer seed for reproducibility (default 123).
#' @param .verbose Logical for progress messages (default TRUE).
#' @param .small.size Integer threshold for straggler absorption (default 3). Clusters smaller 
#'   than this are merged into most connected neighbors.
#'   
#' @return Updated \code{.lm.obj} with \code{$graph$clustering} containing:
#'   \itemize{
#'     \item \code{ids}: Factor of cluster assignments for each landmark
#'     \item \code{mean.exprs}: Matrix of mean expression per cluster (top PCs for RNA, all markers for cyto)
#'     \item \code{pheatmap}: Heatmap object showing cluster profiles
#'   }
#'   
#' @details
#' For RNA data, the heatmap displays top contributors to each principal component (6 genes 
#' per PC: 3 most positive, 3 most negative). For cytometry, all markers are shown.
#' 
#' The fuzzy graph method ("fgraph") uses probabilistic UMAP edges and may produce more 
#' visually coherent clusters in UMAP space. The SNN method uses Jaccard similarity and 
#' often gives more stable results for downstream statistics.
#' 
#' @seealso \code{\link{leiden.cluster}}, \code{\link{get.graph}}
#' 
#' @examples
#' \dontrun{
#' # After graph construction
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
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
#' @export
lm.cluster <-
  function(.lm.obj,
           .cl.method = "snn",
           .cl.resolution.parameter = 0.8,
           .seed = 123,
           .verbose = TRUE,
           .small.size = 3){
    
    # R CMD check appeasement
    cell.pop <- NULL
    
    set.seed(seed = .seed)
    
    if(is.null(.lm.obj$graph)){
      stop("Graph component missing. Run get.graph() before clustering.")
    }
    
    # Validate clustering method
    .cl.method <-
      match.arg(arg = .cl.method,
                choices = c("fgraph","snn"))
    
    .lm.obj$graph$clustering <-
      vector(mode = "list")
    
    # Run Leiden clustering on selected similarity matrix
    .lm.obj$graph$clustering$ids <-
      leiden.cluster(
        .lm.obj = .lm.obj,
        .sim.matrix = if(.cl.method == "fgraph"){
          Matrix::drop0(x = .lm.obj$graph$uwot$fgraph,   # Prune weak fuzzy edges
                        tol = 1/20)
        } else {
          .lm.obj$graph$snn                              # Use SNN graph
        },
        .resolution.parameter = .cl.resolution.parameter * 1e-3,  # Scale for CPM
        .small.size = .small.size,
        .verbose = .verbose)
    
    if(isTRUE(x = .verbose)){
      message("cluster heatmap")
    }
    
    # For RNA: select top genes from PCA rotation (top 3 and bottom 3 per PC)
    if(.lm.obj$assay.type == "RNA"){
      
      top <-
        apply(X = .lm.obj$pca$rotation,
              MARGIN = 2,
              FUN = function(PC.rot){
                
                order(PC.rot,
                      decreasing = TRUE) |>
                  (\(x)
                   c(utils::head(x = x, 
                                 n = 3),
                     utils::tail(x = x,
                                 n = 3))
                  )()
                
              }) |>
        as.vector() |> 
        (\(x)
         rownames(x = .lm.obj$pca$rotation)[x]
        )() |>
        unique()
      
    }
    
    # Compute mean expression per cluster
    .lm.obj$graph$clustering$mean.exprs <-
      (if(.lm.obj$assay.type == "RNA") .lm.obj$scaled.lm[,top] else .lm.obj$lm) |>
      dplyr::as_tibble() |>
      cbind(cell.pop = as.character(x = .lm.obj$graph$clustering$ids)) |>
      dplyr::group_by(cell.pop) |>
      dplyr::summarize_all(.funs = mean) |>
      as.data.frame() |>
      (\(x)
       `rownames<-`(x = x[,colnames(x = x) != "cell.pop"],
                    value = x$cell.pop)
      )() |>
      as.matrix()
    
    # Generate heatmap of cluster profiles
    .lm.obj$graph$clustering$pheatmap <-
      pheatmap::pheatmap(mat = .lm.obj$graph$clustering$mean.exprs,
                         color = grDevices::colorRampPalette(
                           unname(obj =
                                    Color.Palette[1,c(1,6,2)]))(100),
                         kmeans_k = NA,
                         breaks = NA,
                         border_color = NA,
                         scale = "none",
                         angle_col = 90,
                         cluster_rows = if(nrow(x = .lm.obj$graph$clustering$mean.exprs) > 1) TRUE else FALSE,
                         cluster_cols = TRUE,
                         cellwidth = 20,
                         cellheight = 20,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         silent = TRUE)
    
    return(.lm.obj)
  }

#' Graph-based feature discovery for landmarks
#'
#' Identifies the most characteristic features (genes or markers) for each landmark by 
#' determining which features contribute most strongly to the landmark's defining principal 
#' component. This creates landmark-specific feature signatures useful for interpretation.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.landmarks} (PCA required).
#'   
#' @return Updated \code{.lm.obj} with \code{$interact.plot$lm.features} containing:
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
#' Results are stored for use in \code{interactFeatPlot}, which displays feature signatures 
#' when hovering over landmarks in UMAP plots.
#' 
#' @seealso \code{\link{interactFeatPlot}}, \code{\link{get.landmarks}}
#' 
#' @examples
#' \dontrun{
#' # After landmark selection and graph construction
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.lm.features.stats()
#' 
#' # View feature signature for first landmark
#' head(lm.cells$interact.plot$lm.features$res[[1]])
#' 
#' # Use in interactive plotting
#' plotPCA(lm.cells, .hover.stats = "marker")
#' }
#' 
#' @export
get.lm.features.stats <-
  function(
    .lm.obj
  ){
    
    # Identify the top (dominant) PC for each landmark
    top.comp <-
      .lm.obj$pca$embed |>
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
      diag(x = .lm.obj$pca$embed[,top.comp]) |>
      sign()
    
    # Signed feature loadings: rotation × sign for each landmark's top PC
    coefs <-
      ((Matrix::t(x = .lm.obj$pca$rotation[,top.comp]) * top.comp.sign) |>
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
      nrow(x = .lm.obj$lm) |>
      seq_len() |>
      stats::setNames(nm = rownames(x = .lm.obj$lm)) |>
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
    
    .lm.obj$interact.plot$lm.features <-
      list(res = res,
           html = html)
    
    return(.lm.obj)
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


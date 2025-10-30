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

#' Leiden clustering
#'
#' This function performs Leiden clustering on a graph of cells. The graph is represented by an adjacency matrix, where the elements of the matrix are the likelihood that pairs of vertices are adjacent. The function returns a factor vector of cluster assignments.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks. If provided, the function will use the Laplacian Eigenmap embedding from the .lm.obj for initialization. If NULL, the function will compute a 3-dimensional embedding using irlba on the similarity matrix.
#' @param .sim.matrix A square matrix representing the graph of the cells. Elements of the matrix should be a similarity measure of pairs of observations.
#' @param .resolution.parameter A numeric value that determines the granularity of the clustering. Higher values result in more clusters.
#' @param .small.size An integer value that determines the size of clusters that are considered "stragglers". Stragglers are assigned to the cluster they are most connected to. Defaults to 0.5% of the number of cells.
#' @param .verbose A logical value indicating whether to print messages about the clustering process.
#' @param .seed An integer value that sets the seed for reproducibility.
#' @return A factor vector of cluster assignments.
leiden.cluster <-
  function(.lm.obj = NULL,
           .sim.matrix,
           .resolution.parameter,
           .small.size = floor(x = nrow(x = .sim.matrix) / 200),
           .verbose = TRUE,
           .seed = 123) {
    
    .init.embed <- 
      if(!is.null(x = .lm.obj)){
        .lm.obj$graph$uwot$embedding
      } else {
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
      stop("similarity matrix must be square")
    }
    
    Matrix::diag(x = .sim.matrix) <-
      0
    
    g <-
      igraph::graph_from_adjacency_matrix(
        adjmatrix = .sim.matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE,
        add.colnames = NULL,
        add.rownames = NA)
    
    # try to avoid singletons by starting with simple kmeans clustering instead of singletons
    kres <-
      stats::kmeans(x = .init.embed,
                    centers = min(25,
                                  (nrow(x = .init.embed) / 2) |>
                                    ceiling()),
                    nstart = 10,
                    iter.max = 100)
    
    #https://stackoverflow.com/a/30055776
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
    
    if (isTRUE(.verbose)) {
      message("checking for stragglers.")
    }
    
    # clusters smaller than .small.size are candidates to absorb
    tbl <- table(ids)
    stragglers <- names(tbl)[tbl < .small.size]
    keepers    <- setdiff(names(tbl), stragglers)
    
    if (length(stragglers) > 0 && length(keepers) > 0) {
      
      # Ensure row/colnames on .sim.matrix match cell order used in 'ids'
      if (is.null(rownames(.sim.matrix))) {
        # If you have real cell IDs, use those instead of seq_len
        rownames(.sim.matrix) <- colnames(.sim.matrix) <- seq_len(nrow(.sim.matrix))
      }
      
      for (i in stragglers) {
        i.cells <- which(ids == i)
        # Build connectivity afresh for this straggler cluster
        connectivity <- numeric(length(keepers))
        names(connectivity) <- keepers
        
        for (j in keepers) {
          j.cells <- which(ids == j)
          sub.sim.matrix <- .sim.matrix[i.cells, j.cells, drop = FALSE]
          
          if (inherits(sub.sim.matrix, "sparseMatrix")) {
            # mean of block (avoid converting full dense)
            connectivity[j] <- sum(sub.sim.matrix@x) / (nrow(sub.sim.matrix) * ncol(sub.sim.matrix))
          } else {
            connectivity[j] <- mean(sub.sim.matrix)
          }
        }
        
        # Tie-handling: pick first max, or sample among maxima if you prefer
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
    
    # relabeling after merge
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
#' This function converts a matrix of nearest neighbors into a non-symetric sparse matrix representation.
#' 
#' @param .nn.idx A matrix where each row represents a cell and each column represents a nearest neighbor. The elements of the matrix are the indices of the nearest neighbors.
#' @return A non-symetric sparse matrix representation of the nearest neighbors.
#'
get.adj.matrix <-
  function(.nn.idx){
    Matrix::sparseMatrix(i = as.vector(x = .nn.idx),
                         j = nrow(x = .nn.idx) |>
                           seq_len() |>
                           rep(times = ncol(x = .nn.idx)),
                         x = 1,
                         dims = nrow(x = .nn.idx) |>
                           rep(times = 2),
                         repr = "C", 
                         dimnames = list(rownames(x = .nn.idx),
                                         rownames(x = .nn.idx)))
  }

#' Shared nearest neighbors via fast Jaccard index calculation
#'
#' This function calculates the Jaccard index between cells based on their nearest neighbors. The function is optimized for speed and memory efficiency.
#'
#' @param .adj.matrix A non-symetric sparse matrix representation of the nearest neighbors.
#' @param .prune A numeric value that determines the threshold for pruning small values from the Jaccard index matrix.
#' @return A shared nearest neighbors sparse adjacency matrix of Jaccard indices.
fast.jaccard.r <-
  function(
    .adj.matrix,
    .prune = 1/15
  ){
    
    adj.cp <-
      Matrix::crossprod(x = .adj.matrix)
    
    adj.cp@x <-
      1 / (((sum(x = .adj.matrix[,1]) * 2) / adj.cp@x) - 1)
    
    adj.cp <-
      Matrix::drop0(x = adj.cp,
                    tol = .prune) |>
      Matrix::forceSymmetric() |>
      methods::as(Class = "generalMatrix")
    
    return(adj.cp)
    
  }

#' Graph embedding of landmarks
#'
#' This function performs graph-based clustering, dimensionality reduction, and map building of landmark
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks.
#' @param .k An integer value that determines the number of nearest neighbors to consider in the graph.
#' @param .scale A logical value indicating whether to scale the data before embedding. Default is FALSE for assay type RNA (to avoid scaling PCA results) and TRUE for cyto.
#' @param .verbose A logical value indicating whether to print messages about the embedding process.
#' @param .seed An integer value that sets the seed for reproducibility.
#' @param .cl.method A character value that determines the clustering method. Options are "snn" for shared nearest neighbors and "fgraph" for fuzzy graph.
#' @param .cl.resolution.parameter A numeric value that determines the granularity of the clustering. Higher values result in more clusters.
#' @param .small.size An integer value that determines the size of clusters that are considered "stragglers". Stragglers are assigned to the cluster they are most connected to. Defaults to 0.5% of the number of cells.
#' @return The .lm.obj now containing a graph object with, for example, the UMAP embedding, cluster assignments, mean cluster expression values as well as the mean expression heatmap.
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
    
    .cl.method <-
      match.arg(arg = .cl.method,
                choices = c("fgraph","snn"))
    
    set.seed(seed = .seed)
    .lm.obj$graph$uwot <-
      uwot::umap(X = if(!is.null(x = .lm.obj$harmony.obj) || .lm.obj$assay.type == "RNA") .lm.obj$pca$embed else .lm.obj$lm,
                 n_neighbors = .k,
                 n_components = 2,
                 n_epochs = 500,
                 scale = .scale,
                 pca = NULL,
                 verbose = isTRUE(x = .verbose),
                 ret_model = TRUE,
                 batch = TRUE,
                 seed = .seed,
                 n_threads = .lm.obj$n.threads,
                 fast_sgd = FALSE,
                 n_sgd_threads = .lm.obj$n.threads,
                 ret_extra = c("fgraph",
                               "nn"))
    
    colnames(x = .lm.obj$graph$uwot$embedding) <-
      c("umap.1",
        "umap.2")
    
    .lm.obj$graph$adj.matrix <-
      get.adj.matrix(.nn.idx = .lm.obj$graph$uwot$nn$euclidean$idx)
    
    if(isTRUE(x = .verbose)){
      message("getting snn")
    }
    
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
    
    # see https://github.com/jlmelville/uwot/blob/f9e576e97d9df44d48be2cc559412282838dc4a5/R/init.R
    .lm.obj$graph$LE$W.sym <-
      (((.lm.obj$graph$adj.matrix > 0) |
          (Matrix::t(x = .lm.obj$graph$adj.matrix) > 0)) * 1) |>
      Matrix::Matrix(sparse = TRUE)
    
    # Compute as many spectrum components as PCA dims
    target_k <-
      ncol(x = .lm.obj$pca$embed)
    nv <-
      min(target_k, 
          nrow(x = .lm.obj$pca$embed) - 1L) # cannot exceed n-1
    stopifnot(nv >= 2L)
    
    .lm.obj$graph$LE <- 
      c(.lm.obj$graph$LE,
        form_modified_laplacian(A = .lm.obj$graph$LE$W.sym,
                                ret_d = TRUE))
    
    .lm.obj$graph$LE <-
      c(.lm.obj$graph$LE,
        irlba_spectral_tsvd(L = .lm.obj$graph$LE$L, 
                            n = nv))
    
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
      
      # Drop ALL near-zero eigenvalues
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
        
        if(length(x = .lm.obj$graph$LE$nontriv) < 4L){
          
          k <- length(x = .lm.obj$graph$LE$nontriv)
          
        } else {
          
          .lm.obj$graph$LE$elbow <-
            elbow.sec.deriv(x = .lm.obj$graph$LE$values[.lm.obj$graph$LE$nontriv],
                            smooth = TRUE,
                            df = NULL,
                            sort.order = "asc")$index
          
          # Bound k by available non-trivial eigenpairs and target_k
          k <- 
            min(.lm.obj$graph$LE$elbow, 
                target_k, 
                length(x = .lm.obj$graph$LE$nontriv)) |>
            max(2L)
          
        }
        
        # Take the first k non-trivial eigenvectors (ascending spectrum)
        .lm.obj$graph$LE$keep <-
          .lm.obj$graph$LE$nontriv[seq_len(length.out = k)]
        
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
#' This function maps cells to landmarks to derive the probability that each cell is in the neighborhood of a landmark.
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks and get.graph.
#' @param .ref.obj A symphony object with Z_corr field populated.
#' @param .integrate.vars A character vector of batch variables to integrate (column names in metadata).
#' @param .celltype.col.name The column name in the reference object metadata that contains cell type information.
#' @param .cl.ct.to.ign A character value that determines the cluster or cell type to ignore during mapping. Default is NULL. Use this to ignore clusters or cell types that are not relevant for the statistical analysis, such as erythrocytes in PBMC or BMMC, for example, wheresince their presence in such samples is more likely to be associated with technical artifacts during sample processing in the lab rather than biological phenotypes.
#' @param .irrel.par A character vector of irrelevant parameters to exclude from the mapping. Default is NULL.
#' @param .verbose A logical value indicating whether to print messages about the mapping process.
#' @param .seed An integer value that sets the seed for reproducibility.
#' @return The .lm.obj now containing a list that includes the sample-wise sum of cell-landmark neighboring probabilities, celltypes (if available, see "celltyping") and cluster assignments for each cell, and nearest neighboring landmarks for each cell.
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
    
    cell.pop <- id <- value <- ri <- NULL
    
    if(!is.null(x = .ref.obj)){
      
      if(.lm.obj$assay.type != "RNA"){
        stop("celltyping with reference object is only supported for assay type RNA")
      }
      
      if(inherits(x = .ref.obj,
                  what = "list")){
        
        if(is.null(x = .ref.obj$Z_corr)){
          stop(".ref.obj must be a symphony object with Z_corr field populated")
          
        }
        
      } else {
        
        stop(".ref.obj must be a symphony object with Z_corr field populated")
        
      }
      
      if(!(.celltype.col.name %in% colnames(x = .ref.obj$meta_data))){
        
        stop(paste0("Column name '", .celltype.col.name,
                    "' not found in .ref.obj$meta_data."))
        
      }
      
      warning("Using symphony reference object for cell typing and removing previous celltyping from .lm.obj")
      
      if(isTRUE(x = .verbose)){
        message("building annoy index for symphony reference object")
      }
      
      .lm.obj$symphony.obj <-
        .ref.obj
      
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
    
    if(!is.null(x = .cl.ct.to.ign)){
      
      if(length(x = .cl.ct.to.ign) > 1){
        stop("Please provide only one cluster or cell type to ignore. If more than one cluster or cell type should be ignored, please merge them manually using the `celltyping` function.")
      }
      
      if(!(.cl.ct.to.ign %in% levels(x = .lm.obj$graph$clustering$ids)) &
         !((.cl.ct.to.ign %in% levels(x = .lm.obj$graph$celltyping$ids)) |
           ((.cl.ct.to.ign %in% .ref.obj$meta_data[[.celltype.col.name]])))){
        stop("the cluster or cell type to ignore is not present in the clustering, celltyping nor .ref.obj metadata.")
      }
      
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
#' This function performs Leiden clustering on the landmarks.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj.
#' @param .cl.method A character value that determines the clustering method. Options are "snn" for shared nearest neighbors and "fgraph" for fuzzy graph.
#' @param .cl.resolution.parameter A numeric value that determines the granularity of the clustering. Higher values result in more clusters.
#' @param .seed An integer value that sets the seed for reproducibility.
#' @param .verbose A logical value indicating whether to print messages about the clustering process.
#' @param .small.size An integer value that determines the size of clusters that are considered "stragglers". Stragglers are assigned to the cluster they are most connected to.
#' @return The .lm.obj now containing a list that includes the cluster assignments for each landmark, mean cluster expression values as well as the mean expression heatmap.
#' @export
lm.cluster <-
  function(.lm.obj,
           .cl.method = "snn",
           .cl.resolution.parameter = 0.8,
           .seed = 123,
           .verbose = TRUE,
           .small.size = 3){
    
    cell.pop <- NULL
    
    set.seed(seed = .seed)
    
    if(is.null(.lm.obj$graph)){
      stop("please run get.graph() first")
    }
    
    .cl.method <-
      match.arg(arg = .cl.method,
                choices = c("fgraph","snn"))
    
    .lm.obj$graph$clustering <-
      vector(mode = "list")
    
    .lm.obj$graph$clustering$ids <-
      leiden.cluster(
        .lm.obj = .lm.obj,
        .sim.matrix = if(.cl.method == "fgraph"){
          Matrix::drop0(x = .lm.obj$graph$uwot$fgraph,
                        tol = 1/20)
        } else {
          .lm.obj$graph$snn
        },
        .resolution.parameter = .cl.resolution.parameter * 1e-3,
        .small.size = .small.size,
        .verbose = .verbose)
    
    if(isTRUE(x = .verbose)){
      message("cluster heatmap")
    }
    
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
#' This function identifies features that are most characteristic for each landmark based on the PCA rotation matrix.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks.
#' @return A list that includes the feature discovery results for each landmark.
#' @export
get.lm.features.stats <-
  function(
    .lm.obj
  ){
    
    top.comp <-
      .lm.obj$pca$embed |>
      abs() |>
      as.matrix() |>
      (\(x)
       matrixStats::rowRanks(x = x,
                             ties.method = "max") == # using max here to ensure that, if the top comp is a tie, the highest rank is assigned and corresponds to the ncol of the matrix
         ncol(x = x)
      )() |>
      which(arr.ind = TRUE) |>
      (\(x)
       x[order(x = x[,"row"]),"col"]
      )()
    
    top.comp.sign <-
      diag(x = .lm.obj$pca$embed[,top.comp]) |>
      sign()
    
    coefs <-
      ((Matrix::t(x = .lm.obj$pca$rotation[,top.comp]) * top.comp.sign) |>
         Matrix::t())
    
    hits <-
      coefs |>
      (\(x)
       matrixStats::colRanks(x = x,
                             ties.method = "max",
                             preserveShape = TRUE) > # using max here to ensure that, the top comp is a tie, the highest rank is assigned and corresponds to the ncol of the matrix
         (nrow(x = x) - 10)
      )() |>
      which(arr.ind = TRUE) |>
      as.data.frame()
    
    res <-
      nrow(x = .lm.obj$lm) |>
      seq_len() |>
      stats::setNames(nm = rownames(x = .lm.obj$lm)) |>
      lapply(FUN = function(lm.idx){
        
        coefs[hits$row[hits$col == lm.idx],
              lm.idx] |>
          sort(decreasing = TRUE) |>
          as.data.frame() |>
          (\(x)
           `colnames<-`(x = x,
                        value = "feature importance")
          )()
        
      })
    
    html <-
      lapply(X = res,
             FUN = knitr::kable,
             format = "html")
    
    .lm.obj$interact.plot$lm.features <-
      list(res = res,
           html = html)
    
    return(.lm.obj)
  }

annoy_build <-
  getFromNamespace(x = "annoy_build",
                   ns = "uwot")

annoy_search <-
  getFromNamespace(x = "annoy_search",
                   ns = "uwot")

annoy_nn <-
  getFromNamespace(x = "annoy_nn",
                   ns = "uwot")

form_modified_laplacian <-
  getFromNamespace(x = "form_modified_laplacian",
                   ns = "uwot")

irlba_spectral_tsvd <-
  getFromNamespace(x = "irlba_spectral_tsvd",
                   ns = "uwot")


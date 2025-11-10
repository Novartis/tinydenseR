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

#' Summarize gene expression patterns for genes of interest
#' 
#' Generates summaries of gene expression for specified genes of interest (GOI),
#' including detection status (positive/negative), cell counts, and percentages
#' across samples. Useful for identifying marker gene expression patterns across
#' clusters or cell types.
#' 
#' @details
#' For each gene of interest, the function:
#' \enumerate{
#'   \item Determines which cells express the gene (pos.) vs don't express it (neg.)
#'   \item Creates three parallel summary structures:
#'     \itemize{
#'       \item \code{clustering} - summaries by cluster IDs with pos./neg. prefix
#'       \item \code{celltyping} - summaries by cell type labels with pos./neg. prefix (if available)
#'       \item \code{all} - summaries by pos./neg. status only
#'     }
#'   \item Computes cell counts and percentages for each category per sample
#' }
#' 
#' The function operates on mapped cells (after \code{get.map}), allowing you to
#' optionally filter to specific clusters or cell types using \code{.id} parameter.
#' 
#' @param .lm.obj A list object initialized with \code{setup.lm.obj} and processed 
#'   with \code{get.graph} and \code{get.map}. Must contain RNA assay data.
#' @param .goi Character vector of gene names to summarize. All genes must exist 
#'   in \code{.lm.obj$raw.lm} column names.
#' @param .id.idx Integer index of a specific landmark to analyze. If NULL (default), 
#'   uses \code{.id} parameter or all cells. Rarely used directly.
#' @param .id Character vector of cluster or cell type IDs to filter analysis. 
#'   If NULL (default), includes all clusters/cell types. Use this to focus on 
#'   specific populations (e.g., \code{c("cluster.01", "cluster.02")}).
#' @param .id.from Character specifying whether \code{.id} refers to "clustering" 
#'   or "celltyping" identifiers. Default is "clustering".
#' @param .verbose Logical, whether to print progress messages. Default TRUE.
#' 
#' @return A named list with one element per gene in \code{.goi}. Each gene's 
#'   element contains:
#'   \describe{
#'     \item{\code{$clustering}}{List with \code{$ids}, \code{$cell.count}, and 
#'       \code{$cell.perc} for cluster-level summaries with pos./neg. prefixes}
#'     \item{\code{$celltyping}}{Same structure as clustering but for cell types 
#'       (NULL if cell typing not performed)}
#'     \item{\code{$all}}{Same structure but only pos./neg. categories without 
#'       cluster/cell type breakdown}
#'   }
#'   
#'   Cell counts are matrices with samples as rows and categories as columns.
#'   Cell percentages show within-sample proportions.
#' 
#' @examples
#' \dontrun{
#' # Complete workflow with gene expression summary
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.map()
#' 
#' # Summarize expression of marker genes
#' markers <- c("CD4", "CD8A", "CD19")
#' goi_results <- goi.summary(lm.cells, markers)
#' 
#' # Check CD4 expression by cluster
#' goi_results$CD4$clustering$cell.perc
#' 
#' # Focus analysis on specific clusters only
#' goi_results <- goi.summary(lm.cells, markers, 
#'                            .id = c("cluster.01", "cluster.02"))
#' }
#' @export
goi.summary <-
  function(
    .lm.obj,
    .goi,
    .id.idx = NULL,
    .id = NULL,
    .id.from = "clustering",
    .verbose = TRUE
  ){
    
    if(!all(.goi %in% colnames(x = .lm.obj$raw.lm))){
      stop("Gene(s) not found in data: ",
           paste(.goi[!(.goi %in% colnames(x = .lm.obj$raw.lm))],
                 collapse = ", "),
           "\nCheck gene names (case-sensitive) or use colnames(.lm.obj$raw.lm) to see available genes.")
    }
    
    if(.lm.obj$assay.type != "RNA"){
      stop("goi.summary() only supports RNA assay data.\n",
           "Current assay type: ", .lm.obj$assay.type)
    }
    
    if(is.null(x = .lm.obj$map)){
      stop("Cell mapping not found. Run get.graph() and get.map() before calling goi.summary().")
    }
    
    goi <-
      vector(mode = "list",
             length = 0)
    
    # Determine which cell indices to analyze
    # Three modes: (1) all cells, (2) cells in specific clusters/celltypes, (3) specific landmark index
    if(is.null(x = .id.idx)){
      
      if(!is.null(x = .id)){
        
        .id.from <-
          match.arg(arg = .id.from,
                    choices = c("clustering",
                                "celltyping"))
        
        if(!all(.id %in% unique(x = .lm.obj$graph[[.id.from]]$ids))){
          
          stop("Invalid ", .id.from, " ID(s): ",
               paste(.id[!(.id %in% unique(x = .lm.obj$graph[[.id.from]]$ids))],
                     collapse = ", "),
               "\nThese IDs were not found. Use unique(.lm.obj$graph$", .id.from, "$ids) to see valid IDs.")
          
        }
        
        # Filter to cells belonging to specified cluster/celltype IDs
        .id.idx <-
          lapply(X = .lm.obj$map[[.id.from]]$ids,
                 FUN = function(smpl){
                   which(x = smpl %in% .id)
                 })
        
      } else {
        
        # Use all cells (default behavior)
        .id.idx <-
          lapply(X = .lm.obj$map[[.id.from]]$ids,
                 FUN = seq_along)
        
      }
    } else {
      
      # Filter to cells mapped to a specific landmark index
      .id.idx <-
        lapply(X = .lm.obj$map$nearest.landmarks,
               FUN = function(smpl.knn)
                 which(x = smpl.knn[,1] == .id.idx))
      
    }
    
    # Process each sample: extract expression, label cells as pos./neg. for each gene
    goi <-
      seq_along(along.with = .lm.obj$cells) |>
      stats::setNames(nm = names(x = .lm.obj$cells)) |>
      lapply(FUN = function(cells.elem){
        
        # Load expression data for selected cells and genes
        goi.exprs.mat <-
          readRDS(file = .lm.obj$cells[[cells.elem]]) |> 
          (\(x)
           x[.goi,.id.idx[[cells.elem]],drop = FALSE]
          )() |>
          Matrix::t()
        
        # Determine which cells express each gene (expression > 0)
        all.goi.detected <-
          as.matrix(x = goi.exprs.mat > 0)
        
        # Create cluster ID matrix for all cells × genes combinations
        ids <-
          .lm.obj$map$clustering$ids[[cells.elem]][.id.idx[[cells.elem]]] |>
          rep(times = length(x = .goi)) |>
          matrix(byrow = FALSE,
                 nrow = length(x = .id.idx[[cells.elem]]),
                 ncol = length(x = .goi),
                 dimnames = list(names(x = .lm.obj$map$clustering$ids[[cells.elem]][.id.idx[[cells.elem]]]),
                                 .goi))
        
        # Add "pos." prefix for expressing cells, "neg." for non-expressing
        ids[all.goi.detected] <-
          paste0("pos.",
                 ids[all.goi.detected])
        
        ids[!all.goi.detected] <-
          paste0("neg.",
                 ids[!all.goi.detected])
        
        # Repeat same process for cell type labels (if available)
        if(!is.null(x = .lm.obj$map$celltyping)){
          
          ct.ids <-
            .lm.obj$map$celltyping$ids[[cells.elem]][.id.idx[[cells.elem]]] |>
            rep(times = length(x = .goi)) |>
            matrix(byrow = FALSE,
                   nrow = length(x = .id.idx[[cells.elem]]),
                   ncol = length(x = .goi),
                   dimnames = list(names(x = .lm.obj$map$celltyping$ids[[cells.elem]][.id.idx[[cells.elem]]]),
                                   .goi))
          
          ct.ids[all.goi.detected] <-
            paste0("pos.",
                   ct.ids[all.goi.detected])
          
          ct.ids[!all.goi.detected] <-
            paste0("neg.",
                   ct.ids[!all.goi.detected])
          
        } else {
          
          ct.ids <- NULL
          
        }
        
        if(isTRUE(x = .verbose)){
          message(
            paste0("progress: ",
                   round(x = cells.elem * 100 / length(x = .lm.obj$cells),
                         digits = 2),
                   "%")
          )
        }
        
        return(
          list(
            clustering = 
              as.data.frame(x = ids) |>
              as.list(),
            celltyping =
              as.data.frame(x = ct.ids) |>
              as.list(),
            # "all" category: strip pos./neg. prefix to get just expression status
            # Uses "." as delimiter (literal dot, not regex) to split "pos.cluster.01" -> "pos"
            all = 
              strsplit(x = ids,
                       split = ".", 
                       fixed = TRUE) |>
              lapply(FUN = utils::head,
                     n = 1) |>
              unlist(use.names = TRUE) |>
              matrix(nrow = length(x = .id.idx[[cells.elem]]),
                     ncol = length(x = .goi),
                     dimnames = list(names(x = .lm.obj$map$clustering$ids[[cells.elem]][.id.idx[[cells.elem]]]),
                                     .goi)) |>
              as.data.frame() |>
              as.list()
          ))
        
      }) |> 
      # Restructure: from sample-centric to gene-centric organization
      # Input: list by sample, each with clustering/celltyping/all for all genes
      # Output: list by gene, each with clustering/celltyping/all across all samples
      (\(x)
       stats::setNames(object = .goi,
                       nm = .goi) |>
         lapply(FUN = function(goi.elem) {
           list(clustering = list(ids = lapply(X = x,
                                               FUN = function(y){
                                                 y$clustering[[goi.elem]]
                                               })),
                celltyping = list(ids = lapply(X = x,
                                               FUN = function(y){
                                                 y$celltyping[[goi.elem]]
                                               })),
                all = list(ids = lapply(X = x,
                                        FUN = function(y){
                                          y$all[[goi.elem]]
                                        })))
         })   
      )()
    
    
    if(isTRUE(x = .verbose)){
      message("Retrieving cell counts and percentages.")
    }
    
    # For each gene, compute cell counts and percentages across samples
    goi <-
      lapply(X = goi, 
             FUN = function(goi.elem){
               
               # Count cells in each pos./neg. × cluster category per sample
               goi.elem$clustering$cell.count <-
                 seq_along(along.with = goi.elem$clustering$ids) |>
                 stats::setNames(nm = names(x = goi.elem$clustering$ids)) |>
                 lapply(FUN = function(smpl.idx){
                   
                   table(goi.elem$clustering$ids[[smpl.idx]]) |>
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
                 )()
               
               # Replace NA with 0 for categories not present in a sample
               goi.elem$clustering$cell.count[
                 is.na(x = goi.elem$clustering$cell.count)
               ] <- 0
               
               # Convert counts to within-sample percentages
               goi.elem$clustering$cell.perc <-
                 (goi.elem$clustering$cell.count * 100) /
                 Matrix::rowSums(x = goi.elem$clustering$cell.count)
               
               if(!is.null(x = .lm.obj$map$celltyping)){
                 
                 # Same process for cell type labels
                 goi.elem$celltyping$cell.count <-
                   seq_along(along.with = goi.elem$celltyping$ids) |>
                   stats::setNames(nm = names(x = goi.elem$celltyping$ids)) |>
                   lapply(FUN = function(smpl.idx){
                            
                            table(goi.elem$celltyping$ids[[smpl.idx]]) |>
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
                   )()
                 
                 goi.elem$celltyping$cell.count[
                   is.na(x = goi.elem$celltyping$cell.count)
                 ] <- 0
                 
                 goi.elem$celltyping$cell.perc <-
                   (goi.elem$celltyping$cell.count * 100) /
                   Matrix::rowSums(x = goi.elem$celltyping$cell.count)
                 
               }
               
               # Same process for "all" category (pos./neg. only, no cluster/celltype breakdown)
               goi.elem$all$cell.count <-
                 seq_along(along.with = goi.elem$all$ids) |>
                 stats::setNames(nm = names(x = goi.elem$all$ids)) |>
                 lapply(FUN = function(smpl.idx){
                   
                   table(goi.elem$all$ids[[smpl.idx]]) |>
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
                 )()
               
               goi.elem$all$cell.count[
                 is.na(x = goi.elem$all$cell.count)
               ] <- 0
               
               goi.elem$all$cell.perc <-
                 (goi.elem$all$cell.count * 100) /
                 Matrix::rowSums(x = goi.elem$all$cell.count)
               
               return(goi.elem)
               
             })
    
    return(goi)
    
  }


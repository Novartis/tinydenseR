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

#' Manually assign cell type labels to clusters
#'
#' This function allows manual annotation of clusters by mapping one or more
#' clusters to biologically meaningful cell type labels. This is particularly
#' useful in cytometry experiments, where there typically is a larger number of
#' cells and lineage marker-based hierarchical analysis is common.
#'
#' @details
#' Manual cell type annotation is useful when:
#' \itemize{
#'   \item Domain expertise is needed for fine-grained annotations
#'   \item Custom groupings are required for specific analyses
#'   \item Working with cytometry data where marker combinations define cell types
#' }
#' 
#' Cell type labels assigned here will be used in downstream analyses:
#' \itemize{
#'   \item \code{\link{get.map}} propagates labels to all cells via nearest landmarks
#'   \item \code{\link{get.stats}} enables cell-type-level differential abundance testing
#'   \item \code{\link{get.dea}} enables cell-type-specific differential expression analysis
#' }
#'
#' @param .lm.obj A list object initialized with \code{setup.lm.obj} and processed 
#'   with \code{get.landmarks} and \code{get.graph}. Must contain 
#'   \code{.lm.obj$graph$clustering$ids}.
#' @param .celltyping.map A named list mapping cell types to cluster IDs. 
#'   List names are cell type labels (e.g., "CD4.T.cells"), and each element 
#'   is a character vector of cluster IDs (e.g., c("cluster.01", "cluster.02")) 
#'   that belong to that cell type. Requirements:
#'   \itemize{
#'     \item Each cell type name must be unique
#'     \item Each cluster can only map to one cell type
#'     \item All clusters in \code{.lm.obj$graph$clustering$ids} must be mapped
#'   }
#' @examples
#' \dontrun{
#' # After clustering with get.graph()
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#' 
#' # Map clusters to cell types
#' celltype_map <- list(
#'   "CD4.T.cells" = c("cluster.01", "cluster.02"),
#'   "CD8.T.cells" = c("cluster.03"),
#'   "B.cells" = c("cluster.04", "cluster.05")
#' )
#' lm.cells <- celltyping(lm.cells, celltype_map)
#' }
#' @return The \code{.lm.obj} with the following updated fields:
#'   \describe{
#'     \item{\code{$graph$celltyping$ids}}{Factor vector of cell type labels 
#'       for each landmark}
#'     \item{\code{$graph$celltyping$mean.exprs}}{Matrix of mean expression 
#'       values per cell type (rows) by features (columns). For RNA data, 
#'       displays z-scored expression of top PC-loading genes. For cytometry 
#'       data, shows z-scored marker intensities.}
#'     \item{\code{$graph$celltyping$pheatmap}}{Heatmap object (from 
#'       \code{pheatmap} package) visualizing mean expression patterns across 
#'       cell types}
#'   }
#' @seealso 
#'   \code{\link{get.map}} for automatic cell typing with symphony reference objects,
#'   \code{\link{lm.cluster}} for the clustering that produces cluster IDs used here
#' @export
celltyping <-
  function(.lm.obj,
           .celltyping.map){
    
    cell.pop <- NULL
    
    if(names(x = .celltyping.map) |>
       is.null()){
      stop("Cell type names are missing. Each element in .celltyping.map must have a name.\n",
           "Example: list(\"CD4.T.cells\" = c(\"cluster.01\"), \"CD8.T.cells\" = c(\"cluster.02\"))")
    }
    
    if(names(x = .celltyping.map) |>
       duplicated() |>
       any()){
      stop(paste0("Duplicate cell type names detected: ",
                  paste(names(.celltyping.map)[duplicated(names(.celltyping.map))],
                        collapse = ", "),
                  "\nEach cell type must have a unique name."))
    }
    
    cls.in.map <-
      unlist(x = .celltyping.map,
             use.names = FALSE)
    
    if(anyDuplicated(x = cls.in.map)){
      stop(paste0("Cluster(s) mapped to multiple cell types: ",
                  paste(cls.in.map[duplicated(x = cls.in.map)],
                        collapse = ", "),
                  "\nEach cluster can only belong to one cell type."))
    }
    
    if(any(!(unique(x = .lm.obj$graph$clustering$ids) %in%
             cls.in.map))){
      stop(paste0("Every cluster must be mapped to a cell type. Unmapped clusters: ",
                  paste(unique(x = .lm.obj$graph$clustering$ids)[
                    !(unique(x = .lm.obj$graph$clustering$ids) %in%
                        cls.in.map)],
                    collapse = ", "),
                  "\nAdd these to .celltyping.map or merge them with existing clusters."))
      
    }
    
    if(!all(cls.in.map %in%
            unique(x = .lm.obj$graph$clustering$ids))){
      stop(paste0("Invalid cluster ID(s) in .celltyping.map: ",
                  paste(cls.in.map[!(cls.in.map %in%
                                       unique(x = .lm.obj$graph$clustering$ids))],
                        collapse = ", "),
                  "\nThese clusters do not exist in .lm.obj$graph$clustering$ids.",
                  "\nRun lm.cluster() to see available cluster IDs."))
    }
    
    .lm.obj$graph$celltyping <-
      vector(mode = "list")
    
    # Create cell type labels by:
    # 1. Inverting the map: cluster IDs -> cell type names
    # 2. Indexing by cluster IDs to get corresponding cell type for each landmark
    # 3. Converting to factor to maintain consistency with clustering$ids
    .lm.obj$graph$celltyping$ids <-
      .celltyping.map |>
      (\(x)
       stats::setNames(
         object = names(x = x) |>
           rep(times = lapply(X = x,
                              length) |>
                 unlist(use.names = FALSE)),
         nm = unlist(x = x,
                     use.names = FALSE))[
                       as.character(x = .lm.obj$graph$clustering$ids)]
      )() |>
      as.factor()
    
    # For RNA: select top PC-loading genes for heatmap visualization
    # Takes top 3 positive and top 3 negative loadings per PC to capture
    # genes that drive the major axes of variation
    if(.lm.obj$assay.type == "RNA") {
      
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
    
    return(.lm.obj)
    
  }

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

#' Map clusters to celltypes.
#'
#' This function maps clusters to cell types.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks and get.graph.
#' @param .celltyping.map A named list, where names are celltypes and each elemt is a character vector of cluster mapping to each respective celltype. A celltype can be comprised of more than one cluster, but each cluser can only be mapped to one celltype.
#' @returns The .lm.obj with the celltyping field populated with the celltype mapping, mean celltype expression values as well as the mean expression heatmap.
#' @export
celltyping <-
  function(.lm.obj,
           .celltyping.map){
    
    cell.pop <- NULL
    
    if(names(x = .celltyping.map) |>
       is.null()){
      stop("names of .celltyping.map cannot be NULL")
    }
    
    if(names(x = .celltyping.map) |>
       duplicated() |>
       any()){
      stop("names of .celltyping.map cannot be duplicated")
    }
    
    cls.in.map <-
      unlist(x = .celltyping.map,
             use.names = FALSE)
    
    if(anyDuplicated(x = cls.in.map)){
      stop(paste0("following cluster(s) mapping to more than one celltyping: ",
                  paste(cls.in.map[duplicated(x = cls.in.map)],
                        collapse = ", ")))
    }
    
    if(any(!(unique(x = .lm.obj$graph$clustering$ids) %in%
             cls.in.map))){
      stop(paste0("every cluster must be mapped to at least one celltyping: ",
                  paste(unique(x = .lm.obj$graph$clustering$ids)[
                    !(unique(x = .lm.obj$graph$clustering$ids) %in%
                        cls.in.map)],
                    collapse = ", ")))
      
    }
    
    if(!all(cls.in.map %in%
            unique(x = .lm.obj$graph$clustering$ids))){
      stop(paste0("at least one cluster in .celltyping.map is not present in .lm.obj$graph$clustering$ids: ",
                  paste(cls.in.map[!(cls.in.map %in%
                                       unique(x = .lm.obj$graph$clustering$ids))],
                        collapse = ", ")))
    }
    
    .lm.obj$graph$celltyping <-
      vector(mode = "list")
    
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

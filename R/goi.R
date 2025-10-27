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

#' Creates several summaries from a list of genes of interest.
#' 
#' This function creates several summaries typically used in scRNAseq data 
#' analysis from  a list of genes of interest (GOI).
#' 
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with `get.graph` and `get.map`.
#' @param .goi A character vector of genes of interest.
#' @param .id.idx An index of the cluster or celltype. If NULL, all landmarks are used.
#' @param .id The id of the cluster or celltype. If NULL, all landmarks are used.
#' @param .id.from The source of the id. Can be set to "clustering" or "celltyping". Defaults to NULL.
#' @param .verbose Logical indicating whether to print verbose output.
#' @returns The .lm.obj with the goi field populated with the updated cluster/celltype names.
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
      stop(paste0("Following genes could not be found: ",
                  paste(.goi[!(.goi %in% colnames(x = .lm.obj$raw.lm))],
                        collapse = ", ")))
    }
    
    if(.lm.obj$assay.type != "RNA"){
      stop("Currently, goi.summary only works with RNA.")
    }
    
    if(is.null(x = .lm.obj$map)){
      stop("First run tinydenseR `get.graph` and `get.map`.")
    }
    
    goi <-
      vector(mode = "list",
             length = 0)
    
    if(is.null(x = .id.idx)){
      
      if(!is.null(x = .id)){
        
        .id.from <-
          match.arg(arg = .id.from,
                    choices = c("clustering",
                                "celltyping"))
        
        if(!all(.id %in% unique(x = .lm.obj$graph[[.id.from]]$ids))){
          
          stop(paste0(paste0(.id[!(.id %in% unique(x = .lm.obj$graph[[.id.from]]$ids))],
                             collapse = ", "),
                      " not found in ",
                      .id.from))
          
        }
        
        .id.idx <-
          lapply(X = .lm.obj$map[[.id.from]]$ids,
                 FUN = function(smpl){
                   which(x = smpl %in% .id)
                 })
        
      } else {
        
        .id.idx <-
          lapply(X = .lm.obj$map[[.id.from]]$ids,
                 FUN = seq_along)
        
      }
    } else {
      
      .id.idx <-
        lapply(X = .lm.obj$map$nearest.landmarks,
               FUN = function(smpl.knn)
                 which(x = smpl.knn[,1] == .id.idx))
      
    }
    
    goi <-
      seq_along(along.with = .lm.obj$cells) |>
      stats::setNames(nm = names(x = .lm.obj$cells)) |>
      lapply(FUN = function(cells.elem){
        
        goi.exprs.mat <-
          readRDS(file = .lm.obj$cells[[cells.elem]]) |> 
          (\(x)
           x[.goi,.id.idx[[cells.elem]],drop = FALSE]
          )() |>
          Matrix::t()
        
        all.goi.detected <-
          as.matrix(x = goi.exprs.mat > 0)
        
        ids <-
          .lm.obj$map$clustering$ids[[cells.elem]][.id.idx[[cells.elem]]] |>
          rep(times = length(x = .goi)) |>
          matrix(byrow = FALSE,
                 nrow = length(x = .id.idx[[cells.elem]]),
                 ncol = length(x = .goi),
                 dimnames = list(names(x = .lm.obj$map$clustering$ids[[cells.elem]][.id.idx[[cells.elem]]]),
                                 .goi))
        
        ids[all.goi.detected] <-
          paste0("pos.",
                 ids[all.goi.detected])
        
        ids[!all.goi.detected] <-
          paste0("neg.",
                 ids[!all.goi.detected])
        
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
    
    goi <-
      lapply(X = goi, 
             FUN = function(goi.elem){
               
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
               
               goi.elem$clustering$cell.count[
                 is.na(x = goi.elem$clustering$cell.count)
               ] <- 0
               
               goi.elem$clustering$cell.perc <-
                 (goi.elem$clustering$cell.count * 100) /
                 Matrix::rowSums(x = goi.elem$clustering$cell.count)
               
               if(!is.null(x = .lm.obj$map$celltyping)){
                 
                 goi.elem$celltyping$cell.count <-
                   seq_along(along.with = goi.elem$celltyping$ids) |>
                   stats::setNames(nm = names(x = goi.elem$celltyping$ids)) |>
                   lapply(X = ,
                          FUN = function(smpl.idx){
                            
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


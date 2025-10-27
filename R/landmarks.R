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

#' Setup object for landmark identification and modeling
#'
#' This function initializes an object for landmark identification and modeling.
#'
#' @param .cells A named list containing the on-disk location of RDS files containing one expression matrix for each sample. List names correspond to sample names.
#' @param .meta A data frame containing the sample-wise metadata.
#' @param .markers A character vector specifying the markers to use for landmark identification. Default is NULL. Applies only if `assay.type` is `"cyto"` since all features are included for `"RNA"`.
#' @param .harmony.var A character vector specifying the metadata variables to use for Harmony integration. Default is NULL. If provided, these variables must be present in the metadata.
#' @param .assay.type A character string specifying the assay type. Default is `"cyto"`. Must be one of `"cyto"` or `"RNA"`.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @param .seed An integer specifying the .seed for reproducibility. Default is 123.
#' @param .prop.landmarks A numeric value between 0 and 1 specifying the proportion of cells to use as landmarks. Default is 0.1 (10% of total cells), with a maximum of 5000 landmarks.
#' @param .n.threads An integer specifying the number of threads to use. Default is the max allowed by the system.
#' @return A list object initialized with the input parameters.
#' @export
#'
setup.lm.obj <-
  function(
    .cells,
    .meta,
    .markers = NULL,
    .harmony.var = NULL,
    .assay.type = "cyto",
    .verbose = TRUE,
    .seed = 123,
    .prop.landmarks = 0.1,
    .n.threads = if(is.hpc()){
      max(RhpcBLASctl::blas_get_num_procs(),
          RhpcBLASctl::omp_get_num_procs(),
          RhpcBLASctl::omp_get_max_threads(),
          na.rm = TRUE)
    } else {
      parallel::detectCores(logical = TRUE)
    }){
    
    .cell <- NULL
    
    .assay.type <-
      match.arg(arg = .assay.type,
                choices = c("cyto",
                            "RNA"))
    
    if(!is.null(x = .harmony.var)){
      
      if(!inherits(x = .harmony.var,
                   what = "character")){
        stop(".harmony.var must be a character vector")
      }
      
      if(!all(.harmony.var %in% colnames(x = .meta))){
        stop(paste0("The following variables in .harmony.var were not found in the metadata: ",
                    paste0(.harmony.var[!(.harmony.var %in% colnames(x = .meta))],
                           collapse = ", ")))
      }
      
    }
    
    if(!is.null(x = .markers)){
      if(.assay.type == "RNA"){
        stop(".markers argument is only applicable for cytometry data analysis" )
      } else if(length(x = .markers) < 3){
        stop(".markers must contain at least 3 markers")
      }
    }
    
    if(!inherits(x = .meta,
                 what = "data.frame")){
      stop(".meta must be a data.frame object")
    }
    
    if(!all(rownames(x = .meta) ==
            names(x = .cells))){
      stop("names of .cells must be the same as rownames of .meta")
    }
    
    if(!inherits(x = .cells,
                 what = "list")){
      stop(".cells must be a list object")
    }
    
    if((.prop.landmarks < 0) | (.prop.landmarks > 1)){
      stop(".prop.landmarks must be a numeric value between 0 and 1")
    }
    
    .lm.obj.nm <-
      c("cells",
        "lm",
        "scaled.lm",
        "raw.lm",
        "metadata",
        "key",
        "pca",
        "graph",
        "map",
        "assay.type",
        "spec",
        "n.threads",
        "markers",
        "harmony.var",
        "interact.plot",
        "harmony.obj")
    
    .lm.obj <- vector(mode = "list",
                      length = length(x = .lm.obj.nm)) |>
      stats::setNames(nm = .lm.obj.nm)
    
    .lm.obj$n.threads <-
      .n.threads
    
    .lm.obj$assay.type <-
      .assay.type
    
    if(lapply(X = .cells,
              FUN = inherits,
              what = "character") |>
       unlist() |>
       all()){
      
      n.cells <-
        seq_along(along.with = .cells) |>
        stats::setNames(nm = names(x = .cells)) |>
        lapply(FUN = function(cell_elem){
          
          if(isTRUE(x = .verbose)){
            paste0("reading expression matrix: ",
                   names(x = .cells)[cell_elem],
                   "\n",
                   "progress: ",
                   round(x = cell_elem * 100 / length(x = .cells),
                         digits = 2),
                   "%") |>
              message()
          }
          
          .cell <- readRDS(file = .cells[[cell_elem]])
          
          if(!is.null(x = .markers)){
            if(!all(.markers %in% colnames(x = .cell))){
              stop(paste0("\n",
                          .markers,
                          " could not be found in the data. Arg .markers must be one or more of the following: ",
                          paste0(colnames(x = .cell),
                                 collapse = ", ")))
            }
          }
          
          if(rownames(x = .cell) |>
             is.null()){
            
            stop(paste0("Expression matrices in .cells must have rownames. Missing rownames in: ",
                        names(x = .cells)[cell_elem]))
            
          } else if(colnames(x = .cell) |>
                    is.null()){
            
            stop(paste0("Expression matrices in .cells must have colnames. Missing colnames in: ",
                        names(x = .cells)[cell_elem]))
            
          } else {
            
            if(.assay.type == "RNA"){
              return(ncol(x = .cell))
            } else {
              return(nrow(x = .cell))
            }
            
          }
          
          if(isTRUE(x = .verbose)){
            paste0("progress: ",
                   round(x = cell_elem * 100 / length(x = .cells),
                         digits = 2),
                   "%") |>
              message()
          }
          
          
        }) |>
        unlist()
      
      .lm.obj$cells <- .cells
      
    } else {
      
      stop(".cells has to be a named list of characters pointing to the location of the count matrices on disk.")
      
    }
    
    .lm.obj$spec$n.cells <-
      n.cells
    
    if((max(.lm.obj$spec$n.cells) / min(.lm.obj$spec$n.cells)) > 10){
      
      if(any(.lm.obj$spec$n.cells < 1000)){
        warning("A very large difference in the number of cells across samples was detected. In particular for cytometry experiments, it is advisable to remove samples with fewer than 1000 cells.")
      }
      
      warning(paste0("\nThe ratio of the largest to the smallest number of cells in the dataset is greater than 10.\n\nThis can lead to unreliable results.\n\nThe number of cells in the smallest sample is ",
                     min(.lm.obj$spec$n.cells),
                     ".\n\nMake sure to exclude low quality samples and proceed at your own discretion."))
      
    } 
    
    .lm.obj$spec$target.lm.n <-
      pmin(sum(.lm.obj$spec$n.cells) * .prop.landmarks,
           5e3)
    
    .lm.obj$spec$n.perSample <-
      pmin(ceiling(x = .lm.obj$spec$n.cells * .prop.landmarks),
           ceiling(x = .lm.obj$spec$target.lm.n / length(x = .lm.obj$cells)))
    
    .lm.obj$key <-
      seq_along(along.with = .lm.obj$cells) |>
      rep(times = .lm.obj$spec$n.perSample) |>
      (\(x)
       stats::setNames(object = x,
                       nm = names(x = .lm.obj$cells)[x])
      )()
    
    .lm.obj$metadata <-
      .meta
    
    .lm.obj$metadata$n.perSample <- 
      .lm.obj$spec$n.perSample
    
    .lm.obj$metadata$n.cells <-
      .lm.obj$spec$n.cells
    
    .lm.obj$metadata$log10.n.cells <-
      log10(x = .lm.obj$spec$n.cells)
    
    if(.assay.type == "cyto"){
      .lm.obj$markers <- 
        if(is.null(x = .markers)) readRDS(file = .lm.obj$cells[[1]]) |> colnames() else .markers
    }
    
    if(!is.null(x = .harmony.var)){
      .lm.obj$harmony.var <- 
        .harmony.var  
    }
    
    return(.lm.obj)
    
  }

#' Get landmarks
#'
#' This function identifies landmarks from the input expression matrices.
#'
#' @param .lm.obj An object initialized with setup.lm.obj.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @param .seed An integer specifying the .seed for reproducibility. Default is 123.
#' @param .nHVG An integer specifying the number of highly variable genes to use for landmark identification. Default is 5000
#' @param .nPC An integer specifying the number of principal components to use for landmark identification. Default is 30.
#' @param .exc.vdj.mito.ribo.genes.from.hvg A logical specifying whether to exclude VDJ, mitochondrial, and ribosomal genes from the highly variable genes. Default is TRUE.
#' @param .force.in A character vector specifying the genes to force into the landmark set. Default is NULL. Applies only if `assay.type` is `"RNA"` since all features are included for `"cyto"`.
#' @return The .lm.obj now containing the landmarks, a key linking each landmark to a row of .metadata and the PCA object.
#' @export
#'
get.landmarks <-
  function(.lm.obj,
           .verbose = TRUE,
           .seed = 123,
           .nHVG = 5000,
           .nPC = 30,
           .exc.vdj.mito.ribo.genes.from.hvg = TRUE,
           .force.in = NULL){
    
    set.seed(seed = .seed)
    
    if(names(x = .lm.obj) |>
       is.null()){
      stop("Input .lm.obj is not valid. Make sure to initialize it with setup.lm.obj")
    }
    
    if(!identical(
      x = names(x = .lm.obj), 
      y = c("cells",
            "lm",
            "scaled.lm",
            "raw.lm",
            "metadata",
            "key",
            "pca",
            "graph",
            "map",
            "assay.type",
            "spec",
            "n.threads",
            "markers",
            "harmony.var",
            "interact.plot",
            "harmony.obj"))){
      stop("Input .lm.obj is not valid. Make sure to initialize it with setup.lm.obj")
    }
    
    if(any(.nPC > min(.lm.obj$metadata$n.cells))){
      stop(paste0("Number of principal components (.nPC) must be less than or equal to the number of cells in the smallest sample:",
                  min(.lm.obj$metadata$n.cells)))
    }
    
    .lm.obj$raw.lm <-
      seq_along(along.with = .lm.obj$cells) |>
      stats::setNames(nm = names(x = .lm.obj$cells)) |>
      lapply(FUN = function(.cells.idx){
        
        mat <-
          readRDS(file = .lm.obj$cells[[.cells.idx]])
        
        if(.lm.obj$assay.type == "RNA"){
          
          if(!inherits(x = mat,
                       what = "dgCMatrix")){
            stop("RNA expression matrices must be a sparse matrix of class dgCMatrix")
          }
          
          if(isTRUE(x = .exc.vdj.mito.ribo.genes.from.hvg)) {
            
            vdj.mito.ribo <-
              grep(pattern = "^TR[AB][VDJ]\\d+|^IG[KHL][VDJ]\\d+|^RP|^MT-",
                   x = colnames(x = mat),
                   fixed = FALSE,
                   value = TRUE)
            
          } else {
            
            vdj.mito.ribo <-
              NULL
            
          }
          
          mat <-
            Matrix::t(x = mat) |>
            (\(x)
             x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))
            )()
          
          mat@x <-
            log2(x = mat@x + 1)
          
          HVG <-
            sparseMatrixStats::colVars(x = mat[,!colnames(x = mat) %in% vdj.mito.ribo]) |>
            order(decreasing = TRUE) |>
            utils::head(n = .nHVG)
          
          mat <-
            mat[,HVG]
          
        } else {
          
          mat <-
            mat[,.lm.obj$markers]
          
        }
        
        mat.mean <-
          Matrix::colMeans(x = mat)
        
        mat.sd <-
          sparseMatrixStats::colSds(x = mat)
        
        if(.lm.obj$assay.type == "RNA"){
          
          lev.score <-
            irlba::irlba(A = mat,
                         nv = .nPC,
                         center = mat.mean,
                         scale = mat.sd)$u |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        } else {
          
          mat <-
            ((Matrix::t(x = mat) - mat.mean) /
               mat.sd) |>
            Matrix::t()
          
          lev.score <-
            base::svd(x = mat,
                      nu = ncol(x = mat),
                      nv = ncol(x = mat))$u |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        }
        
        lm.sample <-
          sample(x = nrow(x = mat),
                 size = .lm.obj$spec$n.perSample[[.cells.idx]],
                 replace = FALSE,
                 prob = lev.score)
        
        landmarks <-
          if(.lm.obj$assay.type == "RNA"){
            
            readRDS(file = .lm.obj$cells[[.cells.idx]])[,lm.sample] |>
              Matrix::t() |>
              (\(x)
               `rownames<-`(x = x,
                            value = paste0(names(x = .lm.obj$cells)[.cells.idx],
                                           "_",
                                           rownames(x = x)))
              )()
            
          } else {
            readRDS(file = .lm.obj$cells[[.cells.idx]])[lm.sample,.lm.obj$markers]
          }
        
        if(isTRUE(x = .verbose)){
          paste0("progress: ",
                 round(x = .cells.idx * 100 / length(x = .lm.obj$cells),
                       digits = 2),
                 "%") |>
            message()
        }
        
        return(landmarks)
        
      }) |>
      do.call(what = rbind)
    
    if(.lm.obj$assay.type == "RNA"){
      
      .lm.obj$raw.lm <-
        methods::as(object = .lm.obj$raw.lm,
                   Class = "CsparseMatrix")
      
    }
    
    if(isTRUE(x = .verbose)){
      message("getting landmark pca")
    }
    
    if(.lm.obj$assay.type == "RNA"){
      
      if(isTRUE(x = .exc.vdj.mito.ribo.genes.from.hvg)) {
        
        vdj.mito.ribo <-
          grep(pattern = "^TR[AB][VDJ]\\d+|^IG[KHL][VDJ]\\d+|^RP|^MT-",
               x = colnames(x = .lm.obj$raw.lm),
               fixed = FALSE,
               value = TRUE,
               ignore.case = TRUE)
        
      } else {
        
        vdj.mito.ribo <-
          NULL
        
      }
      
      if(!is.null(x = .force.in)){
        if(inherits(x = .force.in,
                    what = "character")){
          if(any(!.force.in %in% colnames(x = .lm.obj$raw.lm))){
            
            warning(paste0(paste0(.force.in[!(.force.in %in% colnames(x = .lm.obj$raw.lm))],
                                  collapse = ", "),
                           "could not be found in the gene expression set"))
            
          }
          
          vdj.mito.ribo <-
            vdj.mito.ribo[!vdj.mito.ribo %in% .force.in]
          
          .force.in.idx <-
            which(x = colnames(x = .lm.obj$raw.lm)[!colnames(x = .lm.obj$raw.lm) %in% vdj.mito.ribo] %in% .force.in)
          
        } else {
          stop(".force.in must be either NULL or a character vector")
        }
      } else {
        .force.in.idx <-
          vector(mode = "integer",
                 length = 0)
      }
      
      .lm.obj$lm <-
        .lm.obj$raw.lm /
        (Matrix::rowSums(x = .lm.obj$raw.lm) /
           mean(x = Matrix::rowSums(x = .lm.obj$raw.lm)))
      
      .lm.obj$lm@x <-
        log2(x = .lm.obj$lm@x + 1)
      
      HVG <-
        sparseMatrixStats::colVars(
          x = .lm.obj$lm[,!colnames(x = .lm.obj$lm) %in% vdj.mito.ribo]) |>
        order(decreasing = TRUE) |>
        utils::head(n = .nHVG) |>
        (\(x)
         c(.force.in.idx,
           x) |>
           unique()
        )() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .lm.obj$lm)[!colnames(x = .lm.obj$lm) %in% vdj.mito.ribo][x])
        )()
      
    } else {
      
      .lm.obj$lm <-
        .lm.obj$raw.lm
      
      HVG <-
        ncol(x = .lm.obj$lm) |>
        seq_len() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .lm.obj$lm)[x])
        )()
      
    }
    
    .lm.obj$lm <-
      if(.lm.obj$assay.type == "RNA") {
        .lm.obj$lm[,colnames(x = .lm.obj$lm)[!colnames(x = .lm.obj$lm) %in% vdj.mito.ribo][HVG]]
      } else {
        .lm.obj$lm[,HVG]
      }
    
    .lm.obj$pca$center <-
      Matrix::colMeans(x = .lm.obj$lm)
    
    .lm.obj$pca$scale <-
      sparseMatrixStats::colSds(x = .lm.obj$lm)
    
    if(.lm.obj$assay.type == "RNA"){
      
      pca_res <-
        irlba::irlba(
          A = .lm.obj$lm,
          nv = .nPC,
          center = .lm.obj$pca$center,
          scale = .lm.obj$pca$scale)
      
    } else {
      
      pca_res <-
        Matrix::t(x = .lm.obj$lm) |>
        (\(x)
         (x - .lm.obj$pca$center) /
           .lm.obj$pca$scale
        )() |>
        Matrix::t() |>
        base::svd()
      
    }
    
    colnames(x = pca_res$u) <- colnames(x = pca_res$v) <- names(x = pca_res$d) <-
      paste0("PC",
             seq_along(along.with = pca_res$d))
    
    rownames(x = pca_res$u) <-
      rownames(x = .lm.obj$lm)
    
    rownames(x = pca_res$v) <-
      colnames(x = .lm.obj$lm)
    
    pca_res$embed <-
      pca_res$u %*%
      (base::diag(x = pca_res$d) |>
         (\(x)
          `dimnames<-`(x= x,
                       value = list(names(x = pca_res$d),
                                    names(x = pca_res$d)))
         )())
    
    .lm.obj$pca <-
      c(pca_res,
        .lm.obj$pca)
    
    if(!is.null(x = .lm.obj$harmony.var)){
      
      set.seed(seed = .seed)
      tmp.harmony <- 
        harmony::RunHarmony(data_mat = .lm.obj$pca$embed, 
                            meta_data = .lm.obj$metadata[.lm.obj$key,],
                            vars_use = .lm.obj$harmony.var, 
                            return_obj = TRUE,
                            # see https://github.com/immunogenomics/harmony/blob/b36bab002c1767af6e665c81f186b40a87870e64/R/ui.R#L191
                            nclust = min(round(x = nrow(x = .lm.obj$pca$embed) / 30), 100) |>
                              max(20),
                            verbose = .verbose)
      
      .lm.obj$harmony.obj <-
        symphony::buildReferenceFromHarmonyObj(
          harmony_obj = tmp.harmony,
          metadata = .lm.obj$metadata[.lm.obj$key,],
          vargenes_means_sds =
            dplyr::tibble(
              symbol = colnames(x = .lm.obj$lm),
              mean = .lm.obj$pca$center,
              stddev = .lm.obj$pca$scale
            ),
          pca_loadings = .lm.obj$pca$v,
          verbose = .verbose,
          do_umap = FALSE,
          save_uwot_path = NULL,
          #umap_min_dist = 0.1,
          seed = .seed)  
      
      if(isTRUE(x = .verbose)){
        message("Returning harmony corrected embedding.")
      }
      
      .lm.obj$pca$embed <-
        Matrix::t(x = .lm.obj$harmony.obj$Z_corr) |> 
        (\(x)
         `rownames<-`(x = x,
                      value = rownames(x = .lm.obj$lm))
         )()
      
    }
    
    # use learned dataset-wide PCA to re-assess sample-wise landmarks leverage scores
    .lm.obj$raw.lm <-
      seq_along(along.with = .lm.obj$cells) |>
      stats::setNames(nm = names(x = .lm.obj$cells)) |> 
      lapply(FUN = function(.cells.idx){
        
        set.seed(seed = .seed)
        
        mat <-
          readRDS(file = .lm.obj$cells[[.cells.idx]])
        
        if(.lm.obj$assay.type == "RNA"){
          
          mat <-
            Matrix::t(x = mat) |>
            (\(x)
             # size factor normalization, taking into consideration size factor of landmarks
             (x / (Matrix::rowSums(x = x) /
                     mean(x = Matrix::rowSums(x = x)))) * 
               (mean(x = Matrix::rowSums(x = x)) / mean(x = Matrix::rowSums(x = .lm.obj$raw.lm)))
            )()
          
          mat@x <-
            log2(x = mat@x + 1)
          
        }
        
        if(!is.null(x = .lm.obj$harmony.obj)){
          
          # Map query
          corr_embed <-
            symphony::mapQuery(exp_query = Matrix::t(x = mat[,colnames(x = .lm.obj$lm)]),
                               metadata_query = matrix(data = 1,
                                                       nrow = nrow(x = mat),
                                                       ncol = 2),
                               ref_obj = .lm.obj$harmony.obj,
                               vars = NULL,
                               verbose = .verbose,
                               do_normalize = FALSE,
                               do_umap = FALSE)$Z
          
          Y <-
            Matrix::t(x = mat[,colnames(x = .lm.obj$lm)]) |>
            (\(x)
             (x - .lm.obj$pca$center) /
               .lm.obj$pca$scale
            )() |>
            Matrix::t()
          
          Yt.X <-
            Matrix::crossprod(x = (corr_embed - Matrix::rowMeans(x = corr_embed)) |>
                                Matrix::t(),
                              y = Y)
          
          res <-
            base::svd(x = Yt.X, 
                      nu = nrow(x = corr_embed),
                      nv = nrow(x = corr_embed))
          
          dimnames(x = res$v) <-
            list(colnames(x = .lm.obj$lm),
                 rownames(x = corr_embed))
          
          #res$v <-
          #  (Matrix::t(x = res$v) * 
          #     (res$u[(abs(x = res$u) |> 
          #               matrixStats::rowRanks()) == ncol(x = res$u)] |>
          #        sign())) |>
          #  Matrix::t()
          
          lev.score <-
            ((Y %*% res$v) %*%
               diag(x = 1 / res$d)) |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        } else {
          
          lev.score <-
            (((((Matrix::t(x = mat[,colnames(x = .lm.obj$lm)]) - .lm.obj$pca$center) /
                  .lm.obj$pca$scale) |>
                 Matrix::t()) %*%
                .lm.obj$pca$v) %*%
               diag(x = 1 / .lm.obj$pca$d)) |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        }
        
        lm.sample <-
          sample(x = nrow(x = mat),
                 size = .lm.obj$spec$n.perSample[[.cells.idx]],
                 replace = FALSE,
                 prob = lev.score)
        
        landmarks <-
          if(.lm.obj$assay.type == "RNA"){
            
            readRDS(file = .lm.obj$cells[[.cells.idx]])[,lm.sample] |>
              Matrix::t() |>
              (\(x)
               `rownames<-`(x = x,
                            value = paste0(names(x = .lm.obj$cells)[.cells.idx],
                                           "_",
                                           rownames(x = x)))
              )()
            
          } else {
            readRDS(file = .lm.obj$cells[[.cells.idx]])[lm.sample,.lm.obj$markers]
          }
        
        if(isTRUE(x = .verbose)){
          paste0("progress: ",
                 round(x = .cells.idx * 100 / length(x = .lm.obj$cells),
                       digits = 2),
                 "%") |>
            message()
        }
        
        return(landmarks)
        
      }) |>
      do.call(what = rbind)
    
    if(isTRUE(x = .verbose)){
      message("getting landmark pca")
    }
    
    if(.lm.obj$assay.type == "RNA"){
      
      .lm.obj$lm <-
        .lm.obj$raw.lm /
        (Matrix::rowSums(x = .lm.obj$raw.lm) /
           mean(x = Matrix::rowSums(x = .lm.obj$raw.lm)))
      
      .lm.obj$lm@x <-
        log2(x = .lm.obj$lm@x + 1)
      
      HVG <-
        sparseMatrixStats::colVars(
          x = .lm.obj$lm[,!colnames(x = .lm.obj$lm) %in% vdj.mito.ribo]) |>
        order(decreasing = TRUE) |>
        utils::head(n = .nHVG) |>
        (\(x)
         c(.force.in.idx,
           x) |>
           unique()
        )() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .lm.obj$lm)[!colnames(x = .lm.obj$lm) %in% vdj.mito.ribo][x])
        )()
      
    } else {
      
      .lm.obj$lm <-
        .lm.obj$raw.lm
      
      HVG <-
        ncol(x = .lm.obj$lm) |>
        seq_len() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .lm.obj$lm)[x])
        )()
      
    }
    
    .lm.obj$lm <-
      if(.lm.obj$assay.type == "RNA") {
        .lm.obj$lm[,colnames(x = .lm.obj$lm)[!colnames(x = .lm.obj$lm) %in% vdj.mito.ribo][HVG]]
      } else {
        .lm.obj$lm[,HVG]
      }
    
    .lm.obj$pca$center <-
      Matrix::colMeans(x = .lm.obj$lm)
    
    .lm.obj$pca$scale <-
      sparseMatrixStats::colSds(x = .lm.obj$lm)
    
    if(.lm.obj$assay.type == "RNA"){
      
      pca_res <-
        irlba::irlba(
          A = .lm.obj$lm,
          nv = .nPC,
          center = .lm.obj$pca$center,
          scale = .lm.obj$pca$scale)
      
    } else {
      
      pca_res <-
        Matrix::t(x = .lm.obj$lm) |>
        (\(x)
         (x - .lm.obj$pca$center) /
           .lm.obj$pca$scale
        )() |>
        Matrix::t() |>
        base::svd()
      
    }
    
    .lm.obj$pca <-
      c(pca_res,
        .lm.obj$pca)
    
    dimnames(x = .lm.obj$pca$u) <-
      list(rownames(x = .lm.obj$lm),
           paste0("PC",
                  seq_along(along.with = .lm.obj$pca$d)))
    
    names(x = .lm.obj$pca$d) <-
      paste0("PC",
             seq_along(along.with = .lm.obj$pca$d))
    
    dimnames(x = .lm.obj$pca$v) <-
      list(colnames(x = .lm.obj$lm),
           paste0("PC",
                  seq_along(along.with = .lm.obj$pca$d)))
    
    .lm.obj$pca$embed <-
      .lm.obj$pca$u %*%
      (base::diag(x = .lm.obj$pca$d) |>
         (\(x)
          `dimnames<-`(x= x,
                       value = list(names(x = .lm.obj$pca$d),
                                    names(x = .lm.obj$pca$d)))
         )())
    
    .lm.obj$pca$sdev <-
      .lm.obj$pca$d /
      sqrt(x = max(1, nrow(x = .lm.obj$lm) - 1))
    
    .lm.obj$pca$rotation <-
      .lm.obj$pca$v
    
    .lm.obj$pca$v <- .lm.obj$pca$d <- .lm.obj$pca$u <-
      NULL
    
    .lm.obj$pca$HVG <-
      names(x = HVG)
    
    if(.lm.obj$assay.type == "RNA"){
      
      .lm.obj$scaled.lm <-
        (.lm.obj$pca$embed %*% Matrix::t(x = .lm.obj$pca$rotation)) |>
        Matrix::t() |>
        (\(x)
         (x - Matrix::rowMeans(x = x)) /
           matrixStats::rowSds(x = x)
        )() |>
        Matrix::t() |>
        (\(x)
         `dimnames<-`(x = x,
                      value = list(rownames(x = .lm.obj$lm),
                                   colnames(x = .lm.obj$lm)))
        )()
      
      .lm.obj$lm <-
        (.lm.obj$pca$embed %*% Matrix::t(x = .lm.obj$pca$rotation)) |>
        Matrix::t() |>
        (\(x)
         (x * .lm.obj$pca$scale) +
           .lm.obj$pca$center)() |>
        Matrix::t() |>
        (\(x)
         `dimnames<-`(x = x,
                      value = list(rownames(x = .lm.obj$lm),
                                   colnames(x = .lm.obj$lm)))
        )()
      
    } else {
      
      .lm.obj$scaled.lm <-
        Matrix::t(x = .lm.obj$lm) |>
        (\(x)
         (x - Matrix::rowMeans(x = x)) /
           matrixStats::rowSds(x = x)
        )() |>
        Matrix::t()
      
    }
    
    if(!is.null(x = .lm.obj$harmony.var)){
      
      warning("PCA embedding and rotations will be replaced by harmony-corrected embedding and roughly approximated rotations.")
      
      set.seed(seed = .seed)
      tmp.harmony <- 
        harmony::RunHarmony(data_mat = .lm.obj$pca$embed, 
                            meta_data = .lm.obj$metadata[.lm.obj$key,],
                            vars_use = .lm.obj$harmony.var, 
                            return_obj = TRUE,
                            # see https://github.com/immunogenomics/harmony/blob/b36bab002c1767af6e665c81f186b40a87870e64/R/ui.R#L191
                            nclust = min(round(x = nrow(x = .lm.obj$pca$embed) / 30), 100) |>
                              max(20),
                            verbose = .verbose)
      
      .lm.obj$harmony.obj <-
        symphony::buildReferenceFromHarmonyObj(
          harmony_obj = tmp.harmony,
          metadata = .lm.obj$metadata[.lm.obj$key,],
          vargenes_means_sds =
            dplyr::tibble(
              symbol = colnames(x = .lm.obj$lm),
              mean = .lm.obj$pca$center,
              stddev = .lm.obj$pca$scale
            ),
          pca_loadings = .lm.obj$pca$rotation,
          verbose = .verbose,
          do_umap = FALSE,
          save_uwot_path = NULL,
          #umap_min_dist = 0.1,
          seed = .seed)  
      
      if(isTRUE(x = .verbose)){
        message("Returning harmony corrected embedding.")
      }
      
      .lm.obj$pca$embed <-
        Matrix::t(x = .lm.obj$harmony.obj$Z_corr)  |> 
        (\(x)
         `rownames<-`(x = x,
                      value = rownames(x = .lm.obj$lm))
        )()
      
      if(isTRUE(x = .verbose)){
        message("Returning approximated rotations.")
      }
      
      # Build design matrix Z with intercept
      Z <-
        stats::model.matrix(object = ~ 1 + ., 
                            data = .lm.obj$metadata[.lm.obj$key, 
                                                    .lm.obj$harmony.var,
                                                    drop = FALSE])
      
      # Residualize each column in .lm.obj$lm (unfortunately best can do is linear)
      Y.resid <- 
        stats::lm.fit(x = Z, # fast base R linear model fit
                      y = as.matrix(x = .lm.obj$lm)) |>  
        stats::residuals()  # same shape as Y
      
      # Scale by stored stddev
      Y.resid <- 
        Matrix::t(x = Y.resid) |> 
        (\(x)
         x / .lm.obj$pca$scale
        )() |>
        Matrix::t()
      
      Yt.X <-
        Matrix::crossprod(x = (.lm.obj$harmony.obj$Z_corr - Matrix::rowMeans(x = .lm.obj$harmony.obj$Z_corr)) |>
                            Matrix::t(),
                          y = Y.resid)
      
      res <-
        base::svd(x = Yt.X, 
                  nu = nrow(x = .lm.obj$harmony.obj$Z_corr),
                  nv = nrow(x = .lm.obj$harmony.obj$Z_corr))
      
      dimnames(x = res$v) <-
        list(colnames(x = .lm.obj$lm),
             rownames(x = .lm.obj$harmony.obj$Z_corr))
      
      res$v <-
        (Matrix::t(x = res$v) * 
           (res$u[(abs(x = res$u) |> 
                     matrixStats::rowRanks(ties.method = "average")) ==
                    ncol(x = res$u)] |>
              sign())) |>
        Matrix::t()
      
      .lm.obj$pca$rotation <-
        res$v
      
    }
    
    return(.lm.obj)
    
  }

#' Check if running on HPC
is.hpc <- function() {
  # Common environment variables for HPC schedulers
  scheduler_vars <- c("SLURM_JOB_ID", "LSB_JOBID", "PBS_JOBID", "SGE_JOB_ID", "OAR_JOB_ID")
  
  # Check if any of the scheduler variables are present
  for (var in scheduler_vars) {
    if (Sys.getenv(var) != "") {
      return(TRUE)
    }
  }
  
  # If no scheduler variables are found, assume it's a local machine
  return(FALSE)
}

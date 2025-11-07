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

#' Get stats
#'
#' Fits a model of outcomes as a function of the sample wise sum of cell-landmark neighboring probabilities
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.map.
#' @param .design A matrix specifying the design matrix to be used.
#' @param .contrasts A matrix specifying the contrasts to be used.
#' @param .block A character vector of length one indicating which column of `.lm.obj$metadata` to be used as blocking factor.
#' @param .verbose A logical indicating whether to print progress messages.
#' @param .seed An integer to set the seed for reproducibility.
#' @return The .lm.obj now containing the stats from fitting a voomLmFit model.
#' @export
get.stats <-
  function(
    .lm.obj,
    .design,
    .contrasts = NULL,
    .block = NULL,
    .verbose = TRUE,
    .seed = 123){
    
    if(nrow(x = .design) != length(x = .lm.obj$cells)){
      stop("Number of rows in design matrix must be equal to the number of samples")
    }
    
    # check if no intercept is present but any of the covariates is a continuous variable
    no.intercept <-
      (!(0L %in% attr(.design,
                      which = "assign")))
    
    num.covariate <-
      !all(as.vector(x = .design) %in% c(0,1))
    
    if(no.intercept & num.covariate){
      warning("No intercept in the model but at least one of the covariates is a continuous variable. This model is not the same as the one with an intercept because it assumes that the continuous variable is centered at 0.")
    }
    
    #	Check design (source: https://rdrr.io/bioc/edgeR/src/R/glmfit.R)
    nlib <- length(x = .lm.obj$cells)
    .design <- as.matrix(x = .design)
    if(nrow(x = .design) != nlib) {
      stop("nrow(design) disagrees with length(x = .lm.obj$cells)")
    }
    ne <- limma::nonEstimable(x = .design)
    if(!is.null(x = ne)) {
      stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n",
                 paste(ne,
                       collapse = " ")))
    }
    
    if(is.null(x = .lm.obj$map)){
      stop("First run get.map")
    }
    
    # create list to hold results
    stats <-
      vector(mode = "list",
             length = 0)
    
    #stats$counts <-
    #  seq_along(along.with = .lm.obj$map$nearest.landmarks) |>
    #  stats::setNames(nm = names(x = .lm.obj$map$nearest.landmarks)) |>
    #  lapply(FUN = function(smpl.idx){
    #    
    #    res <-
    #      table(.lm.obj$map$nearest.landmarks[[smpl.idx]][,1]) |>
    #      (\(x)
    #       stats::setNames(object = as.vector(x = x),
    #                       nm = rownames(x = .lm.obj$lm)[names(x = x) |>
    #                                                       as.integer()])
    #      )() |>
    #      (\(x)
    #       x[match(x = rownames(x = .lm.obj$lm),
    #               table = names(x = x))]
    #      )()
    #    
    #    res[is.na(x = res)] <- 0
    #    
    #    names(x = res) <-
    #      rownames(x = .lm.obj$lm)
    #    
    #    return(res)
    #    
    #  }) |>
    #  dplyr::bind_cols() |>
    #  as.matrix() |>
    #  (\(x)
    #   `rownames<-`(x = x,
    #                value = rownames(x = .lm.obj$lm))
    #  )()
    
    stats$y <-
      log2(x = .lm.obj$map$fdens + 0.5)
    
    if(nrow(x = stats$y) !=
       nrow(x = .lm.obj$lm)){
      
      stats$y <-
        stats$y[match(x = rownames(x = .lm.obj$lm),
                      table = rownames(x = stats$y)),]
      
      stats$y[
        is.na(x = stats$y)
      ] <- min(stats$y,
               na.rm = TRUE)
      
      rownames(x = stats$y) <-
        rownames(x = .lm.obj$lm)
      
    }
    
    if(!(is.null(x = .block))){
      
      if(length(x = .block) != 1){
        stop("Block must be a vector of length 1")
      }
      if(!(.block %in% colnames(x = .lm.obj$metadata))){
        stop(paste0(.block,
                    " not found in metadata"))
      }
      
      if(isTRUE(x = .verbose)){
        message("\nestimating the intra-block correlation")
      }
      
      # https://support.bioconductor.org/p/125489/#125602
      # duplicateCorrelation is more general and is THE ONLY SOLUTION when
      # you want to compare across blocking levels, e.g., comparing diseased
      # and healthy donors when each donor also contributes before/after treatment samples.
      dupcor <- 
        limma::duplicateCorrelation(object = stats$y,
                                    design = .design,
                                    block = .lm.obj$metadata[[.block]])
    }
    
    if(isTRUE(x = .verbose)){
      message("\nfitting linear models")
    }
    
    stats$fit <-
      limma::lmFit(object = stats$y,
                   design = .design,
                   block = if(exists(x = "dupcor")) .lm.obj$metadata[[.block]] else  NULL,
                   correlation = if(exists(x = "dupcor")) dupcor$consensus else NULL)
    
    if(!is.null(x = .contrasts)){
      
      stats$fit <-
        limma::contrasts.fit(fit = stats$fit,
                             contrasts = .contrasts)
      
    }
    
    stats$fit <-
      limma::eBayes(fit = stats$fit,
                    robust = TRUE)
    
    if(isTRUE(x = .verbose)){
      message("\nretrieving density-weighted BH FDR")
    }
    
    stats$fit$density.weighted.bh.fdr <-
      apply(X = stats$fit$p.value,
            MARGIN = 2,
            FUN = function(pvalues){
              
              # source: https://github.com/MarioniLab/miloR/blob/38bbebec6e3dd2cac526a81afa56de1fad7bdd1d/R/graphSpatialFDR.R
              haspval <- !is.na(x = pvalues)
              
              if (!all(haspval)) {
                pvalues <- pvalues[haspval]
              }
              
              # use 1/sum of probabilities as the weighting for the weighted BH adjustment from Cydar
              w <- 1 / log10(x = Matrix::rowSums(x = .lm.obj$map$fdens) + 1)
              w[is.infinite(x = w)] <- 1
              
              # Computing a density-weighted q-value.
              o <- order(pvalues)
              pvalues <- pvalues[o]
              w <- w[o]
              
              adjp <- numeric(length = length(x = o))
              adjp[o] <- rev(x = cummin(x = rev(x = sum(w)*pvalues/cumsum(x = w))))
              adjp <- pmin(adjp, 1)
              
              if (!all(haspval)) {
                refp <- rep(x = NA_real_,
                            length(x = haspval))
                refp[haspval] <- adjp
                adjp <- refp
              }
              return(adjp)
              
            })
    
    if(isTRUE(x = .verbose)){
      message("\nretrieving PCA-weighted q-values")
    }
    
    if(nrow(x = stats$fit$p.value) < 100){
      warning("q-value estimation is not recommended for fewer than 100 tests. Using BH instead.")
    } else if(nrow(x = stats$fit$p.value) < 1000){
      warning("PCA-weighted q-value estimation is not recommended for fewer than 1000 tests. Using standard q-value estimation instead.")
    } else {
      
      # regress design matrix out of PCs to guarantee independence from covariates.
      
      # use limma so block structure can be residualized as needed
      
      # keep all PCs for now since it's cheap to residualize against them
      # and stabilizes the correlation estimates in dupCor and blocked LM
      # but keep only top k later to avoid overfitting pi0
      tX <-
        Matrix::t(x = .lm.obj$pca$embed)
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .lm.obj$metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        # https://support.bioconductor.org/p/125489/#125602
        # duplicateCorrelation is more general and is THE ONLY SOLUTION when
        # you want to compare across blocking levels, e.g., comparing diseased
        # and healthy donors when each donor also contributes before/after treatment samples.
        dupcor.q <- 
          limma::duplicateCorrelation(object = tX,
                                      design = stats$fit$design[.lm.obj$key,],
                                      block = .lm.obj$metadata[[.block]][.lm.obj$key])
      }
      
      fit.q <-
        limma::lmFit(object = tX,
                     design = stats$fit$design[.lm.obj$key,],
                     block = if(exists(x = "dupcor.q")) .lm.obj$metadata[[.block]][.lm.obj$key] else  NULL,
                     correlation = if(exists(x = "dupcor.q")) dupcor.q$consensus else NULL)
      
      # we only need residuals that remove all design effects, so
      # contrasts here are counter-productive bc they actually change
      # the coefficients in the limma object and make it impossible
      # to retrieve accurate fitted values for all covariates.
      #if(!is.null(x = .contrasts)){
      #  
      #  fit.q <-
      #    limma::contrasts.fit(fit = fit.q,
      #                         contrasts = .contrasts)
      #  
      #}
      
      X.resid <- 
        # observed - fitted
        Matrix::t(x = tX -
                    # fitted values:
                    (fit.q$coefficients %*% Matrix::t(x = fit.q$design))) |>
        (\(x)
         # keep only top k PCs that capture most variation
         # to avoid overfitting pi0 later
         x[,1:elbow.sec.deriv(x = .lm.obj$pca$sdev,
                              sort.order = "desc")$index,
           drop=FALSE]
        )()
      
    }
    
    stats$fit$pca.weighted.q <-
      apply(X = stats$fit$p.value,
            MARGIN = 2,
            FUN = function(pvalues){
              
              if(length(x = pvalues) < 100){
                
                res <- 
                  stats::p.adjust(p = pvalues,
                                  method = "BH")
                
                return(res)
                
              } else if(length(x = pvalues) < 1000){
                
                res <- 
                  tryCatch(expr = {
                    
                    suppressWarnings(expr = {
                      qvalue::qvalue(p = pvalues)$qvalues
                    })
                    
                  },
                  error = function(e){
                    
                    warning("q-value estimation failed. Using BH instead.")
                    
                    stats::p.adjust(p = pvalues,
                                    method = "BH")
                    
                  })
                return(res)
                
              } else {
                
                # at very large π₁ (most tests are truly non‑null), small n and
                # correlated tests, π₀ collapses to zero.
                # In such regimes, the covariate‑qvalue regression (and also
                # Storey’s global π₀) can underestimate π₀—sometimes to ~0, 
                # making q-values tiny for almost everything (anti‑conservative).
                # thus, guardrail by global pi0
                
                # pi0 estimates from the qvalue package is unfortunately unstable.
                # Here, we estimate pi0_hat for several λ values and take the median
                # as a conservative summary. This mimics the spirit of the smoother
                # but avoids the fragility of smooth.spline() in extreme cases
                # (like when all most p-values are very small).
                lam <- 
                  seq(from = 0.10, 
                      to = 0.80, 
                      by = 0.05)
                
                # Ensure p-values are finite and within [0,1]
                p.clean <-
                  pvalues[is.finite(x = pvalues)]
                p.clean <- 
                  pmax(p.clean, 0) |>
                  pmin(1)
                
                pi0.grid <- 
                  (outer(X = p.clean, 
                         Y = lam, 
                         FUN = `>`) |>
                     colMeans(na.rm = TRUE)) /
                  (1 - lam)
                
                # Replace non-finite with NA, then take a robust summary
                pi0.grid[!is.finite(x = pi0.grid)] <-
                  NA_real_
                pi0.global <- 
                  stats::median(x = pi0.grid, 
                                na.rm = TRUE)
                # clamp to [0,1]
                pi0.global <-
                  min(pi0.global, 1) |>
                  max(0)  
                
                if(pi0.global < 0.6) {
                  
                  warning("Global pi0 estimate is less than 0.6. Using fixed pi0 floor of 0.6 for PCA-weighted q-value estimation to avoid anti-conservative estimates.")
                  # compute with fixed pi0 floor
                  pi0.obj <-
                    suppressWarnings(expr = {
                      swfdr::lm_pi0(p = pvalues,
                                    X = X.resid)
                    })
                  
                  pi0.obj$pi0 <- 
                    pmin(pi0.obj$pi0, 1) |>
                    pmax(0.6)
                  
                  res <-
                    suppressWarnings(expr = {
                      swfdr::lm_qvalue(p = pvalues, 
                                       X = X.resid, 
                                       pi0 = pi0.obj)$qvalues
                    })
                  
                } else {
                  res <-
                    suppressWarnings(expr = {
                      swfdr::lm_qvalue(p = pvalues,
                                       X = X.resid)$qvalues  
                    })
                  
                }
                
                return(res)
                
              }
              
            }) |>
      (\(x)
       `dimnames<-`(x = x,
                    value = dimnames(x = stats$fit$p.value))
      )()
    
    if(isTRUE(x = .verbose)){
      message("\ngetting traditional stats")
    }
    
    stats$trad <-
      vector(mode = "list",
             length = 2) |>
      stats::setNames(nm = c("clustering",
                             "celltyping"))
    
    if(!(is.null(x = .block))){
      
      if(length(x = .block) != 1){
        stop("Block must be a vector of length 1")
      }
      if(!(.block %in% colnames(x = .lm.obj$metadata))){
        stop(paste0(.block,
                    " not found in metadata"))
      }
      
      if(isTRUE(x = .verbose)){
        message("\nestimating the intra-block correlation for stats from clustering")
      }
      
      cl.dupcor <- 
        log2(x = .lm.obj$map$clustering$cell.perc + 0.5) |>
        Matrix::t() |>
        limma::duplicateCorrelation(design = .design,
                                    block = .lm.obj$metadata[[.block]])
    }
    
    if(isTRUE(x = .verbose)){
      message("\nfitting linear models")
    }
    
    stats$trad$clustering$fit <-
      limma::lmFit(object = log2(x = .lm.obj$map$clustering$cell.perc + 0.5) |>
                     Matrix::t(),
                   design = .design,
                   block = if(exists(x = "cl.dupcor")) .lm.obj$metadata[[.block]] else NULL,
                   correlation = if(exists(x = "cl.dupcor")) cl.dupcor$consensus else NULL)
    
    if(!is.null(x = .contrasts)){
      
      stats$trad$clustering$fit <-
        limma::contrasts.fit(fit = stats$trad$clustering$fit,
                             contrasts = .contrasts)
      
    }
    
    stats$trad$clustering$fit <-
      limma::eBayes(fit = stats$trad$clustering$fit)
    
    stats$trad$clustering$fit$adj.p <-
      apply(X = stats$trad$clustering$fit$p.value,
            MARGIN = 2,
            FUN = stats::p.adjust,
            method = "fdr")
    
    if(!is.null(x = .lm.obj$graph$celltyping)){
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .lm.obj$metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        if(isTRUE(x = .verbose)){
          message("\nestimating the intra-block correlation for stats from celltyping")
        }
        
        ct.dupcor <- 
          log2(x = .lm.obj$map$celltyping$cell.perc + 0.5) |>
          Matrix::t() |>
          limma::duplicateCorrelation(design = .design,
                                      block = .lm.obj$metadata[[.block]])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      stats$trad$celltyping$fit <-
        limma::lmFit(object = log2(x = .lm.obj$map$celltyping$cell.perc + 0.5) |>
                       Matrix::t(),
                     design = .design,
                     block = if(exists(x = "ct.dupcor")) .lm.obj$metadata[[.block]] else NULL,
                     correlation = if(exists(x = "ct.dupcor")) ct.dupcor$consensus else NULL)
      
      if(!is.null(x = .contrasts)){
        
        stats$trad$celltyping$fit <-
          limma::contrasts.fit(fit = stats$trad$celltyping$fit,
                               contrasts = .contrasts)
        
      }
      
      stats$trad$celltyping$fit <-
        limma::eBayes(fit = stats$trad$celltyping$fit)
      
      stats$trad$celltyping$fit$adj.p <-
        apply(X = stats$trad$celltyping$fit$p.value,
              MARGIN = 2,
              FUN = stats::p.adjust,
              method = "fdr")
      
    }
    
    return(stats)
    
  }

#' Differential Expression Analysis
#'
#' This function performs differential expression analysis for each feature
#' using .coef as response. If provided, only cells within .id are used.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with
#' get.landmarks, get.graph, get.map and get.stats.
#' @param .design A matrix specifying the design matrix to be used.
#' @param .contrasts A matrix specifying the contrasts to be used.
#' @param .block A character vector of length one indicating which column of `.lm.obj$metadata` to be used as blocking factor.
#' @param .geneset.ls A list of character vectors specifying the gene sets to be used in differential gene set score analysis (via pseudobulk GSVA). Default is NULL (no gene set scoring is run).
#' @param .id.idx A list of integer vectors specifying the indices of the landmarks to be used. Default is NULL. If provided, cells in each sample that have the selected landmark as first neighbor are included.
#' @param .id A character vector specifying the ids to be used. Default is NULL.
#' @param .id.from A character specifying the source of the ids. Default is NULL.
#' Options are "clustering" and "celltyping".
#' @param .verbose A logical indicating whether to print progress messages. Default is TRUE.
#' @return The .lm.obj now containing the stats from fitting a linear model for each marker in each cluster.
#' @export
#'
get.dea <-
  function(
    .lm.obj,
    .design,
    .contrasts = NULL,
    .block = NULL,
    .geneset.ls = NULL,
    .id.idx = NULL,
    .id = NULL,
    .id.from = NULL,
    .verbose = TRUE
  ) {
    
    if(!is.null(x = .geneset.ls)){
      if(.lm.obj$assay.type != "RNA"){
        stop(".geneset.ls is only supported for RNA assay type")
      } else if(!is.list(x = .geneset.ls)){
        stop(".geneset.ls must be a list of character vectors")
      }
    } 
    
    if(nrow(x = .design) != length(x = .lm.obj$cells)){
      stop("Number of rows in design matrix must be equal to the number of samples")
    }
    
    # check if no intercept is present but any of the covariates is a continuous variable
    no.intercept <-
      !(any(Matrix::colSums(x = .design == 1) == nrow(x = .design)))
    
    num.covariate <-
      !all(as.vector(x = .design) %in% c(0,1))
    
    if(no.intercept & num.covariate){
      warning("No intercept in the model but at least one of the covariates is a continuous variable. This model is not the same as the one with an intercept because it assumes that the continuous variable is centered at 0.")
    }
    
    #	Check design (source: https://rdrr.io/bioc/edgeR/src/R/glmfit.R)
    nlib <- length(x = .lm.obj$cells)
    .design <- as.matrix(x = .design)
    if(nrow(x = .design) != nlib) {
      stop("nrow(design) disagrees with length(x = .lm.obj$cells)")
    }
    ne <- limma::nonEstimable(x = .design)
    if(!is.null(x = ne)) {
      stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n",
                 paste(ne,
                       collapse = " ")))
    }
    
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
          lapply(X = .lm.obj$map$nearest.landmarks,
                 FUN = nrow) |>
          lapply(FUN = seq_len)
        
      }
    } else {
      
      if(!all(.id.idx %in% (nrow(x = .lm.obj$lm) |> seq_len()))) {
        stop(paste0(".id.idx must be a list of integer vectors from  1 to ",
                    nrow(x = .lm.obj$lm)))
      }
      
      .id.idx <-
        lapply(X = .lm.obj$map$nearest.landmarks,
               FUN = function(smpl.knn)
                 which(x = smpl.knn[,1] %in% .id.idx))
      
    }
    
    # number of cells in each pseudobulk
    n.pseudo <-
      lengths(x = .id.idx)
    
    # remove outlier samples with too few cells
    smpl.outlier <-
      (n.pseudo / mean(x = n.pseudo)) < 0.1
    
    # get pseudobulk
    if(.lm.obj$assay.type == "RNA"){
      
      if(any(smpl.outlier)){
        warning(paste0("The following samples were removed since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
                       paste(names(x = smpl.outlier)[smpl.outlier],
                             collapse = "\n")))
      }
      
      .id.idx <-
        .id.idx[!smpl.outlier]
      
      counts <-
        stats::setNames(object = names(x = .lm.obj$cells)[!smpl.outlier],
                        nm = names(x = .lm.obj$cells)[!smpl.outlier]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          res <-
            Matrix::rowSums(x = exprs.mat[,.id.idx[[smpl]],drop = FALSE])
          
          return(res)
          
        }) |>
        do.call(what = cbind)
      
      dge <-
        edgeR::DGEList(counts = counts)
      
      tmp.design <-
        .design[!smpl.outlier,] |>
        (\(x)
         x[,Matrix::colSums(x = x == 0) != nrow(x = x),drop = FALSE]
        )()
      
      keep <-
        edgeR::filterByExpr(y = dge,
                            design = tmp.design)
      
      dge <-
        dge[keep,,keep.lib.sizes=FALSE]
      
      dge <-
        edgeR::calcNormFactors(object = dge,
                               method = "TMM")
      
      v <- limma::voom(counts = dge,
                       design = tmp.design,
                       plot = FALSE)
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .lm.obj$metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        if(isTRUE(x = .verbose)){
          message("\nestimating the intra-block correlation")
        }
        
        # https://support.bioconductor.org/p/125489/#125602
        # duplicateCorrelation is more general and is THE ONLY SOLUTION when
        # you want to compare across blocking levels, e.g., comparing diseased
        # and healthy donors when each donor also contributes before/after treatment samples.
        dupcor <- 
          limma::duplicateCorrelation(object = v,
                                      design = tmp.design,
                                      block = .lm.obj$metadata[[.block]][!smpl.outlier])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      .de <-
        limma::lmFit(object = v,
                     design = tmp.design,
                     block = if(exists(x = "dupcor")) .lm.obj$metadata[[.block]][!smpl.outlier] else  NULL,
                     correlation = if(exists(x = "dupcor")) dupcor$consensus else NULL)
      
      if(!is.null(x = .contrasts)){
        
        tmp.contrasts <-
          .contrasts[colnames(x = tmp.design),,drop = FALSE] |>
          (\(x)
           x[,Matrix::colSums(x = x == 0) != nrow(x = x),drop = FALSE]
          )()
        
        .de <-
          limma::contrasts.fit(fit = .de,
                               contrasts = tmp.contrasts)
        
      }
      
      .de <-
        limma::eBayes(fit = .de,
                      robust = TRUE)
      
      .de$adj.p <-
        apply(X = .de$p.value,
              MARGIN = 2,
              FUN = stats::p.adjust,
              method = "fdr")
      
      if(!is.null(x = .geneset.ls)){
        
        gsva.es <- 
          GSVA::gsvaParam(
            exprData = v$E,
            geneSets = .geneset.ls,
            kcdf = "Gaussian",
            maxDiff = TRUE) |>
          GSVA::gsva()
        
        if(!(is.null(x = .block))){
          
          if(isTRUE(x = .verbose)){
            message("\nestimating the intra-block correlation")
          }
          
          # https://support.bioconductor.org/p/125489/#125602
          # duplicateCorrelation is more general and is THE ONLY SOLUTION when
          # you want to compare across blocking levels, e.g., comparing diseased
          # and healthy donors when each donor also contributes before/after treatment samples.
          gsva.dupcor <- 
            limma::duplicateCorrelation(object = gsva.es,
                                        design = tmp.design,
                                        block = .lm.obj$metadata[[.block]][!smpl.outlier])
        }
        
        gsva.fit <- 
          limma::lmFit(object = gsva.es,
                       design = tmp.design,
                       block = if(exists(x = "gsva.dupcor")) .lm.obj$metadata[[.block]][!smpl.outlier] else  NULL,
                       correlation = if(exists(x = "gsva.dupcor")) gsva.dupcor$consensus else NULL)
        
        if(!is.null(x = .contrasts)){
          
          gsva.fit <-
            limma::contrasts.fit(fit = gsva.fit,
                                 contrasts = tmp.contrasts)
          
        }
        
        gsva.fit <- 
          limma::eBayes(fit = gsva.fit,
                        robust = TRUE)
        
        gsva.fit$adj.p <-
          apply(X = gsva.fit$p.value,
                MARGIN = 2,
                FUN = stats::p.adjust,
                method = "fdr")
        
        .de$geneset <-
          gsva.fit
        
        .de$geneset$E <-
          gsva.es
        
      }
      
    } else {
      
      counts <-
        stats::setNames(object = names(x = .lm.obj$cells),#[!smpl.outlier],
                        nm = names(x = .lm.obj$cells)) |>#[!smpl.outlier]) |>
        lapply(FUN = function(smpl){
          
          if(isTRUE(x = .verbose)){
            which(x = names(x = .lm.obj$cells) == smpl) |> 
              (\(x)
               paste0("progress: ",
                      round(x = x * 100 / length(x = .lm.obj$cells),
                            digits = 2),
                      "%")
              )()|> 
              message()
          }
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[,colnames(x = .lm.obj$lm)]
          
          res <-
            matrixStats::colMedians(x = exprs.mat[.id.idx[[smpl]],,drop = FALSE])
          
          return(res)
          
        }) |>
        do.call(what = rbind)
      
      tmp.design <-
        .design |>#[!smpl.outlier,] |>
        (\(x)
         x[,Matrix::colSums(x = x == 0) != nrow(x = x),drop = FALSE]
        )()
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .lm.obj$metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        if(isTRUE(x = .verbose)){
          message("\nestimating the intra-block correlation")
        }
        
        # https://support.bioconductor.org/p/125489/#125602
        # duplicateCorrelation is more general and is THE ONLY SOLUTION when
        # you want to compare across blocking levels, e.g., comparing diseased
        # and healthy donors when each donor also contributes before/after treatment samples.
        dupcor <- 
          limma::duplicateCorrelation(object = Matrix::t(x = counts),
                                      design = tmp.design,
                                      block = .lm.obj$metadata[[.block]])#[!smpl.outlier])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      .de <-
        limma::lmFit(object = Matrix::t(x = counts),
                     design = tmp.design,
                     block = if(exists(x = "dupcor")) .lm.obj$metadata[[.block]] else  NULL,#[!smpl.outlier] else  NULL,
                     correlation = if(exists(x = "dupcor")) dupcor$consensus else NULL)
      
      if(!is.null(x = .contrasts)){
        
        tmp.contrasts <-
          .contrasts[colnames(x = tmp.design),,drop = FALSE] |>
          (\(x)
           x[,Matrix::colSums(x = x == 0) != nrow(x = x),drop = FALSE]
          )()
        
        .de <-
          limma::contrasts.fit(fit = .de,
                               contrasts = tmp.contrasts)
        
      }
      
      .de <-
        limma::eBayes(fit = .de,
                      robust = TRUE)
      
      .de$adj.p <-
        apply(X = .de$p.value,
              MARGIN = 2,
              FUN = stats::p.adjust,
              method = "fdr")
      
    }
    
    .de$smpl.outlier <-
      smpl.outlier
    
    .de$id.idx <-
      .id.idx
    
    .de$n.pseudo <-
      n.pseudo
    
    return(.de)
    
  }

#' Differential Expression Analysis Comparing Cell Subsets
#'
#' This function performs differential expression analysis for each feature
#' comparing cell subsets.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with
#' get.landmarks, get.graph, get.map and get.stats.
#' @param .geneset.ls A list of character vectors specifying the gene sets to be used in differential gene set score analysis (via pseudobulk GSVA). Default is NULL (no gene set scoring is run).
#' @param .id1 A character vector specifying the ids to be used.
#' @param .id2 A character vector specifying the ids to be used. Default is `..all.other.subsets..`.
#' @param .id.from A character specifying the source of the ids. Options are "clustering" and "celltyping".
#' @param .id1.idx A list of integer vectors specifying the indices of the landmarks to be used for .id1. Default is NULL.
#' @param .id2.idx A list of integer vectors specifying the indices of the landmarks to be used for .id2. Default is NULL.
#' @return The .lm.obj now containing the stats from fitting a linear model for each marker in each cluster.
#' @export
#'
get.marker <-
  function(
    .lm.obj,
    .geneset.ls = NULL,
    .id1.idx = NULL,
    .id2.idx = NULL,
    .id1 = NULL,
    .id2 = "..all.other.subsets..",
    .id.from = "clustering"
  ) {
    
    if(!is.null(x = .geneset.ls)){
      if(.lm.obj$assay.type != "RNA"){
        stop(".geneset.ls is only supported for RNA assay type")
      } else if(!is.list(x = .geneset.ls)){
        stop(".geneset.ls must be a list of character vectors")
      }
    } 
    
    .id.from <-
      match.arg(arg = .id.from,
                choices = c("clustering",
                            "celltyping"))
    
    if(is.null(x = .id1) &&
       is.null(x = .id1.idx)){
      stop("Please provide .id1 or .id1.idx")
    }
    
    if(is.null(x = .id2) &&
       is.null(x = .id2.idx)){
      stop("Please provide .id2 or .id2.idx")
    }
    
    if(is.null(x = .id1.idx)){
      if(!all(.id1 %in% unique(x = .lm.obj$graph[[.id.from]]$ids))){
        
        stop(paste0(paste0(.id1[!(.id1 %in% unique(x = .lm.obj$graph[[.id.from]]$ids))],
                           collapse = ", "),
                    " not found in ",
                    .id.from))
        
      }
    }
    
    if(is.null(x = .id2.idx)){
      if(!all(.id2 %in% c("..all.other.subsets..",
                          unique(x = .lm.obj$graph[[.id.from]]$ids) |>
                          as.character()))){
        
        stop(paste0(paste0(.id2[!(.id2 %in% c("..all.other.subsets..",
                                              unique(x = .lm.obj$graph[[.id.from]]$ids) |>
                                                as.character()))],
                           collapse = ", "),
                    " not found in ",
                    .id.from))
        
      }
    }
    
    if(is.null(x = .id1.idx)){
      
      .id1.idx <-
        lapply(X = .lm.obj$map[[.id.from]]$ids,
               FUN = function(smpl){
                 which(x = smpl %in% .id1)
               })
      
    } else {
      
      if(!all(.id1.idx %in% (nrow(x = .lm.obj$lm) |> seq_len()))) {
        stop(paste0(".id1.idx must be a list of integer vectors from  1 to ",
                    nrow(x = .lm.obj$lm)))
      }
      
      .id1.idx <-
        lapply(X = .lm.obj$map$nearest.landmarks,
               FUN = function(smpl.knn)
                 which(x = smpl.knn[,1] %in% .id1.idx))
      
    }
    
    if("..all.other.subsets.." %in% .id2){
      
      message("using `..all.other.subsets..` for .id2")
      
      .id2 <- "..all.other.subsets.."
      
      .id2.idx <-
        names(x = .lm.obj$map$nearest.landmarks) |>
        stats::setNames(nm = names(x = .lm.obj$map$nearest.landmarks)) |>
        lapply(FUN = function(smpl.nm){
          (nrow(x = .lm.obj$map$nearest.landmarks[[smpl.nm]]) |>
             seq_len())[-.id1.idx[[smpl.nm]]]
        })
      
    } else {
      
      if(is.null(x = .id2.idx)){
        
        .id2.idx <-
          lapply(X = .lm.obj$map[[.id.from]]$ids,
                 FUN = function(smpl){
                   which(x = smpl %in% .id2)
                 })
        
      } else {
        
        if(!all(.id2.idx %in% (nrow(x = .lm.obj$lm) |> seq_len()))) {
          stop(paste0(".id2.idx must be a list of integer vectors from  1 to ",
                      nrow(x = .lm.obj$lm)))
        }
        
        .id2.idx <-
          lapply(X = .lm.obj$map$nearest.landmarks,
                 FUN = function(smpl.knn)
                   which(x = smpl.knn[,1] %in% .id2.idx))
        
      }
      
    }
    
    # number of cells in each .id1 and .id2 pseudobulks
    n.pseudo1 <-
      lengths(x = .id1.idx)
    
    n.pseudo2 <-
      lengths(x = .id2.idx)
    
    # warn if outlier samples with too few cells are detected
    smpl.outlier.1 <-
      (n.pseudo1 /
         mean(x = c(n.pseudo1,
                    n.pseudo2))) < 0.1
    
    smpl.outlier.2 <-
      (n.pseudo2 /
         mean(x = c(n.pseudo1,
                    n.pseudo2))) < 0.1
    
    # get pseudobulk
    if(.lm.obj$assay.type == "RNA"){
      
      if(any(smpl.outlier.1)){
        warning(paste0("The following samples were removed from .id1 since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
                       paste(names(x = smpl.outlier.1)[smpl.outlier.1],
                             collapse = "\n")))
        .id1.idx <-
          .id1.idx[!smpl.outlier.1]
        
      }
      
      if(any(smpl.outlier.2)){
        warning(paste0("The following samples were removed from .id2 since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
                       paste(names(x = smpl.outlier.2)[smpl.outlier.2],
                             collapse = "\n")))
        .id2.idx <-
          .id2.idx[!smpl.outlier.2]
        
      }
      
      counts1 <-
        names(x = .id1.idx) |>
        stats::setNames(nm = names(x = .id1.idx)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          res <-
            tryCatch(
              expr = {
                
                Matrix::rowSums(x = exprs.mat[,.id1.idx[[smpl]],drop = FALSE])
                
              },
              error = function(e) {
                
                stats::setNames(object = rep(x = 0,
                                             times = nrow(x = exprs.mat)),
                                nm = rownames(x = exprs.mat))
                
              })
          
          return(res)
          
        }) |>
        do.call(what = cbind) |>
        (\(x)
         `colnames<-`(x = x,
                      value = paste0("counts1.",
                                     colnames(x = x)))
        )()
      
      counts2 <-
        names(x = .id2.idx) |>
        stats::setNames(nm = names(x = .id2.idx)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          res <-
            tryCatch(
              expr = {
                
                Matrix::rowSums(x = exprs.mat[,.id2.idx[[smpl]],drop = FALSE])
                
              },
              error = function(e) {
                
                stats::setNames(object = rep(x = 0,
                                             times = nrow(x = exprs.mat)),
                                nm = rownames(x = exprs.mat))
                
              })
          
          return(res)
          
        }) |>
        do.call(what = cbind) |>
        (\(x)
         `colnames<-`(x = x,
                      value = paste0("counts2.",
                                     colnames(x = x)))
        )()
      
      dge <-
        cbind(counts1,
              counts2) |>
        edgeR::DGEList()
      
      .ids <-
        colnames(x = dge) |>
        strsplit(split = ".",
                 fixed = TRUE) |>
        lapply(FUN = "[",
               i = 1) |>
        unlist() |>
        gsub(pattern = "counts",
             replacement = ".id",
             fixed = TRUE) |>
        factor(levels = c(".id2",
                          ".id1"))
      
      .pairs <-
        colnames(x = dge) |>
        gsub(pattern = "^counts1.|^counts2.",
             replacement = "",
             fixed = FALSE)
      
      tmp.design <-
        stats::model.matrix(object = ~ .ids + .pairs) |> 
        (\(x)
         `colnames<-`(x = x,
                      value = colnames(x = x) |>
                        gsub(pattern = "^\\.ids|^\\.pairs",
                             replacement = "",
                             fixed = FALSE))
        )() |> 
        (\(x)
         # this is to ensure only true pairs are kept
         x[,Matrix::colSums(x = x) > 1]
        )()
      
      keep <-
        edgeR::filterByExpr(y = dge,
                            design = tmp.design)
      
      dge <-
        dge[keep,,keep.lib.sizes=FALSE]
      
      dge <-
        edgeR::calcNormFactors(object = dge,
                               method = "TMM")
      
      v <-
        limma::voom(counts = dge,
                    design = tmp.design,
                    plot = FALSE)
      
      .de <-
        limma::lmFit(object = v,
                     design = tmp.design)
      
      .de <-
        limma::eBayes(fit = .de,
                      robust = TRUE)
      
      .de$adj.p <-
        apply(X = .de$p.value,
              MARGIN = 2,
              FUN = stats::p.adjust,
              method = "fdr")
      
      if(!is.null(x = .geneset.ls)){
        
        gsva.es <- 
          GSVA::gsvaParam(
            exprData = v$E,
            geneSets = .geneset.ls,
            kcdf = "Gaussian",
            maxDiff = TRUE) |>
          GSVA::gsva()
        
        gsva.fit <- 
          limma::lmFit(object = gsva.es,
                       design = tmp.design)
        
        gsva.fit <- 
          limma::eBayes(fit = gsva.fit,
                        robust = TRUE)
        
        gsva.fit$adj.p <-
          apply(X = gsva.fit$p.value,
                MARGIN = 2,
                FUN = stats::p.adjust,
                method = "fdr")
        
        .de$geneset <-
          gsva.fit
        
        .de$geneset$E <-
          gsva.es
        
      }
      
    } else {
      
      counts1 <-
        names(x = .lm.obj$cells) |>
        stats::setNames(nm = names(x = .lm.obj$cells)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[,colnames(x = .lm.obj$lm)]
          
          res <-
            tryCatch(
              expr = {
                
                matrixStats::colMedians(x = exprs.mat[.id1.idx[[smpl]],,drop = FALSE])
                
              },
              error = function(e) {
                
                stats::setNames(object = rep(x = 0,
                                             times = ncol(x = exprs.mat)),
                                nm = colnames(x = exprs.mat))
                
              })
          
          return(res)
          
        }) |>
        do.call(what = rbind)  |>
        (\(x)
         `rownames<-`(x = x,
                      value = paste0("counts1.",
                                     rownames(x = x)))
        )()
      
      counts2 <-
        names(x = .lm.obj$cells) |>
        stats::setNames(nm = names(x = .lm.obj$cells)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[,colnames(x = .lm.obj$lm)]
          
          res <-
            tryCatch(
              expr = {
                
                matrixStats::colMedians(x = exprs.mat[.id2.idx[[smpl]],,drop = FALSE])
                
              },
              error = function(e) {
                
                stats::setNames(object = rep(x = 0,
                                             times = ncol(x = exprs.mat)),
                                nm = colnames(x = exprs.mat))
                
              })
          
          return(res)
          
        }) |>
        do.call(what = rbind)  |>
        (\(x)
         `rownames<-`(x = x,
                      value = paste0("counts2.",
                                     rownames(x = x)))
        )()
      
      counts <-
        rbind(counts1,
              counts2) |>
        Matrix::t()
      
      .ids <-
        colnames(x = counts) |>
        strsplit(split = ".",
                 fixed = TRUE) |>
        lapply(FUN = "[",
               i = 1) |>
        unlist() |>
        gsub(pattern = "counts",
             replacement = ".id",
             fixed = TRUE) |>
        factor(levels = c(".id2",
                          ".id1"))
      
      .pairs <-
        colnames(x = counts) |>
        gsub(pattern = "^counts1.|^counts2.",
             replacement = "",
             fixed = FALSE)
      
      tmp.design <-
        stats::model.matrix(object = ~ .ids + .pairs) |> 
        (\(x)
         `colnames<-`(x = x,
                      value = colnames(x = x) |>
                        gsub(pattern = "^\\.ids|^\\.pairs",
                             replacement = "",
                             fixed = FALSE))
        )() |> 
        (\(x)
         # this is to ensure only true pairs are kept
         x[,Matrix::colSums(x = x) > 1]
        )()
      
      .de <-
        limma::lmFit(object = counts,
                     design = tmp.design)
      
      .de <-
        limma::eBayes(fit = .de,
                      robust = TRUE)
      
      .de$adj.p <-
        apply(X = .de$p.value,
              MARGIN = 2,
              FUN = stats::p.adjust,
              method = "fdr")
      
    }
    
    .de$smpl.outlier.1 <-
      smpl.outlier.1
    
    .de$smpl.outlier.2 <-
      smpl.outlier.2
    
    .de$id1.idx <-
      .id1.idx
    
    .de$id2.idx <-
      .id2.idx
    
    .de$n.pseudo1 <-
      n.pseudo1
    
    .de$n.pseudo2 <-
      n.pseudo2
    
    return(.de)
    
  }

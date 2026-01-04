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

#' Differential Abundance Testing
#'
#' Performs landmark-based differential abundance (DA) testing using limma on log-transformed fuzzy 
#' densities. Tests which landmarks (cell populations) change in abundance between conditions. More 
#' sensitive than traditional cluster-level testing because landmarks can capture within-cluster 
#' heterogeneity. Uses PCA-weighted q-values that leverage the correlation structure among 
#' landmarks to improve statistical power.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .design Design matrix specifying the experimental design. Rows correspond to samples 
#'   (matching \code{.lm.obj$cells}), columns to coefficients. Create with \code{model.matrix()}.
#' @param .contrasts Optional contrast matrix for specific comparisons. Each column defines one 
#'   contrast. Create with \code{limma::makeContrasts()}. If NULL, tests all coefficients in 
#'   \code{.design}.
#' @param .block Optional character: column name in \code{.lm.obj$metadata} for blocking factor 
#'   (e.g., "Donor", "Batch"). Accounts for within-block correlation using \code{duplicateCorrelation}.
#' @param .verbose Logical: print progress messages? Default TRUE.
#' @param .seed Integer: random seed for reproducibility in blocking correlation estimation. Default 123.
#'   
#' @return List containing differential abundance results:
#'   \describe{
#'     \item{\code{y}}{log2(fuzzy density + 0.5) matrix: landmarks x samples}
#'     \item{\code{fit}}{limma MArrayLM object after eBayes with moderated statistics, containing:}
#'     \item{\code{fit$coefficients}}{Log fold changes (landmarks x coefficients), nested in fit}
#'     \item{\code{fit$p.value}}{P-values (landmarks x coefficients), nested in fit}
#'     \item{\code{fit$density.weighted.bh.fdr}}{Density-weighted BH FDR adjustment (landmarks x coefficients)}
#'     \item{\code{fit$pca.weighted.q}}{PCA-weighted q-values (landmarks x coefficients)}
#'     \item{\code{trad}}{List with traditional cluster-level DA results}
#'     \item{\code{trad$clustering$fit}}{limma fit for cluster percentages with \code{$adj.p} slot}
#'     \item{\code{trad$celltyping$fit}}{limma fit for celltype percentages with \code{$adj.p} slot (if celltyping exists)}
#'   }
#'   
#' @details
#' The landmark-based DA testing workflow:
#' \enumerate{
#'   \item Computes log2(fuzzy density + 0.5) for each landmark across samples
#'   \item Fits linear model: \code{lmFit(y ~ design)}
#'   \item If blocking specified: estimates within-block correlation via \code{duplicateCorrelation}
#'   \item Applies contrasts if provided
#'   \item Performs empirical Bayes moderation: \code{eBayes()}
#'   \item Computes density-weighted FDR and PCA-weighted q-values
#'   \item Performs traditional cluster/celltype-level DA for comparison
#' }
#' 
#' \strong{Advantages over cluster-level testing:}
#' \itemize{
#'   \item Detects subtle shifts within clusters
#'   \item No arbitrary clustering thresholds
#'   \item Continuous representation of cell state space
#'   \item Powered by limma's variance shrinkage
#' }
#' 
#' \strong{Design matrix considerations:}
#' \itemize{
#'   \item Include intercept for standard comparisons
#'   \item No-intercept models with continuous covariates assume zero-centering (warning issued)
#'   \item Must be full rank (checked automatically)
#' }
#' 
#' \strong{PCA-weighted q-value procedure:}
#' 
#' Standard FDR methods assume independent tests, but landmarks are correlated in 
#' the cell state space. To account for this:
#' 
#' \enumerate{
#'   \item \strong{Residualize PCs}: Regress out design matrix from landmark PCA embeddings to 
#'     obtain covariates independent of experimental factors
#'   \item \strong{Estimate pi0}: Use \code{swfdr::lm_pi0} to estimate the proportion of true 
#'     nulls conditional on PC coordinates. Landmarks in similar cell states have correlated 
#'     p-values; PCA captures this spatial structure
#'   \item \strong{Conservative safeguard}: If global pi0 < 0.6, apply pi0 floor of 0.6 to prevent 
#'     anti-conservative estimates in high-signal regimes
#'   \item \strong{Compute q-values}: Use \code{swfdr::lm_qvalue} to calculate spatially-weighted 
#'     FDR that properly accounts for landmark correlation structure
#' }
#' 
#' This approach is more powerful than standard BH adjustment when landmarks exhibit correlated 
#' expression, while maintaining proper FDR control. Falls back to standard q-value 
#' or BH for small test sets (<1000 landmarks).
#' 
#' @seealso \code{\link{get.map}} (required predecessor), \code{\link{plotBeeswarm}} for 
#'   visualization, \code{\link{plotTradStats}} for comparison with cluster-level tests
#'   
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#' 
#' # Simple two-group comparison
#' design <- model.matrix(~ Condition, data = .meta)
#' stats <- get.stats(lm.cells, .design = design)
#' 
#' # Visualize results
#' plotBeeswarm(lm.cells, stats, .coefs = "ConditionB")
#' 
#' # Complex design with contrasts
#' design <- model.matrix(~ 0 + Group, data = .meta)
#' contrasts <- limma::makeContrasts(
#'   TvsC = GroupTreatment - GroupControl,
#'   T1vsT2 = GroupTreatment1 - GroupTreatment2,
#'   levels = design
#' )
#' stats <- get.stats(lm.cells, .design = design, .contrasts = contrasts)
#' 
#' # With blocking for paired samples
#' design <- model.matrix(~ Timepoint, data = .meta)
#' stats <- get.stats(lm.cells, .design = design, .block = "Subject")
#' }
#' 
#' @export
#'
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
              
            }) |>
      (\(x)
       `dimnames<-`(x = x,
                    value = dimnames(x = stats$fit$p.value))
      )()
    
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
              
              if((is.na(x = pvalues) |> mean()) == 1){
                return(rep(x = NA_real_,
                           length.out = length(x = pvalues)))
              }
              
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
                
                # at very large pi1 (most tests are truly non-null), small n and
                # correlated tests, pi0 collapses to zero.
                # In such regimes, the covariate-qvalue regression (and also
                # Storey's global pi0) can underestimate pi0--sometimes to ~0, 
                # making q-values tiny for almost everything (anti-conservative).
                # thus, guardrail by global pi0
                
                # pi0 estimates from the qvalue package is unfortunately unstable.
                # Here, we estimate pi0_hat for several lambda values and take the median
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
#' Performs pseudobulk differential expression (DE) analysis for genes/markers. Aggregates cells 
#' into pseudobulk samples, then uses limma-voom (RNA) or limma (cytometry) to test for DE. Can 
#' restrict analysis to specific clusters/cell types, and optionally perform gene set enrichment 
#' via GSVA.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .design Design matrix specifying experimental design. Rows = samples, columns = coefficients.
#' @param .contrasts Optional contrast matrix for specific comparisons. Create with 
#'   \code{limma::makeContrasts()}. If NULL, tests all \code{.design} coefficients.
#' @param .block Optional character: column name in \code{.lm.obj$metadata} for blocking factor 
#'   (e.g., "Donor"). Accounts for within-block correlation.
#' @param .geneset.ls Optional named list of character vectors defining gene sets for GSVA enrichment 
#'   analysis. Only for RNA data. Example: \code{list("Tcell" = c("CD3D", "CD3E"), "Bcell" = c("CD19", "MS4A1"))}.
#' @param .id.idx Optional list of integer vectors specifying landmark indices to include. If provided, 
#'   only cells with these landmarks as nearest neighbors are aggregated. Useful for custom cell selection.
#' @param .id Optional character vector of cluster/celltype IDs to restrict analysis to. Uses 
#'   \code{.id.from} to determine source.
#' @param .id.from Character: "clustering" or "celltyping". Source of IDs when \code{.id} is specified.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return List containing DE analysis results:
#'   \describe{
#'     \item{\code{fit}}{limma MArrayLM fit object with moderated statistics}
#'     \item{\code{coefficients}}{Log fold change matrix (features x coefficients)}
#'     \item{\code{p.value}}{P-values (features x coefficients)}
#'     \item{\code{adj.p}}{FDR-adjusted p-values (features x coefficients)}
#'     \item{\code{smpl.outlier}}{Logical vector indicating which samples were excluded as outliers}
#'     \item{\code{id.idx}}{List of integer vectors showing which cells were included per sample}
#'     \item{\code{n.pseudo}}{Integer vector of pseudobulk cell counts per sample}
#'     \item{\code{geneset}}{(RNA only, if \code{.geneset.ls} provided) GSVA results with \code{$fit} 
#'       (limma fit), \code{$E} (enrichment scores), \code{$adj.p} (adjusted p-values)}
#'   }
#'   
#' @details
#' The pseudobulk DE workflow:
#' \enumerate{
#'   \item \strong{Cell selection}: If \code{.id} or \code{.id.idx} specified, filter to those cells. 
#'     Otherwise use all cells.
#'   \item \strong{Pseudobulk aggregation}: 
#'     \itemize{
#'       \item RNA: Sum raw counts across cells within each sample
#'       \item Cytometry: Compute column medians of marker expression across cells within each sample
#'     }
#'   \item \strong{Outlier removal}: (RNA only) Exclude samples with <10\% of average pseudobulk cell count
#'   \item \strong{Normalization}:
#'     \itemize{
#'       \item RNA: TMM normalization (\code{edgeR::calcNormFactors}) + voom transformation for variance modeling
#'       \item Cytometry: Data used as-is (assumes pre-transformed input)
#'     }
#'   \item \strong{Linear modeling}: \code{limma::lmFit} with optional blocking for paired/repeated samples
#'   \item \strong{Empirical Bayes}: Variance shrinkage via \code{limma::eBayes(robust = TRUE)}
#'   \item \strong{GSVA} (optional, RNA only): Gene set variation analysis on voom-normalized expression
#' }
#' 
#' \strong{When to use pseudobulk DE:}
#' \itemize{
#'   \item Identifying marker genes for cell types/states
#'   \item Testing treatment effects on gene expression
#'   \item Comparing expression between conditions within specific populations
#' }
#' 
#' \strong{Pseudobulk vs single-cell DE:}
#' Pseudobulk aggregation:
#' \itemize{
#'   \item Respects biological replication (samples as units)
#'   \item Appropriate statistical framework (no pseudoreplication)
#'   \item Better power for moderate effects
#'   \item Recommended for between-sample comparisons
#' }
#' 
#' @seealso \code{\link{get.map}} (required), \code{\link{plotDEA}} for heatmap visualization,
#'   \code{\link{get.marker}} for marker gene identification
#'   
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#' 
#' # DE analysis on all cells
#' design <- model.matrix(~ Condition, data = .meta)
#' dea <- get.dea(lm.cells, .design = design)
#' plotDEA(lm.cells, dea, .coefs = "ConditionB")
#' 
#' # DE within specific cell type
#' dea.tcell <- get.dea(lm.cells, .design = design,
#'                      .id = c("1", "2", "3"),
#'                      .id.from = "clustering")
#' 
#' # With gene set enrichment
#' hallmark.sets <- list(
#'   "INTERFERON_RESPONSE" = c("ISG15", "ISG20", "IFIT1"),
#'   "INFLAMMATORY_RESPONSE" = c("IL1B", "TNF", "CXCL8")
#' )
#' dea.gsva <- get.dea(lm.cells, .design = design,
#'                     .geneset.ls = hallmark.sets)
#' }
#' 
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
        
        .lm.idx <-
          which(x = .lm.obj$graph[[.id.from]]$ids %in% .id)

      } else {
        
        .lm.idx <-
          nrow(x = .lm.obj$lm) |>
          seq_len()
        
      }
    } else {
      
      if(!all(.id.idx %in% (nrow(x = .lm.obj$lm) |> seq_len()))) {
        stop(paste0(".id.idx must be a list of integer vectors from  1 to ",
                    nrow(x = .lm.obj$lm)))
      }
      
      .lm.idx <-
        .id.idx 

    }
    
    # number of cells in each pseudobulk
    n.pseudo <-
      lapply(X = .lm.obj$map$nearest.landmarks,
             FUN = function(smpl){
               sum(smpl[,1,drop = TRUE] %in% .lm.idx)
             }) |>
      unlist()
    
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
      
      counts <-
        stats::setNames(object = names(x = .lm.obj$cells)[!smpl.outlier],
                        nm = names(x = .lm.obj$cells)[!smpl.outlier]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          wcl <- 
            .lm.obj$map$connect.prob[[smpl]][,.lm.idx,drop = FALSE]
          
          # get weighted sum
          wsum <-
            exprs.mat %*% wcl
          
          # get weighted mean using sum of weights and scale by number of cells
          res <-
            (Matrix::rowSums(x = wsum) / sum(wcl)) * (sum(Matrix::rowSums(x = wcl) > 0))

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
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[,colnames(x = .lm.obj$lm)]
          
          wcl <- 
            .lm.obj$map$connect.prob[[smpl]][,.lm.idx,drop=FALSE]
          
          # get weighted sum
          wsum <-
            (Matrix::t(x = exprs.mat) %*% wcl)
          
          # get weighted mean using sum of weights
          res <-
            Matrix::rowSums(x = wsum) / sum(wcl)
          
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

#' Marker Gene/Protein Identification via Pairwise Comparison
#'
#' Identifies marker genes/proteins that distinguish one cell subset from another (or all others) by 
#' performing pseudobulk differential expression between groups. Unlike \code{get.dea} which tests 
#' experimental conditions, this compares cell populations to find defining features. Useful for 
#' characterizing clusters and validating cell type annotations.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .geneset.ls Optional named list of character vectors defining gene sets for GSVA enrichment. 
#'   Only for RNA data.
#' @param .id1.idx Optional list of integer vectors specifying landmark indices for group 1. If 
#'   provided, only cells with these landmarks as nearest neighbors are included.
#' @param .id2.idx Optional list of integer vectors specifying landmark indices for group 2. If 
#'   provided, only cells with these landmarks as nearest neighbors are included.
#' @param .id1 Character vector of cluster/celltype IDs for group 1 (test group). Required if 
#'   \code{.id1.idx} not provided.
#' @param .id2 Character vector of cluster/celltype IDs for group 2 (reference group). Default 
#'   \code{"..all.other.landmarks.."} compares group 1 to all other cells. Can specify specific IDs 
#'   for pairwise comparisons.
#' @param .id.from Character: "clustering" (default) or "celltyping". Source of IDs in \code{.id1} 
#'   and \code{.id2}.
#'   
#' @return limma MArrayLM fit object with moderated statistics and additional slots:
#'   \describe{
#'     \item{\code{coefficients}}{Log fold changes (features x coefficients), nested in fit}
#'     \item{\code{p.value}}{P-values (features x coefficients), nested in fit}
#'     \item{\code{adj.p}}{FDR-adjusted p-values (features x coefficients)}
#'     \item{\code{smpl.outlier.1}}{Logical vector: samples excluded from group 1}
#'     \item{\code{smpl.outlier.2}}{Logical vector: samples excluded from group 2}
#'     \item{\code{id1.idx}}{List of integer vectors: cells included per sample for group 1}
#'     \item{\code{id2.idx}}{List of integer vectors: cells included per sample for group 2}
#'     \item{\code{n.pseudo1}}{Integer vector: pseudobulk cell counts per sample for group 1}
#'     \item{\code{n.pseudo2}}{Integer vector: pseudobulk cell counts per sample for group 2}
#'     \item{\code{geneset}}{(RNA only, if \code{.geneset.ls} provided) GSVA fit object with \code{$adj.p} and \code{$E} (enrichment scores)}
#'   }
#'   
#' @details
#' The marker identification workflow:
#' \enumerate{
#'   \item \strong{Cell selection}: Extract cells belonging to group 1 (\code{.id1}) and group 2 
#'     (\code{.id2})
#'   \item \strong{Pseudobulk aggregation}: 
#'     \itemize{
#'       \item RNA: Sum raw counts across cells within each group per sample
#'       \item Cytometry: Compute column medians of marker expression across cells within each group per sample
#'     }
#'   \item \strong{Outlier removal}: (RNA only) Exclude samples with <10\% of average pseudobulk cell count from either group
#'   \item \strong{Paired comparison}: Fit model \code{~ .ids + .pairs} where \code{.ids} is the group indicator (.id1 vs .id2) and \code{.pairs} controls for sample of origin (blocking factor)
#'   \item \strong{Statistical testing}: limma-voom with TMM normalization (RNA) or limma (cytometry)
#'   \item \strong{FDR correction}: Adjust p-values across all features
#' }
#' 
#' \strong{Common use cases:}
#' \itemize{
#'   \item \strong{One-vs-all}: Default \code{.id2 = "..all.other.landmarks.."} finds markers that 
#'     uniquely identify a cluster
#'   \item \strong{Pairwise}: Specify both \code{.id1} and \code{.id2} to compare two specific 
#'     populations (e.g., CD4 vs CD8 T cells)
#'   \item \strong{Combined groups}: Provide multiple IDs in \code{.id1} to find markers for a 
#'     meta-population
#' }
#' 
#' \strong{Interpreting results:}
#' \itemize{
#'   \item Positive logFC: Higher expression in group 1 (potential markers)
#'   \item Negative logFC: Higher expression in group 2 (anti-markers)
#'   \item Look for both high fold change AND statistical significance
#' }
#' 
#' \strong{Difference from get.dea:}
#' \itemize{
#'   \item \code{get.dea}: Tests experimental conditions (treatment vs control) within populations
#'   \item \code{get.marker}: Tests cell populations (cluster A vs B) within the same experiment
#' }
#' 
#' @seealso \code{\link{get.map}} (required), \code{\link{get.dea}} for condition-based DE,
#'   \code{\link{plotHeatmap}} for visualizing expression patterns
#'   
#' @examples
#' \dontrun{
#' # After mapping and clustering
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#' 
#' # Find markers for cluster 1 vs all others
#' markers.cl1 <- get.marker(lm.cells,
#'                           .id1 = "cluster.1",
#'                           .id.from = "clustering")
#' 
#' # Top upregulated markers (accessing .id1 coefficient)
#' top.markers <- rownames(markers.cl1$coefficients)[
#'   markers.cl1$adj.p[, ".id1"] < 0.05 & 
#'   markers.cl1$coefficients[, ".id1"] > 1
#' ]
#' 
#' # Pairwise comparison: CD4 T cells vs CD8 T cells
#' markers.cd4vcd8 <- get.marker(lm.cells,
#'                               .id1 = c("cluster.1", "cluster.3"),  # CD4 T clusters
#'                               .id2 = c("cluster.2", "cluster.4"),  # CD8 T clusters
#'                               .id.from = "clustering")
#' 
#' # With gene set enrichment
#' hallmark.sets <- list(
#'   "GLYCOLYSIS" = c("HK2", "PFKP", "LDHA"),
#'   "OXIDATIVE_PHOS" = c("ATP5A", "COX5A", "NDUFA4")
#' )
#' markers.gsva <- get.marker(lm.cells,
#'                            .id1 = "cluster.1",
#'                            .geneset.ls = hallmark.sets)
#' }
#' 
#' @export
#'
get.marker <-
  function(
    .lm.obj,
    .geneset.ls = NULL,
    .id1.idx = NULL,
    .id2.idx = NULL,
    .id1 = NULL,
    .id2 = "..all.other.landmarks..",
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
      if(!all(.id2 %in% c("..all.other.landmarks..",
                          unique(x = .lm.obj$graph[[.id.from]]$ids) |>
                          as.character()))){
        
        stop(paste0(paste0(.id2[!(.id2 %in% c("..all.other.landmarks..",
                                              unique(x = .lm.obj$graph[[.id.from]]$ids) |>
                                                as.character()))],
                           collapse = ", "),
                    " not found in ",
                    .id.from))
        
      }
    }
    
    if(is.null(x = .id1.idx)){
      
      .lm1.idx <-
        which(x = .lm.obj$graph[[.id.from]]$ids %in% .id1)
      
    } else {
      
      if(!all(.id1.idx %in% (nrow(x = .lm.obj$lm) |> seq_len()))) {
        stop(paste0(".id1.idx must be a list of integer vectors from  1 to ",
                    nrow(x = .lm.obj$lm)))
      }
      
      .lm1.idx <-
        .id1.idx
      
    }
    
    if("..all.other.landmarks.." %in% .id2){
      
      message("using `..all.other.landmarks..` for .id2")
      
      .id2 <- "..all.other.landmarks.."
      
      .lm2.idx <-
        which(x = !(.lm.obj$graph[[.id.from]]$ids %in% .id1))
      
    } else {
      
      if(is.null(x = .id2.idx)){
        
        .lm2.idx <-
          which(x = .lm.obj$graph[[.id.from]]$ids %in% .id2)
        
      } else {
        
        if(!all(.id2.idx %in% (nrow(x = .lm.obj$lm) |> seq_len()))) {
          stop(paste0(".id2.idx must be a list of integer vectors from  1 to ",
                      nrow(x = .lm.obj$lm)))
        }
        
        .lm2.idx <-
          .id2.idx
        
      }
      
    }
    
    # number of cells in each .id1 and .id2 pseudobulks
    n.pseudo1 <-
      lapply(X = .lm.obj$map$nearest.landmarks,
             FUN = function(smpl){
               sum(smpl[,1,drop = TRUE] %in% .lm1.idx)
             }) |>
      unlist()
    
    n.pseudo2 <-
      lapply(X = .lm.obj$map$nearest.landmarks,
             FUN = function(smpl){
               sum(smpl[,1,drop = TRUE] %in% .lm2.idx)
             }) |>
      unlist()
    
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
      }
      
      if(any(smpl.outlier.2)){
        warning(paste0("The following samples were removed from .id2 since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
                       paste(names(x = smpl.outlier.2)[smpl.outlier.2],
                             collapse = "\n")))
      }
      
      counts1 <-
        stats::setNames(object = names(x = .lm.obj$cells)[!smpl.outlier.1],
                        nm = names(x = .lm.obj$cells)[!smpl.outlier.1]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          res <-
            tryCatch(
              expr = {
                
                wcl <- 
                  .lm.obj$map$connect.prob[[smpl]][,.lm1.idx,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  exprs.mat %*% wcl
                
                # get weighted mean using sum of weights and scale by number of cells
                tmp.res <-
                  (Matrix::rowSums(x = wsum) / sum(wcl)) * (sum(Matrix::rowSums(x = wcl) > 0))
                
                return(tmp.res)
                
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
        stats::setNames(object = names(x = .lm.obj$cells)[!smpl.outlier.2],
                        nm = names(x = .lm.obj$cells)[!smpl.outlier.2]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          res <-
            tryCatch(
              expr = {
                
                wcl <- 
                  .lm.obj$map$connect.prob[[smpl]][,.lm2.idx,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  exprs.mat %*% wcl
                
                # get weighted mean using sum of weights and scale by number of cells
                res <-
                  (Matrix::rowSums(x = wsum) / sum(wcl)) * (sum(Matrix::rowSums(x = wcl) > 0))
                
                #Matrix::rowSums(x = exprs.mat[,.id1.idx[[smpl]],drop = FALSE])
                
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
        stats::setNames(object = names(x = .lm.obj$cells),
                        nm = names(x = .lm.obj$cells)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[,colnames(x = .lm.obj$lm)]
          
          res <-
            tryCatch(
              expr = {
                
                wcl <- 
                  .lm.obj$map$connect.prob[[smpl]][,.lm1.idx,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  (Matrix::t(x = exprs.mat) %*% wcl)
                
                # get weighted mean using sum of weights
                tmp.res <-
                  Matrix::rowSums(x = wsum) / sum(wcl)
                
                return(tmp.res)
                
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
        stats::setNames(object = names(x = .lm.obj$cells),
                        nm = names(x = .lm.obj$cells)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .lm.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[,colnames(x = .lm.obj$lm)]
          
          res <-
            tryCatch(
              expr = {
                
                wcl <- 
                  .lm.obj$map$connect.prob[[smpl]][,.lm2.idx,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  (Matrix::t(x = exprs.mat) %*% wcl)
                
                # get weighted mean using sum of weights
                tmp.res <-
                  Matrix::rowSums(x = wsum) / sum(wcl)
                
                return(tmp.res)
                
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

#' Compute Sample Embedding from Partial Fitted Values
#'
#' Computes a low-dimensional sample embedding by learning PCA axes from the partial 
#' fitted component \eqn{\Delta\hat{Y} = (H_{full} - H_{red})Y}, then projecting
#' nuisance-residualized samples onto those axes. This isolates variation attributable 
#' to an effect of interest (e.g., a contrast or additional covariates) while preserving
#' residual scatter. Exact for OLS; if duplicateCorrelation/blocking is used, the 
#' decomposition is approximate.
#'
#' @param .lm.obj The tinydenseR object for accessing metadata and blocking information.
#' @param .stats.obj A stats object from \code{get.stats()} containing the full model fit.
#'   Returns a modified copy of .stats.obj with the embedding added.
#' @param .term.of.interest Character string: the covariate/term being isolated when using the 
#'   nested model method (\code{.red.stats.obj}). Must match a column name in 
#'   \code{.lm.obj$metadata}. Ignored when using \code{.contrast.of.interest}.
#' @param .red.stats.obj Optional: A stats object from \code{get.stats()} containing a 
#'   reduced model fit (without the effect of interest). If provided, computes 
#'   \eqn{\Delta\hat{Y} = \hat{Y}_{full} - \hat{Y}_{red}} directly.
#' @param .contrast.of.interest Optional: Character string naming the contrast to extract.
#'   Must match a column name in \code{.stats.obj$fit$contrasts}. When specified, uses
#'   the Frisch-Waugh-Lovell (FWL) theorem to compute the true partial fitted component.
#'   The embedding will be stored in \code{.stats.obj$embedding[[.contrast.of.interest]]}.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'
#' @return The modified \code{.stats.obj} with new embedding slot. For FWL method, the slot
#'   is named after the contrast; for nested model method, the slot is named after
#'   \code{.term.of.interest}. Each embedding is a flat list with:
#'   \describe{
#'     \item{\code{sdev}}{Standard deviations of PCs learned from \eqn{\Delta\hat{Y}^\top}}
#'     \item{\code{rotation}}{Loadings matrix (landmarks x PCs): the PCA basis \eqn{V} 
#'       learned from \eqn{\Delta\hat{Y}^\top}}
#'     \item{\code{center}}{Centering vector: column means of \eqn{\Delta\hat{Y}^\top}}
#'     \item{\code{scale}}{Scaling vector (FALSE since scale. = FALSE)}
#'     \item{\code{x}}{PCA scores on \eqn{\Delta\hat{Y}^\top} (samples x PCs): the partial-effect 
#'       embedding from prcomp}
#'     \item{\code{coord}}{Nuisance-residualized projection (samples x PCs): 
#'       \eqn{(E_{red}^\top - center) V} where \eqn{E_{red} = Y - \hat{Y}_{red}}.
#'       This is the main embedding for visualization--samples spread along PCs by their 
#'       alignment with the effect direction, plus residual scatter.}
#'     \item{\code{effect.resid.cor}}{Numeric vector: per-PC correlation between \code{x} 
#'       (partial-effect scores) and \code{coord} (residualized projection). High values 
#'       (>0.8) indicate residuals align well with the effect direction; low values suggest
#'       strong nuisance removal or weak signal.}
#'     \item{\code{delta.Yhat}}{Matrix (landmarks x samples) of partial fitted values
#'       \eqn{\Delta\hat{Y} = (H_{full} - H_{red})Y}}
#'     \item{\code{method}}{Character: "fwl_contrast" or "nested_models"}
#'   }
#'
#' @details
#' \strong{Two Methods for Computing \eqn{\Delta\hat{Y}}:}
#'
#' \strong{Method 1: FWL-based contrast extraction} (when \code{.contrast.of.interest} is provided)
#'
#' For a cell-means model with nuisance covariates, uses the Frisch-Waugh-Lovell theorem:
#' \enumerate{
#'   \item Compute contrast regressor: \eqn{x_c = X_{group} \cdot c}
#'   \item Residualize against nuisance: \eqn{x_{c,\perp} = (I - H_Z) x_c}
#'   \item Residualize Y against nuisance: \eqn{Y_\perp = (I - H_Z) Y}
#'   \item Compute partial coefficient: \eqn{\gamma_g = (x_{c,\perp}' Y_{\perp,g}) / (x_{c,\perp}' x_{c,\perp})}
#'   \item Form \eqn{\Delta\hat{Y} = \gamma_g \otimes x_{c,\perp}}
#' }
#'
#' \strong{Method 2: Nested model comparison} (when \code{.red.stats.obj} is provided)
#'
#' Direct computation: \eqn{\Delta\hat{Y} = \hat{Y}_{full} - \hat{Y}_{red}}
#'
#' This method is useful when you want to isolate the effect of one or more covariates
#' without using contrasts (e.g., comparing \code{~ disease + sex + age} vs \code{~ sex + age}
#' to extract the disease effect).
#'
#' \strong{Embedding via PCA:}
#'
#' The embedding is computed in two steps:
#' \enumerate{
#'   \item Learn PCA basis \eqn{V} from \eqn{\Delta\hat{Y}^\top} (partial fitted values):
#'     \code{pca <- prcomp(t(delta.Yhat), center = TRUE, scale. = FALSE, rank. = r)}
#'   \item Project nuisance-residualized samples onto this basis:
#'     \code{coord <- t(Y - Yhat.red - pca$center) \%*\% pca$rotation}
#' }
#'
#' Here \eqn{E_{red} = Y - \hat{Y}_{red}} are the nuisance-residualized features (what 
#' remains after removing nuisance effects). The loadings \eqn{V} encode the linear 
#' subspace (often rank-1 for a single contrast) that maximizes variance in the 
#' partial-effect matrix \eqn{\Delta\hat{Y}^\top}.
#'
#' Geometrically, this evaluates how the residualized samples align with the effect 
#' direction. Scores increase when residuals align with the partial-effect axis.
#' Scatter along PCs reflects residual alignment with the effect; scatter orthogonal 
#' to PCs is residual-only noise.
#'
#' @seealso \code{\link{get.stats}} for model fitting, \code{\link{plotSampleEmbedding}} for
#'   visualizing the embedding
#'
#' @examples
#' \dontrun{
#' # ===========================================================================
#' # Example 1: Nested model comparison (no contrasts)
#' # ===========================================================================
#' 
#' # Fit full model with disease effect
#' .design <- model.matrix(object = ~ disease + sex + age,
#'                         data = lm.obj$metadata)
#' 
#' stats.res <- tinydenseR::get.stats(
#'     .lm.obj = lm.obj,
#'     .design = .design,
#'     .verbose = TRUE
#' )
#' 
#' # Fit reduced model without disease
#' .red.design <- model.matrix(object = ~ sex + age,
#'                             data = lm.obj$metadata)
#' 
#' red.stats.res <- tinydenseR::get.stats(
#'     .lm.obj = lm.obj,
#'     .design = .red.design,
#'     .verbose = TRUE
#' )
#' 
#' # Extract disease effect embedding (slot named "disease")
#' stats.res <- get.embedding(
#'     .lm.obj = lm.obj,
#'     .stats.obj = stats.res,
#'     .term.of.interest = "disease",
#'     .red.stats.obj = red.stats.res
#' )
#' 
#' # Plot the embedding (colored by disease)
#' plotSampleEmbedding(lm.obj, stats.res, .embedding.slot = "disease")
#' 
#' # ===========================================================================
#' # Example 2: FWL-based contrast extraction (cell-means model)
#' # ===========================================================================
#' 
#' # Fit cell-means model with contrasts
#' design <- model.matrix(~ 0 + Group + Batch + Age, data = lm.obj$metadata)
#' contrasts <- limma::makeContrasts(
#'     TrtVsCtrl = GroupTrt - GroupCtrl,
#'     KO_TrtVsCtrl = GroupKO_Trt - GroupKO_Ctrl,
#'     levels = design
#' )
#' stats <- get.stats(lm.obj, .design = design, .contrasts = contrasts)
#'
#' # Extract each contrast embedding (slot named after contrast)
#' stats <- get.embedding(
#'     .lm.obj = lm.obj,
#'     .stats.obj = stats, 
#'     .contrast.of.interest = "TrtVsCtrl"
#' )
#'
#' stats <- get.embedding(
#'     .lm.obj = lm.obj,
#'     .stats.obj = stats, 
#'     .contrast.of.interest = "KO_TrtVsCtrl"
#' )
#'
#' # Plot the embedding (specify .color.by for metadata column)
#' plotSampleEmbedding(lm.obj, stats, 
#'                     .embedding.slot = "TrtVsCtrl", 
#'                     .color.by = "Group")
#'
#' # Access the embedding components (slot named after contrast)
#' stats$embedding$TrtVsCtrl$coord      # sample coordinates
#' stats$embedding$TrtVsCtrl$sdev       # PC standard deviations
#' stats$embedding$TrtVsCtrl$rotation   # landmark loadings
#' stats$embedding$TrtVsCtrl$delta.Yhat # partial fitted values
#'
#' # Multiple contrasts stored separately
#' names(stats$embedding)  # c("TrtVsCtrl", "KO_TrtVsCtrl")
#' }
#'
#' @export
#'
get.embedding <-
  function(
    .lm.obj,
    .stats.obj,
    .term.of.interest = NULL,
    .red.stats.obj = NULL,
    .contrast.of.interest = NULL,
    .verbose = TRUE
  ){
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    if(is.null(x = .red.stats.obj) &&
       is.null(x = .contrast.of.interest)){
      
      stop("Please provide either '.red.stats.obj' or '.contrast.of.interest'")
      
    }
    
    if(!is.null(x = .red.stats.obj) &&
       !is.null(x = .contrast.of.interest)){
      
      warning("Both '.red.stats.obj' and '.contrast.of.interest' provided. Using '.contrast.of.interest' (FWL method).")
      
    }
    
    # Determine slot name based on method
    if(!is.null(x = .contrast.of.interest)){
      
      # FWL method: use contrast name as slot
      slot.name <- .contrast.of.interest
      
    } else {
      
      # Nested model method: require .term.of.interest
      if(is.null(x = .term.of.interest) ||
         !is.character(x = .term.of.interest) || 
         length(x = .term.of.interest) != 1){
        
        stop("'.term.of.interest' must be a single character string when using nested model method")
        
      }
      if(!(.term.of.interest %in% colnames(x = .lm.obj$metadata))){
        
        stop("'.term.of.interest' must be a column name of '.lm.obj$metadata'")
        
      }
      
      slot.name <- .term.of.interest
      
    }
    
    # Initialize embedding slot if it doesn't exist
    if(is.null(x = .stats.obj$embedding)){
      
      .stats.obj$embedding <- 
        vector(mode = "list",
               length = 0)
      
    }
    
    # -------------------------------------------------------------------------
    # METHOD 1: FWL-based contrast extraction
    # -------------------------------------------------------------------------
    
    if(!is.null(x = .contrast.of.interest)){
      
      if(is.null(x = .stats.obj$fit$contrasts)){
        stop("No contrasts found in '.stats.obj'. Run get.stats() with .contrasts argument.")
      }
      
      if(!(.contrast.of.interest %in% 
           colnames(x = .stats.obj$fit$contrasts))){
        stop(paste0("'",
                    .contrast.of.interest, 
                    "' not found in contrast names: ",
                    paste(colnames(x = .stats.obj$fit$contrasts),
                          collapse = ", ")))
      }
      
      if(isTRUE(x = .verbose)){
        message("\nComputing embedding via FWL decomposition")
        message("  Slot name: ", 
                slot.name)
        message("  Contrast: ",
                .contrast.of.interest)
      }
      
      # Get the contrast vector
      contrast.vec <- 
        .stats.obj$fit$contrasts[, .contrast.of.interest]
      
      # Identify which design columns are involved in contrasts vs nuisance
      contrast.participation <- 
        Matrix::rowSums(x = abs(x = .stats.obj$fit$contrasts)) > 0
      
      .full.design <- 
        .stats.obj$fit$design
      
      group.cols <-
        which(x = contrast.participation)
      
      nuisance.cols <-
        which(x = !contrast.participation)
      
      if(isTRUE(x = .verbose)){
        
        message("  Group columns: ", 
                paste(colnames(x = .full.design)[group.cols], 
                      collapse = ", "))
        
        if(length(nuisance.cols) > 0){
          
          message("  Nuisance columns: ", 
                  paste(colnames(x = .full.design)[nuisance.cols],
                        collapse = ", "))
          
        } else {
          
          message("  No nuisance columns detected")
          
        }
      }
      
      # Get the expression matrix Y (landmarks x samples)
      Y <- 
        .stats.obj$y
      
      G <- 
        nrow(x = Y)
      
      # -----------------------------------------------------------------------
      # Step 1: Compute contrast regressor x_c = X_group %*% c
      # -----------------------------------------------------------------------
      X.group <- 
        .full.design[, group.cols, drop = FALSE]
      
      x_c <- 
        as.vector(x = X.group %*% contrast.vec[group.cols])
      
      # -----------------------------------------------------------------------
      # Step 2 & 3: FWL residualization (if nuisance covariates exist)
      # -----------------------------------------------------------------------
      
      if(length(x = nuisance.cols) > 0){
        
        Z <- 
          .full.design[, nuisance.cols, drop = FALSE]
        
        has.intercept <-
          apply(X = Z, 
                MARGIN = 2, 
                FUN = function(col) all(col == 1)) |>
          any()
        
        if(!has.intercept){
          Z <- cbind("(Intercept)" = 1, Z)
        }
        
        if(isTRUE(x = .verbose)){
          message("\nResidualizing against nuisance covariates")
        }
        
        # Handle blocking if present
        .block <- NULL
        dupcor <- NULL
        
        if(!is.null(x = .stats.obj$fit$block)){
          
          .block <- 
            lapply(X = .lm.obj$metadata,
                   FUN = function(meta.col){
                     if(length(meta.col) != length(.stats.obj$fit$block)) return(FALSE)
                     all(as.character(x = meta.col) == as.character(x = .stats.obj$fit$block))
                   }) |>
            unlist(use.names = TRUE) |>
            which() |>
            names()
          
          if(length(.block) > 0){
            .block <- .block[1]
            
            if(isTRUE(x = .verbose)){
              message("  Blocking variable: ", .block)
            }
            
            dupcor <- 
              limma::duplicateCorrelation(object = Y,
                                          design = Z,
                                          block = .lm.obj$metadata[[.block]])
          }
        }
        
        # FWL Step 1: Residualize x_c against Z
        fit.xc.on.Z <-
          stats::lm.fit(x = Z, 
                        y = x_c)
        
        x_c_perp <-
          stats::residuals(object = fit.xc.on.Z)
        
        # FWL Step 2: Residualize Y against Z
        fit.Y.on.Z <- 
          limma::lmFit(object = Y,
                       design = Z,
                       block = if(!is.null(dupcor)) .lm.obj$metadata[[.block]] else NULL,
                       correlation = if(!is.null(dupcor)) dupcor$consensus else NULL)
        
        # Yhat.red = fitted values from nuisance-only model (used for projection)
        Yhat.red <- 
          Matrix::tcrossprod(x = fit.Y.on.Z$coefficients, 
                             y = Z)
        
        Y_perp <-
          Y - Yhat.red
        
      } else {
        
        # No nuisance covariates: center only
        if(isTRUE(x = .verbose)){
          message("\nNo nuisance covariates: centering only")
        }
        
        x_c_perp <- 
          x_c - mean(x = x_c)
        
        # Yhat.red = grand mean (intercept-only model)
        Yhat.red <-
          matrix(data = Matrix::rowMeans(x = Y),
                 nrow = G,
                 ncol = ncol(x = Y),
                 dimnames = dimnames(x = Y))
        
        Y_perp <-
          Y - Yhat.red
        
      }
      
      # -----------------------------------------------------------------------
      # Step 3: Compute partial regression coefficient gamma_g
      # -----------------------------------------------------------------------
      
      if(isTRUE(x = .verbose)){
        message("Computing partial regression coefficients")
      }
      
      denom <- sum(x_c_perp^2)
      
      if(denom < .Machine$double.eps * 100){
        
        warning("x_c_perp has near-zero variance. Contrast may be aliased with nuisance covariates.")
        
        gamma_g <- 
          rep(x = 0,
              times = G)
        
      } else {
        
        gamma_g <-
          as.vector(x = Y_perp %*% x_c_perp) / denom
        
      }
      
      # -----------------------------------------------------------------------
      # Step 4: Compute delta.Yhat = gamma_g (outer) x_c_perp
      # -----------------------------------------------------------------------
      
      delta.Yhat <- 
        Matrix::tcrossprod(x = matrix(gamma_g, ncol = 1),
                                       y = matrix(x_c_perp, ncol = 1))
      
      rownames(x = delta.Yhat) <- 
        rownames(x = Y)
      
      colnames(x = delta.Yhat) <- 
        colnames(x = Y)
      
      method <- "fwl_contrast"
      
    } else {
      
      # -----------------------------------------------------------------------
      # METHOD 2: Direct nested model comparison
      # -----------------------------------------------------------------------
      
      if(!is.null(x = .stats.obj$fit$contrasts)){
        stop("Contrasts found in '.stats.obj'. That means, coefficients are 're-oriented' via contrasts.fit. Use '.contrast.of.interest' instead.")
      }
      
      if(isTRUE(x = .verbose)){
        message("\nComputing embedding via nested model comparison")
        message("  Slot name: ", slot.name)
        message("  Full model columns: ", paste(colnames(.stats.obj$fit$design), collapse = ", "))
        message("  Reduced model columns: ", paste(colnames(.red.stats.obj$fit$design), collapse = ", "))
      }
      
      # Validate dimensions
      if(nrow(.stats.obj$y) != nrow(.red.stats.obj$y)){
        stop("Number of landmarks differs between full and reduced models")
      }
      if(ncol(.stats.obj$y) != ncol(.red.stats.obj$y)){
        stop("Number of samples differs between full and reduced models")
      }
      
      if(!identical(colnames(.stats.obj$y), colnames(.red.stats.obj$y))){
        warning("Sample names differ between models. Assuming same order.")
      }
      
      # Compute fitted values for each model
      Yhat.full <- 
        Matrix::tcrossprod(x = .stats.obj$fit$coefficients, 
                                      y = .stats.obj$fit$design)
      
      Yhat.red <-
        Matrix::tcrossprod(x = .red.stats.obj$fit$coefficients,
                                     y = .red.stats.obj$fit$design)
      
      # delta.Yhat = (H_full - H_red)Y
      delta.Yhat <- 
        Yhat.full - Yhat.red
      
      rownames(x = delta.Yhat) <- 
        rownames(x = .stats.obj$y)
      
      colnames(x = delta.Yhat) <- 
        colnames(x = .stats.obj$y)
      
      method <- "nested_models"
      
    }
    
    # -------------------------------------------------------------------------
    # Learn axis of maximal of variance along of variable interest from delta.Yhat
    # -------------------------------------------------------------------------
    
    if(isTRUE(x = .verbose)){
      message("\nComputing PCA embedding")
    }
    
    # Determine rank of delta.Yhat
    delta.Yhat.t <- 
      Matrix::t(x = delta.Yhat)
    
    embed.rank <-
      (Matrix::qr(x = delta.Yhat.t)$rank) |> 
      (\(x)
       max(1, 
           min(x, 
               nrow(x = delta.Yhat.t)-1, 
               ncol(x = delta.Yhat.t)))
      )()
    
    
    if(isTRUE(x = .verbose)){
      message("  Rank of delta.Yhat: ", embed.rank)
    }
    
    # Compute PCA
    pca <- 
      stats::prcomp(x = as.matrix(x = delta.Yhat.t), 
                    center = TRUE, 
                    scale. = FALSE,
                    rank. = embed.rank)
    
    #if(isTRUE(x = .verbose)){
    #  var.explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    #  message("  Variance explained by PCs: ", 
    #          paste(paste0("PC", seq_along(var.explained), "=", var.explained, "%"), collapse = ", "))
    #}
    
    # -------------------------------------------------------------------------
    # Get projection of residualized samples onto the partial-effect axis(-es)
    # -------------------------------------------------------------------------
    
    # project E.red.t using delta.Yhat.t's V
    pca$coord <-
      Matrix::t(x = (.stats.obj$y - Yhat.red) - pca$center) %*%
      pca$rotation
    
    # -------------------------------------------------------------------------
    # Compute per-PC correlation: effect scores (x) vs residualized projection (coord)
    # -------------------------------------------------------------------------
    
    effect.resid.cor <- 
      vapply(X = seq_len(ncol(x = pca$rotation)),
             FUN = function(j) {
               stats::cor(x = pca$x[, j], 
                          y = pca$coord[, j])
             },
             FUN.VALUE = numeric(length = 1))
    
    names(x = effect.resid.cor) <- 
      paste0("PC", seq_along(along.with = effect.resid.cor))
    
    if(isTRUE(x = .verbose)){
      message("\nEffect-residual alignment (per-PC correlation):")
      for(j in seq_along(along.with = effect.resid.cor)){
        message(sprintf("  PC%d: %.3f", j, effect.resid.cor[j]))
      }
    }
    
    # -------------------------------------------------------------------------
    # Store results in .stats.obj$embedding[[slot.name]]
    # -------------------------------------------------------------------------
    
    # Flatten PCA result into shallow list
    #pca$coord <- pca$x
    #pca$x <- NULL
    
    .stats.obj$embedding[[slot.name]] <- 
      c(pca,
        list(effect.resid.cor = effect.resid.cor,
             delta.Yhat = as.matrix(x = delta.Yhat),
             method = method)
      )
    
    if(isTRUE(x = .verbose)){
      message("\nEmbedding stored in: .stats.obj$embedding$", slot.name)
    }
    
    return(.stats.obj)
    
  }
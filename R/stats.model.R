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

#' Differential Density Testing
#'
#' Performs landmark-based differential density testing using limma on size factor-normalized and 
#' log-transformed fuzzy densities. Tests which landmarks (cell states) change in density
#' between conditions. More sensitive than traditional cluster-level testing because
#' landmarks can capture within-cluster heterogeneity. Uses PCA-weighted q-values that 
#' leverage the correlation structure among landmarks to improve statistical power.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .design Design matrix specifying the experimental design. Rows correspond to samples 
#'   (matching \code{.tdr.obj$cells}), columns to coefficients. Create with \code{model.matrix()}.
#' @param .contrasts Optional contrast matrix for specific comparisons. Each column defines one 
#'   contrast. Create with \code{limma::makeContrasts()}. If NULL, tests all coefficients in 
#'   \code{.design}.
#' @param .block Optional character: column name in \code{.tdr.obj$metadata} for blocking factor 
#'   (e.g., "Donor", "Batch"). Accounts for within-block correlation using \code{duplicateCorrelation}.
#' @param .model.name Character string naming this model fit (default "default"). Results are stored 
#'   in \code{.tdr.obj$map$lm[[.model.name]]}. Use different names to store multiple model fits 
#'   (e.g., full vs reduced models for nested model comparisons).
#' @param .force.recalc Logical: if TRUE, overwrite existing results in the specified slot 
#'   (default FALSE). If FALSE and slot already exists, an error is thrown.
#' @param .verbose Logical: print progress messages? Default TRUE.
#' @param .seed Integer: random seed for reproducibility in blocking correlation estimation. Default 123.
#'   
#' @return The modified \code{.tdr.obj} with results stored in \code{.tdr.obj$map$lm[[.model.name]]}:
#'   \describe{
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
#' lm.obj <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#' 
#' # Simple two-group comparison (results stored in lm.obj$map$lm$default)
#' design <- model.matrix(~ Condition, data = .meta)
#' lm.obj <- get.lm(lm.obj, .design = design)
#' 
#' # Visualize results
#' plotBeeswarm(lm.obj, .coefs = "ConditionB")
#' 
#' # Complex design with contrasts
#' design <- model.matrix(~ 0 + Group, data = .meta)
#' contrasts <- limma::makeContrasts(
#'   TvsC = GroupTreatment - GroupControl,
#'   T1vsT2 = GroupTreatment1 - GroupTreatment2,
#'   levels = design
#' )
#' lm.obj <- get.lm(lm.obj, .design = design, .contrasts = contrasts,
#'                     .model.name = "full", .force.recalc = TRUE)
#' 
#' # With blocking for paired samples
#' design <- model.matrix(~ Timepoint, data = .meta)
#' lm.obj <- get.lm(lm.obj, .design = design, .block = "Subject",
#'                     .model.name = "timepoint")
#' 
#' # Nested model comparison: fit reduced model
#' red.design <- model.matrix(~ 1, data = .meta)  # intercept only
#' lm.obj <- get.lm(lm.obj, .design = red.design, .model.name = "reduced")
#' 
#' # Access results by model name
#' lm.obj$map$lm$full$fit$coefficients
#' lm.obj$map$lm$reduced$fit$coefficients
#' }
#' 
#' @export
#'
get.lm <-
  function(
    .tdr.obj,
    .design,
    .contrasts = NULL,
    .block = NULL,
    .model.name = "default",
    .force.recalc = FALSE,
    .verbose = TRUE,
    .seed = 123){
    
    # -------------------------------------------------------------------------
    # Validate on-disk cache (if active) before proceeding
    # -------------------------------------------------------------------------
    
    .tdr_cache_validate_quiet(.tdr.obj)
    
    # -------------------------------------------------------------------------
    # Check if slot already exists
    # -------------------------------------------------------------------------
    
    if(!is.null(x = .tdr.obj$map$lm[[.model.name]]) && !isTRUE(x = .force.recalc)){
      stop(paste0("Model '", .model.name, "' already exists in .tdr.obj$map$lm. ",
                  "Use a different .model.name or set .force.recalc = TRUE to overwrite."))
    }
    
    if(nrow(x = .design) != length(x = .tdr.obj$cells)){
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
    nlib <- length(x = .tdr.obj$cells)
    .design <- as.matrix(x = .design)
    if(nrow(x = .design) != nlib) {
      stop("nrow(design) disagrees with length(x = .tdr.obj$cells)")
    }
    ne <- limma::nonEstimable(x = .design)
    if(!is.null(x = ne)) {
      stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n",
                 paste(ne,
                       collapse = " ")))
    }
    
    if(is.null(x = .tdr.obj$map)){
      stop("First run get.map")
    }
    
    # create list to hold results
    stats <-
      vector(mode = "list",
             length = 0)
    
    # Use pre-computed Y from get.map
    Y <-
      .tdr.obj$map$Y
    
    #if(nrow(x = Y) !=
    #   nrow(x = .tdr.obj$landmarks)){
    #  
    #  stats$y <-
    #    stats$y[match(x = rownames(x = .tdr.obj$landmarks),
    #                  table = rownames(x = stats$y)),]
    #  
    #  stats$y[
    #    is.na(x = stats$y)
    #  ] <- min(stats$y,
    #           na.rm = TRUE)
    #  
    #  rownames(x = stats$y) <-
    #    rownames(x = .tdr.obj$landmarks)
    #  
    #}
    
    if(!(is.null(x = .block))){
      
      if(length(x = .block) != 1){
        stop("Block must be a vector of length 1")
      }
      if(!(.block %in% colnames(x = .tdr.obj$metadata))){
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
        limma::duplicateCorrelation(object = Y,
                                    design = .design,
                                    block = .tdr.obj$metadata[[.block]])
    }
    
    if(isTRUE(x = .verbose)){
      message("\nfitting linear models")
    }
    
    stats$fit <-
      limma::lmFit(object = Y,
                   design = .design,
                   block = if(exists(x = "dupcor")) .tdr.obj$metadata[[.block]] else  NULL,
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
              w <- 1 / log10(x = Matrix::rowSums(x = .tdr.obj$map$fdens) + 1)
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
        Matrix::t(x = .tdr.obj$pca$embed)
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .tdr.obj$metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        # https://support.bioconductor.org/p/125489/#125602
        # duplicateCorrelation is more general and is THE ONLY SOLUTION when
        # you want to compare across blocking levels, e.g., comparing diseased
        # and healthy donors when each donor also contributes before/after treatment samples.
        dupcor.q <- 
          limma::duplicateCorrelation(object = tX,
                                      design = stats$fit$design[.tdr.obj$config$key,],
                                      block = .tdr.obj$metadata[[.block]][.tdr.obj$config$key])
      }
      
      fit.q <-
        limma::lmFit(object = tX,
                     design = stats$fit$design[.tdr.obj$config$key,],
                     block = if(exists(x = "dupcor.q")) .tdr.obj$metadata[[.block]][.tdr.obj$config$key] else  NULL,
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
         x[,1:elbow.sec.deriv(x = .tdr.obj$pca$sdev^2,
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
      if(!(.block %in% colnames(x = .tdr.obj$metadata))){
        stop(paste0(.block,
                    " not found in metadata"))
      }
      
      if(isTRUE(x = .verbose)){
        message("\nestimating the intra-block correlation for stats from clustering")
      }
      
      cl.dupcor <- 
        log2(x = .tdr.obj$map$clustering$cell.perc + 0.5) |>
        Matrix::t() |>
        limma::duplicateCorrelation(design = .design,
                                    block = .tdr.obj$metadata[[.block]])
    }
    
    if(isTRUE(x = .verbose)){
      message("\nfitting linear models")
    }
    
    stats$trad$clustering$fit <-
      limma::lmFit(object = log2(x = .tdr.obj$map$clustering$cell.perc + 0.5) |>
                     Matrix::t(),
                   design = .design,
                   block = if(exists(x = "cl.dupcor")) .tdr.obj$metadata[[.block]] else NULL,
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
    
    if(!is.null(x = .tdr.obj$graph$celltyping)){
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .tdr.obj$metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        if(isTRUE(x = .verbose)){
          message("\nestimating the intra-block correlation for stats from celltyping")
        }
        
        ct.dupcor <- 
          log2(x = .tdr.obj$map$celltyping$cell.perc + 0.5) |>
          Matrix::t() |>
          limma::duplicateCorrelation(design = .design,
                                      block = .tdr.obj$metadata[[.block]])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      stats$trad$celltyping$fit <-
        limma::lmFit(object = log2(x = .tdr.obj$map$celltyping$cell.perc + 0.5) |>
                       Matrix::t(),
                     design = .design,
                     block = if(exists(x = "ct.dupcor")) .tdr.obj$metadata[[.block]] else NULL,
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
    
    # -------------------------------------------------------------------------
    # Store results and return modified .tdr.obj
    # -------------------------------------------------------------------------
    
    # Initialize stats slot if needed
    if(is.null(x = .tdr.obj$map$lm)){
      .tdr.obj$map$lm <- list()
    }
    
    # Store results under the model name
    .tdr.obj$map$lm[[.model.name]] <- stats
    
    if(isTRUE(x = .verbose)){
      message("\nResults stored in: .tdr.obj$map$lm$", .model.name)
    }
    
    return(.tdr.obj)
    
  }

#' Pseudobulk Differential Expression Analysis
#'
#' Performs pseudobulk differential expression (DE) analysis for genes/markers. Aggregates cells 
#' into pseudobulk samples using fuzzy landmark membership, then uses limma-voom (RNA) or limma 
#' (cytometry) to test for DE. Can restrict analysis to specific clusters/cell types, and optionally 
#' perform gene set enrichment via GSVA.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .design Design matrix specifying experimental design. Rows = samples, columns = coefficients.
#' @param .contrasts Optional contrast matrix for specific comparisons. Create with 
#'   \code{limma::makeContrasts()}. If NULL, tests all \code{.design} coefficients.
#' @param .block Optional character: column name in \code{.tdr.obj$metadata} for blocking factor 
#'   (e.g., "Donor"). Accounts for within-block correlation.
#' @param .geneset.ls Optional named list of character vectors defining gene sets for GSVA enrichment 
#'   analysis. Only for RNA data. Example: \code{list("Tcell" = c("CD3D", "CD3E"), "Bcell" = c("CD19", "MS4A1"))}.
#' @param .id.idx Optional integer vector specifying landmark indices to include. If provided,
#'   cells are assigned to these landmarks using fuzzy membership confidence (see \code{.label.confidence}).
#'   Only cells whose fuzzy mass ratio to the target landmarks meets the confidence threshold are included.
#' @param .id Optional character vector of cluster/celltype IDs to restrict analysis to. Uses 
#'   \code{.id.from} to determine source.
#' @param .id.from Character: "clustering" or "celltyping". Source of IDs when \code{.id} is specified.
#' @param .model.name Character string naming this model fit (default "default"). Results are stored 
#'   in \code{.tdr.obj$pbDE[[.model.name]][[.population.name]]}.
#' @param .population.name Character string naming this population (default "all"). If \code{.id} is 
#'   specified and \code{.population.name} is not changed from default, it is auto-set to the value 
#'   of \code{.id} (or concatenated with "_" if multiple IDs). Enables storing multiple 
#'   population-specific DE results under the same model.
#' @param .force.recalc Logical: if TRUE, overwrite existing results in the specified slot 
#'   (default FALSE). If FALSE and slot already exists, an error is thrown.
#' @param .verbose Logical: print progress messages? Default TRUE.
#' @param .label.confidence Numeric scalar in \code{[0.5,1]} controlling the minimum posterior
#'   confidence required to assign a cell to a set of target landmarks (used when \code{.id.idx} is
#'   provided). For each cell, the ratio of fuzzy membership mass to the target landmarks vs total
#'   fuzzy mass is computed; cells below this threshold are excluded from aggregation. Default 0.8.
#'   
#' @return The modified \code{.tdr.obj} with results stored in \code{.tdr.obj$pbDE[[.model.name]][[.population.name]]}:
#'   \describe{
#'     \item{\code{fit}}{limma MArrayLM fit object with moderated statistics}
#'     \item{\code{coefficients}}{Log fold change matrix (features x coefficients)}
#'     \item{\code{p.value}}{P-values (features x coefficients)}
#'     \item{\code{adj.p}}{FDR-adjusted p-values (features x coefficients)}
#'     \item{\code{smpl.outlier}}{Logical vector indicating which samples were excluded as outliers}
#'     \item{\code{lm.idx}}{Integer vector of landmark indices used for aggregation}
#'     \item{\code{n.pseudo}}{Integer vector of pseudobulk cell counts per sample}
#'     \item{\code{geneset}}{(RNA only, if \code{.geneset.ls} provided) GSVA results with \code{$fit} 
#'       (limma fit), \code{$E} (enrichment scores), \code{$adj.p} (adjusted p-values)}
#'   }
#'   
#' @details
#' The pseudobulk DE workflow:
#' \enumerate{
#'   \item \strong{Cell selection}: If \code{.id} specified, select cells whose transferred label
#'     matches the requested population. If \code{.id.idx} specified, select cells whose fuzzy
#'     membership confidence to the target landmarks meets \code{.label.confidence}. Otherwise
#'     use all cells.
#'   \item \strong{Pseudobulk aggregation}: For each sample, compute weighted expression of the
#'     selected cells using their fuzzy membership across all landmarks:
#'     \itemize{
#'       \item RNA: Weighted sum of counts, scaled by effective cell count
#'       \item Cytometry: Weighted mean of marker expression
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
#' \strong{Cell-centric aggregation:}
#' Cells are first selected into the population of interest, then their expression is aggregated
#' using fuzzy membership weights across all landmarks. This approach first determines \emph{which}
#' cells belong to the population, then lets the fuzzy graph determine \emph{how} each selected
#' cell's expression is distributed across landmark neighborhoods.
#' 
#' \strong{When to use pseudobulk DE:}
#' \itemize{
#'   \item Testing treatment effects on gene expression within populations
#'   \item Comparing expression between conditions using samples as biological replicates
#'   \item Respects the true experimental design (no pseudoreplication)
#' }
#' 
#' @seealso \code{\link{get.map}} (required), \code{\link{plotPbDE}} for heatmap visualization,
#'   \code{\link{get.marker}} for marker gene identification, \code{\link{get.specDE}} for 
#'   density contrast-coupled expression programs
#'   
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#' 
#' # DE analysis on all cells (stored in $pbDE$default$all)
#' design <- model.matrix(~ Condition, data = .meta)
#' lm.cells <- get.pbDE(lm.cells, .design = design)
#' plotPbDE(lm.cells, .coefs = "ConditionB")
#' 
#' # DE within specific cell type (stored in $pbDE$default$tcells)
#' lm.cells <- get.pbDE(lm.cells, .design = design,
#'                      .id = c("1", "2", "3"),
#'                      .id.from = "clustering",
#'                      .population.name = "tcells")
#' 
#' # With gene set enrichment
#' hallmark.sets <- list(
#'   "INTERFERON_RESPONSE" = c("ISG15", "ISG20", "IFIT1"),
#'   "INFLAMMATORY_RESPONSE" = c("IL1B", "TNF", "CXCL8")
#' )
#' lm.cells <- get.pbDE(lm.cells, .design = design,
#'                      .geneset.ls = hallmark.sets,
#'                      .model.name = "gsva_model")
#' 
#' # Access results
#' lm.cells$pbDE$default$all$adj.p
#' lm.cells$pbDE$default$tcells$coefficients
#' }
#' 
#' @export
#'
get.pbDE <-
  function(
    .tdr.obj,
    .design,
    .contrasts = NULL,
    .block = NULL,
    .geneset.ls = NULL,
    .id.idx = NULL,
    .id = NULL,
    .id.from = NULL,
    .model.name = "default",
    .population.name = "all",
    .force.recalc = FALSE,
    .verbose = TRUE,
    .label.confidence = 0.8
  ) {
    
    i <- j <- landmark <- cell <- label <- x <- confidence <-
      NULL
    
    if(.label.confidence < 0.5 || .label.confidence > 1){
      stop(".label.confidence must be between 0.5 and 1")
    }
    
    if(!is.null(x = .geneset.ls)){
      if(.tdr.obj$config$assay.type != "RNA"){
        stop(".geneset.ls is only supported for RNA assay type")
      } else if(!is.list(x = .geneset.ls)){
        stop(".geneset.ls must be a list of character vectors")
      }
    } 
    
    if(nrow(x = .design) != length(x = .tdr.obj$cells)){
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
    nlib <- length(x = .tdr.obj$cells)
    .design <- as.matrix(x = .design)
    if(nrow(x = .design) != nlib) {
      stop("nrow(design) disagrees with length(x = .tdr.obj$cells)")
    }
    ne <- limma::nonEstimable(x = .design)
    if(!is.null(x = ne)) {
      stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n",
                 paste(ne,
                       collapse = " ")))
    }
    
    # -------------------------------------------------------------------------
    # Auto-set .population.name if .id is specified and .population.name is default
    # -------------------------------------------------------------------------
    
    if(!is.null(x = .id) && .population.name == "all"){
      # Auto-set population name to the specified cluster/celltype ID(s)
      .population.name <- paste(.id, collapse = "_")
      if(isTRUE(x = .verbose)){
        message(sprintf("Auto-setting .population.name to '%s'", .population.name))
      }
    }
    
    if(!is.null(x = .id.idx) && .population.name == "all"){
      # Custom indices require explicit naming
      .population.name <- paste0("idx_", length(x = .id.idx))
      if(isTRUE(x = .verbose)){
        message(sprintf("Auto-setting .population.name to '%s' (custom landmark indices)", .population.name))
      }
    }
    
    # -------------------------------------------------------------------------
    # Validate on-disk cache (if active) before proceeding
    # -------------------------------------------------------------------------
    
    .tdr_cache_validate_quiet(.tdr.obj)
    
    # -------------------------------------------------------------------------
    # Check if slot already exists
    # -------------------------------------------------------------------------
    
    if(!is.null(x = .tdr.obj$pbDE[[.model.name]][[.population.name]]) && !isTRUE(x = .force.recalc)){
      stop(paste0("Results for model '", .model.name, "' and population '", .population.name, 
                  "' already exist in .tdr.obj$pbDE. ",
                  "Use different names or set .force.recalc = TRUE to overwrite."))
    }
    
    # -------------------------------------------------------------------------
    # Determine cell indices for aggregation (cell-centric approach)
    # .lm.idx is a per-sample list of cell indices
    # -------------------------------------------------------------------------
    
    if(is.null(x = .id.idx)){
      
      if(!is.null(x = .id)){
        
        .id.from <-
          match.arg(arg = .id.from,
                    choices = c("clustering",
                                "celltyping"))
        
        if(!all(.id %in% unique(x = .tdr.obj$graph[[.id.from]]$ids))){
          
          stop(paste0(paste0(.id[!(.id %in% unique(x = .tdr.obj$graph[[.id.from]]$ids))],
                             collapse = ", "),
                      " not found in ",
                      .id.from))
          
        }
        
        # Cell-centric: get per-sample cell indices matching the requested population
        .lm.idx <-
          lapply(X = .tdr.obj$map[[.id.from]]$ids,
                 FUN = function(smpl){
                   which(x = smpl %in% .id)
                 })
        
      } else {
        
        # Use all cells per sample
        .lm.idx <-
          stats::setNames(object = .tdr.obj$metadata$n.cells,
                          nm = names(x = .tdr.obj$cells)) |>
          lapply(FUN = seq_len)
        
      }
    } else {
      
      if(!all(.id.idx %in% (nrow(x = .tdr.obj$landmarks) |> seq_len()))) {
        stop(paste0(".id.idx must be an integer vector between 1 and ",
                    nrow(x = .tdr.obj$landmarks)))
      }
      
      # Cell-centric: use .label.confidence to determine which cells
      # confidently belong to the specified landmark set
      tmp.lbl <-
        rep(x = "out",
            times = nrow(x = .tdr.obj$landmarks))
      
      tmp.lbl[.id.idx] <-
        "in"
      
      .lm.idx <-
        lapply(X = .tdr_get_map_slot_all(.tdr.obj, "fuzzy.graph"),
               FUN = function(smpl) {
                 
                 in.and.out <-
                   Matrix::summary(object = smpl) |>
                   dplyr::rename(cell = i,
                                 landmark = j) |>
                   dplyr::mutate(label = tmp.lbl[landmark]) |>
                   dplyr::select(-landmark) |>
                   collapse::fgroup_by(cell,
                                       label,
                                       sort = FALSE) |>
                   collapse::fsum() |>
                   collapse::fgroup_by(cell,
                                       sort = FALSE) |>
                   collapse::fmutate(confidence = x / collapse::fsum(x)) |>
                   collapse::fungroup() |>
                   collapse::fsubset(confidence >= .label.confidence) |>
                   (\(x)
                    {
                      label <-
                        rep(x = "..low.confidence..",
                            times = nrow(x = smpl))
                      
                      label[x$cell] <-
                        x$label
                      
                      label
                      
                   }
                   )()
                 
                 which(x = in.and.out == "in")
                 
               })
      
    }
    
    if(isTRUE(x = .verbose)){
      message(sprintf("\nUsing %s cells for pseudobulk aggregation",
                      paste(lengths(x = .lm.idx), collapse = ", ")))
    }
    
    # number of cells in each pseudobulk
    n.pseudo <-
      lengths(x = .lm.idx)
    
    # remove outlier samples with too few cells
    smpl.outlier <-
      (n.pseudo / mean(x = n.pseudo)) < 0.1
    
    # get pseudobulk
    if(.tdr.obj$config$assay.type == "RNA"){
      
      if(any(smpl.outlier)){
        warning(paste0("The following samples were removed since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
                       paste(names(x = smpl.outlier)[smpl.outlier],
                             collapse = "\n")))
      }
      
      counts <-
        stats::setNames(object = names(x = .tdr.obj$cells)[!smpl.outlier],
                        nm = names(x = .tdr.obj$cells)[!smpl.outlier]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .tdr.obj$cells[[smpl]])[,.lm.idx[[smpl]]]
          
          wcl <- 
            .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", smpl)[.lm.idx[[smpl]],,drop = FALSE]
          
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
        if(!(.block %in% colnames(x = .tdr.obj$metadata))){
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
                                      block = .tdr.obj$metadata[[.block]][!smpl.outlier])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      .de <-
        limma::lmFit(object = v,
                     design = tmp.design,
                     block = if(exists(x = "dupcor")) .tdr.obj$metadata[[.block]][!smpl.outlier] else  NULL,
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
                                        block = .tdr.obj$metadata[[.block]][!smpl.outlier])
        }
        
        gsva.fit <- 
          limma::lmFit(object = gsva.es,
                       design = tmp.design,
                       block = if(exists(x = "gsva.dupcor")) .tdr.obj$metadata[[.block]][!smpl.outlier] else  NULL,
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
        stats::setNames(object = names(x = .tdr.obj$cells),#[!smpl.outlier],
                        nm = names(x = .tdr.obj$cells)) |>#[!smpl.outlier]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .tdr.obj$cells[[smpl]])
          
          exprs.mat <-
            exprs.mat[.lm.idx[[smpl]],colnames(x = .tdr.obj$landmarks)]
          
          wcl <- 
            .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", smpl)[.lm.idx[[smpl]],,drop=FALSE]
          
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
        if(!(.block %in% colnames(x = .tdr.obj$metadata))){
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
                                      block = .tdr.obj$metadata[[.block]])#[!smpl.outlier])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      .de <-
        limma::lmFit(object = Matrix::t(x = counts),
                     design = tmp.design,
                     block = if(exists(x = "dupcor")) .tdr.obj$metadata[[.block]] else  NULL,#[!smpl.outlier] else  NULL,
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
      .lm.idx
    
    .de$n.pseudo <-
      n.pseudo
    
    # -------------------------------------------------------------------------
    # Store results in .tdr.obj$pbDE[[.model.name]][[.population.name]]
    # -------------------------------------------------------------------------
    
    if(is.null(x = .tdr.obj$pbDE)){
      .tdr.obj$pbDE <- list()
    }
    
    if(is.null(x = .tdr.obj$pbDE[[.model.name]])){
      .tdr.obj$pbDE[[.model.name]] <- list()
    }
    
    .tdr.obj$pbDE[[.model.name]][[.population.name]] <- .de
    
    if(isTRUE(x = .verbose)){
      message(sprintf("\nResults stored in .tdr.obj$pbDE$%s$%s", .model.name, .population.name))
    }
    
    return(.tdr.obj)
    
  }

#' Pseudobulk Differential Expression Analysis (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' \code{get.dea()} has been renamed to \code{\link{get.pbDE}()} for clarity. This function

#' is provided for backward compatibility and will be removed in a future version.
#'
#' @inheritParams get.pbDE
#' @return A list containing DE analysis results (legacy format). Use \code{get.pbDE()} instead
#'   which stores results in \code{.tdr.obj$pbDE} and returns the modified object.
#'   
#' @seealso \code{\link{get.pbDE}}
#' @keywords internal
#' @export
get.dea <- function(
    .tdr.obj,
    .design,
    .contrasts = NULL,
    .block = NULL,
    .geneset.ls = NULL,
    .id.idx = NULL,
    .id = NULL,
    .id.from = NULL,
    .verbose = TRUE,
    .label.confidence = 0.8
) {
  .Deprecated("get.pbDE", 
              msg = "get.dea() is deprecated. Use get.pbDE() instead, which stores results in .tdr.obj$pbDE.")
  
  # Call get.pbDE and extract results for backward compatibility
  result <- get.pbDE(
    .tdr.obj = .tdr.obj,
    .design = .design,
    .contrasts = .contrasts,
    .block = .block,
    .geneset.ls = .geneset.ls,
    .id.idx = .id.idx,
    .id = .id,
    .id.from = .id.from,
    .model.name = ".deprecated.call",
    .population.name = "result",
    .force.recalc = TRUE,
    .verbose = .verbose,
    .label.confidence = .label.confidence
  )
  
  # Return just the DE results (legacy behavior)
  return(result$pbDE$.deprecated.call$result)
}

#' Marker Gene/Protein Identification via Pairwise Comparison
#'
#' Identifies marker genes/proteins that distinguish one cell subset from another (or all others) by 
#' performing pseudobulk differential expression between groups. Unlike \code{get.pbDE} which tests 
#' experimental conditions, this compares cell populations to find defining features. Useful for 
#' characterizing clusters and validating cell type annotations.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .geneset.ls Optional named list of character vectors defining gene sets for GSVA enrichment. 
#'   Only for RNA data.
#' @param .id1.idx Optional integer vector specifying landmark indices for group 1. When provided,
#'   the \code{.label.confidence} pipeline determines which cells confidently belong to this 
#'   landmark set based on their fuzzy graph weights.
#' @param .id2.idx Optional integer vector specifying landmark indices for group 2. When provided,
#'   the \code{.label.confidence} pipeline determines which cells confidently belong to this 
#'   landmark set based on their fuzzy graph weights.
#' @param .id1 Character vector of cluster/celltype IDs for group 1 (test group). Required if 
#'   \code{.id1.idx} not provided.
#' @param .id2 Character vector of cluster/celltype IDs for group 2 (reference group). Default 
#'   \code{"..all.other.landmarks.."} compares group 1 to all other cells. Can specify specific IDs 
#'   for pairwise comparisons.
#' @param .id.from Character: "clustering" (default) or "celltyping". Source of IDs in \code{.id1} 
#'   and \code{.id2}.
#' @param .model.name Character: A name for the analysis model (default: "default"). Used for 
#'   organizing multiple analyses.
#' @param .comparison.name Character: A name for this comparison. If NULL (default), auto-generated 
#'   from .id1 and .id2 values (e.g., "cluster.1_vs_all" or "cluster.1_vs_cluster.2").
#' @param .force.recalc Logical: If TRUE, recalculate even if results exist (default: FALSE).
#' @param .label.confidence Numeric (0.5-1): minimum confidence threshold for cell-to-population
#'   assignment when using \code{.id1.idx} or \code{.id2.idx} (default: 0.8). For each cell, the
#'   fraction of its fuzzy graph weight falling on "in" vs "out" landmarks is computed; cells below
#'   this threshold are excluded as low-confidence assignments.
#'   
#' @return The input \code{.tdr.obj} with results stored in 
#'   \code{.tdr.obj$markerDE[[.model.name]][[.comparison.name]]} as a limma MArrayLM fit object 
#'   with moderated statistics and additional slots:
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
#' \strong{Difference from get.pbDE:}
#' \itemize{
#'   \item \code{get.pbDE}: Tests experimental conditions (treatment vs control) within populations
#'   \item \code{get.markerDE}: Tests cell populations (cluster A vs B) within the same experiment
#' }
#' 
#' @seealso \code{\link{get.map}} (required), \code{\link{get.pbDE}} for condition-based DE,
#'   \code{\link{plotHeatmap}} for visualizing expression patterns, \code{\link{plotMarkerDE}}
#'   for visualizing marker results
#'   
#' @examples
#' \dontrun{
#' # After mapping and clustering
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#' 
#' # Find markers for cluster 1 vs all others (stored in .tdr.obj$markerDE)
#' lm.cells <- get.markerDE(lm.cells,
#'                          .id1 = "cluster.1",
#'                          .id.from = "clustering",
#'                          .comparison.name = "cluster1_vs_all")
#' 
#' # Access results
#' marker.results <- lm.cells$markerDE$default$cluster1_vs_all
#' 
#' # Top upregulated markers (accessing .id1 coefficient)
#' top.markers <- rownames(marker.results$coefficients)[
#'   marker.results$adj.p[, ".id1"] < 0.05 & 
#'   marker.results$coefficients[, ".id1"] > 1
#' ]
#' 
#' # Pairwise comparison: CD4 T cells vs CD8 T cells
#' lm.cells <- get.markerDE(lm.cells,
#'                          .id1 = c("cluster.1", "cluster.3"),  # CD4 T clusters
#'                          .id2 = c("cluster.2", "cluster.4"),  # CD8 T clusters
#'                          .id.from = "clustering",
#'                          .comparison.name = "cd4_vs_cd8")
#' 
#' # With gene set enrichment
#' hallmark.sets <- list(
#'   "GLYCOLYSIS" = c("HK2", "PFKP", "LDHA"),
#'   "OXIDATIVE_PHOS" = c("ATP5A", "COX5A", "NDUFA4")
#' )
#' lm.cells <- get.markerDE(lm.cells,
#'                          .id1 = "cluster.1",
#'                          .geneset.ls = hallmark.sets,
#'                          .comparison.name = "cluster1_gsva")
#' }
#' 
#' @export
#'
get.markerDE <-
  function(
    .tdr.obj,
    .geneset.ls = NULL,
    .id1.idx = NULL,
    .id2.idx = NULL,
    .id1 = NULL,
    .id2 = "..all.other.landmarks..",
    .id.from = "clustering",
    .model.name = "default",
    .comparison.name = NULL,
    .force.recalc = FALSE,
    .label.confidence = 0.8
  ) {
    
    i <- j <- landmark <- cell <- label <- x <- confidence <-
      NULL
    
    if(.label.confidence < 0.5 || .label.confidence > 1){
      stop(".label.confidence must be between 0.5 and 1")
    }
    
    if(!is.null(x = .geneset.ls)){
      if(.tdr.obj$config$assay.type != "RNA"){
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
      if(!all(.id1 %in% unique(x = .tdr.obj$graph[[.id.from]]$ids))){
        
        stop(paste0(paste0(.id1[!(.id1 %in% unique(x = .tdr.obj$graph[[.id.from]]$ids))],
                           collapse = ", "),
                    " not found in ",
                    .id.from))
        
      }
    }
    
    if(is.null(x = .id2.idx)){
      if(!all(.id2 %in% c("..all.other.landmarks..",
                          unique(x = .tdr.obj$graph[[.id.from]]$ids) |>
                          as.character()))){
        
        stop(paste0(paste0(.id2[!(.id2 %in% c("..all.other.landmarks..",
                                              unique(x = .tdr.obj$graph[[.id.from]]$ids) |>
                                                as.character()))],
                           collapse = ", "),
                    " not found in ",
                    .id.from))
        
      }
    }
    
    if(is.null(x = .id1.idx)){
      
      # Cell-centric: get per-sample cell indices matching group 1
      .lm1.idx <-
        lapply(X = .tdr.obj$map[[.id.from]]$ids,
               FUN = function(smpl){
                 which(x = smpl %in% .id1)
               })
      
    } else {
      
      if(!all(.id1.idx %in% (nrow(x = .tdr.obj$landmarks) |> seq_len()))) {
        stop(paste0(".id1.idx must be an integer vector between 1 and ",
                    nrow(x = .tdr.obj$landmarks)))
      }
      
      # Cell-centric: use .label.confidence to determine which cells
      # confidently belong to the specified landmark set (group 1)
      tmp.lbl1 <-
        rep(x = "out",
            times = nrow(x = .tdr.obj$landmarks))
      
      tmp.lbl1[.id1.idx] <-
        "in"
      
      .lm1.idx <-
        lapply(X = .tdr_get_map_slot_all(.tdr.obj, "fuzzy.graph"),
               FUN = function(smpl) {
                 
                 in.and.out <-
                   Matrix::summary(object = smpl) |>
                   dplyr::rename(cell = i,
                                 landmark = j) |>
                   dplyr::mutate(label = tmp.lbl1[landmark]) |>
                   dplyr::select(-landmark) |>
                   collapse::fgroup_by(cell,
                                       label,
                                       sort = FALSE) |>
                   collapse::fsum() |>
                   collapse::fgroup_by(cell,
                                       sort = FALSE) |>
                   collapse::fmutate(confidence = x / collapse::fsum(x)) |>
                   collapse::fungroup() |>
                   collapse::fsubset(confidence >= .label.confidence) |>
                   (\(x)
                    {
                      label <-
                        rep(x = "..low.confidence..",
                            times = nrow(x = smpl))
                      
                      label[x$cell] <-
                        x$label
                      
                      label
                      
                   }
                   )()
                 
                 which(x = in.and.out == "in")
                 
               })
      
    }
    
    if("..all.other.landmarks.." %in% .id2){
      
      message("using `..all.other.landmarks..` for .id2")
      
      .id2 <- "..all.other.landmarks.."
      
      # Cell-centric complement: all cells NOT in group 1
      .lm2.idx <-
        seq_along(along.with = .lm1.idx) |>
        stats::setNames(nm = names(x = .lm1.idx)) |>
        lapply(FUN = function(smpl){
          seq_len(length.out = .tdr.obj$config$sampling$n.cells[smpl])[-.lm1.idx[[smpl]]]
        })
      
    } else {
      
      if(is.null(x = .id2.idx)){
        
        # Cell-centric: get per-sample cell indices matching group 2
        .lm2.idx <-
          lapply(X = .tdr.obj$map[[.id.from]]$ids,
                 FUN = function(smpl){
                   which(x = smpl %in% .id2)
                 })
        
      } else {
        
        if(!all(.id2.idx %in% (nrow(x = .tdr.obj$landmarks) |> seq_len()))) {
          stop(paste0(".id2.idx must be an integer vector between 1 and ",
                      nrow(x = .tdr.obj$landmarks)))
        }
        
        # Cell-centric: use .label.confidence to determine which cells
        # confidently belong to the specified landmark set (group 2)
        tmp.lbl2 <-
          rep(x = "out",
              times = nrow(x = .tdr.obj$landmarks))
        
        tmp.lbl2[.id2.idx] <-
          "in"
        
        .lm2.idx <-
          lapply(X = .tdr_get_map_slot_all(.tdr.obj, "fuzzy.graph"),
                 FUN = function(smpl) {
                   
                   in.and.out <-
                     Matrix::summary(object = smpl) |>
                     dplyr::rename(cell = i,
                                   landmark = j) |>
                     dplyr::mutate(label = tmp.lbl2[landmark]) |>
                     dplyr::select(-landmark) |>
                     collapse::fgroup_by(cell,
                                         label,
                                         sort = FALSE) |>
                     collapse::fsum() |>
                     collapse::fgroup_by(cell,
                                         sort = FALSE) |>
                     collapse::fmutate(confidence = x / collapse::fsum(x)) |>
                     collapse::fungroup() |>
                     collapse::fsubset(confidence >= .label.confidence) |>
                     (\(x)
                      {
                        label <-
                          rep(x = "..low.confidence..",
                              times = nrow(x = smpl))
                        
                        label[x$cell] <-
                          x$label
                        
                        label
                        
                     }
                     )()
                   
                   which(x = in.and.out == "in")
                   
                 })
        
      }
      
    }
    
    # Auto-generate comparison name if not provided
    if(is.null(x = .comparison.name)){
      id1.str <- paste0(.id1, collapse = "_")
      if(identical(.id2, "..all.other.landmarks..")){
        id2.str <- "all"
      } else {
        id2.str <- paste0(.id2, collapse = "_")
      }
      .comparison.name <- paste0(id1.str, "_vs_", id2.str)
    }
    
    # Initialize storage if needed
    if(is.null(x = .tdr.obj$markerDE)){
      .tdr.obj$markerDE <- list()
    }
    if(is.null(x = .tdr.obj$markerDE[[.model.name]])){
      .tdr.obj$markerDE[[.model.name]] <- list()
    }
    
    # Check for existing results
    if(!.force.recalc && !is.null(x = .tdr.obj$markerDE[[.model.name]][[.comparison.name]])){
      message(sprintf("Results already exist at .tdr.obj$markerDE$%s$%s. Use .force.recalc = TRUE to recalculate.", 
                      .model.name, .comparison.name))
      return(.tdr.obj)
    }
    
    # number of cells in each .id1 and .id2 pseudobulks
    n.pseudo1 <-
      lengths(x = .lm1.idx)
    
    n.pseudo2 <-
      lengths(x = .lm2.idx)
    
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
    if(.tdr.obj$config$assay.type == "RNA"){
      
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
        stats::setNames(object = names(x = .tdr.obj$cells)[!smpl.outlier.1],
                        nm = names(x = .tdr.obj$cells)[!smpl.outlier.1]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .tdr.obj$cells[[smpl]])
          
          res <-
            tryCatch(
              expr = {
                
                cell.idx <- .lm1.idx[[smpl]]
                
                wcl <- 
                  .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", smpl)[cell.idx,,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  exprs.mat[,cell.idx] %*% wcl
                
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
        stats::setNames(object = names(x = .tdr.obj$cells)[!smpl.outlier.2],
                        nm = names(x = .tdr.obj$cells)[!smpl.outlier.2]) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .tdr.obj$cells[[smpl]])
          
          res <-
            tryCatch(
              expr = {
                
                cell.idx <- .lm2.idx[[smpl]]
                
                wcl <- 
                  .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", smpl)[cell.idx,,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  exprs.mat[,cell.idx] %*% wcl
                
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
        stats::setNames(object = names(x = .tdr.obj$cells),
                        nm = names(x = .tdr.obj$cells)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .tdr.obj$cells[[smpl]])
          
          cols <-
            colnames(x = .tdr.obj$landmarks)
          
          res <-
            tryCatch(
              expr = {
                
                cell.idx <- .lm1.idx[[smpl]]
                
                wcl <- 
                  .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", smpl)[cell.idx,,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  (Matrix::t(x = exprs.mat[cell.idx, cols]) %*% wcl)
                
                # get weighted mean using sum of weights
                tmp.res <-
                  Matrix::rowSums(x = wsum) / sum(wcl)
                
                return(tmp.res)
                
              },
              error = function(e) {
                
                stats::setNames(object = rep(x = 0,
                                             times = length(x = cols)),
                                nm = cols)
                
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
        stats::setNames(object = names(x = .tdr.obj$cells),
                        nm = names(x = .tdr.obj$cells)) |>
        lapply(FUN = function(smpl){
          
          exprs.mat <-
            readRDS(file = .tdr.obj$cells[[smpl]])
          
          cols <-
            colnames(x = .tdr.obj$landmarks)
          
          res <-
            tryCatch(
              expr = {
                
                cell.idx <- .lm2.idx[[smpl]]
                
                wcl <- 
                  .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", smpl)[cell.idx,,drop = FALSE]
                
                # get weighted sum
                wsum <-
                  (Matrix::t(x = exprs.mat[cell.idx, cols]) %*% wcl)
                
                # get weighted mean using sum of weights
                tmp.res <-
                  Matrix::rowSums(x = wsum) / sum(wcl)
                
                return(tmp.res)
                
              },
              error = function(e) {
                
                stats::setNames(object = rep(x = 0,
                                             times = length(x = cols)),
                                nm = cols)
                
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
      .lm1.idx
    
    .de$id2.idx <-
      .lm2.idx
    
    .de$n.pseudo1 <-
      n.pseudo1
    
    .de$n.pseudo2 <-
      n.pseudo2
    
    # Store results in .tdr.obj
    .tdr.obj$markerDE[[.model.name]][[.comparison.name]] <- .de
    
    message(sprintf("Marker DE results stored at .tdr.obj$markerDE$%s$%s", .model.name, .comparison.name))
    
    return(.tdr.obj)
    
  }

#' Deprecated: Use get.markerDE instead
#'
#' \code{get.marker()} has been renamed to \code{\link{get.markerDE}()} for API consistency. 
#' This function is deprecated and will be removed in a future release.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .geneset.ls Optional named list of character vectors for GSVA.
#' @param .id1.idx Optional landmark indices for group 1.
#' @param .id2.idx Optional landmark indices for group 2.
#' @param .id1 Cluster/celltype IDs for group 1.
#' @param .id2 Reference group IDs. Default \code{"..all.other.landmarks.."}.
#' @param .id.from "clustering" or "celltyping".
#' @param .label.confidence Numeric (0.5-1): minimum confidence for cell assignment (default: 0.8).
#'
#' @return A .tdr.obj with results stored in .tdr.obj$markerDE$.deprecated.call$result.
#' @seealso \code{\link{get.markerDE}}
#' @export
get.marker <- function(
    .tdr.obj,
    .geneset.ls = NULL,
    .id1.idx = NULL,
    .id2.idx = NULL,
    .id1 = NULL,
    .id2 = "..all.other.landmarks..",
    .id.from = "clustering",
    .label.confidence = 0.8
) {
  .Deprecated("get.markerDE", 
              msg = "get.marker() is deprecated. Use get.markerDE() instead, which stores results in .tdr.obj$markerDE.")
  
  # Call get.markerDE and return the .tdr.obj directly
  return(get.markerDE(
    .tdr.obj = .tdr.obj,
    .geneset.ls = .geneset.ls,
    .id1.idx = .id1.idx,
    .id2.idx = .id2.idx,
    .id1 = .id1,
    .id2 = .id2,
    .id.from = .id.from,
    .model.name = ".deprecated.call",
    .comparison.name = "result",
    .force.recalc = TRUE,
    .label.confidence = .label.confidence
  ))
}

# =============================================================================
# Helper function: Compute unsupervised embeddings (pca and traj)
# =============================================================================

#' @noRd
#' @keywords internal
.compute.unsupervised.embeddings <-
  function(
    .tdr.obj,
    .n.eigs = 20,
    .n.pcs = 20,
    .verbose = TRUE,
    .ret.trajectory = FALSE,
    .traj.dist.metric = "cosine"
  ){
    
    # Initialize embedding slot if not present
    if(is.null(x = .tdr.obj$map$embedding)){
      .tdr.obj$map$embedding <- list()
    }
    
    # -------------------------------------------------------------------------
    # Compute PCA if not already present
    # -------------------------------------------------------------------------
    
    if(is.null(x = .tdr.obj$map$embedding$pca)){
      
      if(isTRUE(x = .verbose)){
        message("\nComputing sample-level PCA on landmark densities")
      }
      
      # Compute PCA via irlba (samples as rows, using .tdr.obj$map$Y)
      pca.res <-
        Matrix::t(x = .tdr.obj$map$Y) |>
        (\(x)
         irlba::prcomp_irlba(
           x = x,
           center = Matrix::colMeans(x = x),
           scale. = FALSE,
           n = min(c(.n.pcs, nrow(x = x) - 1, ncol(x = x)))
         )
        )()
      
      # Store in format consistent with get.landmarks():
      # coord = embedding coordinates (samples x PCs), 
      # rotation = loadings (landmarks x PCs),
      # center, scale, sdev
      .tdr.obj$map$embedding$pca <- list(
        coord = pca.res$x,
        rotation = pca.res$rotation,
        center = pca.res$center,
        scale = pca.res$scale,
        sdev = pca.res$sdev
      )
      
      # Set column names for coord
      colnames(x = .tdr.obj$map$embedding$pca$coord) <- 
        paste0("PC", seq_len(length.out = ncol(x = .tdr.obj$map$embedding$pca$coord)))
      
      if(isTRUE(x = .verbose)){
        var.explained <- round(100 * pca.res$sdev^2 / sum(pca.res$sdev^2), 1)
        message("  PCA variance explained: ", 
                paste(paste0("PC", seq_along(var.explained), "=", var.explained, "%"), collapse = ", "))
      }
      
    } else {
      
      if(isTRUE(x = .verbose)){
        message("\nPCA embedding already exists; skipping recomputation")
      }
      
    }
    
    # -------------------------------------------------------------------------
    # Compute diffusion map trajectory (traj) if not already present
    # -------------------------------------------------------------------------
    
    if(isTRUE(x = .ret.trajectory)){
      
      if(is.null(x = .tdr.obj$map$embedding$traj)){
        
        if(isTRUE(x = .verbose)){
          message("\nComputing diffusion map trajectory embedding")
        }
        
        # Check if destiny is available
        if(!requireNamespace("destiny", quietly = TRUE)){
          warning("Package 'destiny' is required for trajectory embedding. ",
                  "Install with BiocManager::install('destiny'). ",
                  "Skipping trajectory computation.",
                  call. = FALSE)
          return(.tdr.obj)
        }
        
        # Compute diffusion map on transposed Y (samples as rows)
        Matrix.warnDeprecatedCoerce <-
          # destiny relies on Deprecated Matrix package function: see https://github.com/theislab/destiny/issues/61
          getOption(x = "Matrix.warnDeprecatedCoerce")
        
        options(Matrix.warnDeprecatedCoerce = 0)
        
        traj.res <- 
          Matrix::t(x = .tdr.obj$map$Y) |>
          (\(x)
           destiny::DiffusionMap(
             data = x,
             n_eigs = min(c(.n.eigs, nrow(x = x) - 2, ncol(x = x))), 
             n_pcs = NA, # suppress PCA inside destiny
             distance = .traj.dist.metric
           )
          )()
        
        options(Matrix.warnDeprecatedCoerce = Matrix.warnDeprecatedCoerce)
        
        # Extract eigenvalues and compute cell embeddings
        stdev <- traj.res@eigenvalues
        stdev[stdev < 0] <- 0
        
        D <- 
          abs(x = stdev) |>
          sqrt() |>
          diag()
        
        traj.coord <- 
          traj.res@eigenvectors %*% D
        
        # Set column names
        colnames(x = traj.coord) <- 
          paste0("DC", seq_len(length.out = ncol(x = traj.coord)))
        
        rownames(x = traj.coord) <- 
          rownames(x = .tdr.obj$metadata)
        
        # Store in format similar to pca slot
        .tdr.obj$map$embedding$traj <- list(
          coord = traj.coord,
          eigenvalues = traj.res@eigenvalues,
          eigenvectors = traj.res@eigenvectors,
          sdev = abs(x = stdev) |> sqrt()
        )
        
        if(isTRUE(x = .verbose)){
          message("  Diffusion map computed with ", ncol(x = traj.coord), " components")
        }
        
      } else {
        
        if(isTRUE(x = .verbose)){
          message("\nTrajectory embedding already exists; skipping recomputation")
        }
        
      }
    }
    
    return(.tdr.obj)
    
  }

#' Compute Sample Embedding from Partial Fitted Values
#'
#' Computes embeddings for sample-level visualization. When called without supervised
#' arguments, computes unsupervised embeddings (PCA and diffusion map trajectory) on
#' landmark densities stored in \code{.tdr.obj$map$Y}. When called with supervised arguments 
#' (\code{.contrast.of.interest} or \code{.red.model}), additionally computes a 
#' partial-effect PCA (pePC) that isolates variation attributable to a specific effect. 
#' Exact for OLS; if duplicateCorrelation/blocking is used, the decomposition is approximate.
#'
#' @param .tdr.obj The tinydenseR object processed through \code{get.map()}. Contains
#'   \code{$map$Y} (log2-transformed densities) used for unsupervised embeddings. Statistical
#'   model fits should be stored in \code{$map$lm} via \code{get.lm()}.
#' @param .full.model Character string naming the full model in \code{.tdr.obj$map$lm}
#'   (default "default"). Required only when computing supervised embeddings (pePC).
#' @param .term.of.interest Character string: the covariate/term being isolated when using the 
#'   nested model method (\code{.red.model}). Must match a column name in 
#'   \code{.tdr.obj$metadata}. Ignored when using \code{.contrast.of.interest}.
#' @param .red.model Optional character string naming the reduced model in \code{.tdr.obj$map$lm}
#'   (without the effect of interest). If provided, computes 
#'   \eqn{\Delta\hat{Y} = \hat{Y}_{full} - \hat{Y}_{red}} directly. Make sure to construct
#'   nested models by "dropping terms" (so reduced is a strict subset of full)!
#' @param .contrast.of.interest Optional: Character string naming the contrast to extract.
#'   Must match a column name in the full model's contrasts. When specified, uses
#'   the Frisch-Waugh-Lovell (FWL) theorem to compute the true partial fitted component.
#' @param .n.eigs Integer: number of eigenvectors for diffusion map (default 20).
#' @param .n.pcs Integer: number of PCs for unsupervised PCA and diffusion map (default 20).
#' @param .ret.trajectory Logical: whether to compute diffusion map trajectory embedding (default FALSE).
#' @param .traj.dist.metric Character: distance metric for diffusion map (default "cosine").
#' @param .verbose Logical: print progress messages? Default TRUE.
#'
#' @return Modified \code{.tdr.obj} with embeddings stored in \code{.tdr.obj$map$embedding}:
#'   \describe{
#'     \item{\code{pca}}{Unsupervised PCA on log-transformed landmark densities:
#'       \itemize{
#'         \item \code{coord}: Sample coordinates (samples x PCs)
#'         \item \code{rotation}: Loadings (landmarks x PCs)
#'         \item \code{center}, \code{scale}, \code{sdev}: Centering, scaling, and standard deviations
#'       }}
#'     \item{\code{traj}}{Diffusion map trajectory embedding:
#'       \itemize{
#'         \item \code{coord}: Sample coordinates (samples x DCs)
#'         \item \code{eigenvalues}, \code{eigenvectors}: Raw diffusion map components
#'         \item \code{sdev}: sqrt(abs(eigenvalues))
#'       }}
#'     \item{\code{pePC}}{Supervised partial-effect embeddings (when pePC args provided):
#'       \itemize{
#'         \item \code{<contrast>} or \code{<term>}: Named list containing:
#'           \itemize{
#'             \item \code{coord}: Nuisance-residualized projection (samples x pePCs)
#'             \item \code{x}: PCA scores on partial fitted values
#'             \item \code{rotation}: Landmark loadings from partial-effect PCA
#'             \item \code{effect.resid.cor}: Per-PC correlation diagnostic
#'             \item \code{delta.Yhat}: Partial fitted values matrix
#'             \item \code{method}: "fwl_contrast" or "nested_models"
#'           }
#'       }}
#'   }
#'
#' @details
#' \strong{Unsupervised Embeddings (always computed):}
#
# \code{pca}: Standard PCA on log2-transformed landmark densities (log2(fdens + 0.5)),
# centered and scaled.
#
# \code{traj}: Diffusion map trajectory computed via \code{destiny::DiffusionMap} on
# the log2-transformed landmark density matrix. Useful for visualizing continuous trajectories.
#
# \strong{Supervised Embedding (pePC, when .contrast.of.interest or .red.model provided):}
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
#' \strong{Method 2: Nested model comparison} (when \code{.red.model} is provided)
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
#' @seealso \code{\link{get.lm}} for model fitting, \code{\link{plotSampleEmbedding}} for
#'   visualizing the embedding
#'
#' @examples
#' \dontrun{
#' # ===========================================================================
#' # Example 0: Unsupervised embeddings only (before get.lm)
#' # ===========================================================================
#' 
#' # After get.map, compute unsupervised embeddings
#' lm.obj <- get.embedding(.tdr.obj = lm.obj)
#' 
#' # Access unsupervised embeddings
#' lm.obj$map$embedding$pca$coord   # sample PCA coordinates
#' lm.obj$map$embedding$traj$coord  # sample diffusion map coordinates
#' 
#' # Plot unsupervised embeddings
#' plotSampleEmbedding(lm.obj, .embedding = "pca", .color.by = "Condition")
#' plotSampleEmbedding(lm.obj, .embedding = "traj", .color.by = "Timepoint")
#' 
#' # ===========================================================================
#' # Example 1: Nested model comparison (no contrasts)
#' # ===========================================================================
#' 
#' # Fit full model with disease effect (stored as "full")
#' .design <- model.matrix(object = ~ disease + sex + age,
#'                         data = lm.obj$metadata)
#' lm.obj <- get.lm(.tdr.obj = lm.obj, .design = .design, 
#'                     .model.name = "full")
#' 
#' # Fit reduced model without disease (stored as "reduced")
#' .red.design <- model.matrix(object = ~ sex + age,
#'                             data = lm.obj$metadata)
#' lm.obj <- get.lm(.tdr.obj = lm.obj, .design = .red.design, 
#'                     .model.name = "reduced")
#' 
#' # Extract disease effect embedding (all stored in lm.obj)
#' lm.obj <- get.embedding(
#'     .tdr.obj = lm.obj,
#'     .full.model = "full",
#'     .term.of.interest = "disease",
#'     .red.model = "reduced"
#' )
#' 
#' # Plot the embedding (colored by disease)
#' plotSampleEmbedding(lm.obj, .embedding = "pePC", .sup.embed.slot = "disease")
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
#' lm.obj <- get.lm(lm.obj, .design = design, .contrasts = contrasts)
#'
#' # Extract each contrast embedding
#' lm.obj <- get.embedding(lm.obj, .contrast.of.interest = "TrtVsCtrl")
#' lm.obj <- get.embedding(lm.obj, .contrast.of.interest = "KO_TrtVsCtrl")
#'
#' # Plot the embedding (specify .color.by for metadata column)
#' plotSampleEmbedding(lm.obj, .embedding = "pePC",
#'                     .sup.embed.slot = "TrtVsCtrl", .color.by = "Group")
#'
#' # Access the embedding components (all in lm.obj$map$embedding)
#' lm.obj$map$embedding$pca$coord              # sample PCA coordinates
#' lm.obj$map$embedding$traj$coord             # sample trajectory coordinates
#' lm.obj$map$embedding$pePC$TrtVsCtrl$coord   # sample pePC coordinates
#' lm.obj$map$embedding$pePC$TrtVsCtrl$rotation # landmark loadings
#' lm.obj$map$embedding$pePC$TrtVsCtrl$delta.Yhat # partial fitted values
#'
#' # Multiple embeddings stored together
#' names(lm.obj$map$embedding)       # c("pca", "traj", "pePC")
#' names(lm.obj$map$embedding$pePC)  # c("TrtVsCtrl", "KO_TrtVsCtrl")
#' }
#'
#' @export
#'
get.embedding <-
  function(
    .tdr.obj,
    .full.model = "default",
    .term.of.interest = NULL,
    .red.model = NULL,
    .contrast.of.interest = NULL,
    .n.eigs = 20,
    .n.pcs = 20,
    .ret.trajectory = FALSE,
    .traj.dist.metric = "cosine",
    .verbose = TRUE
  ){
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    # Check that .tdr.obj has map$Y
    if(is.null(x = .tdr.obj$map$Y)){
      stop("'.tdr.obj$map$Y' not found. Run get.map() first.")
    }
    
    # Determine if this is unsupervised-only mode (no pePC args provided)
    unsupervised.only <- 
      is.null(x = .red.model) && is.null(x = .contrast.of.interest)
    
    # If supervised args provided, require the full model to exist
    if(!unsupervised.only){
      
      if(is.null(x = .tdr.obj$map$lm[[.full.model]])){
        stop(paste0("Model '", .full.model, "' not found in .tdr.obj$map$lm. ",
                    "Run get.lm() first, or omit '.contrast.of.interest'/'.red.model' for unsupervised-only."))
      }
      
    }
    
    if(!unsupervised.only &&
       !is.null(x = .red.model) &&
       !is.null(x = .contrast.of.interest)){
      
      warning("Both '.red.model' and '.contrast.of.interest' provided. Using '.contrast.of.interest' (FWL method).")
      
    }
    
    # Reference the stats object for convenience
    .stats.obj <- .tdr.obj$map$lm[[.full.model]]
    
    # -------------------------------------------------------------------------
    # Compute unsupervised embeddings (pca, traj) - stored in .tdr.obj$map$embedding
    # -------------------------------------------------------------------------
    
    .tdr.obj <- 
      .compute.unsupervised.embeddings(
        .tdr.obj = .tdr.obj,
        .n.eigs = .n.eigs,
        .n.pcs = .n.pcs,
        .ret.trajectory = .ret.trajectory,
        .traj.dist.metric = .traj.dist.metric,
        .verbose = .verbose
      )
    
    # If no supervised args provided, return .tdr.obj with unsupervised embeddings
    if(unsupervised.only){
      
      if(isTRUE(x = .verbose)){
        message("\nNo supervised embedding args provided. Returning .tdr.obj with unsupervised embeddings.")
      }
      
      return(.tdr.obj)
      
    }
    
    # -------------------------------------------------------------------------
    # Determine slot name for supervised embedding
    # -------------------------------------------------------------------------
    
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
      if(!(.term.of.interest %in% colnames(x = .tdr.obj$metadata))){
        
        stop("'.term.of.interest' must be a column name of '.tdr.obj$metadata'")
        
      }
      
      slot.name <- .term.of.interest
      
    }
    
    # -------------------------------------------------------------------------
    # METHOD 1: FWL-based contrast extraction
    # -------------------------------------------------------------------------
    
    if(!is.null(x = .contrast.of.interest)){
      
      if(is.null(x = .stats.obj$fit$contrasts)){
        stop("No contrasts found in '.stats.obj'. Run get.lm() with .contrasts argument.")
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
      
      # Get the expression matrix Y (landmarks x samples) from .tdr.obj$map$Y
      Y <- 
        .tdr.obj$map$Y
      
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
            lapply(X = .tdr.obj$metadata,
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
                                          block = .tdr.obj$metadata[[.block]])
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
                       block = if(!is.null(dupcor)) .tdr.obj$metadata[[.block]] else NULL,
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
      
      # Store nuisance design Z for downstream use (get.lm.spectral.dea)
      # If no nuisance, Z is intercept-only
      if(!exists("Z") || is.null(Z)){
        Z <- matrix(1, nrow = ncol(x = Y), ncol = 1, 
                    dimnames = list(colnames(Y), "(Intercept)"))
      }
      nuisance.design <- Z
      
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
        stop("Contrasts found in full model '", .full.model, "'. That means, coefficients are 're-oriented' via contrasts.fit. Use '.contrast.of.interest' instead.")
      }
      
      # Validate reduced model exists
      if(is.null(x = .tdr.obj$map$lm[[.red.model]])){
        stop(paste0("Reduced model '", .red.model, "' not found in .tdr.obj$map$lm. ",
                    "Run get.lm() with .model.name = '", .red.model, "' first."))
      }
      
      .red.stats.obj <- .tdr.obj$map$lm[[.red.model]]
      
      if(isTRUE(x = .verbose)){
        message("\nComputing embedding via nested model comparison")
        message("  Slot name: ", slot.name)
        message("  Full model ('", .full.model, "') columns: ", paste(colnames(.stats.obj$fit$design), collapse = ", "))
        message("  Reduced model ('", .red.model, "') columns: ", paste(colnames(.red.stats.obj$fit$design), collapse = ", "))
      }
      
      # Get Y from .tdr.obj$map$Y for dimension validation
      Y <- .tdr.obj$map$Y
      
      # Validate dimensions against model fits
      if(nrow(x = .stats.obj$fit$coefficients) != nrow(x = Y)){
        stop("Number of landmarks in full model does not match .tdr.obj$map$Y")
      }
      if(nrow(x = .red.stats.obj$fit$coefficients) != nrow(x = Y)){
        stop("Number of landmarks in reduced model does not match .tdr.obj$map$Y")
      }
      if(ncol(x = .stats.obj$fit$design) != ncol(x = .red.stats.obj$fit$design) ||
         nrow(x = .stats.obj$fit$design) != nrow(x = .red.stats.obj$fit$design)){
        # Different sample counts is ok if designs have different # of samples
        # But same samples should have same design row count
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
        rownames(x = Y)
      
      colnames(x = delta.Yhat) <- 
        colnames(x = Y)
      
      # Store nuisance design Z for downstream use (get.lm.spectral.dea)
      # For nested models, nuisance = reduced model design
      nuisance.design <- .red.stats.obj$fit$design
      
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
    
    # Compute rank analytically (O(1)) instead of expensive QR decomposition
    if(method == "fwl_contrast"){
      # Contrasts are rank-1 by construction (outer product gamma_g (x) x_c_perp)
      embed.rank <- 1L
    } else {
      # Nested models: rank = df_full - df_red, capped at matrix dimensions
      embed.rank <- 
        min(
          ncol(x = .stats.obj$fit$design) - ncol(x = .red.stats.obj$fit$design),
          nrow(x = delta.Yhat.t) - 1L,
          ncol(x = delta.Yhat.t)
        )
    }
    embed.rank <- 
      max(1L,
          embed.rank)
    
    
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
      Matrix::t(x = (.tdr.obj$map$Y - Yhat.red) - pca$center) %*%
      pca$rotation
    
    colnames(x = pca$coord) <- 
      paste0("pePC", 
             seq_len(length.out = ncol(x = pca$coord)))
    
    # partial-effect principal component variance as a fraction of TOTAL variance
    perc.tot.var.exp <-
      (100 * ((pca$sdev[1:ncol(x = pca$coord)])^2) /
         (matrixStats::rowVars(x = .tdr.obj$map$Y) |>
            sum())) |>
      stats::setNames(nm = colnames(x = pca$coord))
    
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
    # Store results in .tdr.obj$map$embedding$pePC[[slot.name]]
    # -------------------------------------------------------------------------
    
    # Initialize pePC list if needed
    if(is.null(x = .tdr.obj$map$embedding$pePC)){
      .tdr.obj$map$embedding$pePC <- list()
    }
    
    .tdr.obj$map$embedding$pePC[[slot.name]] <- 
      c(pca,
        list(perc.tot.var.exp = perc.tot.var.exp,
             effect.resid.cor = effect.resid.cor,
             delta.Yhat = as.matrix(x = delta.Yhat),
             nuisance.design = nuisance.design,
             method = method)
      )
    
    if(isTRUE(x = .verbose)){
      message("\nSupervised embedding stored in: .tdr.obj$map$embedding$pePC$", slot.name)
      message("Unsupervised embeddings stored in: .tdr.obj$map$embedding$pca, .tdr.obj$map$embedding$traj")
    }
    
    # Return .tdr.obj with all embeddings
    return(.tdr.obj)
    
  }


# =============================================================================
# Internal helpers shared by get.specDE, get.nmfDE, get.plsDE
# =============================================================================

# -----------------------------------------------------------------------------
# .validate.DE.inputs
# Validates common preconditions for all *DE methods.
# Returns coef.mat (the coefficient matrix from the fitted model).
# -----------------------------------------------------------------------------

.validate.DE.inputs <-
  function(.tdr.obj,
           .model.name,
           .coef.col) {
    
    # Validate on-disk cache (if active) before proceeding
    .tdr_cache_validate_quiet(.tdr.obj)
    
    if (is.null(x = .tdr.obj$map$lm[[.model.name]])) {
      stop("Model '", .model.name, "' not found. Run get.lm() first.")
    }
    
    coef.mat <-
      .tdr.obj$map$lm[[.model.name]]$fit$coefficients
    
    if (!(.coef.col %in% colnames(x = coef.mat))) {
      stop("Coefficient '", .coef.col, "' not found in model coefficients.\n",
           "Available: ", paste(colnames(x = coef.mat), collapse = ", "))
    }
    
    if (is.null(x = .tdr.obj$graph$snn)) {
      stop("SNN graph not found. Run get.graph() first.")
    }
    
    if (!Matrix::isSymmetric(object = .tdr.obj$graph$snn)) {
      stop("SNN graph not symmetric.")
    }
    
    if (is.null(x = .tdr.obj$raw.landmarks)) {
      stop("Raw landmarks not found. Run get.landmarks() first.")
    }
    
    if (is.null(x = .tdr.obj$pca$embed)) {
      stop("PCA embedding not found. Run get.landmarks() first.")
    }
    
    return(coef.mat)
    
  }


# -----------------------------------------------------------------------------
# .build.P
# Builds the random-walk normalized transition matrix from SNN, with optional
# degree regularization and lazy walk.  Returns P (sparse matrix).
# -----------------------------------------------------------------------------

.build.P <-
  function(SNN,
           .degree.reg = FALSE,
           .tau.mult = 1,
           .lazy.alpha = 1,
           .verbose = TRUE) {
    
    # Validate lazy.alpha
    if (!is.numeric(x = .lazy.alpha) || .lazy.alpha <= 0 || .lazy.alpha > 1) {
      stop(".lazy.alpha must be in (0, 1]")
    }
    
    # Check for disconnected nodes
    d <-
      Matrix::rowSums(x = SNN)
    
    if (any(d == 0)) {
      n.disconnected <- sum(d == 0)
      stop("Graph contains ", n.disconnected, " disconnected node(s) (degree = 0). ",
           "Cannot compute random-walk matrix. Check graph construction or increase k.")
    }
    
    # Degree regularization: add self-loops W_tau = W + tau*I, then row-normalize
    if (isTRUE(x = .degree.reg)) {
      tau <- .tau.mult * mean(x = d)
      W.tau <- SNN + tau * Matrix::Diagonal(n = nrow(x = SNN))
      d.tau <- Matrix::rowSums(x = W.tau)
      D.inv <- Matrix::Diagonal(x = 1 / d.tau)
      P <- D.inv %*% W.tau
      if (isTRUE(x = .verbose)) {
        message("  Degree regularization: tau = ", round(x = tau, digits = 2),
                " (tau.mult = ", .tau.mult, ")")
      }
    } else {
      D.inv <- Matrix::Diagonal(x = 1 / d)
      P <- D.inv %*% SNN
    }
    
    # Set dimnames
    P <-
      `dimnames<-`(x = P,
                   value = dimnames(x = SNN))
    
    # Lazy random walk: P_lazy = (1 - alpha) * I + alpha * P
    if (.lazy.alpha < 1) {
      I.mat <- Matrix::Diagonal(n = nrow(x = P))
      P <- (1 - .lazy.alpha) * I.mat + .lazy.alpha * P
      if (isTRUE(x = .verbose)) {
        message("  Lazy walk: alpha = ", .lazy.alpha)
      }
    }
    
    return(P)
    
  }


# -----------------------------------------------------------------------------
# .prepare.X
# Filters, normalizes (size-factor + log2 for RNA), and centers the expression
# matrix.  Returns a list with:
#   $X  — pre-centering (sparse; useful for sparse-aware loadings and nmfDE mass)
#   $Xc — column-centered (sparse; needed for M.local and DE-semantic loadings)
# -----------------------------------------------------------------------------

.prepare.X <-
  function(.tdr.obj,
           .min.prop = 0.005,
           .verbose = TRUE) {
    
    if (.tdr.obj$config$assay.type == "RNA") {
      
      # Filter genes: detected in at least min.prop of landmarks
      X <-
        .tdr.obj$raw.landmarks |>
        (\(x)
         x[, Matrix::colSums(x = x > 0) > (nrow(x = x) * .min.prop)]
        )()
      
      if (isTRUE(x = .verbose)) {
        message("  Genes after filtering (>", .min.prop * 100, "% detection): ",
                ncol(x = X), " / ", ncol(x = .tdr.obj$raw.landmarks))
      }
      
      # Size factor normalization
      X <-
        X |>
        Matrix::t() |>
        (\(x)
         x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))
        )() |>
        Matrix::t()
      
      # Log2 transform (in-place for sparse matrix efficiency)
      X@x <-
        log2(x = X@x + 1)
      
    } else {
      
      # Cytometry: use raw landmarks directly
      X <-
        .tdr.obj$raw.landmarks
      
    }
    
    # Center expression (genes/markers)
    Xc <-
      X |>
      (\(x)
       Matrix::t(x = x) - Matrix::colMeans(x = x)
      )() |>
      Matrix::t()
    
    return(list(X = X, Xc = Xc))
    
  }


# -----------------------------------------------------------------------------
# .compute.Sk
# Laplacian smoothness for each column of a score matrix.
# Uses the normalized Laplacian: L_sym = I - D^{-1/2} W D^{-1/2}.
# Smoothness_k = 1 - energy_k / 2, where energy_k is the Dirichlet energy
# of the centered score vector against L_sym.
# Returns: named numeric vector of smoothness values.
# -----------------------------------------------------------------------------

.compute.Sk <-
  function(scores,
           SNN) {
    
    D.vec <-
      Matrix::rowSums(x = SNN)
    
    D.inv.sqrt <-
      Matrix::Diagonal(x = 1 / sqrt(x = D.vec))
    
    L.sym <-
      Matrix::Diagonal(n = nrow(x = SNN)) -
      (D.inv.sqrt %*% SNN %*% D.inv.sqrt)
    
    smoothness <-
      seq_len(length.out = ncol(x = scores)) |>
      sapply(FUN = function(k) {
        
        s <-
          scores[, k] |>
          (\(x)
           x - mean(x = x)
          )()
        
        # Laplacian Dirichlet energy
        energy <-
          as.numeric(x = t(x = s) %*% L.sym %*% s) /
          sum(s^2)
        
        # Smoothness = 1 - energy (capped at 1 for normalized Laplacian)
        1 - pmin(pmax(energy, 0), 2) / 2
      })
    
    names(x = smoothness) <-
      colnames(x = scores)
    
    return(smoothness)
    
  }


# -----------------------------------------------------------------------------
# .compute.loadings
# Compute feature loadings for each score column against an expression matrix.
#   .method = "pearson"  — Pearson correlation (scale-invariant, sparse path)
#   .method = "ols"      — OLS regression beta (preserves magnitude, sparse path)
#   .method = "spearman" — Spearman rank correlation (robust, calls .sparse.spearman)
# X.sparse: pre-centering sparse expression matrix (exploits structural zeros).
# scores:   matrix (landmarks x K).
# Returns: matrix (genes x K).
# -----------------------------------------------------------------------------

.compute.loadings <-
  function(scores,
           X.sparse,
           .method = c("pearson", "ols", "spearman")) {
    
    .method <-
      match.arg(arg = .method,
                choices = c("pearson", "ols", "spearman"))
    
    n <- nrow(x = scores)
    n.genes <- ncol(x = X.sparse)
    comp.names <- colnames(x = scores)
    ncomp <- ncol(x = scores)
    
    if (.method == "pearson") {
      
      # Gene-level sufficient statistics (sparse-efficient: X^2 stays sparse)
      col.sums <-
        Matrix::colSums(x = X.sparse) |>
        as.numeric()
      
      col.sum.sq <-
        Matrix::colSums(x = X.sparse^2) |>
        as.numeric()
      
      gene.sd <-
        sqrt(x = pmax(col.sum.sq - (col.sums^2) / n, 0) /
               (n - 1))
      
      loadings <-
        seq_len(length.out = ncomp) |>
        stats::setNames(nm = comp.names) |>
        lapply(FUN = function(k) {
          
          score.k <-
            scores[, k] - mean(x = scores[, k])
          
          sd.score <-
            sqrt(x = sum(score.k^2) / (n - 1))
          
          if (sd.score < .Machine$double.eps) return(rep(x = NA_real_, times = n.genes))
          
          # crossprod on sparse X: touches only nonzero entries
          num <-
            (Matrix::crossprod(x = score.k,
                               y = X.sparse) |>
               as.numeric()) / (n - 1)
          
          r <-
            num / (sd.score * gene.sd)
          
          # Genes with zero variance: set to 0 (no information)
          r[gene.sd < .Machine$double.eps] <-
            0
          
          return(r)
          
        }) |>
        do.call(what = cbind)
      
    } else if (.method == "spearman") {
      
      loadings <-
        seq_len(length.out = ncomp) |>
        stats::setNames(nm = comp.names) |>
        lapply(FUN = function(k) {
          .sparse.spearman(s = scores[, k],
                           X.sparse = X.sparse)
        }) |>
        do.call(what = cbind)
      
    } else {
      
      # OLS regression: with centered score, crossprod on sparse X gives same
      # result as on centered Xc (centering term vanishes when sum(score.k) = 0)
      loadings <-
        seq_len(length.out = ncomp) |>
        stats::setNames(nm = comp.names) |>
        lapply(FUN = function(k) {
          
          score.k <-
            scores[, k] - mean(x = scores[, k])
          
          denom <-
            sum(score.k^2)
          
          if (denom < .Machine$double.eps) return(rep(x = NA_real_, times = n.genes))
          
          beta <-
            (Matrix::crossprod(x = score.k,
                               y = X.sparse) |>
               as.numeric()) /
            denom
          
          return(beta)
          
        }) |>
        do.call(what = cbind)
      
    }
    
    rownames(x = loadings) <-
      colnames(x = X.sparse)
    
    return(loadings)
    
  }


#' Spectral Differential Expression (specDE)
#'
#' Decomposes gene/marker expression into latent programs that covary with a
#' density contrast field across the cell graph. specDE identifies spectral modes
#' of a graph-conditioned differential expression operator, providing continuous
#' program scores per landmark and program-specific gene/marker loadings.
#'
#' @details
#' specDE answers: "What are the dominant expression programs that covary with a 
#' density contrast across the graph?" The method:
#' \enumerate{
#'   \item Extracts the density contrast vector Y (coefficients from \code{get.lm})
#'   \item Prepares expression matrix X (size-factor normalized + log2 for RNA; raw for cyto)
#'   \item Builds degree-regularized, lazy random-walk graph P from SNN
#'   \item Computes Y-weighted expression interaction: \code{diag(Y) \%*\% X}
#'   \item Graph-smooths: \code{M.local = P \%*\% diag(Y) \%*\% X}, then centers
#'   \item Decomposes M.local via SVD into latent programs
#'   \item Orders components by \code{abs(cor(Y, score))} (specDE1 = most Y-aligned)
#'   \item Computes loadings via regression of X on each component
#'   \item Computes Laplacian smoothness (Sk) for each component
#' }
#'
#' \strong{Diagnostic metrics (Ak, Vk, Sk):}
#'
#' \strong{Ak (Y-alignment):} Measures how strongly component scores covary with density 
#' contrast Y. High Ak indicates strong coupling to the tested contrast; low Ak indicates
#' weak coupling. Components are ordered by Ak (specDE1 = most Y-aligned). Note: high Ak
#' measures association with Y, not mechanism---it does not imply DA is the cause.
#'
#' \strong{Vk (variance explained):} Proportion of M.local variance captured by the
#' component. High Vk indicates a dominant mode; low Vk indicates a minor/specific pattern.
#' Vk is not a relevance score---high Vk does not imply biological importance or Y-coupling.
#'
#' \strong{Sk (graph smoothness):} Derived from the normalized Laplacian. High Sk indicates
#' large-scale, graph-smooth structure (low-frequency); low Sk indicates localized, rapidly
#' varying structure (high-frequency). Sk measures spatial scale, not biological mechanism---
#' both DA and DE programs can be smooth. A smooth (high Sk) component can still represent
#' DE-only biology if its association with Y is weak or non-monotonic.
#'
#' \strong{Joint interpretation (DA vs DE):}
#' \itemize{
#'   \item DA-dominated (typical): high Ak, high Sk, moderate-high Vk
#'   \item DE-dominated (functional): moderate/low Ak, any Sk, often lower Vk
#'   \item Structural/constraint: low Ak, high Vk, often high Sk
#' }
#' No single metric determines DA vs DE. Interpretation requires their joint pattern plus
#' visualization (e.g., Y vs score scatter via \code{plotSpecDE}).
#'
#' @param .tdr.obj A tinydenseR object after \code{get.lm()}.
#' @param .coef.col Character: column name in model coefficients to use as density contrast Y.
#'   Must be a valid column in \code{.tdr.obj$map$lm[[.model.name]]$fit$coefficients}.
#' @param .model.name Character: name of the fitted model to use (default "default").
#' @param .nv Integer: number of specDE components to compute. Defaults to 
#'   \code{ncol(.tdr.obj$pca$embed)}, matching the number of PCs from \code{get.landmarks}.
#' @param .min.prop Numeric: for RNA, minimum proportion of landmarks where a gene must be 
#'   detected (>0) to be included. Default 0.005 (0.5 percent).
#' @param .store.M Logical: if TRUE, store the M.local matrix in output. Default FALSE
#'   (saves memory; M.local can be large for RNA).
#' @param .degree.reg Logical: if TRUE, apply degree regularization by adding
#'   tau * I (self-loops) to W before row-normalizing: W_tau = W + tau*I, then
#'   P = D(W_tau)^(-1) W_tau. This stabilizes diffusion on graphs with heterogeneous
#'   degree (e.g., small k). Default FALSE (standard random walk).
#' @param .tau.mult Numeric: multiplier for tau when .degree.reg = TRUE.
#'   tau = .tau.mult * mean(degree). Default 1. Use smaller values (e.g., 0.1)
#'   to avoid over-damping on well-behaved graphs.
#' @param .lazy.alpha Numeric in (0, 1]: mixing parameter for lazy random walk.
#'   P_lazy = (1 - alpha) * I + alpha * P. Default 1 (standard walk).
#'   Values < 1 damp oscillations on sparse/irregular graphs.
#' @param .loading.method Character: method for computing feature loadings.
#'   \code{"pearson"} computes Pearson correlation between component scores
#'   and expression (scale-invariant, fast sparse-BLAS path).
#'   \code{"ols"} computes OLS regression coefficients (preserves magnitude;
#'   scale-dependent). \code{"spearman"} computes Spearman rank correlation
#'   (robust to outliers and nonlinearity; slower). Default \code{"pearson"}.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'
#' @return The modified \code{.tdr.obj} with results stored in \code{.tdr.obj$specDE[[.coef.col]]}:
#'   \describe{
#'     \item{scores}{Matrix (landmarks x components): specDE scores per landmark, ordered by Y-alignment}
#'     \item{loadings}{Matrix (genes/markers x components): gene/marker loadings per component}
#'     \item{sdev}{Numeric vector: standard deviation per component}
#'     \item{var.explained}{Numeric vector (Vk): proportion of variance explained per component}
#'     \item{Y.alignment}{Numeric vector (Ak): absolute correlation with Y per component}
#'     \item{smoothness}{Numeric vector (Sk): Laplacian smoothness per component (higher = smoother)}
#'     \item{Y}{Numeric vector: the centered density contrast used}
#'     \item{M.local}{Matrix (optional): the density-weighted expression matrix (if .store.M = TRUE)}
#'     \item{params}{List: parameters used (.model.name, .coef.col, .nv, .min.prop)}
#'   }
#'
#' @seealso \code{\link{get.lm}} (required predecessor)
#'
#' @examples
#' \dontrun{
#' # After fitting linear model
#' lm.obj <- get.lm(lm.obj, .design = design)
#' 
#' # Run specDE for "Infection" coefficient
#' lm.obj <- get.specDE(lm.obj, .coef.col = "Infection")
#' 
#' # Access results
#' lm.obj$specDE$Infection$scores[, "specDE1"]
#' lm.obj$specDE$Infection$loadings[, "specDE1"]
#' 
#' # Diagnostic table
#' data.frame(
#'   component = colnames(lm.obj$specDE$Infection$scores),
#'   Ak = lm.obj$specDE$Infection$Y.alignment,
#'   Vk = lm.obj$specDE$Infection$var.explained,
#'   Sk = lm.obj$specDE$Infection$smoothness
#' )
#' }
#'
#' @export
get.specDE <-
  function(.tdr.obj,
           .coef.col,
           .model.name = "default",
           .nv = NULL,
           .min.prop = 0.005,
           .store.M = FALSE,
           .degree.reg = FALSE,
           .tau.mult = 1,
           .lazy.alpha = 1,
           .loading.method = c("pearson", "ols", "spearman"),
           .verbose = TRUE) {
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    .loading.method <-
      match.arg(arg = .loading.method,
                choices = c("pearson",
                            "ols",
                            "spearman"))
    
    coef.mat <-
      .validate.DE.inputs(.tdr.obj = .tdr.obj,
                          .model.name = .model.name,
                          .coef.col = .coef.col)
    
    # Default .nv to number of PCs from get.landmarks
    if (is.null(x = .nv)) {
      .nv <-
        ncol(x = .tdr.obj$pca$embed)
    }
    
    if (!is.numeric(x = .nv) || length(x = .nv) != 1 || .nv < 1) {
      stop(".nv must be a positive integer.")
    }
    
    .nv <-
      as.integer(x = .nv)
    
    # -------------------------------------------------------------------------
    # Extract Y: density contrast vector (centered)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("\n=== specDE: Spectral Differential Expression ===")
      message("Coefficient: ", .coef.col)
    }
    
    Y <-
      coef.mat[, .coef.col] |>
      (\(x)
       x - mean(x = x)
      )()
    
    # -------------------------------------------------------------------------
    # Prepare expression matrix
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Preparing expression matrix...")
    }
    
    prep <-
      .prepare.X(.tdr.obj = .tdr.obj,
                 .min.prop = .min.prop,
                 .verbose = .verbose)
    
    X  <- prep$X
    Xc <- prep$Xc
    
    # -------------------------------------------------------------------------
    # Build random-walk normalized graph P
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Building random-walk normalized graph...")
    }
    
    SNN <-
      .tdr.obj$graph$snn
    
    P <-
      .build.P(SNN = SNN,
               .degree.reg = .degree.reg,
               .tau.mult = .tau.mult,
               .lazy.alpha = .lazy.alpha,
               .verbose = .verbose)
    
    # -------------------------------------------------------------------------
    # Y-X interaction and graph smoothing
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing density-weighted expression (M.local)...")
    }
    
    # Y-weighted expression: diag(Y) %*% Xc
    YX.interaction <-
      (Matrix::Diagonal(x = Y) %*% Xc) |>
      (\(x)
       `dimnames<-`(x = x,
                    value = dimnames(x = Xc))
      )()
    
    # Graph-smooth and center
    M.local <-
      (P %*% YX.interaction) |>
      (\(x)
       Matrix::t(x = x) - Matrix::colMeans(x = x)
      )() |>
      Matrix::t()
    
    n.landmarks <-
      nrow(x = M.local)
    
    # -------------------------------------------------------------------------
    # SVD decomposition (following get.landmarks pattern)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing specDE decomposition...")
    }
    
    if (.tdr.obj$config$assay.type == "RNA") {
      
      # RNA: use truncated SVD (irlba) for efficiency
      svd.res <-
        irlba::irlba(A = M.local,
                     nv = .nv,
                     nu = .nv)
      
    } else {
      
      # Cytometry: use full SVD then subset
      svd.res <-
        as.matrix(x = M.local) |>
        base::svd()
      
      .nv <-
        min(.nv, length(x = svd.res$d))
      
      svd.res$d <-
        svd.res$d[seq_len(length.out = .nv)]
      
      svd.res$u <-
        svd.res$u[, seq_len(length.out = .nv), drop = FALSE]
      
      svd.res$v <-
        svd.res$v[, seq_len(length.out = .nv), drop = FALSE]
      
    }
    
    if (isTRUE(x = .verbose)) {
      message("  -> ", .nv, " components retained")
    }
    
    # -------------------------------------------------------------------------
    # Extract scores with sign correction
    # -------------------------------------------------------------------------
    
    # Standard deviation
    sdev <-
      svd.res$d / sqrt(x = max(1, n.landmarks - 1))
    
    # Direction scores: sign correction based on P %*% Y
    dir.scores <-
      (svd.res$u %*% base::diag(x = svd.res$d)) * 
      (as.numeric(x = P %*% Y) |> sign())
    
    # -------------------------------------------------------------------------
    # Order by Y-alignment
    # -------------------------------------------------------------------------
    
    Y.cor <-
      stats::cor(x = Y,
                 y = dir.scores)[1, ]
    
    Y.cor.ord <-
      abs(x = Y.cor) |>
      order(decreasing = TRUE)
    
    # Reorder scores and apply sign so positive = aligned with Y
    scores <-
      (Matrix::t(x = dir.scores[, Y.cor.ord, drop = FALSE]) *
         sign(x = Y.cor[Y.cor.ord])) |>
      Matrix::t() |>
      as.matrix() |>
      (\(x)
       `dimnames<-`(x = x,
                    value = list(rownames(x = Xc),
                                 paste0("specDE", 
                                        seq_len(length.out = ncol(x = x)))))
      )()
    
    # Reorder sdev
    sdev <-
      sdev[Y.cor.ord]
    
    names(x = sdev) <-
      colnames(x = scores)
    
    # Y-alignment (absolute correlation, reordered)
    Y.alignment <-
      abs(x = Y.cor[Y.cor.ord])
    
    names(x = Y.alignment) <-
      colnames(x = scores)
    
    if (isTRUE(x = .verbose)) {
      message("Y-alignment (top 5): ",
              paste(sprintf("%.3f", utils::head(x = Y.alignment, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Compute loadings
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing gene/marker loadings (", .loading.method, ")...")
    }
    
    loadings <-
      .compute.loadings(scores = scores,
                        X.sparse = X,
                        .method = .loading.method)
    
    # -------------------------------------------------------------------------
    # Variance explained
    # -------------------------------------------------------------------------
    
    # Because M.local is centered by columns, compute the sample-variance-consistent Frobenius energy cheaply and exactly
    M.local.var <-
      (Matrix::colSums(x = M.local^2) |>
         sum()) / 
      (n.landmarks - 1)
    
    var.explained <-
      (sdev^2) / M.local.var
    
    names(x = var.explained) <-
      colnames(x = scores)
    
    if (isTRUE(x = .verbose)) {
      message("Variance explained (top 5): ",
              paste(sprintf("%.1f%%", utils::head(x = var.explained * 100, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Laplacian smoothness (Sk)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing Laplacian smoothness...")
    }
    
    smoothness <-
      .compute.Sk(scores = scores,
                  SNN = SNN)
    
    if (isTRUE(x = .verbose)) {
      message("Smoothness (top 5): ",
              paste(sprintf("%.3f", utils::head(x = smoothness, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Store results
    # -------------------------------------------------------------------------
    
    if (is.null(x = .tdr.obj$specDE)) {
      .tdr.obj$specDE <-
        list()
    }
    
    .tdr.obj$specDE[[.coef.col]] <-
      list(
        scores = scores,
        loadings = loadings,
        sdev = sdev,
        var.explained = var.explained,
        Y.alignment = Y.alignment,
        smoothness = smoothness,
        Y = Y,
        params = list(
          model.name = .model.name,
          coef.col = .coef.col,
          nv = .nv,
          min.prop = .min.prop,
          degree.reg = .degree.reg,
          tau.mult = .tau.mult,
          lazy.alpha = .lazy.alpha,
          loading.method = .loading.method
        )
      )
    
    if (isTRUE(x = .store.M)) {
      .tdr.obj$specDE[[.coef.col]]$M.local <-
        M.local
    }
    
    if (isTRUE(x = .verbose)) {
      message("\nResults stored in: .tdr.obj$specDE$", .coef.col)
      message("  $scores       : ", nrow(x = scores), " landmarks x ", ncol(x = scores), " components")
      message("  $loadings     : ", nrow(x = loadings), " features x ", ncol(x = loadings), " components (", .loading.method, ")")
      message("  $var.explained: Vk (proportion of variance)")
      message("  $Y.alignment  : Ak (|cor(Y, score)|)")
      message("  $smoothness   : Sk (Laplacian smoothness)")
    }
    
    return(.tdr.obj)
    
  }


#' Mass-Diffusion NMF for Density-Contrast Differential Expression
#'
#' Decomposes a nonnegative, graph-smoothed, density-weighted expression matrix
#' into additive mass programs via NMF. Unlike \code{get.specDE} (SVD on signed,
#' double-centered M.local), nmfDE operates on nonnegative mass matrices and
#' supports a literal "transported mass" interpretation.
#'
#' @details
#' nmfDE answers: "What are the dominant nonneg expression programs that carry
#' density-contrast mass across the graph?" The method:
#' \enumerate{
#'   \item Extracts the density contrast vector Y (coefficients from \code{get.lm})
#'   \item Prepares nonneg expression matrix X0 (size-factor normalized + log2 for RNA;
#'         raw for cyto; \strong{not centered})
#'   \item Builds degree-regularized, lazy random-walk graph P from SNN
#'   \item Splits Y into Y+ = max(Y,0), Y- = max(-Y,0)
#'   \item Computes mass matrices: M+ = P \%*\% diag(Y+) \%*\% X0,
#'         M- = P \%*\% diag(Y-) \%*\% X0
#'   \item Balances blocks so ||M+||_F ~ ||M-||_F (prevents trivial block separation)
#'   \item Concatenates A = \[M+ | M-\] and factorizes via \code{RcppML::nmf}: A ~ w d h
#'   \item Splits h into h+ and h- (gene templates for positive/negative contrast mass)
#'   \item Scales gene templates by d (RcppML normalizes w/h; d carries component strength)
#'   \item Derives signed gene loadings (d*h+ - d*h-) and net mass per component
#'   \item Defines polarity via cor(Y, W) (stable; does not depend on net mass)
#'   \item Orders components by |cor(Y, W * polarity)| (nmfDE1 = most Y-aligned)
#'   \item Computes loadings via regression of X0 on W (mass semantics) and
#'         of centered X on signed scores (DE semantics)
#'   \item Computes reconstruction quality via Frobenius identity (memory-efficient)
#'   \item Computes Laplacian smoothness (Sk) for each component
#' }
#'
#' \strong{Diagnostic metrics (Ak, Rk, Sk):}
#'
#' \strong{Ak (Y-alignment):} Measures how strongly component signed scores covary with
#' density contrast Y. High Ak indicates strong coupling to the tested contrast.
#' Components are ordered by Ak (nmfDE1 = most Y-aligned).
#'
#' \strong{Rk (reconstruction):} Global reconstruction quality: 1 - ||A - wdh||_F^2 / ||A||_F^2.
#' A single scalar for the full model (not per-component, since NMF components are
#' not orthogonal and variance cannot be cleanly partitioned).
#'
#' \strong{Sk (graph smoothness):} Derived from the normalized Laplacian applied to W columns.
#' High Sk indicates large-scale, graph-smooth structure. Note: because NMF components are
#' not orthogonal, Sk values across components are not independent.
#'
#' @param .tdr.obj A tinydenseR object after \code{get.lm()}.
#' @param .coef.col Character: column name in model coefficients to use as density contrast Y.
#'   Must be a valid column in \code{.tdr.obj$map$lm[[.model.name]]$fit$coefficients}.
#' @param .model.name Character: name of the fitted model to use (default "default").
#' @param .k Integer: number of NMF components. Defaults to
#'   \code{ncol(.tdr.obj$pca$embed)}, matching the number of PCs from \code{get.landmarks}.
#' @param .min.prop Numeric: for RNA, minimum proportion of landmarks where a gene must be
#'   detected (>0) to be included. Default 0.005 (0.5 percent).
#' @param .store.M Logical: if TRUE, store M.pos and M.neg in output. Default FALSE
#'   (saves memory; these can be large for RNA).
#' @param .degree.reg Logical: if TRUE, apply degree regularization by adding
#'   tau * I (self-loops) to W before row-normalizing. Default FALSE.
#' @param .tau.mult Numeric: multiplier for tau when .degree.reg = TRUE.
#'   tau = .tau.mult * mean(degree). Default 1.
#' @param .lazy.alpha Numeric in (0, 1]: mixing parameter for lazy random walk.
#'   P_lazy = (1 - alpha) * I + alpha * P. Default 1 (standard walk).
#' @param .loading.method Character: method for computing feature loadings
#'   (both mass-semantic and DE-semantic variants). \code{"pearson"} computes
#'   Pearson correlation (scale-invariant). \code{"ols"} computes OLS regression
#'   coefficients (preserves magnitude). \code{"spearman"} computes Spearman
#'   rank correlation (robust to outliers). Default \code{"pearson"}.
#' @param .L1 Numeric vector of length 2: L1/LASSO penalties for w and h in
#'   \code{RcppML::nmf}. Values in \[0, 1\]. Default c(0, 0).
#'   Increase for sparser, more interpretable programs.
#' @param .tol Numeric: convergence tolerance for \code{RcppML::nmf}. Default 1e-4.
#' @param .maxit Integer: maximum iterations for \code{RcppML::nmf}. Default 100.
#' @param .seed Integer: random seed for \code{RcppML::nmf} (passed as \code{seed} argument,
#'   not \code{set.seed}). Required for reproducibility since NMF is non-convex.
#'   Default 42.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'
#' @return The modified \code{.tdr.obj} with results stored in \code{.tdr.obj$nmfDE[[.coef.col]]}:
#'   \describe{
#'     \item{scores}{Matrix (landmarks x K): nonneg NMF activations (W), ordered by Y-alignment}
#'     \item{signed.scores}{Matrix (landmarks x K): W * polarity, oriented so positive = aligned with Y}
#'     \item{d}{Numeric vector (K): scaling diagonal from NMF}
#'     \item{h.pos}{Matrix (K x G): nonneg gene templates for positive-contrast mass}
#'     \item{h.neg}{Matrix (K x G): nonneg gene templates for negative-contrast mass}
#'     \item{gene.loadings.signed}{Matrix (G x K): derived signed gene loadings (d*h+ - d*h-),
#'       scaled by component strength}
#'     \item{loadings.mass}{Matrix (G x K): regression of nonneg X0 on W (mass semantics)}
#'     \item{loadings.de}{Matrix (G x K): regression of centered X on signed scores (DE semantics)}
#'     \item{polarity}{Integer vector (K): +1/-1 per component (sign convention)}
#'     \item{Y.alignment}{Numeric vector (K): |cor(Y, signed_score)| per component}
#'     \item{reconstruction}{Numeric scalar: 1 - ||A - wdh||_F^2 / ||A||_F^2}
#'     \item{smoothness}{Numeric vector (K): Laplacian smoothness per W column}
#'     \item{Y}{Numeric vector: the centered density contrast used}
#'     \item{block.scale}{Numeric vector of length 2: Frobenius norms of M+ and M- before balancing}
#'     \item{M.pos}{Matrix (optional): nonneg positive-contrast mass (if .store.M = TRUE)}
#'     \item{M.neg}{Matrix (optional): nonneg negative-contrast mass (if .store.M = TRUE)}
#'     \item{params}{List: parameters used}
#'   }
#'
#' @seealso \code{\link{get.lm}} (required predecessor), \code{\link{get.specDE}} (SVD variant),
#'   \code{\link{plotNmfDE}} (visualization), \code{\link{plotNmfDEHeatmap}} (heatmap)
#'
#' @examples
#' \dontrun{
#' # After fitting linear model
#' lm.obj <- get.lm(lm.obj, .design = design)
#'
#' # Run nmfDE for "Infection" coefficient
#' lm.obj <- get.nmfDE(lm.obj, .coef.col = "Infection")
#'
#' # Access results
#' lm.obj$nmfDE$Infection$scores[, "nmfDE1"]          # nonneg activations
#' lm.obj$nmfDE$Infection$signed.scores[, "nmfDE1"]    # oriented scores
#' lm.obj$nmfDE$Infection$gene.loadings.signed[, "nmfDE1"]  # h+ - h-
#'
#' # Diagnostic table
#' data.frame(
#'   component = colnames(lm.obj$nmfDE$Infection$scores),
#'   Ak = lm.obj$nmfDE$Infection$Y.alignment,
#'   Sk = lm.obj$nmfDE$Infection$smoothness
#' )
#' }
#'
#' @export
#'
get.nmfDE <-
  function(.tdr.obj,
           .coef.col,
           .model.name = "default",
           .k = NULL,
           .min.prop = 0.005,
           .store.M = FALSE,
           .degree.reg = FALSE,
           .tau.mult = 1,
           .lazy.alpha = 1,
           .loading.method = c("pearson", "ols", "spearman"),
           .L1 = c(0, 0),
           .tol = 1e-4,
           .maxit = 100L,
           .seed = 42L,
           .verbose = TRUE) {
    
    # -------------------------------------------------------------------------
    # Dependency check
    # -------------------------------------------------------------------------
    
    if (!requireNamespace("RcppML", quietly = TRUE)) {
      stop("Package 'RcppML' is required for get.nmfDE(). ",
           "Install it with: install.packages('RcppML')")
    }
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    .loading.method <-
      match.arg(arg = .loading.method,
                choices = c("pearson",
                            "ols",
                            "spearman"))
    
    coef.mat <-
      .validate.DE.inputs(.tdr.obj = .tdr.obj,
                          .model.name = .model.name,
                          .coef.col = .coef.col)
    
    # Default .k to number of PCs from get.landmarks
    if (is.null(x = .k)) {
      .k <-
        ncol(x = .tdr.obj$pca$embed)
    }
    
    if (!is.numeric(x = .k) || length(x = .k) != 1 || .k < 1) {
      stop(".k must be a positive integer.")
    }
    
    .k <-
      as.integer(x = .k)
    
    if (length(x = .L1) != 2 || any(.L1 < 0) || any(.L1 > 1)) {
      stop(".L1 must be a numeric vector of length 2 with values in [0, 1].")
    }
    
    # -------------------------------------------------------------------------
    # Extract Y: density contrast vector (centered)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("\n=== nmfDE: Mass-Diffusion NMF for Density-Contrast DE ===")
      message("Coefficient: ", .coef.col)
    }
    
    Y <-
      coef.mat[, .coef.col] |>
      (\(x)
       x - mean(x = x)
      )()
    
    # -------------------------------------------------------------------------
    # Prepare expression matrix
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Preparing expression matrix...")
    }
    
    prep <-
      .prepare.X(.tdr.obj = .tdr.obj,
                 .min.prop = .min.prop,
                 .verbose = .verbose)
    
    X0 <- prep$X
    
    n.landmarks <-
      nrow(x = X0)
    
    n.genes <-
      ncol(x = X0)
    
    gene.names <-
      colnames(x = X0)
    
    # -------------------------------------------------------------------------
    # Build random-walk normalized graph P
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Building random-walk normalized graph...")
    }
    
    SNN <-
      .tdr.obj$graph$snn
    
    P <-
      .build.P(SNN = SNN,
               .degree.reg = .degree.reg,
               .tau.mult = .tau.mult,
               .lazy.alpha = .lazy.alpha,
               .verbose = .verbose)
    
    # -------------------------------------------------------------------------
    # Split Y into positive/negative components and build mass matrices
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing nonneg mass matrices M+ and M-...")
    }
    
    Y.pos <-
      pmax(Y, 0)
    
    Y.neg <-
      pmax(-Y, 0)
    
    # M+ = P %*% diag(Y+) %*% X0: mass flowing through positive-contrast landmarks
    M.pos <-
      P %*% (Matrix::Diagonal(x = Y.pos) %*% X0)
    
    # M- = P %*% diag(Y-) %*% X0: mass flowing through negative-contrast landmarks
    M.neg <-
      P %*% (Matrix::Diagonal(x = Y.neg) %*% X0)
    
    # -------------------------------------------------------------------------
    # Block balancing: scale so ||M+||_F ~ ||M-||_F
    # Prevents NMF from trivially separating blocks by magnitude
    # -------------------------------------------------------------------------
    
    frob.pos <-
      sqrt(x = sum(M.pos^2))
    
    frob.neg <-
      sqrt(x = sum(M.neg^2))
    
    if (isTRUE(x = .verbose)) {
      message("  ||M+||_F = ", sprintf("%.2f", frob.pos),
              ", ||M-||_F = ", sprintf("%.2f", frob.neg))
    }
    
    # Scale both blocks to geometric mean of their Frobenius norms
    frob.geomean <-
      sqrt(x = frob.pos * frob.neg)
    
    if (frob.pos > .Machine$double.eps) {
      M.pos.balanced <-
        M.pos * (frob.geomean / frob.pos)
    } else {
      M.pos.balanced <-
        M.pos
    }
    
    if (frob.neg > .Machine$double.eps) {
      M.neg.balanced <-
        M.neg * (frob.geomean / frob.neg)
    } else {
      M.neg.balanced <-
        M.neg
    }
    
    if (isTRUE(x = .verbose)) {
      message("  Block-balanced to geometric mean: ", sprintf("%.2f", frob.geomean))
    }
    
    # -------------------------------------------------------------------------
    # Concatenate A = [M+_balanced | M-_balanced] and run NMF
    # -------------------------------------------------------------------------
    
    A <-
      cbind(M.pos.balanced, M.neg.balanced)
    
    # RcppML::nmf expects features-by-samples (columns = observations),
    # but our A is landmarks-by-genes (rows = landmarks).
    # Transpose so that landmarks are "samples" (columns) and genes are "features" (rows).
    A.t <-
      Matrix::t(x = A) |>
      methods::as(Class = "CsparseMatrix")
    
    if (isTRUE(x = .verbose)) {
      message("Running NMF (k = ", .k, ", tol = ", .tol,
              ", maxit = ", .maxit, ", seed = ", .seed, ")...")
      if (any(.L1 > 0)) {
        message("  L1 penalties: w = ", .L1[1], ", h = ", .L1[2])
      }
    }
    
    nmf.res <-
      RcppML::nmf(A = A.t,
                  k = .k,
                  tol = .tol,
                  maxit = .maxit,
                  verbose = FALSE,
                  L1 = .L1,
                  seed = .seed,
                  diag = TRUE,
                  nonneg = TRUE)
    
    if (isTRUE(x = .verbose)) {
      message("  NMF converged at tol = ", sprintf("%.2e", nmf.res$tol),
              " after ", nmf.res$iter, " iterations")
    }
    
    # RcppML::nmf returns:  A.t ~ w %*% diag(d) %*% h
    #   w: (2G x K) — feature factor matrix (rows = genes in [pos|neg] block, cols = components)
    #   d: (K)      — scaling diagonal
    #   h: (K x L)  — sample factor matrix (rows = components, cols = landmarks)
    #
    # To get landmark scores (L x K) and gene templates (K x G):
    #   W = t(h)         — landmark activations (L x K), nonneg
    #   H_full = t(w)    — gene templates in concatenated space (K x 2G)
    #   H_pos = H_full[, 1:G]       — positive-block gene templates
    #   H_neg = H_full[, (G+1):2G]  — negative-block gene templates
    
    W <-
      t(x = nmf.res$h)
    
    nmf.d <-
      nmf.res$d
    
    H.full <-
      t(x = nmf.res$w)
    
    H.pos <-
      H.full[, seq_len(length.out = n.genes), drop = FALSE]
    
    H.neg <-
      H.full[, n.genes + seq_len(length.out = n.genes), drop = FALSE]
    
    colnames(x = H.pos) <-
      gene.names
    
    colnames(x = H.neg) <-
      gene.names
    
    if (isTRUE(x = .verbose)) {
      message("  W: ", nrow(x = W), " landmarks x ", ncol(x = W), " components")
      message("  H+, H-: ", nrow(x = H.pos), " components x ", ncol(x = H.pos), " features each")
    }
    
    # -------------------------------------------------------------------------
    # Scale templates by d to get effective component strength
    # RcppML::nmf(diag=TRUE) normalizes w/h to sum to 1, pushing all scale
    # into d. Without scaling, net.mass and gene.loadings.signed would treat
    # a dominant and negligible component equally.
    # -------------------------------------------------------------------------
    
    H.pos.scaled <-
      as.matrix(x = Matrix::Diagonal(x = nmf.d) %*% H.pos)
    
    H.neg.scaled <-
      as.matrix(x = Matrix::Diagonal(x = nmf.d) %*% H.neg)
    
    # -------------------------------------------------------------------------
    # Derive polarity and signed scores
    # -------------------------------------------------------------------------
    
    # Signed gene loadings: (d * h+) - (d * h-) (derived, not a factor of the nonneg model)
    H.signed <-
      H.pos.scaled - H.neg.scaled
    
    # Net mass per component: ||d_k * h+_k||_1 - ||d_k * h-_k||_1
    # Positive net mass -> component carries more positive-contrast expression
    net.mass <-
      rowSums(x = H.pos.scaled) - rowSums(x = H.neg.scaled)
    
    # Polarity: use direct cor(Y, W) -- more stable than cor(Y, W*net.mass)
    # When net.mass approx 0 the product W*net.mass -> 0 and correlation becomes noise.
    # cor(Y, W) is always well-defined since W columns are nonneg with variance > 0.
    Y.cor.W <-
      stats::cor(x = Y,
                 y = W,
                 method = "spearman")[1, ]
    
    polarity <-
      sign(x = Y.cor.W) |>
      as.integer()
    
    # Handle zero correlation (default to +1)
    polarity[polarity == 0L] <-
      1L
    
    # Signed scores: W * polarity (oriented by Y-alignment)
    # net.mass is stored as a separate descriptor, not baked into scores
    signed.scores <-
      t(x = t(x = W) * polarity)
    
    # -------------------------------------------------------------------------
    # Order by Y-alignment (|cor(Y, signed.scores)|, descending)
    # -------------------------------------------------------------------------
    
    Y.alignment.raw <-
      abs(x = Y.cor.W)
    
    comp.ord <-
      order(Y.alignment.raw, decreasing = TRUE)
    
    # Reorder all per-component outputs
    W <-
      W[, comp.ord, drop = FALSE]
    
    signed.scores <-
      signed.scores[, comp.ord, drop = FALSE]
    
    nmf.d <-
      nmf.d[comp.ord]
    
    H.pos <-
      H.pos[comp.ord, , drop = FALSE]
    
    H.neg <-
      H.neg[comp.ord, , drop = FALSE]
    
    H.signed <-
      H.signed[comp.ord, , drop = FALSE]
    
    polarity <-
      polarity[comp.ord]
    
    net.mass <-
      net.mass[comp.ord]
    
    Y.alignment <-
      Y.alignment.raw[comp.ord]
    
    # Set names
    comp.names <-
      paste0("nmfDE", seq_len(length.out = .k))
    
    colnames(x = W) <-
      comp.names
    
    rownames(x = W) <-
      rownames(x = X0)
    
    colnames(x = signed.scores) <-
      comp.names
    
    rownames(x = signed.scores) <-
      rownames(x = X0)
    
    names(x = nmf.d) <-
      comp.names
    
    rownames(x = H.pos) <-
      comp.names
    
    rownames(x = H.neg) <-
      comp.names
    
    rownames(x = H.signed) <-
      comp.names
    
    names(x = polarity) <-
      comp.names
    
    names(x = net.mass) <-
      comp.names
    
    names(x = Y.alignment) <-
      comp.names
    
    if (isTRUE(x = .verbose)) {
      message("Y-alignment (top 5): ",
              paste(sprintf("%.3f", utils::head(x = Y.alignment, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Compute loadings (two semantic variants)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing gene/marker loadings (", .loading.method, ")...")
    }
    
    # Variant 1: Mass semantics -- loadings of nonneg X0 on W
    # Interpretation: genes whose expression increases with component mass activation
    loadings.mass <-
      .compute.loadings(scores = W,
                        X.sparse = X0,
                        .method = .loading.method)
    
    # Variant 2: DE semantics -- loadings of centered X on signed scores
    # Interpretation: genes associated with net positive vs negative mass along contrast
    loadings.de <-
      .compute.loadings(scores = signed.scores,
                        X.sparse = X0,
                        .method = .loading.method)
    
    # -------------------------------------------------------------------------
    # Reconstruction quality (global, not per-component)
    # Uses Frobenius identity to avoid materializing the dense (2G x L) A-hat:
    #   ||A - A-hat||^2 = ||A||^2 - 2 tr(A^T A-hat) + ||A-hat||^2
    # All terms computed via KxK intermediaries + one (LxK) product.
    # -------------------------------------------------------------------------
    
    # ||A.t||^2_F -- efficient for sparse dgCMatrix (sum of squared nonzero slots)
    A.t.frob2 <-
      sum(A.t@x^2)
    
    # tr(A.t^T A-hat) = tr(A.t^T w D h) = sum(Q * t(Dh))
    #   where Q = A.t^T w (LxK), Dh = diag(d) h (KxL)
    Q <-
      Matrix::crossprod(x = A.t, y = nmf.res$w)  # L x K
    
    Dh <-
      nmf.res$d * nmf.res$h  # K x L (d recycled over rows)
    
    trAtAhat <-
      sum(Q * Matrix::t(x = Dh))  # sum of element-wise product, L x K
    
    # ||A-hat||^2_F = sum_{ij} ((wDh)_{ij})^2 via KxK Hadamard:
    #   = sum( (d_i d_j w^Tw_{ij}) * (hh^T)_{ij} )
    WtW <-
      Matrix::crossprod(x = nmf.res$w)  # K x K
    
    DWtWD <-
      (nmf.res$d %o% nmf.res$d) * as.matrix(x = WtW)  # K x K
    
    HHt <-
      Matrix::tcrossprod(x = nmf.res$h)  # K x K
    
    Ahat.frob2 <-
      sum(DWtWD * as.matrix(x = HHt))  # trace via Hadamard product
    
    # ||A.t - A-hat||^2_F = ||A.t||^2 - 2 tr(A^T A-hat) + ||A-hat||^2
    resid.frob2 <-
      A.t.frob2 - 2 * trAtAhat + Ahat.frob2
    
    reconstruction <-
      if (A.t.frob2 > .Machine$double.eps) {
        1 - max(resid.frob2, 0) / A.t.frob2
      } else {
        NA_real_
      }
    
    if (isTRUE(x = .verbose)) {
      message("Reconstruction quality (1 - ||residual||^2/||A||^2): ",
              sprintf("%.4f", reconstruction))
    }
    
    # -------------------------------------------------------------------------
    # Laplacian smoothness (Sk) -- computed on W columns
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing Laplacian smoothness...")
    }
    
    smoothness <-
      .compute.Sk(scores = W,
                  SNN = SNN)
    
    if (isTRUE(x = .verbose)) {
      message("Smoothness (top 5): ",
              paste(sprintf("%.3f", utils::head(x = smoothness, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Store results
    # -------------------------------------------------------------------------
    
    if (is.null(x = .tdr.obj$nmfDE)) {
      .tdr.obj$nmfDE <-
        list()
    }
    
    .tdr.obj$nmfDE[[.coef.col]] <-
      list(
        scores = W,
        signed.scores = signed.scores,
        d = nmf.d,
        h.pos = H.pos,
        h.neg = H.neg,
        gene.loadings.signed = t(x = H.signed),
        loadings.mass = loadings.mass,
        loadings.de = loadings.de,
        polarity = polarity,
        net.mass = net.mass,
        Y.alignment = Y.alignment,
        reconstruction = reconstruction,
        smoothness = smoothness,
        Y = Y,
        block.scale = c(pos = frob.pos, neg = frob.neg),
        params = list(
          model.name = .model.name,
          coef.col = .coef.col,
          k = .k,
          min.prop = .min.prop,
          degree.reg = .degree.reg,
          tau.mult = .tau.mult,
          lazy.alpha = .lazy.alpha,
          loading.method = .loading.method,
          L1 = .L1,
          tol = .tol,
          maxit = .maxit,
          seed = .seed
        )
      )
    
    if (isTRUE(x = .store.M)) {
      .tdr.obj$nmfDE[[.coef.col]]$M.pos <-
        M.pos
      .tdr.obj$nmfDE[[.coef.col]]$M.neg <-
        M.neg
    }
    
    if (isTRUE(x = .verbose)) {
      message("\nResults stored in: .tdr.obj$nmfDE$", .coef.col)
      message("  $scores              : ", nrow(x = W), " landmarks x ", .k, " components (nonneg W)")
      message("  $signed.scores       : ", nrow(x = W), " landmarks x ", .k, " components (oriented)")
      message("  $d                   : scaling diagonal (", .k, ")")
      message("  $h.pos, $h.neg       : ", .k, " x ", n.genes, " gene templates (nonneg)")
      message("  $gene.loadings.signed: ", n.genes, " x ", .k, " (h+ - h-, derived)")
      message("  $loadings.mass       : ", n.genes, " x ", .k, " (", .loading.method, ", X0 ~ W)")
      message("  $loadings.de         : ", n.genes, " x ", .k, " (", .loading.method, ", Xc ~ signed.scores)")
      message("  $Y.alignment         : Ak (|cor(Y, signed.score)|)")
      message("  $reconstruction      : Rk (global, 1 - ||resid||^2/||A||^2)")
      message("  $smoothness          : Sk (Laplacian smoothness on W)")
    }
    
    return(.tdr.obj)
    
  }


#' Graph-Diffused, Density Contrast-Aligned PLS for Differential Expression
#'
#' Decomposes an expression interaction matrix M.local via NIPALS PLS1
#' against the density-contrast vector Y. plsDE identifies features that drive
#' the density contrast — including population markers (DA), differentially
#' expressed genes (DE), and their mixture — by maximizing covariance between
#' graph-smoothed expression and Y.
#'
#' @details
#' plsDE is an \strong{interpretive decomposition}, not a formal hypothesis test.
#' It answers: "What features — through expression patterns, population identity,
#' or both — explain the density contrast?" Loadings reflect this combined signal:
#' markers of depleted populations carry negative loadings (high expression where
#' density is low), DE genes carry signed loadings tracking their direction of
#' change, and mixed DA/DE features appear alongside both.
#'
#' When \code{.YX.interaction = TRUE} (default), the interaction term diag(Y)
#' bakes the density contrast into the data matrix:
#' M.local = P \%*\% diag(Y) \%*\% Xc. This gives plsDE comprehensive scope:
#' it captures DA markers, DE genes, and their interplay in a single
#' decomposition. The interaction amplifies expression signal in landmarks
#' with strong density contrast, making the method sensitive even to subtle DE
#' in the absence of DA.
#'
#' When \code{.YX.interaction = FALSE}, M.local = P \%*\% Xc (graph-smoothed
#' centered expression without Y-weighting). Y appears only on the response
#' side, so loadings reflect only features whose graph-smoothed expression
#' independently covaries with Y. This mode is less comprehensive — it misses
#' DA markers unless their expression happens to correlate with Y through
#' smoothed space — but serves as a useful diagnostic comparator: features
#' that rank high in both modes are robust.
#'
#' The method:
#' \enumerate{
#'   \item Extracts the density contrast vector Y (coefficients from \code{get.lm})
#'   \item Prepares centered expression matrix Xc (size-factor normalized + log2 for RNA;
#'         raw for cyto; then centered)
#'   \item Builds random-walk normalized graph P from SNN (with optional degree
#'         regularization and lazy walk)
#'   \item Constructs M.local: if \code{.YX.interaction = TRUE},
#'         M.local = P \%*\% diag(Y) \%*\% Xc (density-weighted); if FALSE,
#'         M.local = P \%*\% Xc (graph-smoothed only). Both are centered by columns.
#'   \item Runs NIPALS PLS1: iteratively finds feature weights w maximizing
#'         cov(M.local w, Y), with deflation of both M.local and Y
#'   \item Applies sign convention: positive scores = aligned with Y
#'   \item Computes feature loadings (Pearson r, OLS regression, or Spearman rank
#'         correlation)
#'   \item Computes Laplacian smoothness (Sk) for each component
#' }
#'
#' \strong{Diagnostic metrics (Ak, Sk):}
#'
#' \strong{Ak (Y-alignment):} Measures how strongly component scores covary with
#' density contrast Y. Ak is high by construction for the leading components:
#' NIPALS PLS1 maximizes covariance with Y at every deflation step. Components
#' are ordered by extraction (plsDE1 = highest covariance with Y).
#'
#' \strong{Sk (graph smoothness):} Derived from the normalized Laplacian applied to
#' score vectors. High Sk indicates large-scale, graph-smooth structure.
#'
#' \strong{About Y appearing on both sides (when .YX.interaction = TRUE):}
#' Y appears in both the data matrix (via diag(Y) in M.local) and the PLS objective
#' (maximize covariance with Y). This is intentional: it gives plsDE comprehensive
#' sensitivity to features driving density changes. The design is not circular
#' because the expression matrix Xc sits between: the product is large only when
#' features exist whose Y-weighted, graph-smoothed expression genuinely covaries
#' with Y. With permuted Y, Ak collapses. Setting \code{.YX.interaction = FALSE}
#' removes Y from the data matrix entirely, providing a diagnostic comparator.
#'
#' @note plsDE is designed as an exploratory tool to help interpret density
#'   changes in terms of the features driving them. It does not provide
#'   gene-level p-values or formal multiple testing correction. For rigorous
#'   differential expression testing with p-values following field standards,
#'   use \code{\link{get.pbDE}} (pseudobulk DE via edgeR/limma) or
#'   \code{\link{get.markerDE}} (marker-level DE). plsDE complements these
#'   methods by providing a multivariate, graph-aware decomposition that
#'   captures joint DA/DE patterns not visible in gene-by-gene tests.
#'
#' @param .tdr.obj A tinydenseR object after \code{get.lm()}.
#' @param .coef.col Character: column name in model coefficients to use as density contrast Y.
#'   Must be a valid column in \code{.tdr.obj$map$lm[[.model.name]]$fit$coefficients}.
#' @param .model.name Character: name of the fitted model to use (default "default").
#' @param .ncomp Integer: number of PLS components. Defaults to
#'   \code{ncol(.tdr.obj$pca$embed)}, matching the number of PCs from \code{get.landmarks}.
#' @param .min.prop Numeric: for RNA, minimum proportion of landmarks where a gene must be
#'   detected (>0) to be included. Default 0.005 (0.5 percent).
#' @param .store.M Logical: if TRUE, store M.local in output. Default FALSE
#'   (saves memory; can be large for RNA).
#' @param .degree.reg Logical: if TRUE, apply degree regularization by adding
#'   tau * I (self-loops) to SNN before row-normalizing. Default FALSE.
#' @param .tau.mult Numeric: multiplier for tau when .degree.reg = TRUE.
#'   tau = .tau.mult * mean(degree). Default 1.
#' @param .lazy.alpha Numeric in (0, 1]: mixing parameter for lazy random walk.
#'   P_lazy = (1 - alpha) * I + alpha * P. Default 1 (standard walk).
#' @param .YX.interaction Logical: if TRUE (default), construct
#'   M.local = P diag(Y) Xc (Y-weighted interaction). Loadings capture features
#'   driving the density contrast through both differential abundance (population
#'   markers) and differential expression. If FALSE, construct M.local = P Xc
#'   (graph-smoothed expression only; Y appears only on the response side).
#'   The FALSE mode is less comprehensive — it misses DA markers unless their
#'   expression independently covaries with Y.
#' @param .loading.method Character: method for computing feature loadings. 
#'   \code{"pearson"} computes Pearson correlation between component scores
#'   and centered expression (scale-invariant, fast sparse-BLAS path).
#'   \code{"ols"} computes OLS regression coefficients of centered Xc on component
#'   scores (same sparse-BLAS path; preserves magnitude information but is
#'   scale-dependent, so high-variance features rank higher). \code{"spearman"}
#'   computes Spearman rank correlation via a sparse-aware implementation
#'   (robust to outliers and nonlinearity; slower, O(nnz * log(m_g) * K)).
#' @param .verbose Logical: print progress messages? Default TRUE.
#'
#' @return The modified \code{.tdr.obj} with results stored in \code{.tdr.obj$plsDE[[.coef.col]]}:
#'   \describe{
#'     \item{scores}{Matrix (landmarks x K): PLS scores (oriented so positive = aligned with Y)}
#'     \item{feature.weights}{Matrix (genes x K): PLS feature weight vectors w (unit-norm)}
#'     \item{x.loadings}{Matrix (genes x K): PLS deflation loadings p}
#'     \item{y.loadings}{Numeric vector (K): PLS Y-loadings q (scalar per component)}
#'     \item{loadings}{Matrix (genes x K): feature loadings per component (Pearson r,
#'       OLS regression, or Spearman rank correlation, depending on
#'       \code{.loading.method}); positive = upregulated with Y,
#'       negative = downregulated}
#'     \item{Y.alignment}{Numeric vector (K): |cor(Y, score_k)| per component}
#'     \item{smoothness}{Numeric vector (K): Laplacian smoothness per score vector}
#'     \item{Y}{Numeric vector: the centered density contrast used}
#'     \item{M.local}{Matrix (optional): the interaction matrix (if .store.M = TRUE)}
#'     \item{params}{List: parameters used}
#'   }
#'
#' @seealso \code{\link{get.lm}} (required predecessor),
#'   \code{\link{get.specDE}} (specDE: SVD-based peer method),
#'   \code{\link{get.nmfDE}} (nmfDE: NMF-based peer method),
#'   \code{\link{plotPlsDE}} (visualization), \code{\link{plotPlsDEHeatmap}} (heatmap)
#'
#' @examples
#' \dontrun{
#' # After fitting linear model
#' lm.obj <- get.lm(lm.obj, .design = design)
#'
#' # Run plsDE for "Infection" coefficient
#' lm.obj <- get.plsDE(lm.obj, .coef.col = "Infection")
#'
#' # Access results
#' lm.obj$plsDE$Infection$scores[, "plsDE1"]      # PLS scores (Y-aligned)
#' lm.obj$plsDE$Infection$loadings[, "plsDE1"]     # gene loadings (Pearson, OLS or Spearman)
#'
#' # Diagnostic table
#' data.frame(
#'   component = colnames(lm.obj$plsDE$Infection$scores),
#'   Ak = lm.obj$plsDE$Infection$Y.alignment,
#'   Sk = lm.obj$plsDE$Infection$smoothness
#' )
#' }
#'
#' @export
#'
get.plsDE <-
  function(.tdr.obj,
           .coef.col,
           .model.name = "default",
           .ncomp = NULL,
           .min.prop = 0.005,
           .store.M = FALSE,
           .degree.reg = FALSE,
           .tau.mult = 1,
           .lazy.alpha = 1,
           .YX.interaction = TRUE,
           .loading.method = c("pearson", "ols", "spearman"),
           .verbose = TRUE) {
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    if (!is.logical(x = .YX.interaction) || length(x = .YX.interaction) != 1) {
      stop(".YX.interaction must be TRUE or FALSE.")
    }

    .loading.method <-
      match.arg(arg = .loading.method,
                choices = c("pearson",
                            "ols",
                            "spearman"))
    
    coef.mat <-
      .validate.DE.inputs(.tdr.obj = .tdr.obj,
                          .model.name = .model.name,
                          .coef.col = .coef.col)
    
    # Default .ncomp to number of PCs from get.landmarks
    if (is.null(x = .ncomp)) {
      .ncomp <-
        ncol(x = .tdr.obj$pca$embed)
    }
    
    if (!is.numeric(x = .ncomp) || length(x = .ncomp) != 1 || .ncomp < 1) {
      stop(".ncomp must be a positive integer.")
    }
    
    .ncomp <-
      as.integer(x = .ncomp)
    
    # -------------------------------------------------------------------------
    # Extract Y: density contrast vector (centered)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("\n=== plsDE: Graph-Diffused, Density Contrast-Aligned PLS for Differential Expression ===")
      message("Coefficient: ", .coef.col)
    }
    
    Y <-
      coef.mat[, .coef.col] |>
      (\(x)
       x - mean(x = x)
      )()
    
    # -------------------------------------------------------------------------
    # Prepare expression matrix
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Preparing expression matrix...")
    }
    
    prep <-
      .prepare.X(.tdr.obj = .tdr.obj,
                 .min.prop = .min.prop,
                 .verbose = .verbose)
    
    X  <- prep$X
    Xc <- prep$Xc
    
    n.landmarks <-
      nrow(x = Xc)
    
    n.genes <-
      ncol(x = Xc)
    
    gene.names <-
      colnames(x = Xc)
    
    # -------------------------------------------------------------------------
    # Build random-walk normalized graph P
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Building random-walk normalized graph...")
    }
    
    SNN <-
      .tdr.obj$graph$snn
    
    P <-
      .build.P(SNN = SNN,
               .degree.reg = .degree.reg,
               .tau.mult = .tau.mult,
               .lazy.alpha = .lazy.alpha,
               .verbose = .verbose)
    
    # -------------------------------------------------------------------------
    # Construct M.local
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .YX.interaction)) {
      
      # M.local = P %*% diag(Y) %*% Xc (density-weighted interaction)
      if (isTRUE(x = .verbose)) {
        message("Computing density-weighted expression (M.local = P diag(Y) Xc)...")
      }
      
      YX.term <-
        (Matrix::Diagonal(x = Y) %*% Xc) |>
        (\(x)
         `dimnames<-`(x = x,
                      value = dimnames(x = Xc))
        )()
      
    } else {
      
      # M.local = P %*% Xc (graph-smoothed expression, no Y-weighting)
      if (isTRUE(x = .verbose)) {
        message("Computing graph-smoothed expression (M.local = P Xc)...")
      }
      
      YX.term <- Xc
      
    }
    
    # Graph-smooth and center
    M.local <-
      (P %*% YX.term) |>
      (\(x)
       Matrix::t(x = x) - Matrix::colMeans(x = x)
      )() |>
      Matrix::t()
    
    # -------------------------------------------------------------------------
    # NIPALS PLS1: M.local vs Y
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Running NIPALS PLS1 (", .ncomp, " components)...")
    }
    
    # Dense copy for deflation (deflation breaks sparsity on first step)
    Z.work <-
      as.matrix(x = M.local)
    
    Y.work <-
      Y
    
    # Pre-allocate storage
    W.mat <-
      matrix(data = 0,
             nrow = n.genes,
             ncol = .ncomp)
    
    T.mat <-
      matrix(data = 0,
             nrow = n.landmarks,
             ncol = .ncomp)
    
    P.mat <-
      matrix(data = 0,
             nrow = n.genes,
             ncol = .ncomp)
    
    Q.vec <-
      numeric(length = .ncomp)
    
    for (k in seq_len(length.out = .ncomp)) {
      
      # Gene weight: w = Z'Y / ||Z'Y||
      ZtY <-
        crossprod(x = Z.work,
                  y = Y.work) |>
        as.numeric()
      
      w.norm <-
        sqrt(x = sum(ZtY^2))
      
      if (w.norm < .Machine$double.eps) {
        if (isTRUE(x = .verbose)) {
          message("  Component ", k, ": Z'Y norm ~ 0, stopping early.")
        }
        .ncomp <- k - 1L
        break
      }
      
      w <-
        ZtY / w.norm
      
      # Score: t = Z w
      t.score <-
        (Z.work %*% w) |>
        as.numeric()
      
      # Deflation loading: p = Z't / (t't)
      tt <-
        sum(t.score^2)
      
      p.load <-
        (crossprod(x = Z.work,
                   y = t.score) |>
           as.numeric()) / tt
      
      # Y loading: q = Y't / (t't)
      q.load <-
        sum(Y.work * t.score) / tt
      
      # Deflate Z and Y
      Z.work <-
        Z.work - tcrossprod(x = t.score,
                            y = p.load)
      
      Y.work <-
        Y.work - t.score * q.load
      
      # Store
      W.mat[, k] <-
        w
      
      T.mat[, k] <-
        t.score
      
      P.mat[, k] <-
        p.load
      
      Q.vec[k] <-
        q.load
      
      if (isTRUE(x = .verbose)) {
        Ak <-
          abs(x = stats::cor(x = Y, y = t.score))
        message(sprintf("  plsDE%d: cov = %.2f, Ak = %.4f, q = %.4f, resid.var(Y) = %.4f",
                        k,
                        abs(x = sum(t.score * Y)),
                        Ak,
                        q.load,
                        stats::var(x = Y.work) / stats::var(x = Y)))
      }
    }
    
    # Trim if stopped early
    if (.ncomp < ncol(x = W.mat)) {
      W.mat <- W.mat[, seq_len(length.out = .ncomp), drop = FALSE]
      T.mat <- T.mat[, seq_len(length.out = .ncomp), drop = FALSE]
      P.mat <- P.mat[, seq_len(length.out = .ncomp), drop = FALSE]
      Q.vec <- Q.vec[seq_len(length.out = .ncomp)]
    }
    
    # -------------------------------------------------------------------------
    # Sign convention: positive scores = aligned with Y
    # -------------------------------------------------------------------------
    
    Y.cor <-
      stats::cor(x = Y,
                 y = T.mat)[1, ]
    
    sign.flip <-
      sign(x = Y.cor)
    
    # Handle zero correlation (default to +1)
    sign.flip[sign.flip == 0] <-
      1
    
    # Flip scores, weights, and loadings
    T.mat <-
      t(x = t(x = T.mat) * sign.flip)
    
    W.mat <-
      t(x = t(x = W.mat) * sign.flip)
    
    P.mat <-
      t(x = t(x = P.mat) * sign.flip)
    
    Q.vec <-
      Q.vec * sign.flip
    
    Y.alignment <-
      abs(x = Y.cor)
    
    # -------------------------------------------------------------------------
    # Name everything
    # -------------------------------------------------------------------------
    
    comp.names <-
      paste0("plsDE", seq_len(length.out = .ncomp))
    
    colnames(x = W.mat) <-
      colnames(x = T.mat) <-
      colnames(x = P.mat) <-
      names(x = Q.vec) <-
      names(x = Y.alignment) <-
      comp.names
    
    rownames(x = W.mat) <-
      rownames(x = P.mat) <-
      gene.names
    
    rownames(x = T.mat) <-
      rownames(x = Xc)
    
    if (isTRUE(x = .verbose)) {
      message("Y-alignment (top 5): ",
              paste(sprintf("%.3f", utils::head(x = Y.alignment, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Compute loadings
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing gene/marker loadings (", .loading.method, ")...")
    }
    
    loadings <-
      .compute.loadings(scores = T.mat,
                        X.sparse = X,
                        .method = .loading.method)
    
    # -------------------------------------------------------------------------
    # Laplacian smoothness (Sk)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing Laplacian smoothness...")
    }
    
    smoothness <-
      .compute.Sk(scores = T.mat,
                  SNN = SNN)
    
    if (isTRUE(x = .verbose)) {
      message("Smoothness (top 5): ",
              paste(sprintf("%.3f", utils::head(x = smoothness, n = 5)), collapse = ", "))
    }
    
    # -------------------------------------------------------------------------
    # Store results
    # -------------------------------------------------------------------------
    
    if (is.null(x = .tdr.obj$plsDE)) {
      .tdr.obj$plsDE <-
        list()
    }
    
    .tdr.obj$plsDE[[.coef.col]] <-
      list(
        scores = T.mat,
        gene.weights = W.mat,
        x.loadings = P.mat,
        y.loadings = Q.vec,
        loadings = loadings,
        Y.alignment = Y.alignment,
        smoothness = smoothness,
        Y = Y,
        params = list(
          model.name = .model.name,
          coef.col = .coef.col,
          ncomp = .ncomp,
          min.prop = .min.prop,
          degree.reg = .degree.reg,
          tau.mult = .tau.mult,
          lazy.alpha = .lazy.alpha,
          YX.interaction = .YX.interaction,
          loading.method = .loading.method
        )
      )
    
    if (isTRUE(x = .store.M)) {
      .tdr.obj$plsDE[[.coef.col]]$M.local <-
        M.local
    }
    
    if (isTRUE(x = .verbose)) {
      message("\nResults stored in: .tdr.obj$plsDE$", .coef.col)
      message("  $scores       : ", n.landmarks, " landmarks x ", .ncomp, " components")
      message("  $gene.weights : ", n.genes, " features x ", .ncomp, " components (PLS w)")
      message("  $loadings     : ", n.genes, " features x ", .ncomp, " components (", .loading.method, ")")
      message("  $Y.alignment  : Ak (|cor(Y, score)|)")
      message("  $smoothness   : Sk (Laplacian smoothness)")
    }
    
    return(.tdr.obj)
    
  }


# =============================================================================
# Internal: Sparse-aware Spearman correlation
# =============================================================================

# Spearman rank correlation between a dense score vector and each column of a
# sparse matrix.  Exploits the fact that after log2(x+1) transformation,
# structural zeros are the column minimum.  Their midranks are known
# analytically, so only the nonzero entries need sorting.
#
# Args:
#   s        - numeric vector (length n): dense score (e.g. PLS component)
#   X.sparse - dgCMatrix (n x p): sparse expression matrix where structural
#              zeros are the column-minimum value (pre-centering)
#
# Value:
#   numeric vector (length p): Spearman rho for each column vs s

.sparse.spearman <-
  function(s, X.sparse) {
    
    if(is.na(x = X.sparse) |>
       any()) {
      stop("NAs found in normalized expression matrix")
    }
    
    X.sparse <-
      Matrix::drop0(x = X.sparse)
    
    if(min(X.sparse) < 0) {
      stop("X.sparse contains zeros")
    }
    
    n <- length(x = s)
    mid <- (n + 1) / 2
    
    # Rank the score vector once -- O(n log n)
    r.s <- rank(x = s)
    r.s.c <- r.s - mid                        # centered ranks
    ss.rs <- sum(r.s.c^2)                      # score-side denominator (scalar)
    
    # CSC slot access
    p.ptr <- X.sparse@p
    x.vals <- X.sparse@x
    i.idx <- X.sparse@i                        # 0-based row indices
    
    n.genes <- ncol(x = X.sparse)
    rho <- numeric(length = n.genes)
    
    for (g in seq_len(length.out = n.genes)) {
      
      idx.start <- p.ptr[g] + 1L               # 1-based start in x.vals / i.idx
      idx.end <- p.ptr[g + 1L]                 # 1-based end (inclusive)
      m.g <- idx.end - idx.start + 1L          # number of nonzeros
      z.g <- n - m.g                            # number of structural zeros
      
      if (m.g == 0L) next                       # all-zero gene: rho stays 0
      
      # Midrank for the zero group
      r0 <- (1 + z.g) / 2
      
      # Ranks of nonzero values within themselves, shifted to full-vector position
      nz.idx <- idx.start:idx.end
      nz.ranks <- rank(x = x.vals[nz.idx]) + z.g     # z.g+1 .. n
      
      # Centered score ranks for the nonzero rows
      nz.rows <- i.idx[nz.idx] + 1L            # 1-based row indices
      nz.r.s.c <- r.s.c[nz.rows]
      
      # Numerator:
      #   sum_all r.s.c * r.g.c
      # = sum_nnz r.s.c * (nz.rank - mid) + sum_zero r.s.c * (r0 - mid)
      # Since sum_all r.s.c = 0:  sum_zero r.s.c = -sum_nnz r.s.c
      sum.nz.rs <- sum(nz.r.s.c)
      num <- sum(nz.r.s.c * (nz.ranks - mid)) +
        (r0 - mid) * (-sum.nz.rs)
      
      # Denominator (gene-side sum of squared centered ranks)
      ss.rg <- z.g * (r0 - mid)^2 +
        sum((nz.ranks - mid)^2)
      
      denom <- sqrt(x = ss.rs * ss.rg)
      
      if (denom > 0) rho[g] <- num / denom
    }
    
    rho
  }

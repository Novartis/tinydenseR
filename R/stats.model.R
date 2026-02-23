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
        lapply(X = .tdr.obj$map$fuzzy.graph,
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
            .tdr.obj$map$fuzzy.graph[[smpl]][.lm.idx[[smpl]],,drop = FALSE]
          
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
            .tdr.obj$map$fuzzy.graph[[smpl]][.lm.idx[[smpl]],,drop=FALSE]
          
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
        lapply(X = .tdr.obj$map$fuzzy.graph,
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
        Map(f = function(n, idx1){
          setdiff(x = seq_len(length.out = n),
                  y = idx1)
        },
        n = .tdr.obj$metadata$n.cells,
        idx1 = .lm1.idx)
      
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
          lapply(X = .tdr.obj$map$fuzzy.graph,
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
                  .tdr.obj$map$fuzzy.graph[[smpl]][cell.idx,,drop = FALSE]
                
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
                  .tdr.obj$map$fuzzy.graph[[smpl]][cell.idx,,drop = FALSE]
                
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
                  .tdr.obj$map$fuzzy.graph[[smpl]][cell.idx,,drop = FALSE]
                
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
                  .tdr.obj$map$fuzzy.graph[[smpl]][cell.idx,,drop = FALSE]
                
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
    .dist.metric = "cosine",
    .verbose = TRUE
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
          distance = .dist.metric
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
#' @param .dist.metric Character: distance metric for diffusion map (default "cosine").
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
    .dist.metric = "cosine",
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
        .dist.metric = .dist.metric,
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
           .verbose = TRUE) {
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
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
    
    if(!Matrix::isSymmetric(object = .tdr.obj$graph$snn)){
      stop("SNN graph not symmetric.")
    }
    
    if (is.null(x = .tdr.obj$raw.landmarks)) {
      stop("Raw landmarks not found. Run get.landmarks() first.")
    }
    
    if (is.null(x = .tdr.obj$pca$embed)) {
      stop("PCA embedding not found. Run get.landmarks() first.")
    }
    
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
    # Prepare expression matrix X
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Preparing expression matrix...")
    }
    
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
    X <-
      X |>
      (\(x)
       Matrix::t(x = x) - Matrix::colMeans(x = x)
      )() |>
      Matrix::t()
    
    # -------------------------------------------------------------------------
    # Build random-walk normalized graph P (with optional degree regularization + lazy walk)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Building random-walk normalized graph...")
    }
    
    # Validate lazy.alpha
    if (!is.numeric(x = .lazy.alpha) || .lazy.alpha <= 0 || .lazy.alpha > 1) {
      stop(".lazy.alpha must be in (0, 1]")
    }
    
    SNN <-
      .tdr.obj$graph$snn
    
    # Check for disconnected nodes
    d <-
      Matrix::rowSums(x = SNN)
    
    if (any(d == 0)) {
      n.disconnected <- sum(d == 0)
      stop("Graph contains ", n.disconnected, " disconnected node(s) (degree = 0). ",
           "Cannot compute random-walk matrix. Check graph construction or increase k.")
    }
    
    # Degree regularization: add self-loops W_tau = W + tau*I, then row-normalize
    # This keeps P as a proper Markov operator (rowSums = 1) while stabilizing low-degree nodes
    if (isTRUE(x = .degree.reg)) {
      tau <- .tau.mult * mean(x = d)
      W.tau <- SNN + tau * Matrix::Diagonal(n = nrow(x = SNN))
      d.tau <- Matrix::rowSums(x = W.tau)  # = d + tau
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
    
    # -------------------------------------------------------------------------
    # Y-X interaction and graph smoothing
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing density-weighted expression (M.local)...")
    }
    
    # Y-weighted expression: diag(Y) %*% X
    YX.interaction <-
      (Matrix::Diagonal(x = Y) %*% X) |>
      (\(x)
       `dimnames<-`(x = x,
                    value = dimnames(x = X))
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
                    value = list(rownames(x = X),
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
    # Compute loadings via regression
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Computing gene/marker loadings...")
    }
    
    loadings <-
      seq_len(length.out = ncol(x = scores)) |>
      stats::setNames(nm = colnames(x = scores)) |>
      lapply(FUN = function(k) {
        
        score.k <-
          scores[, k] - mean(x = scores[, k])
        
        denom <-
          Matrix::crossprod(x = score.k) |>
          as.numeric()
        
        if (denom < .Machine$double.eps) return(rep(x = NA_real_, times = ncol(X)))
        
        beta <-
          (Matrix::crossprod(x = score.k, 
                             y = X) |> 
             as.numeric()) / 
          denom
        
        return(beta)
        
      }) |>
      do.call(what = cbind)
    
    rownames(x = loadings) <-
      colnames(x = X)
    
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
    
    # Normalized Laplacian: L_sym = I - D^{-1/2} W D^{-1/2}
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
          lazy.alpha = .lazy.alpha
        )
      )
    
    if (isTRUE(x = .store.M)) {
      .tdr.obj$specDE[[.coef.col]]$M.local <-
        M.local
    }
    
    if (isTRUE(x = .verbose)) {
      message("\nResults stored in: .tdr.obj$specDE$", .coef.col)
      message("  $scores       : ", nrow(x = scores), " landmarks x ", ncol(x = scores), " components")
      message("  $loadings     : ", nrow(x = loadings), " features x ", ncol(x = loadings), " components")
      message("  $var.explained: Vk (proportion of variance)")
      message("  $Y.alignment  : Ak (|cor(Y, score)|)")
      message("  $smoothness   : Sk (Laplacian smoothness)")
    }
    
    return(.tdr.obj)
    
  }

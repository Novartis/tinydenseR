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
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object processed through \code{get.map()}.
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
#' @param ... Additional arguments passed to methods.
#' @export
#'
get.lm <- function(x, ...) UseMethod("get.lm")

#' @rdname get.lm
#' @export
get.lm.TDRObj <-
  function(
    x,
    .design,
    .contrasts = NULL,
    .block = NULL,
    .model.name = "default",
    .force.recalc = FALSE,
    .verbose = TRUE,
    ...){
    .tdr.obj <- x
    
    # -------------------------------------------------------------------------
    # Validate on-disk cache (if active) before proceeding
    # -------------------------------------------------------------------------
    
    .tdr_cache_validate_quiet(.tdr.obj)
    
    # -------------------------------------------------------------------------
    # Check if slot already exists
    # -------------------------------------------------------------------------
    
    if(!is.null(x = .tdr.obj@results$lm[[.model.name]]) && !isTRUE(x = .force.recalc)){
      stop(paste0("Model '", .model.name, "' already exists in .tdr.obj$map$lm. ",
                  "Use a different .model.name or set .force.recalc = TRUE to overwrite."))
    }
    
    if(nrow(x = .design) != length(x = .tdr.obj@cells)){
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
    nlib <- length(x = .tdr.obj@cells)
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
    
    if(is.null(x = .tdr.obj@density$norm)){
      stop("First run get.map")
    }
    
    # create list to hold results
    stats <-
      vector(mode = "list",
             length = 0)
    
    # Use pre-computed log.norm from get.map
    Y <-
      .tdr.obj@density$log.norm
    
    #if(nrow(x = Y) !=
    #   nrow(x = .tdr.obj@assay$expr)){
    #  
    #  stats$y <-
    #    stats$y[match(x = rownames(x = .tdr.obj@assay$expr),
    #                  table = rownames(x = stats$y)),]
    #  
    #  stats$y[
    #    is.na(x = stats$y)
    #  ] <- min(stats$y,
    #           na.rm = TRUE)
    #  
    #  rownames(x = stats$y) <-
    #    rownames(x = .tdr.obj@assay$expr)
    #  
    #}
    
    if(!(is.null(x = .block))){
      
      if(length(x = .block) != 1){
        stop("Block must be a vector of length 1")
      }
      if(!(.block %in% colnames(x = .tdr.obj@metadata))){
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
                                    block = .tdr.obj@metadata[[.block]])
    }
    
    if(isTRUE(x = .verbose)){
      message("\nfitting linear models")
    }
    
    stats$fit <-
      limma::lmFit(object = Y,
                   design = .design,
                   block = if(exists(x = "dupcor")) .tdr.obj@metadata[[.block]] else  NULL,
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
              w <- 1 / log10(x = Matrix::rowSums(x = .tdr.obj@density$norm) + 1)
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
        Matrix::t(x = .tdr.obj@landmark.embed$pca$coord)
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .tdr.obj@metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        # https://support.bioconductor.org/p/125489/#125602
        # duplicateCorrelation is more general and is THE ONLY SOLUTION when
        # you want to compare across blocking levels, e.g., comparing diseased
        # and healthy donors when each donor also contributes before/after treatment samples.
        dupcor.q <- 
          limma::duplicateCorrelation(object = tX,
                                      design = stats$fit$design[.tdr.obj@config$key,],
                                      block = .tdr.obj@metadata[[.block]][.tdr.obj@config$key])
      }
      
      fit.q <-
        limma::lmFit(object = tX,
                     design = stats$fit$design[.tdr.obj@config$key,],
                     block = if(exists(x = "dupcor.q")) .tdr.obj@metadata[[.block]][.tdr.obj@config$key] else  NULL,
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
         x[,1:elbow.sec.deriv(x = .tdr.obj@landmark.embed$pca$sdev^2,
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
    
    # Number of unique clusters/celltypes — skip composition limma fit when
    # only 1 level exists (nothing to compare, and limma needs >= 2 rows)
    n.clusters <-
      length(x = unique(x = .tdr.obj@landmark.annot$clustering$ids))
    
    n.celltypes <-
      if (!is.null(x = .tdr.obj@landmark.annot$celltyping$ids)) {
        length(x = unique(x = .tdr.obj@landmark.annot$celltyping$ids))
      } else {
        0L
      }

    if (n.clusters <= 1L) {
      warning("Only ", n.clusters, " cluster detected. ",
              "Skipping limma fit on cluster composition (nothing to compare). ",
              "Density-level analysis (get.plsD, get.dea, etc.) is unaffected.")
    }

    if (n.clusters > 1L) {

    if(!(is.null(x = .block))){
      
      if(length(x = .block) != 1){
        stop("Block must be a vector of length 1")
      }
      if(!(.block %in% colnames(x = .tdr.obj@metadata))){
        stop(paste0(.block,
                    " not found in metadata"))
      }
      
      if(isTRUE(x = .verbose)){
        message("\nestimating the intra-block correlation for stats from clustering")
      }
      
      cl.dupcor <- 
        log2(x = .tdr.obj@density$composition$clustering$cell.perc + 0.5) |>
        Matrix::t() |>
        limma::duplicateCorrelation(design = .design,
                                    block = .tdr.obj@metadata[[.block]])
    }
    
    if(isTRUE(x = .verbose)){
      message("\nfitting linear models")
    }
    
    stats$trad$clustering$fit <-
      limma::lmFit(object = log2(x = .tdr.obj@density$composition$clustering$cell.perc + 0.5) |>
                     Matrix::t(),
                   design = .design,
                   block = if(exists(x = "cl.dupcor")) .tdr.obj@metadata[[.block]] else NULL,
                   correlation = if(exists(x = "cl.dupcor")) cl.dupcor$consensus else NULL)
    
    if(!is.null(x = .contrasts)){
      
      stats$trad$clustering$fit <-
        limma::contrasts.fit(fit = stats$trad$clustering$fit,
                             contrasts = .contrasts)
      
    }
    
    stats$trad$clustering$fit <-
      limma::eBayes(fit = stats$trad$clustering$fit, robust = TRUE)
    
    stats$trad$clustering$fit$adj.p <-
      apply(X = stats$trad$clustering$fit$p.value,
            MARGIN = 2,
            FUN = stats::p.adjust,
            method = "fdr")
    
    } # end if (n.clusters > 1L)

    if(!is.null(x = .tdr.obj@landmark.annot$celltyping$ids) &&
       n.celltypes > 1L){
      
      if(!(is.null(x = .block))){
        
        if(length(x = .block) != 1){
          stop("Block must be a vector of length 1")
        }
        if(!(.block %in% colnames(x = .tdr.obj@metadata))){
          stop(paste0(.block,
                      " not found in metadata"))
        }
        
        if(isTRUE(x = .verbose)){
          message("\nestimating the intra-block correlation for stats from celltyping")
        }
        
        ct.dupcor <- 
          log2(x = .tdr.obj@density$composition$celltyping$cell.perc + 0.5) |>
          Matrix::t() |>
          limma::duplicateCorrelation(design = .design,
                                      block = .tdr.obj@metadata[[.block]])
      }
      
      if(isTRUE(x = .verbose)){
        message("\nfitting linear models")
      }
      
      stats$trad$celltyping$fit <-
        limma::lmFit(object = log2(x = .tdr.obj@density$composition$celltyping$cell.perc + 0.5) |>
                       Matrix::t(),
                     design = .design,
                     block = if(exists(x = "ct.dupcor")) .tdr.obj@metadata[[.block]] else NULL,
                     correlation = if(exists(x = "ct.dupcor")) ct.dupcor$consensus else NULL)
      
      if(!is.null(x = .contrasts)){
        
        stats$trad$celltyping$fit <-
          limma::contrasts.fit(fit = stats$trad$celltyping$fit,
                               contrasts = .contrasts)
        
      }
      
      stats$trad$celltyping$fit <-
        limma::eBayes(fit = stats$trad$celltyping$fit, robust = TRUE)
      
      stats$trad$celltyping$fit$adj.p <-
        apply(X = stats$trad$celltyping$fit$p.value,
              MARGIN = 2,
              FUN = stats::p.adjust,
              method = "fdr")
      
    }

    if (!is.null(x = .tdr.obj@landmark.annot$celltyping$ids) &&
        n.celltypes <= 1L) {
      warning("Only ", n.celltypes, " cell type detected. ",
              "Skipping limma fit on celltype composition (nothing to compare). ",
              "Density-level analysis (get.plsD, get.dea, etc.) is unaffected.")
    }
    
    # -------------------------------------------------------------------------
    # Store results and return modified .tdr.obj
    # -------------------------------------------------------------------------
    
    # Initialize stats slot if needed
    if(is.null(x = .tdr.obj@results$lm)){
      .tdr.obj@results$lm <- list()
    }
    
    # Store results under the model name
    .tdr.obj@results$lm[[.model.name]] <- stats
    
    if(isTRUE(x = .verbose)){
      message("\nResults stored in: .tdr.obj$map$lm$", .model.name)
    }
    
    return(.tdr.obj)
    
  }

# =============================================================================
# Internal helpers for pseudobulk aggregation (shared by design and marker modes)
# =============================================================================

#' Resolve cell indices from labels or raw indices
#'
#' Given either label-based (\code{.id}) or index-based (\code{.id.idx}) cell
#' specification, returns a named list of per-sample integer vectors of cell
#' indices.
#'
#' @param .tdr.obj A TDRObj.
#' @param .id Character vector of cluster/celltype IDs, or NULL.
#' @param .id.idx Integer vector of landmark indices, or NULL.
#' @param .id.from "clustering" or "celltyping".
#' @param .label.confidence Numeric (0-1) for fuzzy confidence thresholding.
#' @return Named list (one element per sample) of integer vectors.
#' @keywords internal
#' @noRd
.tdr_resolve_cell_idx <- function(.tdr.obj,
                                   .id = NULL,
                                   .id.idx = NULL,
                                   .id.from = "clustering",
                                   .label.confidence = 0.5) {

  if (!is.null(x = .id.idx)) {

    n.landmarks <- nrow(x = .tdr.obj@assay$expr)
    if (!is.null(x = n.landmarks) &&
        !all(.id.idx %in% seq_len(length.out = n.landmarks))) {
      stop(paste0(".id.idx must be an integer vector between 1 and ", n.landmarks))
    }

    tmp.lbl <- rep(x = "out", times = nrow(x = .tdr.obj@assay$expr))
    tmp.lbl[.id.idx] <- "in"

    lapply(
      X = .tdr_get_map_slot_all(.tdr.obj, "fuzzy.graphs"),
      FUN = function(smpl) {
        in.and.out <-
          .tdr_transfer_labels(
            .method = "fuzzy",
            .label.confidence = .label.confidence,
            .n.cells = nrow(x = smpl),
            .cell.names = seq_len(length.out = nrow(x = smpl)),
            .fgraph = smpl,
            .landmark.labels = tmp.lbl
          )
        which(x = in.and.out == "in")
      }
    )

  } else if (!is.null(x = .id)) {

    .id.from <- match.arg(arg = .id.from, choices = c("clustering", "celltyping"))

    if (!all(.id %in% unique(x = .tdr.obj@landmark.annot[[.id.from]]$ids))) {
      stop(paste0(
        paste0(.id[!(.id %in% unique(x = .tdr.obj@landmark.annot[[.id.from]]$ids))],
               collapse = ", "),
        " not found in ", .id.from
      ))
    }

    lapply(
      X = .tdr_get_map_slot_all(.tdr.obj, .id.from),
      FUN = function(smpl) {
        which(x = smpl %in% .id)
      }
    )

  } else {

    # All cells per sample
    stats::setNames(
      object = .tdr.obj@metadata$n.cells,
      nm = names(x = .tdr.obj@cells)
    ) |> lapply(FUN = seq_len)

  }
}

#' Pseudobulk aggregation for RNA data
#'
#' Computes fuzzy-weighted pseudobulk expression for each sample (RNA).
#'
#' @param .source Raw data object (or NULL for file backend).
#' @param .tdr.obj A TDRObj.
#' @param .samples Character vector of sample names.
#' @param .cell.idx.list Named list of per-sample cell index vectors.
#' @param .robust If TRUE, wrap per-sample aggregation in tryCatch (return
#'   zero vector on error).
#' @return A genes x samples matrix of pseudobulk counts.
#' @keywords internal
#' @noRd
.tdr_pseudobulk_aggregate_rna <- function(.source, .tdr.obj,
                                           .samples, .cell.idx.list,
                                           .robust = FALSE) {
  stats::setNames(object = .samples, nm = .samples) |>
    lapply(FUN = function(smpl) {
      exprs.mat <-
        .get_sample_matrix(.source, .tdr.obj, smpl) |>
        methods::as(Class = "dgCMatrix")

      .agg_fn <- function() {
        cell.idx <- .cell.idx.list[[smpl]]
        wcl <-
          .tdr_get_map_slot(.tdr.obj, "fuzzy.graphs", smpl)[cell.idx, , drop = FALSE]
        wsum <- exprs.mat[, cell.idx] %*% wcl
        (Matrix::rowSums(x = wsum) / sum(wcl)) *
          (sum(Matrix::rowSums(x = wcl) > 0))
      }

      if (isTRUE(x = .robust)) {
        tryCatch(
          expr = .agg_fn(),
          error = function(e) {
            stats::setNames(
              object = rep(x = 0, times = nrow(x = exprs.mat)),
              nm = rownames(x = exprs.mat)
            )
          }
        )
      } else {
        .agg_fn()
      }
    }) |>
    do.call(what = cbind)
}

#' Pseudobulk aggregation for cytometry data
#'
#' Computes fuzzy-weighted pseudobulk expression for each sample (cytometry).
#'
#' @param .source Raw data object (or NULL for file backend).
#' @param .tdr.obj A TDRObj.
#' @param .samples Character vector of sample names.
#' @param .cell.idx.list Named list of per-sample cell index vectors.
#' @param .robust If TRUE, wrap per-sample aggregation in tryCatch (return
#'   zero vector on error).
#' @return A samples x features matrix of pseudobulk expression.
#' @keywords internal
#' @noRd
.tdr_pseudobulk_aggregate_cyto <- function(.source, .tdr.obj,
                                            .samples, .cell.idx.list,
                                            .robust = FALSE) {
  cols <- colnames(x = .tdr.obj@assay$expr)

  stats::setNames(object = .samples, nm = .samples) |>
    lapply(FUN = function(smpl) {
      exprs.mat <- .get_sample_matrix(.source, .tdr.obj, smpl)

      .agg_fn <- function() {
        cell.idx <- .cell.idx.list[[smpl]]
        wcl <-
          .tdr_get_map_slot(.tdr.obj, "fuzzy.graphs", smpl)[cell.idx, , drop = FALSE]
        wsum <- Matrix::t(x = exprs.mat[cell.idx, cols]) %*% wcl
        Matrix::rowSums(x = wsum) / sum(wcl)
      }

      if (isTRUE(x = .robust)) {
        tryCatch(
          expr = .agg_fn(),
          error = function(e) {
            stats::setNames(
              object = rep(x = 0, times = length(x = cols)),
              nm = cols
            )
          }
        )
      } else {
        .agg_fn()
      }
    }) |>
    do.call(what = rbind)
}

#' Pseudobulk Differential Expression Analysis
#'
#' Performs pseudobulk differential expression (DE) analysis for genes/markers with two modes.
#'
#' \strong{Design mode} (\code{.mode = "design"}): Aggregates cells into pseudobulk samples
#' using fuzzy landmark membership, then uses limma-voom (RNA) or limma (cytometry) to test for
#' DE across experimental conditions. Uses a user-supplied design matrix with samples as replicates.
#'
#' \strong{Marker mode} (\code{.mode = "marker"}): Identifies marker genes/proteins distinguishing
#' one cell population from another (or all others) via a within-sample paired comparison. Unlike
#' design mode which tests experimental conditions, marker mode compares cell populations to find
#' defining features.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object processed through \code{get.map()}.
#' @param .source The raw data object for non-file backends. \code{NULL} (default) for
#'   the files backend; otherwise a Seurat, SingleCellExperiment, or anndataR AnnData object.
#'   Used by \code{.get_sample_matrix()} to retrieve per-sample expression matrices.
#' @param .mode Character: analysis mode. One of \code{"design"} or \code{"marker"}.
#'   If \code{NULL} (default), auto-detected from arguments:
#'   \itemize{
#'     \item \code{.design} provided \eqn{\Rightarrow} \code{"design"}
#'     \item \code{.id}/\code{.id.idx} provided without \code{.design} \eqn{\Rightarrow} \code{"marker"}
#'   }
#' @param .design Design matrix specifying experimental design (design mode only). Rows = samples,
#'   columns = coefficients.
#' @param .contrasts Optional contrast matrix for specific comparisons (design mode only). Create
#'   with \code{limma::makeContrasts()}. If NULL, tests all \code{.design} coefficients.
#' @param .block Optional character: column name in \code{.tdr.obj$metadata} for blocking factor
#'   (design mode only, e.g., "Donor"). Accounts for within-block correlation.
#' @param .geneset.ls Optional named list of character vectors defining gene sets for GSVA enrichment
#'   analysis. Only for RNA data. Example: \code{list("Tcell" = c("CD3D", "CD3E"), "Bcell" = c("CD19", "MS4A1"))}.
#' @param .id Optional character vector of cluster/celltype IDs. In design mode, restricts
#'   analysis to cells matching these IDs. In marker mode, defines group 1 (test group).
#' @param .id.idx Optional integer vector specifying landmark indices. In design mode, restricts
#'   analysis to cells confidently assigned to these landmarks. In marker mode, defines group 1
#'   landmark indices.
#' @param .id2 Character vector of cluster/celltype IDs for group 2 (marker mode only). Default
#'   \code{"..all.other.landmarks.."} compares group 1 to all other cells. Can specify specific IDs
#'   for pairwise comparisons.
#' @param .id2.idx Optional integer vector specifying landmark indices for group 2 (marker mode
#'   only). When provided, takes priority over \code{.id2}.
#' @param .id.from Character: \code{"clustering"} or \code{"celltyping"}. Source of IDs in
#'   \code{.id} and \code{.id2}. Default \code{NULL} (resolved to \code{"clustering"} when needed).
#' @param .model.name Character string naming this model fit (default \code{"default"}).
#' @param .result.name Character string naming this result. In design mode defaults to \code{"all"};
#'   in marker mode auto-generated from \code{.id} and \code{.id2}. Used as the storage key:
#'   \code{$pbDE[[.model.name]][[.result.name]]} (design) or
#'   \code{$markerDE[[.model.name]][[.result.name]]} (marker).
#' @param .population.name \code{NULL}. Deprecated alias for \code{.result.name} (design mode).
#' @param .comparison.name \code{NULL}. Deprecated alias for \code{.result.name} (marker mode).
#' @param .force.recalc Logical: if TRUE, overwrite existing results in the specified slot
#'   (default FALSE).
#' @param .verbose Logical: print progress messages? Default TRUE.
#' @param .label.confidence Numeric scalar in \code{[0,1]} controlling the minimum posterior
#'   confidence required to assign a cell to a set of target landmarks (used when \code{.id.idx} or
#'   \code{.id2.idx} is provided). Default 0.5.
#'
#' @return The modified \code{.tdr.obj} with results stored depending on mode:
#'   \subsection{Design mode}{Results in \code{.tdr.obj$pbDE[[.model.name]][[.result.name]]}:
#'     \describe{
#'       \item{\code{coefficients}}{Log fold change matrix (features x coefficients)}
#'       \item{\code{p.value}}{P-values (features x coefficients)}
#'       \item{\code{adj.p}}{FDR-adjusted p-values (features x coefficients)}
#'       \item{\code{smpl.outlier}}{Logical vector indicating outlier samples}
#'       \item{\code{id.idx}}{Per-sample list of cell indices used}
#'       \item{\code{n.pseudo}}{Integer vector of pseudobulk cell counts per sample}
#'       \item{\code{geneset}}{(RNA + \code{.geneset.ls}) GSVA results}
#'     }
#'   }
#'   \subsection{Marker mode}{Results in \code{.tdr.obj$markerDE[[.model.name]][[.result.name]]}:
#'     \describe{
#'       \item{\code{coefficients}}{Log fold changes (features x coefficients)}
#'       \item{\code{p.value}}{P-values (features x coefficients)}
#'       \item{\code{adj.p}}{FDR-adjusted p-values (features x coefficients)}
#'       \item{\code{smpl.outlier.1}}{Logical: samples excluded from group 1}
#'       \item{\code{smpl.outlier.2}}{Logical: samples excluded from group 2}
#'       \item{\code{id1.idx}}{Per-sample cell indices for group 1}
#'       \item{\code{id2.idx}}{Per-sample cell indices for group 2}
#'       \item{\code{n.pseudo1}}{Pseudobulk cell counts per sample for group 1}
#'       \item{\code{n.pseudo2}}{Pseudobulk cell counts per sample for group 2}
#'       \item{\code{geneset}}{(RNA + \code{.geneset.ls}) GSVA results}
#'     }
#'   }
#'
#' @details
#' \subsection{Design mode}{
#' Tests for DE across experimental conditions using a user-supplied design matrix:
#' \enumerate{
#'   \item \strong{Cell selection}: If \code{.id} specified, select matching cells. If
#'     \code{.id.idx} specified, use fuzzy confidence thresholding. Otherwise use all cells.
#'   \item \strong{Pseudobulk aggregation}: Fuzzy-weighted expression per sample.
#'   \item \strong{Outlier removal}: (RNA only) Exclude samples with <10\% of average cell count.
#'   \item \strong{Normalization}: TMM + voom (RNA) or as-is (cytometry).
#'   \item \strong{Linear modeling}: \code{limma::lmFit} with optional blocking.
#'   \item \strong{Empirical Bayes}: \code{limma::eBayes(robust = TRUE)}.
#'   \item \strong{GSVA} (optional, RNA only): Gene set variation analysis.
#' }
#' }
#'
#' \subsection{Marker mode}{
#' Compares cell populations via a within-sample paired design (\code{~ .ids + .pairs}):
#' \enumerate{
#'   \item \strong{Cell selection}: Extract cells for group 1 (\code{.id}) and group 2 (\code{.id2}).
#'   \item \strong{Pseudobulk aggregation}: Independent aggregation per group per sample.
#'   \item \strong{Outlier removal}: (RNA only) Cross-group outlier flagging.
#'   \item \strong{Paired comparison}: \code{limma::lmFit} with sample-as-pair blocking.
#' }
#' Positive logFC = higher in group 1, negative = higher in group 2.
#' }
#'
#' @seealso \code{\link{get.map}} (required), \code{\link{plotPbDE}}, \code{\link{plotMarkerDE}},
#'   \code{\link{get.plsD}}
#'
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |> get.map()
#'
#' # --- Design mode (default when .design is provided) ---
#' design <- model.matrix(~ Condition, data = .meta)
#' lm.cells <- get.pbDE(lm.cells, .design = design)
#' plotPbDE(lm.cells, .coefs = "ConditionB")
#'
#' # DE within specific cell type
#' lm.cells <- get.pbDE(lm.cells, .design = design,
#'                      .id = c("1", "2", "3"),
#'                      .id.from = "clustering",
#'                      .result.name = "tcells")
#'
#' # --- Marker mode (auto-detected when .id is provided without .design) ---
#' lm.cells <- get.pbDE(lm.cells, .mode = "marker",
#'                      .id = "cluster.3",
#'                      .result.name = "cluster3_markers")
#' plotMarkerDE(lm.cells, .comparison.name = "cluster3_markers")
#'
#' # Pairwise comparison
#' lm.cells <- get.pbDE(lm.cells, .mode = "marker",
#'                      .id = c("cluster.1", "cluster.3"),
#'                      .id2 = c("cluster.2", "cluster.4"),
#'                      .result.name = "cd4_vs_cd8")
#'
#' # Access results
#' lm.cells$pbDE$default$all$adj.p
#' lm.cells$markerDE$default$cluster3_markers$coefficients
#' }
#'
#' @param ... Additional arguments passed to methods.
#' @export
#'
get.pbDE <- function(x, ...) UseMethod("get.pbDE")

#' @rdname get.pbDE
#' @export
get.pbDE.TDRObj <-
  function(
    x,
    .source = NULL,
    .mode = NULL,
    .design = NULL,
    .contrasts = NULL,
    .block = NULL,
    .geneset.ls = NULL,
    .id = NULL,
    .id.idx = NULL,
    .id2 = "..all.other.landmarks..",
    .id2.idx = NULL,
    .id.from = NULL,
    .model.name = "default",
    .result.name = NULL,
    .population.name = NULL,
    .comparison.name = NULL,
    .force.recalc = FALSE,
    .verbose = TRUE,
    .label.confidence = 0.5,
    ...
  ) {
    .tdr.obj <- x

    .tdr_validate_label_confidence(.label.confidence)

    # R CMD check appeasement
    i <- j <- landmark <- cell <- label <- x <- confidence <- NULL

    # =========================================================================
    # Resolve deprecated aliases for .result.name
    # =========================================================================

    if (!is.null(x = .population.name)) {
      warning("`.population.name` is deprecated. Use `.result.name` instead.",
              call. = FALSE)
      if (is.null(x = .result.name)) .result.name <- .population.name
    }
    if (!is.null(x = .comparison.name)) {
      warning("`.comparison.name` is deprecated. Use `.result.name` instead.",
              call. = FALSE)
      if (is.null(x = .result.name)) .result.name <- .comparison.name
    }

    # =========================================================================
    # Mode detection & validation
    # =========================================================================

    .has.design <- !is.null(x = .design)
    .has.id     <- !is.null(x = .id) || !is.null(x = .id.idx)
    .has.id2    <- !is.null(x = .id2.idx) ||
      !identical(x = .id2, y = "..all.other.landmarks..")

    if (is.null(x = .mode)) {
      if (.has.design && !.has.id2 && is.null(x = .id2.idx)) {
        .mode <- "design"
      } else if (.has.design && (.has.id2 || !is.null(x = .id2.idx))) {
        stop(
          "Conflicting arguments: .design and .id2/.id2.idx cannot both be provided.\n",
          "  - For pseudobulk DE comparing samples according to a design matrix, ",
          "provide .design (and optionally .id/.id.idx to subset to a population).\n",
          "  - For marker identification comparing cell populations within samples, ",
          "provide .id/.id.idx (group 1) and .id2/.id2.idx (group 2) without .design.",
          call. = FALSE
        )
      } else if (.has.id) {
        .mode <- "marker"
      } else if (!.has.design && !.has.id) {
        stop(
          "Cannot determine analysis mode. Provide either:\n",
          "  - .design (for pseudobulk DE comparing samples across conditions), or\n",
          "  - .id/.id.idx (for marker identification comparing cell populations).",
          call. = FALSE
        )
      }
    }

    .mode <- match.arg(arg = .mode, choices = c("design", "marker"))

    if (.mode == "design" && !.has.design) {
      stop("In design mode, .design is required.", call. = FALSE)
    }
    if (.mode == "marker" && .has.design) {
      stop(
        ".design is not used in marker mode.\n",
        "  - For pseudobulk DE comparing samples, use .mode = 'design'.\n",
        "  - For marker identification, do not provide .design.",
        call. = FALSE
      )
    }
    if (.mode == "marker" && is.null(x = .id) && is.null(x = .id.idx)) {
      stop("In marker mode, .id or .id.idx (group 1) is required.", call. = FALSE)
    }
    if (.mode == "design" && !is.null(x = .id2.idx)) {
      warning(".id2.idx is ignored in design mode.", call. = FALSE)
    }
    if (.mode == "design" &&
        !identical(x = .id2, y = "..all.other.landmarks..")) {
      warning(".id2 is ignored in design mode.", call. = FALSE)
    }
    if (.mode == "marker" && !is.null(x = .contrasts)) {
      warning(".contrasts is not used in marker mode and will be ignored.",
              call. = FALSE)
    }
    if (.mode == "marker" && !is.null(x = .block)) {
      warning(".block is not used in marker mode and will be ignored.",
              call. = FALSE)
    }

    # Validate shared arguments
    if (!is.null(x = .geneset.ls)) {
      if (.tdr.obj@config$assay.type != "RNA") {
        stop(".geneset.ls is only supported for RNA assay type")
      } else if (!is.list(x = .geneset.ls)) {
        stop(".geneset.ls must be a list of character vectors")
      }
    }

    # Resolve .id.from default
    if (is.null(x = .id.from) && (!is.null(x = .id) || .mode == "marker")) {
      .id.from <- "clustering"
    }

    # Validate on-disk cache
    .tdr_cache_validate_quiet(.tdr.obj)

    # =========================================================================
    # DESIGN MODE
    # =========================================================================
    if (.mode == "design") {

      .design <- as.matrix(x = .design)

      if (nrow(x = .design) != length(x = .tdr.obj@cells)) {
        stop("Number of rows in design matrix must be equal to the number of samples")
      }

      # Check design rank (source: https://rdrr.io/bioc/edgeR/src/R/glmfit.R)
      ne <- limma::nonEstimable(x = .design)
      if (!is.null(x = ne)) {
        stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n",
                   paste(ne, collapse = " ")))
      }

      # Check intercept + continuous covariate
      no.intercept <-
        !(any(Matrix::colSums(x = .design == 1) == nrow(x = .design)))
      num.covariate <-
        !all(as.vector(x = .design) %in% c(0, 1))
      if (no.intercept & num.covariate) {
        warning("No intercept in the model but at least one of the covariates is a continuous variable. This model is not the same as the one with an intercept because it assumes that the continuous variable is centered at 0.")
      }

      # --- Auto-set .result.name ---
      if (is.null(x = .result.name)) {
        if (!is.null(x = .id)) {
          .result.name <- paste(.id, collapse = "_")
          if (isTRUE(x = .verbose)) {
            message(sprintf("Auto-setting .result.name to '%s'", .result.name))
          }
        } else if (!is.null(x = .id.idx)) {
          .result.name <- paste0("idx_", length(x = .id.idx))
          if (isTRUE(x = .verbose)) {
            message(sprintf("Auto-setting .result.name to '%s' (custom landmark indices)",
                            .result.name))
          }
        } else {
          .result.name <- "all"
        }
      }

      # Check if slot already exists
      if (!is.null(x = .tdr.obj@results$pb[[.model.name]][[.result.name]]) &&
          !isTRUE(x = .force.recalc)) {
        stop(paste0("Results for model '", .model.name, "' and result '", .result.name,
                    "' already exist in .tdr.obj$pbDE. ",
                    "Use different names or set .force.recalc = TRUE to overwrite."))
      }

      # --- Resolve cell indices ---
      .lm.idx <- .tdr_resolve_cell_idx(
        .tdr.obj, .id = .id, .id.idx = .id.idx,
        .id.from = if (!is.null(x = .id.from)) .id.from else "clustering",
        .label.confidence = .label.confidence
      )

      if (isTRUE(x = .verbose)) {
        message(sprintf("\nUsing %s cells for pseudobulk aggregation",
                        paste(lengths(x = .lm.idx), collapse = ", ")))
      }

      n.pseudo <- lengths(x = .lm.idx)
      smpl.outlier <- (n.pseudo / mean(x = n.pseudo)) < 0.1

      # --- RNA pseudobulk aggregation ---
      if (.tdr.obj@config$assay.type == "RNA") {

        if (any(smpl.outlier)) {
          warning(paste0(
            "The following samples were removed since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
            paste(names(x = smpl.outlier)[smpl.outlier], collapse = "\n")
          ))
        }

        counts <- .tdr_pseudobulk_aggregate_rna(
          .source, .tdr.obj,
          .samples = names(x = .tdr.obj@cells)[!smpl.outlier],
          .cell.idx.list = .lm.idx,
          .robust = FALSE
        )

        dge <- edgeR::DGEList(counts = counts)

        tmp.design <-
          .design[!smpl.outlier, ] |>
          (\(x) x[, Matrix::colSums(x = x == 0) != nrow(x = x), drop = FALSE])()

        keep <- edgeR::filterByExpr(y = dge, design = tmp.design)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        dge <- edgeR::calcNormFactors(object = dge, method = "TMM")
        v <- limma::voom(counts = dge, design = tmp.design, plot = FALSE)

        if (!is.null(x = .block)) {
          if (length(x = .block) != 1) stop("Block must be a vector of length 1")
          if (!(.block %in% colnames(x = .tdr.obj@metadata))) {
            stop(paste0(.block, " not found in metadata"))
          }
          if (isTRUE(x = .verbose)) message("\nestimating the intra-block correlation")
          dupcor <-
            limma::duplicateCorrelation(
              object = v, design = tmp.design,
              block = .tdr.obj@metadata[[.block]][!smpl.outlier]
            )
        }

        if (isTRUE(x = .verbose)) message("\nfitting linear models")

        .de <-
          limma::lmFit(
            object = v, design = tmp.design,
            block = if (exists(x = "dupcor")) .tdr.obj@metadata[[.block]][!smpl.outlier] else NULL,
            correlation = if (exists(x = "dupcor")) dupcor$consensus else NULL
          )

        if (!is.null(x = .contrasts)) {
          tmp.contrasts <-
            .contrasts[colnames(x = tmp.design), , drop = FALSE] |>
            (\(x) x[, Matrix::colSums(x = x == 0) != nrow(x = x), drop = FALSE])()
          .de <- limma::contrasts.fit(fit = .de, contrasts = tmp.contrasts)
        }

        .de <- limma::eBayes(fit = .de, robust = TRUE)
        .de$adj.p <- apply(X = .de$p.value, MARGIN = 2,
                           FUN = stats::p.adjust, method = "fdr")

        # GSVA (optional)
        if (!is.null(x = .geneset.ls)) {
          gsva.es <-
            GSVA::gsvaParam(exprData = v$E, geneSets = .geneset.ls,
                            kcdf = "Gaussian", maxDiff = TRUE) |>
            GSVA::gsva()

          if (!is.null(x = .block)) {
            if (isTRUE(x = .verbose)) message("\nestimating the intra-block correlation")
            gsva.dupcor <-
              limma::duplicateCorrelation(
                object = gsva.es, design = tmp.design,
                block = .tdr.obj@metadata[[.block]][!smpl.outlier]
              )
          }

          gsva.fit <-
            limma::lmFit(
              object = gsva.es, design = tmp.design,
              block = if (exists(x = "gsva.dupcor")) .tdr.obj@metadata[[.block]][!smpl.outlier] else NULL,
              correlation = if (exists(x = "gsva.dupcor")) gsva.dupcor$consensus else NULL
            )

          if (!is.null(x = .contrasts)) {
            gsva.fit <- limma::contrasts.fit(fit = gsva.fit, contrasts = tmp.contrasts)
          }

          gsva.fit <- limma::eBayes(fit = gsva.fit, robust = TRUE)
          gsva.fit$adj.p <- apply(X = gsva.fit$p.value, MARGIN = 2,
                                  FUN = stats::p.adjust, method = "fdr")
          .de$geneset <- gsva.fit
          .de$geneset$E <- gsva.es
        }

      } else {
        # --- Cytometry pseudobulk aggregation ---

        counts <- .tdr_pseudobulk_aggregate_cyto(
          .source, .tdr.obj,
          .samples = names(x = .tdr.obj@cells),
          .cell.idx.list = .lm.idx,
          .robust = FALSE
        )

        tmp.design <-
          .design |>
          (\(x) x[, Matrix::colSums(x = x == 0) != nrow(x = x), drop = FALSE])()

        if (!is.null(x = .block)) {
          if (length(x = .block) != 1) stop("Block must be a vector of length 1")
          if (!(.block %in% colnames(x = .tdr.obj@metadata))) {
            stop(paste0(.block, " not found in metadata"))
          }
          if (isTRUE(x = .verbose)) message("\nestimating the intra-block correlation")
          dupcor <-
            limma::duplicateCorrelation(
              object = Matrix::t(x = counts), design = tmp.design,
              block = .tdr.obj@metadata[[.block]]
            )
        }

        if (isTRUE(x = .verbose)) message("\nfitting linear models")

        .de <-
          limma::lmFit(
            object = Matrix::t(x = counts), design = tmp.design,
            block = if (exists(x = "dupcor")) .tdr.obj@metadata[[.block]] else NULL,
            correlation = if (exists(x = "dupcor")) dupcor$consensus else NULL
          )

        if (!is.null(x = .contrasts)) {
          tmp.contrasts <-
            .contrasts[colnames(x = tmp.design), , drop = FALSE] |>
            (\(x) x[, Matrix::colSums(x = x == 0) != nrow(x = x), drop = FALSE])()
          .de <- limma::contrasts.fit(fit = .de, contrasts = tmp.contrasts)
        }

        .de <- limma::eBayes(fit = .de, robust = TRUE)
        .de$adj.p <- apply(X = .de$p.value, MARGIN = 2,
                           FUN = stats::p.adjust, method = "fdr")
      }

      .de$smpl.outlier <- smpl.outlier
      .de$id.idx <- .lm.idx
      .de$n.pseudo <- n.pseudo

      # Store results in @results$pb
      if (is.null(x = .tdr.obj@results$pb)) .tdr.obj@results$pb <- list()
      if (is.null(x = .tdr.obj@results$pb[[.model.name]])) {
        .tdr.obj@results$pb[[.model.name]] <- list()
      }
      .tdr.obj@results$pb[[.model.name]][[.result.name]] <- .de

      if (isTRUE(x = .verbose)) {
        message(sprintf("\nResults stored in .tdr.obj$pbDE$%s$%s",
                        .model.name, .result.name))
      }

      return(.tdr.obj)
    }

    # =========================================================================
    # MARKER MODE
    # =========================================================================
    if (.mode == "marker") {

      if (is.null(x = .id.from)) .id.from <- "clustering"
      .id.from <- match.arg(arg = .id.from, choices = c("clustering", "celltyping"))

      # --- Resolve group 1 cell indices ---
      .lm1.idx <- .tdr_resolve_cell_idx(
        .tdr.obj, .id = .id, .id.idx = .id.idx,
        .id.from = .id.from,
        .label.confidence = .label.confidence
      )

      # --- Resolve group 2 cell indices ---
      if (!is.null(x = .id2.idx)) {

        # .id2.idx takes priority over .id2
        .lm2.idx <- .tdr_resolve_cell_idx(
          .tdr.obj, .id = NULL, .id.idx = .id2.idx,
          .id.from = .id.from,
          .label.confidence = .label.confidence
        )

      } else if ("..all.other.landmarks.." %in% .id2) {

        if (isTRUE(x = .verbose)) message("using `..all.other.landmarks..` for .id2")

        .id2 <- "..all.other.landmarks.."

        # Complement: all cells NOT in group 1
        .lm2.idx <-
          seq_along(along.with = .lm1.idx) |>
          stats::setNames(nm = names(x = .lm1.idx)) |>
          lapply(FUN = function(smpl) {
            seq_len(
              length.out = .tdr.obj@config$sampling$n.cells[smpl]
            )[-.lm1.idx[[smpl]]]
          })

      } else {

        # Validate .id2 labels
        if (!all(.id2 %in% c(
          unique(x = .tdr.obj@landmark.annot[[.id.from]]$ids) |> as.character()
        ))) {
          stop(paste0(
            paste0(.id2[!(.id2 %in% (
              unique(x = .tdr.obj@landmark.annot[[.id.from]]$ids) |>
                as.character()
            ))], collapse = ", "),
            " not found in ", .id.from
          ))
        }

        .lm2.idx <-
          lapply(
            X = .tdr_get_map_slot_all(.tdr.obj, .id.from),
            FUN = function(smpl) which(x = smpl %in% .id2)
          )
      }

      # --- Auto-generate .result.name ---
      if (is.null(x = .result.name)) {
        id1.str <- if (!is.null(x = .id)) paste0(.id, collapse = "_") else
          paste0("idx", length(x = .id.idx))
        if (identical(x = .id2, y = "..all.other.landmarks..")) {
          id2.str <- "all"
        } else if (!is.null(x = .id2.idx)) {
          id2.str <- paste0("idx", length(x = .id2.idx))
        } else {
          id2.str <- paste0(.id2, collapse = "_")
        }
        .result.name <- paste0(id1.str, "_vs_", id2.str)
      }

      # Initialize storage
      if (is.null(x = .tdr.obj@results$marker)) .tdr.obj@results$marker <- list()
      if (is.null(x = .tdr.obj@results$marker[[.model.name]])) {
        .tdr.obj@results$marker[[.model.name]] <- list()
      }

      # Check for existing results
      if (!isTRUE(x = .force.recalc) &&
          !is.null(x = .tdr.obj@results$marker[[.model.name]][[.result.name]])) {
        message(sprintf(
          "Results already exist at .tdr.obj$markerDE$%s$%s. Use .force.recalc = TRUE to recalculate.",
          .model.name, .result.name
        ))
        return(.tdr.obj)
      }

      # --- Pseudobulk cell counts and outlier detection ---
      n.pseudo1 <- lengths(x = .lm1.idx)
      n.pseudo2 <- lengths(x = .lm2.idx)

      smpl.outlier.1 <-
        (n.pseudo1 / mean(x = c(n.pseudo1, n.pseudo2))) < 0.1
      smpl.outlier.2 <-
        (n.pseudo2 / mean(x = c(n.pseudo1, n.pseudo2))) < 0.1

      # --- RNA marker aggregation ---
      if (.tdr.obj@config$assay.type == "RNA") {

        if (any(smpl.outlier.1)) {
          warning(paste0(
            "The following samples were removed from group 1 since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
            paste(names(x = smpl.outlier.1)[smpl.outlier.1], collapse = "\n")
          ))
        }
        if (any(smpl.outlier.2)) {
          warning(paste0(
            "The following samples were removed from group 2 since they have less than a tenth of the average number of cells in pseudobulks, which can lead to misleading results:\n",
            paste(names(x = smpl.outlier.2)[smpl.outlier.2], collapse = "\n")
          ))
        }

        counts1 <- .tdr_pseudobulk_aggregate_rna(
          .source, .tdr.obj,
          .samples = names(x = .tdr.obj@cells)[!smpl.outlier.1],
          .cell.idx.list = .lm1.idx,
          .robust = TRUE
        )
        colnames(counts1) <- paste0("counts1.", colnames(x = counts1))

        counts2 <- .tdr_pseudobulk_aggregate_rna(
          .source, .tdr.obj,
          .samples = names(x = .tdr.obj@cells)[!smpl.outlier.2],
          .cell.idx.list = .lm2.idx,
          .robust = TRUE
        )
        colnames(counts2) <- paste0("counts2.", colnames(x = counts2))

        dge <- cbind(counts1, counts2) |> edgeR::DGEList()

        .ids <-
          colnames(x = dge) |>
          strsplit(split = ".", fixed = TRUE) |>
          lapply(FUN = "[", i = 1) |>
          unlist() |>
          gsub(pattern = "counts", replacement = ".id", fixed = TRUE) |>
          factor(levels = c(".id2", ".id1"))

        .pairs <-
          colnames(x = dge) |>
          gsub(pattern = "^counts1.|^counts2.", replacement = "", fixed = FALSE)

        tmp.design <-
          stats::model.matrix(object = ~ .ids + .pairs) |>
          (\(x) `colnames<-`(
            x = x,
            value = colnames(x = x) |>
              gsub(pattern = "^\\.ids|^\\.pairs", replacement = "", fixed = FALSE)
          ))() |>
          (\(x) x[, Matrix::colSums(x = x) > 1])()

        keep <- edgeR::filterByExpr(y = dge, design = tmp.design)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        dge <- edgeR::calcNormFactors(object = dge, method = "TMM")
        v <- limma::voom(counts = dge, design = tmp.design, plot = FALSE)

        .de <- limma::lmFit(object = v, design = tmp.design)
        .de <- limma::eBayes(fit = .de, robust = TRUE)
        .de$adj.p <- apply(X = .de$p.value, MARGIN = 2,
                           FUN = stats::p.adjust, method = "fdr")

        # GSVA (optional)
        if (!is.null(x = .geneset.ls)) {
          gsva.es <-
            GSVA::gsvaParam(exprData = v$E, geneSets = .geneset.ls,
                            kcdf = "Gaussian", maxDiff = TRUE) |>
            GSVA::gsva()

          gsva.fit <- limma::lmFit(object = gsva.es, design = tmp.design)
          gsva.fit <- limma::eBayes(fit = gsva.fit, robust = TRUE)
          gsva.fit$adj.p <- apply(X = gsva.fit$p.value, MARGIN = 2,
                                  FUN = stats::p.adjust, method = "fdr")
          .de$geneset <- gsva.fit
          .de$geneset$E <- gsva.es
        }

      } else {
        # --- Cytometry marker aggregation ---

        counts1 <- .tdr_pseudobulk_aggregate_cyto(
          .source, .tdr.obj,
          .samples = names(x = .tdr.obj@cells),
          .cell.idx.list = .lm1.idx,
          .robust = TRUE
        )
        rownames(counts1) <- paste0("counts1.", rownames(x = counts1))

        counts2 <- .tdr_pseudobulk_aggregate_cyto(
          .source, .tdr.obj,
          .samples = names(x = .tdr.obj@cells),
          .cell.idx.list = .lm2.idx,
          .robust = TRUE
        )
        rownames(counts2) <- paste0("counts2.", rownames(x = counts2))

        counts <- rbind(counts1, counts2) |> Matrix::t()

        .ids <-
          colnames(x = counts) |>
          strsplit(split = ".", fixed = TRUE) |>
          lapply(FUN = "[", i = 1) |>
          unlist() |>
          gsub(pattern = "counts", replacement = ".id", fixed = TRUE) |>
          factor(levels = c(".id2", ".id1"))

        .pairs <-
          colnames(x = counts) |>
          gsub(pattern = "^counts1.|^counts2.", replacement = "", fixed = FALSE)

        tmp.design <-
          stats::model.matrix(object = ~ .ids + .pairs) |>
          (\(x) `colnames<-`(
            x = x,
            value = colnames(x = x) |>
              gsub(pattern = "^\\.ids|^\\.pairs", replacement = "", fixed = FALSE)
          ))() |>
          (\(x) x[, Matrix::colSums(x = x) > 1])()

        .de <- limma::lmFit(object = counts, design = tmp.design)
        .de <- limma::eBayes(fit = .de, robust = TRUE)
        .de$adj.p <- apply(X = .de$p.value, MARGIN = 2,
                           FUN = stats::p.adjust, method = "fdr")
      }

      .de$smpl.outlier.1 <- smpl.outlier.1
      .de$smpl.outlier.2 <- smpl.outlier.2
      .de$id1.idx <- .lm1.idx
      .de$id2.idx <- .lm2.idx
      .de$n.pseudo1 <- n.pseudo1
      .de$n.pseudo2 <- n.pseudo2

      # Store results in @results$marker
      .tdr.obj@results$marker[[.model.name]][[.result.name]] <- .de

      if (isTRUE(x = .verbose)) {
        message(sprintf(
          "\nMarker DE results stored in .tdr.obj$markerDE$%s$%s",
          .model.name, .result.name
        ))
      }

      return(.tdr.obj)
    }
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
    .label.confidence = 0.5
) {
  .tdr_validate_label_confidence(.label.confidence)

  .Deprecated("get.pbDE", 
              msg = "get.dea() is deprecated. Use get.pbDE() instead, which stores results in .tdr.obj$pbDE.")
  
  # Call get.pbDE and extract results for backward compatibility
  result <- get.pbDE(
    x = .tdr.obj,
    .mode = "design",
    .design = .design,
    .contrasts = .contrasts,
    .block = .block,
    .geneset.ls = .geneset.ls,
    .id.idx = .id.idx,
    .id = .id,
    .id.from = .id.from,
    .model.name = ".deprecated.call",
    .result.name = "result",
    .force.recalc = TRUE,
    .verbose = .verbose,
    .label.confidence = .label.confidence
  )
  
  # Return just the DE results (legacy behavior)
  return(result$pbDE$.deprecated.call$result)
}

#' Marker Gene/Protein Identification (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' \code{get.markerDE()} is deprecated. Use \code{\link{get.pbDE}(.mode = "marker")} instead.
#' This function dispatches to \code{get.pbDE} with \code{.mode = "marker"}.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object processed through \code{get.map()}.
#' @param .source The raw data object (or NULL for file backend).
#' @param .geneset.ls Optional named list of character vectors for GSVA. RNA only.
#' @param .id1.idx Optional integer vector of landmark indices for group 1.
#' @param .id2.idx Optional integer vector of landmark indices for group 2.
#' @param .id1 Character vector of cluster/celltype IDs for group 1.
#' @param .id2 Character vector of IDs for group 2. Default \code{"..all.other.landmarks.."}.
#' @param .id.from Character: "clustering" or "celltyping".
#' @param .model.name Character: model name (default "default").
#' @param .comparison.name Character: comparison name, passed as \code{.result.name}.
#' @param .force.recalc Logical: recalculate even if results exist? Default FALSE.
#' @param .label.confidence Numeric (0-1): confidence threshold. Default 0.5.
#' @param ... Additional arguments passed to \code{get.pbDE}.
#'
#' @return The modified \code{x} with results stored in \code{$markerDE}.
#' @seealso \code{\link{get.pbDE}}
#' @keywords internal
#' @export
get.markerDE <- function(x, ...) UseMethod("get.markerDE")

#' @rdname get.markerDE
#' @export
get.markerDE.TDRObj <-
  function(
    x,
    .source = NULL,
    .geneset.ls = NULL,
    .id1.idx = NULL,
    .id2.idx = NULL,
    .id1 = NULL,
    .id2 = "..all.other.landmarks..",
    .id.from = "clustering",
    .model.name = "default",
    .comparison.name = NULL,
    .force.recalc = FALSE,
    .label.confidence = 0.5,
    ...
  ) {
    .Deprecated("get.pbDE",
                msg = paste0(
                  "get.markerDE() is soft-deprecated. ",
                  "Use get.pbDE(x, .mode = 'marker', .id = ..., .id2 = ...) instead.\n",
                  "get.markerDE() will be removed in a future release."
                ))

    get.pbDE.TDRObj(
      x = x,
      .source = .source,
      .mode = "marker",
      .geneset.ls = .geneset.ls,
      .id = .id1,
      .id.idx = .id1.idx,
      .id2 = .id2,
      .id2.idx = .id2.idx,
      .id.from = .id.from,
      .model.name = .model.name,
      .result.name = .comparison.name,
      .force.recalc = .force.recalc,
      .verbose = TRUE,
      .label.confidence = .label.confidence,
      ...
    )
  }

#' Deprecated: Use get.pbDE(.mode = "marker") instead
#'
#' \code{get.marker()} is deprecated. Use \code{\link{get.pbDE}(.mode = "marker")} instead.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .geneset.ls Optional named list of character vectors for GSVA.
#' @param .id1.idx Optional landmark indices for group 1.
#' @param .id2.idx Optional landmark indices for group 2.
#' @param .id1 Cluster/celltype IDs for group 1.
#' @param .id2 Reference group IDs. Default \code{"..all.other.landmarks.."}.
#' @param .id.from "clustering" or "celltyping".
#' @param .label.confidence Numeric (0-1): minimum confidence for cell assignment (default: 0.5).
#'
#' @return A .tdr.obj with results stored in .tdr.obj$markerDE.
#' @seealso \code{\link{get.pbDE}}
#' @export
get.marker <- function(
    .tdr.obj,
    .geneset.ls = NULL,
    .id1.idx = NULL,
    .id2.idx = NULL,
    .id1 = NULL,
    .id2 = "..all.other.landmarks..",
    .id.from = "clustering",
    .label.confidence = 0.5
) {
  .tdr_validate_label_confidence(.label.confidence)

  .Deprecated("get.pbDE",
              msg = "get.marker() is deprecated. Use get.pbDE(.mode = 'marker') instead.")

  # Dispatch directly to get.pbDE marker mode
  return(get.pbDE(
    x = .tdr.obj,
    .mode = "marker",
    .geneset.ls = .geneset.ls,
    .id = .id1,
    .id.idx = .id1.idx,
    .id2 = .id2,
    .id2.idx = .id2.idx,
    .id.from = .id.from,
    .model.name = ".deprecated.call",
    .result.name = "result",
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
    .traj.dist.metric = "cosine",
    .seed = 123
  ){
    
    # Initialize embedding slot if not present
    if(is.null(x = .tdr.obj@sample.embed)){
      .tdr.obj@sample.embed <- list()
    }
    
    # -------------------------------------------------------------------------
    # Compute PCA if not already present
    # -------------------------------------------------------------------------
    
    if(is.null(x = .tdr.obj@sample.embed$pca)){
      
      if(isTRUE(x = .verbose)){
        message("\nComputing sample-level PCA on landmark densities")
      }
      
      # Compute PCA via irlba (samples as rows, using .tdr.obj@density$log.norm)
      set.seed(seed = .seed)
      pca.res <-
        Matrix::t(x = .tdr.obj@density$log.norm) |>
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
      .tdr.obj@sample.embed$pca <- list(
        coord = pca.res$x,
        rotation = pca.res$rotation,
        center = pca.res$center,
        scale = pca.res$scale,
        sdev = pca.res$sdev
      )
      
      # Set column names for coord
      colnames(x = .tdr.obj@sample.embed$pca$coord) <- 
        paste0("PC", seq_len(length.out = ncol(x = .tdr.obj@sample.embed$pca$coord)))
      
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
      
      if(is.null(x = .tdr.obj@sample.embed$traj)){
        
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
        set.seed(seed = .seed)
        Matrix.warnDeprecatedCoerce <-
          # destiny relies on Deprecated Matrix package function: see https://github.com/theislab/destiny/issues/61
          getOption(x = "Matrix.warnDeprecatedCoerce")
        
        options(Matrix.warnDeprecatedCoerce = 0)
        
        traj.res <- 
          Matrix::t(x = .tdr.obj@density$log.norm) |>
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
          rownames(x = .tdr.obj@metadata)
        
        # Store in format similar to pca slot
        .tdr.obj@sample.embed$traj <- list(
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
#' landmark densities stored in \code{.tdr.obj@density$log.norm}. When called with supervised arguments 
#' (\code{.contrast.of.interest} or \code{.red.model}), additionally computes a 
#' partial-effect PCA (pePC) that isolates variation attributable to a specific effect. 
#' Exact for OLS; if duplicateCorrelation/blocking is used, the decomposition is approximate.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object processed through \code{get.map()}. Contains
#'   \code{@density$log.norm} (log2-transformed densities) used for unsupervised embeddings. Statistical
#'   model fits should be stored in \code{@results$lm} via \code{get.lm()}.
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
#' @param .seed Integer: random seed for reproducibility of sample-level PCA
#'   (via \code{irlba::prcomp_irlba}) and diffusion map computation. Default 123.
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
# \code{pca}: Standard PCA on log2-transformed landmark densities (log2(norm + 0.5)),
# centered and scaled.
#
# \code{traj}: Diffusion map trajectory computed via \code{destiny::DiffusionMap} on
# the log2-transformed landmark density matrix (\code{@density$log.norm}). Useful for
# visualizing continuous trajectories.
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
#' direction(s). For rank-1 effects (single contrast or two-level factor), scores
#' increase when residuals align with the unique partial-effect axis. For rank > 1
#' effects, each pePC axis captures a variance-maximizing direction within the
#' multi-dimensional effect subspace (see below).
#'
#' \subsection{Multi-level terms (rank > 1 effects)}{
#'
#' When the dropped term has rank > 1 (e.g., a factor with more than two levels),
#' the nested model method (\code{.red.model}) produces \eqn{\Delta\hat{Y}} with
#' rank \eqn{r = p_{full} - p_{red} > 1}. PCA then learns \eqn{r} orthogonal axes
#' that maximize variance within this effect subspace.
#'
#' \strong{Key property:} these axes are rotations that maximize explained variance.
#' They do \emph{not} correspond to individual level-vs-baseline comparisons. For
#' example, with a three-level factor (Baseline, D1, D7), the nested model pePC
#' produces two axes (pePC1, pePC2) that span the same subspace as the D1-vs-Baseline
#' and D7-vs-Baseline effects, but pePC1 might capture "shared activation" while
#' pePC2 captures the difference between D1 and D7.
#'
#' \strong{The whole-term embedding is correct and parameterization-invariant:}
#' \eqn{\Delta\hat{Y}} depends only on the column spaces of the full and reduced
#' design matrices, not on the choice of factor coding (treatment, sum, Helmert) or
#' reference level. The pePC subspace faithfully represents the complete effect of the
#' dropped term.
#'
#' \strong{For per-level interpretability, use explicit contrasts instead.}
#' The FWL-based method (\code{.contrast.of.interest}) produces rank-1 embeddings
#' that isolate specific comparisons (e.g., D1 vs Baseline, D7 vs Baseline). This
#' requires fitting a cell-means model with \code{limma::makeContrasts}. See
#' Example 3 below.
#'
#' \strong{Practical guidance:}
#' \itemize{
#'   \item Two-level factors or single contrasts: both methods are equivalent and
#'     produce a single, interpretable pePC axis.
#'   \item Multi-level factors, whole-term view: use nested models. Useful for
#'     omnibus visualization ("does this factor matter at all?") and for downstream
#'     analyses that operate on the full effect subspace.
#'   \item Multi-level factors, per-level view: use explicit contrasts. Each call to
#'     \code{get.embedding()} with a different \code{.contrast.of.interest} produces
#'     a separate rank-1 embedding with a clear biological interpretation.
#' }
#' }
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
#'
#' # ===========================================================================
#' # Example 3: Multi-level factor — per-level embeddings via explicit contrasts
#' # ===========================================================================
#'
#' # Suppose Timepoint has three levels: Baseline, D1, D7.
#' #
#' # Option A (nested models): whole-term embedding with rank 2.
#' # pePC1 and pePC2 are variance-maximizing rotations within the
#' # Timepoint effect subspace — they do NOT correspond to individual
#' # level-vs-baseline comparisons.
#'
#' full.design  <- model.matrix(~ Timepoint + Batch, data = lm.obj$metadata)
#' red.design   <- model.matrix(~ Batch, data = lm.obj$metadata)
#' lm.obj <- get.lm(lm.obj, .design = full.design,  .model.name = "tp_full")
#' lm.obj <- get.lm(lm.obj, .design = red.design,   .model.name = "tp_red")
#'
#' lm.obj <- get.embedding(lm.obj,
#'     .full.model = "tp_full", .red.model = "tp_red",
#'     .term.of.interest = "Timepoint")
#' # -> produces pePC1, pePC2: omnibus Timepoint effect
#'
#' # Option B (explicit contrasts): one rank-1 embedding per comparison.
#' # Each pePC axis directly captures a specific level-vs-baseline effect.
#'
#' cm.design <- model.matrix(~ 0 + Timepoint + Batch, data = lm.obj$metadata)
#' tp.contrasts <- limma::makeContrasts(
#'     D1vsBaseline = TimepointD1  - TimepointBaseline,
#'     D7vsBaseline = TimepointD7  - TimepointBaseline,
#'     levels = cm.design
#' )
#' lm.obj <- get.lm(lm.obj, .design = cm.design,
#'                   .contrasts = tp.contrasts,
#'                   .model.name = "tp_contrasts")
#'
#' lm.obj <- get.embedding(lm.obj,
#'     .full.model = "tp_contrasts",
#'     .contrast.of.interest = "D1vsBaseline")
#' lm.obj <- get.embedding(lm.obj,
#'     .full.model = "tp_contrasts",
#'     .contrast.of.interest = "D7vsBaseline")
#'
#' # Each embedding is rank 1 with a clear biological meaning:
#' plotSampleEmbedding(lm.obj, .embedding = "pePC",
#'     .sup.embed.slot = "D1vsBaseline", .color.by = "Timepoint")
#' plotSampleEmbedding(lm.obj, .embedding = "pePC",
#'     .sup.embed.slot = "D7vsBaseline", .color.by = "Timepoint")
#' }
#'
#' @param ... Additional arguments passed to methods.
#' @export
#'
get.embedding <- function(x, ...) UseMethod("get.embedding")

#' @rdname get.embedding
#' @export
get.embedding.TDRObj <-
  function(
    x,
    .full.model = "default",
    .term.of.interest = NULL,
    .red.model = NULL,
    .contrast.of.interest = NULL,
    .n.eigs = 20,
    .n.pcs = 20,
    .ret.trajectory = FALSE,
    .traj.dist.metric = "cosine",
    .seed = 123,
    .verbose = TRUE,
    ...
  ){
    .tdr.obj <- x
    
    # Protect caller's RNG state; restore on exit (safe in fresh sessions)
    withr::local_preserve_seed()
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    # Check that .tdr.obj has @density$log.norm
    if(is.null(x = .tdr.obj@density$log.norm)){
      stop("'.tdr.obj@density$log.norm' not found. Run get.map() first.")
    }
    
    # Determine if this is unsupervised-only mode (no pePC args provided)
    unsupervised.only <- 
      is.null(x = .red.model) && is.null(x = .contrast.of.interest)
    
    # If supervised args provided, require the full model to exist
    if(!unsupervised.only){
      
      if(is.null(x = .tdr.obj@results$lm[[.full.model]])){
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
    .stats.obj <- .tdr.obj@results$lm[[.full.model]]
    
    # -------------------------------------------------------------------------
    # Compute unsupervised embeddings (pca, traj) - stored in .tdr.obj@sample.embed
    # -------------------------------------------------------------------------
    
    .tdr.obj <- 
      .compute.unsupervised.embeddings(
        .tdr.obj = .tdr.obj,
        .n.eigs = .n.eigs,
        .n.pcs = .n.pcs,
        .ret.trajectory = .ret.trajectory,
        .traj.dist.metric = .traj.dist.metric,
        .seed = .seed,
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
      if(!(.term.of.interest %in% colnames(x = .tdr.obj@metadata))){
        
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
      
      # Get the log-normalized density matrix (landmarks x samples)
      Y <- 
        .tdr.obj@density$log.norm
      
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
            lapply(X = .tdr.obj@metadata,
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
                                          block = .tdr.obj@metadata[[.block]])
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
                       block = if(!is.null(dupcor)) .tdr.obj@metadata[[.block]] else NULL,
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
      if(is.null(x = .tdr.obj@results$lm[[.red.model]])){
        stop(paste0("Reduced model '", .red.model, "' not found in .tdr.obj$map$lm. ",
                    "Run get.lm() with .model.name = '", .red.model, "' first."))
      }
      
      .red.stats.obj <- .tdr.obj@results$lm[[.red.model]]
      
      if(isTRUE(x = .verbose)){
        message("\nComputing embedding via nested model comparison")
        message("  Slot name: ", slot.name)
        message("  Full model ('", .full.model, "') columns: ", paste(colnames(.stats.obj$fit$design), collapse = ", "))
        message("  Reduced model ('", .red.model, "') columns: ", paste(colnames(.red.stats.obj$fit$design), collapse = ", "))
      }
      
      # Get log.norm from .tdr.obj@density$log.norm for dimension validation
      Y <- .tdr.obj@density$log.norm
      
      # Validate dimensions against model fits
      if(nrow(x = .stats.obj$fit$coefficients) != nrow(x = Y)){
        stop("Number of landmarks in full model does not match .tdr.obj@density$log.norm")
      }
      if(nrow(x = .red.stats.obj$fit$coefficients) != nrow(x = Y)){
        stop("Number of landmarks in reduced model does not match .tdr.obj@density$log.norm")
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
    
    if(embed.rank > 1L && method == "nested_models" && isTRUE(x = .verbose)){
      message(
        "\n  Note: Term '", slot.name, "' has rank ", embed.rank,
        ". pePC axes are variance-maximizing\n",
        "  rotations within the full effect subspace -- they do NOT correspond\n",
        "  to individual level-vs-baseline comparisons. For per-level embeddings,\n",
        "  use get.embedding() with .contrast.of.interest and explicit contrasts.\n",
        "  See ?get.embedding section 'Multi-level terms' for details."
      )
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
      Matrix::t(x = (.tdr.obj@density$log.norm - Yhat.red) - pca$center) %*%
      pca$rotation
    
    colnames(x = pca$coord) <- 
      paste0("pePC", 
             seq_len(length.out = ncol(x = pca$coord)))
    
    # partial-effect principal component variance as a fraction of TOTAL variance
    perc.tot.var.exp <-
      (100 * ((pca$sdev[1:ncol(x = pca$coord)])^2) /
         (matrixStats::rowVars(x = .tdr.obj@density$log.norm) |>
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
    # Store results in .tdr.obj@sample.embed$pepc[[slot.name]]
    # -------------------------------------------------------------------------
    
    # Initialize pePC list if needed
    if(is.null(x = .tdr.obj@sample.embed$pepc)){
      .tdr.obj@sample.embed$pepc <- list()
    }
    
    .tdr.obj@sample.embed$pepc[[slot.name]] <- 
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
# Internal helpers shared by get.plsD
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
    
    if (is.null(x = .tdr.obj@results$lm[[.model.name]])) {
      stop("Model '", .model.name, "' not found. Run get.lm() first.")
    }
    
    coef.mat <-
      .tdr.obj@results$lm[[.model.name]]$fit$coefficients
    
    if (!(.coef.col %in% colnames(x = coef.mat))) {
      stop("Coefficient '", .coef.col, "' not found in model coefficients.\n",
           "Available: ", paste(colnames(x = coef.mat), collapse = ", "))
    }
    
    if (is.null(x = .tdr.obj@graphs$snn)) {
      stop("SNN graph not found. Run get.graph() first.")
    }
    
    if (!Matrix::isSymmetric(object = .tdr.obj@graphs$snn)) {
      stop("SNN graph not symmetric.")
    }
    
    if (is.null(x = .tdr.obj@assay$raw)) {
      stop("Raw landmarks not found. Run get.landmarks() first.")
    }
    
    if (is.null(x = .tdr.obj@landmark.embed$pca$coord)) {
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
# .build.nuisance.Z
#
# Extracts the nuisance portion of the design matrix, expanded to landmark
# level via @config$key.  Optionally augments Z with dummy-coded blocking
# columns when fit$block is non-NULL (i.e. get.lm() was called with .block).
#
# Block-awareness rationale:
#   When get.lm() uses limma::duplicateCorrelation with a blocking variable
#   (e.g. donor/subject), the density contrast Y is estimated via GLS that
#   accounts for within-block correlation.  However, the blocking factor does
#   NOT appear in the design matrix — it enters only through the covariance
#   structure.  Without explicit handling, the expression matrix X retains
#   between-block mean shifts that can contaminate PLS directions.
#
#   The fix is to include the blocking factor as additional nuisance columns in
#   Z so that the OLS projection (I - H_Z) X removes block-level expression
#   means.  This preserves the implicit-operator architecture and does not
#   require GLS on X (which would be both computationally infeasible at
#   landmark scale and statistically inappropriate for observed data).
#
# Returns: expanded nuisance matrix Z (n_landmarks x q_z), or NULL if no
#          nuisance columns exist after accounting for design + block.
# -----------------------------------------------------------------------------

.build.nuisance.Z <-
  function(.tdr.obj,
           .model.name,
           .coef.col,
           .verbose = TRUE) {

    fit <-
      .tdr.obj@results$lm[[.model.name]]$fit

    design <-
      fit$design                       # original design (unchanged by contrasts.fit)

    if (!is.null(x = fit$contrasts)) {
      # Contrasts present: group = columns participating in any contrast
      contrast.participation <-
        Matrix::rowSums(x = abs(x = fit$contrasts)) > 0
      nuisance.cols <-
        which(x = !contrast.participation)
    } else {
      # No contrasts: group = .coef.col column
      coef.col.idx <-
        which(x = colnames(x = design) == .coef.col)
      if (length(x = coef.col.idx) == 0) {
        stop("'.coef.col' not found in design columns: ",
             paste(colnames(x = design), collapse = ", "))
      }
      nuisance.cols <-
        setdiff(x = seq_len(length.out = ncol(x = design)),
                y = coef.col.idx)
    }

    # -------------------------------------------------------------------------
    # Build sample-level nuisance matrix from design columns
    # -------------------------------------------------------------------------

    if (length(x = nuisance.cols) > 0) {
      Z.sample <-
        design[, nuisance.cols, drop = FALSE]
    } else {
      Z.sample <- NULL
    }

    # -------------------------------------------------------------------------
    # Detect and append blocking variable (if used in get.lm)
    #
    # fit$block is a character/factor vector of length n_samples stored by
    # limma::lmFit when block= is supplied.  We reverse-lookup the metadata
    # column that matches it (same approach as get.embedding's pePC path).
    # The block factor is dummy-encoded with treatment contrasts; the intercept
    # (already guaranteed in Z) absorbs the reference level.
    # -------------------------------------------------------------------------

    block.cols.sample <- NULL

    if (!is.null(x = fit$block)) {

      # Reverse-lookup: find the metadata column matching fit$block
      block.name <-
        lapply(X = .tdr.obj@metadata,
               FUN = function(meta.col) {
                 if (length(x = meta.col) != length(x = fit$block)) return(FALSE)
                 all(as.character(x = meta.col) == as.character(x = fit$block))
               }) |>
        unlist(use.names = TRUE) |>
        which() |>
        names()

      if (length(x = block.name) > 0) {
        block.name <- block.name[1]

        block.factor <-
          factor(x = fit$block)

        n.levels <- nlevels(x = block.factor)

        if (n.levels > 1) {
          # Treatment-coded dummies (drops reference = first level)
          # With intercept in Z, this gives full-rank parameterization
          block.dummies <-
            stats::model.matrix(object = ~ block.factor)[, -1, drop = FALSE]

          # Clean column names: "block.factorB" -> "block:B"
          colnames(x = block.dummies) <-
            gsub(pattern = "^block\\.factor",
                 replacement = paste0("block:", block.name, ":"),
                 x = colnames(x = block.dummies))

          block.cols.sample <- block.dummies

          if (isTRUE(x = .verbose)) {
            message("  Including blocking variable '", block.name,
                    "' (", n.levels, " levels) in nuisance design.")
          }
        }
      }
    }

    # -------------------------------------------------------------------------
    # Combine design-based nuisance + block dummies at sample level
    # -------------------------------------------------------------------------

    if (is.null(x = Z.sample) && is.null(x = block.cols.sample)) {
      return(NULL)
    }

    if (!is.null(x = Z.sample) && !is.null(x = block.cols.sample)) {
      Z.sample <- cbind(Z.sample, block.cols.sample)
    } else if (is.null(x = Z.sample)) {
      Z.sample <- block.cols.sample
    }

    # -------------------------------------------------------------------------
    # Expand to landmark level via @config$key
    # -------------------------------------------------------------------------

    Z <-
      Z.sample[.tdr.obj@config$key, , drop = FALSE]

    # Ensure intercept is present (needed for proper centering)
    has.intercept <-
      apply(X = Z,
            MARGIN = 2,
            FUN = function(col) all(col == 1)) |>
      any()

    if (!has.intercept) {
      Z <- cbind("(Intercept)" = 1, Z)
    }

    return(Z)
  }


# -----------------------------------------------------------------------------
# .build.resid.operators
# Precomputes the implicit residualization operators:
#   B = (Z'Z)^{-1} Z' X          (q_z x p, dense)
#   muX.resid = muX - zbar' B    (length p)
#   ZtZ = Z'Z                     (q_z x q_z, for loading SD correction)
#
# Uses pivoted QR if Z'Z is near-singular.
# Returns a list; NULL Z input returns NULL.
# -----------------------------------------------------------------------------

.build.resid.operators <-
  function(Z,
           X.sparse,
           muX,
           .verbose = TRUE) {

    if (is.null(x = Z)) return(NULL)

    q.z <- ncol(x = Z)
    n   <- nrow(x = Z)

    ZtZ <-
      crossprod(x = Z)                              # q_z x q_z

    # Check conditioning
    rc <- rcond(x = ZtZ)
    if (rc < .Machine$double.eps * 100) {
      # Rank-deficient: reduce Z via pivoted QR
      qr.Z <- qr(x = Z)
      keep  <- seq_len(length.out = qr.Z$rank)
      Z     <- qr.Q(qr = qr.Z)[, keep, drop = FALSE]
      ZtZ   <- crossprod(x = Z)
      q.z   <- ncol(x = Z)
      if (isTRUE(x = .verbose)) {
        warning("Nuisance design is rank-deficient after landmark expansion. ",
                "Reduced to rank ", q.z, ".")
      }
    }

    ZtZ.inv <-
      solve(a = ZtZ)                                # q_z x q_z

    ZtX <-
      as.matrix(
        x = Matrix::crossprod(x = Z, y = X.sparse)  # q_z x p
      )

    B <-
      ZtZ.inv %*% ZtX                               # q_z x p

    z.bar <-
      colMeans(x = Z)                                # length q_z

    muX.resid <-
      muX - as.numeric(x = crossprod(x = B, y = z.bar))  # length p

    list(B         = B,
         Z         = Z,
         ZtZ       = ZtZ,
         z.bar     = z.bar,
         muX.resid = muX.resid)
  }


# -----------------------------------------------------------------------------
# .prepare.X
# Filters, normalizes (size-factor + log2 for RNA), and optionally centers the
# expression matrix.  Returns a list with:
#   $X   — pre-centering (sparse; useful for sparse-aware loadings)
#   $muX — column means of X (numeric; always returned for implicit-centering)
#   $Xc  — column-centered (dense; only when .center = TRUE)
# -----------------------------------------------------------------------------

.prepare.X <-
  function(.tdr.obj,
           .min.prop = 0.005,
           .center = TRUE,
           .verbose = TRUE) {
    
    if (.tdr.obj@config$assay.type == "RNA") {
      
      # Filter genes: detected in at least min.prop of landmarks
      X <-
        .tdr.obj@assay$raw |>
        (\(x)
         x[, Matrix::colSums(x = x > 0) > (nrow(x = x) * .min.prop)]
        )()
      
      if (isTRUE(x = .verbose)) {
        message("  Genes after filtering (>", .min.prop * 100, "% detection): ",
                ncol(x = X), " / ", ncol(x = .tdr.obj@assay$raw))
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
        .tdr.obj@assay$raw
      
    }
    
    # Column means (sparse-aware, O(nnz))
    muX <-
      Matrix::colMeans(x = X)
    
    # Center expression (genes/markers) — only when requested
    if (isTRUE(x = .center)) {
      Xc <-
        X |>
        (\(x)
         Matrix::t(x = x) - muX
        )() |>
        Matrix::t()
    } else {
      Xc <- NULL
    }
    
    return(list(X = X, muX = muX, Xc = Xc))
    
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
           .method = c("pearson", "ols", "spearman"),
           Y0.vec = NULL,
           post.sd.vec = NULL,
           muX = NULL,
           resid.ops = NULL) {
    
    .method <-
      match.arg(arg = .method,
                choices = c("pearson", "ols", "spearman"))
    
    .do.resid <- !is.null(x = resid.ops)

    muX.eff <-
      if (.do.resid) resid.ops$muX.resid else muX

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

      if (.do.resid) {
        # Adjust sufficient statistics for Xtilde = X - Z B
        # colSums(Xtilde) = colSums(X) - n * zbar' B
        col.sums <- col.sums -
          n * as.numeric(x = crossprod(x = resid.ops$B,
                                       y = resid.ops$z.bar))

        # colSumSq(Xtilde) = colSumSq(X) - diag(B' Z'Z B)
        R <- Matrix::chol(x = resid.ops$ZtZ)
        RB <- R %*% resid.ops$B                          # q_z x p
        col.sum.sq <- col.sum.sq -
          colSums(x = RB^2)

        gene.sd <-
          sqrt(x = pmax(col.sum.sq - (col.sums^2) / n, 0) /
                 (n - 1))
      }
      
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
          
          # Residualization correction to numerator
          if (.do.resid) {
            Zt.score <-
              as.numeric(x = crossprod(x = resid.ops$Z, y = score.k))
            num <- num -
              as.numeric(x = crossprod(x = resid.ops$B, y = Zt.score)) /
              (n - 1)
          }

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
          
          num.ols <-
            Matrix::crossprod(x = score.k,
                              y = X.sparse) |>
            as.numeric()

          if (.do.resid) {
            Zt.score <-
              as.numeric(x = crossprod(x = resid.ops$Z, y = score.k))
            num.ols <- num.ols -
              as.numeric(x = crossprod(x = resid.ops$B, y = Zt.score))
          }

          beta <- num.ols / denom
          
          return(beta)
          
        }) |>
        do.call(what = cbind)
      
    }
    
    # -----------------------------------------------------------------------
    # Per-gene soft concordance weights
    #
    # c_jk = clip( A_jk / N_jk, 0, 1 ) where
    #   N_jk = t_kc' Xc  (full loading numerator; sum(t_kc) = 0 so muX drops out)
    #   A_jk = (s_k * t_kc)' Xc  (soft-concordant portion)
    #   s_ki = pnorm( Y0_i * sign(t_kci) / post.sd_i )
    #   post.sd_i defaults to 0 (hard threshold) if post.sd.vec is NULL
    #
    # Not defined for spearman (NA).
    # Not computed when Y0.vec is NULL.
    # -----------------------------------------------------------------------

    concordance.weights <-
      if (is.null(x = Y0.vec) || .method == "spearman") {
        
        matrix(data = NA_real_,
               nrow = n.genes,
               ncol = ncomp,
               dimnames = list(colnames(x = X.sparse), comp.names))
        
      } else {
        
        sapply(
          X   = seq_len(length.out = ncomp),
          FUN = function(k) {
            
            score.k <- scores[, k] - mean(x = scores[, k])
            
            # Full loading numerator: t_kc' X (muX term = 0 since sum(t_kc)=0)
            full.num <-
              as.numeric(x = Matrix::crossprod(x = X.sparse, y = score.k))
            
            if (.do.resid) {
              Zt.score <-
                as.numeric(x = crossprod(x = resid.ops$Z, y = score.k))
              full.num <- full.num -
                as.numeric(x = crossprod(x = resid.ops$B, y = Zt.score))
            }
            
            # Soft concordance indicator: pnorm( Y0 * sign(t_kc) / post.sd )
            # When t_kc = 0: sign returns 0 -> argument = 0 -> pnorm(0) = 0.5 (neutral)
            # When post.sd.vec is NULL: use hard threshold (indicator of sign agreement)
            if (is.null(x = post.sd.vec)) {
              s.vec <- as.numeric(x = sign(x = score.k) == sign(x = Y0.vec) &
                                    Y0.vec != 0)
            } else {
              s.vec <- stats::pnorm(
                q = Y0.vec * sign(x = score.k) / post.sd.vec
              )
              # sign(0) = 0 -> argument = 0 -> pnorm(0) = 0.5; already correct
            }
            
            # Concordant numerator: (s * t_kc)' Xc
            # Note: sum(s * t_kc) != 0 in general, so muX correction IS needed
            s.t <- s.vec * score.k
            conc.num <-
              as.numeric(x = Matrix::crossprod(x = X.sparse, y = s.t))

            if (.do.resid) {
              Zt.st <-
                as.numeric(x = crossprod(x = resid.ops$Z, y = s.t))
              conc.num <- conc.num -
                as.numeric(x = crossprod(x = resid.ops$B, y = Zt.st))
            }

            conc.num <- conc.num - muX.eff * sum(s.t)
            
            # Concordance fraction
            near.zero <-
              abs(x = full.num) <
              (.Machine$double.eps * max(abs(x = full.num), 1))
            
            cw <-
              ifelse(
                test = near.zero,
                yes  = 0.5,
                no   = conc.num / full.num
              )
            
            # Clip to [0, 1]
            pmax(0, pmin(1, cw))
            
          }
        )
        
      }

    if (!is.null(x = Y0.vec) && .method != "spearman") {
      rownames(x = concordance.weights) <- colnames(x = X.sparse)
      colnames(x = concordance.weights) <- comp.names
    }

    rownames(x = loadings) <-
      colnames(x = X.sparse)
    
    return(
      list(
        raw.loadings         = loadings,
        concordance.weights  = concordance.weights
      )
    )
    
  }


#' Graph-Diffused, Density Contrast-Aligned PLS Decomposition (plsD)
#'
#' Decomposes an expression interaction matrix M.local via NIPALS PLS1
#' against the density-contrast vector Y. plsD generates candidate features that drive
#' the density contrast — including population markers (DA), differentially
#' expressed genes (DE), and their mixture — by maximizing covariance between
#' graph-smoothed expression and Y.
#'
#' @details
#' plsD is an \strong{interpretive decomposition}, not a formal hypothesis test.
#' It answers: "What features — through expression patterns, population identity,
#' or both — explain the density contrast?" Loadings reflect combined signal, like:
#' markers of depleted populations carry negative loadings (high expression where
#' density is low), DE genes carry signed loadings tracking their direction of
#' change, and mixed DA/DE features appear alongside both.
#'
#' When \code{.YX.interaction = TRUE} (default), the interaction term diag(Y)
#' bakes the density contrast into the data matrix:
#' M.local = P \%*\% diag(Y) \%*\% Xc. This gives plsD comprehensive scope:
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
#' smoothed space. In the presence of strongly bimodal density contrasts and/or,
#' datasets with a small number of extreme landmarks, setting 
#' \code{.YX.interaction = FALSE} can help distinguish features whose expression 
#' genuinely covaries with Y from those that are strong but purely geometric
#' counterweights.
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
#'         Optionally, residualize Xc against nuisance covariates before constructing M.local.
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
#' are ordered by extraction (plsD1 = highest covariance with Y).
#'
#' \strong{Sk (graph smoothness):} Derived from the normalized Laplacian applied to
#' score vectors. High Sk indicates large-scale, graph-smooth structure.
#'
#' \strong{About Y appearing on both sides (when .YX.interaction = TRUE):}
#' Y appears in both the data matrix (via diag(Y) in M.local) and the PLS objective
#' (maximize covariance with Y). This is intentional: it gives plsD comprehensive
#' sensitivity to features driving density changes. The design is not circular
#' because the expression matrix Xc sits between: the product is large only when
#' features exist whose Y-weighted, graph-smoothed expression genuinely covaries
#' with Y. With permuted Y, Ak collapses. Setting \code{.YX.interaction = FALSE}
#' removes Y from the data matrix entirely, providing a diagnostic comparator.
#'
#' \strong{Nuisance residualization (.residualize = TRUE):}
#' When the experimental design includes nuisance covariates (e.g. batch, sex,
#' age), the density contrast Y from \code{get.lm} is already adjusted for these.
#' However, the expression matrix X is not. Setting \code{.residualize = TRUE}
#' projects X onto the orthogonal complement of the nuisance design, so that
#' both sides of the PLS objective are free of nuisance effects. This improves
#' interpretational symmetry: loadings reflect features whose expression covaries
#' with the \emph{adjusted} density contrast, with nuisance-driven expression
#' variation removed. The residualization uses the Frisch-Waugh-Lovell (FWL)
#' decomposition: X is replaced by (I - H_Z) X, where H_Z is the hat matrix of
#' the nuisance design expanded to landmark level via \code{@config$key}. Sparsity
#' is preserved via implicit operators. When no nuisance covariates exist (e.g.
#' a simple two-group design), \code{.residualize = TRUE} is equivalent to
#' \code{.residualize = FALSE}.
#'
#' \strong{Blocking variable handling:}
#' When \code{get.lm()} was called with a \code{.block} argument (invoking
#' \code{limma::duplicateCorrelation}), the blocking factor does NOT appear in
#' the design matrix — it enters only through the GLS covariance structure for
#' the density contrast Y. If \code{.residualize = TRUE}, the blocking variable
#' is automatically detected from \code{fit$block} and appended to the nuisance
#' design Z as dummy-coded columns. This ensures that block-level mean expression
#' shifts (e.g. donor-specific expression baselines) are removed from X, preventing
#' them from contaminating PLS directions via spurious covariance with the
#' block-correlated residual structure in Y.
#'
#' \strong{Structural score constraints and the balancing effect:}
#' The score vectors computed by plsD are structurally mean-zero across landmarks. This
#' is a mathematical consequence of column-centering \code{M.local}: for any weight vector
#' \code{w}, the score \code{t = M.local w} satisfies \code{sum(t) = 0}. As a result, a
#' dominant positive pole (landmarks with large positive scores) must be balanced by
#' compensatory negative scores somewhere on the manifold.
#'
#' In datasets where a small number of biologically extreme landmarks simultaneously have
#' large |Y_c| and unusual expression profiles (large row-wise expression norm), these
#' landmarks can dominate the first PLS component. The compensatory negative region may
#' then appear coherent and strong, mimicking a genuine opposing biological program.
#'
#' \strong{Distinguishing genuine signal from structural counterweight}: the
#' \code{plotPlsD} scatter panel colors landmarks by their \emph{raw} (uncentered)
#' density contrast coefficient. If landmarks with large negative scores show near-zero
#' or positive raw coefficients, they are geometric counterweights rather than genuinely
#' depleted populations. Conversely (depending on contrast direction and component),
#' large-magnitude positive-score landmarks with near-zero raw Y may also reflect
#' structural balance.
#'
#' \strong{Pre-analysis diagnostic}: if unexpected strong balancing is observed, compute
#' the row-wise interaction norm distribution before running plsD:
#' \preformatted{
#' Yc <- coef.mat[, .coef.col]; Yc <- Yc - mean(Yc)
#' # Xc: size-factor-normalized, log2-transformed, column-centered expression
#' row_norms <- sqrt(rowSums((Yc * as.matrix(Xc))^2))
#' quantile(row_norms)
#' }
#' A heavily right-skewed distribution (max >> Q75, or top 10 landmarks holding >20% of
#' Frobenius norm squared) indicates dominant leverage; reduce \code{.lazy.alpha} or
#' winsorize \code{Yc} before proceeding.
#'
#' @note plsD is designed as an exploratory tool to help interpret density
#'   changes in terms of the features driving them. It does not provide
#'   gene-level p-values or formal multiple testing correction. For rigorous
#'   differential expression testing with p-values following field standards,
#'   use \code{\link{get.pbDE}} (pseudobulk DE via edgeR/limma, both design and marker modes).
#'   plsD complements these
#'   methods by providing a multivariate, graph-aware decomposition that
#'   captures joint DA/DE patterns not visible in gene-by-gene tests.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object after \code{get.lm()}.
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
#'   Values < 1 reduce the spatial propagation range of expression scores on the graph,
#'   making the decomposition more local. This is particularly useful when the density
#'   contrast is dominated by a single biologically extreme population: reducing
#'   \code{.lazy.alpha} limits how far the structural score counterweight (see Details)
#'   propagates to distant manifold regions. Try 0.5 as a first step when strong
#'   unexpected balancing is observed. Values < 1 also incidentally damp oscillations
#'   on sparse/irregular graphs.
#' @param .YX.interaction Logical: if TRUE (default), construct
#'   M.local = P diag(Y) Xc (Y-weighted interaction). Loadings capture candidate features
#'   driving the density contrast through both differential abundance (population
#'   markers) and differential expression. If FALSE, construct M.local = P Xc
#'   (graph-smoothed expression only; Y appears only on the response side).
#'   The FALSE mode is less comprehensive — it misses DA markers unless their
#'   expression independently covaries with Y.
#'   Note: Y is centered before forming diag(Y), so the interaction weights landmarks
#'   by deviation from the mean coefficient, not from zero. Landmarks with below-average
#'   (but still positive) raw coefficients receive a negative weight, contributing to
#'   the structural opposite pole in the decomposition.
#' @param .loading.method Character: method for computing feature loadings. 
#'   \code{"pearson"} computes Pearson correlation between component scores
#'   and centered expression (scale-invariant, fast sparse-BLAS path).
#'   \code{"ols"} computes OLS regression coefficients of centered Xc on component
#'   scores (same sparse-BLAS path; preserves magnitude information but is
#'   scale-dependent, so high-variance features rank higher). \code{"spearman"}
#'   computes Spearman rank correlation via a sparse-aware implementation
#'   (robust to outliers and nonlinearity; slower, O(nnz * log(m_g) * K)).
#' @param .residualize Logical: if TRUE, project out nuisance covariates from
#'   the expression matrix before decomposition. Default FALSE. Nuisance columns
#'   are identified automatically: when contrasts are present, design columns
#'   that do not participate in any contrast are nuisance; when no contrasts are
#'   used, all design columns except \code{.coef.col} are nuisance. If
#'   \code{get.lm()} was called with \code{.block}, the blocking variable is
#'   additionally included in the nuisance design (dummy-coded). The design
#'   is expanded from sample level to landmark level via \code{@config$key}.
#'   Sparsity of the expression matrix is preserved via implicit operators (the
#'   residualized matrix is never materialized). Not supported with
#'   \code{.loading.method = "spearman"}.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'
#' @return The modified \code{.tdr.obj} with results stored in \code{.tdr.obj$plsD[[.coef.col]]}:
#'   \describe{
#'     \item{scores}{Matrix (landmarks x K): PLS scores (oriented so positive = aligned with Y)}
#'     \item{feature.weights}{Matrix (genes x K): PLS feature weight vectors w (unit-norm)}
#'     \item{x.loadings}{Matrix (genes x K): PLS deflation loadings p}
#'     \item{y.loadings}{Numeric vector (K): PLS Y-loadings q (scalar per component)}
#'     \item{raw.loadings}{Matrix (genes x K): raw feature loadings per component (Pearson r,
#'       OLS regression, or Spearman rank correlation, depending on
#'       \code{.loading.method}); positive = upregulated with Y,
#'       negative = downregulated}
#'     \item{loadings}{Matrix (genes x K): concordance-filtered feature loadings,
#'       defined as \code{raw.loadings * concordance.weights} for Pearson/OLS.
#'       For \code{.loading.method = "spearman"}, \code{concordance.weights}
#'       are \code{NA} and \code{loadings = raw.loadings}. Positive = upregulated
#'       with Y, negative = downregulated.}
#'     \item{concordance.weights}{Matrix (genes/markers x K): per-gene soft
#'       concordance weight \eqn{c_{jk} \in [0,1]} for each component. For gene
#'       \eqn{j} and component \eqn{k}, \eqn{c_{jk}} is the fraction of the OLS
#'       loading numerator (\eqn{t_{k}^\top X_{c,j}}) attributable to landmarks
#'       where the PLS score and the raw density contrast coefficient agree in sign,
#'       weighted by the posterior probability that each landmark's coefficient is
#'       genuinely nonzero and concordant:
#'       \eqn{s_{k,i} = \Phi(Y_{0,i} \cdot \operatorname{sign}(t_{k,c,i}) /
#'       \hat\sigma_i)}, where \eqn{\hat\sigma_i = \sqrt{s^{2,\text{post}}_i}
#'       \cdot \text{stdev.unscaled}_{i,c}} is the limma eBayes posterior standard
#'       deviation. \eqn{c_{jk} = 1}: loading arises entirely from concordant,
#'       statistically confident landmarks; \eqn{c_{jk} = 0}: loading arises
#'       entirely from discordant or high-uncertainty landmarks (structural
#'       counterweight); \eqn{c_{jk} = 0.5}: neutral (near-zero loading or
#'       near-zero/highly-uncertain coefficient). \code{NA} for
#'       \code{.loading.method = "spearman"} (rank-based numerators are not
#'       additively decomposable). Compare \code{raw.loadings} to
#'       \code{loadings} to assess attenuation of structural balancing effects.}
#'     \item{Y.alignment}{Numeric vector (K): |cor(Y, score_k)| per component}
#'     \item{smoothness}{Numeric vector (K): Laplacian smoothness per score vector}
#'     \item{Y}{Numeric vector: the density contrast used, centered to mean zero (Y - mean(Y)).
#'       Note: this differs from the raw coefficient by a constant shift. Use the raw
#'       coefficient from the model fit for scientific interpretation of absolute density changes.}
#'     \item{Y.mean}{Numeric: mean of the density contrast used (mean(Y))}
#'     \item{M.local}{Matrix (optional): the interaction matrix (if .store.M = TRUE)}
#'     \item{params}{List: parameters used, including:
#'       \code{params$residualize} (logical: whether residualization was applied) and
#'       \code{params$nuisance.cols} (character: names of nuisance design columns,
#'       NULL if not residualized)}
#'   }
#'
#' @seealso \code{\link{get.lm}} (required predecessor),
#'   \code{\link{plotPlsD}} (visualization), \code{\link{plotPlsDHeatmap}} (heatmap)
#'
#' @examples
#' \dontrun{
#' # After fitting linear model
#' lm.obj <- get.lm(lm.obj, .design = design)
#'
#' # Run plsD for "Infection" coefficient
#' lm.obj <- get.plsD(lm.obj, .coef.col = "Infection")
#'
#' # Access results
#' lm.obj$plsD$Infection$scores[, "plsD1"]      # PLS scores (Y-aligned)
#' lm.obj$plsD$Infection$loadings[, "plsD1"]      # concordance-filtered gene loadings
#' lm.obj$plsD$Infection$raw.loadings[, "plsD1"]  # raw gene loadings before concordance filtering
#'
#' # Diagnostic table
#' data.frame(
#'   component = colnames(lm.obj$plsD$Infection$scores),
#'   Ak = lm.obj$plsD$Infection$Y.alignment,
#'   Sk = lm.obj$plsD$Infection$smoothness
#' )
#' }
#'
#' @param ... Additional arguments passed to methods.
#' @export
#'
get.plsD <- function(x, ...) UseMethod("get.plsD")

#' @rdname get.plsD
#' @export
get.plsD.TDRObj <-
  function(x,
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
           .residualize = FALSE,
           .verbose = TRUE,
           ...) {
    .tdr.obj <- x
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    if (!is.logical(x = .YX.interaction) || length(x = .YX.interaction) != 1) {
      stop(".YX.interaction must be TRUE or FALSE.")
    }

    if (!is.logical(x = .residualize) || length(x = .residualize) != 1) {
      stop(".residualize must be TRUE or FALSE.")
    }

    .loading.method <-
      match.arg(arg = .loading.method,
                choices = c("pearson",
                            "ols",
                            "spearman"))

    if (isTRUE(x = .residualize) &&
        .loading.method == "spearman") {
      stop(".residualize = TRUE is not supported with .loading.method = 'spearman'.\n",
           "Rank-based correlations cannot be computed via the implicit ",
           "residualization operator.\n",
           "Use .loading.method = 'pearson' or 'ols' instead.")
    }
    
    coef.mat <-
      .validate.DE.inputs(.tdr.obj = .tdr.obj,
                          .model.name = .model.name,
                          .coef.col = .coef.col)
    
    # Default .ncomp to number of PCs from get.landmarks
    if (is.null(x = .ncomp)) {
      .ncomp <-
        ncol(x = .tdr.obj@landmark.embed$pca$coord)
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
      message("\n=== plsD: Graph-Diffused, Density Contrast-Aligned PLS Decomposition ===")
      message("Coefficient: ", .coef.col)
    }
    
    Y <-
      coef.mat[, .coef.col] |>
      (\(x)
       x - mean(x = x)
      )()
    
    # -------------------------------------------------------------------------
    # Prepare expression matrix (sparse only; no dense Xc allocated)
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Preparing expression matrix...")
    }
    
    prep <-
      .prepare.X(.tdr.obj = .tdr.obj,
                 .min.prop = .min.prop,
                 .center = FALSE,
                 .verbose = .verbose)
    
    X   <- prep$X     # sparse  n x p
    muX <- prep$muX   # numeric length p
    
    n.landmarks <-
      nrow(x = X)
    
    n.genes <-
      ncol(x = X)
    
    gene.names <-
      colnames(x = X)
    
    # -------------------------------------------------------------------------
    # Nuisance residualization (optional)
    # -------------------------------------------------------------------------

    resid.ops <- NULL
    .do.resid <- FALSE
    nuisance.col.names <- NULL

    if (isTRUE(x = .residualize)) {

      if (is.null(x = .tdr.obj@config$key)) {
        stop("@config$key not found. Cannot expand design to landmark level.")
      }

      Z.nuisance <-
        .build.nuisance.Z(.tdr.obj = .tdr.obj,
                          .model.name = .model.name,
                          .coef.col = .coef.col,
                          .verbose = .verbose)

      if (!is.null(x = Z.nuisance)) {

        # Preserve original column names before .build.resid.operators,
        # which may QR-reduce Z and lose them
        nuisance.col.names <- colnames(x = Z.nuisance)

        if (isTRUE(x = .verbose)) {
          message("Residualizing expression against ",
                  ncol(x = Z.nuisance),
                  " nuisance column(s): ",
                  paste(nuisance.col.names, collapse = ", "))
        }

        resid.ops <-
          .build.resid.operators(Z = Z.nuisance,
                                X.sparse = X,
                                muX = muX,
                                .verbose = .verbose)

        .do.resid <- TRUE

      } else {
        if (isTRUE(x = .verbose)) {
          message("No nuisance covariates found; skipping residualization.")
        }
      }
    }

    # Effective column means (residualized or original)
    muX.eff <-
      if (.do.resid) resid.ops$muX.resid else muX

    # -------------------------------------------------------------------------
    # Build random-walk normalized graph P
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Building random-walk normalized graph...")
    }
    
    SNN <-
      .tdr.obj@graphs$snn
    
    P.graph <-
      .build.P(SNN = SNN,
               .degree.reg = .degree.reg,
               .tau.mult = .tau.mult,
               .lazy.alpha = .lazy.alpha,
               .verbose = .verbose)
    
    # -------------------------------------------------------------------------
    # Implicit M.local operators  (no dense n x p matrix ever formed)
    #
    # M_local = columnCenter( P.graph  diag(Y)  (X - 1 muX') )   (YX.interaction)
    # M_local = columnCenter( P.graph  (X - 1 muX') )             (no interaction)
    #
    # We precompute muM (column means of M_local) from vectors alone,
    # then define M_mv(w) = M_local %*% w   and
    #             Mt_mv(v) = t(M_local) %*% v
    # touching X and P.graph only through sparse-dense matvec products.
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      if (isTRUE(x = .YX.interaction)) {
        message("Setting up sparse-implicit M.local operators (P diag(Y) Xc)...")
      } else {
        message("Setting up sparse-implicit M.local operators (P Xc)...")
      }
    }
    
    # --- precompute column means of M.local (all O(nnz) or O(n)) ---
    #
    # muM_j = (1/n) * [P.graph' 1]' diag(Y?) (X - 1 muX')
    #       = (1/n) * { (qw)' X - sum(qw) muX }
    #
    # where  q = P.graph' 1  (column sums of P.graph, length n)
    #        qw = q * Y      (if YX.interaction)  or  q  (otherwise)
    
    q.col <-
      as.numeric(
        x = Matrix::crossprod(
          x = P.graph,
          y = rep(x = 1, times = n.landmarks)
        )
      )
    
    if (isTRUE(x = .YX.interaction)) {
      qw <- q.col * Y
    } else {
      qw <- q.col
    }
    
    muM <-
      (1 / n.landmarks) *
      (as.numeric(x = Matrix::crossprod(x = X, y = qw)) -
         sum(qw) * muX)

    # Correct muM for residualization: X -> X - Z B
    if (.do.resid) {
      Zt.qw <-
        as.numeric(x = crossprod(x = resid.ops$Z, y = qw))   # length q_z
      muM <- muM -
        (1 / n.landmarks) *
        (as.numeric(x = crossprod(x = resid.ops$B, y = Zt.qw)) +
           sum(qw) * (muX.eff - muX))
    }
    
    # --- M_mv(w): M_local %*% w  (returns dense length n) ---
    #
    # = P.graph  diag(Y?)  (X w - (muX'w) 1)  -  (muM'w) 1
    #
    # Cost: one sparse X %*% w  +  one sparse P.graph %*% vec  =  O(nnz(X) + nnz(P))
    
    M_mv <- function(w) {
      w <- as.numeric(x = w)
      # sparse X %*% dense w  ->  dense length n
      xw <- as.numeric(x = X %*% w)
      # residualization correction: (X - Z B) w = X w - Z(Bw)
      if (.do.resid) {
        xw <- xw -
          as.numeric(x = resid.ops$Z %*% (resid.ops$B %*% w))
      }
      # implicit centering: Xc w = Xtilde w - (muX.eff'w) 1
      xcw <- xw - sum(muX.eff * w)
      # density weighting (element-wise, length n)
      if (isTRUE(x = .YX.interaction)) xcw <- Y * xcw
      # graph smoothing (sparse P.graph %*% dense vec)
      pxcw <- as.numeric(x = P.graph %*% xcw)
      # center M.local (subtract column-mean contribution)
      pxcw - sum(muM * w)
    }
    
    # --- Mt_mv(v): t(M_local) %*% v  (returns dense length p) ---
    #
    # = (X - 1 muX')' diag(Y?) P.graph' v  -  muM sum(v)
    # = X' (Y? * P.graph' v) - muX sum(Y? * P.graph' v)  -  muM sum(v)
    #
    # Cost: one sparse crossprod(X, vec)  +  one sparse crossprod(P.graph, vec)
    
    Mt_mv <- function(v) {
      v <- as.numeric(x = v)
      # sparse P.graph^T %*% dense v  ->  dense length n
      ptv <- as.numeric(
        x = Matrix::crossprod(x = P.graph, y = v)
      )
      # density weighting
      if (isTRUE(x = .YX.interaction)) ptv <- Y * ptv
      # sparse X^T %*% result  ->  dense length p
      xtv <- as.numeric(
        x = Matrix::crossprod(x = X, y = ptv)
      )
      # residualization correction: (X - Z B)' v = X' v - B' (Z' v)
      if (.do.resid) {
        xtv <- xtv -
          as.numeric(x = crossprod(x = resid.ops$B,
                                   y = crossprod(x = resid.ops$Z, y = ptv)))
      }
      # implicit centering of X (transpose side)
      xctv <- xtv - muX.eff * sum(ptv)
      # center M.local (transpose side)
      xctv - muM * sum(v)
    }
    
    # -------------------------------------------------------------------------
    # NIPALS PLS1 with implicit deflation
    #
    # Instead of  Z <- Z - t p'  (which densifies on step 1), we accumulate
    # T  (n x k)  and  Pload  (p x k)  and substitute:
    #   M_k  w = M_0 w - T_{<k} (Pload_{<k}' w)
    #   M_k' v = M_0'v - Pload_{<k} (T_{<k}' v)
    # Only Y is deflated in-place (a length-n vector).
    # -------------------------------------------------------------------------
    
    if (isTRUE(x = .verbose)) {
      message("Running NIPALS PLS1 (", .ncomp, " components, sparse-implicit)...")
    }
    
    eps <- .Machine$double.eps
    
    Y.work <-
      Y
    
    # Pre-allocate storage (small dense matrices)
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
    
    nfit <- 0L
    
    for (k in seq_len(length.out = .ncomp)) {
      
      # Gene weight:  w proportional to M_k' Y.work
      #   M_k' Y.work = M_0' Y.work - Pload_{<k} (T_{<k}' Y.work)
      cvec <-
        Mt_mv(v = Y.work)
      
      if (k > 1L) {
        cvec <- cvec -
          P.mat[, 1:(k - 1), drop = FALSE] %*%
          crossprod(
            x = T.mat[, 1:(k - 1), drop = FALSE],
            y = Y.work
          )
        cvec <- as.numeric(x = cvec)
      }
      
      w.norm <-
        sqrt(x = sum(cvec^2))
      
      if (!is.finite(x = w.norm) || w.norm < eps) {
        if (isTRUE(x = .verbose)) {
          message("  Component ", k, ": ||M_k' Y|| ~ 0, stopping early.")
        }
        break
      }
      
      w <-
        cvec / w.norm
      
      # Score:  t = M_k w = M_0 w - T_{<k} (Pload_{<k}' w)
      t.score <-
        M_mv(w = w)
      
      if (k > 1L) {
        t.score <- t.score -
          T.mat[, 1:(k - 1), drop = FALSE] %*%
          crossprod(
            x = P.mat[, 1:(k - 1), drop = FALSE],
            y = w
          )
        t.score <- as.numeric(x = t.score)
      }
      
      tt <-
        sum(t.score^2)
      
      if (!is.finite(x = tt) || tt < eps) {
        if (isTRUE(x = .verbose)) {
          message("  Component ", k, ": t't ~ 0, stopping early.")
        }
        break
      }
      
      # X-loading:  p = M_k' t / (t't)
      pvec <-
        Mt_mv(v = t.score)
      
      if (k > 1L) {
        pvec <- pvec -
          P.mat[, 1:(k - 1), drop = FALSE] %*%
          crossprod(
            x = T.mat[, 1:(k - 1), drop = FALSE],
            y = t.score
          )
        pvec <- as.numeric(x = pvec)
      }
      pvec <- pvec / tt
      
      # Y loading (scalar):  q = Y.work' t / (t't)
      q.load <-
        sum(Y.work * t.score) / tt
      
      # Deflate Y only (length-n vector, negligible memory)
      Y.work <-
        Y.work - t.score * q.load
      
      # Store
      W.mat[, k] <-
        w
      
      T.mat[, k] <-
        t.score
      
      P.mat[, k] <-
        pvec
      
      Q.vec[k] <-
        q.load
      
      nfit <- k
      
      if (isTRUE(x = .verbose)) {
        Ak <-
          abs(x = stats::cor(x = Y, y = t.score))
        message(sprintf("  plsD%d: cov = %.2f, Ak = %.4f, q = %.4f, resid.var(Y) = %.4f",
                        k,
                        abs(x = sum(t.score * Y)),
                        Ak,
                        q.load,
                        stats::var(x = Y.work) / stats::var(x = Y)))
      }
    }
    
    # Trim if stopped early
    .ncomp <- nfit
    if (.ncomp < ncol(x = W.mat)) {
      W.mat <- W.mat[, seq_len(length.out = .ncomp), drop = FALSE]
      T.mat <- T.mat[, seq_len(length.out = .ncomp), drop = FALSE]
      P.mat <- P.mat[, seq_len(length.out = .ncomp), drop = FALSE]
      Q.vec <- Q.vec[seq_len(length.out = .ncomp)]
    }
    
    if (.ncomp == 0L) {
      stop("plsD: no components could be extracted (M' Y ~ 0 on first step).")
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
      paste0("plsD", seq_len(length.out = .ncomp))
    
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
      rownames(x = X)
    
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
    
    # Compute posterior SD from limma eBayes slots
    .fit <- .tdr.obj@results$lm[[.model.name]]$fit
    post.sd.vec <-
      sqrt(x = .fit$s2.post) *
      .fit$stdev.unscaled[, .coef.col]

    .loadings.res <-
      .compute.loadings(
        scores      = T.mat,
        X.sparse    = X,
        .method     = .loading.method,
        Y0.vec      = coef.mat[, .coef.col],
        post.sd.vec = post.sd.vec,
        muX         = muX,
        resid.ops   = resid.ops
      )

    raw.loadings        <- .loadings.res$raw.loadings
    concordance.weights <- .loadings.res$concordance.weights
    loadings            <-
      if (all(is.na(x = concordance.weights))) {
        raw.loadings
      } else {
        raw.loadings * concordance.weights
      }
    
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
    
    if (is.null(x = .tdr.obj@results$pls)) {
      .tdr.obj@results$pls <-
        list()
    }
    
    .tdr.obj@results$pls[[.coef.col]] <-
      list(
        coord = T.mat,
        gene.weights = W.mat,
        x.loadings = P.mat,
        y.loadings = Q.vec,
        raw.loadings = raw.loadings,
        loadings = loadings,
        concordance.weights = concordance.weights,
        Y.alignment = Y.alignment,
        smoothness = smoothness,
        Y = Y,
        Y.mean = mean(x = coef.mat[, .coef.col]),
        params = list(
          model.name = .model.name,
          coef.col = .coef.col,
          ncomp = .ncomp,
          min.prop = .min.prop,
          degree.reg = .degree.reg,
          tau.mult = .tau.mult,
          lazy.alpha = .lazy.alpha,
          YX.interaction = .YX.interaction,
          loading.method = .loading.method,
          residualize = .residualize,
          nuisance.cols = nuisance.col.names
        )
      )
    
    if (isTRUE(x = .store.M)) {
      # Materialise M.local on demand (not kept in memory during PLS)
      if (isTRUE(x = .verbose)) {
        message("Materialising M.local for storage (.store.M = TRUE)...")
      }
      Xc.tmp <-
        X |>
        (\(x) {
          # Apply residualization if active
          if (.do.resid) {
            # Materialize Xtilde = X - Z B (this WILL be dense)
            x <- x - resid.ops$Z %*% resid.ops$B
          }
          Matrix::t(x = x) - muX.eff
        })() |>
        Matrix::t()
      if (isTRUE(x = .YX.interaction)) {
        YX.term <-
          Matrix::Diagonal(x = Y) %*% Xc.tmp
      } else {
        YX.term <- Xc.tmp
      }
      M.local <-
        (P.graph %*% YX.term) |>
        (\(x)
         Matrix::t(x = x) - Matrix::colMeans(x = x)
        )() |>
        Matrix::t()
      .tdr.obj@results$pls[[.coef.col]]$M.local <-
        M.local
      rm(Xc.tmp, YX.term, M.local)
    }
    
    if (isTRUE(x = .verbose)) {
      message("\nResults stored in: .tdr.obj$plsD$", .coef.col)
      message("  $scores       : ", n.landmarks, " landmarks x ", .ncomp, " components")
      message("  $gene.weights : ", n.genes, " features x ", .ncomp, " components (PLS w)")
          message("  $raw.loadings : ", n.genes, " features x ", .ncomp, " components (", .loading.method, ")")
          message("  $loadings     : ", n.genes, " features x ", .ncomp,
            " components (concordance-filtered; raw for spearman)")
      message("  $concordance.weights : ", n.genes, " features x ", .ncomp,
              " components (soft concordance, NA if spearman)")
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

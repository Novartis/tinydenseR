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

#' Initialize tinydenseR object for landmark-based analysis
#'
#' Creates and validates the main tinydenseR object structure that will hold
#' expression data, metadata, and analysis results throughout the workflow.
#' This is the required first step before landmark identification.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Validates input data structure and compatibility
#'   \item Calculates the number of landmarks to sample per sample (max 5000 total)
#'   \item Creates a "key" vector mapping landmarks to samples
#'   \item Initializes empty slots for downstream analyses (PCA, graph, mapping)
#'   \item Performs quality checks (warns if sample sizes vary >10-fold)
#' }
#' 
#' The landmark sampling strategy aims for proportional representation across samples
#' while capping total landmarks at 5000 for computational efficiency. Smaller samples
#' may have relatively more landmarks sampled to ensure adequate representation.
#' 
#' @param .cells A named list of file paths (character strings) pointing to RDS files,
#'   each containing one expression matrix per sample. For RNA: sparse matrix (dgCMatrix)
#'   with genes as rows and cells as columns. For cytometry: matrix with cells as rows and
#' markers as columns. List names become sample identifiers.
#' @param .meta A data.frame with sample-level metadata. Rownames must match names in
#'   \code{.cells} exactly. Used for batch correction and downstream modeling.
#' @param .markers Character vector of marker names to use for landmark identification.
#'   Only applicable for \code{.assay.type = "cyto"}. If NULL, uses all markers from
#'   the first sample. Ignored for RNA data (uses HVG selection instead).
#' @param .harmony.var Character vector of column names from \code{.meta} to use for
#'   Harmony batch correction. If NULL, no batch correction is performed.
#' @param .assay.type Character string: "cyto" for cytometry (default) or "RNA" for
#'   scRNA-seq. Determines normalization strategy and feature selection approach.
#' @param .verbose Logical, whether to print progress messages. Default TRUE.
#' @param .seed Integer for random seed to ensure reproducibility. Default 123.
#' @param .prop.landmarks Numeric between 0 and 1 specifying proportion of cells to
#'   sample as landmarks. Default 0.1 (10%). Total landmarks capped at about 5000 regardless.
#' @param .n.threads Integer for parallel processing. Default automatically detects
#'   maximum available threads (using BLAS settings on HPC, or \code{detectCores()} locally).
#' 
#' @return A list object with initialized structure containing:
#'   \describe{
#'     \item{\code{$cells}}{Input file paths}
#'     \item{\code{$metadata}}{Sample metadata with added cell count columns}
#'     \item{\code{$key}}{Named vector mapping future landmarks to samples}
#'     \item{\code{$spec}}{Specifications including \code{n.cells}, \code{n.perSample}}
#'     \item{\code{$assay.type}}{Assay type ("cyto" or "RNA")}
#'     \item{\code{$markers}}{Marker list (cytometry only)}
#'     \item{\code{$harmony.var}}{Batch variables (if provided)}
#'     \item{Empty slots}{lm, scaled.lm, raw.lm, pca, graph, map, etc. - populated by downstream functions}
#'   }
#' 
#' @examples
#' \dontrun{
#' # Prepare data (from README workflow)
#' .meta <- get.meta(.obj = sim_trajectory, .sample.var = "Sample")
#' .cells <- get.cells(.exprs = sim_trajectory, 
#'                     .meta = .meta,
#'                     .sample.var = "Sample")
#' 
#' # Basic setup for RNA data
#' lm.cells <- setup.lm.obj(
#'   .cells = .cells,
#'   .meta = .meta,
#'   .assay.type = "RNA",
#'   .prop.landmarks = 0.15
#' )
#' 
#' # Cytometry workflow with marker selection
#' lm.cells <- setup.lm.obj(
#'   .cells = .cells,
#'   .meta = .meta,
#'   .markers = c("CD3", "CD4", "CD8", "CD19"),
#'   .assay.type = "cyto"
#' )
#' 
#' # RNA workflow with batch correction
#' lm.cells <- setup.lm.obj(
#'   .cells = .cells,
#'   .meta = .meta,
#'   .harmony.var = "batch",
#'   .assay.type = "RNA"
#' )
#' }
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
        stop(".harmony.var must be a character vector.")
      }
      
      if(!all(.harmony.var %in% colnames(x = .meta))){
        stop("Variables not found in metadata: ",
             paste(.harmony.var[!(.harmony.var %in% colnames(x = .meta))],
                   collapse = ", "),
             "\nCheck column names in .meta with colnames(.meta).")
      }
      
    }
    
    if(!is.null(x = .markers)){
      if(.assay.type == "RNA"){
        stop(".markers argument only applies to cytometry data.\n",
             "For RNA data, feature selection uses highly variable genes (HVG) automatically.")
      } else if(length(x = .markers) < 3){
        stop(".markers must contain at least 3 markers for meaningful dimensionality reduction.")
      }
    }
    
    if(!inherits(x = .meta,
                 what = "data.frame")){
      stop(".meta must be a data.frame object")
    }
    
    if(!all(rownames(x = .meta) ==
            names(x = .cells))){
      stop("Sample names mismatch between .cells and .meta.\n",
           "names(.cells) must exactly match rownames(.meta) in the same order.")
    }
    
    if(!inherits(x = .cells,
                 what = "list")){
      stop(".cells must be a list object")
    }
    
    if((.prop.landmarks < 0) | (.prop.landmarks > 1)){
      stop(".prop.landmarks must be between 0 and 1 (e.g., 0.1 for 10% of cells).\n",
           "Current value: ", .prop.landmarks)
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
    
    # Quality check: warn if sample sizes are highly imbalanced (>10-fold difference)
    if((max(.lm.obj$spec$n.cells) / min(.lm.obj$spec$n.cells)) > 10){
      
      warning("Sample size imbalance detected: largest/smallest ratio > 10.\n",
              "Smallest sample has ", min(.lm.obj$spec$n.cells), " cells.\n",
              "Consider removing low-quality samples or downsampling large samples for more balanced representation.")
      
      if(any(.lm.obj$spec$n.cells < 1000)){
        warning("Large variation in sample sizes detected. For cytometry, samples with <1000 cells may be unreliable.")
      }
      
    } 
    
    # Calculate target number of landmarks: 10% of total cells, capped at 5000
    .lm.obj$spec$target.lm.n <-
      pmin(sum(.lm.obj$spec$n.cells) * .prop.landmarks,
           5e3)
    
    # Allocate landmarks per sample: proportional to sample size, but capped
    .lm.obj$spec$n.perSample <-
      pmin(ceiling(x = .lm.obj$spec$n.cells * .prop.landmarks),
           ceiling(x = .lm.obj$spec$target.lm.n / length(x = .lm.obj$cells)))
    
    # Create key vector: maps each future landmark to its sample
    # (will be used after landmark selection to link back to metadata)
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

#' Identify landmarks via leverage score sampling
#'
#' Selects representative "landmark" cells from the full dataset using a two-pass
#' leverage score sampling strategy. Landmarks capture the representative cells
#' while being computationally tractable for downstream graph construction and mapping.
#'
#' @details
#' **Two-pass landmark selection algorithm:**
#' 
#' \strong{Pass 1 - Initial sampling:}
#' \enumerate{
#'   \item For RNA: normalize, log-transform, select top HVGs per sample
#'   \item For cytometry: use specified markers
#'   \item Compute sample-specific PCA
#'   \item Calculate leverage scores (sum of squared PC loadings per cell)
#'   \item Sample landmarks proportionally to leverage scores
#' }
#' 
#' \strong{Pass 2 - Refinement:}
#' \enumerate{
#'   \item Pool landmarks from all samples
#'   \item Compute dataset-wide PCA on pooled landmarks
#'   \item Project ALL cells onto this shared PC space
#'   \item Recalculate leverage scores using shared PCA
#'   \item Resample landmarks with improved scores (final set)
#' }
#' 
#' This two-pass approach ensures landmarks are representative of global
#' (not just sample-specific) variation patterns. Leverage score sampling
#' prioritizes cells in high-variance regions while maintaining diversity.
#' 
#' **Optional Harmony integration:**
#' If \code{.harmony.var} was specified in \code{setup.lm.obj}, performs
#' batch correction on landmark PC embeddings. This creates a reference
#' object for mapping query cells in a batch-corrected space.
#' 
#' @param .lm.obj Object initialized with \code{setup.lm.obj}.
#' @param .verbose Logical, print progress messages. Default TRUE.
#' @param .seed Integer for reproducibility. Default 123.
#' @param .nHVG Integer, number of highly variable genes to select for RNA data.
#'   Default 5000. Ignored for cytometry. Higher values capture more variation
#'   but increase computation time.
#' @param .nPC Integer, number of principal components for dimensionality reduction.
#'   Default 30. Must be less than the number of cells in smallest sample.
#' @param .exc.vdj.mito.ribo.genes.from.hvg Logical, whether to exclude V(D)J
#'   recombination genes, mitochondrial genes (MT-), and ribosomal genes (RP)
#'   from HVG selection (RNA only). Default TRUE. Recommended to avoid
#'   technical/biological noise dominating variation.
#' @param .force.in Character vector of gene names to force into the feature set
#'   regardless of variance (RNA only). Useful for known markers. Default NULL.
#' 
#' @return Updated \code{.lm.obj} with populated fields:
#'   \describe{
#'     \item{\code{$raw.lm}}{Raw counts matrix for landmarks (landmarks × features)}
#'     \item{\code{$lm}}{Processed landmark expression on selected features (landmarks × features):
#'       \itemize{
#'         \item RNA: PCA-reconstructed denoised expression (log2-scale after library size normalization)
#'         \item Cytometry: Original marker values on selected markers
#'       }
#'     }
#'     \item{\code{$scaled.lm}}{Z-scored landmark expression (landmarks × features, for visualization/heatmaps)}
#'     \item{\code{$pca}}{List containing PCA results:
#'       \itemize{
#'         \item \code{$embed} - PC coordinates for landmarks (landmarks × PCs)
#'         \item \code{$rotation} - Feature loadings (features × PCs)
#'         \item \code{$center} - Feature means (length = # features)
#'         \item \code{$scale} - Feature standard deviations (length = # features)
#'         \item \code{$sdev} - Standard deviations of PCs (length = # PCs)
#'         \item \code{$HVG} - Selected feature names (character vector)
#'       }
#'       If Harmony used: \code{$embed} and \code{$rotation} are Harmony-corrected/approximated
#'     }
#'     \item{\code{$harmony.obj}}{Symphony reference object (if \code{.harmony.var} specified), used for batch-corrected mapping of query cells}
#'   }
#' 
#' @examples
#' \dontrun{
#' # Typical workflow (from README)
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500, .nPC = 3)
#' 
#' # RNA with more PCs and custom HVGs
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nPC = 50, .nHVG = 3000)
#' 
#' # Force specific markers into feature set
#' lm.cells <- get.landmarks(lm.cells, 
#'                           .force.in = c("CD3D", "CD4", "CD8A"))
#' }
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
      stop("Invalid .lm.obj: missing structure.\n",
           "Initialize with setup.lm.obj() first.")
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
      stop("Invalid .lm.obj structure.\n",
           "Ensure it was created with setup.lm.obj() and not manually modified.")
    }
    
    if(any(.nPC > min(.lm.obj$metadata$n.cells))){
      stop("Too many PCs requested: .nPC (", .nPC, ") exceeds smallest sample size (",
           min(.lm.obj$metadata$n.cells), ").\n",
           "Reduce .nPC or remove small samples.")
    }
    
    # PASS 1: Sample-wise landmark selection using local PCA
    # For each sample independently, compute PCA and select landmarks via leverage scores
    .lm.obj$raw.lm <-
      seq_along(along.with = .lm.obj$cells) |>
      stats::setNames(nm = names(x = .lm.obj$cells)) |>
      lapply(FUN = function(.cells.idx){
        
        mat <-
          readRDS(file = .lm.obj$cells[[.cells.idx]])
        
        if(.lm.obj$assay.type == "RNA"){
          
          if(!inherits(x = mat,
                       what = "dgCMatrix")){
            stop("RNA expression matrices must be sparse (dgCMatrix class).\n")
          }
          
          # Identify genes to exclude from HVG selection (VDJ, mitochondrial, ribosomal)
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
          
          # Normalize and log-transform RNA data
          mat <-
            Matrix::t(x = mat) |>
            (\(x)
             # Size factor normalization: scale each cell by its total counts
             x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))
            )()
          
          mat@x <-
            log2(x = mat@x + 1)
          
          # Select top HVGs by variance
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
        
        # Compute leverage scores: sum of squared PC loadings per cell
        # Cells with high leverage contribute most to defining PC space
        if(.lm.obj$assay.type == "RNA"){
          
          # For RNA: use efficient truncated SVD (irlba)
          lev.score <-
            irlba::irlba(A = mat,
                         nv = .nPC,
                         center = mat.mean,
                         scale = mat.sd)$u |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        } else {
          
          # For cytometry: scale then use full SVD (smaller feature space)
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
        
        # Sample landmarks proportionally to leverage scores
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
    
    # Compute dataset-wide HVG selection and PCA on pooled Pass 1 landmarks
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
            
            warning("Genes not found in expression data: ",
                    paste(.force.in[!(.force.in %in% colnames(x = .lm.obj$raw.lm))],
                          collapse = ", "))
            
          }
          
          # Ensure forced-in genes are not excluded
          vdj.mito.ribo <-
            vdj.mito.ribo[!vdj.mito.ribo %in% .force.in]
          
          # Get indices of forced-in genes for later merging with HVGs
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
      
      # Normalize and log-transform pooled landmarks
      .lm.obj$lm <-
        .lm.obj$raw.lm /
        (Matrix::rowSums(x = .lm.obj$raw.lm) /
           mean(x = Matrix::rowSums(x = .lm.obj$raw.lm)))
      
      .lm.obj$lm@x <-
        log2(x = .lm.obj$lm@x + 1)
      
      # Select HVGs across all pooled landmarks, merge with forced-in genes
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
    
    # Optional: Perform Harmony batch correction on landmark embeddings
    if(!is.null(x = .lm.obj$harmony.var)){
      
      set.seed(seed = .seed)
      # Run Harmony to correct for batch effects in PC space
      # nclust parameter: adaptive based on landmark count (min 20, max 100)
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
    
    # PASS 2: Refine landmark selection using dataset-wide PCA
    # Now that we have initial landmarks from all samples, compute global PCA
    # and reselect landmarks based on leverage scores in this shared space
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
        
        # Calculate leverage scores in shared PC space (Harmony-corrected if applicable)
        if(!is.null(x = .lm.obj$harmony.obj)){
          
          # Map cells to Harmony-corrected space using Symphony
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
          
          # Standardize features for SVD
          Y <-
            Matrix::t(x = mat[,colnames(x = .lm.obj$lm)]) |>
            (\(x)
             (x - .lm.obj$pca$center) /
               .lm.obj$pca$scale
            )() |>
            Matrix::t()
          
          # Compute approximate feature loadings in Harmony space via SVD
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
          
          # Compute leverage scores using Harmony-corrected loadings
          lev.score <-
            ((Y %*% res$v) %*%
               diag(x = 1 / res$d)) |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        } else {
          
          # Without Harmony: project cells onto global PCA and compute leverage scores
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
        
        # Sample landmarks proportionally to leverage scores
        lm.sample <-
          sample(x = nrow(x = mat),
                 size = .lm.obj$spec$n.perSample[[.cells.idx]],
                 replace = FALSE,
                 prob = lev.score)
        
        # Extract raw counts for sampled landmarks
        landmarks <-
          if(.lm.obj$assay.type == "RNA"){
            
            # For RNA: reload to get raw counts, add sample prefix to cell IDs
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
    
    # Recompute PCA on refined Pass 2 landmarks (final feature set and embedding)
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
    
    # Reconstruct landmark expression matrices from PCA (denoised versions)
    if(.lm.obj$assay.type == "RNA"){
      
      # Z-scored version: for visualization/heatmaps
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
      
      # Log-normalized version: back-transformed from PCA, de-noised
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
      
      # Approximate feature loadings for Harmony-corrected space
      # Since Harmony operates in PC space (not feature space), we need to
      # back-calculate approximate loadings via SVD of corrected_PCs vs residualized_features
      # This is an approximation: Harmony corrections are non-linear, so perfect
      # feature loadings don't exist. But this gives reasonable loadings for interpretation.
      
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

#' Detect if running on High-Performance Computing cluster (internal)
#' 
#' Checks environment variables to determine if code is executing on an HPC
#' system with a job scheduler (SLURM, LSF, PBS, SGE, OAR). Used to optimize
#' thread detection strategy.
#' 
#' @return Logical: TRUE if HPC scheduler detected, FALSE otherwise
#' @keywords internal
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

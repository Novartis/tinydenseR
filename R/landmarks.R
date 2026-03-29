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

# Internal helper for progress display (not exported)
# Combines single-line overwriting progress bar with stage headers
.show_progress <- function(current, total, stage_label = NULL, item_label = NULL, 

                           start_time = NULL, width = 30) {
  pct <- current / total
  filled <- floor(pct * width)
  bar <- paste0("[", strrep("=", filled), 
                if(filled < width) ">" else "", 
                strrep(" ", max(0, width - filled - 1)), "]")
  
  # Build progress line
  line <- sprintf("\r%s %5.1f%% (%d/%d)", bar, pct * 100, current, total)
  

  # Add item label if provided (truncate if too long)
  if(!is.null(item_label)) {
    max_label <- 40
    if(nchar(item_label) > max_label) {
      item_label <- paste0(substr(item_label, 1, max_label - 3), "...")
    }
    line <- paste0(line, " ", item_label)
  }
  
  # Pad to clear previous content
  line <- sprintf("%-80s", line)
  
  cat(line)
  utils::flush.console()
  
  # Print completion with timing
  if(current == total) {
    if(!is.null(start_time)) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat(sprintf("\r%s 100.0%% (%d/%d) done in %.1fs%s\n", 
                  paste0("[", strrep("=", width), "]"),
                  total, total, elapsed, strrep(" ", 30)))
    } else {
      cat("\n")
    }
  }
}

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
#' while capping total landmarks at 5000 for computational efficiency. Large samples are
#' capped to prevent domination and ensure adequate representation.
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
#' @param .celltype.vec Optional named character vector mapping cell IDs to cell type labels.
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
#'     \item{\code{$config}}{Configuration parameters:
#'       \describe{
#'         \item{\code{$key}}{Named vector mapping future landmarks to samples}
#'         \item{\code{$sampling}}{Sampling parameters including \code{n.cells}, \code{n.perSample}}
#'         \item{\code{$assay.type}}{Assay type ("cyto" or "RNA")}
#'         \item{\code{$markers}}{Marker list (cytometry only)}
#'         \item{\code{$n.threads}}{Number of threads for parallel processing}
#'       }
#'     }
#'     \item{\code{$integration}}{Integration/batch correction results:
#'       \describe{
#'         \item{\code{$harmony.var}}{Batch variables (if provided)}
#'         \item{\code{$harmony.obj}}{Symphony reference object (populated by \code{get.landmarks})}
#'       }
#'     }
#'     \item{Empty slots}{landmarks, scaled.landmarks, raw.landmarks, pca, graph, map, etc. - populated by downstream functions}
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
#' tdr.obj <- setup.tdr.obj(
#'   .cells = .cells,
#'   .meta = .meta,
#'   .assay.type = "RNA",
#'   .prop.landmarks = 0.15
#' )
#' 
#' # Cytometry workflow with marker selection
#' tdr.obj <- setup.tdr.obj(
#'   .cells = .cells,
#'   .meta = .meta,
#'   .markers = c("CD3", "CD4", "CD8", "CD19"),
#'   .assay.type = "cyto"
#' )
#' 
#' # RNA workflow with batch correction
#' tdr.obj <- setup.tdr.obj(
#'   .cells = .cells,
#'   .meta = .meta,
#'   .harmony.var = "batch",
#'   .assay.type = "RNA"
#' )
#' }
#' @import Matrix
#' @export
#'
setup.tdr.obj <-
  function(
    .cells,
    .meta,
    .markers = NULL,
    .harmony.var = NULL,
    .assay.type = "cyto",
    .celltype.vec = NULL,
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
    
    .tdr.obj <- TDRObj(
      config = list(
        key = NULL,
        sampling = NULL,
        assay.type = .assay.type,
        markers = NULL,
        n.threads = .n.threads
      ),
      integration = list(
        harmony.var = NULL,
        harmony.obj = NULL
      )
    )
    
    if(lapply(X = .cells,
              FUN = inherits,
              what = "character") |>
       unlist() |>
       all()){
      
      if(isTRUE(x = .verbose)){
        message("-> Reading ", length(.cells), " expression matrices...")
        .read_start <- Sys.time()
      }
      
      n.cells <-
        seq_along(along.with = .cells) |>
        stats::setNames(nm = names(x = .cells)) |>
        lapply(FUN = function(cell_elem){
          
          .cell <- readRDS(file = .cells[[cell_elem]])
          
          if(isTRUE(x = .verbose)){
            .show_progress(current = cell_elem, 
                           total = length(.cells),
                           item_label = names(x = .cells)[cell_elem],
                           start_time = .read_start)
          }
          
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
          
          
        }) |>
        unlist()
      
      .tdr.obj@cells <- .cells
      
    } else {
      
      stop(".cells has to be a named list of characters pointing to the location of the count matrices on disk.")
      
    }
    
    .tdr.obj@config$sampling$n.cells <-
      n.cells
    
    # Validate .celltype.vec if provided
    if (!is.null(.celltype.vec)) {
      if (!is.character(.celltype.vec) || is.null(names(.celltype.vec))) {
        stop(".celltype.vec must be a named character vector ",
             "(names = cell IDs, values = cell type labels).")
      }
      expected_n <- sum(n.cells)
      if (length(.celltype.vec) != expected_n) {
        stop(".celltype.vec length (", length(.celltype.vec),
             ") does not match total cell count (", expected_n, ").")
      }
    }
    
    # Quality check: warn if sample sizes are highly imbalanced (>10-fold difference)
    if((max(.tdr.obj@config$sampling$n.cells) / min(.tdr.obj@config$sampling$n.cells)) > 10){
      
      warning("Sample size imbalance detected: largest/smallest ratio > 10.\n",
              "Smallest sample has ", min(.tdr.obj@config$sampling$n.cells), " cells.\n",
              "Consider removing low-quality samples.")
      
      if(any(.tdr.obj@config$sampling$n.cells < 1000)){
        warning("Large variation in sample sizes detected. For cytometry, samples with <1000 cells may be unreliable.")
      }
      
    } 
    
    # Calculate target number of landmarks: 10% of total cells, capped at 5000
    .tdr.obj@config$sampling$target.lm.n <-
      pmin(sum(.tdr.obj@config$sampling$n.cells) * .prop.landmarks,
           5e3)
    
    # Allocate landmarks per sample: proportional to sample size, but capped
    .tdr.obj@config$sampling$n.perSample <-
      pmin(ceiling(x = .tdr.obj@config$sampling$n.cells * .prop.landmarks),
           ceiling(x = .tdr.obj@config$sampling$target.lm.n / length(x = .tdr.obj@cells)))
    
    # Create key vector: maps each future landmark to its sample
    # (will be used after landmark selection to link back to metadata)
    .tdr.obj@config$key <-
      seq_along(along.with = .tdr.obj@cells) |>
      rep(times = .tdr.obj@config$sampling$n.perSample) |>
      (\(x)
       stats::setNames(object = x,
                       nm = names(x = .tdr.obj@cells)[x])
      )()
    
    .tdr.obj@metadata <-
      .meta
    
    .tdr.obj@metadata$n.perSample <- 
      .tdr.obj@config$sampling$n.perSample
    
    .tdr.obj@metadata$n.cells <-
      .tdr.obj@config$sampling$n.cells
    
    .tdr.obj@metadata$log10.n.cells <-
      log10(x = .tdr.obj@config$sampling$n.cells)
    
    if(.assay.type == "cyto"){
      .tdr.obj@config$markers <- 
        if(is.null(x = .markers)) readRDS(file = .tdr.obj@cells[[1]]) |> colnames() else .markers
    }
    
    if(!is.null(x = .harmony.var)){
      .tdr.obj@integration$harmony.var <- 
        .harmony.var  
    }
    
    if (!is.null(.celltype.vec)) {
      .tdr.obj@config$celltype.vec <- .celltype.vec
    }
    
    return(.tdr.obj)
    
  }

#' @rdname setup.tdr.obj
#' @description `setup.lm.obj()` is deprecated; use `setup.tdr.obj()` instead.
#' @param ... Arguments passed to `setup.tdr.obj()`.
#' @export
setup.lm.obj <- function(...) {
  .Deprecated("setup.tdr.obj", 
              msg = "setup.lm.obj() is deprecated. Use setup.tdr.obj() instead.")
  setup.tdr.obj(...)
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
#' If \code{.harmony.var} was specified in \code{setup.tdr.obj}, performs
#' batch correction on landmark PC embeddings. This creates a reference
#' object for mapping query cells in a batch-corrected space.
#' 
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object initialized with \code{setup.tdr.obj}.
#' @param .source The raw data object for non-file backends. \code{NULL} (default) for 
#'   the files backend; otherwise a Seurat, SingleCellExperiment, or anndataR AnnData object. 
#'   Used by \code{.get_sample_matrix()} to retrieve per-sample expression matrices.
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
#' @return Updated \code{.tdr.obj} with populated fields:
#'   \describe{
#'     \item{\code{$raw.landmarks}}{Raw counts matrix for landmarks (landmarks × features)}
#'     \item{\code{$landmarks}}{Processed landmark expression on selected features (landmarks × features):
#'       \itemize{
#'         \item RNA: PCA-reconstructed denoised expression (log2-scale after library size normalization)
#'         \item Cytometry: Original marker values on selected markers
#'       }
#'     }
#'     \item{\code{$scaled.landmarks}}{Z-scored landmark expression (landmarks × features, for visualization/heatmaps)}
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
#'     \item{\code{$integration$harmony.obj}}{Symphony reference object (if \code{.harmony.var} specified), used for batch-corrected mapping of query cells}
#'   }
#' 
#' @examples
#' \dontrun{
#' # Typical workflow (from README)
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500, .nPC = 3)
#' 
#' # RNA with more PCs and custom HVGs
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nPC = 50, .nHVG = 3000)
#' 
#' # Force specific markers into feature set
#' lm.cells <- get.landmarks(lm.cells, 
#'                           .force.in = c("CD3D", "CD4", "CD8A"))
#' }
#' @param ... Additional arguments passed to methods.
#' @export
#'
get.landmarks <- function(x, ...) UseMethod("get.landmarks")

#' @rdname get.landmarks
#' @export
get.landmarks.TDRObj <-
  function(x,
           .source = NULL,
           .verbose = TRUE,
           .seed = 123,
           .nHVG = 5000,
           .nPC = 30,
           .exc.vdj.mito.ribo.genes.from.hvg = TRUE,
           .force.in = NULL,
           ...){
    .tdr.obj <- x
    
    set.seed(seed = .seed)
    
    if(!inherits(x = .tdr.obj, what = "TDRObj")){
      # Fallback: accept legacy lists with correct structure
      if(is.list(.tdr.obj) && !is.null(names(.tdr.obj))){
        expected_names <- c("cells","landmarks","scaled.landmarks","raw.landmarks",
                            "metadata","config","integration","pca","graph","map",
                            "specDE","pbDE","markerDE","interact.plot")  # legacy slot names
        if(identical(x = names(x = .tdr.obj), y = expected_names)){
          .tdr.obj <- methods::as(.tdr.obj, "TDRObj")
        } else {
          stop("Invalid .tdr.obj structure.\n",
               "Ensure it was created with setup.tdr.obj() and not manually modified.")
        }
      } else {
        stop("Invalid .tdr.obj: missing structure.\n",
             "Initialize with setup.tdr.obj() first.")
      }
    }
    
    if(any(.nPC > min(.tdr.obj@metadata$n.cells))){
      stop("Too many PCs requested: .nPC (", .nPC, ") exceeds smallest sample size (",
           min(.tdr.obj@metadata$n.cells), ").\n",
           "Reduce .nPC or remove small samples.")
    }
    
    # PASS 1: Sample-wise landmark selection using local PCA
    # For each sample independently, compute PCA and select landmarks via leverage scores
    if(isTRUE(x = .verbose)){
      message("-> Pass 1: Selecting landmarks from ", length(.tdr.obj@cells), " samples...")
      .pass1_start <- Sys.time()
    }
    
    .tdr.obj@assay$raw <-
      seq_along(along.with = .tdr.obj@cells) |>
      stats::setNames(nm = names(x = .tdr.obj@cells)) |>
      lapply(FUN = function(.cells.idx){
        
        mat <-
          .get_sample_matrix(.source, .tdr.obj, .cells.idx)
        
        if(.tdr.obj@config$assay.type == "RNA"){
          
          if(inherits(x = mat,
                       what = c("IterableMatrix", "DelayedMatrix"))){
            mat <-
              methods::as(object = mat,
                         Class = "dgCMatrix")
          } else if(!inherits(x = mat,
                       what = "dgCMatrix")){
            stop("RNA expression matrices must be a dgCMatrix, IterableMatrix, or DelayedMatrix.\n")
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
            mat[,.tdr.obj@config$markers]
          
        }
        
        mat.mean <-
          Matrix::colMeans(x = mat)
        
        mat.sd <-
          sparseMatrixStats::colSds(x = mat)
        
        # Compute leverage scores: sum of squared PC loadings per cell
        # Cells with high leverage contribute most to defining PC space
        if(.tdr.obj@config$assay.type == "RNA"){
          
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
                 size = .tdr.obj@config$sampling$n.perSample[[.cells.idx]],
                 replace = FALSE,
                 prob = lev.score)
        
        landmarks <-
          if(.tdr.obj@config$assay.type == "RNA"){
            
            .get_sample_matrix(.source, .tdr.obj, .cells.idx)[,lm.sample] |>
              methods::as(Class = "dgCMatrix") |>
              Matrix::t() |>
              (\(x)
               `rownames<-`(x = x,
                            value = paste0(names(x = .tdr.obj@cells)[.cells.idx],
                                           "_",
                                           rownames(x = x)))
              )()
            
          } else {
            .get_sample_matrix(.source, .tdr.obj, .cells.idx)[lm.sample,.tdr.obj@config$markers] |>
              (\(x)
               `rownames<-`(x = x,
                            value = paste0(names(x = .tdr.obj@cells)[.cells.idx],
                                           "_",
                                           rownames(x = x)))
              )()
          }
        
        if(isTRUE(x = .verbose)){
          .show_progress(current = .cells.idx, 
                         total = length(.tdr.obj@cells),
                         item_label = names(x = .tdr.obj@cells)[.cells.idx],
                         start_time = .pass1_start)
        }
        
        return(landmarks)
        
      }) |>
      do.call(what = rbind)
    
    if(.tdr.obj@config$assay.type == "RNA"){
      
      .tdr.obj@assay$raw <-
        methods::as(object = .tdr.obj@assay$raw,
                   Class = "dgCMatrix")
      
    }
    
    if(isTRUE(x = .verbose)){
      message("-> Computing landmark PCA...")
    }
    
    # Compute dataset-wide HVG selection and PCA on pooled Pass 1 landmarks
    if(.tdr.obj@config$assay.type == "RNA"){
      
      if(isTRUE(x = .exc.vdj.mito.ribo.genes.from.hvg)) {
        
        vdj.mito.ribo <-
          grep(pattern = "^TR[AB][VDJ]\\d+|^IG[KHL][VDJ]\\d+|^RP|^MT-",
               x = colnames(x = .tdr.obj@assay$raw),
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
          if(any(!.force.in %in% colnames(x = .tdr.obj@assay$raw))){
            
            warning("Genes not found in expression data: ",
                    paste(.force.in[!(.force.in %in% colnames(x = .tdr.obj@assay$raw))],
                          collapse = ", "))
            
          }
          
          # Ensure forced-in genes are not excluded
          vdj.mito.ribo <-
            vdj.mito.ribo[!vdj.mito.ribo %in% .force.in]
          
          # Get indices of forced-in genes for later merging with HVGs
          .force.in.idx <-
            which(x = colnames(x = .tdr.obj@assay$raw)[!colnames(x = .tdr.obj@assay$raw) %in% vdj.mito.ribo] %in% .force.in)
          
        } else {
          stop(".force.in must be either NULL or a character vector")
        }
      } else {
        .force.in.idx <-
          vector(mode = "integer",
                 length = 0)
      }
      
      # Normalize and log-transform pooled landmarks
      .tdr.obj@assay$expr <-
        .tdr.obj@assay$raw /
        (Matrix::rowSums(x = .tdr.obj@assay$raw) /
           mean(x = Matrix::rowSums(x = .tdr.obj@assay$raw)))
      
      .tdr.obj@assay$expr@x <-
        log2(x = .tdr.obj@assay$expr@x + 1)
      
      # Select HVGs across all pooled landmarks, merge with forced-in genes
      HVG <-
        sparseMatrixStats::colVars(
          x = .tdr.obj@assay$expr[,!colnames(x = .tdr.obj@assay$expr) %in% vdj.mito.ribo]) |>
        order(decreasing = TRUE) |>
        utils::head(n = .nHVG) |>
        (\(x)
         c(.force.in.idx,
           x) |>
           unique()
        )() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .tdr.obj@assay$expr)[!colnames(x = .tdr.obj@assay$expr) %in% vdj.mito.ribo][x])
        )()
      
    } else {
      
      .tdr.obj@assay$expr <-
        .tdr.obj@assay$raw
      
      HVG <-
        ncol(x = .tdr.obj@assay$expr) |>
        seq_len() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .tdr.obj@assay$expr)[x])
        )()
      
    }
    
    .tdr.obj@assay$expr <-
      if(.tdr.obj@config$assay.type == "RNA") {
        .tdr.obj@assay$expr[,colnames(x = .tdr.obj@assay$expr)[!colnames(x = .tdr.obj@assay$expr) %in% vdj.mito.ribo][HVG]]
      } else {
        .tdr.obj@assay$expr[,HVG]
      }
    
    pca_res1 <-
      vector(mode = "list")
    
    pca_res1$center <-
      Matrix::colMeans(x = .tdr.obj@assay$expr)
    
    pca_res1$scale <-
      sparseMatrixStats::colSds(x = .tdr.obj@assay$expr)
    
    if(.tdr.obj@config$assay.type == "RNA"){
      
      pca_res1 <-
        c(pca_res1,
          irlba::irlba(
          A = .tdr.obj@assay$expr,
          nv = .nPC,
          center = pca_res1$center,
          scale = pca_res1$scale))
      
    } else {
      
      pca_res1 <-
        c(pca_res1,
          Matrix::t(x = .tdr.obj@assay$expr) |>
        (\(x)
         (x - pca_res1$center) /
           pca_res1$scale
        )() |>
        Matrix::t() |>
        base::svd())
      
    }
    
    colnames(x = pca_res1$u) <- colnames(x = pca_res1$v) <- names(x = pca_res1$d) <-
      paste0("PC",
             seq_along(along.with = pca_res1$d))
    
    rownames(x = pca_res1$u) <-
      rownames(x = .tdr.obj@assay$expr)
    
    rownames(x = pca_res1$v) <-
      colnames(x = .tdr.obj@assay$expr)
    
    pca_res1$embed <-
      pca_res1$u %*%
      (base::diag(x = pca_res1$d) |>
         (\(x)
          `dimnames<-`(x= x,
                       value = list(names(x = pca_res1$d),
                                    names(x = pca_res1$d)))
         )())
    
#    .tdr.obj$pca <-
#      c(pca_res,
#        .tdr.obj$pca)
    
    # Optional: Perform Harmony batch correction on landmark embeddings
    if(!is.null(x = .tdr.obj@integration$harmony.var)){
      
      set.seed(seed = .seed)
      # Run Harmony to correct for batch effects in PC space
      # nclust parameter: adaptive based on landmark count (min 20, max 100)
      tmp.harmony <- 
        harmony::RunHarmony(data_mat = pca_res1$embed, 
                            meta_data = .tdr.obj@metadata[.tdr.obj@config$key,],
                            vars_use = .tdr.obj@integration$harmony.var, 
                            return_obj = TRUE,
                            # see https://github.com/immunogenomics/harmony/blob/b36bab002c1767af6e665c81f186b40a87870e64/R/ui.R#L191
                            nclust = min(round(x = nrow(x = pca_res1$embed) / 30), 100) |>
                              max(20),
                            verbose = .verbose)
      
      .tdr.obj@integration$harmony.obj <-
        symphony::buildReferenceFromHarmonyObj(
          harmony_obj = tmp.harmony,
          metadata = .tdr.obj@metadata[.tdr.obj@config$key,],
          vargenes_means_sds =
            dplyr::tibble(
              symbol = colnames(x = .tdr.obj@assay$expr),
              mean = pca_res1$center,
              stddev = pca_res1$scale
            ),
          pca_loadings = pca_res1$v,
          verbose = .verbose,
          do_umap = FALSE,
          save_uwot_path = NULL,
          #umap_min_dist = 0.1,
          seed = .seed)  
      
      if(isTRUE(x = .verbose)){
        message("Returning harmony corrected embedding.")
      }
      
      pca_res1$embed <-
        Matrix::t(x = .tdr.obj@integration$harmony.obj$Z_corr) |> 
        (\(x)
         `rownames<-`(x = x,
                      value = rownames(x = .tdr.obj@assay$expr))
         )()
      
    }
    
    # PASS 2: Refine landmark selection using dataset-wide PCA
    # Now that we have initial landmarks from all samples, compute global PCA
    # and reselect landmarks based on leverage scores in this shared space
    if(isTRUE(x = .verbose)){
      message("-> Pass 2: Refining landmarks from ", length(.tdr.obj@cells), " samples...")
      .pass2_start <- Sys.time()
    }
    
    .tdr.obj@assay$raw <-
      seq_along(along.with = .tdr.obj@cells) |>
      stats::setNames(nm = names(x = .tdr.obj@cells)) |> 
      lapply(FUN = function(.cells.idx){
        
        set.seed(seed = .seed)
        
        mat <-
          .get_sample_matrix(.source, .tdr.obj, .cells.idx)
        
        if(.tdr.obj@config$assay.type == "RNA"){
          
          mat <-
            methods::as(object = mat,
                       Class = "dgCMatrix") |>
            Matrix::t() |>
            (\(x)
             # size factor normalization, taking into consideration size factor of landmarks
             (x / (Matrix::rowSums(x = x) /
                     mean(x = Matrix::rowSums(x = x)))) * 
               (mean(x = Matrix::rowSums(x = x)) / mean(x = Matrix::rowSums(x = .tdr.obj@assay$raw)))
            )()
          
          mat@x <-
            log2(x = mat@x + 1)
          
        }
        
        # Calculate leverage scores in shared PC space (Harmony-corrected if applicable)
        if(!is.null(x = .tdr.obj@integration$harmony.obj)){
          
          # Map cells to Harmony-corrected space using Symphony
          corr_embed <-
            symphony::mapQuery(exp_query = Matrix::t(x = mat[,colnames(x = .tdr.obj@assay$expr)]),
                               metadata_query = matrix(data = 1,
                                                       nrow = nrow(x = mat),
                                                       ncol = 2),
                               ref_obj = .tdr.obj@integration$harmony.obj,
                               vars = NULL,
                               verbose = .verbose,
                               do_normalize = FALSE,
                               do_umap = FALSE)$Z
          
          # Standardize features for SVD
          Y <-
            Matrix::t(x = mat[,colnames(x = .tdr.obj@assay$expr)]) |>
            (\(x)
             (x - pca_res1$center) /
               pca_res1$scale
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
            list(colnames(x = .tdr.obj@assay$expr),
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
            (((((Matrix::t(x = mat[,colnames(x = .tdr.obj@assay$expr)]) - pca_res1$center) /
                  pca_res1$scale) |>
                 Matrix::t()) %*%
                pca_res1$v) %*%
               diag(x = 1 / pca_res1$d)) |>
            (\(x)
             Matrix::rowSums(x = x ^ 2)
            )()
          
        }
        
        # Sample landmarks proportionally to leverage scores
        lm.sample <-
          sample(x = nrow(x = mat),
                 size = .tdr.obj@config$sampling$n.perSample[[.cells.idx]],
                 replace = FALSE,
                 prob = lev.score)
        
        # Extract raw counts for sampled landmarks
        landmarks <-
          if(.tdr.obj@config$assay.type == "RNA"){
            
            # For RNA: reload to get raw counts, add sample prefix to cell IDs
            .get_sample_matrix(.source, .tdr.obj, .cells.idx)[,lm.sample] |>
              methods::as(Class = "dgCMatrix") |>
              Matrix::t() |>
              (\(x)
               `rownames<-`(x = x,
                            value = paste0(names(x = .tdr.obj@cells)[.cells.idx],
                                           "_",
                                           rownames(x = x)))
              )()
            
          } else {
            .get_sample_matrix(.source, .tdr.obj, .cells.idx)[lm.sample,.tdr.obj@config$markers] |>
              (\(x)
               `rownames<-`(x = x,
                            value = paste0(names(x = .tdr.obj@cells)[.cells.idx],
                                           "_",
                                           rownames(x = x)))
              )()
          }
        
        if(isTRUE(x = .verbose)){
          .show_progress(current = .cells.idx, 
                         total = length(.tdr.obj@cells),
                         item_label = names(x = .tdr.obj@cells)[.cells.idx],
                         start_time = .pass2_start)
        }
        
        return(landmarks)
        
      }) |>
      do.call(what = rbind)
    
    if(isTRUE(x = .verbose)){
      message("-> Recomputing landmark PCA...")
    }
    
    # Recompute PCA on refined Pass 2 landmarks (final feature set and embedding)
    if(.tdr.obj@config$assay.type == "RNA"){
      
      .tdr.obj@assay$expr <-
        .tdr.obj@assay$raw /
        (Matrix::rowSums(x = .tdr.obj@assay$raw) /
           mean(x = Matrix::rowSums(x = .tdr.obj@assay$raw)))
      
      .tdr.obj@assay$expr@x <-
        log2(x = .tdr.obj@assay$expr@x + 1)
      
      HVG <-
        sparseMatrixStats::colVars(
          x = .tdr.obj@assay$expr[,!colnames(x = .tdr.obj@assay$expr) %in% vdj.mito.ribo]) |>
        order(decreasing = TRUE) |>
        utils::head(n = .nHVG) |>
        (\(x)
         c(.force.in.idx,
           x) |>
           unique()
        )() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .tdr.obj@assay$expr)[!colnames(x = .tdr.obj@assay$expr) %in% vdj.mito.ribo][x])
        )()
      
    } else {
      
      .tdr.obj@assay$expr <-
        .tdr.obj@assay$raw
      
      HVG <-
        ncol(x = .tdr.obj@assay$expr) |>
        seq_len() |>
        (\(x)
         `names<-`(x = x,
                   value = colnames(x = .tdr.obj@assay$expr)[x])
        )()
      
    }
    
    .tdr.obj@assay$expr <-
      if(.tdr.obj@config$assay.type == "RNA") {
        .tdr.obj@assay$expr[,colnames(x = .tdr.obj@assay$expr)[!colnames(x = .tdr.obj@assay$expr) %in% vdj.mito.ribo][HVG]]
      } else {
        .tdr.obj@assay$expr[,HVG]
      }
    
    pca_res2 <-
      vector(mode = "list")
    
    pca_res2$center <-
      Matrix::colMeans(x = .tdr.obj@assay$expr)
    
    pca_res2$scale <-
      sparseMatrixStats::colSds(x = .tdr.obj@assay$expr)
    
    if(.tdr.obj@config$assay.type == "RNA"){
      
      pca_res2 <-
        c(pca_res2,
          irlba::irlba(
          A = .tdr.obj@assay$expr,
          nv = .nPC,
          center = pca_res2$center,
          scale = pca_res2$scale))
      
    } else {
      
      pca_res2 <-
        c(pca_res2,
          Matrix::t(x = .tdr.obj@assay$expr) |>
        (\(x)
         (x - pca_res2$center) /
           pca_res2$scale
        )() |>
        Matrix::t() |>
        base::svd())
      
    }
    
    #.tdr.obj$pca <-
    #  c(pca_res,
    #    .tdr.obj$pca)
    
    dimnames(x = pca_res2$u) <-
      list(rownames(x = .tdr.obj@assay$expr),
           paste0("PC",
                  seq_along(along.with = pca_res2$d)))
    
    names(x = pca_res2$d) <-
      paste0("PC",
             seq_along(along.with = pca_res2$d))
    
    dimnames(x = pca_res2$v) <-
      list(colnames(x = .tdr.obj@assay$expr),
           paste0("PC",
                  seq_along(along.with = pca_res2$d)))
    
    pca_res2$coord <-
      pca_res2$u %*%
      (base::diag(x = pca_res2$d) |>
         (\(x)
          `dimnames<-`(x= x,
                       value = list(names(x = pca_res2$d),
                                    names(x = pca_res2$d)))
         )())
    
    pca_res2$sdev <-
      pca_res2$d /
      sqrt(x = max(1, nrow(x = .tdr.obj@assay$expr) - 1))
    
    pca_res2$rotation <-
      pca_res2$v
    
    #pca_res2$v <- pca_res2$d <- pca_res2$u <-
    #  NULL
    
    pca_res2$HVG <-
      names(x = HVG)
    
    # Reconstruct landmark expression matrices from PCA (denoised versions)
    if(.tdr.obj@config$assay.type == "RNA"){
      
      # Z-scored version: for visualization/heatmaps
      .tdr.obj@assay$scaled <-
        (pca_res2$coord %*% Matrix::t(x = pca_res2$rotation)) |>
        Matrix::t() |>
        (\(x)
         (x - Matrix::rowMeans(x = x)) /
           matrixStats::rowSds(x = x)
        )() |>
        Matrix::t() |>
        (\(x)
         `dimnames<-`(x = x,
                      value = list(rownames(x = .tdr.obj@assay$expr),
                                   colnames(x = .tdr.obj@assay$expr)))
        )()
      
      # Log-normalized version: back-transformed from PCA, de-noised
      .tdr.obj@assay$expr <-
        (pca_res2$coord %*% Matrix::t(x = pca_res2$rotation)) |>
        Matrix::t() |>
        (\(x)
         (x * pca_res2$scale) +
           pca_res2$center)() |>
        Matrix::t() |>
        (\(x)
         `dimnames<-`(x = x,
                      value = list(rownames(x = .tdr.obj@assay$expr),
                                   colnames(x = .tdr.obj@assay$expr)))
        )()
      
    } else {
      
      .tdr.obj@assay$scaled <-
        Matrix::t(x = .tdr.obj@assay$expr) |>
        (\(x)
         (x - Matrix::rowMeans(x = x)) /
           matrixStats::rowSds(x = x)
        )() |>
        Matrix::t()
      
    }
    
    if(!is.null(x = .tdr.obj@integration$harmony.var)){
      
      warning("PCA embedding and rotations will be replaced by harmony-corrected embedding and roughly approximated rotations.")
      
      set.seed(seed = .seed)
      tmp.harmony <- 
        harmony::RunHarmony(data_mat = pca_res2$coord, 
                            meta_data = .tdr.obj@metadata[.tdr.obj@config$key,],
                            vars_use = .tdr.obj@integration$harmony.var, 
                            return_obj = TRUE,
                            # see https://github.com/immunogenomics/harmony/blob/b36bab002c1767af6e665c81f186b40a87870e64/R/ui.R#L191
                            nclust = min(round(x = nrow(x = pca_res2$coord) / 30), 100) |>
                              max(20),
                            verbose = .verbose)
      
      .tdr.obj@integration$harmony.obj <-
        symphony::buildReferenceFromHarmonyObj(
          harmony_obj = tmp.harmony,
          metadata = .tdr.obj@metadata[.tdr.obj@config$key,],
          vargenes_means_sds =
            dplyr::tibble(
              symbol = colnames(x = .tdr.obj@assay$expr),
              mean = pca_res2$center,
              stddev = pca_res2$scale
            ),
          pca_loadings = pca_res2$rotation,
          verbose = .verbose,
          do_umap = FALSE,
          save_uwot_path = NULL,
          #umap_min_dist = 0.1,
          seed = .seed)  
      
      if(isTRUE(x = .verbose)){
        message("Returning harmony corrected embedding.")
      }
      
      pca_res2$coord <-
        Matrix::t(x = .tdr.obj@integration$harmony.obj$Z_corr)  |> 
        (\(x)
         `rownames<-`(x = x,
                      value = rownames(x = .tdr.obj@assay$expr))
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
                            data = .tdr.obj@metadata[.tdr.obj@config$key, 
                                                    .tdr.obj@integration$harmony.var,
                                                    drop = FALSE])
      
      # Residualize each column in .tdr.obj$landmarks (unfortunately best can do is linear)
      Y.resid <- 
        stats::lm.fit(x = Z, # fast base R linear model fit
                      y = as.matrix(x = .tdr.obj@assay$expr)) |>  
        stats::residuals()  # same shape as Y
      
      # Scale by stored stddev
      Y.resid <- 
        Matrix::t(x = Y.resid) |> 
        (\(x)
         x / pca_res2$scale
        )() |>
        Matrix::t()
      
      Yt.X <-
        Matrix::crossprod(x = (.tdr.obj@integration$harmony.obj$Z_corr - Matrix::rowMeans(x = .tdr.obj@integration$harmony.obj$Z_corr)) |>
                            Matrix::t(),
                          y = Y.resid)
      
      res <-
        base::svd(x = Yt.X, 
                  nu = nrow(x = .tdr.obj@integration$harmony.obj$Z_corr),
                  nv = nrow(x = .tdr.obj@integration$harmony.obj$Z_corr))
      
      dimnames(x = res$v) <-
        list(colnames(x = .tdr.obj@assay$expr),
             rownames(x = .tdr.obj@integration$harmony.obj$Z_corr))
      
      res$v <-
        (Matrix::t(x = res$v) * 
           (res$u[(abs(x = res$u) |> 
                     matrixStats::rowRanks(ties.method = "average")) ==
                    ncol(x = res$u)] |>
              sign())) |>
        Matrix::t()
      
      pca_res2$rotation <-
        res$v
      
    }
    
    .tdr.obj@landmark.embed$pca <-
      pca_res2
    
    return(.tdr.obj)
    
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

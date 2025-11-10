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

# ============================================================================
# Data Preprocessing Functions
# ============================================================================
# This file contains functions for converting various single-cell data formats
# into standardized formats compatible with tinydenseR workflows

#' Create .cells Object from Count Matrices
#'
#' Converts a list of count matrices into tinydenseR's standardized .cells format. Each matrix 
#' represents one sample and is saved as a temporary RDS file for memory-efficient processing. 
#' Use this when you already have count data loaded in R memory.
#'
#' @param .count.mat.list Named list of count matrices, one per sample. Names must match 
#'   \code{rownames(.meta)}. Each matrix can be dense, sparse (dgCMatrix), or similar matrix-like 
#'   object with features as rows and cells as columns.
#' @param .meta Data frame with sample-level metadata. Rownames must match \code{names(.count.mat.list)}. 
#'   Should include experimental variables like Condition, Replicate, etc.
#' @param .compress Logical: compress RDS files? Default FALSE for faster I/O. Set TRUE to save disk 
#'   space for large datasets.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Named list of file paths to temporary RDS files, one per sample. Structure suitable for 
#'   \code{setup.lm.obj(.cells = ...)}.
#'   
#' @details
#' This function creates temporary RDS files containing each sample's count matrix. These files 
#' persist for the R session and allow tinydenseR to process data without loading all samples 
#' into memory simultaneously.
#' 
#' All matrices in \code{.count.mat.list} should have:
#' \itemize{
#'   \item Same feature names (rownames) across samples
#'   \item Cell barcodes as column names
#'   \item Raw or normalized counts (specify via \code{setup.lm.obj(.assay.type = ...)})
#' }
#' 
#' @seealso \code{\link{get.cells}} for automatic format detection, 
#'   \code{\link{get.cells.SCE}} for SingleCellExperiment input,
#'   \code{\link{setup.lm.obj}} for next step in workflow
#'   
#' @examples
#' \dontrun{
#' # Load example trajectory data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta
#' sim_trajectory <- trajectory_data$SCE
#'   
#' # Prepare sample-level metadata
#' sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#' 
#' # Extract count matrices for two samples
#' count.matrices <- list(
#'   A_R1 = SingleCellExperiment::counts(sim_trajectory)[, 
#'            sim_trajectory$Sample == "A_R1"],
#'   B_R1 = SingleCellExperiment::counts(sim_trajectory)[, 
#'            sim_trajectory$Sample == "B_R1"]
#' )
#' 
#' # Create .cells object
#' cells <- get.cells.list.mat(.count.mat.list = count.matrices,
#'                             .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ])
#' 
#' # Use in tinydenseR workflow
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                        .assay.type = "RNA")
#' }
#' 
#' @export
#'
get.cells.list.mat <-
  function(.count.mat.list,
           .meta,
           .compress = FALSE,
           .verbose = TRUE){
    
    # Validate input types
    if(!inherits(x = .count.mat.list,
                 what = "list")){
      stop(".count.mat.list must be a list object")
    }
    
    if(is.null(x = names(x = .count.mat.list))){
      stop("names of .count.mat.list cannot be NULL")
    }
    
    # Ensure list names match metadata rownames
    if(!all(rownames(x = .meta) ==
            names(x = .count.mat.list))){
      stop("names of .count.mat.list must be the same as rownames of .meta")
    }
    
    if(!inherits(x = .compress,
                 what = "logical")){
      stop(".compress must be a logical")
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    if(.verbose){
      cat("Creating .cells object from", length(x = .count.mat.list), "count matrices...\n")
    }
    
    # Create temporary RDS file for each sample matrix
    .cells <-
      .count.mat.list[rownames(x = .meta)] |>
      lapply(FUN = function(mat){
        
        # Generate unique temporary file path
        uri <-
          tempfile(fileext = ".RDS")
        
        # Save matrix to disk for memory-efficient processing
        saveRDS(object = mat,
                file = uri,
                compress = .compress)
        
        return(uri)
        
      })
    
    if(.verbose){
      cat("Successfully created .cells object with", length(x = .cells), "temporary RDS files\n")
    }
    
    return(.cells)
    
  }

#' Create .cells Object from Seurat v5 Object
#'
#' Converts a Seurat v5 object with layered data into tinydenseR's .cells format. Handles Seurat's 
#' v5 assay structure where data is split across multiple layers. Each sample's cells are extracted, 
#' combined across layers, and saved as temporary RDS files.
#'
#' @param .seurat.obj Seurat v5 object containing layered RNA data. Must use Assay5 class.
#' @param .meta Data frame with sample-level metadata. Rownames must be sample IDs matching values 
#'   in \code{.seurat.obj@meta.data[[.sample.var]]}. Required for filtering valid samples.
#' @param .sample.var Character: column name in \code{.seurat.obj@meta.data} identifying which sample 
#'   each cell belongs to (e.g., "sample_id", "Sample", "orig.ident").
#' @param .assay Character: assay name to extract. Default "RNA".
#' @param .layer.pattern Character: pattern to match layer names. Default "counts" matches layers 
#'   containing count data (e.g., "counts.Sample1", "counts.Sample2"). Uses \code{fixed = TRUE} matching.
#' @param .min.cells.per.sample Integer: minimum cells required per sample. Samples with fewer cells 
#'   are excluded. Default 10.
#' @param .compress Logical: compress RDS files? Default FALSE for faster I/O.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Named list of file paths to temporary RDS files, one per sample. Only includes samples 
#'   meeting minimum cell threshold and present in \code{.meta}.
#'   
#' @details
#' Seurat v5 uses a layered assay structure where data for different samples may be stored in 
#' separate layers. This function:
#' \enumerate{
#'   \item Identifies all layers matching \code{.layer.pattern}
#'   \item Finds common genes across all layers
#'   \item Extracts cells for each sample across matching layers
#'   \item Combines layer data into single sparse matrix per sample
#'   \item Saves as temporary RDS files for memory-efficient processing
#' }
#' 
#' Only samples appearing in both the Seurat object and \code{.meta} rownames are processed.
#' 
#' @seealso \code{\link{get.cells}} for automatic format detection,
#'   \code{\link{get.cells.Seurat}} for Seurat v4 objects,
#'   \code{\link{setup.lm.obj}} for next workflow step
#'   
#' @examples
#' \dontrun{
#' # Load example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta |>
#'   dplyr::select(Condition, Replicate, Sample) |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#'
#' # Create Seurat v5 object from SCE data
#' # Extract counts for two samples
#' counts.A_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
#'                  trajectory_data$SCE$Sample == "A_R1"]
#' counts.B_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
#'                  trajectory_data$SCE$Sample == "B_R1"]
#' combined.counts <- cbind(counts.A_R1, counts.B_R1)
#' 
#' # Build cell metadata
#' cell.meta <- data.frame(
#'   sample_id = c(rep("A_R1", ncol(counts.A_R1)),
#'                 rep("B_R1", ncol(counts.B_R1))),
#'   row.names = colnames(combined.counts)
#' )
#' 
#' # Create Seurat v5 object
#' seurat.obj <- CreateSeuratObject(counts = combined.counts,
#'                                  meta.data = cell.meta)
#' seurat.obj[["RNA"]] <- as(seurat.obj[["RNA"]], "Assay5")
#' seurat.obj <- JoinLayers(seurat.obj)
#' 
#' # Convert to .cells format
#' cells <- get.cells.Seurat5(.seurat.obj = seurat.obj,
#'                            .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                            .sample.var = "sample_id")
#' 
#' # Use in tinydenseR workflow
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                        .assay.type = "RNA")
#' }
#' 
#' @export
#'
get.cells.Seurat5 <-
  function(.seurat.obj,
           .meta,
           .sample.var,
           .assay = "RNA",
           .layer.pattern = "counts",
           .min.cells.per.sample = 10,
           .compress = FALSE,
           .verbose = TRUE){
    
    if(!inherits(x = .seurat.obj,
                 what = "Seurat")){
      stop(".seurat.obj must be a Seurat object")
    }
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!(.sample.var %in% colnames(x = .seurat.obj@meta.data))){
      stop(paste0(".sample.var '", .sample.var, "' not found in Seurat object metadata"))
    }
    
    if(!(.assay %in% names(x = .seurat.obj@assays))){
      stop(paste0("Assay '", .assay, "' not found in Seurat object"))
    }
    
    if(!inherits(x = .layer.pattern,
                 what = "character")){
      stop(".layer.pattern must be a character string")
    }
    
    if(!inherits(x = .min.cells.per.sample,
                 what = "numeric") ||
       .min.cells.per.sample < 1){
      stop(".min.cells.per.sample must be a positive integer")
    }
    
    if(!inherits(x = .compress,
                 what = "logical")){
      stop(".compress must be a logical")
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    # Get available layers matching the pattern
    .available.layers <-
      names(x = .seurat.obj@assays[[.assay]]@layers) |>
      grep(pattern = .layer.pattern,
           value = TRUE,
           fixed = TRUE)
    
    if(length(x = .available.layers) == 0){
      stop(paste0("No layers found matching pattern '", .layer.pattern, "' in assay '", .assay, "'"))
    }
    
    if(.verbose){
      cat("Found", length(x = .available.layers), "layers matching pattern:", paste0(.available.layers, collapse = ", "), "\n")
    }
    
    # Get intersect of genes across all layers
    .genes.names <-
      .available.layers |>
      lapply(FUN = function(layer){
        
        SeuratObject::LayerData(object = .seurat.obj,
                                assay = .assay,
                                layer = layer) |>
          rownames()
        
      }) |>
      Reduce(f = intersect)
    
    if(.verbose){
      cat("Using", length(x = .genes.names), "genes common across all layers\n")
    }
    
    # Get unique samples
    .unique.samples <-
      unique(x = .seurat.obj@meta.data[[.sample.var]])
    
    # Filter samples: must have minimum cells AND appear in .meta
    .sample.cell.counts <-
      table(.seurat.obj@meta.data[[.sample.var]])
    
    .valid.samples <-
      names(x = .sample.cell.counts)[.sample.cell.counts >= .min.cells.per.sample] |> 
      (\(x)
       x[x %in% rownames(x = .meta)]
       )()
    
    if(.verbose){
      cat("Found", length(x = .unique.samples), "total samples,", length(x = .valid.samples), "with >=", .min.cells.per.sample, "cells\n")
    }
    
    if(length(x = .valid.samples) == 0){
      stop("No samples meet the minimum cell count requirement")
    }
    
    # Create .cells object
    .cells <-
      .valid.samples |>
      (\(x)
       `names<-`(x = x, value = x)
      )() |>
      lapply(FUN = function(smpl){
        
        if(.verbose){
          .progress <-
            (which(x = .valid.samples == smpl) * 100 / length(x = .valid.samples)) |>
            round(digits = 2)
          cat("Processing sample", smpl, "(", .progress, "%)\n")
        }
        
        # Create temporary file for this sample
        uri <-
          tempfile(fileext = ".RDS")
        
        # Get cells belonging to this sample
        .cells.in.smpl <-
          .seurat.obj@meta.data[[.sample.var]] == smpl
        
        if(sum(.cells.in.smpl) == 0){
          return(NULL)
        }
        
        # Extract and combine count data from all layers for this sample
        .sample.counts <-
          .available.layers |>
          lapply(FUN = function(layer){
            
            # Extract the count matrix for the sample from this layer
            SeuratObject::LayerData(object = .seurat.obj,
                                    assay = .assay,
                                    layer = layer)[.genes.names,
                                                   .cells.in.smpl,
                                                   drop = FALSE]
            
          }) |>
          do.call(what = cbind)
        
        # Convert to dgCMatrix and save
        .sample.counts |>
          methods::as(Class = "dgCMatrix") |>
          (\(x)
           saveRDS(object = x,
                   file = uri,
                   compress = .compress)
          )()
        
        return(uri)
        
      }) |>
      (\(x)
       x[lengths(x = x) > 0]
      )()
    
    if(.verbose){
      cat("Successfully created .cells object with", length(x = .cells), "samples\n")
    }
    
    return(.cells)
    
  }

#' Create .cells Object from Seurat v4 Object
#'
#' Converts a Seurat v4 object into tinydenseR's .cells format. Unlike v5, Seurat v4 stores all 
#' data in a single slot per assay. Each sample's cells are extracted and saved as temporary RDS files.
#'
#' @param .seurat.obj Seurat v4 object containing expression data. Must use Assay class (not Assay5).
#' @param .meta Data frame with sample-level metadata. Rownames must be sample IDs matching values 
#'   in \code{.seurat.obj@meta.data[[.sample.var]]}. Required for filtering valid samples.
#' @param .sample.var Character: column name in \code{.seurat.obj@meta.data} identifying which sample 
#'   each cell belongs to (e.g., "sample_id", "Sample", "orig.ident").
#' @param .assay Character: assay name to extract. Default "RNA".
#' @param .slot Character: data slot to extract. Default "counts". Can also be "data" (normalized) 
#'   or "scale.data" (scaled).
#' @param .min.cells.per.sample Integer: minimum cells required per sample. Default 10.
#' @param .compress Logical: compress RDS files? Default FALSE for faster I/O.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Named list of file paths to temporary RDS files, one per sample. Only includes samples 
#'   meeting minimum cell threshold and present in \code{.meta}.
#'   
#' @details
#' Seurat v4 uses a simpler assay structure than v5, with all data stored in a single matrix per slot. 
#' This function extracts cells sample-by-sample and converts to sparse dgCMatrix format for 
#' memory-efficient storage.
#' 
#' Compatible with Seurat versions < 5.0.0 or when using legacy Assay class in v5.
#' 
#' @seealso \code{\link{get.cells}} for automatic format detection,
#'   \code{\link{get.cells.Seurat5}} for Seurat v5 objects,
#'   \code{\link{setup.lm.obj}} for next workflow step
#'   
#' @examples
#' \dontrun{
#' # Load example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta |>
#'   dplyr::select(Condition, Replicate, Sample) |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#'
#' # Create Seurat v4 object
#' counts.A_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
#'                  trajectory_data$SCE$Sample == "A_R1"]
#' counts.B_R1 <- SingleCellExperiment::counts(trajectory_data$SCE)[, 
#'                  trajectory_data$SCE$Sample == "B_R1"]
#' combined.counts <- cbind(counts.A_R1, counts.B_R1)
#' 
#' cell.meta <- data.frame(
#'   sample_id = c(rep("A_R1", ncol(counts.A_R1)),
#'                 rep("B_R1", ncol(counts.B_R1))),
#'   row.names = colnames(combined.counts)
#' )
#' 
#' seurat.obj <- CreateSeuratObject(counts = combined.counts,
#'                                  meta.data = cell.meta)
#' 
#' # Convert to .cells format
#' cells <- get.cells.Seurat(.seurat.obj = seurat.obj,
#'                           .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                           .sample.var = "sample_id")
#' 
#' # Use in tinydenseR workflow
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                        .assay.type = "RNA")
#' }
#' 
#' @export
#'
get.cells.Seurat <-
  function(.seurat.obj,
           .meta,
           .sample.var,
           .assay = "RNA",
           .slot = "counts",
           .min.cells.per.sample = 10,
           .compress = FALSE,
           .verbose = TRUE){
    
    if(!inherits(x = .seurat.obj,
                 what = "Seurat")){
      stop(".seurat.obj must be a Seurat object")
    }
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!(.sample.var %in% colnames(x = .seurat.obj@meta.data))){
      stop(paste0(".sample.var '", .sample.var, "' not found in Seurat object metadata"))
    }
    
    if(!(.assay %in% names(x = .seurat.obj@assays))){
      stop(paste0("Assay '", .assay, "' not found in Seurat object"))
    }
    
    if(!inherits(x = .slot,
                 what = "character")){
      stop(".slot must be a character string")
    }
    
    if(!inherits(x = .min.cells.per.sample,
                 what = "numeric") ||
       .min.cells.per.sample < 1){
      stop(".min.cells.per.sample must be a positive integer")
    }
    
    if(!inherits(x = .compress,
                 what = "logical")){
      stop(".compress must be a logical")
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    # Extract count matrix using GetAssayData for v4 compatibility
    .count.matrix <-
      SeuratObject::GetAssayData(object = .seurat.obj,
                                 assay = .assay,
                                 slot = .slot)
    
    # Extract metadata  
    .meta.data <-
      .seurat.obj@meta.data
    
    if(.verbose){
      cat("Extracting data from Seurat v4 object with", ncol(x = .count.matrix), "cells and", nrow(x = .count.matrix), "features\n")
    }
    
    # Get unique samples
    .unique.samples <-
      unique(x = .meta.data[[.sample.var]])
    
    # Filter samples by minimum cell count
    .sample.cell.counts <-
      table(.meta.data[[.sample.var]])
    
    .valid.samples <-
      names(x = .sample.cell.counts)[.sample.cell.counts >= .min.cells.per.sample] |> 
      (\(x)
       x[x %in% rownames(x = .meta)]
      )()
    
    if(.verbose){
      cat("Found", length(x = .unique.samples), "total samples,", length(x = .valid.samples), "with >=", .min.cells.per.sample, "cells\n")
    }
    
    if(length(x = .valid.samples) == 0){
      stop("No samples meet the minimum cell count requirement")
    }
    
    # Create .cells object
    .cells <-
      .valid.samples |>
      (\(x)
       `names<-`(x = x, value = x)
      )() |>
      lapply(FUN = function(smpl){
        
        if(.verbose){
          .progress <-
            (which(x = .valid.samples == smpl) * 100 / length(x = .valid.samples)) |>
            round(digits = 2)
          cat("Processing sample", smpl, "(", .progress, "%)\n")
        }
        
        # Create temporary file for this sample
        uri <-
          tempfile(fileext = ".RDS")
        
        # Get cells belonging to this sample
        .cells.in.smpl <-
          .meta.data[[.sample.var]] == smpl
        
        if(sum(.cells.in.smpl) == 0){
          return(NULL)
        }
        
        # Extract count matrix for this sample
        .sample.counts <-
          .count.matrix[, .cells.in.smpl, drop = FALSE]
        
        # Convert to dgCMatrix and save
        .sample.counts |>
          methods::as(Class = "dgCMatrix") |>
          (\(x)
           saveRDS(object = x,
                   file = uri,
                   compress = .compress)
          )()
        
        return(uri)
        
      }) |>
      (\(x)
       x[lengths(x = x) > 0]
      )()
    
    if(.verbose){
      cat("Successfully created .cells object with", length(x = .cells), "samples\n")
    }
    
    return(.cells)
    
  }

#' Create .cells Object from SingleCellExperiment Object
#'
#' Converts a SingleCellExperiment (SCE) object into tinydenseR's .cells format. SCE is the standard 
#' Bioconductor container for single-cell data. Each sample's cells are extracted and saved as 
#' temporary RDS files.
#'
#' @param .sce.obj SingleCellExperiment object containing expression data in assays.
#' @param .meta Data frame with sample-level metadata. Rownames must be sample IDs matching values 
#'   in \code{.sce.obj@colData[[.sample.var]]}. Required for filtering valid samples.
#' @param .sample.var Character: column name in \code{colData(.sce.obj)} identifying which sample 
#'   each cell belongs to (e.g., "sample_id", "Sample").
#' @param .assay Character: assay name to extract. Default "counts". Can also be "logcounts", 
#'   "normcounts", or any assay in \code{assayNames(.sce.obj)}.
#' @param .min.cells.per.sample Integer: minimum cells required per sample. Default 10.
#' @param .compress Logical: compress RDS files? Default FALSE for faster I/O.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Named list of file paths to temporary RDS files, one per sample. Only includes samples 
#'   meeting minimum cell threshold and present in \code{.meta}.
#'   
#' @details
#' SingleCellExperiment is the Bioconductor standard for single-cell data, used by many analysis 
#' packages. This function extracts cells sample-by-sample and converts to sparse dgCMatrix format 
#' for memory-efficient storage.
#' 
#' The function verifies that \code{colnames(.sce.obj)} exist, as these are required for proper 
#' cell tracking.
#' 
#' @seealso \code{\link{get.cells}} for automatic format detection,
#'   \code{\link{get.cells.Seurat}} for Seurat objects,
#'   \code{\link{setup.lm.obj}} for next workflow step
#'   
#' @examples
#' \dontrun{
#' # Load example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta |>
#'   dplyr::select(Condition, Replicate, Sample) |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#' sim_trajectory <- trajectory_data$SCE
#'
#' # Convert SCE directly to .cells format
#' cells <- get.cells.SCE(.sce.obj = sim_trajectory,
#'                        .meta = sim_trajectory.meta,
#'                        .sample.var = "Sample")
#'
#' # Use in tinydenseR workflow
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta,
#'                        .assay.type = "RNA")
#' }
#' 
#' @export
#'
get.cells.SCE <-
  function(.sce.obj,
           .meta,
           .sample.var,
           .assay = "counts",
           .min.cells.per.sample = 10,
           .compress = FALSE,
           .verbose = TRUE){
    
    if(!inherits(x = .sce.obj,
                 what = "SingleCellExperiment")){
      stop(".sce.obj must be a SingleCellExperiment object")
    }
    
    if(colnames(x = .sce.obj) |>
       is.null()){
      stop("colnames of .sce.obj cannot be NULL")
    }
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!(.sample.var %in% colnames(x = SummarizedExperiment::colData(.sce.obj)))){
      stop(paste0(".sample.var '", .sample.var, "' not found in SingleCellExperiment colData"))
    }
    
    if(!(.assay %in% SummarizedExperiment::assayNames(.sce.obj))){
      stop(paste0("Assay '", .assay, "' not found in SingleCellExperiment object"))
    }
    
    if(!inherits(x = .min.cells.per.sample,
                 what = "numeric") ||
       .min.cells.per.sample < 1){
      stop(".min.cells.per.sample must be a positive integer")
    }
    
    if(!inherits(x = .compress,
                 what = "logical")){
      stop(".compress must be a logical")
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    # Extract count matrix
    .count.matrix <-
      SummarizedExperiment::assay(x = .sce.obj, 
                                  i = .assay,
                                  withDimnames = TRUE)
    
    # Extract metadata  
    .coldata <-
      SummarizedExperiment::colData(.sce.obj) |>
      as.data.frame()
    
    if(.verbose){
      cat("Extracting data from SingleCellExperiment object with", ncol(x = .count.matrix), "cells and", nrow(x = .count.matrix), "features\n")
    }
    
    # Get unique samples
    .unique.samples <-
      unique(x = .coldata[[.sample.var]])
    
    # Filter samples by minimum cell count
    .sample.cell.counts <-
      table(.coldata[[.sample.var]])
    
    .valid.samples <-
      names(x = .sample.cell.counts)[.sample.cell.counts >= .min.cells.per.sample] |> 
      (\(x)
       x[x %in% rownames(x = .meta)]
      )()
    
    if(.verbose){
      cat("Found", length(x = .unique.samples), "total samples,", length(x = .valid.samples), "with >=", .min.cells.per.sample, "cells\n")
    }
    
    if(length(x = .valid.samples) == 0){
      stop("No samples meet the minimum cell count requirement")
    }
    
    # Create .cells object
    .cells <-
      .valid.samples |>
      (\(x)
       `names<-`(x = x, value = x)
      )() |>
      lapply(FUN = function(smpl){
        
        if(.verbose){
          .progress <-
            (which(x = .valid.samples == smpl) * 100 / length(x = .valid.samples)) |>
            round(digits = 2)
          cat("Processing sample", smpl, "(", .progress, "%)\n")
        }
        
        # Create temporary file for this sample
        uri <-
          tempfile(fileext = ".RDS")
        
        # Get cells belonging to this sample
        .cells.in.smpl <-
          .coldata[[.sample.var]] == smpl
        
        if(sum(.cells.in.smpl) == 0){
          return(NULL)
        }
        
        # Extract count matrix for this sample
        .sample.counts <-
          .count.matrix[, .cells.in.smpl, drop = FALSE]
        
        # Convert to dgCMatrix and save
        .sample.counts |>
          methods::as(Class = "dgCMatrix") |>
          (\(x)
           saveRDS(object = x,
                   file = uri,
                   compress = .compress)
          )()
        
        return(uri)
        
      }) |>
      (\(x)
       x[lengths(x = x) > 0]
      )()
    
    if(.verbose){
      cat("Successfully created .cells object with", length(x = .cells), "samples\n")
    }
    
    return(.cells)
    
  }

#' Create .cells Object with Automatic Format Detection
#'
#' Convenience wrapper that automatically detects input format and calls the appropriate conversion 
#' function. Supports Seurat (v4/v5), SingleCellExperiment, and list of count matrices. Simplifies 
#' workflow by eliminating need to know which specific function to call.
#'
#' @param .exprs Expression data in one of three formats:
#'   \itemize{
#'     \item Named list of count matrices (calls \code{get.cells.list.mat})
#'     \item Seurat object v4 or v5 (calls \code{get.cells.Seurat} or \code{get.cells.Seurat5})
#'     \item SingleCellExperiment object (calls \code{get.cells.SCE})
#'   }
#' @param .meta Data frame with sample-level metadata. Required for list input. Rownames must match 
#'   sample IDs. For Seurat/SCE, rownames must match values in \code{.sample.var} column.
#' @param .sample.var Character: metadata column identifying samples. Required for Seurat/SCE, 
#'   ignored for list input.
#' @param .assay Character: assay name. Default "RNA" for Seurat, "counts" for SCE. Ignored for 
#'   list input.
#' @param .layer.pattern Character: layer name pattern for Seurat v5. Default "counts". Ignored 
#'   for other formats.
#' @param .slot Character: data slot for Seurat v4. Default "counts". Ignored for other formats.
#' @param .min.cells.per.sample Integer: minimum cells per sample. Default 10. Ignored for list input.
#' @param .compress Logical: compress RDS files? Default FALSE.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Named list of file paths to temporary RDS files, one per sample.
#'   
#' @details
#' This function inspects the class of \code{.exprs} to determine format:
#' \itemize{
#'   \item \code{Seurat}: Checks for Assay5 class to distinguish v4 from v5
#'   \item \code{SingleCellExperiment}: Maps "RNA" assay to "counts" if needed
#'   \item \code{list}: Validates names match \code{.meta} rownames
#' }
#' 
#' For Seurat v5, automatically handles layered data structure. For SCE, uses standard 
#' \code{assay()} accessor.
#' 
#' @seealso \code{\link{get.cells.list.mat}}, \code{\link{get.cells.Seurat}}, 
#'   \code{\link{get.cells.Seurat5}}, \code{\link{get.cells.SCE}}, \code{\link{setup.lm.obj}}
#'   
#' @examples
#' \dontrun{
#' # Example 1: List of count matrices
#' count.matrices <- list(A_R1 = counts.A_R1, B_R1 = counts.B_R1)
#' cells <- get.cells(.exprs = count.matrices,
#'                    .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ])
#' 
#' # Example 2: Seurat object (auto-detects v4 vs v5)
#' cells <- get.cells(.exprs = seurat.obj,
#'                    .meta = sim_trajectory.meta,
#'                    .sample.var = "sample_id")
#' 
#' # Example 3: SingleCellExperiment object
#' cells <- get.cells(.exprs = sce.obj,
#'                    .meta = sim_trajectory.meta,
#'                    .sample.var = "Sample",
#'                    .assay = "counts")
#' }
#' 
#' @export
#'
get.cells <-
  function(.exprs,
           .meta,
           .sample.var = NULL,
           .assay = "RNA",
           .layer.pattern = "counts",
           .slot = "counts",
           .min.cells.per.sample = 10,
           .compress = FALSE,
           .verbose = TRUE){
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    # Detect input type and dispatch to appropriate function
    if(inherits(x = .exprs,
                what = "Seurat")){
      
      if(.verbose){
        cat("Detected Seurat object...\n")
      }
      
      if(is.null(x = .sample.var)){
        stop(".sample.var must be provided for Seurat input")
      }
      
      # Check if this is a Seurat v5 object by looking for Assay5 class
      .is.seurat5 <-
        inherits(x = .exprs@assays[[.assay]],
                 what = "Assay5")
      
      if(.is.seurat5){
        
        if(.verbose){
          cat("Using get.cells.Seurat5 for v5 object...\n")
        }
        
        return(get.cells.Seurat5(.seurat.obj = .exprs,
                                 .sample.var = .sample.var,
                                 .meta = .meta,
                                 .assay = .assay,
                                 .layer.pattern = .layer.pattern,
                                 .min.cells.per.sample = .min.cells.per.sample,
                                 .compress = .compress,
                                 .verbose = .verbose))
        
      } else {
        
        if(.verbose){
          cat("Using get.cells.Seurat for v4 object...\n")
        }
        
        return(get.cells.Seurat(.seurat.obj = .exprs,
                                .sample.var = .sample.var,
                                .meta = .meta,
                                .assay = .assay,
                                .slot = .slot,
                                .min.cells.per.sample = .min.cells.per.sample,
                                .compress = .compress,
                                .verbose = .verbose))
        
      }
      
    } else if(inherits(x = .exprs,
                       what = "SingleCellExperiment")){
      
      if(.verbose){
        cat("Detected SingleCellExperiment object, using get.cells.SCE...\n")
      }
      
      if(is.null(x = .sample.var)){
        stop(".sample.var must be provided for SingleCellExperiment input")
      }
      
      # For SCE, default assay should be "counts" if user didn't specify
      .sce.assay <-
        if(.assay == "RNA") "counts" else .assay
      
      return(get.cells.SCE(.sce.obj = .exprs,
                           .sample.var = .sample.var,
                           .meta = .meta,
                           .assay = .sce.assay,
                           .min.cells.per.sample = .min.cells.per.sample,
                           .compress = .compress,
                           .verbose = .verbose))
      
    } else if(inherits(x = .exprs,
                       what = "list")){
      
      if(.verbose){
        cat("Detected list of matrices, using get.cells.list.mat...\n")
      }
      
      if(is.null(x = .meta)){
        stop(".meta must be provided for list input")
      }
      
      return(get.cells.list.mat(.count.mat.list = .exprs,
                                .meta = .meta,
                                .compress = .compress,
                                .verbose = .verbose))
      
    } else {
      
      stop("Input type not supported. .exprs must be either a Seurat object, a SingleCellExperiment object, or a named list of count matrices")
      
    }
    
  }

# ============================================================================
# Metadata Extraction Functions
# ============================================================================

#' Extract Sample-Level Metadata from Seurat v5 Object
#'
#' Extracts sample-wise metadata from a Seurat v5 object by identifying which metadata columns have 
#' consistent values within each sample (sample-level variables). Cell-level variables are excluded 
#' with a warning. Returns a data frame suitable for \code{setup.lm.obj(.meta = ...)}.
#'
#' @param .seurat.obj Seurat v5 object with cell-level metadata in \code{@meta.data}.
#' @param .sample.var Character: column name in \code{.seurat.obj@meta.data} identifying samples.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Data frame with one row per sample, containing only sample-level metadata columns. 
#'   Rownames are sample IDs from \code{.sample.var}.
#'   
#' @details
#' This function automatically distinguishes sample-level from cell-level metadata:
#' \itemize{
#'   \item \strong{Sample-level}: Variables with identical values for all cells within each sample 
#'     (e.g., Condition, Batch, Treatment). These are kept.
#'   \item \strong{Cell-level}: Variables that vary between cells in the same sample (e.g., nCount_RNA, 
#'     percent.mt, cluster). These are excluded with a warning.
#' }
#' 
#' The \code{.sample.var} column is always included. Only distinct rows are returned (one per sample).
#' 
#' Useful when you've stored experimental design information in the Seurat object and want to 
#' extract it for tinydenseR analysis.
#' 
#' @seealso \code{\link{get.meta}} for automatic format detection, 
#'   \code{\link{get.meta.Seurat}} for Seurat v4,
#'   \code{\link{get.cells.Seurat5}} for extracting count data
#'   
#' @examples
#' \dontrun{
#' # Assuming Seurat object has sample metadata
#' meta <- get.meta.Seurat5(.seurat.obj = seurat.obj,
#'                          .sample.var = "sample_id")
#' 
#' # Use with get.cells
#' cells <- get.cells.Seurat5(.seurat.obj = seurat.obj,
#'                            .meta = meta,
#'                            .sample.var = "sample_id")
#' }
#' 
#' @export
#'
get.meta.Seurat5 <-
  function(.seurat.obj,
           .sample.var,
           .verbose = TRUE){
    
    if(!inherits(x = .seurat.obj,
                 what = "Seurat")){
      stop(".seurat.obj must be a Seurat object")
    }
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!(.sample.var %in% colnames(x = .seurat.obj@meta.data))){
      stop(paste0(".sample.var '", .sample.var, "' not found in Seurat object metadata"))
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    if(.verbose){
      cat("Extracting metadata from Seurat v5 object...\n")
    }
    
    # Identify sample-level columns (where all values within each sample are identical)
    .sample.level.cols <-
      .seurat.obj@meta.data |>
      dplyr::group_by(!!rlang::sym(.sample.var)) |>
      dplyr::summarize(dplyr::across(.cols = dplyr::everything(),
                                     .fns = ~ length(x = unique(x = .x)) == 1),
                       .groups = "drop") |> 
      (\(x)
       dplyr::select(.data = x[,!colnames(x = x) %in% .sample.var,
                               drop = FALSE],
                     dplyr::where(fn = all))
      )() |>
      colnames()
    
    # Always include the sample variable itself
    .sample.level.cols <-
      c(.sample.var, .sample.level.cols)
    
    # Identify excluded columns  
    .all.cols <-
      colnames(x = .seurat.obj@meta.data)
    
    .excluded.cols <-
      .all.cols[!.all.cols %in% .sample.level.cols]
    
    if(length(x = .excluded.cols) > 0 && .verbose){
      warning(paste0("Excluded ", length(x = .excluded.cols), 
                     " cell-level columns from .meta: ",
                     paste0(.excluded.cols, collapse = ", ")))
    }
    
    # Extract sample-wise metadata for sample-level columns only
    .sample.meta <-
      .seurat.obj@meta.data[, .sample.level.cols, drop = FALSE] |>
      dplyr::distinct() |>
      as.data.frame() |>
      (\(x)
       `rownames<-`(x = x, value = x[[.sample.var]])
      )()
    
    if(.verbose){
      cat("Extracted metadata for", nrow(x = .sample.meta), "samples with", ncol(x = .sample.meta), "sample-level variables\n")
    }
    
    return(.sample.meta)
    
  }

#' Extract Sample-Level Metadata from Seurat v4 Object
#'
#' Extracts sample-wise metadata from a Seurat v4 object by identifying which metadata columns have 
#' consistent values within each sample. Identical functionality to \code{get.meta.Seurat5} but for 
#' v4 objects. Returns a data frame suitable for \code{setup.lm.obj(.meta = ...)}.
#'
#' @param .seurat.obj Seurat v4 object with cell-level metadata in \code{@meta.data}.
#' @param .sample.var Character: column name in \code{.seurat.obj@meta.data} identifying samples.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Data frame with one row per sample, containing only sample-level metadata columns. 
#'   Rownames are sample IDs from \code{.sample.var}.
#'   
#' @details
#' Automatically distinguishes sample-level (consistent within samples) from cell-level (varying 
#' within samples) metadata. Only sample-level columns are retained. Cell-level columns are 
#' excluded with a warning if \code{.verbose = TRUE}.
#' 
#' The algorithm groups by \code{.sample.var} and tests each column for uniqueness within groups. 
#' Columns with multiple unique values within any sample are classified as cell-level.
#' 
#' @seealso \code{\link{get.meta}} for automatic format detection, 
#'   \code{\link{get.meta.Seurat5}} for Seurat v5,
#'   \code{\link{get.cells.Seurat}} for extracting count data
#'   
#' @examples
#' \dontrun{
#' # Extract sample metadata from Seurat v4 object
#' meta <- get.meta.Seurat(.seurat.obj = seurat.obj,
#'                         .sample.var = "sample_id")
#' 
#' # Use with get.cells
#' cells <- get.cells.Seurat(.seurat.obj = seurat.obj,
#'                           .meta = meta,
#'                           .sample.var = "sample_id")
#' }
#' 
#' @export
#'
get.meta.Seurat <-
  function(.seurat.obj,
           .sample.var,
           .verbose = TRUE){
    
    if(!inherits(x = .seurat.obj,
                 what = "Seurat")){
      stop(".seurat.obj must be a Seurat object")
    }
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!(.sample.var %in% colnames(x = .seurat.obj@meta.data))){
      stop(paste0(".sample.var '", .sample.var, "' not found in Seurat object metadata"))
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    if(.verbose){
      cat("Extracting metadata from Seurat v4 object...\n")
    }
    
    # Identify sample-level columns (where all values within each sample are identical)
    .sample.level.cols <-
      .seurat.obj@meta.data |>
      dplyr::group_by(!!rlang::sym(.sample.var)) |>
      dplyr::summarize(dplyr::across(.cols = dplyr::everything(),
                                     .fns = ~ length(x = unique(x = .x)) == 1),
                       .groups = "drop") |> 
      (\(x)
       dplyr::select(.data = x[,!colnames(x = x) %in% .sample.var,
                               drop = FALSE],
                     dplyr::where(fn = all))
      )() |>
      colnames()
    
    # Always include the sample variable itself
    .sample.level.cols <-
      c(.sample.var, .sample.level.cols)
    
    # Identify excluded columns  
    .all.cols <-
      colnames(x = .seurat.obj@meta.data)
    
    .excluded.cols <-
      .all.cols[!.all.cols %in% .sample.level.cols]
    
    if(length(x = .excluded.cols) > 0 && .verbose){
      warning(paste0("Excluded ", length(x = .excluded.cols), 
                     " cell-level columns from .meta: ",
                     paste0(.excluded.cols, collapse = ", ")))
    }
    
    # Extract sample-wise metadata for sample-level columns only
    .sample.meta <-
      .seurat.obj@meta.data[, .sample.level.cols, drop = FALSE] |>
      dplyr::distinct() |>
      as.data.frame() |>
      (\(x)
       `rownames<-`(x = x, value = x[[.sample.var]])
      )()
    
    if(.verbose){
      cat("Extracted metadata for", nrow(x = .sample.meta), "samples with", ncol(x = .sample.meta), "sample-level variables\n")
    }
    
    return(.sample.meta)
    
  }

#' Extract Sample-Level Metadata from SingleCellExperiment Object
#'
#' Extracts sample-wise metadata from a SingleCellExperiment object by identifying which colData 
#' columns have consistent values within each sample. Identical logic to Seurat metadata extraction 
#' but adapted for SCE's \code{colData} structure. Returns a data frame suitable for 
#' \code{setup.lm.obj(.meta = ...)}.
#'
#' @param .sce.obj SingleCellExperiment object with cell-level metadata in \code{colData}.
#' @param .sample.var Character: column name in \code{colData(.sce.obj)} identifying samples.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Data frame with one row per sample, containing only sample-level metadata columns. 
#'   Rownames are sample IDs from \code{.sample.var}.
#'   
#' @details
#' SingleCellExperiment stores cell-level metadata in \code{colData}. This function:
#' \enumerate{
#'   \item Converts \code{colData} to data frame
#'   \item Groups by \code{.sample.var}
#'   \item Identifies columns with unique values within each sample (sample-level)
#'   \item Excludes varying columns (cell-level) with warning
#'   \item Returns one row per sample
#' }
#' 
#' Common sample-level variables include experimental conditions, batches, donors. Common cell-level 
#' variables include QC metrics, clusters, cell types.
#' 
#' @seealso \code{\link{get.meta}} for automatic format detection, 
#'   \code{\link{get.meta.Seurat}} for Seurat objects,
#'   \code{\link{get.cells.SCE}} for extracting count data
#'   
#' @examples
#' \dontrun{
#' # Extract sample metadata from SCE object
#' meta <- get.meta.SCE(.sce.obj = sce.obj,
#'                      .sample.var = "sample_id")
#' 
#' # Use with get.cells
#' cells <- get.cells.SCE(.sce.obj = sce.obj,
#'                        .meta = meta,
#'                        .sample.var = "sample_id")
#' }
#' 
#' @export
#'
get.meta.SCE <-
  function(.sce.obj,
           .sample.var,
           .verbose = TRUE){
    
    if(!inherits(x = .sce.obj,
                 what = "SingleCellExperiment")){
      stop(".sce.obj must be a SingleCellExperiment object")
    }
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!(.sample.var %in% colnames(x = SummarizedExperiment::colData(.sce.obj)))){
      stop(paste0(".sample.var '", .sample.var, "' not found in SingleCellExperiment colData"))
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    if(.verbose){
      cat("Extracting metadata from SingleCellExperiment object...\n")
    }
    
    # Convert colData to data frame for analysis
    .coldata.df <-
      SummarizedExperiment::colData(.sce.obj) |>
      as.data.frame()
    
    # Identify sample-level columns (where all values within each sample are identical)
    .sample.level.cols <-
      .coldata.df |>
      dplyr::group_by(!!rlang::sym(.sample.var)) |>
      dplyr::summarize(dplyr::across(.cols = dplyr::everything(),
                                     .fns = ~ length(x = unique(x = .x)) == 1),
                       .groups = "drop") |> 
      (\(x)
       dplyr::select(.data = x[,!colnames(x = x) %in% .sample.var,
                               drop = FALSE],
                     dplyr::where(fn = all))
      )() |>
      colnames()
    
    # Always include the sample variable itself
    .sample.level.cols <-
      c(.sample.var, .sample.level.cols)
    
    # Identify excluded columns  
    .all.cols <-
      colnames(x = .coldata.df)
    
    .excluded.cols <-
      .all.cols[!.all.cols %in% .sample.level.cols]
    
    if(length(x = .excluded.cols) > 0 && .verbose){
      warning(paste0("Excluded ", length(x = .excluded.cols), 
                     " cell-level columns from .meta: ",
                     paste0(.excluded.cols, collapse = ", ")))
    }
    
    # Extract sample-wise metadata for sample-level columns only
    .sample.meta <-
      .coldata.df[, .sample.level.cols, drop = FALSE] |>
      dplyr::distinct() |>
      as.data.frame() |>
      (\(x)
       `rownames<-`(x = x, value = x[[.sample.var]])
      )()
    
    if(.verbose){
      cat("Extracted metadata for", nrow(x = .sample.meta), "samples with", ncol(x = .sample.meta), "sample-level variables\n")
    }
    
    return(.sample.meta)
    
  }

#' Extract Sample-Level Metadata with Automatic Format Detection
#'
#' Convenience wrapper that automatically detects input format (Seurat v4/v5 or SingleCellExperiment) 
#' and calls the appropriate metadata extraction function. Simplifies workflow by eliminating need 
#' to know object version.
#'
#' @param .obj Object containing cell-level metadata. Can be:
#'   \itemize{
#'     \item Seurat object (v4 or v5) - calls \code{get.meta.Seurat} or \code{get.meta.Seurat5}
#'     \item SingleCellExperiment object - calls \code{get.meta.SCE}
#'   }
#' @param .sample.var Character: column name in object metadata identifying samples.
#' @param .verbose Logical: print progress messages? Default TRUE.
#'   
#' @return Data frame with one row per sample, containing only sample-level metadata columns. 
#'   Rownames are sample IDs.
#'   
#' @details
#' This function inspects the class of \code{.obj} to determine format:
#' \itemize{
#'   \item \code{Seurat}: Checks first assay for Assay5 class to distinguish v4 from v5
#'   \item \code{SingleCellExperiment}: Uses standard \code{colData} accessor
#' }
#' 
#' Automatically filters to sample-level metadata by testing which columns have consistent values 
#' within each sample. Cell-level columns (varying within samples) are excluded with a warning.
#' 
#' @seealso \code{\link{get.meta.Seurat}}, \code{\link{get.meta.Seurat5}}, 
#'   \code{\link{get.meta.SCE}}, \code{\link{get.cells}} for extracting count data
#'   
#' @examples
#' \dontrun{
#' # Example 1: Seurat object (auto-detects v4 vs v5)
#' meta <- get.meta(.obj = seurat.obj,
#'                  .sample.var = "sample_id")
#' 
#' # Example 2: SingleCellExperiment object
#' meta <- get.meta(.obj = sce.obj,
#'                  .sample.var = "Sample")
#' 
#' # Use with get.cells
#' cells <- get.cells(.exprs = seurat.obj,
#'                    .meta = meta,
#'                    .sample.var = "sample_id")
#' }
#' 
#' @export
#'
get.meta <-
  function(.obj,
           .sample.var,
           .verbose = TRUE){
    
    if(!inherits(x = .sample.var,
                 what = "character") ||
       length(x = .sample.var) != 1){
      stop(".sample.var must be a single character string")
    }
    
    if(!inherits(x = .verbose,
                 what = "logical")){
      stop(".verbose must be a logical")
    }
    
    # Detect input type and dispatch to appropriate function
    if(inherits(x = .obj,
                what = "Seurat")){
      
      if(.verbose){
        cat("Detected Seurat object...\n")
      }
      
      # Check if this is a Seurat v5 object by looking for Assay5 class
      # Use the first available assay to check
      .first.assay <-
        names(x = .obj@assays)[1]
      
      .is.seurat5 <-
        inherits(x = .obj@assays[[.first.assay]],
                 what = "Assay5")
      
      if(.is.seurat5){
        
        if(.verbose){
          cat("Using get.meta.Seurat5 for v5 object...\n")
        }
        
        return(get.meta.Seurat5(.seurat.obj = .obj,
                                .sample.var = .sample.var,
                                .verbose = .verbose))
        
      } else {
        
        if(.verbose){
          cat("Using get.meta.Seurat for v4 object...\n")
        }
        
        return(get.meta.Seurat(.seurat.obj = .obj,
                               .sample.var = .sample.var,
                               .verbose = .verbose))
        
      }
      
    } else if(inherits(x = .obj,
                       what = "SingleCellExperiment")){
      
      if(.verbose){
        cat("Detected SingleCellExperiment object, using get.meta.SCE...\n")
      }
      
      return(get.meta.SCE(.sce.obj = .obj,
                          .sample.var = .sample.var,
                          .verbose = .verbose))
      
    } else {
      
      stop("Input type not supported. .obj must be either a Seurat object or a SingleCellExperiment object")
      
    }
    
  }
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

#' Create .cells object from count matrices
#'
#' This function facilitates the creation of a .cells object suitable for setup.lm.obj from count matrices already loaded into memory.
#'
#' @param .count.mat.list A named list containing count matrices for each sample. List names correspond to sample names and must match rownames of .meta. Each element can be a matrix, dgCMatrix, or similar matrix-like object.
#' @param .meta A data frame containing the sample-wise metadata. Rownames must match names of .count.mat.list.
#' @param .compress A logical specifying whether to compress the RDS files. Default is FALSE for faster I/O.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A named list containing the on-disk location of RDS files containing one expression matrix for each sample, suitable for use with setup.lm.obj.
#' @examples
#' \dontrun{
#' # Fetch example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta
#' sim_trajectory <- trajectory_data$SCE
#'   
#' # Process metadata for tinydenseR format
#' sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#' 
#' # Create .cells object using SCE method for subset
#' .min.meta <- sim_trajectory.meta[c("A_R1", "B_R1"), ]
#' .min.cells <- get.cells.SCE(.sce.obj = sim_trajectory,
#'                             .sample.var = "Sample") |>
#'   (\(x) x[.min.meta$Sample])()
#' 
#' # Create .cells object
#' cells <- get.cells.list.mat(.count.mat.list = count.matrices,
#'                    .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ])
#' 
#' # Use with setup.lm.obj
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                        .assay.type = "RNA")
#' }
#' @export
#'
get.cells.list.mat <-
  function(.count.mat.list,
           .meta,
           .compress = FALSE,
           .verbose = TRUE){
    
    if(!inherits(x = .count.mat.list,
                 what = "list")){
      stop(".count.mat.list must be a list object")
    }
    
    if(is.null(x = names(x = .count.mat.list))){
      stop("names of .count.mat.list cannot be NULL")
    }
    
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
    
    .cells <-
      .count.mat.list[rownames(x = .meta)] |>
      lapply(FUN = function(mat){
        
        # Create temporary file for this matrix
        uri <-
          tempfile(fileext = ".RDS")
        
        # Save matrix as RDS file
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

#' Create .cells object from Seurat v5 object
#'
#' This function facilitates the creation of a .cells object suitable for setup.lm.obj from a Seurat v5 object with layered data.
#'
#' @param .seurat.obj A Seurat v5 object containing layered RNA data.
#' @param .sample.var A character string specifying the metadata column to use for sample identification. This column will be used to split cells into samples.
#' @param .meta A data frame containing the sample-wise metadata. Row names must be sample names matching elements of the object metadata specified by `.sample.var`. Required for filtering valid samples.
#' @param .assay A character string specifying the assay to use. Default is "RNA".
#' @param .layer.pattern A character string specifying the pattern to match layer names. Default is "counts" to match layers containing count data.
#' @param .min.cells.per.sample An integer specifying the minimum number of cells required per sample. Samples with fewer cells will be excluded. Default is 10.
#' @param .compress A logical specifying whether to compress the RDS files. Default is FALSE for faster I/O.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A named list containing the on-disk location of RDS files containing one expression matrix for each sample, suitable for use with setup.lm.obj.
#' @examples
#' \dontrun{
#' # Fetch example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta
#' sim_trajectory <- trajectory_data$SCE
#'
#' # Process metadata for tinydenseR format
#' sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#'
#' # Use the SCE object directly for creating Seurat v5 object
#' # Extract count data from SCE for subset of samples
#' sim_trajectory.counts.A_R1 <- 
#'  SingleCellExperiment::counts(sim_trajectory)[, sim_trajectory@colData$Sample == "A_R1"]
#' sim_trajectory.counts.B_R1 <- 
#'  SingleCellExperiment::counts(sim_trajectory)[, sim_trajectory@colData$Sample == "B_R1"]
#'
#' # Create a Seurat v5 object from the example data
#' # First, combine the matrices
#' combined.counts <- cbind(sim_trajectory.counts.A_R1, sim_trajectory.counts.B_R1)
#' 
#' # Create cell metadata
#' cell.meta <- data.frame(
#'   sample_id = c(rep("A_R1", ncol(sim_trajectory.counts.A_R1)),
#'                 rep("B_R1", ncol(sim_trajectory.counts.B_R1))),
#'   row.names = colnames(combined.counts)
#' )
#' 
#' # Create Seurat object
#' seurat.obj <- CreateSeuratObject(counts = combined.counts,
#'                                  meta.data = cell.meta)
#' 
#' # Convert to v5 assay structure if needed
#' if(packageVersion("Seurat") >= "5.0.0"){
#'   seurat.obj[["RNA"]] <- as(seurat.obj[["RNA"]], "Assay5")
#'   seurat.obj <- JoinLayers(seurat.obj)
#' }
#' 
#' # Create .cells object
#' cells <- get.cells.Seurat5(.seurat.obj = seurat.obj,
#'                            .sample.var = "sample_id")
#' 
#' # Use with setup.lm.obj
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                        .assay.type = "RNA")
#' }
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
    
    # Filter samples by minimum cell count
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

#' Create .cells object from Seurat v4 object
#'
#' This function facilitates the creation of a .cells object suitable for setup.lm.obj from a Seurat v4 object.
#'
#' @param .seurat.obj A Seurat v4 object containing expression data.
#' @param .sample.var A character string specifying the metadata column to use for sample identification. This column will be used to split cells into samples.
#' @param .meta A data frame containing the sample-wise metadata. Row names must be sample names matching elements of the object metadata specified by `.sample.var`. Required for filtering valid samples.
#' @param .assay A character string specifying the assay to use. Default is "RNA".
#' @param .slot A character string specifying the data slot to use. Default is "counts".
#' @param .min.cells.per.sample An integer specifying the minimum number of cells required per sample. Samples with fewer cells will be excluded. Default is 10.
#' @param .compress A logical specifying whether to compress the RDS files. Default is FALSE for faster I/O.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A named list containing the on-disk location of RDS files containing one expression matrix for each sample, suitable for use with setup.lm.obj.
#' @examples
#' \dontrun{
#' # Fetch example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta
#' sim_trajectory <- trajectory_data$SCE
#'
#' # Process metadata for tinydenseR format
#' sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#'
#' # Extract count data from SCE for subset of samples
#' sim_trajectory.counts.A_R1 <-
#'  SingleCellExperiment::counts(sim_trajectory)[, sim_trajectory@colData$Sample == "A_R1"]
#' sim_trajectory.counts.B_R1 <-
#'  SingleCellExperiment::counts(sim_trajectory)[, sim_trajectory@colData$Sample == "B_R1"]
#'
#' # Create a Seurat v4 object from the example data
#' # First, combine the matrices
#' combined.counts <- cbind(sim_trajectory.counts.A_R1, sim_trajectory.counts.B_R1)
#' 
#' # Create cell metadata
#' cell.meta <- data.frame(
#'   sample_id = c(rep("A_R1", ncol(sim_trajectory.counts.A_R1)),
#'                 rep("B_R1", ncol(sim_trajectory.counts.B_R1))),
#'   row.names = colnames(combined.counts)
#' )
#' 
#' # Create Seurat object (will be v4 format by default in older Seurat versions)
#' seurat.obj <- CreateSeuratObject(counts = combined.counts,
#'                                  meta.data = cell.meta)
#' 
#' # Create .cells object
#' cells <- get.cells.Seurat(.seurat.obj = seurat.obj,
#'                           .sample.var = "sample_id")
#' 
#' # Use with setup.lm.obj
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ],
#'                        .assay.type = "RNA")
#' }
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

#' Create .cells object from SingleCellExperiment object
#'
#' This function facilitates the creation of a .cells object suitable for setup.lm.obj from a SingleCellExperiment object.
#'
#' @param .sce.obj A SingleCellExperiment object containing expression data.
#' @param .sample.var A character string specifying the metadata column to use for sample identification. This column will be used to split cells into samples.
#' @param .meta A data frame containing the sample-wise metadata. Row names must be sample names matching elements of the object metadata specified by `.sample.var`. Required for filtering valid samples.
#' @param .assay A character string specifying the assay to use. Default is "counts".
#' @param .min.cells.per.sample An integer specifying the minimum number of cells required per sample. Samples with fewer cells will be excluded. Default is 10.
#' @param .compress A logical specifying whether to compress the RDS files. Default is FALSE for faster I/O.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A named list containing the on-disk location of RDS files containing one expression matrix for each sample, suitable for use with setup.lm.obj.
#' @examples
#' \dontrun{
#' # Fetch example data
#' trajectory_data <- fetch_trajectory_data()
#' sim_trajectory.meta <- trajectory_data$meta
#' sim_trajectory <- trajectory_data$SCE
#'
#' # Process metadata for tinydenseR format
#' sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
#'   unique()
#' rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample
#'
#' # Use the SCE object directly
#' cells <- get.cells.SCE(.sce.obj = sim_trajectory,
#'                        .sample.var = "Sample")
#'
#' # Use with setup.lm.obj
#' lm.obj <- setup.lm.obj(.cells = cells,
#'                        .meta = sim_trajectory.meta,
#'                        .assay.type = "RNA")
#' }
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

#' Create .cells object automatically detecting input type
#'
#' This function automatically detects the input type and calls the appropriate function to create a .cells object suitable for setup.lm.obj.
#'
#' @param .exprs The expression data. Can be either:
#'   - A named list of count matrices (for get.cells.list.mat)
#'   - A Seurat object (for get.cells.Seurat or get.cells.Seurat5)
#'   - A SingleCellExperiment object (for get.cells.SCE)
#' @param .meta A data frame containing the sample-wise metadata. Row names must be sample names matching the names of .exprs if .exprs is a list or elements of the object metadata specified by `.sample.var` if .exprs is a Seurat/SCE object. Required for list input.
#' @param .sample.var A character string specifying the metadata column to use for sample identification. Required for Seurat/SCE input, ignored for list input.
#' @param .assay A character string specifying the assay to use. Default is "RNA" for Seurat, "counts" for SCE. Applies only to Seurat/SCE input.
#' @param .layer.pattern A character string specifying the pattern to match layer names. Default is "counts". Applies only to Seurat v5 input.
#' @param .slot A character string specifying the data slot to use. Default is "counts". Applies only to Seurat v4 input.
#' @param .min.cells.per.sample An integer specifying the minimum number of cells required per sample. Default is 10. Applies only to Seurat/SCE input.
#' @param .compress A logical specifying whether to compress the RDS files. Default is FALSE for faster I/O.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A named list containing the on-disk location of RDS files containing one expression matrix for each sample, suitable for use with setup.lm.obj.
#' @examples
#' \dontrun{
#' # Example 1: Using a named list of count matrices
#' count.matrices <- list(
#'   A_R1 = sim_trajectory.counts.A_R1,
#'   B_R1 = sim_trajectory.counts.B_R1
#' )
#' cells <- get.cells(.exprs = count.matrices,
#'                    .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ])
#' 
#' # Example 2: Using a Seurat object
#' cells <- get.cells(.exprs = seurat.obj,
#'                    .sample.var = "sample_id")
#' 
#' # Example 3: Using a SingleCellExperiment object
#' cells <- get.cells(.exprs = sce.obj,
#'                    .sample.var = "sample_id",
#'                    .assay = "counts")
#' }
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

#' Extract metadata from Seurat v5 object
#'
#' This function extracts sample-wise metadata from a Seurat v5 object.
#'
#' @param .seurat.obj A Seurat v5 object containing cell metadata.
#' @param .sample.var A character string specifying the metadata column to use for sample identification.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A data frame containing sample-wise metadata with sample names as rownames.
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat v5 object with sample metadata
#' meta <- get.meta.Seurat5(.seurat.obj = seurat.obj,
#'                          .sample.var = "sample_id")
#' }
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

#' Extract metadata from Seurat v4 object
#'
#' This function extracts sample-wise metadata from a Seurat v4 object.
#'
#' @param .seurat.obj A Seurat v4 object containing cell metadata.
#' @param .sample.var A character string specifying the metadata column to use for sample identification.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A data frame containing sample-wise metadata with sample names as rownames.
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat v4 object with sample metadata
#' meta <- get.meta.Seurat(.seurat.obj = seurat.obj,
#'                         .sample.var = "sample_id")
#' }
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

#' Extract metadata from SingleCellExperiment object
#'
#' This function extracts sample-wise metadata from a SingleCellExperiment object.
#'
#' @param .sce.obj A SingleCellExperiment object containing cell metadata.
#' @param .sample.var A character string specifying the metadata column to use for sample identification.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A data frame containing sample-wise metadata with sample names as rownames.
#' @examples
#' \dontrun{
#' # Assuming you have a SingleCellExperiment object with sample metadata
#' meta <- get.meta.SCE(.sce.obj = sce.obj,
#'                      .sample.var = "sample_id")
#' }
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

#' Extract metadata automatically detecting input type
#'
#' This function automatically detects the input type and calls the appropriate function to extract sample-wise metadata.
#'
#' @param .obj The object containing metadata. Can be either:
#'   - A Seurat object (for get.meta.Seurat or get.meta.Seurat5)
#'   - A SingleCellExperiment object (for get.meta.SCE)
#' @param .sample.var A character string specifying the metadata column to use for sample identification.
#' @param .verbose A logical specifying whether to print progress messages. Default is TRUE.
#' @return A data frame containing sample-wise metadata with sample names as rownames.
#' @examples
#' \dontrun{
#' # Example 1: Using a Seurat object
#' meta <- get.meta(.obj = seurat.obj,
#'                  .sample.var = "sample_id")
#' 
#' # Example 2: Using a SingleCellExperiment object
#' meta <- get.meta(.obj = sce.obj,
#'                  .sample.var = "sample_id")
#' }
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
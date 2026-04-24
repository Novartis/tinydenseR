#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

#' Internal accessor for per-sample expression matrices
#'
#' Dispatches to the appropriate backend to retrieve a single sample's
#' expression matrix.  For the default \code{"files"} backend the function
#' simply calls \code{\link[base]{readRDS}}; other backends slice into a
#' live in-memory object.
#'
#' @param .source  The raw data object.
#'   \code{NULL} for the files backend; otherwise a Seurat,
#'   SingleCellExperiment, or anndataR AnnData object.
#' @param .tdr.obj A \code{\linkS4class{TDRObj}}.
#' @param sample_idx Scalar index (integer or character name) into
#'   \code{.tdr.obj@@cells}.
#'
#' @return A matrix (dense or sparse) with the same orientation and class
#'   as the files backend would return.
#'
#' @keywords internal
#' @noRd
.get_sample_matrix <- function(.source, .tdr.obj, sample_idx) {

  backend <- .tdr.obj@config$backend
  if (is.null(backend)) backend <- "files"

  # For seurat/sce backends, fall back to stored source reference when

  # .source is NULL (e.g., child TDRObj from get.subset)
  if (is.null(.source) && backend %in% c("seurat", "sce")) {
    .source <- .tdr.obj@config$source.env$source
    if (is.null(.source)) {
      stop("No .source object available for '", backend, "' backend. ",
           "Pass the source object explicitly or ensure the child was ",
           "created via the container dispatch (get.subset.Seurat / ",
           "get.subset.SingleCellExperiment).",
           call. = FALSE)
    }
  }

  result <- switch(backend,
    "files" = {
      entry <- .tdr.obj@cells[[sample_idx]]
      if (is.list(entry)) {
        # Subset TDRObj: entry is list(path = ..., idx = ...)
        mat <- readRDS(file = entry$path)
        if (.tdr.obj@config$assay.type == "RNA") {
          mat[, entry$idx, drop = FALSE]
        } else {
          mat[entry$idx, , drop = FALSE]
        }
      } else {
        readRDS(file = entry)
      }
    },
    "seurat" = {
      col_idx <- .tdr.obj@cells[[sample_idx]]
      SeuratObject::LayerData(.source,
                              assay  = .tdr.obj@config$source.assay,
                              layer  = .tdr.obj@config$source.layer
                              )[, col_idx, drop = FALSE]
    },
    "sce" = {
      col_idx <- .tdr.obj@cells[[sample_idx]]
      SummarizedExperiment::assay(.source,
                                  .tdr.obj@config$source.assay
                                  )[, col_idx, drop = FALSE]
    },
    "matrix" = {
      col_idx <- .tdr.obj@cells[[sample_idx]]
      src <- .tdr.obj@config$source.env$mat
      if (.tdr.obj@config$assay.type == "RNA") {
        src[, col_idx, drop = FALSE]
      } else {
        src[col_idx, , drop = FALSE]
      }
    },
    "cyto" = {
      cs <- .tdr.obj@config$source.env$cs
      # Handle both integer index and character sample name
      if (is.numeric(sample_idx)) {
        sample_name <- names(.tdr.obj@cells)[sample_idx]
      } else {
        sample_name <- sample_idx
      }
      mat <- flowCore::exprs(cs[[sample_name]])
      # strip $PnN names attribute from colnames
      colnames(mat) <- unname(colnames(mat))
      if (is.null(rownames(mat))) {
        rownames(mat) <- paste0("event_", seq_len(nrow(mat)))
      }
      # Apply cell filter for subset TDRObjs
      cell_idx <- .tdr.obj@cells[[sample_idx]]
      if (length(cell_idx) < nrow(mat)) {
        mat <- mat[cell_idx, , drop = FALSE]
      }
      mat
    },
    stop("Unknown backend: ", backend)
  )

  # Orientation assertion for cyto data

  if (.tdr.obj@config$assay.type == "cyto" &&
      !is.null(.tdr.obj@config$markers)) {
    if (!all(.tdr.obj@config$markers %in% colnames(result))) {
      stop("Matrix from backend '", backend,
           "' does not contain expected marker columns.")
    }
  }

  return(result)
}

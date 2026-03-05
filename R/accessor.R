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

  result <- switch(backend,
    "files" = {
      readRDS(file = .tdr.obj@cells[[sample_idx]])
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
    "h5ad" = {
      col_idx <- .tdr.obj@cells[[sample_idx]]
      mat <- .source$X[col_idx, , drop = FALSE]
      if (.tdr.obj@config$assay.type == "RNA") {
        mat <- Matrix::t(mat)
      }
      if (!inherits(mat, "dgCMatrix") && methods::is(mat, "Matrix")) {
        mat <- methods::as(mat, "dgCMatrix")
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

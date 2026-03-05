#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

#' Run the full tinydenseR pipeline
#'
#' S3 generic that executes the complete tinydenseR analysis pipeline:
#' landmark selection, graph construction, optional cell typing, and mapping.
#'
#' @param x An object to run the pipeline on. See methods for supported types.
#' @param ... Additional arguments passed to methods.
#'
#' @return Depends on the method; see individual method documentation.
#'
#' @export
RunTDR <- function(x, ...) UseMethod("RunTDR")

#' @describeIn RunTDR Run the pipeline on an existing TDRObj
#'
#' Executes \code{get.landmarks → get.graph → celltyping → get.map}
#' on a pre-built \code{\linkS4class{TDRObj}}.
#'
#' @param x A \code{\linkS4class{TDRObj}} (created by \code{\link{setup.tdr.obj}}).
#' @param .celltype.vec Named character vector of per-cell labels
#'   (names = cell IDs, values = cell type labels). Length must equal
#'   \code{sum(x@@config$sampling$n.cells)}. If \code{NULL}, no cell
#'   typing is performed (unless one is already stored in the object).
#' @param .celltype.vec.overwrite Logical. If \code{TRUE}, replace an
#'   existing \code{.celltype.vec} stored in the TDRObj config.
#' @param .verbose Logical. Print progress messages.
#' @param .seed Integer. Random seed for reproducibility.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return The updated \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.default <- function(x,
                           .celltype.vec = NULL,
                           .celltype.vec.overwrite = FALSE,
                           .verbose = TRUE,
                           .seed = 123,
                           ...) {

  stopifnot(is.TDRObj(x))

  # --- .celltype.vec handling ---
  if (!is.null(.celltype.vec)) {

    if (!is.null(x@config$celltype.vec) && !isTRUE(.celltype.vec.overwrite)) {
      stop("A .celltype.vec already exists in TDRObj config. Set ",
           ".celltype.vec.overwrite = TRUE to replace it.")
    }

    if (!is.null(x@config$celltype.vec) && isTRUE(.celltype.vec.overwrite)) {
      x@config$celltype.vec <- .celltype.vec
    }

    if (is.null(x@config$celltype.vec)) {
      # Validate
      if (!is.character(.celltype.vec) || is.null(names(.celltype.vec))) {
        stop(".celltype.vec must be a named character vector ",
             "(names = cell IDs, values = cell type labels).")
      }
      expected_n <- sum(x@config$sampling$n.cells)
      if (length(.celltype.vec) != expected_n) {
        stop(".celltype.vec length (", length(.celltype.vec),
             ") does not match total cell count (", expected_n, ").")
      }
      x@config$celltype.vec <- .celltype.vec
    }

  }

  eff_ct <- x@config$celltype.vec

  # --- Pipeline ---
  x <- get.landmarks(x, .source = NULL, .seed = .seed,
                     .verbose = .verbose, ...)
  x <- get.graph(x, .seed = .seed, .verbose = .verbose, ...)

  if (!is.null(eff_ct)) {
    x <- celltyping(x, .celltyping.map = eff_ct, .verbose = .verbose)
  }

  x <- get.map(x, .source = NULL, .seed = .seed,
               .verbose = .verbose, ...)

  return(x)
}

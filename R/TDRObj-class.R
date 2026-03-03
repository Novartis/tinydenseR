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

#' TDRObj: S4 class for tinydenseR analysis objects
#'
#' The core data structure for tinydenseR landmark-based analysis. Contains
#' expression data references, metadata, dimensionality reduction results,
#' graph structure, and differential expression results.
#'
#' @slot cells list. Named list of file paths to expression matrices.
#' @slot landmarks ANY. Processed landmark expression matrix (NULL initially).
#' @slot scaled.landmarks ANY. Z-scored landmark expression matrix (NULL initially).
#' @slot raw.landmarks ANY. Raw counts matrix for landmarks (NULL initially).
#' @slot metadata data.frame. Sample-level metadata.
#' @slot config list. Configuration parameters (key, sampling, assay.type, markers, n.threads).
#' @slot integration list. Integration/batch correction results (harmony.var, harmony.obj).
#' @slot pca list. PCA results (embed, rotation, center, scale, sdev, HVG, u, v, d).
#' @slot graph list. Graph structure (uwot, adj.matrix, snn, LE, clustering, celltyping).
#' @slot map list. Mapping results (fdens, Y, clustering, celltyping, nearest.landmarks, fuzzy.graph, cl.ct.to.ign, lm).
#' @slot specDE ANY. Spectral DE results (NULL initially).
#' @slot pbDE ANY. Pseudo-bulk DE results (NULL initially).
#' @slot markerDE ANY. Marker DE results (NULL initially).
#' @slot interact.plot ANY. Interactive plot data (NULL initially).
#' @slot symphony.obj ANY. Symphony reference object for cell typing (NULL initially).
#' @slot nmfDE ANY. NMF-based DE results (NULL initially).
#' @slot plsDE ANY. PLS-based DE results (NULL initially).
#'
#' @name TDRObj-class
#' @rdname TDRObj-class
#' @exportClass TDRObj
#' @importFrom methods setClass setMethod setReplaceMethod setAs new slot slot<- slotNames show validObject is
setClass(
  Class = "TDRObj",
  slots = c(
    cells            = "list",
    landmarks        = "ANY",
    scaled.landmarks = "ANY",
    raw.landmarks    = "ANY",
    metadata         = "data.frame",
    config           = "list",
    integration      = "list",
    pca              = "list",
    graph            = "list",
    map              = "list",
    specDE           = "ANY",
    pbDE             = "ANY",
    markerDE         = "ANY",
    interact.plot    = "ANY",
    symphony.obj     = "ANY",
    nmfDE            = "ANY",
    plsDE            = "ANY"
  ),
  prototype = list(
    cells            = list(),
    landmarks        = NULL,
    scaled.landmarks = NULL,
    raw.landmarks    = NULL,
    metadata         = data.frame(),
    config           = list(),
    integration      = list(),
    pca              = list(),
    graph            = list(),
    map              = list(),
    specDE           = NULL,
    pbDE             = NULL,
    markerDE         = NULL,
    interact.plot    = NULL,
    symphony.obj     = NULL,
    nmfDE            = NULL,
    plsDE            = NULL
  ),
  validity = function(object) {
    errors <- character()
    if (!is.data.frame(object@metadata)) {
      errors <- c(errors, "metadata must be a data.frame")
    }
    if (!is.list(object@cells)) {
      errors <- c(errors, "cells must be a list")
    }
    if (!is.list(object@config)) {
      errors <- c(errors, "config must be a list")
    }
    if (length(object@config) > 0 && !("assay.type" %in% names(object@config))) {
      errors <- c(errors, "config must contain 'assay.type'")
    }
    if (length(errors) == 0) TRUE else errors
  }
)

#' Construct a TDRObj
#'
#' @param ... Named arguments corresponding to TDRObj slots.
#' @return A \code{TDRObj} object.
#' @export
TDRObj <- function(...) {
  new("TDRObj", ...)
}

#' Check if an object is a TDRObj
#'
#' @param x An object to test.
#' @return Logical: TRUE if \code{x} is a TDRObj, FALSE otherwise.
#' @export
is.TDRObj <- function(x) {
  inherits(x, "TDRObj")
}

# --- Backward-compatible $ and $<- methods ---
# These shims ensure existing code using .tdr.obj$slot continues to work
# during the migration from list-based to S4-based objects.

#' @rdname TDRObj-class
setMethod("$", "TDRObj", function(x, name) {
  slot(x, name)
})

#' @rdname TDRObj-class
setReplaceMethod("$", "TDRObj", function(x, name, value) {
  slot(x, name) <- value
  x
})

# --- names method for compatibility with names(.tdr.obj) calls ---

#' @rdname TDRObj-class
setMethod("names", "TDRObj", function(x) {
  slotNames(x)
})

# --- show method ---

#' @rdname TDRObj-class
setMethod("show", "TDRObj", function(object) {
  cat("TDRObj\n")
  cat("  Samples:", length(object@cells), "\n")
  backend <- if ("backend" %in% names(object@config)) object@config$backend else "files"
  assay <- if ("assay.type" %in% names(object@config)) object@config$assay.type else "unknown"
  cat("  Backend:", backend, "\n")
  cat("  Assay type:", assay, "\n")
  n_lm <- if (!is.null(object@landmarks)) nrow(object@landmarks) else 0L
  cat("  Landmarks:", n_lm, "\n")
  cat("  Graph computed:", length(object@graph) > 0, "\n")
  cat("  Map computed:", length(object@map) > 0, "\n")
})

# --- Coercion from list ---

setAs("list", "TDRObj", function(from) {
  # Only pass elements that are valid slot names
  valid <- intersect(names(from), slotNames("TDRObj"))
  do.call(new, c("TDRObj", from[valid]))
})

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
#' @slot cells list. Named list of per-sample file paths to expression matrices.
#' @slot metadata data.frame. Sample-level metadata.
#' @slot config list. Run parameters: key, sampling, assay.type, markers, n.threads.
#' @slot integration list. Trained projection models and batch variables
#'   (harmony.var, harmony.obj, symphony.obj, umap.model).
#' @slot assay list. Landmark expression layers (L x features matrices):
#'   raw (raw counts), expr (normalized/log expression), scaled (Z-scored).
#' @slot landmark.embed list. Landmark-space coordinate matrices; each entry
#'   has a $coord slot. Contains pca, le, and umap sub-lists.
#' @slot landmark.annot list. Per-landmark categorical annotations (factor, length L).
#'   Contains clustering and celltyping sub-lists, each with an $ids factor.
#' @slot graphs list. Landmark-landmark connectivity matrices: adj.matrix, snn, fgraph.
#' @slot density list. fdens-centric sample-level analytics (L x N matrices):
#'   fdens, Y, composition (clustering/celltyping cell counts/percentages),
#'   ignored, .cache.
#' @slot sample.embed list. Sample-level embeddings (N x k matrices), each with $coord.
#'   Contains pca, traj, and pepc sub-lists.
#' @slot cellmap list. Per-cell, per-sample cached lists: cluster.ids, celltype.ids,
#'   nearest.lm, fuzzy.graph.
#' @slot results list. All statistical outputs: lm, pb, marker, spec, nmf, pls,
#'   clustering, celltyping, features.
#'
#' @param x A TDRObj object.
#' @param name A character string naming the slot to access.
#' @param value The value to assign to the slot.
#' @param object A TDRObj object (used in show method).
#'
#' @name TDRObj-class
#' @rdname TDRObj-class
#' @exportClass TDRObj
#' @importFrom methods as setClass setMethod setReplaceMethod setAs new slot slot<- slotNames show validObject is
setClass(
  Class = "TDRObj",
  slots = c(
    cells          = "list",
    metadata       = "data.frame",
    config         = "list",
    integration    = "list",
    assay          = "list",
    landmark.embed = "list",
    landmark.annot = "list",
    graphs         = "list",
    density        = "list",
    sample.embed   = "list",
    cellmap        = "list",
    results        = "list"
  ),
  prototype = list(
    cells          = list(),
    metadata       = data.frame(),
    config         = list(),
    integration    = list(),
    assay          = list(),
    landmark.embed = list(),
    landmark.annot = list(),
    graphs         = list(),
    density        = list(),
    sample.embed   = list(),
    cellmap        = list(),
    results        = list()
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
#'   Legacy (17-slot) slot names are automatically mapped to new 13-slot paths.
#' @return A \code{TDRObj} object.
#' @export
TDRObj <- function(...) {
  args <- list(...)
  new_slots <- c("cells", "metadata", "config", "integration",
                 "assay", "landmark.embed", "landmark.annot",
                 "graphs", "density", "sample.embed", "cellmap",
                 "results")

  # ── Translate legacy slot names to new structure ──────────────────────────
  # raw.landmarks / landmarks / scaled.landmarks → assay$raw/expr/scaled
  for (nm in c("raw.landmarks", "landmarks", "scaled.landmarks")) {
    if (nm %in% names(args)) {
      if (is.null(args$assay)) args$assay <- list()
      sub <- switch(nm,
        raw.landmarks    = "raw",
        landmarks        = "expr",
        scaled.landmarks = "scaled")
      args$assay[[sub]] <- args[[nm]]
      args[[nm]] <- NULL
    }
  }

  # pca → landmark.embed$pca
  if ("pca" %in% names(args)) {
    if (is.null(args$landmark.embed)) args$landmark.embed <- list()
    pca_val <- args$pca
    # Handle old pca$embed → pca$coord rename
    if (!is.null(pca_val$embed) && is.null(pca_val$coord)) {
      pca_val$coord <- pca_val$embed
      pca_val$embed <- NULL
    }
    args$landmark.embed$pca <- pca_val
    args$pca <- NULL
  }

  # symphony.obj → integration$symphony.obj
  if ("symphony.obj" %in% names(args)) {
    if (is.null(args$integration)) args$integration <- list()
    args$integration$symphony.obj <- args$symphony.obj
    args$symphony.obj <- NULL
  }

  # graph → graphs (simple rename; or decompose if it has clustering/celltyping)
  if ("graph" %in% names(args)) {
    g <- args$graph
    args$graph <- NULL

    # Decompose @graph into new slots
    if (!is.null(g$uwot)) {
      if (is.null(args$integration)) args$integration <- list()
      args$integration$umap.model <- g$uwot
      if (!is.null(args$landmark.embed)) {
        args$landmark.embed$umap <- list(coord = g$uwot$embedding)
      } else {
        args$landmark.embed <- list(umap = list(coord = g$uwot$embedding))
      }
      if (is.null(args$graphs)) args$graphs <- list()
      if (!is.null(g$uwot$fgraph)) args$graphs$fgraph <- g$uwot$fgraph
    }
    if (!is.null(g$adj.matrix)) {
      if (is.null(args$graphs)) args$graphs <- list()
      args$graphs$adj.matrix <- g$adj.matrix
    }
    if (!is.null(g$snn)) {
      if (is.null(args$graphs)) args$graphs <- list()
      args$graphs$snn <- g$snn
    }
    if (!is.null(g$LE)) {
      if (is.null(args$landmark.embed)) args$landmark.embed <- list()
      le <- g$LE
      if (!is.null(le$embed) && is.null(le$coord)) {
        le$coord <- le$embed; le$embed <- NULL
      }
      args$landmark.embed$le <- le
    }
    if (!is.null(g$clustering)) {
      cl <- g$clustering
      if (!is.null(cl$ids)) {
        if (is.null(args$landmark.annot)) args$landmark.annot <- list()
        args$landmark.annot$clustering <- list(ids = cl$ids)
      }
      rest <- cl[setdiff(names(cl), "ids")]
      if (length(rest) > 0) {
        if (is.null(args$results)) args$results <- list()
        args$results$clustering <- rest
      }
    }
    if (!is.null(g$celltyping)) {
      ct <- g$celltyping
      if (!is.null(ct$ids)) {
        if (is.null(args$landmark.annot)) args$landmark.annot <- list()
        args$landmark.annot$celltyping <- list(ids = ct$ids)
      }
      rest <- ct[setdiff(names(ct), "ids")]
      if (length(rest) > 0) {
        if (is.null(args$results)) args$results <- list()
        args$results$celltyping <- rest
      }
    }
  }

  # map → density + cellmap (decompose each sub-element)
  if ("map" %in% names(args)) {
    m <- args$map
    args$map <- NULL
    if (is.null(args$density)) args$density <- list()
    if (is.null(args$cellmap)) args$cellmap <- list()

    if (!is.null(m$fdens))   args$density$fdens   <- m$fdens
    if (!is.null(m$Y))       args$density$Y       <- m$Y
    if (!is.null(m$.cache))  args$density$.cache  <- m$.cache

    # Only store .cache=NULL explicitly if it was passed
    if (is.element(".cache", names(m)) && is.null(m$.cache)) {
      args$density$.cache <- NULL
    }

    if (!is.null(m$cl.ct.to.ign)) args$density$ignored <- m$cl.ct.to.ign

    if (!is.null(m$clustering)) {
      cl <- m$clustering
      if (!is.null(cl$ids))        args$cellmap$cluster.ids <- cl$ids
      if (is.null(args$density$composition)) args$density$composition <- list()
      args$density$composition$clustering <- cl[setdiff(names(cl), "ids")]
    }
    if (!is.null(m$celltyping)) {
      ct <- m$celltyping
      if (!is.null(ct$ids))        args$cellmap$celltype.ids <- ct$ids
      if (is.null(args$density$composition)) args$density$composition <- list()
      args$density$composition$celltyping <- ct[setdiff(names(ct), "ids")]
    }
    if (!is.null(m$nearest.landmarks)) args$cellmap$nearest.lm  <- m$nearest.landmarks
    if (!is.null(m$fuzzy.graph))       args$cellmap$fuzzy.graph <- m$fuzzy.graph
    if (!is.null(m$lm)) {
      if (is.null(args$results)) args$results <- list()
      args$results$lm <- m$lm
    }
    if (!is.null(m$embedding)) {
      emb <- m$embedding
      if (is.null(args$sample.embed)) args$sample.embed <- list()
      if (!is.null(emb$pca))  args$sample.embed$pca  <- emb$pca
      if (!is.null(emb$traj)) args$sample.embed$traj <- emb$traj
      if (!is.null(emb$pePC)) args$sample.embed$pepc <- emb$pePC
    }
  }

  # specDE / pbDE / markerDE / nmfDE / plsDE / interact.plot → results$...
  for (pair in list(
    c("specDE",       "spec"),
    c("pbDE",         "pb"),
    c("markerDE",     "marker"),
    c("nmfDE",        "nmf"),
    c("plsDE",        "pls"),
    c("interact.plot","features")
  )) {
    old_nm <- pair[1]; new_nm <- pair[2]
    if (old_nm %in% names(args)) {
      if (is.null(args$results)) args$results <- list()
      args$results[[new_nm]] <- args[[old_nm]]
      args[[old_nm]] <- NULL
    }
  }

  # Drop any remaining unrecognised slot names
  args <- args[names(args) %in% new_slots]

  do.call(new, c("TDRObj", args))
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
# These shims map old 17-slot names to new 13-slot paths so legacy code
# continues to work during the migration period. After all phases are
# confirmed passing, the shim can be simplified to the 12 new slot names.

# New slot names (direct pass-through)
.tdrobj_new_slots <- c("cells", "metadata", "config", "integration",
                        "assay", "landmark.embed", "landmark.annot",
                        "graphs", "density", "sample.embed", "cellmap",
                        "results")

#' @rdname TDRObj-class
setMethod("$", "TDRObj", function(x, name) {
  if (name %in% .tdrobj_new_slots) {
    return(slot(x, name))
  }
  # Legacy name shims (old 17-slot names → new paths)
  switch(name,
    "landmarks"        = slot(x, "assay")$expr,
    "scaled.landmarks" = slot(x, "assay")$scaled,
    "raw.landmarks"    = slot(x, "assay")$raw,
    "pca"              = slot(x, "landmark.embed")$pca,
    "graph"            = slot(x, "graphs"),
    "map"              = {
      d <- slot(x, "density")
      comp <- d$composition
      d$composition <- NULL
      if (!is.null(comp)) c(d, comp) else d
    },
    "specDE"           = slot(x, "results")$spec,
    "pbDE"             = slot(x, "results")$pb,
    "markerDE"         = slot(x, "results")$marker,
    "nmfDE"            = slot(x, "results")$nmf,
    "plsDE"            = slot(x, "results")$pls,
    "interact.plot"    = slot(x, "results")$features,
    "symphony.obj"     = slot(x, "integration")$symphony.obj,
    stop("Unknown slot name: ", name)
  )
})

#' @rdname TDRObj-class
setReplaceMethod("$", "TDRObj", function(x, name, value) {
  if (name %in% .tdrobj_new_slots) {
    slot(x, name) <- value
    return(x)
  }
  # Legacy name shims for write
  switch(name,
    "landmarks"        = { s <- slot(x, "assay"); s$expr <- value; slot(x, "assay") <- s },
    "scaled.landmarks" = { s <- slot(x, "assay"); s$scaled <- value; slot(x, "assay") <- s },
    "raw.landmarks"    = { s <- slot(x, "assay"); s$raw <- value; slot(x, "assay") <- s },
    "pca"              = { s <- slot(x, "landmark.embed"); s$pca <- value; slot(x, "landmark.embed") <- s },
    "graph"            = { slot(x, "graphs") <- value },
    "map"              = { slot(x, "density") <- value },
    "specDE"           = { s <- slot(x, "results"); s$spec <- value; slot(x, "results") <- s },
    "pbDE"             = { s <- slot(x, "results"); s$pb <- value; slot(x, "results") <- s },
    "markerDE"         = { s <- slot(x, "results"); s$marker <- value; slot(x, "results") <- s },
    "nmfDE"            = { s <- slot(x, "results"); s$nmf <- value; slot(x, "results") <- s },
    "plsDE"            = { s <- slot(x, "results"); s$pls <- value; slot(x, "results") <- s },
    "interact.plot"    = { s <- slot(x, "results"); s$features <- value; slot(x, "results") <- s },
    "symphony.obj"     = { s <- slot(x, "integration"); s$symphony.obj <- value; slot(x, "integration") <- s },
    stop("Unknown slot name: ", name)
  )
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
  n_lm <- if (length(object@assay) > 0 && !is.null(object@assay$expr)) nrow(object@assay$expr) else 0L
  cat("  Landmarks:", n_lm, "\n")
  cat("  Graph computed:", length(object@graphs) > 0, "\n")
  cat("  Density/map computed:", !is.null(object@density$fdens), "\n")
})

# --- Coercion from list ---

setAs("list", "TDRObj", function(from) {
  # Map old 17-slot names to new 13-slot names where possible
  name_map <- c(
    landmarks        = "assay",
    scaled.landmarks = "assay",
    raw.landmarks    = "assay",
    pca              = "landmark.embed",
    graph            = "graphs",
    map              = "density",
    specDE           = "results",
    pbDE             = "results",
    markerDE         = "results",
    nmfDE            = "results",
    plsDE            = "results",
    interact.plot    = "results",
    symphony.obj     = "integration"
  )
  new_slots <- slotNames("TDRObj")
  # Only pass elements that are valid new slot names
  valid <- intersect(names(from), new_slots)
  do.call(new, c("TDRObj", from[valid]))
})

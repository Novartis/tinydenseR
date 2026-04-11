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

#' Manually assign cell type labels to clusters
#'
#' This function allows manual annotation of clusters by mapping one or more
#' clusters to biologically meaningful cell type labels. This is particularly
#' useful in cytometry experiments, where there typically is a larger number of
#' cells and lineage marker-based hierarchical analysis is common.
#'
#' \code{celltyping()} can be called at any point in the pipeline — before or
#' after \code{\link{get.map}}.
#' When called \strong{after} \code{get.map()}, it automatically refreshes all
#' celltyping-dependent downstream slots (cell-level IDs, composition matrices,
#' and the summary heatmap) without re-reading expression data or re-running
#' the UMAP transform.  Any existing \code{get.lm()} results that depend on
#' celltyping are invalidated with a warning so the user can re-run
#' \code{get.lm()}.
#'
#' @details
#' Manual cell type annotation is useful when:
#' \itemize{
#'   \item Domain expertise is needed for fine-grained annotations
#'   \item Custom groupings are required for specific analyses
#'   \item Working with cytometry data where marker combinations define cell types
#' }
#' 
#' Cell type labels assigned here will be used in downstream analyses:
#' \itemize{
#'   \item \code{\link{get.map}} propagates labels to all cells via nearest landmarks
#'   \item \code{\link{get.lm}} enables cell-type-level differential density testing in traditional analysis
#'   \item \code{\link{get.pbDE}} enables cell-type-specific pseudobulk differential expression analysis
#' }
#'
#' @param .celltyping.map Cell type assignments, supplied in one of two
#'   mutually exclusive formats:
#'   \describe{
#'     \item{Mode A — cluster map (named \code{list})}{
#'       List names are cell-type labels (e.g., \code{"CD4.T.cells"}), and
#'       each element is a character vector of cluster IDs
#'       (e.g., \code{c("cluster.01", "cluster.02")}) that belong to that
#'       cell type.  Every cluster in
#'       \code{.tdr.obj$landmark.annot$clustering$ids} must appear in exactly
#'       one element.
#'     }
#'     \item{Mode B — per-cell labels (named \code{character} vector)}{
#'       A named character vector where \code{names()} are the original cell
#'       IDs (before the \code{paste0(sample, "_", ...)} prefix added by
#'       tinydenseR) and values are cell-type labels.  The vector must cover
#'       **all** cells across all samples; landmark-only labels are extracted
#'       automatically.  Duplicate names are not allowed.
#'     }
#'   }
#' @param .verbose Logical; if \code{TRUE}, print progress messages when
#'   refreshing downstream slots (default \code{TRUE}).
#' @examples
#' \dontrun{
#' # After clustering with get.graph()
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#' 
#' # Map clusters to cell types (before get.map — traditional order)
#' celltype_map <- list(
#'   "CD4.T.cells" = c("cluster.01", "cluster.02"),
#'   "CD8.T.cells" = c("cluster.03"),
#'   "B.cells" = c("cluster.04", "cluster.05"),
#'   "unknown" = c("cluster.06")
#' )
#' lm.cells <- celltyping(lm.cells, celltype_map)
#'
#' # Late celltyping — after get.map (new supported workflow)
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.map() |>
#'   celltyping(celltype_map)
#'
#' # Mode B — per-cell labels (named character vector)
#' # cell_labels is a named character vector: names = original cell IDs,
#' # values = cell-type labels, covering ALL cells across all samples.
#' lm.cells <- celltyping(lm.cells, cell_labels)
#' }
#' @return The \code{.tdr.obj} with the following updated fields:
#'   \describe{
#'     \item{\code{$landmark.annot$celltyping$ids}}{Factor vector of cell type labels 
#'       for each landmark}
#'   }
#'   If \code{get.map()} has already been run, the following slots are also
#'   refreshed (no re-mapping required):
#'   \describe{
#'     \item{\code{$cellmap$celltyping$ids}}{Per-sample list of cell-level celltype
#'       assignments, re-derived from the existing fuzzy graph}
#'     \item{\code{$density$composition$celltyping$cell.count}}{Samples x cell
#'       types count matrix}
#'     \item{\code{$density$composition$celltyping$cell.perc}}{Samples x cell
#'       types percentage matrix}
#'   }
#'   If \code{get.lm()} has been run, any existing
#'   \code{$results$lm[[model]]$trad$celltyping} fits are invalidated (set to
#'   \code{NULL}) with a warning.
#' @seealso 
#'   \code{\link{get.map}} for automatic cell typing with symphony reference objects,
#'   \code{\link{lm.cluster}} for the clustering that produces cluster IDs used here,
#'   \code{\link{import_cell_annotations}} for automatic multi-column ingestion of
#'     cell-level annotations from metadata
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object.
#' @param ... Additional arguments passed to methods.
#' @export
celltyping <- function(x, ...) UseMethod("celltyping")

#' @rdname celltyping
#' @export
celltyping.TDRObj <-
  function(x,
           .celltyping.map,
           .verbose = TRUE,
           .name = NULL,
           ...){
    .tdr.obj <- x
    
    cell.pop <- NULL
    
    # Default name based on mode
    if (is.null(.name)) {
      .name <- if (is.list(.celltyping.map)) {
        paste0("manual.", format(Sys.time(), "%Y%m%d.%H%M%S"))
      } else {
        paste0("labels.", format(Sys.time(), "%Y%m%d.%H%M%S"))
      }
    }
    if (.name == "ids") {
      stop("'.name' cannot be 'ids' (reserved for the active solution).")
    }
    
    # ── Mode dispatch ────────────────────────────────────────────────
    if (is.list(.celltyping.map)) {
      # → Mode A: cluster-to-label mapping (existing behavior)
      .tdr.obj <- .celltyping_mode_a(.tdr.obj, .celltyping.map)
    } else if (is.character(.celltyping.map) && !is.null(names(.celltyping.map))) {
      # → Mode B: per-cell labels (new path)
      .tdr.obj <- .celltyping_mode_b(.tdr.obj, .celltyping.map)
    } else {
      stop("`.celltyping.map` must be either:\n",
           "  - A named list (cluster->label mapping), or\n",
           "  - A named character vector (cell->label mapping).")
    }
    
    # Store as named solution
    .tdr.obj@landmark.annot$celltyping[[.name]] <-
      .tdr.obj@landmark.annot$celltyping$ids
    
    # Refresh all celltyping-dependent slots (shared by both modes)
    .tdr.obj <- .refresh_celltyping(.tdr.obj, .verbose = .verbose)
    
    return(.tdr.obj)
    
  }

# ──────────────────────────────────────────────────────────────────────
# Mode A: cluster-to-label mapping (original behavior)
# ──────────────────────────────────────────────────────────────────────
.celltyping_mode_a <- function(.tdr.obj, .celltyping.map) {
  
  if(names(x = .celltyping.map) |>
     is.null()){
    stop("Cell type names are missing. Each element in .celltyping.map must have a name.\n",
         "Example: list(\"CD4.T.cells\" = c(\"cluster.01\"), \"CD8.T.cells\" = c(\"cluster.02\"))")
  }
  
  if(names(x = .celltyping.map) |>
     duplicated() |>
     any()){
    stop(paste0("Duplicate cell type names detected: ",
                paste(names(.celltyping.map)[duplicated(names(.celltyping.map))],
                      collapse = ", "),
                "\nEach cell type must have a unique name."))
  }
  
  cls.in.map <-
    unlist(x = .celltyping.map,
           use.names = FALSE)
  
  if(anyDuplicated(x = cls.in.map)){
    stop(paste0("Cluster(s) mapped to multiple cell types: ",
                paste(cls.in.map[duplicated(x = cls.in.map)],
                      collapse = ", "),
                "\nEach cluster can only belong to one cell type."))
  }
  
  if(any(!(unique(x = .tdr.obj@landmark.annot$clustering$ids) %in%
           cls.in.map))){
    stop(paste0("Every cluster must be mapped to a cell type. Unmapped clusters: ",
                paste(unique(x = .tdr.obj@landmark.annot$clustering$ids)[
                  !(unique(x = .tdr.obj@landmark.annot$clustering$ids) %in%
                      cls.in.map)],
                  collapse = ", "),
                "\nAdd these to .celltyping.map or merge them with existing clusters."))
    
  }
  
  if(!all(cls.in.map %in%
          unique(x = .tdr.obj@landmark.annot$clustering$ids))){
    stop(paste0("Invalid cluster ID(s) in .celltyping.map: ",
                paste(cls.in.map[!(cls.in.map %in%
                                     unique(x = .tdr.obj@landmark.annot$clustering$ids))],
                      collapse = ", "),
                "\nThese clusters do not exist in .tdr.obj$landmark.annot$clustering$ids.",
                "\nRun lm.cluster() to see available cluster IDs."))
  }
  
  # Create cell type labels by:
  # 1. Inverting the map: cluster IDs -> cell type names
  # 2. Indexing by cluster IDs to get corresponding cell type for each landmark
  # 3. Converting to factor to maintain consistency with clustering$ids
  .tdr.obj@landmark.annot$celltyping$ids <-
    .celltyping.map |>
    (\(x)
     stats::setNames(
       object = names(x = x) |>
         rep(times = lapply(X = x,
                            length) |>
               unlist(use.names = FALSE)),
       nm = unlist(x = x,
                   use.names = FALSE))[
                     as.character(x = .tdr.obj@landmark.annot$clustering$ids)]
    )() |>
    as.factor()
  
  # Store provenance
  
  .tdr.obj
}

# ──────────────────────────────────────────────────────────────────────
# Mode B: per-cell labels (new)
# ──────────────────────────────────────────────────────────────────────
.celltyping_mode_b <- function(.tdr.obj, .celltyping.map) {

  # --- Validate: length matches total cell count ---
  expected_n <- sum(.tdr.obj@config$sampling$n.cells)
  if (length(.celltyping.map) != expected_n) {
    stop("Length of .celltyping.map (", length(.celltyping.map),
         ") does not match total cell count (", expected_n, ").")
  }

  # --- Get landmark row names from the assay matrix ---
  lm_names <- rownames(.tdr.obj@assay$expr)

  # config$key names tell us which sample each landmark belongs to
  sample_names <- names(.tdr.obj@cells)
  n_cells      <- .tdr.obj@config$sampling$n.cells
  lm_sample    <- names(.tdr.obj@config$key)

  labels <- character(length(lm_names))

  if (anyDuplicated(names(.celltyping.map))) {
    # ── Non-unique cell IDs (e.g. cytoset event_N) ──
    # Positional split: the vec MUST be ordered by sample.
    ct_splits <- split(
      .celltyping.map,
      rep(sample_names, times = n_cells)
    )[sample_names]

    # Validate: no duplicate cell IDs within any sample
    for (sn in sample_names) {
      chunk_names <- names(ct_splits[[sn]])
      if (anyDuplicated(chunk_names)) {
        dupes <- unique(chunk_names[duplicated(chunk_names)])
        stop("Duplicate cell IDs in sample '", sn, "': ",
             paste(head(dupes, 5), collapse = ", "),
             if (length(dupes) > 5)
               paste0(" ... (", length(dupes), " total)"),
             "\nEach cell must appear exactly once within a sample.")
      }
    }

    for (sn in sample_names) {
      idx <- which(lm_sample == sn)
      if (length(idx) == 0L) next
      stripped <- substring(lm_names[idx], nchar(sn) + 2L)
      m <- match(stripped, names(ct_splits[[sn]]))
      if (anyNA(m)) {
        missing <- lm_names[idx][is.na(m)]
        stop("Could not find labels for ", length(missing),
             " landmark(s) in sample '", sn, "'.\n",
             "First few: ", paste(head(missing, 5), collapse = ", "), "\n",
             "Ensure names(.celltyping.map) contain the original cell IDs.")
      }
      labels[idx] <- ct_splits[[sn]][m]
    }

  } else {
    # ── Globally unique cell IDs (Seurat barcodes, matrix colnames, ──
    # ── files-backend rownames): direct match against full map.     ──
    for (sn in sample_names) {
      idx <- which(lm_sample == sn)
      if (length(idx) == 0L) next
      stripped <- substring(lm_names[idx], nchar(sn) + 2L)
      m <- match(stripped, names(.celltyping.map))
      if (anyNA(m)) {
        missing <- lm_names[idx][is.na(m)]
        stop("Could not find labels for ", length(missing),
             " landmark(s) in sample '", sn, "'.\n",
             "First few: ", paste(head(missing, 5), collapse = ", "), "\n",
             "Ensure names(.celltyping.map) contain the original cell IDs.")
      }
      labels[idx] <- .celltyping.map[m]
    }
  }

  # --- Warn if all landmarks get the same label ---
  if (length(unique(labels)) == 1L) {
    warning("All landmarks received the same cell-type label ('",
            labels[1],
            "'). This is likely a user error.", call. = FALSE)
  }

  # --- Assign ---
  .tdr.obj@landmark.annot$celltyping$ids <- factor(labels)
  names(.tdr.obj@landmark.annot$celltyping$ids) <- lm_names

  .tdr.obj
}

# ──────────────────────────────────────────────────────────────────────
# Private helpers for annotation refresh (parameterised by .annot.type)
# ──────────────────────────────────────────────────────────────────────

# Mapping from annotation type to internal slot names
.annot_slot_map <- function(.annot.type) {
  switch(.annot.type,
    clustering = list(
      slot      = "clustering",
      comp_slot = "clustering"
    ),
    celltyping = list(
      slot      = "celltyping",
      comp_slot = "celltyping"
    ),
    stop("Unknown .annot.type: ", .annot.type)
  )
}

#' Re-derive cell-level annotation IDs from fuzzy graph + landmark labels
#'
#' Reads the cached (on-disk or in-memory) fuzzy graph for each sample and
#' applies the current landmark-level annotation via the same
#' weighted-voting logic used in \code{get.map()}.  No expression data is
#' re-read and no UMAP transform is repeated.
#'
#' @param .tdr.obj A TDRObj where \code{get.map()} has already been run.
#' @param .annot.type Character: \code{"celltyping"} (default) or
#'   \code{"clustering"}.
#' @param .verbose Logical; print progress messages.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.relabel_cellmap <- function(.tdr.obj, .annot.type = "celltyping", .verbose = FALSE) {
  
  # R CMD check appeasement
  cell <- landmark <- label <- x <- confidence <- i <- j <- NULL
  
  smap <- .annot_slot_map(.annot.type)
  
  annot_ids <- .tdr.obj@landmark.annot[[.annot.type]]$ids
  .label.confidence <- if (!is.null(.tdr.obj@config$label.confidence)) {
    .tdr.obj@config$label.confidence
  } else {
    0.5
  }

  .tdr_validate_label_confidence(.label.confidence)
  
  # Detect on-disk mode: check if any entry in the slot is a path string
  .is_path_string <- function(val) {
    is.character(val) && length(val) == 1L && !is.null(attr(val, "schema_v"))
  }
  existing_entries <- .tdr.obj@cellmap[[smap$slot]]$ids
  if (is.null(existing_entries)) {
    # Check fuzzy.graphs to determine on-disk status
    existing_entries <- .tdr.obj@cellmap$fuzzy.graphs
  }
  is_cached <- !is.null(existing_entries) &&
    length(existing_entries) > 0L &&
    .is_path_string(existing_entries[[1]])
  cache_root <- .tdr.obj@config$.cache.root
  
  sample_names <- names(.tdr.obj@cells)
  
  new_ids <- stats::setNames(
    vector("list", length(sample_names)),
    sample_names
  )
  
  for (sn in sample_names) {
    
    if (isTRUE(.verbose)) {
      message("-> Relabeling ", .annot.type, " for sample: ", sn)
    }
    
    # Read fuzzy graph (from disk cache or in-memory)
    fgraph <- .tdr_get_map_slot(.tdr.obj, "fuzzy.graphs", sn)
    
    n_cells <- nrow(fgraph)
    
    new_ids[[sn]] <-
      .tdr_transfer_labels(
        .method = "fuzzy",
        .label.confidence = .label.confidence,
        .n.cells = n_cells,
        .cell.names = names(.tdr_get_map_slot(.tdr.obj, "clustering", sn)),
        .fgraph = fgraph,
        .landmark.labels = annot_ids
      )
    
    # Write to on-disk cache if active
    if (is_cached && !is.null(cache_root)) {
      path_str <- .tdr_cache_write(
        object = new_ids[[sn]],
        cache_dir = cache_root,
        slot_name = smap$slot,
        sample_name = sn
      )
      .tdr.obj@cellmap[[smap$slot]]$ids[[sn]] <- path_str
    }
  }
  
  # Store in-memory (data directly) if not on-disk
  if (!is_cached) {
    .tdr.obj@cellmap[[smap$slot]]$ids <- new_ids
  }
  
  .tdr.obj
}

#' Recompute composition matrices from cell-level annotation IDs
#'
#' Aggregates cell-level assignments into samples x populations
#' count and percentage matrices.
#'
#' @param .tdr.obj A TDRObj with freshly relabeled annotation IDs.
#' @param .annot.type Character: \code{"celltyping"} (default) or
#'   \code{"clustering"}.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.recompute_composition <- function(.tdr.obj, .annot.type = "celltyping") {
  
  smap <- .annot_slot_map(.annot.type)
  sample_names <- names(.tdr.obj@cells)
  
  counts <- lapply(
    X = stats::setNames(sample_names, sample_names),
    FUN = function(sn) {
      ids <- .tdr_get_map_slot(.tdr.obj, smap$slot, sn)
      tbl <- table(ids)
      stats::setNames(as.vector(tbl), names(tbl))
    }
  )
  
  .tdr.obj@density$composition[[smap$comp_slot]]$cell.count <-
    counts |>
    dplyr::bind_rows(.id = "sample") |>
    as.data.frame() |>
    (\(x)
     `rownames<-`(x = as.matrix(x = x[, colnames(x = x) != "sample"]),
                  value = x$sample)
    )() |>
    (\(x)
     x[, colnames(x = x) |>
         sort()]
    )()
  
  .tdr.obj@density$composition[[smap$comp_slot]]$cell.count[
    is.na(x = .tdr.obj@density$composition[[smap$comp_slot]]$cell.count)
  ] <- 0
  
  .tdr.obj@density$composition[[smap$comp_slot]]$cell.perc <-
    (.tdr.obj@density$composition[[smap$comp_slot]]$cell.count * 100) /
    Matrix::rowSums(x = .tdr.obj@density$composition[[smap$comp_slot]]$cell.count)
  
  .tdr.obj
}

#' Invalidate stale traditional analysis fits
#'
#' Sets \code{@results$lm[[model]]$trad[[.annot.type]]} to \code{NULL}
#' for every stored model and emits a warning.
#'
#' @param .tdr.obj A TDRObj.
#' @param .annot.type Character: \code{"celltyping"} (default) or
#'   \code{"clustering"}.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.invalidate_trad <- function(.tdr.obj, .annot.type = "celltyping") {
  
  model_names <- names(.tdr.obj@results$lm)
  
  for (mn in model_names) {
    if (!is.null(.tdr.obj@results$lm[[mn]]$trad[[.annot.type]])) {
      .tdr.obj@results$lm[[mn]]$trad[[.annot.type]] <- NULL
      warning("Invalidated trad$", .annot.type, " fit for model '", mn,
              "'. Re-run get.lm() to recompute.",
              call. = FALSE)
    }
  }
  
  .tdr.obj
}

#' Warn about potentially stale pseudobulk / marker DE results
#'
#' Checks whether any \code{get.pbDE()} results (design or marker mode)
#' were computed with \code{.id.from} matching \code{.annot.type} and
#' emits a warning if so.
#'
#' @param .tdr.obj A TDRObj.
#' @param .annot.type Character: \code{"celltyping"} (default) or
#'   \code{"clustering"}.
#' @return The (unmodified) \code{.tdr.obj}, invisibly.
#' @keywords internal
.warn_stale_de_results <- function(.tdr.obj, .annot.type = "celltyping") {
  
  # Check pseudobulk DE results
  if (!is.null(.tdr.obj@results$pb)) {
    for (mn in names(.tdr.obj@results$pb)) {
      for (pn in names(.tdr.obj@results$pb[[mn]])) {
        res <- .tdr.obj@results$pb[[mn]][[pn]]
        if (!is.null(res$.id.from) && res$.id.from == .annot.type) {
          warning("Pseudobulk DE results for model '", mn,
                  "', population '", pn,
                  "' were computed with old ", .annot.type, " labels and may be stale.",
                  call. = FALSE)
        }
      }
    }
  }
  
  # Check marker DE results
  if (!is.null(.tdr.obj@results$marker)) {
    for (mn in names(.tdr.obj@results$marker)) {
      for (cn in names(.tdr.obj@results$marker[[mn]])) {
        res <- .tdr.obj@results$marker[[mn]][[cn]]
        if (!is.null(res$.id.from) && res$.id.from == .annot.type) {
          warning("Marker DE results for model '", mn,
                  "', comparison '", cn,
                  "' were computed with old ", .annot.type, " labels and may be stale.",
                  call. = FALSE)
        }
      }
    }
  }
  
  invisible(.tdr.obj)
}

#' Refresh all celltyping-dependent slots
#'
#' Inspects the TDRObj and updates every populated celltyping-dependent slot.
#' This is called at the end of \code{celltyping()} so that late-bound
#' celltyping (after \code{get.map()}) produces a fully consistent object.
#'
#' Slots that have not yet been computed (e.g. \code{get.map()} has not run)
#' are simply skipped — the object remains valid at whatever pipeline stage
#' it is in.
#'
#' @param .tdr.obj A TDRObj with \code{@landmark.annot$celltyping$ids}
#'   already set to the new values.
#' @param .verbose Logical; print progress messages.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.refresh_celltyping <- function(.tdr.obj, .verbose = FALSE) {
  
  # --- Cell-level IDs & composition (only if get.map() has run) ---
  if (!is.null(.tdr.obj@density$fdens)) {
    
    if (isTRUE(.verbose)) {
      message("-> get.map() results detected; refreshing cell-level celltype IDs and composition...")
    }
    
    .tdr.obj <- .relabel_cellmap(.tdr.obj, .annot.type = "celltyping", .verbose = .verbose)
    .tdr.obj <- .recompute_composition(.tdr.obj, .annot.type = "celltyping")
    
    if (isTRUE(.verbose)) {
      message("-> Cell-level celltype IDs and composition matrices refreshed.")
    }
  }
  
  # --- trad celltyping fit (only if get.lm() was run) ---
  if (!is.null(.tdr.obj@results$lm)) {
    .tdr.obj <- .invalidate_trad(.tdr.obj, .annot.type = "celltyping")
  }
  
  # --- Warn about potentially stale pbDE / markerDE results ---
  .warn_stale_de_results(.tdr.obj, .annot.type = "celltyping")
  
  .tdr.obj
}

#' Set active celltyping solution
#'
#' Switches the active celltyping to a previously stored solution and
#' refreshes all celltyping-dependent downstream slots.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, SingleCellExperiment, or HDF5AnnData
#'   (anndataR) object.
#' @param .column.name Character: name of the stored celltyping solution to activate.
#' @param .verbose Logical (default TRUE).
#' @param ... Additional arguments passed to methods.
#' @return Updated object with the selected celltyping as active \code{$ids}.
#' @export
set_active_celltyping <- function(x, ...) UseMethod("set_active_celltyping")

#' @rdname set_active_celltyping
#' @export
set_active_celltyping.TDRObj <-
  function(x,
           .column.name,
           .verbose = TRUE,
           ...) {
    .tdr.obj <- x
    
    if (.column.name == "ids") {
      stop("'.column.name' cannot be 'ids'. Specify a named solution.")
    }
    
    stored <- .tdr.obj@landmark.annot$celltyping[[.column.name]]
    if (is.null(stored)) {
      avail <- setdiff(names(.tdr.obj@landmark.annot$celltyping), "ids")
      stop("Celltyping solution '", .column.name, "' not found.\n",
           "Available solutions: ", paste(avail, collapse = ", "))
    }
    
    .tdr.obj@landmark.annot$celltyping$ids <- stored
    .tdr.obj <- .refresh_celltyping(.tdr.obj, .verbose = .verbose)
    
    return(.tdr.obj)
  }


# ──────────────────────────────────────────────────────────────────────
# import_cell_annotations: multi-column cell-level annotation ingestion
# ──────────────────────────────────────────────────────────────────────

#' Import multiple cell-level annotation columns as landmark-level celltyping solutions
#'
#' Scans a cell-level metadata data.frame for categorical columns and
#' imports each as a named landmark-level celltyping solution via
#' \code{\link{celltyping}} Mode B.  After importing, users can switch
#' between solutions with \code{\link{set_active_celltyping}}.
#'
#' This function is designed to be called after \code{\link{get.graph}}
#' and \strong{before} \code{\link{get.map}} in the pipeline.  Because
#' \code{get.map()} has not yet run, \code{.refresh_celltyping()} inside
#' each \code{celltyping()} call is a no-op, making the loop efficient.
#'
#' @section Column selection criteria:
#' A column in \code{.cell.meta} is imported if it is:
#' \itemize{
#'   \item \code{character} or \code{factor} type (not numeric, logical,
#'     integer, etc.)
#'   \item Not the \code{.sample.var} column
#'   \item Has more than 1 unique non-NA value
#'   \item Has fewer unique non-NA values than total cells (excludes
#'     ID-like columns)
#'   \item Column name is not \code{"ids"} (reserved by the multi-solution
#'     system)
#' }
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, or SingleCellExperiment
#'   object.
#' @param .cell.meta A \code{data.frame} with rownames set to cell IDs
#'   (original IDs before the tinydenseR sample prefix).  Must contain
#'   at least the \code{.sample.var} column.
#' @param .sample.var Character(1).  Column name in \code{.cell.meta}
#'   identifying sample membership.
#' @param .celltype.vec Character(1) or \code{NULL}.  If non-NULL, the
#'   name of a column in \code{.cell.meta} that should be set as the
#'   active \code{$ids} after import.
#' @param .verbose Logical (default \code{TRUE}).  If \code{TRUE}, prints
#'   a message listing the imported columns.
#' @param ... Additional arguments passed to methods.
#'
#' @return The updated object with imported celltyping solutions stored in
#'   \code{@landmark.annot$celltyping$<column_name>}.
#'
#' @seealso \code{\link{celltyping}}, \code{\link{set_active_celltyping}},
#'   \code{\link{list_celltyping_solutions}}
#'
#' @examples
#' \dontrun{
#' tdr.obj <- setup.tdr.obj(.cells, .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#'
#' # Import all categorical cell-level columns
#' tdr.obj <- import_cell_annotations(tdr.obj,
#'   .cell.meta = cell_metadata,
#'   .sample.var = "sample_id",
#'   .celltype.vec = "cell_type_l1"
#' )
#'
#' # Check what was imported
#' list_celltyping_solutions(tdr.obj)
#'
#' # Switch to a different annotation
#' tdr.obj <- set_active_celltyping(tdr.obj, "cell_type_l2")
#' }
#'
#' @export
import_cell_annotations <- function(x, ...) UseMethod("import_cell_annotations")

#' @rdname import_cell_annotations
#' @export
import_cell_annotations.TDRObj <-
  function(x,
           .cell.meta,
           .sample.var,
           .celltype.vec = NULL,
           .verbose = TRUE,
           ...) {
    .tdr.obj <- x

    # --- Validate inputs ---
    if (!inherits(.cell.meta, "data.frame")) {
      stop(".cell.meta must be a data.frame.")
    }
    rn <- rownames(.cell.meta)
    if (is.null(rn) ||
        identical(rn, as.character(seq_len(nrow(.cell.meta))))) {
      stop(".cell.meta must have meaningful rownames (cell IDs).")
    }
    if (!is.character(.sample.var) || length(.sample.var) != 1L) {
      stop(".sample.var must be a single character string.")
    }
    if (!(.sample.var %in% colnames(.cell.meta))) {
      stop(".sample.var '", .sample.var,
           "' not found in .cell.meta.")
    }

    # --- Identify cell-level categorical columns ---
    # Cell-level = NOT constant within every sample (inverse of get.meta logic)
    sample_ids <- .cell.meta[[.sample.var]]
    n_cells <- nrow(.cell.meta)
    all_cols <- colnames(.cell.meta)

    # Exclude .sample.var and reserved "ids"
    candidate_cols <- setdiff(all_cols, c(.sample.var, "ids"))

    # Filter to character / factor columns
    candidate_cols <- candidate_cols[
      vapply(candidate_cols, function(col) {
        is.character(.cell.meta[[col]]) || is.factor(.cell.meta[[col]])
      }, logical(1))
    ]

    if (length(candidate_cols) == 0L) {
      if (isTRUE(.verbose)) {
        message("import_cell_annotations: no categorical columns found; ",
                "returning object unchanged.")
      }
      return(.tdr.obj)
    }

    # Filter: >1 unique non-NA value AND <n_cells unique values (not ID-like)
    # AND varies within at least one sample (cell-level, not sample-level)
    qualifying <- vapply(candidate_cols, function(col) {
      vals <- .cell.meta[[col]]
      uvals <- unique(vals[!is.na(vals)])
      n_unique <- length(uvals)
      if (n_unique <= 1L || n_unique >= n_cells) return(FALSE)
      # Check if cell-level: varies within at least one sample
      any(vapply(
        split(vals, sample_ids),
        function(chunk) length(unique(chunk[!is.na(chunk)])) > 1L,
        logical(1)
      ))
    }, logical(1))

    qualifying_cols <- candidate_cols[qualifying]

    if (length(qualifying_cols) == 0L) {
      if (isTRUE(.verbose)) {
        message("import_cell_annotations: no qualifying cell-level ",
                "categorical columns found; returning object unchanged.")
      }
      return(.tdr.obj)
    }

    if (isTRUE(.verbose)) {
      message("import_cell_annotations: importing ",
              length(qualifying_cols), " column(s): ",
              paste(qualifying_cols, collapse = ", "))
    }

    # --- Save current $ids state to restore if .celltype.vec is NULL ---
    ids_before <- .tdr.obj@landmark.annot$celltyping$ids

    # --- Import each column via celltyping Mode B ---
    for (col in qualifying_cols) {
      vec <- stats::setNames(
        as.character(.cell.meta[[col]]),
        rownames(.cell.meta)
      )
      .tdr.obj <- celltyping.TDRObj(
        .tdr.obj,
        .celltyping.map = vec,
        .name = col,
        .verbose = FALSE
      )
    }

    # --- Activate requested column or restore prior $ids ---
    if (!is.null(.celltype.vec) && .celltype.vec %in% qualifying_cols) {
      .tdr.obj <- set_active_celltyping.TDRObj(
        .tdr.obj,
        .column.name = .celltype.vec,
        .verbose = .verbose
      )
    } else {
      # Restore prior $ids state (even if NULL)
      .tdr.obj@landmark.annot$celltyping$ids <- ids_before
      if (!is.null(.celltype.vec) && isTRUE(.verbose)) {
        message("import_cell_annotations: .celltype.vec '",
                .celltype.vec,
                "' was not among the imported columns; ",
                "active $ids unchanged.")
      }
    }

    return(.tdr.obj)
  }


# ──────────────────────────────────────────────────────────────────────
# list_celltyping_solutions: utility to list stored solutions
# ──────────────────────────────────────────────────────────────────────

#' List stored celltyping solutions
#'
#' Returns the names of all stored celltyping solutions (excluding the
#' active \code{$ids} slot).
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, or SingleCellExperiment
#'   object.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of solution names.
#'
#' @seealso \code{\link{set_active_celltyping}}, \code{\link{import_cell_annotations}}
#'
#' @examples
#' \dontrun{
#' list_celltyping_solutions(tdr.obj)
#' # [1] "cell_type_l1" "cell_type_l2" "azimuth_pred"
#' }
#'
#' @export
list_celltyping_solutions <- function(x, ...) UseMethod("list_celltyping_solutions")

#' @rdname list_celltyping_solutions
#' @export
list_celltyping_solutions.TDRObj <- function(x, ...) {
  nms <- names(x@landmark.annot$celltyping)
  if (is.null(nms)) return(character(0))
  setdiff(nms, "ids")
}

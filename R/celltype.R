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
#'     \item{\code{$landmark.annot$celltyping$map}}{Mode A only: the
#'       \code{.celltyping.map} list (stored for provenance).
#'       \code{NULL} in Mode B.}
#'     \item{\code{$landmark.annot$celltyping$mode}}{\code{"cluster_map"}
#'       (Mode A) or \code{"cell_labels"} (Mode B)}
#'     \item{\code{$results$celltyping$median.exprs}}{Matrix of median expression 
#'       values per cell type (rows) by features (columns). For RNA data, 
#'       displays z-scored expression of top PC-loading genes. For cytometry 
#'       data, shows z-scored marker intensities.}
#'     \item{\code{$results$celltyping$pheatmap}}{Heatmap object (from 
#'       \code{pheatmap} package) visualizing median expression patterns across 
#'       cell types}
#'   }
#'   If \code{get.map()} has already been run, the following slots are also
#'   refreshed (no re-mapping required):
#'   \describe{
#'     \item{\code{$cellmap$celltype.ids}}{Per-sample list of cell-level celltype
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
#'   \code{\link{lm.cluster}} for the clustering that produces cluster IDs used here
#' @param x Object to operate on (TDRObj, Seurat, or SingleCellExperiment).
#' @param ... Additional arguments passed to methods.
#' @export
celltyping <- function(x, ...) UseMethod("celltyping")

#' @rdname celltyping
#' @export
celltyping.TDRObj <-
  function(x,
           .celltyping.map,
           .verbose = TRUE,
           ...){
    .tdr.obj <- x
    
    cell.pop <- NULL
    
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
  
  .tdr.obj@landmark.annot$celltyping <-
    vector(mode = "list")
  
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
  .tdr.obj@landmark.annot$celltyping$map  <- .celltyping.map
  .tdr.obj@landmark.annot$celltyping$mode <- "cluster_map"
  
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

  # --- Split .celltyping.map by sample using n.cells counts ---
  sample_names <- names(.tdr.obj@cells)
  n_cells      <- .tdr.obj@config$sampling$n.cells
  ct_splits    <- split(
    .celltyping.map,
    rep(sample_names, times = n_cells)
  )[sample_names]              # preserve sample order

  # --- Validate: no duplicate cell IDs within any sample ---
  for (sn in sample_names) {
    chunk_names <- names(ct_splits[[sn]])
    if (anyDuplicated(chunk_names)) {
      dupes <- unique(chunk_names[duplicated(chunk_names)])
      stop("Duplicate cell IDs in sample '", sn, "': ",
           paste(head(dupes, 5), collapse = ", "),
           if (length(dupes) > 5) paste0(" ... (", length(dupes), " total)"),
           "\nEach cell must appear exactly once within a sample.")
    }
  }

  # --- Match sample-by-sample (safe for non-unique cell IDs) ---
  labels <- character(length(lm_names))

  # config$key names tell us which sample each landmark belongs to
  lm_sample <- names(.tdr.obj@config$key)

  for (sn in sample_names) {
    idx    <- which(lm_sample == sn)
    if (length(idx) == 0L) next
    prefix <- paste0("^", sn, "_")
    stripped <- sub(prefix, "", lm_names[idx])
    ct_sample <- ct_splits[[sn]]
    m <- match(stripped, names(ct_sample))
    if (anyNA(m)) {
      missing <- lm_names[idx][is.na(m)]
      stop("Could not find labels for ", length(missing),
           " landmark(s) in sample '", sn, "'.\n",
           "First few: ", paste(head(missing, 5), collapse = ", "), "\n",
           "Ensure names(.celltyping.map) contain the original cell IDs.")
    }
    labels[idx] <- ct_sample[m]
  }

  # --- Warn if all landmarks get the same label ---
  if (length(unique(labels)) == 1L) {
    warning("All landmarks received the same cell-type label ('",
            labels[1],
            "'). This is likely a user error.", call. = FALSE)
  }

  # --- Assign ---
  .tdr.obj@landmark.annot$celltyping <- list()
  .tdr.obj@landmark.annot$celltyping$ids <- factor(labels)
  names(.tdr.obj@landmark.annot$celltyping$ids) <- lm_names

  # Store provenance (no cluster map in Mode B)
  .tdr.obj@landmark.annot$celltyping$map  <- NULL
  .tdr.obj@landmark.annot$celltyping$mode <- "cell_labels"

  .tdr.obj
}

# ──────────────────────────────────────────────────────────────────────
# Private helpers for celltyping refresh
# ──────────────────────────────────────────────────────────────────────

#' Recompute landmark-level celltyping summaries (median.exprs + pheatmap)
#'
#' Derives the \code{@results$celltyping$median.exprs} matrix and
#' \code{@results$celltyping$pheatmap} heatmap object from the current
#' \code{@landmark.annot$celltyping$ids}.  This is a pure function of
#' landmark expression data and celltype labels — no cell-level data is
#' needed.
#'
#' @param .tdr.obj A TDRObj with \code{@landmark.annot$celltyping$ids}
#'   already set.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.recompute_celltyping_summaries <- function(.tdr.obj) {
  
  cell.pop <- NULL
  
  # For RNA: select top PC-loading genes for heatmap visualization
  # Takes top 3 positive and top 3 negative loadings per PC to capture
  # genes that drive the major axes of variation
  if(.tdr.obj@config$assay.type == "RNA") {
    
    top <-
      apply(X = .tdr.obj@landmark.embed$pca$rotation,
            MARGIN = 2,
            FUN = function(PC.rot){
              
              order(PC.rot,
                    decreasing = TRUE) |>
                (\(x)
                 c(utils::head(x = x, 
                               n = 3),
                   utils::tail(x = x,
                               n = 3))
                )()
              
            }) |>
      as.vector() |> 
      (\(x)
       rownames(x = .tdr.obj@landmark.embed$pca$rotation)[x]
      )() |>
      unique()
    
  }
  
  .tdr.obj@results$celltyping$median.exprs <-
    (if(.tdr.obj@config$assay.type == "RNA") .tdr.obj@assay$scaled[,top] else .tdr.obj@assay$expr) |>
    dplyr::as_tibble() |>
    cbind(cell.pop = as.character(x = .tdr.obj@landmark.annot$celltyping$ids)) |>
    dplyr::group_by(cell.pop) |>
    dplyr::summarize_all(.funs = stats::median) |>
    as.data.frame() |>
    (\(x)
     `rownames<-`(x = x[,colnames(x = x) != "cell.pop"],
                  value = x$cell.pop)
    )() |>
    as.matrix()
  
  .tdr.obj@results$celltyping$pheatmap <-
    pheatmap::pheatmap(mat = .tdr.obj@results$celltyping$median.exprs,
                       color = grDevices::colorRampPalette(
                         unname(obj =
                                  Color.Palette[1,c(1,6,2)]))(100),
                       kmeans_k = NA,
                       breaks = NA,
                       border_color = NA,
                       scale = "none",
                       angle_col = 90,
                       cluster_rows = nrow(.tdr.obj@results$celltyping$median.exprs) > 1,
                       cluster_cols = ncol(.tdr.obj@results$celltyping$median.exprs) > 1,
                       cellwidth = 20,
                       cellheight = 20,
                       treeheight_row = 20,
                       treeheight_col = 20,
                       silent = TRUE)
  
  .tdr.obj
}

#' Re-derive cell-level celltype IDs from fuzzy graph + landmark labels
#'
#' Reads the cached (on-disk or in-memory) fuzzy graph for each sample and
#' applies the current \code{@landmark.annot$celltyping$ids} via the same
#' weighted-voting logic used in \code{get.map()}.  No expression data is
#' re-read and no UMAP transform is repeated.
#'
#' @param .tdr.obj A TDRObj where \code{get.map()} has already been run.
#' @param .verbose Logical; print progress messages.
#' @return The modified \code{.tdr.obj} with updated
#'   \code{@cellmap$celltype.ids} (in-memory) or
#'   \code{@density$.cache$manifests$celltyping.ids} (on-disk).
#' @keywords internal
.relabel_cellmap <- function(.tdr.obj, .verbose = FALSE) {
  
  # R CMD check appeasement
  cell <- landmark <- label <- x <- confidence <- i <- j <- NULL
  
  ct_ids <- .tdr.obj@landmark.annot$celltyping$ids
  .label.confidence <- if (!is.null(.tdr.obj@config$label.confidence)) {
    .tdr.obj@config$label.confidence
  } else {
    0.8
  }
  
  cache <- .tdr.obj@density$.cache
  is_cached <- !is.null(cache) && isTRUE(cache$on.disk)
  
  sample_names <- names(.tdr.obj@cells)
  
  new_celltype_ids <- stats::setNames(
    vector("list", length(sample_names)),
    sample_names
  )
  
  for (sn in sample_names) {
    
    if (isTRUE(.verbose)) {
      message("-> Relabeling celltypes for sample: ", sn)
    }
    
    # Read fuzzy graph (from disk cache or in-memory)
    fgraph <- .tdr_get_map_slot(.tdr.obj, "fuzzy.graph", sn)
    
    n_cells <- nrow(fgraph)
    
    # Apply the same weighted-voting logic as get.map():
    # 1. Expand sparse fgraph into (cell, landmark, weight) triplets
    # 2. Map landmarks to celltype labels
    # 3. Sum weights per (cell, label), compute confidence
    # 4. Assign label if confidence >= threshold, else "..low.confidence.."
    new_celltype_ids[[sn]] <-
      Matrix::summary(object = fgraph) |>
      dplyr::rename(cell = i,
                    landmark = j) |>
      dplyr::mutate(label = as.character(x = ct_ids[landmark])) |>
      dplyr::select(-landmark) |>
      collapse::fgroup_by(cell,
                          label,
                          sort = FALSE) |>
      collapse::fsum() |>
      collapse::fgroup_by(cell,
                          sort = FALSE) |>
      collapse::fmutate(confidence = x / collapse::fsum(x)) |>
      collapse::fungroup() |>
      collapse::fsubset(confidence >= .label.confidence) |>
      (\(x)
       {
         lbl <-
           rep(x = "..low.confidence..",
               times = n_cells)
         
         lbl[x$cell] <-
           x$label
         
         # Recover original cell names from clustering IDs (same ordering)
         cell_names <- names(.tdr_get_map_slot(.tdr.obj, "clustering.ids", sn))
         
         stats::setNames(object = lbl,
                         nm = cell_names)
       }
      )()
    
    # Write to on-disk cache if active
    if (is_cached) {
      cache_record <- .tdr_cache_write(
        object = new_celltype_ids[[sn]],
        cache_dir = cache$root,
        slot_name = "celltyping.ids",
        sample_name = sn
      )
      .tdr.obj@density$.cache$manifests$celltyping.ids[[sn]] <- cache_record
    }
  }
  
  # Store in-memory (NULL if on-disk caching, populated otherwise)
  if (!is_cached) {
    .tdr.obj@cellmap$celltype.ids <- new_celltype_ids
  }
  
  .tdr.obj
}

#' Recompute celltyping composition matrices from cell-level IDs
#'
#' Aggregates cell-level celltype assignments into samples x celltypes
#' count and percentage matrices, replacing
#' \code{@density$composition$celltyping$cell.count} and
#' \code{@density$composition$celltyping$cell.perc}.
#'
#' @param .tdr.obj A TDRObj with freshly relabeled celltype IDs.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.recompute_composition_celltyping <- function(.tdr.obj) {
  
  sample_names <- names(.tdr.obj@cells)
  
  ct_counts <- lapply(
    X = stats::setNames(sample_names, sample_names),
    FUN = function(sn) {
      ids <- .tdr_get_map_slot(.tdr.obj, "celltyping.ids", sn)
      tbl <- table(ids)
      stats::setNames(as.vector(tbl), names(tbl))
    }
  )
  
  # Determine which celltypes to exclude (the previously ignored population)
  .cl.ct.to.ign <- .tdr.obj@density$ignored
  
  .tdr.obj@density$composition$celltyping$cell.count <-
    ct_counts |>
    dplyr::bind_rows(.id = "sample") |>
    as.data.frame() |>
    (\(x)
     `rownames<-`(x = as.matrix(x = x[, colnames(x = x) != "sample"]),
                  value = x$sample)
    )() |>
    (\(x)
     x[, colnames(x = x) |>
         sort()]
    )() |>
    (\(x)
     x[, !(colnames(x = x) %in% .cl.ct.to.ign)]
    )()
  
  .tdr.obj@density$composition$celltyping$cell.count[
    is.na(x = .tdr.obj@density$composition$celltyping$cell.count)
  ] <- 0
  
  .tdr.obj@density$composition$celltyping$cell.perc <-
    (.tdr.obj@density$composition$celltyping$cell.count * 100) /
    Matrix::rowSums(x = .tdr.obj@density$composition$celltyping$cell.count)
  
  .tdr.obj
}

#' Invalidate stale trad celltyping fits in get.lm() results
#'
#' Sets \code{@results$lm[[model]]$trad$celltyping} to \code{NULL} for every
#' stored model and emits a warning so the user knows to re-run
#' \code{get.lm()}.
#'
#' @param .tdr.obj A TDRObj.
#' @return The modified \code{.tdr.obj}.
#' @keywords internal
.invalidate_trad_celltyping <- function(.tdr.obj) {
  
  model_names <- names(.tdr.obj@results$lm)
  
  for (mn in model_names) {
    if (!is.null(.tdr.obj@results$lm[[mn]]$trad$celltyping)) {
      .tdr.obj@results$lm[[mn]]$trad$celltyping <- NULL
      warning("Invalidated trad$celltyping fit for model '", mn,
              "'. Re-run get.lm() to recompute with updated celltype labels.",
              call. = FALSE)
    }
  }
  
  .tdr.obj
}

#' Warn about potentially stale pseudobulk / marker DE results
#'
#' Checks whether any \code{get.pbDE()} or \code{get.markerDE()} results
#' were computed with \code{.id.from = "celltyping"} and emits a warning
#' if so, since those results now reflect the old celltype labels.
#'
#' @param .tdr.obj A TDRObj.
#' @return The (unmodified) \code{.tdr.obj}, invisibly.
#' @keywords internal
.warn_stale_de_results <- function(.tdr.obj) {
  
  # Check pseudobulk DE results
  if (!is.null(.tdr.obj@results$pb)) {
    for (mn in names(.tdr.obj@results$pb)) {
      for (pn in names(.tdr.obj@results$pb[[mn]])) {
        res <- .tdr.obj@results$pb[[mn]][[pn]]
        if (!is.null(res$.id.from) && res$.id.from == "celltyping") {
          warning("Pseudobulk DE results for model '", mn,
                  "', population '", pn,
                  "' were computed with old celltyping labels and may be stale.",
                  call. = FALSE)
        }
      }
    }
  }
  
  # Check marker DE results
  if (!is.null(.tdr.obj@results$marker)) {
    for (mn in names(.tdr.obj@results$marker)) {
      res <- .tdr.obj@results$marker[[mn]]
      if (!is.null(res$.id.from) && res$.id.from == "celltyping") {
        warning("Marker DE results for model '", mn,
                "' were computed with old celltyping labels and may be stale.",
                call. = FALSE)
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
  
  # --- P2, P3: median expression & pheatmap (always recompute — cheap) ---
  .tdr.obj <- .recompute_celltyping_summaries(.tdr.obj)
  
  # --- P4–P7: cell-level IDs & composition (only if get.map() has run) ---
  if (!is.null(.tdr.obj@density$fdens)) {
    
    if (isTRUE(.verbose)) {
      message("-> get.map() results detected; refreshing cell-level celltype IDs and composition...")
    }
    
    .tdr.obj <- .relabel_cellmap(.tdr.obj, .verbose = .verbose)
    .tdr.obj <- .recompute_composition_celltyping(.tdr.obj)
    
    if (isTRUE(.verbose)) {
      message("-> Cell-level celltype IDs and composition matrices refreshed.")
    }
  }
  
  # --- P8: trad celltyping fit (only if get.lm() was run) ---
  if (!is.null(.tdr.obj@results$lm)) {
    .tdr.obj <- .invalidate_trad_celltyping(.tdr.obj)
  }
  
  # --- Warn about potentially stale pbDE / markerDE results ---
  .warn_stale_de_results(.tdr.obj)
  
  .tdr.obj
}

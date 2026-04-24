#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

#' Create a hierarchical subset TDRObj
#'
#' Given a parent \code{TDRObj} that has been processed through
#' \code{\link{get.map}}, creates a new \emph{child} \code{TDRObj} containing
#' only the cells matching the specified selection criteria.
#' The child references the parent's expression data \strong{without copying}:
#' for file-backed data, the same RDS paths are reused with an index filter;
#' for in-memory backends (Seurat, SCE, matrix, cytoset), the index vectors are
#' intersected.
#'
#' The child object is returned with empty analysis slots (assay, landmarks,
#' graphs, density, cellmap, results).
#' The user must re-run the full pipeline
#' (\code{get.landmarks} \eqn{\to} \code{get.graph} \eqn{\to}
#'  \code{get.map} \eqn{\to} \code{get.lm} \eqn{\to} \ldots)
#' on the child to obtain subset-specific results.
#'
#' Cell selection uses the same machinery as \code{\link{get.pbDE}}:
#' \itemize{
#'   \item \code{.id}: select cells whose per-cell label (from
#'     \code{@cellmap\$clustering\$ids} or \code{@cellmap\$celltyping\$ids})
#'     matches the specified identifiers.
#'   \item \code{.id.idx}: select cells via fuzzy confidence-thresholded
#'     voting from the specified landmark indices (using the UMAP fuzzy
#'     simplicial set stored in \code{@cellmap\$fuzzy.graphs}).
#' }
#'
#' Nested subsetting is supported: calling \code{get.subset()} on a child
#' object creates a grandchild that still references the original expression
#' data without additional copies.
#'
#' @section Session-scoped:
#' Subset objects are session-scoped.
#' Serialization via \code{\link[base]{saveRDS}} is not supported for child
#' TDRObjs that reference in-memory backends (seurat, sce, matrix, cyto).
#' For the files backend, the child can be saved as long as the parent's RDS
#' files remain on disk.
#'
#' @param x A \code{\linkS4class{TDRObj}}, Seurat, or SingleCellExperiment
#'   object that has been processed through \code{\link{get.map}}.
#' @param .source The raw data object for non-file backends (\code{NULL} for
#'   the files backend).
#'   When calling on a bare \code{TDRObj}, pass the original Seurat or SCE
#'   object here for seurat/sce backends.
#'   When calling via a dispatch wrapper (\code{get.subset.Seurat}), this is
#'   filled automatically.
#' @param .id Optional character vector of cluster or celltype labels to keep.
#'   Exactly one of \code{.id} or \code{.id.idx} must be provided.
#' @param .id.idx Optional integer vector of landmark indices.
#'   Cells are selected via fuzzy confidence-thresholded voting using the
#'   stored UMAP fuzzy graphs.
#' @param .id.from Character: \code{"clustering"} (default) or
#'   \code{"celltyping"}.  Source of labels when \code{.id} is used.
#' @param .label.confidence Numeric scalar in \code{[0,1]}.
#'   Minimum confidence for fuzzy label assignment.
#'   Defaults to the parent's stored threshold
#'   (\code{@config\$label.confidence}), or 0.5 if not set.
#' @param .prop.landmarks Numeric between 0 and 1 specifying proportion of
#'   cells to use as landmarks when the child pipeline runs
#'   \code{get.landmarks()}.  Default 0.1 (10 percent).
#' @param .min.cells.per.sample Integer.  Samples with fewer qualifying cells
#'   than this threshold are excluded from the child (default 10).
#' @param .verbose Logical: print progress? Default \code{TRUE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A new \code{\linkS4class{TDRObj}} with:
#' \describe{
#'   \item{\code{@cells}}{Modified per-sample cell references (no expression
#'     data copied).}
#'   \item{\code{@metadata}}{Filtered to retained samples, with updated
#'     \code{n.cells}.}
#'   \item{\code{@config}}{Inherited parameters (backend, assay.type,
#'     markers, harmony.var) plus provenance in
#'     \code{@config\$.subset.provenance}.}
#'   \item{All other slots}{Empty -- ready for \code{get.landmarks()}.}
#' }
#'
#' @seealso \code{\link{get.map}} (required predecessor),
#'   \code{\link{get.pbDE}} (uses the same \code{.id}/\code{.id.idx}
#'   selection semantics)
#'
#' @examples
#' \dontrun{
#' # After running the full pipeline on a parent population
#' parent <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph() |>
#'   celltyping(.celltyping.map = list(CD4T = c("1","2"), CD8T = "3")) |>
#'   get.map()
#'
#' # Create subset of CD8 T cells
#' child <- get.subset(parent, .id = "CD8T", .id.from = "celltyping")
#'
#' # Run full pipeline on the subset (new landmarks, graph, etc.)
#' child <- child |>
#'   get.landmarks() |> get.graph() |>
#'   celltyping(.celltyping.map = list(Tem = "1", Tn = "2")) |>
#'   get.map() |>
#'   get.lm(.design = design)
#'
#' # Nested subsetting
#' grandchild <- get.subset(child, .id = "Tem", .id.from = "celltyping")
#' grandchild <- grandchild |> get.landmarks() |> get.graph() |> get.map()
#' }
#'
#' @export
get.subset <- function(x, ...) UseMethod("get.subset")

#' @rdname get.subset
#' @export
get.subset.TDRObj <-
  function(
    x,
    .source = NULL,
    .id = NULL,
    .id.idx = NULL,
    .id.from = "clustering",
    .label.confidence = NULL,
    .prop.landmarks = 0.1,
    .min.cells.per.sample = 10,
    .verbose = TRUE,
    ...
  ) {
    .tdr.obj <- x

    # -----------------------------------------------------------------------
    # Input validation
    # -----------------------------------------------------------------------

    if (is.null(.id) && is.null(.id.idx)) {
      stop("Either .id or .id.idx must be provided.", call. = FALSE)
    }
    if (!is.null(.id) && !is.null(.id.idx)) {
      stop("Provide .id or .id.idx, not both.", call. = FALSE)
    }

    # Require get.map() to have been run
    has.fuzzy  <- !is.null(.tdr.obj@cellmap$fuzzy.graphs)
    has.labels <- !is.null(.tdr.obj@cellmap$clustering$ids) ||
                  !is.null(.tdr.obj@cellmap$celltyping$ids)

    if (!has.fuzzy && !has.labels) {
      stop("get.map() must be run on the parent TDRObj before subsetting.\n",
           "The cellmap (fuzzy.graphs, clustering/celltyping labels) is required ",
           "to resolve cell selections.",
           call. = FALSE)
    }

    if (!is.null(.id.idx) && !has.fuzzy) {
      stop("get.map() must be run on the parent TDRObj before subsetting ",
           "with .id.idx.\n",
           "@cellmap$fuzzy.graphs is required for fuzzy label transfer.",
           call. = FALSE)
    }

    .id.from <- match.arg(arg = .id.from,
                          choices = c("clustering", "celltyping"))

    # Default .label.confidence from parent
    if (is.null(.label.confidence)) {
      .label.confidence <- .tdr.obj@config$label.confidence
      if (is.null(.label.confidence)) .label.confidence <- 0.5
    }
    .tdr_validate_label_confidence(.label.confidence)

    if (!is.numeric(.min.cells.per.sample) ||
        length(.min.cells.per.sample) != 1 ||
        .min.cells.per.sample < 1) {
      stop(".min.cells.per.sample must be a positive integer.",
           call. = FALSE)
    }

    # -----------------------------------------------------------------------
    # Resolve cell indices (reuses get.pbDE machinery)
    # -----------------------------------------------------------------------

    cell.idx.list <- .tdr_resolve_cell_idx(
      .tdr.obj,
      .id = .id,
      .id.idx = .id.idx,
      .id.from = .id.from,
      .label.confidence = .label.confidence
    )

    # -----------------------------------------------------------------------
    # Determine backend
    # -----------------------------------------------------------------------

    backend <- .tdr.obj@config$backend
    if (is.null(backend)) backend <- "files"

    # -----------------------------------------------------------------------
    # Build child @cells (zero-copy)
    # -----------------------------------------------------------------------

    new.cells    <- list()
    new.n.cells  <- integer(0)

    for (sn in names(.tdr.obj@cells)) {

      idx <- cell.idx.list[[sn]]
      if (is.null(idx) || length(idx) < .min.cells.per.sample) next

      if (backend == "files") {
        # Per Q5(b): store list(path, idx) in @cells
        parent.entry <- .tdr.obj@cells[[sn]]
        if (is.list(parent.entry)) {
          # Nested subset: compose with existing filter
          new.cells[[sn]] <- list(
            path = parent.entry$path,
            idx  = parent.entry$idx[idx]
          )
        } else {
          # First-level subset from a plain path string
          new.cells[[sn]] <- list(
            path = parent.entry,
            idx  = idx
          )
        }
      } else {
        # Index-based backends (seurat, sce, matrix, cyto):
        # Per Q9(a): replace @cells indices directly
        parent.idx <- .tdr.obj@cells[[sn]]
        new.cells[[sn]] <- parent.idx[idx]
      }

      new.n.cells[sn] <- length(idx)
    }

    if (length(new.cells) == 0) {
      stop("No samples have >= ", .min.cells.per.sample,
           " qualifying cells after subsetting.",
           call. = FALSE)
    }

    # -----------------------------------------------------------------------
    # Build child @metadata
    # -----------------------------------------------------------------------

    new.meta <- .tdr.obj@metadata[names(new.cells), , drop = FALSE]
    new.meta$n.cells      <- new.n.cells[rownames(new.meta)]
    new.meta$log10.n.cells <- log10(new.meta$n.cells)
    # n.perSample will be set from config$sampling after construction
    new.meta$n.perSample  <- NULL

    # -----------------------------------------------------------------------
    # Build child @config
    # -----------------------------------------------------------------------

    new.config <- list(
      assay.type = .tdr.obj@config$assay.type,
      backend    = .tdr.obj@config$backend,
      markers    = .tdr.obj@config$markers,
      n.threads  = .tdr.obj@config$n.threads
    )

    # Inherit Harmony variable (per Q2-a: child re-runs Harmony on its own
    # landmarks, but needs the batch variable name)
    if (!is.null(.tdr.obj@integration$harmony.var)) {
      new.config$harmony.var <- .tdr.obj@integration$harmony.var
    }

    # Inherit source environment for backends that store it
    if (!is.null(.tdr.obj@config$source.env)) {
      new.config$source.env <- .tdr.obj@config$source.env
    }

    # For seurat/sce backends: store .source so descendant pipeline steps
    # can access expression data without the original container
    if (backend %in% c("seurat", "sce") && !is.null(.source)) {
      if (is.null(new.config$source.env)) {
        new.config$source.env <- new.env(parent = emptyenv())
      }
      new.config$source.env$source <- .source
    }

    # Seurat/SCE assay/layer configuration
    if (!is.null(.tdr.obj@config$source.assay)) {
      new.config$source.assay <- .tdr.obj@config$source.assay
    }
    if (!is.null(.tdr.obj@config$source.layer)) {
      new.config$source.layer <- .tdr.obj@config$source.layer
    }

    # Validate .prop.landmarks
    if (.prop.landmarks < 0 || .prop.landmarks > 1) {
      stop(".prop.landmarks must be between 0 and 1 (e.g., 0.1 for 10% of cells).\n",
           "Current value: ", .prop.landmarks,
           call. = FALSE)
    }

    # Sampling configuration needed by get.landmarks
    n.cells.vec <- new.n.cells[names(new.cells)]
    target.lm.n <- pmin(sum(n.cells.vec) * .prop.landmarks, 5e3)
    n.perSample <- pmin(
      ceiling(n.cells.vec * .prop.landmarks),
      ceiling(target.lm.n / length(new.cells))
    )
    new.config$sampling <- list(
      n.cells      = n.cells.vec,
      target.lm.n  = target.lm.n,
      n.perSample  = n.perSample
    )

    # Key vector (maps each future landmark to its sample index)
    new.config$key <-
      seq_along(new.cells) |>
      rep(times = n.perSample) |>
      stats::setNames(nm = names(new.cells)[
        seq_along(new.cells) |> rep(times = n.perSample)
      ])

    # -----------------------------------------------------------------------
    # Provenance
    # -----------------------------------------------------------------------

    provenance.entry <- list(
      selection.type    = if (!is.null(.id)) "id" else "id.idx",
      selection.values  = if (!is.null(.id)) .id else .id.idx,
      selection.from    = .id.from,
      label.confidence  = .label.confidence,
      parent.n.landmarks = nrow(.tdr.obj@assay$expr),
      parent.n.samples  = length(.tdr.obj@cells),
      created.at        = Sys.time()
    )

    # Append to existing chain for nested subsetting
    new.config$.subset.provenance <-
      c(.tdr.obj@config$.subset.provenance,
        list(provenance.entry))

    # -----------------------------------------------------------------------
    # Construct child TDRObj
    # -----------------------------------------------------------------------

    # Set n.perSample in metadata (expected by get.landmarks)
    new.meta$n.perSample <- n.perSample[rownames(new.meta)]

    child <- TDRObj(
      cells    = new.cells,
      metadata = new.meta,
      config   = new.config
    )

    if (isTRUE(.verbose)) {
      message(
        sprintf("\nSubset TDRObj created: %d samples (from %d parent samples)",
                length(new.cells), length(.tdr.obj@cells))
      )
      message(
        sprintf("  Cell counts: min=%d, median=%d, max=%d",
                min(new.n.cells), as.integer(stats::median(new.n.cells)),
                max(new.n.cells))
      )
      message("  Run get.landmarks() to start the subset analysis pipeline.")
    }

    return(child)
  }

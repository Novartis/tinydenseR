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
#' \code{RunTDR} is a convenience wrapper around the explicit step-by-step
#' workflow (\code{get.meta → get.cells → setup.tdr.obj → get.landmarks →
#' get.graph → get.map → get.embedding}).
#' For equivalent inputs, seed, and arguments the two paths produce
#' \strong{numerically identical} results.
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
#' Executes \code{get.landmarks → get.graph → celltyping → get.map → get.embedding}
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

  # --- Argument routing ---
  dots <- list(...)

  .lm.formals  <- setdiff(names(formals(get.landmarks.TDRObj)), c("x", "..."))
  .gr.formals  <- setdiff(names(formals(get.graph.TDRObj)), c("x", "..."))
  .map.formals <- setdiff(names(formals(get.map.TDRObj)), c("x", "..."))
  .emb.formals <- setdiff(names(formals(get.embedding.TDRObj)), c("x", "..."))

  lm_args  <- dots[names(dots) %in% .lm.formals]
  gr_args  <- dots[names(dots) %in% .gr.formals]
  map_args <- dots[names(dots) %in% .map.formals]
  emb_args <- dots[names(dots) %in% .emb.formals]

  all_known <- Reduce(union, list(.lm.formals, .gr.formals, .map.formals, .emb.formals))
  orphans <- setdiff(names(dots), all_known)
  if (length(orphans) > 0) {
    warning("Unknown arguments will be ignored: ",
            paste(orphans, collapse = ", "))
  }

  # --- Pipeline ---
  x <- do.call(get.landmarks, c(list(x, .source = NULL, .seed = .seed, .verbose = .verbose), lm_args))
  x <- do.call(get.graph, c(list(x, .seed = .seed, .verbose = .verbose), gr_args))

  if (!is.null(eff_ct)) {
    x <- celltyping(x, .celltyping.map = eff_ct, .verbose = .verbose)
  }

  x <- do.call(get.map, c(list(x, .source = NULL, .seed = .seed, .verbose = .verbose), map_args))
  x <- do.call(get.embedding, c(list(x, .seed = .seed, .verbose = .verbose), emb_args))

  return(x)
}

#' @describeIn RunTDR Run the pipeline directly on a Seurat object
#'
#' Builds a \code{\linkS4class{TDRObj}} from a Seurat object and executes
#' the full pipeline. The finished TDRObj is stored in
#' \code{SeuratObject::Misc(x, slot = "tdr.obj")}.
#'
#' All categorical cell-level columns in \code{x@@meta.data} are
#' automatically imported as named celltyping solutions via
#' \code{\link{import_cell_annotations}}.  The \code{.celltype.vec}
#' column (if specified) is set as the active annotation.
#'
#' @param x A Seurat object.
#' @param .sample.var Character(1). Column name in \code{x@@meta.data}
#'   identifying sample membership.
#' @param .assay Character(1). Name of the Seurat assay to use.
#' @param .layer Character(1). Layer within the assay (e.g. \code{"counts"}).
#' @param .assay.type Character. \code{"RNA"} or \code{"cyto"}.
#' @param .harmony.var Character vector of batch variable column names
#'   in sample-level metadata, or \code{NULL}.
#' @param .markers Character vector of marker names (required for cyto).
#' @param .celltype.vec Character(1). Column name in \code{x@@meta.data}
#'   containing per-cell type labels, or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for a sample to be
#'   included.
#' @param .verbose Logical. Print progress messages.
#' @param .seed Integer. Random seed.
#' @param .prop.landmarks Numeric in (0, 1]. Proportion of cells as landmarks.
#' @param .n.threads Integer. Number of threads.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return The Seurat object \code{x} with the TDRObj stored in
#'   \code{Misc(x, slot = "tdr.obj")}.
#'
#' @export
RunTDR.Seurat <- function(x,
                          .sample.var,
                          .assay = "RNA",
                          .layer = "counts",
                          .assay.type = "RNA",
                          .harmony.var = NULL,
                          .markers = NULL,
                          .celltype.vec = NULL,
                          .min.cells.per.sample = 10,
                          .verbose = TRUE,
                          .seed = 123,
                          .prop.landmarks = 0.1,
                          .n.threads = if (is.hpc()) {
                            max(RhpcBLASctl::blas_get_num_procs(),
                                RhpcBLASctl::omp_get_num_procs(),
                                RhpcBLASctl::omp_get_max_threads(),
                                na.rm = TRUE)
                          } else {
                            parallel::detectCores(logical = TRUE)
                          },
                          ...) {

  # --- Input validation ---
  if (!inherits(x, "Seurat")) {
    stop("x must be a Seurat object.")
  }

  if (!is.character(.sample.var) || length(.sample.var) != 1) {
    stop(".sample.var must be a single character string.")
  }

  if (!(.sample.var %in% colnames(x@meta.data))) {
    stop(".sample.var '", .sample.var, "' not found in x@meta.data.")
  }

  if (!(.assay %in% names(x@assays))) {
    stop(".assay '", .assay, "' not found in x@assays.")
  }

  # --- Extract sample-level metadata ---
  .meta <- get.meta(.obj = x, .sample.var = .sample.var,
                    .verbose = .verbose)

  # --- Build @cells as named list of sorted integer index vectors ---
  sample_ids <- x@meta.data[[.sample.var]]
  counts_tbl <- table(sample_ids)
  valid <- names(counts_tbl)[
    counts_tbl >= .min.cells.per.sample &
    names(counts_tbl) %in% rownames(.meta)
  ]
  .cells <- lapply(
    stats::setNames(valid, valid),
    function(s) sort(which(sample_ids == s))
  )

  # Align sample ordering with get.meta() for parity with explicit workflow.
  # table() sorts names alphabetically, but get.meta() preserves first-
  # occurrence order.  The explicit workflow feeds get.meta()-ordered .cells

  # to setup.tdr.obj → get.landmarks, so landmark sampling (which depends on
  # processing order through the RNG stream) must see the same ordering here.
  .meta_order <- intersect(rownames(.meta), names(.cells))
  .cells <- .cells[.meta_order]
  .meta  <- .meta[.meta_order, , drop = FALSE]

  # --- .celltype.vec handling ---
  ct_vec <- NULL
  if (!is.null(.celltype.vec)) {
    if (!is.character(.celltype.vec) || length(.celltype.vec) != 1) {
      stop(".celltype.vec must be a single character string ",
           "(column name in x@meta.data).")
    }
    if (!(.celltype.vec %in% colnames(x@meta.data))) {
      stop(".celltype.vec '", .celltype.vec,
           "' not found in x@meta.data.")
    }
    valid_cells <- unlist(.cells)
    ct_vec <- stats::setNames(
      as.character(x@meta.data[[.celltype.vec]][valid_cells]),
      rownames(x@meta.data)[valid_cells]
    )
  }

  # --- Build TDRObj ---
  tdr.obj <- .setup_tdr_from_seurat(
    .cells = .cells,
    .meta = .meta,
    .assay = .assay,
    .layer = .layer,
    .assay.type = .assay.type,
    .markers = .markers,
    .harmony.var = .harmony.var,
    .celltype.vec = ct_vec,
    .prop.landmarks = .prop.landmarks,
    .seed = .seed,
    .n.threads = .n.threads,
    .verbose = .verbose
  )

  # --- Argument routing ---
  dots <- list(...)

  .lm.formals  <- setdiff(names(formals(get.landmarks.TDRObj)), c("x", "..."))
  .gr.formals  <- setdiff(names(formals(get.graph.TDRObj)), c("x", "..."))
  .map.formals <- setdiff(names(formals(get.map.TDRObj)), c("x", "..."))
  .emb.formals <- setdiff(names(formals(get.embedding.TDRObj)), c("x", "..."))

  lm_args  <- dots[names(dots) %in% .lm.formals]
  gr_args  <- dots[names(dots) %in% .gr.formals]
  map_args <- dots[names(dots) %in% .map.formals]
  emb_args <- dots[names(dots) %in% .emb.formals]

  all_known <- Reduce(union, list(.lm.formals, .gr.formals, .map.formals, .emb.formals))
  orphans <- setdiff(names(dots), all_known)
  if (length(orphans) > 0) {
    warning("Unknown arguments will be ignored: ",
            paste(orphans, collapse = ", "))
  }

  # --- Pipeline ---
  tdr.obj <- do.call(get.landmarks.TDRObj, c(list(tdr.obj, .source = x, .seed = .seed, .verbose = .verbose), lm_args))
  tdr.obj <- do.call(get.graph.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), gr_args))

  # --- Import all categorical cell-level annotations ---
  valid_cells <- unlist(.cells)
  .cell.meta.full <- x@meta.data[valid_cells, , drop = FALSE]
  tdr.obj <- import_cell_annotations(tdr.obj,
                                     .cell.meta = .cell.meta.full,
                                     .sample.var = .sample.var,
                                     .celltype.vec = .celltype.vec,
                                     .verbose = .verbose)

  tdr.obj <- do.call(get.map.TDRObj, c(list(tdr.obj, .source = x, .seed = .seed, .verbose = .verbose), map_args))
  tdr.obj <- do.call(get.embedding.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), emb_args))

  # --- Store TDRObj in Seurat Misc slot ---
  SeuratObject::Misc(x, slot = "tdr.obj") <- tdr.obj

  return(x)
}


# ──────────────────────────────────────────────────────────────────────
# Internal helper: build TDRObj from Seurat-derived inputs
# ──────────────────────────────────────────────────────────────────────

#' Build a TDRObj from Seurat-derived cell lists and metadata
#'
#' @keywords internal
#' @noRd
.setup_tdr_from_seurat <- function(.cells,
                                   .meta,
                                   .assay,
                                   .layer,
                                   .assay.type,
                                   .markers,
                                   .harmony.var,
                                   .celltype.vec,
                                   .prop.landmarks,
                                   .seed,
                                   .n.threads,
                                   .verbose) {

  .assay.type <- match.arg(arg = .assay.type,
                           choices = c("cyto", "RNA"))

  # --- Harmony var validation ---
  if (!is.null(.harmony.var)) {
    if (!inherits(.harmony.var, "character")) {
      stop(".harmony.var must be a character vector.")
    }
    if (!all(.harmony.var %in% colnames(.meta))) {
      stop("Variables not found in metadata: ",
           paste(.harmony.var[!(.harmony.var %in% colnames(.meta))],
                 collapse = ", "),
           "\nCheck column names in .meta with colnames(.meta).")
    }
  }

  # --- Markers validation ---
  if (!is.null(.markers)) {
    if (.assay.type == "RNA") {
      stop(".markers argument only applies to cytometry data.\n",
           "For RNA data, feature selection uses highly variable genes (HVG) automatically.")
    } else if (length(.markers) < 3) {
      stop(".markers must contain at least 3 markers for meaningful dimensionality reduction.")
    }
  }

  if (.assay.type == "cyto" && is.null(.markers)) {
    stop("For cyto assay.type with Seurat backend, .markers must be provided.")
  }

  # --- Prop landmarks validation ---
  if ((.prop.landmarks < 0) | (.prop.landmarks > 1)) {
    stop(".prop.landmarks must be between 0 and 1 (e.g., 0.1 for 10% of cells).\n",
         "Current value: ", .prop.landmarks)
  }

  # --- Create TDRObj ---
  .tdr.obj <- TDRObj(
    config = list(
      key = NULL,
      sampling = NULL,
      assay.type = .assay.type,
      markers = NULL,
      n.threads = .n.threads
    ),
    integration = list(
      harmony.var = NULL,
      harmony.obj = NULL
    )
  )

  .tdr.obj@cells <- .cells

  # --- n.cells from integer index vector lengths (no readRDS needed) ---
  n.cells <- lengths(.cells)

  .tdr.obj@config$sampling$n.cells <- n.cells

  # Quality check: warn if sample sizes are highly imbalanced (>10-fold difference)
  if ((max(.tdr.obj@config$sampling$n.cells) /
       min(.tdr.obj@config$sampling$n.cells)) > 10) {

    warning("Sample size imbalance detected: largest/smallest ratio > 10.\n",
            "Smallest sample has ", min(.tdr.obj@config$sampling$n.cells),
            " cells.\n",
            "Consider removing low-quality samples.")

    if (any(.tdr.obj@config$sampling$n.cells < 1000)) {
      warning("Large variation in sample sizes detected. ",
              "For cytometry, samples with <1000 cells may be unreliable.")
    }
  }

  # Calculate target number of landmarks: prop of total cells, capped at 5000
  .tdr.obj@config$sampling$target.lm.n <-
    pmin(sum(.tdr.obj@config$sampling$n.cells) * .prop.landmarks,
         5e3)

  # Allocate landmarks per sample: proportional to sample size, but capped
  .tdr.obj@config$sampling$n.perSample <-
    pmin(ceiling(x = .tdr.obj@config$sampling$n.cells * .prop.landmarks),
         ceiling(x = .tdr.obj@config$sampling$target.lm.n /
                     length(x = .tdr.obj@cells)))

  # Create key vector: maps each future landmark to its sample
  .tdr.obj@config$key <-
    seq_along(along.with = .tdr.obj@cells) |>
    rep(times = .tdr.obj@config$sampling$n.perSample) |>
    (\(x)
     stats::setNames(object = x,
                     nm = names(.tdr.obj@cells)[x])
    )()

  .tdr.obj@metadata <- .meta

  .tdr.obj@metadata$n.perSample <-
    .tdr.obj@config$sampling$n.perSample

  .tdr.obj@metadata$n.cells <-
    .tdr.obj@config$sampling$n.cells

  .tdr.obj@metadata$log10.n.cells <-
    log10(x = .tdr.obj@config$sampling$n.cells)

  # --- Markers ---
  if (.assay.type == "cyto") {
    .tdr.obj@config$markers <- .markers
  }

  # --- Harmony ---
  if (!is.null(.harmony.var)) {
    .tdr.obj@integration$harmony.var <- .harmony.var
  }

  # --- Seurat backend config ---
  .tdr.obj@config$backend <- "seurat"
  .tdr.obj@config$source.assay <- .assay
  .tdr.obj@config$source.layer <- .layer

  # --- Cell type vector ---
  if (!is.null(.celltype.vec)) {
    .tdr.obj@config$celltype.vec <- .celltype.vec
  }

  return(.tdr.obj)
}


#' @describeIn RunTDR Run the pipeline directly on a SingleCellExperiment
#'
#' Builds a \code{\linkS4class{TDRObj}} from a \code{SingleCellExperiment}
#' and executes the full pipeline. The finished TDRObj is stored in
#' \code{S4Vectors::metadata(x)$tdr.obj}.
#'
#' If the assay is a \code{DelayedMatrix} (e.g., HDF5-backed), it is
#' converted to a BPCells on-disk \code{IterableMatrix} and routed
#' through \code{.run_tdr_matrix()} for efficient lazy access. If the
#' assay is an in-memory matrix (e.g., \code{dgCMatrix}), the existing
#' SCE backend path is used.
#'
#' @param x A \code{SingleCellExperiment} object.
#' @param .sample.var Character(1). Column in \code{colData(x)} identifying
#'   sample membership.
#' @param .assay Character(1). Name of the assay in \code{assayNames(x)}.
#' @param .assay.type Character. \code{"RNA"} or \code{"cyto"}.
#' @param .bpcells.dir Character(1) or \code{NULL}. Directory path for
#'   the BPCells on-disk matrix (used only when the assay is a
#'   \code{DelayedMatrix}). If \code{NULL} (default), uses a temporary
#'   directory (\code{tempdir()}) that is cleaned up on session end.
#'   If a path is given and already contains a valid BPCells matrix,
#'   the conversion step is skipped (cache hit).
#' @param .harmony.var Character vector of batch variable column names
#'   in sample-level metadata, or \code{NULL}.
#' @param .markers Character vector of marker names (required for cyto).
#' @param .celltype.vec Character(1). Column name in \code{colData(x)}
#'   containing per-cell type labels, or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for a sample to be
#'   included.
#' @param .verbose Logical. Print progress messages.
#' @param .seed Integer. Random seed.
#' @param .prop.landmarks Numeric in (0, 1]. Proportion of cells as landmarks.
#' @param .n.threads Integer. Number of threads.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return The \code{SingleCellExperiment} \code{x} with the TDRObj stored
#'   in \code{S4Vectors::metadata(x)$tdr.obj}.
#'
#' @export
RunTDR.SingleCellExperiment <- function(x,
                                        .sample.var,
                                        .assay = "counts",
                                        .assay.type = "RNA",
                                        .bpcells.dir = NULL,
                                        .harmony.var = NULL,
                                        .markers = NULL,
                                        .celltype.vec = NULL,
                                        .min.cells.per.sample = 10,
                                        .verbose = TRUE,
                                        .seed = 123,
                                        .prop.landmarks = 0.1,
                                        .n.threads = if (is.hpc()) {
                                          max(RhpcBLASctl::blas_get_num_procs(),
                                              RhpcBLASctl::omp_get_num_procs(),
                                              RhpcBLASctl::omp_get_max_threads(),
                                              na.rm = TRUE)
                                        } else {
                                          parallel::detectCores(logical = TRUE)
                                        },
                                        ...) {

  # --- Input validation ---

  if (!inherits(x, "SingleCellExperiment")) {
    stop("x must be a SingleCellExperiment object.")
  }

  if (!is.character(.sample.var) || length(.sample.var) != 1) {
    stop(".sample.var must be a single character string.")
  }

  if (!(.sample.var %in% colnames(SummarizedExperiment::colData(x)))) {
    stop(".sample.var '", .sample.var,
         "' not found in colData(x).")
  }

  if (!(.assay %in% SummarizedExperiment::assayNames(x))) {
    stop(".assay '", .assay, "' not found in assayNames(x).")
  }

  # --- Check if assay is a DelayedMatrix ---
  is_delayed <- methods::is(
    SummarizedExperiment::assay(x, .assay), "DelayedMatrix")

  if (is_delayed) {

    # ══════════════════════════════════════════════════════════════════
    # BPCells path for DelayedMatrix-backed SCE
    # ══════════════════════════════════════════════════════════════════

    if (!requireNamespace("BPCells", quietly = TRUE)) {
      stop("Package 'BPCells' is required for DelayedMatrix SCE support. ",
           "Install it with: BiocManager::install('BPCells')", call. = FALSE)
    }

    if (isTRUE(.verbose)) {
      cat("DelayedMatrix assay detected; converting to BPCells...\n")
    }

    # Extract the DelayedMatrix
    mat <- SummarizedExperiment::assay(x, .assay)

    # Convert to BPCells on-disk format
    bp_mat <- .delayed_to_bpcells(mat, .bpcells.dir, .verbose)

    # Preserve dimnames if lost during conversion
    if (!is.null(rownames(mat)) && is.null(rownames(bp_mat))) {
      rownames(bp_mat) <- rownames(mat)
    }
    if (!is.null(colnames(mat)) && is.null(colnames(bp_mat))) {
      colnames(bp_mat) <- colnames(mat)
    }

    # Extract cell metadata: colData is a DataFrame, convert to data.frame
    cell_meta <- as.data.frame(SummarizedExperiment::colData(x))

    # Delegate to .run_tdr_matrix (the proven IterableMatrix path)
    tdr.obj <- .run_tdr_matrix(
      x               = bp_mat,
      .cell.meta      = cell_meta,
      .sample.var     = .sample.var,
      .assay.type     = .assay.type,
      .harmony.var    = .harmony.var,
      .markers        = .markers,
      .celltype.vec   = .celltype.vec,
      .min.cells.per.sample = .min.cells.per.sample,
      .verbose        = .verbose,
      .seed           = .seed,
      .prop.landmarks = .prop.landmarks,
      .n.threads      = .n.threads,
      ...
    )

    # Store TDRObj in SCE metadata
    S4Vectors::metadata(x)$tdr.obj <- tdr.obj
    return(x)

  } else {

    # ══════════════════════════════════════════════════════════════════
    # Existing in-memory path for dgCMatrix / non-delayed assays
    # ══════════════════════════════════════════════════════════════════

    # --- Extract sample-level metadata ---
    .meta <- get.meta(.obj = x, .sample.var = .sample.var,
                      .verbose = .verbose)

    # --- Build @cells as named list of sorted integer index vectors ---
    sample_ids <- SummarizedExperiment::colData(x)[[.sample.var]]
    counts_tbl <- table(sample_ids)
    valid <- names(counts_tbl)[
      counts_tbl >= .min.cells.per.sample &
      names(counts_tbl) %in% rownames(.meta)
    ]
    .cells <- lapply(
      stats::setNames(valid, valid),
      function(s) sort(which(sample_ids == s))
    )

    # Align sample ordering with get.meta() for parity with explicit workflow
    .meta_order <- intersect(rownames(.meta), names(.cells))
    .cells <- .cells[.meta_order]
    .meta  <- .meta[.meta_order, , drop = FALSE]

    # --- .celltype.vec handling ---
    ct_vec <- NULL
    if (!is.null(.celltype.vec)) {
      if (!is.character(.celltype.vec) || length(.celltype.vec) != 1) {
        stop(".celltype.vec must be a single character string ",
             "(column name in colData(x)).")
      }
      if (!(.celltype.vec %in%
            colnames(SummarizedExperiment::colData(x)))) {
        stop(".celltype.vec '", .celltype.vec,
             "' not found in colData(x).")
      }
      valid_cells <- unlist(.cells)
      ct_vec <- stats::setNames(
        as.character(
          SummarizedExperiment::colData(x)[[.celltype.vec]][valid_cells]),
        colnames(x)[valid_cells]
      )
    }

    # --- Build TDRObj ---
    tdr.obj <- .setup_tdr_from_sce(
      .cells = .cells,
      .meta = .meta,
      .assay = .assay,
      .assay.type = .assay.type,
      .markers = .markers,
      .harmony.var = .harmony.var,
      .celltype.vec = ct_vec,
      .prop.landmarks = .prop.landmarks,
      .seed = .seed,
      .n.threads = .n.threads,
      .verbose = .verbose
    )

    # --- Argument routing ---
    dots <- list(...)

    .lm.formals  <- setdiff(names(formals(get.landmarks.TDRObj)), c("x", "..."))
    .gr.formals  <- setdiff(names(formals(get.graph.TDRObj)), c("x", "..."))
    .map.formals <- setdiff(names(formals(get.map.TDRObj)), c("x", "..."))
    .emb.formals <- setdiff(names(formals(get.embedding.TDRObj)), c("x", "..."))

    lm_args  <- dots[names(dots) %in% .lm.formals]
    gr_args  <- dots[names(dots) %in% .gr.formals]
    map_args <- dots[names(dots) %in% .map.formals]
    emb_args <- dots[names(dots) %in% .emb.formals]

    all_known <- Reduce(union, list(.lm.formals, .gr.formals, .map.formals, .emb.formals))
    orphans <- setdiff(names(dots), all_known)
    if (length(orphans) > 0) {
      warning("Unknown arguments will be ignored: ",
              paste(orphans, collapse = ", "))
    }

    # --- Pipeline ---
    tdr.obj <- do.call(get.landmarks.TDRObj, c(list(tdr.obj, .source = x, .seed = .seed, .verbose = .verbose), lm_args))
    tdr.obj <- do.call(get.graph.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), gr_args))

    # --- Import all categorical cell-level annotations ---
    valid_cells <- unlist(.cells)
    .cell.meta.full <- as.data.frame(
      SummarizedExperiment::colData(x)
    )[valid_cells, , drop = FALSE]
    tdr.obj <- import_cell_annotations(tdr.obj,
                                       .cell.meta = .cell.meta.full,
                                       .sample.var = .sample.var,
                                       .celltype.vec = .celltype.vec,
                                       .verbose = .verbose)

    tdr.obj <- do.call(get.map.TDRObj, c(list(tdr.obj, .source = x, .seed = .seed, .verbose = .verbose), map_args))
    tdr.obj <- do.call(get.embedding.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), emb_args))

    # --- Store TDRObj in SCE metadata ---
    S4Vectors::metadata(x)$tdr.obj <- tdr.obj

    return(x)
  }
}


# ──────────────────────────────────────────────────────────────────────
# Internal helper: build TDRObj from SCE-derived inputs
# ──────────────────────────────────────────────────────────────────────

#' Build a TDRObj from SingleCellExperiment-derived cell lists and metadata
#'
#' @keywords internal
#' @noRd
.setup_tdr_from_sce <- function(.cells,
                                .meta,
                                .assay,
                                .assay.type,
                                .markers,
                                .harmony.var,
                                .celltype.vec,
                                .prop.landmarks,
                                .seed,
                                .n.threads,
                                .verbose) {

  .assay.type <- match.arg(arg = .assay.type,
                           choices = c("cyto", "RNA"))

  # --- Harmony var validation ---
  if (!is.null(.harmony.var)) {
    if (!inherits(.harmony.var, "character")) {
      stop(".harmony.var must be a character vector.")
    }
    if (!all(.harmony.var %in% colnames(.meta))) {
      stop("Variables not found in metadata: ",
           paste(.harmony.var[!(.harmony.var %in% colnames(.meta))],
                 collapse = ", "),
           "\nCheck column names in .meta with colnames(.meta).")
    }
  }

  # --- Markers validation ---
  if (!is.null(.markers)) {
    if (.assay.type == "RNA") {
      stop(".markers argument only applies to cytometry data.\n",
           "For RNA data, feature selection uses highly variable genes (HVG) automatically.")
    } else if (length(.markers) < 3) {
      stop(".markers must contain at least 3 markers for meaningful dimensionality reduction.")
    }
  }

  if (.assay.type == "cyto" && is.null(.markers)) {
    stop("For cyto assay.type with SCE backend, .markers must be provided.")
  }

  # --- Prop landmarks validation ---
  if ((.prop.landmarks < 0) | (.prop.landmarks > 1)) {
    stop(".prop.landmarks must be between 0 and 1 (e.g., 0.1 for 10% of cells).\n",
         "Current value: ", .prop.landmarks)
  }

  # --- Create TDRObj ---
  .tdr.obj <- TDRObj(
    config = list(
      key = NULL,
      sampling = NULL,
      assay.type = .assay.type,
      markers = NULL,
      n.threads = .n.threads
    ),
    integration = list(
      harmony.var = NULL,
      harmony.obj = NULL
    )
  )

  .tdr.obj@cells <- .cells

  # --- n.cells from integer index vector lengths ---
  n.cells <- lengths(.cells)

  .tdr.obj@config$sampling$n.cells <- n.cells

  # Quality check: warn if sample sizes are highly imbalanced (>10-fold difference)
  if ((max(.tdr.obj@config$sampling$n.cells) /
       min(.tdr.obj@config$sampling$n.cells)) > 10) {

    warning("Sample size imbalance detected: largest/smallest ratio > 10.\n",
            "Smallest sample has ", min(.tdr.obj@config$sampling$n.cells),
            " cells.\n",
            "Consider removing low-quality samples.")

    if (any(.tdr.obj@config$sampling$n.cells < 1000)) {
      warning("Large variation in sample sizes detected. ",
              "For cytometry, samples with <1000 cells may be unreliable.")
    }
  }

  # Calculate target number of landmarks: prop of total cells, capped at 5000
  .tdr.obj@config$sampling$target.lm.n <-
    pmin(sum(.tdr.obj@config$sampling$n.cells) * .prop.landmarks,
         5e3)

  # Allocate landmarks per sample: proportional to sample size, but capped
  .tdr.obj@config$sampling$n.perSample <-
    pmin(ceiling(x = .tdr.obj@config$sampling$n.cells * .prop.landmarks),
         ceiling(x = .tdr.obj@config$sampling$target.lm.n /
                     length(x = .tdr.obj@cells)))

  # Create key vector: maps each future landmark to its sample
  .tdr.obj@config$key <-
    seq_along(along.with = .tdr.obj@cells) |>
    rep(times = .tdr.obj@config$sampling$n.perSample) |>
    (\(x)
     stats::setNames(object = x,
                     nm = names(.tdr.obj@cells)[x])
    )()

  .tdr.obj@metadata <- .meta

  .tdr.obj@metadata$n.perSample <-
    .tdr.obj@config$sampling$n.perSample

  .tdr.obj@metadata$n.cells <-
    .tdr.obj@config$sampling$n.cells

  .tdr.obj@metadata$log10.n.cells <-
    log10(x = .tdr.obj@config$sampling$n.cells)

  # --- Markers ---
  if (.assay.type == "cyto") {
    .tdr.obj@config$markers <- .markers
  }

  # --- Harmony ---
  if (!is.null(.harmony.var)) {
    .tdr.obj@integration$harmony.var <- .harmony.var
  }

  # --- SCE backend config ---
  .tdr.obj@config$backend <- "sce"
  .tdr.obj@config$source.assay <- .assay

  # --- Cell type vector ---
  if (!is.null(.celltype.vec)) {
    .tdr.obj@config$celltype.vec <- .celltype.vec
  }

  return(.tdr.obj)
}


# ──────────────────────────────────────────────────────────────────────
# Internal: rhdf5-based H5AD metadata readers
# ──────────────────────────────────────────────────────────────────────

#' Detect AnnData format version from root HDF5 attributes
#'
#' @param path Character(1). Path to the h5ad file.
#' @return A list with \code{modern} (logical) and
#'   \code{encoding_version} (character or NA).
#' @keywords internal
#' @noRd
.h5ad_detect_version <- function(path) {
  attrs <- rhdf5::h5readAttributes(path, "/")
  rhdf5::H5close()
  if ("encoding-type" %in% names(attrs)) {
    return(list(modern = TRUE,
                encoding_version = attrs[["encoding-version"]]))
  }
  list(modern = FALSE, encoding_version = NA_character_)
}

#' Read obs metadata from an h5ad file as a data.frame via rhdf5
#'
#' Handles both pre-0.8.0 (flat \code{__categories}) and modern (nested
#' group) AnnData categorical encodings.  Categorical codes of \code{-1}
#' are mapped to \code{NA}.
#'
#' @param path Character(1). Path to the h5ad file.
#' @return A \code{data.frame} with rownames set to the cell index.
#' @keywords internal
#' @noRd
.h5ad_read_obs <- function(path) {
  ver    <- .h5ad_detect_version(path)
  ls_all <- rhdf5::h5ls(path, recursive = TRUE)
  rhdf5::H5close()

  obs_entries <- ls_all[ls_all$group == "/obs", ]

  if (!ver$modern) {
    ## -- Pre-0.8.0: __categories pattern ----------------------------------
    cat_names  <- ls_all[ls_all$group == "/obs/__categories", "name"]
    col_names  <- obs_entries$name[obs_entries$name != "__categories"]
    index_key  <- "index"  # always "index" in pre-0.8.0

    result <- list()
    for (col in col_names) {
      if (col == index_key) {
        result[[col]] <- rhdf5::h5read(path, paste0("/obs/", col))
      } else if (col %in% cat_names) {
        codes <- rhdf5::h5read(path, paste0("/obs/", col))
        cats  <- rhdf5::h5read(path, paste0("/obs/__categories/", col))
        vals  <- rep(NA_character_, length(codes))
        valid <- codes >= 0L
        vals[valid] <- cats[codes[valid] + 1L]
        result[[col]] <- factor(vals, levels = cats)
      } else {
        result[[col]] <- rhdf5::h5read(path, paste0("/obs/", col))
      }
      rhdf5::H5close()
    }

    df <- as.data.frame(result, stringsAsFactors = FALSE,
                         check.names = FALSE)
    if (index_key %in% names(df)) {
      rownames(df) <- df[[index_key]]
      df[[index_key]] <- NULL
    }

  } else {
    ## -- Modern (>=0.8.0): nested group categoricals ----------------------
    result    <- list()
    index_key <- NULL

    for (i in seq_len(nrow(obs_entries))) {
      col   <- obs_entries$name[i]
      otype <- obs_entries$otype[i]

      if (col == "_index" || col == "index") {
        index_key <- col
        result[[col]] <- rhdf5::h5read(path, paste0("/obs/", col))
      } else if (otype == "H5I_GROUP") {
        codes <- rhdf5::h5read(path, paste0("/obs/", col, "/codes"))
        cats  <- rhdf5::h5read(path, paste0("/obs/", col, "/categories"))
        vals  <- rep(NA_character_, length(codes))
        valid <- codes >= 0L
        vals[valid] <- cats[codes[valid] + 1L]
        result[[col]] <- factor(vals, levels = cats)
      } else {
        result[[col]] <- rhdf5::h5read(path, paste0("/obs/", col))
      }
      rhdf5::H5close()
    }

    df <- as.data.frame(result, stringsAsFactors = FALSE,
                         check.names = FALSE)
    if (!is.null(index_key) && index_key %in% names(df)) {
      rownames(df) <- df[[index_key]]
      df[[index_key]] <- NULL
    }
  }

  df
}

#' Read var (gene) names from an h5ad file via rhdf5
#'
#' @param path Character(1). Path to the h5ad file.
#' @return Character vector of gene/feature names.
#' @keywords internal
#' @noRd
.h5ad_read_var_names <- function(path) {
  ver <- .h5ad_detect_version(path)
  rhdf5::H5close()
  if (!ver$modern) {
    nms <- rhdf5::h5read(path, "/var/index")
  } else {
    nms <- tryCatch(
      rhdf5::h5read(path, "/var/_index"),
      error = function(e) rhdf5::h5read(path, "/var/index")
    )
  }
  rhdf5::H5close()
  as.character(nms)
}

#' Read obs (cell) names from an h5ad file via rhdf5
#'
#' @param path Character(1). Path to the h5ad file.
#' @return Character vector of cell IDs / barcodes.
#' @keywords internal
#' @noRd
.h5ad_read_obs_names <- function(path) {
  ver <- .h5ad_detect_version(path)
  rhdf5::H5close()
  if (!ver$modern) {
    nms <- rhdf5::h5read(path, "/obs/index")
  } else {
    nms <- tryCatch(
      rhdf5::h5read(path, "/obs/_index"),
      error = function(e) rhdf5::h5read(path, "/obs/index")
    )
  }
  rhdf5::H5close()
  as.character(nms)
}


# ──────────────────────────────────────────────────────────────────────
# Internal: resolve which HDF5 group holds the count matrix in an h5ad
# ──────────────────────────────────────────────────────────────────────

#' Determine the HDF5 group containing counts in an h5ad file
#'
#' Probes first for \code{/layers/counts}, then falls back to \code{/X}.
#' Uses \code{BPCells::open_matrix_anndata_hdf5} to test group validity
#' (no data is read — BPCells only checks group structure).
#'
#' @param file_path Character(1). Path to the h5ad file.
#' @param .h5ad.group Character(1) or \code{NULL}. If non-NULL, used
#'   directly without auto-detection.
#' @return Character(1) — the resolved HDF5 group path.
#' @keywords internal
#' @noRd
.h5ad_resolve_counts_group <- function(file_path, .h5ad.group = NULL) {
  if (!is.null(.h5ad.group)) {
    # User override — validate it actually works
    tryCatch(
      {
        BPCells::open_matrix_anndata_hdf5(file_path, group = .h5ad.group)
        return(.h5ad.group)
      },
      error = function(e) {
        stop("User-specified .h5ad.group '", .h5ad.group,
             "' could not be opened by BPCells in '", file_path, "'.\n",
             "Original error: ", conditionMessage(e), call. = FALSE)
      }
    )
  }

  # Auto-detect: try /layers/counts first, then /X
  for (group in c("/layers/counts", "/X")) {
    res <- tryCatch(
      {
        BPCells::open_matrix_anndata_hdf5(file_path, group = group)
        group
      },
      error = function(e) NULL
    )
    if (!is.null(res)) return(res)
  }

  stop("Could not find a valid count matrix in '", file_path,
       "'. Checked /layers/counts and /X. ",
       "Specify the correct group with .h5ad.group.", call. = FALSE)
}


#' @describeIn RunTDR Run the pipeline on an HDF5AnnData object
#'
#' Builds a \code{\linkS4class{TDRObj}} from an HDF5-backed AnnData
#' object. The expression matrix is converted to a BPCells on-disk
#' directory for efficient lazy access, then the proven
#' \code{IterableMatrix} pipeline is used for all downstream steps.
#'
#' Metadata (obs) is read via \code{anndataR}; the expression matrix
#' is read and converted via \code{BPCells}.
#'
#' @param x An \code{HDF5AnnData} object (created via
#'   \code{anndataR::read_h5ad(..., backend = "HDF5AnnData")}).
#' @param .sample.var Character(1). Column in \code{x$obs} identifying
#'   sample membership.
#' @param .assay.type Character. \code{"RNA"} or \code{"cyto"}.
#' @param .h5ad.group Character(1) or \code{NULL}. HDF5 group path
#'   containing the count matrix (e.g. \code{"/layers/counts"} or
#'   \code{"/X"}). If \code{NULL} (default), auto-detects by probing
#'   \code{/layers/counts} first, then falling back to \code{/X}.
#' @param .bpcells.dir Character(1) or \code{NULL}. Directory path for
#'   the BPCells on-disk matrix. If \code{NULL} (default), uses a
#'   temporary directory (\code{tempdir()}) that is cleaned up on session
#'   end. If a path is given and already contains a valid BPCells
#'   matrix, the conversion step is skipped (cache hit).
#' @param .harmony.var Character vector of batch variable column names
#'   in sample-level metadata, or \code{NULL}.
#' @param .markers Character vector of marker names (required for cyto).
#' @param .celltype.vec Character(1). Column name in \code{x$obs}
#'   containing per-cell type labels, or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for a sample to be
#'   included.
#' @param .verbose Logical. Print progress messages.
#' @param .seed Integer. Random seed.
#' @param .prop.landmarks Numeric in (0, 1]. Proportion of cells as landmarks.
#' @param .n.threads Integer. Number of threads.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return A \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.HDF5AnnData <- function(x,
                        .sample.var,
                        .assay.type = "RNA",
                        .h5ad.group = NULL,
                        .bpcells.dir = NULL,
                        .harmony.var = NULL,
                        .markers = NULL,
                        .celltype.vec = NULL,
                        .min.cells.per.sample = 10,
                        .verbose = TRUE,
                        .seed = 123,
                        .prop.landmarks = 0.1,
                        .n.threads = if (is.hpc()) {
                          max(RhpcBLASctl::blas_get_num_procs(),
                              RhpcBLASctl::omp_get_num_procs(),
                              RhpcBLASctl::omp_get_max_threads(),
                              na.rm = TRUE)
                        } else {
                          parallel::detectCores(logical = TRUE)
                        },
                        ...) {

  # --- Input validation ---
  if (!is.character(.sample.var) || length(.sample.var) != 1) {
    stop(".sample.var must be a single character string.")
  }

  if (!(.sample.var %in% colnames(x$obs))) {
    stop(".sample.var '", .sample.var, "' not found in x$obs.")
  }

  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for HDF5AnnData support. ",
         "Install it with: BiocManager::install('BPCells')", call. = FALSE)
  }

  # --- Get backing file path ---
  # HDF5AnnData stores the HDF5 handle in a private R6 field;
  # extract the filename via rhdf5 (which anndataR already depends on)
  h5ad_path <- tryCatch(
    {
      h5obj <- x$.__enclos_env__$private$.h5obj
      rhdf5::H5Fget_name(h5obj)
    },
    error = function(e) NULL
  )
  if (is.null(h5ad_path) || !file.exists(h5ad_path)) {
    stop("Could not determine h5ad file path from the HDF5AnnData object. ",
         "Ensure x was created with anndataR::read_h5ad().", call. = FALSE)
  }

  if (isTRUE(.verbose)) {
    cat("h5ad file: ", h5ad_path, "\n")
  }

  # --- Resolve HDF5 group ---
  group <- .h5ad_resolve_counts_group(h5ad_path, .h5ad.group)

  if (isTRUE(.verbose)) {
    cat("Using HDF5 group: ", group, "\n")
  }

  # --- Convert to BPCells on-disk format ---
  if (is.null(.bpcells.dir)) {
    .bpcells.dir <- file.path(tempdir(), paste0("bpcells_",
                              tools::file_path_sans_ext(basename(h5ad_path))))
  }

  # Check for cache hit: if the dir exists and contains a valid BPCells matrix
  cache_hit <- FALSE
  if (dir.exists(.bpcells.dir)) {
    bp_mat <- tryCatch(
      {
        BPCells::open_matrix_dir(.bpcells.dir)
      },
      error = function(e) NULL
    )
    if (!is.null(bp_mat)) {
      cache_hit <- TRUE
      if (isTRUE(.verbose)) {
        cat("BPCells cache hit: reusing ", .bpcells.dir, "\n")
      }
    }
  }

  if (!cache_hit) {
    if (isTRUE(.verbose)) {
      cat("Converting h5ad matrix to BPCells on-disk format...\n")
    }

    # Open h5ad matrix lazily via BPCells (nearly instant, no data read)
    h5_mat <- BPCells::open_matrix_anndata_hdf5(h5ad_path, group = group)

    # Streaming write to BPCells on-disk directory
    bp_mat <- BPCells::write_matrix_dir(mat = h5_mat, dir = .bpcells.dir)

    if (isTRUE(.verbose)) {
      cat("BPCells matrix written to: ", .bpcells.dir, "\n")
    }

    # Re-open from on-disk dir for consistent state
    bp_mat <- BPCells::open_matrix_dir(.bpcells.dir)
  }

  # --- Apply gene names from anndataR var metadata ---
  gene_names <- rownames(x$var)
  if (is.null(gene_names)) {
    gene_names <- x$var_names
  }
  if (!is.null(gene_names) && length(gene_names) == nrow(bp_mat)) {
    rownames(bp_mat) <- gene_names
  }

  # --- Apply cell names (barcodes) ---
  cell_names <- rownames(x$obs)
  if (is.null(cell_names)) {
    cell_names <- x$obs_names
  }
  if (!is.null(cell_names) && length(cell_names) == ncol(bp_mat)) {
    colnames(bp_mat) <- cell_names
  }

  # --- Build cell metadata from obs ---
  cell_meta <- as.data.frame(x$obs)

  # --- Delegate to .run_tdr_matrix (the proven IterableMatrix path) ---
  .run_tdr_matrix(
    x               = bp_mat,
    .cell.meta      = cell_meta,
    .sample.var     = .sample.var,
    .assay.type     = .assay.type,
    .harmony.var    = .harmony.var,
    .markers        = .markers,
    .celltype.vec   = .celltype.vec,
    .min.cells.per.sample = .min.cells.per.sample,
    .verbose        = .verbose,
    .seed           = .seed,
    .prop.landmarks = .prop.landmarks,
    .n.threads      = .n.threads,
    ...
  )
}


# ======================================================================
# RunTDR – character (h5ad file path) method
# ======================================================================

#' @describeIn RunTDR Run the pipeline directly from an h5ad file path
#'
#' Reads metadata and gene/cell names from the h5ad file using
#' \code{rhdf5}, opens the expression matrix with \code{BPCells}, and
#' delegates to the internal \code{IterableMatrix} pipeline.
#'
#' This method supports \strong{all} H5AD format versions, including
#' pre-0.8.0 files generated by older versions of Python anndata that
#' are not supported by the \code{anndataR} package.
#'
#' @param x Character(1).
#'   Path to a \code{.h5ad} file.
#' @param .sample.var Character(1). Column in \code{obs} identifying
#'   sample membership.
#' @param .assay.type Character. \code{"RNA"} or \code{"cyto"}.
#' @param .h5ad.group Character(1) or \code{NULL}. HDF5 group path
#'   containing the count matrix (e.g. \code{"/layers/counts"} or
#'   \code{"/X"}). If \code{NULL} (default), auto-detects by probing
#'   \code{/layers/counts} first, then falling back to \code{/X}.
#' @param .bpcells.dir Character(1) or \code{NULL}. Directory path for
#'   the BPCells on-disk matrix. If \code{NULL} (default), uses a
#'   temporary directory (\code{tempdir()}) that is cleaned up on session
#'   end. If a path is given and already contains a valid BPCells
#'   matrix, the conversion step is skipped (cache hit).
#' @param .harmony.var Character vector of batch variable column names
#'   in sample-level metadata, or \code{NULL}.
#' @param .markers Character vector of marker names (required for cyto).
#' @param .celltype.vec Character(1). Column name in \code{obs}
#'   containing per-cell type labels, or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for a sample to be
#'   included.
#' @param .verbose Logical. Print progress messages.
#' @param .seed Integer. Random seed.
#' @param .prop.landmarks Numeric in (0, 1]. Proportion of cells as landmarks.
#' @param .n.threads Integer. Number of threads.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return A \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.character <- function(x,
                             .sample.var,
                             .assay.type = "RNA",
                             .h5ad.group = NULL,
                             .bpcells.dir = NULL,
                             .harmony.var = NULL,
                             .markers = NULL,
                             .celltype.vec = NULL,
                             .min.cells.per.sample = 10,
                             .verbose = TRUE,
                             .seed = 123,
                             .prop.landmarks = 0.1,
                             .n.threads = if (is.hpc()) {
                               max(RhpcBLASctl::blas_get_num_procs(),
                                   RhpcBLASctl::omp_get_num_procs(),
                                   RhpcBLASctl::omp_get_max_threads(),
                                   na.rm = TRUE)
                             } else {
                               parallel::detectCores(logical = TRUE)
                             },
                             ...) {

  # --- Input validation ---
  if (length(x) != 1L) {
    stop("x must be a single file path (length-1 character string).",
         call. = FALSE)
  }
  if (!file.exists(x)) {
    stop("File not found: '", x, "'.", call. = FALSE)
  }
  if (!grepl("\\.h5ad$", x, ignore.case = TRUE)) {
    stop("x does not have an .h5ad extension: '", x, "'.", call. = FALSE)
  }
  if (!is.character(.sample.var) || length(.sample.var) != 1) {
    stop(".sample.var must be a single character string.", call. = FALSE)
  }

  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required for h5ad file path support. ",
         "Install it with: BiocManager::install('rhdf5')", call. = FALSE)
  }
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for h5ad file path support. ",
         "Install it with: BiocManager::install('BPCells')", call. = FALSE)
  }

  h5ad_path <- normalizePath(x, mustWork = TRUE)

  if (isTRUE(.verbose)) {
    cat("h5ad file: ", h5ad_path, "\n")
    ver <- .h5ad_detect_version(h5ad_path)
    cat("H5AD format: ",
        if (ver$modern) paste0("modern (", ver$encoding_version, ")")
        else "pre-0.8.0 (legacy)",
        "\n")
  }

  # --- Read cell metadata from obs ---
  if (isTRUE(.verbose)) cat("Reading obs metadata via rhdf5...\n")
  cell_meta <- .h5ad_read_obs(h5ad_path)

  if (!(.sample.var %in% colnames(cell_meta))) {
    stop(".sample.var '", .sample.var, "' not found in h5ad obs. ",
         "Available columns: ",
         paste(head(colnames(cell_meta), 20), collapse = ", "),
         call. = FALSE)
  }

  # --- Resolve HDF5 group ---
  group <- .h5ad_resolve_counts_group(h5ad_path, .h5ad.group)

  if (isTRUE(.verbose)) {
    cat("Using HDF5 group: ", group, "\n")
  }

  # --- Convert to BPCells on-disk format ---
  if (is.null(.bpcells.dir)) {
    .bpcells.dir <- file.path(tempdir(), paste0("bpcells_",
                              tools::file_path_sans_ext(basename(h5ad_path))))
  }

  cache_hit <- FALSE
  if (dir.exists(.bpcells.dir)) {
    bp_mat <- tryCatch(BPCells::open_matrix_dir(.bpcells.dir),
                       error = function(e) NULL)
    if (!is.null(bp_mat)) {
      cache_hit <- TRUE
      if (isTRUE(.verbose)) {
        cat("BPCells cache hit: reusing ", .bpcells.dir, "\n")
      }
    }
  }

  if (!cache_hit) {
    if (isTRUE(.verbose)) {
      cat("Converting h5ad matrix to BPCells on-disk format...\n")
    }
    h5_mat <- BPCells::open_matrix_anndata_hdf5(h5ad_path, group = group)
    bp_mat <- BPCells::write_matrix_dir(mat = h5_mat, dir = .bpcells.dir)
    if (isTRUE(.verbose)) {
      cat("BPCells matrix written to: ", .bpcells.dir, "\n")
    }
    bp_mat <- BPCells::open_matrix_dir(.bpcells.dir)
  }

  # --- Apply gene names ---
  gene_names <- .h5ad_read_var_names(h5ad_path)
  if (length(gene_names) == nrow(bp_mat)) {
    rownames(bp_mat) <- gene_names
  }

  # --- Apply cell names ---
  cell_names <- .h5ad_read_obs_names(h5ad_path)
  if (length(cell_names) == ncol(bp_mat)) {
    colnames(bp_mat) <- cell_names
  }

  # --- Delegate to .run_tdr_matrix ---
  .run_tdr_matrix(
    x               = bp_mat,
    .cell.meta      = cell_meta,
    .sample.var     = .sample.var,
    .assay.type     = .assay.type,
    .harmony.var    = .harmony.var,
    .markers        = .markers,
    .celltype.vec   = .celltype.vec,
    .min.cells.per.sample = .min.cells.per.sample,
    .verbose        = .verbose,
    .seed           = .seed,
    .prop.landmarks = .prop.landmarks,
    .n.threads      = .n.threads,
    ...
  )
}


# ======================================================================
# RunTDR – flow cytometry methods (cytoset / flowSet)
# ======================================================================

#' @describeIn RunTDR Run the pipeline on a flowWorkspace cytoset
#'
#' Builds a \code{\linkS4class{TDRObj}} from a \code{cytoset} object
#' (one FCS sample per \code{cytoframe}) and executes the full pipeline.
#'
#' @param x A \code{flowWorkspace::cytoset} object.
#' @param .sample.var Character(1). Column name in \code{pData(x)}
#'   identifying sample membership.
#' @param .assay.type Character. Must be \code{"cyto"} (only valid type for
#'   flow cytometry data).
#' @param .harmony.var Character vector of batch variable column names
#'   in pData, or \code{NULL}.
#' @param .markers Character vector of marker/channel names. If \code{NULL},
#'   defaults to all channels in the cytoset.
#' @param .celltype.vec Named character vector of per-cell type labels
#'   (names = cell IDs, values = labels), or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for inclusion.
#' @param .verbose Logical.
#' @param .seed Integer.
#' @param .prop.landmarks Numeric in (0, 1].
#' @param .n.threads Integer.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return The updated \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.cytoset <- function(x,
                           .sample.var,
                           .assay.type = "cyto",
                           .harmony.var = NULL,
                           .markers = NULL,
                           .celltype.vec = NULL,
                           .min.cells.per.sample = 10,
                           .verbose = TRUE,
                           .seed = 123,
                           .prop.landmarks = 0.1,
                           .n.threads = if (is.hpc()) {
                             max(RhpcBLASctl::blas_get_num_procs(),
                                 RhpcBLASctl::omp_get_num_procs(),
                                 RhpcBLASctl::omp_get_max_threads(),
                                 na.rm = TRUE)
                           } else {
                             parallel::detectCores(logical = TRUE)
                           },
                           ...) {
  .RunTDR_flow_common(
    x = x, .sample.var = .sample.var, .assay.type = .assay.type,
    .harmony.var = .harmony.var, .markers = .markers,
    .celltype.vec = .celltype.vec,
    .min.cells.per.sample = .min.cells.per.sample,
    .verbose = .verbose, .seed = .seed,
    .prop.landmarks = .prop.landmarks, .n.threads = .n.threads, ...
  )
}

#' @describeIn RunTDR Run the pipeline on a flowCore flowSet
#'
#' Builds a \code{\linkS4class{TDRObj}} from a \code{flowSet} object
#' (one FCS sample per \code{flowFrame}) and executes the full pipeline.
#' Requires only \pkg{flowCore} (not \pkg{flowWorkspace}).
#'
#' @param x A \code{flowCore::flowSet} object.
#' @param .sample.var Character(1). Column name in \code{pData(x)}
#'   identifying sample membership.
#' @param .assay.type Character. Must be \code{"cyto"} (only valid type for
#'   flow cytometry data).
#' @param .harmony.var Character vector of batch variable column names
#'   in pData, or \code{NULL}.
#' @param .markers Character vector of marker/channel names. If \code{NULL},
#'   defaults to all channels in the flowSet.
#' @param .celltype.vec Named character vector of per-cell type labels
#'   (names = cell IDs, values = labels), or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for inclusion.
#' @param .verbose Logical.
#' @param .seed Integer.
#' @param .prop.landmarks Numeric in (0, 1].
#' @param .n.threads Integer.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return The updated \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.flowSet <- function(x,
                           .sample.var,
                           .assay.type = "cyto",
                           .harmony.var = NULL,
                           .markers = NULL,
                           .celltype.vec = NULL,
                           .min.cells.per.sample = 10,
                           .verbose = TRUE,
                           .seed = 123,
                           .prop.landmarks = 0.1,
                           .n.threads = if (is.hpc()) {
                             max(RhpcBLASctl::blas_get_num_procs(),
                                 RhpcBLASctl::omp_get_num_procs(),
                                 RhpcBLASctl::omp_get_max_threads(),
                                 na.rm = TRUE)
                           } else {
                             parallel::detectCores(logical = TRUE)
                           },
                           ...) {
  .RunTDR_flow_common(
    x = x, .sample.var = .sample.var, .assay.type = .assay.type,
    .harmony.var = .harmony.var, .markers = .markers,
    .celltype.vec = .celltype.vec,
    .min.cells.per.sample = .min.cells.per.sample,
    .verbose = .verbose, .seed = .seed,
    .prop.landmarks = .prop.landmarks, .n.threads = .n.threads, ...
  )
}


# ──────────────────────────────────────────────────────────────────────
# Internal helper: shared logic for RunTDR.cytoset / RunTDR.flowSet
# ──────────────────────────────────────────────────────────────────────

#' Shared pipeline for cytoset and flowSet inputs
#'
#' @keywords internal
#' @noRd
.RunTDR_flow_common <- function(x,
                                .sample.var,
                                .assay.type = "cyto",
                                .harmony.var = NULL,
                                .markers = NULL,
                                .celltype.vec = NULL,
                                .min.cells.per.sample = 10,
                                .verbose = TRUE,
                                .seed = 123,
                                .prop.landmarks = 0.1,
                                .n.threads = if (is.hpc()) {
                                  max(RhpcBLASctl::blas_get_num_procs(),
                                      RhpcBLASctl::omp_get_num_procs(),
                                      RhpcBLASctl::omp_get_max_threads(),
                                      na.rm = TRUE)
                                } else {
                                  parallel::detectCores(logical = TRUE)
                                },
                                ...) {

  # --- Guard: flowCore always required ---
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required for flow cytometry input. ",
         "Install it with: BiocManager::install('flowCore')", call. = FALSE)
  }

  # --- Guard: flowWorkspace required only for cytoset ---
  is_cytoset <- inherits(x, "cytoset")
  if (is_cytoset && !requireNamespace("flowWorkspace", quietly = TRUE)) {
    stop("Package 'flowWorkspace' is required for cytoset input. ",
         "Install it with: BiocManager::install('flowWorkspace')", call. = FALSE)
  }

  # --- Validate input ---
  if (!inherits(x, "flowSet")) {
    stop("x must be a flowCore::flowSet or flowWorkspace::cytoset object.")
  }

  if (!identical(.assay.type, "cyto")) {
    stop(".assay.type must be 'cyto' for flow cytometry input. ",
         "RNA data is not supported with flowSet/cytoset.")
  }

  if (!is.character(.sample.var) || length(.sample.var) != 1) {
    stop(".sample.var must be a single character string.")
  }

  # pData and sampleNames are Biobase generics re-exported by flowCore.
  # flowCore registers the flowSet method; flowWorkspace registers the

  # cytoset method (loaded above via requireNamespace when needed).
  if (!(.sample.var %in% colnames(flowCore::pData(x)))) {
    stop(".sample.var '", .sample.var,
         "' not found in pData(x).")
  }

  # --- Extract sample metadata ---
  sample_meta <- as.data.frame(flowCore::pData(x))
  rownames(sample_meta) <- sample_meta[[.sample.var]]
  sample_meta <- sample_meta[, !colnames(sample_meta) %in% .sample.var,
                             drop = FALSE]

  # --- Build .cells ---
  snames <- flowCore::sampleNames(x)
  .cells <- lapply(
    stats::setNames(snames, snames),
    function(s) seq_len(nrow(x[[s]]))
  )

  # --- Filter by .min.cells.per.sample ---
  keep <- vapply(.cells, length, integer(1)) >= .min.cells.per.sample
  .cells <- .cells[keep]
  sample_meta <- sample_meta[names(.cells), , drop = FALSE]

  if (length(.cells) == 0) {
    stop("No samples have >= ", .min.cells.per.sample,
         " cells. Check .min.cells.per.sample.")
  }

  # --- Resolve markers ---
  all_channels <-
    flowCore::parameters(object = x[[flowCore::sampleNames(x)[1]]])$name |>
    unname()
  if (is.null(.markers)) {
    .markers <- all_channels
  } else {
    bad <- setdiff(.markers, all_channels)
    if (length(bad) > 0) {
      stop("Markers not found in channels: ",
           paste(bad, collapse = ", "))
    }
  }

  # --- Build source_env ---
  source_env <- new.env(parent = emptyenv())
  source_env$cs <- x
  lockBinding(sym = as.name("cs"), env = source_env)

  # --- .celltype.vec handling ---
  ct_vec <- NULL
  if (!is.null(.celltype.vec)) {
    if (!is.character(.celltype.vec) || is.null(names(.celltype.vec))) {
      stop(".celltype.vec must be a named character vector ",
           "(names = cell IDs, values = cell type labels).")
    }
    ct_vec <- .celltype.vec
  }

  # --- Build TDRObj ---
  tdr.obj <- .setup_tdr_from_cytoset(
    .cells = .cells,
    .source.env = source_env,
    .meta = sample_meta,
    .markers = .markers,
    .harmony.var = .harmony.var,
    .celltype.vec = ct_vec,
    .prop.landmarks = .prop.landmarks,
    .seed = .seed,
    .n.threads = .n.threads,
    .verbose = .verbose
  )

  # --- Argument routing ---
  dots <- list(...)

  .lm.formals  <- setdiff(names(formals(get.landmarks.TDRObj)), c("x", "..."))
  .gr.formals  <- setdiff(names(formals(get.graph.TDRObj)), c("x", "..."))
  .map.formals <- setdiff(names(formals(get.map.TDRObj)), c("x", "..."))
  .emb.formals <- setdiff(names(formals(get.embedding.TDRObj)), c("x", "..."))

  lm_args  <- dots[names(dots) %in% .lm.formals]
  gr_args  <- dots[names(dots) %in% .gr.formals]
  map_args <- dots[names(dots) %in% .map.formals]
  emb_args <- dots[names(dots) %in% .emb.formals]

  all_known <- Reduce(union, list(.lm.formals, .gr.formals, .map.formals, .emb.formals))
  orphans <- setdiff(names(dots), all_known)
  if (length(orphans) > 0) {
    warning("Unknown arguments will be ignored: ",
            paste(orphans, collapse = ", "))
  }

  # --- Pipeline ---
  tdr.obj <- do.call(get.landmarks.TDRObj, c(list(tdr.obj, .source = NULL, .seed = .seed, .verbose = .verbose), lm_args))
  tdr.obj <- do.call(get.graph.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), gr_args))

  if (!is.null(tdr.obj@config$celltype.vec)) {
    tdr.obj <- celltyping(tdr.obj,
                          .celltyping.map = tdr.obj@config$celltype.vec,
                          .verbose = .verbose)
  }

  tdr.obj <- do.call(get.map.TDRObj, c(list(tdr.obj, .source = NULL, .seed = .seed, .verbose = .verbose), map_args))
  tdr.obj <- do.call(get.embedding.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), emb_args))

  return(tdr.obj)
}


# ──────────────────────────────────────────────────────────────────────
# Internal helper: build TDRObj from flow-cytometry-derived inputs
# ──────────────────────────────────────────────────────────────────────

#' Build a TDRObj from cytoset/flowSet-derived cell lists and metadata
#'
#' @keywords internal
#' @noRd
.setup_tdr_from_cytoset <- function(.cells,
                                    .source.env,
                                    .meta,
                                    .markers,
                                    .harmony.var,
                                    .celltype.vec,
                                    .prop.landmarks,
                                    .seed,
                                    .n.threads,
                                    .verbose) {

  # --- Harmony var validation ---
  if (!is.null(.harmony.var)) {
    if (!inherits(.harmony.var, "character")) {
      stop(".harmony.var must be a character vector.")
    }
    if (!all(.harmony.var %in% colnames(.meta))) {
      stop("Variables not found in metadata: ",
           paste(.harmony.var[!(.harmony.var %in% colnames(.meta))],
                 collapse = ", "),
           "\nCheck column names in .meta with colnames(.meta).")
    }
  }

  # --- Markers validation ---
  if (is.null(.markers) || length(.markers) < 3) {
    stop(".markers must contain at least 3 markers for meaningful dimensionality reduction.")
  }

  # --- Prop landmarks validation ---
  if ((.prop.landmarks < 0) | (.prop.landmarks > 1)) {
    stop(".prop.landmarks must be between 0 and 1 (e.g., 0.1 for 10% of cells).\n",
         "Current value: ", .prop.landmarks)
  }

  # --- Create TDRObj ---
  .tdr.obj <- TDRObj(
    config = list(
      key = NULL,
      sampling = NULL,
      assay.type = "cyto",
      markers = NULL,
      n.threads = .n.threads
    ),
    integration = list(
      harmony.var = NULL,
      harmony.obj = NULL
    )
  )

  .tdr.obj@cells <- .cells

  # --- n.cells from integer index vector lengths ---
  n.cells <- lengths(.cells)

  .tdr.obj@config$sampling$n.cells <- n.cells

  # Quality check: warn if sample sizes are highly imbalanced (>10-fold difference)
  if ((max(.tdr.obj@config$sampling$n.cells) /
       min(.tdr.obj@config$sampling$n.cells)) > 10) {

    warning("Sample size imbalance detected: largest/smallest ratio > 10.\n",
            "Smallest sample has ", min(.tdr.obj@config$sampling$n.cells),
            " cells.\n",
            "Consider removing low-quality samples.")

    if (any(.tdr.obj@config$sampling$n.cells < 1000)) {
      warning("Large variation in sample sizes detected. ",
              "For cytometry, samples with <1000 cells may be unreliable.")
    }
  }

  # Calculate target number of landmarks: prop of total cells, capped at 5000
  .tdr.obj@config$sampling$target.lm.n <-
    pmin(sum(.tdr.obj@config$sampling$n.cells) * .prop.landmarks,
         5e3)

  # Allocate landmarks per sample: proportional to sample size, but capped
  .tdr.obj@config$sampling$n.perSample <-
    pmin(ceiling(x = .tdr.obj@config$sampling$n.cells * .prop.landmarks),
         ceiling(x = .tdr.obj@config$sampling$target.lm.n /
                     length(x = .tdr.obj@cells)))

  # Create key vector: maps each future landmark to its sample
  .tdr.obj@config$key <-
    seq_along(along.with = .tdr.obj@cells) |>
    rep(times = .tdr.obj@config$sampling$n.perSample) |>
    (\(x)
     stats::setNames(object = x,
                     nm = names(.tdr.obj@cells)[x])
    )()

  .tdr.obj@metadata <- .meta

  .tdr.obj@metadata$n.perSample <-
    .tdr.obj@config$sampling$n.perSample

  .tdr.obj@metadata$n.cells <-
    .tdr.obj@config$sampling$n.cells

  .tdr.obj@metadata$log10.n.cells <-
    log10(x = .tdr.obj@config$sampling$n.cells)

  # --- Markers ---
  .tdr.obj@config$markers <- .markers

  # --- Harmony ---
  if (!is.null(.harmony.var)) {
    .tdr.obj@integration$harmony.var <- .harmony.var
  }

  # --- Cytometry backend: store locked env for .get_sample_matrix ---
  .tdr.obj@config$backend <- "cyto"
  .tdr.obj@config$source.env <- .source.env

  # --- Cell type vector ---
  if (!is.null(.celltype.vec)) {
    .tdr.obj@config$celltype.vec <- .celltype.vec
  }

  return(.tdr.obj)
}


# ======================================================================
# RunTDR – matrix methods (dgCMatrix, DelayedMatrix, IterableMatrix)
# ======================================================================

#' @describeIn RunTDR Run the pipeline on a sparse matrix (dgCMatrix)
#'
#' Builds a \code{\linkS4class{TDRObj}} from a \code{dgCMatrix} and
#' per-cell metadata, then executes the full pipeline.
#'
#' @param x A \code{dgCMatrix} (features x cells for RNA,
#'   cells x features for cyto).
#' @param .cell.meta A \code{data.frame} of per-cell metadata.
#'   Must have one row per cell with rownames matching cell IDs
#'   in \code{x} (\code{colnames} for RNA, \code{rownames} for cyto).
#' @param .sample.var Character(1). Column name in \code{.cell.meta}
#'   identifying sample membership.
#' @param .assay.type Character. \code{"RNA"} or \code{"cyto"}.
#' @param .harmony.var Character vector of batch variable column names
#'   in sample-level metadata, or \code{NULL}.
#' @param .markers Character vector of marker names (required for cyto).
#' @param .celltype.vec Character(1). Column name in \code{.cell.meta}
#'   containing per-cell type labels, or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for a sample to be
#'   included.
#' @param .verbose Logical. Print progress messages.
#' @param .seed Integer. Random seed.
#' @param .prop.landmarks Numeric in (0, 1]. Proportion of cells as landmarks.
#' @param .n.threads Integer. Number of threads.
#' @param ... Additional arguments passed to pipeline functions.
#'
#' @return A \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.dgCMatrix <- function(x, .cell.meta, ...) {
  dots <- list(...)
  if (!is.null(dots$.assay.type) && dots$.assay.type != "RNA") {
    stop("dgCMatrix input is only supported for .assay.type = 'RNA'.\n",
         "For cytometry data, supply a dense matrix instead.")
  }
  .run_tdr_matrix(x, .cell.meta, ...)
}


# ──────────────────────────────────────────────────────────────────────
# Internal: convert a DelayedMatrix to BPCells on-disk format
# ──────────────────────────────────────────────────────────────────────

#' Convert a DelayedMatrix to a BPCells on-disk IterableMatrix
#'
#' Converts a \code{DelayedMatrix} to a BPCells \code{IterableMatrix}
#' without ever fully materializing the matrix in R memory.
#' Three strategies are tried in order:
#' (1) if the seed is already a BPCells \code{IterableMatrix}, return it
#'     directly; (2) if the seed is an \code{HDF5ArraySeed}, open it via
#'     \code{BPCells::open_matrix_hdf5()} for zero-copy streaming;
#' (3) otherwise, iterate over column chunks of size \code{.chunk.size},
#'     coercing each to \code{dgCMatrix} then \code{IterableMatrix},
#'     lazy-\code{cbind}, and write to disk. Peak memory = one chunk.
#' If \code{.bpcells.dir} already contains a valid BPCells matrix,
#' the conversion is skipped (cache hit).
#'
#' @param mat A \code{DelayedMatrix}.
#' @param .bpcells.dir Character(1) or \code{NULL}. Directory for the
#'   BPCells on-disk matrix. If \code{NULL}, uses a temporary directory.
#' @param .verbose Logical. Print progress messages.
#' @param .chunk.size Integer. Number of columns per chunk for the
#'   chunked fallback path. Default 5000.
#' @return A BPCells \code{IterableMatrix} opened from the on-disk
#'   directory, or the seed directly if already BPCells-backed.
#' @keywords internal
#' @noRd
.delayed_to_bpcells <- function(mat, .bpcells.dir = NULL, .verbose = TRUE,
                                .chunk.size = 5000L) {

  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for DelayedMatrix \u2192 BPCells conversion. ",
         "Install it with: BiocManager::install('BPCells')", call. = FALSE)
  }

  # --- Check if the DelayedMatrix is already BPCells-backed ---
  seed <- DelayedArray::seed(mat)
  if (methods::is(seed, "IterableMatrix")) {
    if (isTRUE(.verbose)) {
      cat("DelayedMatrix is already BPCells-backed; skipping conversion.\n")
    }
    return(seed)
  }

  # --- Default directory ---
  if (is.null(.bpcells.dir)) {
    .bpcells.dir <- file.path(tempdir(), paste0("bpcells_delayed_",
                              nrow(mat), "x", ncol(mat)))
  }

  # --- Cache check ---
  if (dir.exists(.bpcells.dir)) {
    bp_mat <- tryCatch(
      BPCells::open_matrix_dir(.bpcells.dir),
      error = function(e) NULL
    )
    if (!is.null(bp_mat)) {
      if (isTRUE(.verbose)) {
        cat("BPCells cache hit: reusing ", .bpcells.dir, "\n")
      }
      return(bp_mat)
    }
  }

  # --- HDF5-backed: try zero-copy streaming via BPCells, fall back to chunked ---
  if (requireNamespace("HDF5Array", quietly = TRUE) &&
      methods::is(seed, "HDF5ArraySeed")) {
    h5_path <- HDF5Array::path(seed)
    h5_group <- seed@name
    if (isTRUE(.verbose)) {
      cat("HDF5-backed DelayedMatrix detected; attempting BPCells streaming...\n")
    }
    bp_mat <- tryCatch({
      bp <- BPCells::open_matrix_hdf5(path = h5_path, group = h5_group)
      BPCells::write_matrix_dir(mat = bp, dir = .bpcells.dir)
      BPCells::open_matrix_dir(.bpcells.dir)
    }, error = function(e) NULL)

    if (!is.null(bp_mat)) {
      if (isTRUE(.verbose)) {
        cat("BPCells matrix written to: ", .bpcells.dir, "\n")
      }
      return(bp_mat)
    }
    if (isTRUE(.verbose)) {
      cat("BPCells cannot read this HDF5 layout; falling back to chunked conversion...\n")
    }
  }

  # --- Chunked fallback for non-HDF5 and incompatible-HDF5 backends ---
  if (isTRUE(.verbose)) {
    cat("Converting DelayedMatrix to BPCells on-disk format (chunked)...\n")
  }
  n <- ncol(mat)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / .chunk.size))
  bp_parts <- lapply(chunks, function(idx) {
    methods::as(methods::as(mat[, idx], "dgCMatrix"), "IterableMatrix")
  })
  bp_full <- do.call(cbind, bp_parts)
  BPCells::write_matrix_dir(mat = bp_full, dir = .bpcells.dir)

  if (isTRUE(.verbose)) {
    cat("BPCells matrix written to: ", .bpcells.dir, "\n")
  }

  # Re-open from on-disk dir for consistent state
  BPCells::open_matrix_dir(.bpcells.dir)
}


#' @describeIn RunTDR Run the pipeline on a DelayedMatrix
#'
#' Converts the \code{DelayedMatrix} to a BPCells on-disk
#' \code{IterableMatrix} for efficient lazy access, then delegates
#' to the proven \code{IterableMatrix} pipeline via
#' \code{.run_tdr_matrix()}.
#'
#' @param x A \code{DelayedMatrix} (features \ifelse{html}{\out{&times;}}{\eqn{\times}} cells for RNA).
#' @param .cell.meta A \code{data.frame} of per-cell metadata.
#'   Must have one row per cell with rownames matching cell IDs
#'   in \code{x} (\code{colnames} for RNA, \code{rownames} for cyto).
#' @param .bpcells.dir Character(1) or \code{NULL}. Directory path for
#'   the BPCells on-disk matrix. If \code{NULL} (default), uses a
#'   temporary directory (\code{tempdir()}) that is cleaned up on session
#'   end. If a path is given and already contains a valid BPCells
#'   matrix, the conversion step is skipped (cache hit).
#' @param ... Additional arguments passed to \code{.run_tdr_matrix}
#'   (e.g. \code{.sample.var}, \code{.assay.type}, \code{.harmony.var},
#'   \code{.markers}, \code{.celltype.vec}, \code{.min.cells.per.sample},
#'   \code{.verbose}, \code{.seed}, \code{.prop.landmarks},
#'   \code{.n.threads}).
#'
#' @return A \code{\linkS4class{TDRObj}}.
#'
#' @export
RunTDR.DelayedMatrix <- function(x, .cell.meta, .bpcells.dir = NULL, ...) {
  dots <- list(...)
  if (!is.null(dots$.assay.type) && dots$.assay.type != "RNA") {
    stop("DelayedMatrix input is only supported for .assay.type = 'RNA'.\n",
         "For cytometry data, supply a dense matrix instead.")
  }

  .verbose <- if (!is.null(dots$.verbose)) dots$.verbose else TRUE

  # --- Convert to BPCells on-disk format ---
  bp_mat <- .delayed_to_bpcells(x, .bpcells.dir, .verbose)

  # --- Preserve dimnames if lost during conversion ---
  if (!is.null(rownames(x)) && is.null(rownames(bp_mat))) {
    rownames(bp_mat) <- rownames(x)
  }
  if (!is.null(colnames(x)) && is.null(colnames(bp_mat))) {
    colnames(bp_mat) <- colnames(x)
  }

  # --- Delegate to .run_tdr_matrix (the proven IterableMatrix path) ---
  .run_tdr_matrix(bp_mat, .cell.meta, ...)
}


#' @describeIn RunTDR Run the pipeline on a BPCells IterableMatrix
#'
#' @export
RunTDR.IterableMatrix <- function(x, .cell.meta, ...) {
  dots <- list(...)
  if (!is.null(dots$.assay.type) && dots$.assay.type != "RNA") {
    stop("IterableMatrix input is only supported for .assay.type = 'RNA'.\n",
         "For cytometry data, supply a dense matrix instead.")
  }
  .run_tdr_matrix(x, .cell.meta, ...)
}


# ──────────────────────────────────────────────────────────────────────
# Internal: shared implementation for all matrix backends
# ──────────────────────────────────────────────────────────────────────

#' Run TDR pipeline on a bare matrix + cell metadata
#'
#' @keywords internal
#' @noRd
.run_tdr_matrix <- function(x,
                            .cell.meta,
                            .sample.var,
                            .assay.type = "RNA",
                            .harmony.var = NULL,
                            .markers = NULL,
                            .celltype.vec = NULL,
                            .min.cells.per.sample = 10,
                            .verbose = TRUE,
                            .seed = 123,
                            .prop.landmarks = 0.1,
                            .n.threads = if (is.hpc()) {
                              max(RhpcBLASctl::blas_get_num_procs(),
                                  RhpcBLASctl::omp_get_num_procs(),
                                  RhpcBLASctl::omp_get_max_threads(),
                                  na.rm = TRUE)
                            } else {
                              parallel::detectCores(logical = TRUE)
                            },
                            ...) {

  # --- Input validation ---
  .assay.type <- match.arg(arg = .assay.type,
                           choices = c("cyto", "RNA"))

  if (!inherits(.cell.meta, "data.frame")) {
    stop(".cell.meta must be a data.frame.")
  }

  if (!is.character(.sample.var) || length(.sample.var) != 1) {
    stop(".sample.var must be a single character string.")
  }

  if (!(.sample.var %in% colnames(.cell.meta))) {
    stop(".sample.var '", .sample.var, "' not found in .cell.meta.")
  }

  # Determine cell IDs based on assay type orientation
  if (.assay.type == "RNA") {
    cell_ids <- colnames(x)
    if (is.null(cell_ids)) {
      stop("Matrix must have colnames (cell IDs) for RNA assay type.")
    }
  } else {
    cell_ids <- rownames(x)
    if (is.null(cell_ids)) {
      stop("Matrix must have rownames (cell IDs) for cyto assay type.")
    }
  }

  # Validate .cell.meta rownames match matrix cell IDs
  if (is.null(rownames(.cell.meta))) {
    stop(".cell.meta must have rownames matching cell IDs in the matrix.")
  }

  shared <- intersect(rownames(.cell.meta), cell_ids)
  if (length(shared) == 0) {
    stop("No overlap between rownames(.cell.meta) and cell IDs in the matrix.")
  }
  if (length(shared) < length(cell_ids)) {
    warning(length(cell_ids) - length(shared),
            " cells in matrix not found in .cell.meta; they will be dropped.")
  }

  # Subset to shared cells
  if (.assay.type == "RNA") {
    x <- x[, shared, drop = FALSE]
  } else {
    x <- x[shared, , drop = FALSE]
  }
  .cell.meta <- .cell.meta[shared, , drop = FALSE]

  # --- Derive sample-level metadata ---
  sample_ids <- .cell.meta[[.sample.var]]
  counts_tbl <- table(sample_ids)
  valid <- names(counts_tbl)[counts_tbl >= .min.cells.per.sample]

  if (length(valid) == 0) {
    stop("No samples have >= ", .min.cells.per.sample,
         " cells. Check .sample.var or lower .min.cells.per.sample.")
  }

  # Identify sample-level columns (constant within each sample)
  sample_level_cols <- .cell.meta |>
    dplyr::group_by(!!rlang::sym(.sample.var)) |>
    dplyr::summarize(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = ~ length(unique(.x)) == 1
      ),
      .groups = "drop"
    ) |>
    (\(df) {
      non_sv <- df[, !colnames(df) %in% .sample.var, drop = FALSE]
      colnames(dplyr::select(non_sv, dplyr::where(all)))
    })()

  sample_level_cols <- c(.sample.var, sample_level_cols)

  sample_meta <- .cell.meta[, sample_level_cols, drop = FALSE] |>
    dplyr::distinct() |>
    as.data.frame()
  rownames(sample_meta) <- sample_meta[[.sample.var]]

  # Align sample ordering with first-occurrence order (from dplyr::distinct)
  # for parity with the explicit get.meta() → get.cells() → setup.tdr.obj()
  # workflow.  intersect() preserves the order of its first argument, so put
  # the first-occurrence rownames first.
  sample_meta <- sample_meta[intersect(rownames(sample_meta), valid),
                             , drop = FALSE]

  # --- Build index-based .cells + locked source env ---
  cell_sample <- .cell.meta[[.sample.var]]
  valid_mask <- cell_sample %in% rownames(sample_meta)

  if (.assay.type == "RNA") {
    x <- x[, valid_mask, drop = FALSE]
  } else {
    x <- x[valid_mask, , drop = FALSE]
  }
  .cell.meta <- .cell.meta[valid_mask, , drop = FALSE]
  cell_sample <- cell_sample[valid_mask]

  # Store source matrix in a locked environment to prevent duplication
  source_env <- new.env(parent = emptyenv())
  source_env$mat <- x
  lockBinding(sym = as.name("mat"), env = source_env)

  # Build .cells as named list of integer index vectors (like Seurat/SCE)
  sample_names <- rownames(sample_meta)
  if (.assay.type == "RNA") {
    .cells <- lapply(
      stats::setNames(sample_names, sample_names),
      function(s) sort(which(cell_sample == s))
    )
  } else {
    .cells <- lapply(
      stats::setNames(sample_names, sample_names),
      function(s) sort(which(cell_sample == s))
    )
  }

  # --- .celltype.vec handling ---
  ct_vec <- NULL
  if (!is.null(.celltype.vec)) {
    if (!is.character(.celltype.vec) || length(.celltype.vec) != 1) {
      stop(".celltype.vec must be a single character string ",
           "(column name in .cell.meta).")
    }
    if (!(.celltype.vec %in% colnames(.cell.meta))) {
      stop(".celltype.vec '", .celltype.vec,
           "' not found in .cell.meta.")
    }
    valid_cells <- if (.assay.type == "RNA") {
      colnames(x)
    } else {
      rownames(x)
    }
    ct_vec <- stats::setNames(
      as.character(.cell.meta[valid_cells, .celltype.vec]),
      valid_cells
    )
  }

  # --- Build TDRObj ---
  tdr.obj <- .setup_tdr_from_matrix(
    .cells = .cells,
    .source.env = source_env,
    .meta = sample_meta,
    .assay.type = .assay.type,
    .markers = .markers,
    .harmony.var = .harmony.var,
    .celltype.vec = ct_vec,
    .prop.landmarks = .prop.landmarks,
    .seed = .seed,
    .n.threads = .n.threads,
    .verbose = .verbose
  )

  # --- Argument routing ---
  dots <- list(...)

  .lm.formals  <- setdiff(names(formals(get.landmarks.TDRObj)), c("x", "..."))
  .gr.formals  <- setdiff(names(formals(get.graph.TDRObj)), c("x", "..."))
  .map.formals <- setdiff(names(formals(get.map.TDRObj)), c("x", "..."))
  .emb.formals <- setdiff(names(formals(get.embedding.TDRObj)), c("x", "..."))

  lm_args  <- dots[names(dots) %in% .lm.formals]
  gr_args  <- dots[names(dots) %in% .gr.formals]
  map_args <- dots[names(dots) %in% .map.formals]
  emb_args <- dots[names(dots) %in% .emb.formals]

  all_known <- Reduce(union, list(.lm.formals, .gr.formals, .map.formals, .emb.formals))
  orphans <- setdiff(names(dots), all_known)
  if (length(orphans) > 0) {
    warning("Unknown arguments will be ignored: ",
            paste(orphans, collapse = ", "))
  }

  # --- Pipeline ---
  tdr.obj <- do.call(get.landmarks.TDRObj, c(list(tdr.obj, .source = NULL, .seed = .seed, .verbose = .verbose), lm_args))
  tdr.obj <- do.call(get.graph.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), gr_args))

  # --- Import all categorical cell-level annotations ---
  tdr.obj <- import_cell_annotations(tdr.obj,
                                     .cell.meta = .cell.meta,
                                     .sample.var = .sample.var,
                                     .celltype.vec = .celltype.vec,
                                     .verbose = .verbose)

  tdr.obj <- do.call(get.map.TDRObj, c(list(tdr.obj, .source = NULL, .seed = .seed, .verbose = .verbose), map_args))
  tdr.obj <- do.call(get.embedding.TDRObj, c(list(tdr.obj, .seed = .seed, .verbose = .verbose), emb_args))

  return(tdr.obj)
}


# ──────────────────────────────────────────────────────────────────────
# Internal helper: build TDRObj from matrix-derived inputs
# ──────────────────────────────────────────────────────────────────────

#' Build a TDRObj from matrix-derived cell lists and metadata
#'
#' @keywords internal
#' @noRd
.setup_tdr_from_matrix <- function(.cells,
                                   .source.env,
                                   .meta,
                                   .assay.type,
                                   .markers,
                                   .harmony.var,
                                   .celltype.vec,
                                   .prop.landmarks,
                                   .seed,
                                   .n.threads,
                                   .verbose) {

  .assay.type <- match.arg(arg = .assay.type,
                           choices = c("cyto", "RNA"))

  # --- Harmony var validation ---
  if (!is.null(.harmony.var)) {
    if (!inherits(.harmony.var, "character")) {
      stop(".harmony.var must be a character vector.")
    }
    if (!all(.harmony.var %in% colnames(.meta))) {
      stop("Variables not found in metadata: ",
           paste(.harmony.var[!(.harmony.var %in% colnames(.meta))],
                 collapse = ", "),
           "\nCheck column names in .meta with colnames(.meta).")
    }
  }

  # --- Markers validation ---
  if (!is.null(.markers)) {
    if (.assay.type == "RNA") {
      stop(".markers argument only applies to cytometry data.\n",
           "For RNA data, feature selection uses highly variable genes (HVG) automatically.")
    } else if (length(.markers) < 3) {
      stop(".markers must contain at least 3 markers for meaningful dimensionality reduction.")
    }
  }

  if (.assay.type == "cyto" && is.null(.markers)) {
    stop("For cyto assay.type with matrix backend, .markers must be provided.")
  }

  # --- Prop landmarks validation ---
  if ((.prop.landmarks < 0) | (.prop.landmarks > 1)) {
    stop(".prop.landmarks must be between 0 and 1 (e.g., 0.1 for 10% of cells).\n",
         "Current value: ", .prop.landmarks)
  }

  # --- Create TDRObj ---
  .tdr.obj <- TDRObj(
    config = list(
      key = NULL,
      sampling = NULL,
      assay.type = .assay.type,
      markers = NULL,
      n.threads = .n.threads
    ),
    integration = list(
      harmony.var = NULL,
      harmony.obj = NULL
    )
  )

  .tdr.obj@cells <- .cells

  # --- n.cells from integer index vector lengths ---
  n.cells <- lengths(.cells)

  .tdr.obj@config$sampling$n.cells <- n.cells

  # Quality check: warn if sample sizes are highly imbalanced (>10-fold difference)
  if ((max(.tdr.obj@config$sampling$n.cells) /
       min(.tdr.obj@config$sampling$n.cells)) > 10) {

    warning("Sample size imbalance detected: largest/smallest ratio > 10.\n",
            "Smallest sample has ", min(.tdr.obj@config$sampling$n.cells),
            " cells.\n",
            "Consider removing low-quality samples.")

    if (any(.tdr.obj@config$sampling$n.cells < 1000)) {
      warning("Large variation in sample sizes detected. ",
              "For cytometry, samples with <1000 cells may be unreliable.")
    }
  }

  # Calculate target number of landmarks: prop of total cells, capped at 5000
  .tdr.obj@config$sampling$target.lm.n <-
    pmin(sum(.tdr.obj@config$sampling$n.cells) * .prop.landmarks,
         5e3)

  # Allocate landmarks per sample: proportional to sample size, but capped
  .tdr.obj@config$sampling$n.perSample <-
    pmin(ceiling(x = .tdr.obj@config$sampling$n.cells * .prop.landmarks),
         ceiling(x = .tdr.obj@config$sampling$target.lm.n /
                     length(x = .tdr.obj@cells)))

  # Create key vector: maps each future landmark to its sample
  .tdr.obj@config$key <-
    seq_along(along.with = .tdr.obj@cells) |>
    rep(times = .tdr.obj@config$sampling$n.perSample) |>
    (\(x)
     stats::setNames(object = x,
                     nm = names(.tdr.obj@cells)[x])
    )()

  .tdr.obj@metadata <- .meta

  .tdr.obj@metadata$n.perSample <-
    .tdr.obj@config$sampling$n.perSample

  .tdr.obj@metadata$n.cells <-
    .tdr.obj@config$sampling$n.cells

  .tdr.obj@metadata$log10.n.cells <-
    log10(x = .tdr.obj@config$sampling$n.cells)

  # --- Markers ---
  if (.assay.type == "cyto") {
    .tdr.obj@config$markers <- .markers
  }

  # --- Harmony ---
  if (!is.null(.harmony.var)) {
    .tdr.obj@integration$harmony.var <- .harmony.var
  }

  # --- Matrix backend: store locked env for .get_sample_matrix ---
  .tdr.obj@config$backend <- "matrix"
  .tdr.obj@config$source.env <- .source.env

  # --- Cell type vector ---
  if (!is.null(.celltype.vec)) {
    .tdr.obj@config$celltype.vec <- .celltype.vec
  }

  return(.tdr.obj)
}


# ======================================================================
# GetTDR – extractor generic + methods
# ======================================================================

#' Extract a TDRObj from a container object
#'
#' S3 generic that retrieves a \code{\linkS4class{TDRObj}} from the
#' object in which \code{\link{RunTDR}} stored it.
#'
#' @param x An object that may contain a TDRObj.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\linkS4class{TDRObj}}.
#'
#' @export
GetTDR <- function(x, ...) UseMethod("GetTDR")

#' @describeIn GetTDR Default method – returns a TDRObj as-is, errors otherwise
#' @export
GetTDR.default <- function(x, ...) {
  if (is.TDRObj(x)) return(x)
  stop("Cannot extract TDRObj from object of class '",
       paste(class(x), collapse = "/"), "'.")
}

#' @describeIn GetTDR Extract TDRObj from a Seurat object's Misc slot
#' @export
GetTDR.Seurat <- function(x, ...) {
  tdr <- SeuratObject::Misc(x, slot = "tdr.obj")
  if (is.null(tdr))
    stop("No TDRObj found in Seurat Misc slot 'tdr.obj'. ",
         "Run RunTDR() first.")
  tdr
}

#' @describeIn GetTDR Extract TDRObj from a SingleCellExperiment's metadata
#' @export
GetTDR.SingleCellExperiment <- function(x, ...) {
  tdr <- S4Vectors::metadata(x)$tdr.obj
  if (is.null(tdr))
    stop("No TDRObj found in SCE metadata slot 'tdr.obj'. ",
         "Run RunTDR() first.")
  tdr
}

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

#' @describeIn RunTDR Run the pipeline directly on a Seurat object
#'
#' Builds a \code{\linkS4class{TDRObj}} from a Seurat object and executes
#' the full pipeline. The finished TDRObj is stored in
#' \code{SeuratObject::Misc(x, slot = "tdr.obj")}.
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
  .meta <- .meta[names(.cells), , drop = FALSE]

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

  # --- Pipeline ---
  tdr.obj <- get.landmarks(tdr.obj, .source = x, .seed = .seed,
                           .verbose = .verbose, ...)
  tdr.obj <- get.graph(tdr.obj, .seed = .seed,
                       .verbose = .verbose, ...)

  if (!is.null(tdr.obj@config$celltype.vec)) {
    tdr.obj <- celltyping(tdr.obj,
                          .celltyping.map = tdr.obj@config$celltype.vec,
                          .verbose = .verbose)
  }

  tdr.obj <- get.map(tdr.obj, .source = x, .seed = .seed,
                     .verbose = .verbose, ...)

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
#' @param x A \code{SingleCellExperiment} object.
#' @param .sample.var Character(1). Column in \code{colData(x)} identifying
#'   sample membership.
#' @param .assay Character(1). Name of the assay in \code{assayNames(x)}.
#' @param .assay.type Character. \code{"RNA"} or \code{"cyto"}.
#' @param .harmony.var Character vector of batch variable column names
#'   in sample-level metadata, or \code{NULL}.
#' @param .markers Character vector of marker names (required for cyto).
#' @param .celltype.vec Character(1). Column name in \code{colData(x)}
#'   containing per-cell type labels, or \code{NULL}.
#' @param .min.cells.per.sample Integer. Minimum cells for a sample to be
#'   included.
#' @param .optimize.hdf5 Logical. If \code{TRUE} and the assay is
#'   HDF5-backed with poor contiguity, reorder the SCE columns by sample
#'   for faster I/O.
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
                                        .harmony.var = NULL,
                                        .markers = NULL,
                                        .celltype.vec = NULL,
                                        .min.cells.per.sample = 10,
                                        .optimize.hdf5 = FALSE,
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
  .meta <- .meta[names(.cells), , drop = FALSE]

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

  # --- HDF5 contiguity check ---
  is_delayed <- methods::is(
    SummarizedExperiment::assay(x, .assay), "DelayedMatrix")

  if (is_delayed) {
    contiguity <- mean(vapply(.cells, function(idx) {
      if (length(idx) <= 1) return(1)
      sum(diff(idx) == 1L) / (length(idx) - 1L)
    }, numeric(1)))

    if (contiguity < 0.80) {
      if (isFALSE(.optimize.hdf5)) {
        warning("Cells in HDF5-backed assay are not grouped by ",
                "sample. Consider .optimize.hdf5 = TRUE for ",
                "~10-50x faster I/O.", call. = FALSE)
      } else {
        # Reorder SCE by sample for contiguous access
        new_order <- unlist(.cells, use.names = FALSE)
        x <- x[, new_order]
        # Rebuild .cells as contiguous ranges
        sample_ids_new <- SummarizedExperiment::colData(x)[[.sample.var]]
        .cells <- lapply(
          stats::setNames(valid, valid),
          function(s) sort(which(sample_ids_new == s))
        )
        if (!is.null(ct_vec)) {
          # Rebuild ct_vec for reordered cells
          ct_vec <- stats::setNames(
            as.character(
              SummarizedExperiment::colData(x)[[.celltype.vec]]),
            colnames(x)
          )[unlist(.cells, use.names = FALSE)]
        }
      }
    }
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

  # --- Pipeline ---
  tdr.obj <- get.landmarks(tdr.obj, .source = x, .seed = .seed,
                           .verbose = .verbose, ...)
  tdr.obj <- get.graph(tdr.obj, .seed = .seed,
                       .verbose = .verbose, ...)

  if (!is.null(tdr.obj@config$celltype.vec)) {
    tdr.obj <- celltyping(tdr.obj,
                          .celltyping.map = tdr.obj@config$celltype.vec,
                          .verbose = .verbose)
  }

  tdr.obj <- get.map(tdr.obj, .source = x, .seed = .seed,
                     .verbose = .verbose, ...)

  # --- Store TDRObj in SCE metadata ---
  S4Vectors::metadata(x)$tdr.obj <- tdr.obj

  return(x)
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

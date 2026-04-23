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

#' Convert an object to SummarizedExperiment
#'
#' S3 generic for converting objects to
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#'
#' @param x An object to convert.
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @export
as.SummarizedExperiment <- function(x, ...) UseMethod("as.SummarizedExperiment")

#' Convert a TDRObj to SummarizedExperiment
#'
#' Converts a tinydenseR \code{\linkS4class{TDRObj}} into a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#'
#' Rows represent tinydenseR **landmarks** (not genes/proteins).
#' The assays stored are:
#' \describe{
#'   \item{\code{normcounts}}{Size-factor-normalized fuzzy density
#'     (landmarks × samples), from \code{@density$fdens}.}
#'   \item{\code{logcounts}}{log2(normcounts + 0.5), used by
#'     \code{\link{get.lm}} for linear modeling, from \code{@density$Y}.}
#'   \item{\code{counts}}{(When available) Raw fuzzy graph density sums
#'     before size-factor normalization, reconstructed from
#'     \code{@cellmap$fuzzy.graphs}.}
#' }
#'
#' \code{rowData} contains all stored clustering and celltyping solutions.
#' \code{colData} contains sample-level metadata.
#' The full TDRObj is preserved in \code{metadata(se)$tdr.obj} for
#' downstream tinydenseR analysis via \code{\link{GetTDR}}.
#'
#' @param x A \code{\linkS4class{TDRObj}} (must have \code{\link{get.map}}
#'   completed).
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#'
#' @seealso \code{\link{GetTDR}}, \code{\link{RunTDR}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Convert TDRObj to SummarizedExperiment
#' se <- as.SummarizedExperiment(tdr.obj)
#'
#' # Access density data
#' SummarizedExperiment::assay(se, "logcounts")[1:5, 1:5]
#'
#' # Round-trip: extract TDRObj for full tinydenseR analysis
#' tdr <- GetTDR(se)
#' }
as.SummarizedExperiment.TDRObj <- function(x, ...) {

  # ── Fail-fast validation ──────────────────────────────────────────────
  stopifnot(is.TDRObj(x))

  if (is.null(x@density$fdens)) {
    stop("get.map() must be run before converting to SummarizedExperiment.")
  }
  if (is.null(x@density$Y)) {
    stop("get.map() must be run before converting to SummarizedExperiment.")
  }
  if (is.null(x@landmark.annot$clustering$ids)) {
    stop("get.map() must be run before converting to SummarizedExperiment.")
  }
  if (nrow(x@metadata) != ncol(x@density$fdens)) {
    stop("Mismatch: nrow(@metadata) [", nrow(x@metadata),
         "] != ncol(@density$fdens) [", ncol(x@density$fdens), "].")
  }

  # ── Assays ────────────────────────────────────────────────────────────
  normcounts <- as.matrix(x@density$fdens)
  logcounts  <- as.matrix(x@density$Y)

  assay_list <- list(normcounts = normcounts, logcounts = logcounts)

  # Best-effort reconstruction of raw counts from fuzzy graphs
  counts_mat <- .reconstruct_counts(x)
  if (!is.null(counts_mat)) {
    assay_list <- c(list(counts = counts_mat), assay_list)
  }

  # ── rowData ───────────────────────────────────────────────────────────
  row_data <- S4Vectors::DataFrame(
    row.names = rownames(x@density$fdens)
  )

  # Clustering solutions
  cl_names <- setdiff(names(x@landmark.annot$clustering), "ids")
  for (nm in cl_names) {
    row_data[[paste0("clustering_", nm)]] <- x@landmark.annot$clustering[[nm]]
  }
  if (!is.null(x@landmark.annot$clustering$ids)) {
    row_data[["clustering_active"]] <- x@landmark.annot$clustering$ids
  }

  # Celltyping solutions
  ct_names <- setdiff(names(x@landmark.annot$celltyping), "ids")
  for (nm in ct_names) {
    row_data[[paste0("celltyping_", nm)]] <- x@landmark.annot$celltyping[[nm]]
  }
  if (!is.null(x@landmark.annot$celltyping$ids)) {
    row_data[["celltyping_active"]] <- x@landmark.annot$celltyping$ids
  }

  # ── colData ───────────────────────────────────────────────────────────
  col_data <- S4Vectors::DataFrame(
    x@metadata,
    row.names = rownames(x@metadata)
  )

  # ── metadata ──────────────────────────────────────────────────────────
  meta_list <- list(
    tdr.obj             = x,
    tinydenseR_version  = utils::packageVersion("tinydenseR"),
    conversion_date     = Sys.Date()
  )

  # ── Assembly ──────────────────────────────────────────────────────────
  SummarizedExperiment::SummarizedExperiment(
    assays   = assay_list,
    rowData  = row_data,
    colData  = col_data,
    metadata = meta_list
  )
}


# ── Internal helper: reconstruct raw counts from fuzzy graphs ────────────
.reconstruct_counts <- function(x) {
  fg_list <- x@cellmap$fuzzy.graphs
  if (is.null(fg_list) || length(fg_list) == 0L) return(NULL)

  lm_names <- rownames(x@density$fdens)
  n_lm     <- nrow(x@density$fdens)
  smpl_names <- names(fg_list)
  cols <- vector("list", length(smpl_names))
  names(cols) <- smpl_names


  for (sn in smpl_names) {
    fg <- tryCatch(
      {
        rec <- fg_list[[sn]]
        if (is.character(rec)) .tdr_cache_read(rec) else rec
      },
      error = function(e) NULL
    )
    if (is.null(fg)) {
      message("Fuzzy graph cache unavailable for some samples; ",
              "'counts' assay omitted.")
      return(NULL)
    }

    # Determine orientation based on dimnames set at lm.graph.embed.R
    if (identical(rownames(fg), lm_names)) {
      cols[[sn]] <- Matrix::rowSums(fg)
    } else {
      cols[[sn]] <- Matrix::colSums(fg)
    }
  }

  counts_mat <- do.call(cbind, cols)
  rownames(counts_mat) <- lm_names
  as.matrix(counts_mat)
}

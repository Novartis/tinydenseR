#####
# Copyright 2026 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

#' Access fuzzy density matrices from a TDRObj
#'
#' Returns one of the three fuzzy density layers stored in
#' \code{@density} after \code{\link{get.map}} has been run.
#'
#' @param x A \code{\linkS4class{TDRObj}} (or Seurat / SCE / HDF5AnnData
#'   wrapping one via \code{\link{GetTDR}}).
#' @param .which Character(1). Which density layer to return:
#'   \describe{
#'     \item{\code{"raw"}}{Pre-size-factor-normalization fuzzy density sums
#'       (landmarks × samples). Each entry is \eqn{\sum_c F_j[c,l]}, the
#'       total fuzzy edge weight connecting sample \eqn{j}'s cells to
#'       landmark \eqn{l}.}
#'     \item{\code{"norm"}}{Size-factor-normalized fuzzy density
#'       (landmarks × samples). Equals \code{raw / size.factors}
#'       (column-wise). This is the primary analytical layer.
#'       Aliases: \code{"fdens"} (deprecated).}
#'     \item{\code{"log.norm"}}{Log2-transformed normalized density:
#'       \code{log2(norm + 0.5)}. Used by \code{get.lm()} and
#'       \code{get.embedding()}.
#'       Aliases: \code{"Y"} (deprecated).}
#'   }
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric matrix (landmarks × samples), or \code{NULL} if the
#'   requested layer is not available.
#'
#' @details
#' The three density layers form a deterministic chain computed by
#' \code{\link{get.map}}:
#'
#' \deqn{\mathrm{raw}[l,j] = \sum_{c \in j} F_j[c,l]}
#'
#' where \eqn{F_j} is the fuzzy graph (UMAP edge weights) for sample
#' \eqn{j}, \eqn{c} indexes cells, and \eqn{l} indexes landmarks.
#'
#' \deqn{\mathrm{norm} = t(t(\mathrm{raw}) \;/\; \mathrm{size.factors})}
#'
#' with \eqn{\mathrm{size.factors}[j] = n_j / \bar{n}} (cells per sample
#' divided by the mean cell count across samples).
#'
#' \deqn{\mathrm{log.norm} = \log_2(\mathrm{norm} + 0.5)}
#'
#' The pseudo-count of 0.5 stabilizes landmarks with zero density.
#'
#' @section Deprecated aliases:
#' \code{"fdens"} maps to \code{"norm"} and \code{"Y"} maps to
#' \code{"log.norm"}. These aliases are retained for backward
#' compatibility with code written against tinydenseR < 0.0.3.
#'
#' @seealso \code{\link{get.map}} (producer), \code{\link{get.lm}},
#'   \code{\link{get.embedding}}, \code{\link{as.SummarizedExperiment.TDRObj}}
#'
#' @export
#' @examples
#' \dontrun{
#' # After running the pipeline:
#' raw_dens  <- get.density(tdr, "raw")       # pre-normalization
#' norm_dens <- get.density(tdr, "norm")       # size-factor normalized
#' log_dens  <- get.density(tdr, "log.norm")   # log2(norm + 0.5)
#'
#' # Verify invariants:
#' all.equal(norm_dens, t(t(raw_dens) / tdr@density$size.factors))
#' all.equal(log_dens, log2(norm_dens + 0.5))
#' }
get.density <- function(x, ...) UseMethod("get.density")

#' @rdname get.density
#' @export
get.density.TDRObj <- function(x, .which = c("raw", "norm", "log.norm",
                                              "fdens", "Y"), ...) {
  .which <- match.arg(.which)

  # Map deprecated aliases
  .which <- switch(.which,
    "fdens" = "norm",
    "Y"     = "log.norm",
    .which
  )

  switch(.which,
    "raw"      = x@density$raw,
    "norm"     = x@density$norm,
    "log.norm" = x@density$log.norm
  )
}

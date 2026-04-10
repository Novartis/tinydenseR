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

# ======================================================================
# SetTDR generic + methods
# ======================================================================

#' Store a TDRObj inside a container object
#'
#' @param x      The container (Seurat, SCE, or TDRObj).
#' @param tdr    A \code{\linkS4class{TDRObj}} to store.
#' @param ...    Additional arguments (currently unused).
#' @return The container with the TDRObj stored.
#' @export
SetTDR <- function(x, tdr, ...) UseMethod("SetTDR")

#' @describeIn SetTDR Default: if x is a TDRObj, returns tdr directly
#' @export
SetTDR.default <- function(x, tdr, ...) {
  if (is.TDRObj(x)) return(tdr)
  stop("Cannot store TDRObj in object of class '",
       paste(class(x), collapse = "/"), "'.")
}

#' @describeIn SetTDR Store TDRObj in Seurat Misc slot
#' @export
SetTDR.Seurat <- function(x, tdr, ...) {
  SeuratObject::Misc(x, slot = "tdr.obj") <- tdr
  x
}

#' @describeIn SetTDR Store TDRObj in SCE metadata
#' @export
SetTDR.SingleCellExperiment <- function(x, tdr, ...) {
  S4Vectors::metadata(x)$tdr.obj <- tdr
  x
}

# ======================================================================
# get.graph Seurat / SCE wrappers
# ======================================================================

#' @export
get.graph.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.graph.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
get.graph.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.graph.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# celltyping Seurat / SCE wrappers
# ======================================================================

#' @export
celltyping.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- celltyping.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
celltyping.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- celltyping.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# lm.cluster Seurat / SCE wrappers
# ======================================================================

#' @export
lm.cluster.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- lm.cluster.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
lm.cluster.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- lm.cluster.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# recluster Seurat / SCE wrappers
# ======================================================================

#' @export
recluster.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- recluster.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
recluster.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- recluster.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# set_active_clustering Seurat / SCE wrappers
# ======================================================================

#' @export
set_active_clustering.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- set_active_clustering.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
set_active_clustering.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- set_active_clustering.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# set_active_celltyping Seurat / SCE wrappers
# ======================================================================

#' @export
set_active_celltyping.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- set_active_celltyping.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
set_active_celltyping.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- set_active_celltyping.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.features Seurat / SCE wrappers
# ======================================================================

#' @method get.features Seurat
#' @export
get.features.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.features.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @method get.features SingleCellExperiment
#' @export
get.features.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.features.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.landmarks Seurat / SCE wrappers (Tier 2)
# ======================================================================

#' @export
get.landmarks.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.landmarks.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

#' @export
get.landmarks.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.landmarks.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.map Seurat / SCE wrappers (Tier 2)
# ======================================================================

#' @export
get.map.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.map.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

#' @export
get.map.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.map.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.pbDE Seurat / SCE wrappers (Tier 2)
# ======================================================================

#' @export
get.pbDE.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.pbDE.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

#' @export
get.pbDE.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.pbDE.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.markerDE Seurat / SCE wrappers (Tier 2)
# ======================================================================

#' @export
get.markerDE.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.markerDE.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

#' @export
get.markerDE.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.markerDE.TDRObj(tdr, .source = x, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# goi.summary Seurat / SCE wrappers (Tier 2, special: no SetTDR)
# ======================================================================

#' @export
goi.summary.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  goi.summary.TDRObj(tdr, .source = x, ...)
}

#' @export
goi.summary.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  goi.summary.TDRObj(tdr, .source = x, ...)
}

# ======================================================================
# get.lm Seurat / SCE wrappers (Tier 1)
# ======================================================================

#' @export
get.lm.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.lm.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
get.lm.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.lm.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.embedding Seurat / SCE wrappers (Tier 1)
# ======================================================================

#' @export
get.embedding.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.embedding.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
get.embedding.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.embedding.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# get.plsD Seurat / SCE wrappers (Tier 1)
# ======================================================================

#' @export
get.plsD.Seurat <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.plsD.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

#' @export
get.plsD.SingleCellExperiment <- function(x, ...) {
  tdr <- GetTDR(x)
  tdr <- get.plsD.TDRObj(tdr, ...)
  SetTDR(x, tdr)
}

# ======================================================================
# Plot function Seurat wrappers (Phase 3)
# ======================================================================

#' @export
plotPCA.Seurat <- function(x, ...) plotPCA.TDRObj(GetTDR(x), ...)

#' @export
plotUMAP.Seurat <- function(x, ...) plotUMAP.TDRObj(GetTDR(x), ...)

#' @export
plotBeeswarm.Seurat <- function(x, ...) plotBeeswarm.TDRObj(GetTDR(x), ...)

#' @export
plot2Markers.Seurat <- function(x, ...) plot2Markers.TDRObj(GetTDR(x), ...)

#' @export
plotSamplePCA.Seurat <- function(x, ...) plotSamplePCA.TDRObj(GetTDR(x), ...)

#' @export
plotSampleEmbedding.Seurat <- function(x, ...) plotSampleEmbedding.TDRObj(GetTDR(x), ...)

#' @export
plotTradStats.Seurat <- function(x, ...) plotTradStats.TDRObj(GetTDR(x), ...)

#' @export
plotTradPerc.Seurat <- function(x, ...) plotTradPerc.TDRObj(GetTDR(x), ...)

#' @export
plotDensity.Seurat <- function(x, ...) plotDensity.TDRObj(GetTDR(x), ...)

#' @export
plotPbDE.Seurat <- function(x, ...) plotPbDE.TDRObj(GetTDR(x), ...)

#' @export
plotDEA.Seurat <- function(x, ...) plotDEA.TDRObj(GetTDR(x), ...)

#' @export
plotMarkerDE.Seurat <- function(x, ...) plotMarkerDE.TDRObj(GetTDR(x), ...)

#' @export
plotHeatmap.Seurat <- function(x, ...) plotHeatmap.TDRObj(GetTDR(x), ...)

#' @export
plotPlsD.Seurat <- function(x, ...) plotPlsD.TDRObj(GetTDR(x), ...)

#' @export
plotPlsDHeatmap.Seurat <- function(x, ...) plotPlsDHeatmap.TDRObj(GetTDR(x), ...)

# ======================================================================
# Plot function SingleCellExperiment wrappers (Phase 3)
# ======================================================================

#' @export
plotPCA.SingleCellExperiment <- function(x, ...) plotPCA.TDRObj(GetTDR(x), ...)

#' @export
plotUMAP.SingleCellExperiment <- function(x, ...) plotUMAP.TDRObj(GetTDR(x), ...)

#' @export
plotBeeswarm.SingleCellExperiment <- function(x, ...) plotBeeswarm.TDRObj(GetTDR(x), ...)

#' @export
plot2Markers.SingleCellExperiment <- function(x, ...) plot2Markers.TDRObj(GetTDR(x), ...)

#' @export
plotSamplePCA.SingleCellExperiment <- function(x, ...) plotSamplePCA.TDRObj(GetTDR(x), ...)

#' @export
plotSampleEmbedding.SingleCellExperiment <- function(x, ...) plotSampleEmbedding.TDRObj(GetTDR(x), ...)

#' @export
plotTradStats.SingleCellExperiment <- function(x, ...) plotTradStats.TDRObj(GetTDR(x), ...)

#' @export
plotTradPerc.SingleCellExperiment <- function(x, ...) plotTradPerc.TDRObj(GetTDR(x), ...)

#' @export
plotDensity.SingleCellExperiment <- function(x, ...) plotDensity.TDRObj(GetTDR(x), ...)

#' @export
plotPbDE.SingleCellExperiment <- function(x, ...) plotPbDE.TDRObj(GetTDR(x), ...)

#' @export
plotDEA.SingleCellExperiment <- function(x, ...) plotDEA.TDRObj(GetTDR(x), ...)

#' @export
plotMarkerDE.SingleCellExperiment <- function(x, ...) plotMarkerDE.TDRObj(GetTDR(x), ...)

#' @export
plotHeatmap.SingleCellExperiment <- function(x, ...) plotHeatmap.TDRObj(GetTDR(x), ...)

#' @export
plotPlsD.SingleCellExperiment <- function(x, ...) plotPlsD.TDRObj(GetTDR(x), ...)

#' @export
plotPlsDHeatmap.SingleCellExperiment <- function(x, ...) plotPlsDHeatmap.TDRObj(GetTDR(x), ...)

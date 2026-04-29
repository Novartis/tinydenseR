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

library(tinydenseR)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(BPCells)
library(ggpubr)
library(ggh4x)
library(rstatix)
library(scales)
library(symphony)

# set wd to local dir
script.path <-
  grep(pattern = "^--file=", 
       x = commandArgs(), 
       value = TRUE,
       fixed = FALSE) |>
  gsub(pattern = "^--file=",
       replacement = "", 
       fixed = FALSE) |>
  (\(x)
   if(length(x = x) == 0) rstudioapi::getSourceEditorContext()$path else x
  )()

script.path |>
  dirname() |>
  setwd()

# create sub-folder for results
rd <- file.path(dirname(script.path), "results", "head_to_head_Seu")

if(!dir.exists(paths = rd)) dir.create(path = rd, recursive = TRUE)

# create sub-folder to store data
dd <- file.path(dirname(script.path), "derived_data")

if(!dir.exists(paths = dd)) dir.create(path = dd, recursive = TRUE)

# ── Required external data ──
# This script requires a merged Seurat v5 object with BPCells backend,
# saved as "derived_data/seurat_object.RDS" relative to this script.
#
# To create this object, follow:
#   1. https://github.com/satijalab/seurat/blob/e35abd442520808a20025e589f861620ddc315af/vignettes/seurat5_bpcells_interaction_vignette.Rmd
#   2. https://github.com/satijalab/seurat/blob/30f82df52159ac5f0feb80b149698abbd876b779/vignettes/COVID_SCTMapping.Rmd
#      (up to line 83)
# Then save the resulting object:
#   saveRDS(object, file.path("derived_data", "seurat_object.RDS"))

# load merged Seurat objects
object <-
  readRDS(file = file.path(dd, "seurat_object.RDS"))

# make sure to have the ref.umap reduction
Seurat::DimPlot(object = object, 
                reduction = "ref.umap", 
                group.by = "predicted.celltype.l2", 
                alpha = 0.1, 
                label = TRUE, 
                split.by = "publication",
                ncol = 3,
                label.size = 3) + 
  Seurat::NoLegend()

# use donor_id and disease to split samples
object$donor_id_disease <-
  paste0(object$donor_id,
         "_",
         object$disease) |>
  make.names()

# get cell type proportions per sample
df_comp <- 
  as.data.frame.matrix(x = table(object$donor_id_disease, 
                                 object$predicted.celltype.l2))

# select samples with > 1000 cells
select.donors <- 
  rownames(x = df_comp)[table(object$donor_id_disease) > 1000]
df_comp <- 
  df_comp[select.donors, ]

# convert to relative abundance
df_comp_relative <-
  sweep(x = df_comp, 
        MARGIN = 1, 
        STATS = rowSums(x = df_comp), 
        FUN = "/")

# add disease info
df_disease <-
  as.data.frame.matrix(x = table(object$donor_id_disease, 
                                 object$disease))[select.donors, ]

df_comp_relative$disease <-
  "other"

df_comp_relative$disease[df_disease$normal != 0] <-
  "normal"

df_comp_relative$disease[df_disease$`COVID-19` != 0] <- 
  "COVID-19"

df_comp_relative$disease <-
  factor(x = df_comp_relative$disease, 
         levels = c("normal", "COVID-19", "other"))

# keep only normal and COVID-19 samples
df_comp_relative <-
  df_comp_relative[df_comp_relative$disease %in% c("normal", "COVID-19"), ]

# statistical testing
COVID.stat.percentages <-
  tidyr::pivot_longer(data = df_comp_relative,
                      cols = -disease,
                      names_to = "cell_type",
                      values_to = "relative_abundance") |> 
  dplyr::mutate(cell_type = gsub(pattern = "Proliferating",
                                 replacement = "prolif.",
                                 x = cell_type,
                                 fixed = TRUE) |>
                  gsub(pattern = "intermediate",
                       replacement = "int.",
                       fixed = TRUE) |>
                  gsub(pattern = "_",
                       replacement = " ",
                       fixed = TRUE) |>
                  gsub(pattern = "bright",
                       replacement = "++",
                       fixed = TRUE)) |>
  dplyr::group_by(cell_type) |>
  rstatix::t_test(formula = relative_abundance ~ disease) |>
  rstatix::adjust_pvalue(method = "fdr") |>
  rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1)) |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "disease") |>
  dplyr::mutate(y.position = y.position*100)

# format relative abundances
COVID.PBMC.Seurat.perc.data <-
  dplyr::mutate(.data = df_comp_relative,
                donor_id_disease = rownames(x = df_comp_relative)) |>
  tidyr::pivot_longer(cols = -c(disease,donor_id_disease),
                      names_to = "cell_type",
                      values_to = "relative_abundance") |> 
  dplyr::mutate(cell_type = gsub(pattern = "Proliferating",
                                 replacement = "prolif.",
                                 x = cell_type,
                                 fixed = TRUE) |>
                  gsub(pattern = "intermediate",
                       replacement = "int.",
                       fixed = TRUE) |>
                  gsub(pattern = "_",
                       replacement = " ",
                       fixed = TRUE) |>
                  gsub(pattern = "bright",
                       replacement = "++",
                       fixed = TRUE))

# plot relative abundances
(p <-
    COVID.PBMC.Seurat.perc.data |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = disease, 
                                          y = relative_abundance*100)) +
     ggplot2::theme_bw() +
     ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                        vjust = 1,
                                                        angle = 30),
                    axis.title.x = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(hjust = 0.5)) +
     ggplot2::facet_wrap(facets = ~ cell_type, 
                         ncol = unique(x = x$cell_type) |>
                           length() |>
                           sqrt() |>
                           ceiling(),
                         scales = "free_y") +
     ggplot2::labs(title = "Seurat: COVID PBMC",
                   y = "% of cells") +
     ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
     ggplot2::geom_boxplot(outlier.shape = NA) +
     ggplot2::geom_point(size = I(x = 0.1),
                         position = ggplot2::position_jitter(width = 0.2,
                                                             height = 0,
                                                             seed = 123)) +
     ggpubr::stat_pvalue_manual(data = COVID.stat.percentages, 
                                label = "p.adj.signif",
                                label.size = I(x = 3)) +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 0.75,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 0.75,
                                                  units = "in"))
  )()); ggplot2::ggsave(
    plot = p,
    filename = file.path(rd,
                         "COVID.PBMC.Seurat.perc.png"),
    width = 7,
    height = 7,
    dpi = 300); rm(p)

# Build symphony reference to use with tinydenseR
reference <-
  symphony::buildReference(
    exp_ref = reference[["SCT"]]$counts, # ONLY SCT FOUND IN THIS DATASET!
    metadata_ref = object@meta.data,
    verbose = TRUE,
    save_uwot_path = file.path(dd, "symphony_ref_uwot_model"), # file path to save uwot UMAP model
    seed = 123
  )

# Use the new tinydenseR API to create .meta and .cells objects
# get.meta.Seurat5 extracts sample-level metadata from Seurat v5 objects
# get.cells.Seurat5 handles multiple BPCells layers, gene intersection, and cell filtering

# Create .meta object using new API
# Note: We need to add the metadata columns we want to keep to the Seurat object first
.meta <-
  tinydenseR::get.meta.Seurat5(
    .seurat.obj = object,
    .sample.var = "donor_id_disease",
    .verbose = TRUE
  )

# Create .cells object using new API
# This automatically handles:
# - Multiple BPCells layers matching "counts" pattern
# - Gene intersection across layers
# - Filtering samples with >= 1000 cells
.cells <-
  tinydenseR::get.cells.Seurat5(
    .seurat.obj = object,
    .meta = .meta,
    .sample.var = "donor_id_disease",
    .assay = "RNA",
    .layer.pattern = "counts",
    .layer.name.source = "publication",  # Use publication to name layers
    .min.cells.per.sample = 1000,  # Filter to samples with > 1000 cells
    .compress = FALSE,
    .verbose = TRUE,
    .dest.path = dd
  )

# Subset .meta to match .cells (samples that passed filtering)
.meta <- 
  .meta[names(x = .cells), ] |>
  droplevels()

covid.lm.cells <-
  tinydenseR::setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .harmony.var = "publication",
    .assay.type = "RNA",
    .verbose = TRUE
  ) |>
  tinydenseR::get.landmarks(
    .verbose = TRUE)

covid.lm.cells <-
  tinydenseR::get.graph(
    x = covid.lm.cells,
    .cl.resolution.parameter = 2,
    .verbose = TRUE)

tinydenseR::plotUMAP(
  x = covid.lm.cells,
  .panel.size = 2
)

tinydenseR::plotUMAP(
  x = covid.lm.cells,
  .feature = covid.lm.cells$metadata$publication[covid.lm.cells$config$key],
  .panel.size = 2
)

tinydenseR::plotPCA(
  x = covid.lm.cells,
  .feature = covid.lm.cells$metadata$publication[covid.lm.cells$config$key],
  .panel.size = 2
)

covid.lm.cells <-
  tinydenseR::get.map(
    x = covid.lm.cells,
    .ref.obj = reference,
    .celltype.col.name = "celltype.l2",
    .verbose = TRUE) |>
  get.embedding()

tinydenseR::plotUMAP(
  x = covid.lm.cells,
  .feature = covid.lm.cells$landmark.annot$celltyping$ids,
  .panel.size = 2
)

tinydenseR::plotSampleEmbedding(
  x = covid.lm.cells,
  .embedding = "pca",
  .color.by = "disease",
  .point.size = 1
)

tinydenseR::plotSampleEmbedding(
  x = covid.lm.cells,
  .embedding = "pca",
  .color.by = "publication",
  .point.size = 1
)

tinydenseR::plotSampleEmbedding(
  x = covid.lm.cells,
  .embedding = "pca",
  .color.by = "sex",
  .point.size = 1
)

COVID.tdr.stat.percentages <-
  as.data.frame(x = covid.lm.cells$map$celltyping$cell.perc[,colnames(x = covid.lm.cells$map$celltyping$cell.perc) != "..low.confidence.."]) |>
  cbind(disease = covid.lm.cells$metadata$disease) |>
  dplyr::filter(disease %in% c("normal",
                               "COVID-19")) |>
  droplevels() |>
  dplyr::mutate(disease = factor(x = disease,
                                 levels = c("normal", "COVID-19"))) |>
  tidyr::pivot_longer(cols = -disease,
                      names_to = "cell_type",
                      values_to = "relative_abundance") |> 
  dplyr::mutate(cell_type = gsub(pattern = "Proliferating",
                                 replacement = "prolif.",
                                 x = cell_type,
                                 fixed = TRUE) |>
                  gsub(pattern = "intermediate",
                       replacement = "int.",
                       fixed = TRUE) |>
                  gsub(pattern = "_",
                       replacement = " ",
                       fixed = TRUE) |>
                  gsub(pattern = "bright",
                       replacement = "++",
                       fixed = TRUE)) |>
  dplyr::group_by(cell_type) |>
  rstatix::t_test(formula = relative_abundance ~ disease) |>
  rstatix::adjust_pvalue(method = "fdr") |>
  rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1)) |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "disease")

COVID.PBMC.perc.data <- 
  as.data.frame(x = covid.lm.cells$map$celltyping$cell.perc[,colnames(x = covid.lm.cells$map$celltyping$cell.perc) != "..low.confidence.."]) |>
  cbind(disease = covid.lm.cells$metadata$disease,
        donor_id_disease = rownames(x = covid.lm.cells$map$celltyping$cell.perc)) |>
  dplyr::filter(disease %in% c("normal",
                               "COVID-19")) |>
  dplyr::filter(disease %in% c("normal",
                               "COVID-19")) |>
  droplevels() |>
  dplyr::mutate(disease = factor(x = disease,
                                 levels = c("normal", "COVID-19"))) |>
  tidyr::pivot_longer(cols = -c(disease, donor_id_disease),
                      names_to = "cell_type",
                      values_to = "relative_abundance") |> 
  dplyr::mutate(cell_type = gsub(pattern = "Proliferating",
                                 replacement = "prolif.",
                                 x = cell_type,
                                 fixed = TRUE) |>
                  gsub(pattern = "intermediate",
                       replacement = "int.",
                       fixed = TRUE) |>
                  gsub(pattern = "_",
                       replacement = " ",
                       fixed = TRUE) |>
                  gsub(pattern = "bright",
                       replacement = "++",
                       fixed = TRUE))

COVID.PBMC.perc <-
  COVID.PBMC.perc.data |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = disease, 
                                          y = relative_abundance)) +
     ggplot2::theme_bw() +
     ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                        vjust = 1,
                                                        angle = 30),
                    axis.title.x = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(hjust = 0.5)) +
     ggplot2::facet_wrap(facets = ~ cell_type, 
                         ncol = unique(x = x$cell_type) |>
                           length() |>
                           sqrt() |>
                           ceiling(),
                         scales = "free_y") +
     ggplot2::labs(title = "tinydenseR: COVID PBMC",
                   y = "% of cells") +
     ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
     ggplot2::geom_boxplot(outlier.shape = NA) +
     ggplot2::geom_point(size = I(x = 0.1),
                         position = ggplot2::position_jitter(width = 0.2,
                                                             height = 0,
                                                             seed = 123)) +
     ggpubr::stat_pvalue_manual(data = COVID.tdr.stat.percentages, 
                                label = "p.adj.signif",
                                label.size = I(x = 3)) +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 0.75,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 0.75,
                                                  units = "in"))
  )()

COVID.PBMC.perc

ggplot2::ggsave(plot = COVID.PBMC.perc,
                filename = file.path(rd,
                                     "COVID.PBMC.perc.png"),
                width = 7,
                height = 7,
                dpi = 300)

(p <-
  list(tinydenseR = COVID.PBMC.perc.data,
     Seurat = dplyr::mutate(.data = COVID.PBMC.Seurat.perc.data,
                            relative_abundance = relative_abundance*100)) |>
  dplyr::bind_rows(.id = "method") |>
  #dplyr::mutate(relative_abundance = log1p(x = relative_abundance)) |>
  tidyr::pivot_wider(names_from = method,
                     values_from = relative_abundance) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = Seurat,
                                          y = tinydenseR)) +
     ggplot2::theme_bw() +
     ggplot2::theme(axis.text = ggplot2::element_text(size = I(x = 6)),
                    plot.title = ggplot2::element_text(hjust = 0.5)) +
     ggplot2::facet_wrap(facets = ~ cell_type,
                         ncol = unique(x = x$cell_type) |>
                           length() |>
                           sqrt() |>
                           ceiling(),
                         scales = "free") +
     ggplot2::labs(title = "COVID PBMC: Seurat vs tinydenseR",
                   x = "Seurat (% of cells)",
                   y = "tinydenseR (% of cells)") +
     scale_x_continuous(transform = "log1p") +
     scale_y_continuous(transform = "log1p") +
     ggplot2::geom_point(size = I(x = 0.1)) +
     ggplot2::geom_abline(slope = 1,
                          intercept = 0,
                          linetype = "dashed",
                          color = "red") +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 0.75,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 0.75,
                                                  units = "in"))
  )()); ggplot2::ggsave(
    plot = p,
    filename = file.path(rd,
                                             "COVID.PBMC.1to1.perc.png"),
                        width = 7,
                        height = 8,
                        dpi = 300); rm(p)

sessionInfo()
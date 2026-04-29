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

library(ggrepel)
library(patchwork)
library(tinydenseR)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggh4x)
library(Matrix)
library(flowCore)
library(flowWorkspace)

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
rd <- file.path(dirname(script.path), "results", "sim_flow_DE")

if(!dir.exists(paths = rd)) dir.create(path = rd, recursive = TRUE)

set.seed(seed = 42)

# Simulate DE data and write FCS files
sim_data <- tinydenseR::simulate_DE_data()
final_data_DE <- sim_data$cell_meta

DE_res_DE <-
  final_data_DE |>
  dplyr::group_by(Sample, Treatment, Batch, SD_Shift, CellType) |>
  dplyr::summarize(n = dplyr::n()) |> 
  dplyr::mutate(prop = n /sum(n)) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |> 
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$SD_Shift) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(setting) {
            lm(formula = log2(x = prop) ~ Treatment + Batch,
               data = x[x$SD_Shift == setting,]) |> 
              (\(mod)
               c(coef = unname(obj = coef(object = mod)["TreatmentActivation"]),
                 confint(object = mod)["TreatmentActivation",])
              )()
          })
  )() |>
  dplyr::bind_rows(.id = "id") 

DE_res_DE |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = coef,
                                          y = id)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
     ggplot2::geom_point(size = I(x = 1)) +
     ggplot2::geom_vline(xintercept = 0,
                         linetype = "dotted",
                         color = "red") +
     ggplot2::geom_errorbar(mapping = ggplot2::aes(xmin = `2.5 %`,
                                                   xmax = `97.5 %`), 
                            width = 0.2) +
     ggplot2::labs(title = "DEA",
                   x = "log2FC in %",
                   y = "") +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 1,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 1,
                                                  units = "in"))
  )()

# 0.5SD
.setting.meta.0.5 <-
  sim_data$sample_meta |>
  dplyr::filter(SD_Shift == "0.5SD")

cs.DE.0.5 <-
  flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(.setting.meta.0.5$fcs_path,
                            .setting.meta.0.5$Sample))

flowWorkspace::pData(cs.DE.0.5)$Sample <-
  flowWorkspace::sampleNames(cs.DE.0.5)

flowWorkspace::pData(cs.DE.0.5)$Treatment <-
  .setting.meta.0.5$Treatment[match(flowWorkspace::sampleNames(cs.DE.0.5),
                                    paste0(.setting.meta.0.5$Sample,
                                           ".fcs"))]

flowWorkspace::pData(cs.DE.0.5)$Batch <-
  .setting.meta.0.5$Batch[match(flowWorkspace::sampleNames(cs.DE.0.5),
                                paste0(.setting.meta.0.5$Sample,
                                       ".fcs"))]

set.seed(seed = 123)
lm.cells.DE.0.5 <-
  tinydenseR::RunTDR(
    x = cs.DE.0.5,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

# Unsupervised sample embedding
(smpl.pca.DE.0.5 <- 
    tinydenseR::plotSampleEmbedding(
      x = lm.cells.DE.0.5,
      .embedding = "pca",
      .color.by = "Treatment",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(2,1)],
      .panel.size = 1.5,
      .point.size = 3) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom"))

(p <- tinydenseR::plotSampleEmbedding(
  x = lm.cells.DE.0.5,
  .embedding = "pca",
  .color.by = "Batch",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(3,4)],
  .panel.size = 1.5,
  .point.size = 3
) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "left")); ggplot2::ggsave(plot = p,
                                                               filename = file.path(rd,
                                                                                    "sample_pca_DE_0.5_Batch.png"),
                                                               width = 3.5, 
                                                               height = 2.5,
                                                               dpi = 300,
                                                               bg = "white"); rm(p)

# supervised analysis
.meta.DE.0.5 <-
  lm.cells.DE.0.5@metadata |>
  dplyr::mutate(Treatment = factor(x = Treatment, levels = c("Baseline", "Activation")))

.design.0.5 <-
  model.matrix(object = ~ Treatment + Batch,
               data = .meta.DE.0.5) |> 
  (\(x)
   `colnames<-`(x = x,
                value = colnames(x = x) |>
                  gsub(pattern = "^Treatment|^Batch",
                       replacement = "",
                       fixed = FALSE))
  )()

lm.cells.DE.0.5 <-
  tinydenseR::get.lm(
    x = lm.cells.DE.0.5,
    .design = .design.0.5)

# Create reduced models to embed samples quantitatively along the variables of interest
noTreatment.design.0.5 <- 
  model.matrix(object = ~ Batch,
               data = lm.cells.DE.0.5@metadata)

noBatch.design.0.5 <- 
  model.matrix(object = ~ Treatment,
               data = lm.cells.DE.0.5@metadata)

lm.cells.DE.0.5 <- 
  tinydenseR::get.lm(
    x = lm.cells.DE.0.5,
    .design = noTreatment.design.0.5,
    .model.name = "noTreatment",
    .verbose = TRUE
  ) |>
  tinydenseR::get.lm(
    .design = noBatch.design.0.5,
    .model.name = "noBatch",
    .verbose = TRUE 
  )

# update stats results to get sample embedding
lm.cells.DE.0.5 <-
  tinydenseR::get.embedding(
    x = lm.cells.DE.0.5,
    .full.model = "default",
    .red.model = "noTreatment",
    .term.of.interest = "Treatment",
    .verbose = FALSE 
  ) |>
  tinydenseR::get.embedding(
    .full.model = "default",
    .red.model = "noBatch",
    .term.of.interest = "Batch",
    .verbose = FALSE 
  )

# Embed samples based on differences along the variable of interest
(smpl.pePC.DE.0.5 <- 
    tinydenseR::plotSampleEmbedding(
      x = lm.cells.DE.0.5,
      .embedding = "pePC",
      .sup.embed.slot = "Treatment",
      .color.by = "Treatment",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(2,1)],
      .panel.size = 1.5,
      .point.size = 3) +
    ggplot2::labs(title = "pePC") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom"))

(p <-
    ((smpl.pca.DE.0.5 | 
        smpl.pePC.DE.0.5) +
       patchwork::plot_layout(guides = "collect") &
       ggplot2::theme(legend.position = "bottom",
                      legend.justification = "center")));ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "sample_embed_DE_0.5.png"),
                        width = 4.5,
                        height = 3,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

(p <- tinydenseR::plotSampleEmbedding(
  x = lm.cells.DE.0.5,
  .embedding = "pePC",
  .sup.embed.slot = "Batch",
  .color.by = "Batch",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(3,4)],
  .panel.size = 1.5,
  .point.size = 2
) +
    ggplot2::labs(title = "pePC: Batch") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none")); ggplot2::ggsave(plot = p,
                                                               filename = file.path(rd,
                                                                                    "sample_pePC_DE_0.5.Batch.png"),
                                                               width = 2.5, 
                                                               height = 2.5,
                                                               dpi = 300,
                                                               bg = "white"); rm(p)

(p <- lapply(X = names(lm.cells.DE.0.5@cells),
             FUN = function(s) flowCore::exprs(cs.DE.0.5[[s]])) |>
    do.call(what = rbind) |>
    (\(x)
     (((Matrix::t(x = x[,lm.cells.DE.0.5@landmark.embed$pca$HVG]) - lm.cells.DE.0.5@landmark.embed$pca$center) /
         lm.cells.DE.0.5@landmark.embed$pca$scale) |>
         Matrix::t()) %*%
       lm.cells.DE.0.5@landmark.embed$pca$rotation
    )() |>
    as.matrix() |>
    as.data.frame() |>
    (\(x)
     dplyr::mutate(
       .data = x,
       Treatment = final_data_DE |>
         dplyr::filter(Sample %in% gsub(pattern = "\\.fcs$", replacement = "", names(x = lm.cells.DE.0.5@cells))) |>
         dplyr::pull(Treatment))
    )() |>
    (\(x)
     ggplot2::ggplot(data = x,
                     mapping = ggplot2::aes(x = PC1,
                                            y = PC2)) +
       ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
       ggplot2::theme_bw() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      legend.position = "none") +
       ggplot2::labs(title = "ground truth") +
       ggplot2::stat_bin_hex(bins = 128) + 
       ggplot2::scale_fill_viridis_c(trans = "log") +
       ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                                 units = "in"),
                               rows = grid::unit(x = 2,
                                                 units = "in"))
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd,
                                               "ground_truth_DE_0.5.png"),
                          width = 5, 
                          height = 3,
                          dpi = 300); rm(p)

# similar to above, but now a simple heatmap of ground truth mean marker expression in target vs other
(p <- 
  flowCore::fsApply(
  x = lm.cells.DE.0.5$config$source.env$cs,
  FUN = flowCore::exprs) |>
  as.data.frame() |>
  cbind(final_data_DE[(final_data_DE$SD_Shift == "0.5SD"),]) |>
  dplyr::filter(CellType == "target") |>
  dplyr::group_by(Sample, CellType) |>
  dplyr::summarize(
    dplyr::across(dplyr::starts_with(match = "Marker"), mean),
    Treatment = unique(x = Treatment),
    .groups = "drop") |>
  tidyr::pivot_longer(cols = dplyr::starts_with("Marker"),
                      names_to = "name",
                      values_to = "mean\nexpr") |>
  ggplot2::ggplot(mapping = ggplot2::aes(x = Treatment,
                                         y = name,
                                         fill = `mean\nexpr`)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        hjust = 0.5,
                                                        vjust = 1),
                   axis.title = ggplot2::element_blank(),
                   legend.position = "right",
                  axis.text.y = ggplot2::element_text()) +
    ggplot2::labs() +
    ggplot2::scale_y_discrete(position = "right", limits = rev) +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::geom_raster() +
    ggh4x::force_panelsizes(cols = grid::unit(x = 0.5,
                                              units = "in"),
                            rows = grid::unit(x = 1.5,
                                              units = "in"))); ggplot2::ggsave(
       filename = file.path(rd,
                            "ground_truth_markers_DE_0.5.png"),
       width = 2, 
       height = 2.5,
       dpi = 300); rm(p)

tinydenseR::plotPCA(x = lm.cells.DE.0.5,
                    .feature = lm.cells.DE.0.5@metadata$Treatment[lm.cells.DE.0.5@config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(DE.0.5.dens <- tinydenseR::plotPCA(x = lm.cells.DE.0.5,
                                    .feature = lm.cells.DE.0.5@results$lm$default$fit$coefficients[,"Activation"],
                                    .plot.title = "Activation vs Baseline",
                                    .color.label = "density\nlog2(+0.5)FC",
                                    .panel.size = 1.5,
                                    .point.size = 1,
                                    .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(DE.0.5.q <- tinydenseR::plotPCA(
  x = lm.cells.DE.0.5,
  .feature =
    ifelse(
      test = lm.cells.DE.0.5@results$lm$default$fit$coefficients[,"Activation"] < 0,
      yes = "lower density",
      no = "higher density") |>
    ifelse(
      test = lm.cells.DE.0.5@results$lm$default$fit$pca.weighted.q[,"Activation"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("lower density",
                      "higher density",
                      "not sig.")),
  .plot.title = "Activation vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2,6)],
  .point.size = 1,
  .panel.size = 1.5,
  .legend.position = "bottom")   +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(cl.DE.0.5 <- tinydenseR::plotPCA(
  x = lm.cells.DE.0.5,
  .feature = lm.cells.DE.0.5@landmark.annot$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 1.5,
  .legend.position = "bottom") |> 
    (\(x)
     x +
       ggplot2::guides(color = ggplot2::guide_legend(ncol = 2,
                                                     override.aes = list(size = I(x = 5)))) +
       ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
       ggrepel::geom_text_repel(data = x$data |>
                            dplyr::group_by(feature) |>
                            dplyr::summarize(PC1 = mean(x = PC1),
                                             PC2 = mean(x = PC2),
                                             .groups = "drop"),
                          mapping = ggplot2::aes(label = feature),
                          size = I(x = 4),
                          color = "black",
                          seed = 123, min.segment.length = -Inf)
    )())

(p <- tinydenseR::plotTradStats(
  x = lm.cells.DE.0.5,
  .model.name = "default")); ggplot2::ggsave(plot = p,
                                             filename = file.path(rd, 
                                                                  "TradStats_DE_0.5.png"),
                                             width = 7, 
                                             height = 5,
                                             dpi = 300); rm(p)

stat.test.percentages.DE.0.5 <-
  lm.cells.DE.0.5@density$composition$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.0.5@metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DE.0.5@results$lm$default$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = lm.cells.DE.0.5@results$lm$default$trad$clustering$fit$adj.p[name,"Activation"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "treatment")

(p <- tinydenseR::plotTradPerc(
  x = lm.cells.DE.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.0.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))); ggplot2::ggsave(plot = p,
                                                                                            filename = file.path(rd, 
                                                                                                                 "TradPerc_DE_0.5.png"),
                                                                                            width = 3.5, 
                                                                                            height = 3,
                                                                                            dpi = 300); rm(p)

(p <- tinydenseR::plotDensity(
  x = lm.cells.DE.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density")); ggplot2::ggsave(
      plot = p,
                                                                      filename = file.path(rd, 
                                                                                           "density_DE_0.5.png"),
                                                                      width = 3.5, 
                                                                      height = 3,
                                                                      dpi = 300); rm(p)

(DE.0.5.bees <- tinydenseR::plotBeeswarm(
  x = lm.cells.DE.0.5,
  .model.name = "default",
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE,
  .q.from = "pca.weighted.q",
  .legend.position = "bottom") +
    ggplot2::geom_vline(xintercept = -1,
                        color = "red",
                        linetype = "dashed"))

lm.cells.DE.0.5 <-
  tinydenseR::get.plsD(
    x = lm.cells.DE.0.5, 
    .coef.col = "Activation",
    .model.name = "default",
    .verbose = TRUE
  )

(p <-
    tinydenseR::plotPlsD(
      x = lm.cells.DE.0.5,
      .coef.col = "Activation",
      .point.size = 2
    )); rm(p)

(DE.0.5.plsD1 <-
    tinydenseR::plotPlsD(
      x = lm.cells.DE.0.5,
      .coef.col = "Activation",
      .plsD.dim = 1,
      .embed = "pca",
      .panel.size = 1.5
    )[[1]])

(p <-
    ((DE.0.5.dens +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                        title.hjust = 0.5)) +
        ggplot2::labs(title = "") +
        ggplot2::theme(legend.position = "bottom",
                       legend.margin = ggplot2::margin(t = -0.1, 
                                                       unit = "in"))) | 
       (DE.0.5.plsD1 +
          ggplot2::labs(title = "") +
          ggplot2::theme(legend.margin = ggplot2::margin(t = 0.1, 
                                                         unit = "in")))) +
    patchwork::plot_annotation(title = "Density contrast: Activation") &
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(hjust = 0.5,
                              margin = ggplot2::margin(t = -0.1, 
                                                       unit = "in")))); ggplot2::ggsave(
                                                         plot = p,
                                                         filename = file.path(rd,
                                                                              "q_DE_0.5.png"),
                                                         width = 5.5,
                                                         height = 3.5,
                                                         units = "in",
                                                         dpi = 300,
                                                         bg = "white");rm(p)

(p <-
    tinydenseR::plotPlsDHeatmap(
      x = lm.cells.DE.0.5,
      .coef.col = "Activation",
      .plsD.dim = 1,
      .order.by = "plsD.dim",
      .panel.height = 2,
      .feature.font.size = I(x = 8)
    )); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "sim.DE.0.5.plsD1.hm.png"),
                        width = 6.5, 
                        height = 4,
                        dpi = 300,
                        bg = "white"); rm(p) 

lm.cells.DE.0.5 <-
  tinydenseR::get.pbDE(
    x = lm.cells.DE.0.5,
    .design = .design.0.5,
    .id = "cluster.4"
  )

(markerDE.DE.0.5 <- 
    tinydenseR::plotPbDE(
      x = lm.cells.DE.0.5,
      .population.name = "cluster.4",
      .coefs = "Activation",
      .row.space.scaler = 0.065,
      .col.space.scaler = 0.2,
      .order.by = "none") +
    ggplot2::labs(title = "cluster.4: Activation",
                  subtitle = "pseudobulk") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                       vjust = 1,
                                                       angle = 30),
                   axis.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank()) +
    ggplot2::guides(fill = guide_colorbar(position = "bottom", 
                                          title.position = "top",
                                          title.hjust = 0.5)) +
    ggplot2::coord_flip())

(p <-
    (cl.DE.0.5 | 
        markerDE.DE.0.5));ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "DEA_DE_0.5.png"),
                        width = 5,
                        height = 3.75,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

# 1SD
.setting.meta.1 <-
  sim_data$sample_meta |>
  dplyr::filter(SD_Shift == "1SD")

cs.DE.1 <-
  flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(.setting.meta.1$fcs_path,
                            .setting.meta.1$Sample))

flowWorkspace::pData(cs.DE.1)$Sample <-
  flowWorkspace::sampleNames(cs.DE.1)

flowWorkspace::pData(cs.DE.1)$Treatment <-
  .setting.meta.1$Treatment[match(flowWorkspace::sampleNames(cs.DE.1),
                                    paste0(.setting.meta.1$Sample,
                                           ".fcs"))]

flowWorkspace::pData(cs.DE.1)$Batch <-
  .setting.meta.1$Batch[match(flowWorkspace::sampleNames(cs.DE.1),
                                paste0(.setting.meta.1$Sample,
                                       ".fcs"))]

set.seed(seed = 123)
lm.cells.DE.1 <-
  tinydenseR::RunTDR(
    x = cs.DE.1,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

# Unsupervised sample embedding
(smpl.pca.DE.1 <- 
    tinydenseR::plotSampleEmbedding(
      x = lm.cells.DE.1,
      .embedding = "pca",
      .color.by = "Treatment",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(2,1)],
      .panel.size = 1.5,
      .point.size = 3) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom"))

(p <- tinydenseR::plotSampleEmbedding(
  x = lm.cells.DE.1,
  .embedding = "pca",
  .color.by = "Batch",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(3,4)],
  .panel.size = 1.5,
  .point.size = 3
) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "left")); ggplot2::ggsave(plot = p,
                                                               filename = file.path(rd,
                                                                                    "sample_pca_DE_1_Batch.png"),
                                                               width = 3.5, 
                                                               height = 2.5,
                                                               dpi = 300,
                                                               bg = "white"); rm(p)

# supervised analysis
.meta.DE.1 <-
  lm.cells.DE.1@metadata |>
  dplyr::mutate(Treatment = factor(x = Treatment, levels = c("Baseline", "Activation")))

.design.1 <-
  model.matrix(object = ~ Treatment + Batch,
               data = .meta.DE.1) |> 
  (\(x)
   `colnames<-`(x = x,
                value = colnames(x = x) |>
                  gsub(pattern = "^Treatment|^Batch",
                       replacement = "",
                       fixed = FALSE))
  )()

lm.cells.DE.1 <-
  tinydenseR::get.lm(
    x = lm.cells.DE.1,
    .design = .design.1)

# Create reduced models to embed samples quantitatively along the variables of interest
noTreatment.design.1 <- 
  model.matrix(object = ~ Batch,
               data = lm.cells.DE.1@metadata)

noBatch.design.1 <- 
  model.matrix(object = ~ Treatment,
               data = lm.cells.DE.1@metadata)

lm.cells.DE.1 <- 
  tinydenseR::get.lm(
    x = lm.cells.DE.1,
    .design = noTreatment.design.1,
    .model.name = "noTreatment",
    .verbose = TRUE
  ) |>
  tinydenseR::get.lm(
    .design = noBatch.design.1,
    .model.name = "noBatch",
    .verbose = TRUE 
  )

# update stats results to get sample embedding
lm.cells.DE.1 <-
  tinydenseR::get.embedding(
    x = lm.cells.DE.1,
    .full.model = "default",
    .red.model = "noTreatment",
    .term.of.interest = "Treatment",
    .verbose = FALSE 
  ) |>
  tinydenseR::get.embedding(
    .full.model = "default",
    .red.model = "noBatch",
    .term.of.interest = "Batch",
    .verbose = FALSE 
  )

# Embed samples based on differences along the variable of interest
(smpl.pePC.DE.1 <- 
    tinydenseR::plotSampleEmbedding(
      x = lm.cells.DE.1,
      .embedding = "pePC",
      .sup.embed.slot = "Treatment",
      .color.by = "Treatment",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(2,1)],
      .panel.size = 1.5,
      .point.size = 3) +
    ggplot2::labs(title = "pePC") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom"))

(p <-
    ((smpl.pca.DE.1 | 
        smpl.pePC.DE.1) +
       patchwork::plot_layout(guides = "collect") &
       ggplot2::theme(legend.position = "bottom",
                      legend.justification = "center")));ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "sample_embed_DE_1.png"),
                        width = 4.5,
                        height = 3,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

(p <- tinydenseR::plotSampleEmbedding(
  x = lm.cells.DE.1,
  .embedding = "pePC",
  .sup.embed.slot = "Batch",
  .color.by = "Batch",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(3,4)],
  .panel.size = 1.5,
  .point.size = 2
) +
    ggplot2::labs(title = "pePC: Batch") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none")); ggplot2::ggsave(plot = p,
                                                               filename = file.path(rd,
                                                                                    "sample_pePC_DE_1.Batch.png"),
                                                               width = 2.5, 
                                                               height = 2.5,
                                                               dpi = 300,
                                                               bg = "white"); rm(p)

(p <- lapply(X = names(lm.cells.DE.1@cells),
             FUN = function(s) flowCore::exprs(cs.DE.1[[s]])) |>
    do.call(what = rbind) |>
    (\(x)
     (((Matrix::t(x = x[,lm.cells.DE.1@landmark.embed$pca$HVG]) - lm.cells.DE.1@landmark.embed$pca$center) /
         lm.cells.DE.1@landmark.embed$pca$scale) |>
         Matrix::t()) %*%
       lm.cells.DE.1@landmark.embed$pca$rotation
    )() |>
    as.matrix() |>
    as.data.frame() |>
    (\(x)
     dplyr::mutate(
       .data = x,
       Treatment = final_data_DE |>
         dplyr::filter(Sample %in% gsub(pattern = "\\.fcs$", replacement = "", names(x = lm.cells.DE.1@cells))) |>
         dplyr::pull(Treatment))
    )() |>
    (\(x)
     ggplot2::ggplot(data = x,
                     mapping = ggplot2::aes(x = PC1,
                                            y = PC2)) +
       ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
       ggplot2::theme_bw() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      legend.position = "none") +
       ggplot2::labs(title = "ground truth") +
       ggplot2::stat_bin_hex(bins = 128) + 
       ggplot2::scale_fill_viridis_c(trans = "log") +
       ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                                 units = "in"),
                               rows = grid::unit(x = 2,
                                                 units = "in"))
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd,
                                               "ground_truth_DE_1.png"),
                          width = 5, 
                          height = 3,
                          dpi = 300); rm(p)

# similar to above, but now a simple heatmap of ground truth mean marker expression in target vs other
(p <- 
  flowCore::fsApply(
  x = lm.cells.DE.1$config$source.env$cs,
  FUN = flowCore::exprs) |>
  as.data.frame() |>
  cbind(final_data_DE[(final_data_DE$SD_Shift == "1SD"),]) |>
  dplyr::filter(CellType == "target") |>
  dplyr::group_by(Sample, CellType) |>
  dplyr::summarize(
    dplyr::across(dplyr::starts_with(match = "Marker"), mean),
    Treatment = unique(x = Treatment),
    .groups = "drop") |>
  tidyr::pivot_longer(cols = dplyr::starts_with("Marker"),
                      names_to = "name",
                      values_to = "mean\nexpr") |>
  ggplot2::ggplot(mapping = ggplot2::aes(x = Treatment,
                                         y = name,
                                         fill = `mean\nexpr`)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        hjust = 0.5,
                                                        vjust = 1),
                   axis.title = ggplot2::element_blank(),
                   legend.position = "right",
                  axis.text.y = ggplot2::element_text()) +
    ggplot2::labs() +
    ggplot2::scale_y_discrete(position = "right", limits = rev) +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::geom_raster() +
    ggh4x::force_panelsizes(cols = grid::unit(x = 0.5,
                                              units = "in"),
                            rows = grid::unit(x = 1.5,
                                              units = "in"))); ggplot2::ggsave(
       filename = file.path(rd,
                            "ground_truth_markers_DE_1.png"),
       width = 2, 
       height = 2.5,
       dpi = 300); rm(p)


tinydenseR::plotPCA(x = lm.cells.DE.1,
                    .feature = lm.cells.DE.1@metadata$Treatment[lm.cells.DE.1@config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(DE.1.dens <- tinydenseR::plotPCA(x = lm.cells.DE.1,
                                    .feature = lm.cells.DE.1@results$lm$default$fit$coefficients[,"Activation"],
                                    .plot.title = "Activation vs Baseline",
                                    .color.label = "density\nlog2(+0.5)FC",
                                    .panel.size = 1.5,
                                    .point.size = 1,
                                    .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(DE.1.q <- tinydenseR::plotPCA(
  x = lm.cells.DE.1,
  .feature =
    ifelse(
      test = lm.cells.DE.1@results$lm$default$fit$coefficients[,"Activation"] < 0,
      yes = "lower density",
      no = "higher density") |>
    ifelse(
      test = lm.cells.DE.1@results$lm$default$fit$pca.weighted.q[,"Activation"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("lower density",
                      "higher density",
                      "not sig.")),
  .plot.title = "Activation vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2,6)],
  .point.size = 1,
  .panel.size = 1.5,
  .legend.position = "bottom")   +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(cl.DE.1 <- tinydenseR::plotPCA(
  x = lm.cells.DE.1,
  .feature = lm.cells.DE.1@landmark.annot$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 1.5,
  .legend.position = "bottom") |> 
    (\(x)
     x +
       ggplot2::guides(color = ggplot2::guide_legend(ncol = 2,
                                                     override.aes = list(size = I(x = 5)))) +
       ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
       ggrepel::geom_text_repel(data = x$data |>
                            dplyr::group_by(feature) |>
                            dplyr::summarize(PC1 = mean(x = PC1),
                                             PC2 = mean(x = PC2),
                                             .groups = "drop"),
                          mapping = ggplot2::aes(label = feature),
                          size = I(x = 4),
                          color = "black",
                          seed = 123)
    )())

(p <- tinydenseR::plotTradStats(
  x = lm.cells.DE.1,
  .model.name = "default")); ggplot2::ggsave(plot = p,
                                             filename = file.path(rd, 
                                                                  "TradStats_DE_1.png"),
                                             width = 7, 
                                             height = 5,
                                             dpi = 300); rm(p)

stat.test.percentages.DE.1 <-
  lm.cells.DE.1@density$composition$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.1@metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DE.1@results$lm$default$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = lm.cells.DE.1@results$lm$default$trad$clustering$fit$adj.p[name,"Activation"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "treatment")

(p <- tinydenseR::plotTradPerc(
  x = lm.cells.DE.1,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.1, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))); ggplot2::ggsave(plot = p,
                                                                                            filename = file.path(rd, 
                                                                                                                 "TradPerc_DE_1.png"),
                                                                                            width = 3.5, 
                                                                                            height = 3,
                                                                                            dpi = 300); rm(p)

(p <- tinydenseR::plotDensity(
  x = lm.cells.DE.1,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density")); ggplot2::ggsave(
      plot = p,
                                                                      filename = file.path(rd, 
                                                                                           "density_DE_1.png"),
                                                                      width = 3.5, 
                                                                      height = 3,
                                                                      dpi = 300); rm(p)

(DE.1.bees <- tinydenseR::plotBeeswarm(
  x = lm.cells.DE.1,
  .model.name = "default",
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE,
  .q.from = "pca.weighted.q",
  .legend.position = "bottom") +
    ggplot2::geom_vline(xintercept = -1,
                        color = "red",
                        linetype = "dashed"))

lm.cells.DE.1 <-
  tinydenseR::get.plsD(
    x = lm.cells.DE.1, 
    .coef.col = "Activation",
    .model.name = "default",
    .verbose = TRUE
  )

(p <-
    tinydenseR::plotPlsD(
      x = lm.cells.DE.1,
      .coef.col = "Activation",
      .point.size = 2
    )); rm(p)

(DE.1.plsD1 <-
    tinydenseR::plotPlsD(
      x = lm.cells.DE.1,
      .coef.col = "Activation",
      .plsD.dim = 1,
      .embed = "pca",
      .panel.size = 1.5
    )[[1]])

(p <-
    ((DE.1.dens +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                        title.hjust = 0.5)) +
        ggplot2::labs(title = "") +
        ggplot2::theme(legend.position = "bottom",
                       legend.margin = ggplot2::margin(t = -0.1, 
                                                       unit = "in"))) | 
       (DE.1.plsD1 +
          ggplot2::labs(title = "") +
          ggplot2::theme(legend.margin = ggplot2::margin(t = 0.1, 
                                                         unit = "in")))) +
    patchwork::plot_annotation(title = "Density contrast: Activation") &
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(hjust = 0.5,
                              margin = ggplot2::margin(t = -0.1, 
                                                       unit = "in")))); ggplot2::ggsave(
                                                         plot = p,
                                                         filename = file.path(rd,
                                                                              "q_DE_1.png"),
                                                         width = 5.5,
                                                         height = 3.5,
                                                         units = "in",
                                                         dpi = 300,
                                                         bg = "white");rm(p)

(p <-
    tinydenseR::plotPlsDHeatmap(
      x = lm.cells.DE.1,
      .coef.col = "Activation",
      .plsD.dim = 1,
      .order.by = "plsD.dim",
      .panel.height = 2,
      .feature.font.size = I(x = 8)
    )); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "sim.DE.1.plsD1.hm.png"),
                        width = 6.5, 
                        height = 4,
                        dpi = 300,
                        bg = "white"); rm(p) 

lm.cells.DE.1 <-
  tinydenseR::get.pbDE(
    x = lm.cells.DE.1,
    .design = .design.1,
    .id = "cluster.3"
  )

(markerDE.DE.1 <- 
    tinydenseR::plotPbDE(
      x = lm.cells.DE.1,
      .population.name = "cluster.3",
      .coefs = "Activation",
      .row.space.scaler = 0.065,
      .col.space.scaler = 0.2,
      .order.by = "none") +
    ggplot2::labs(title = "cluster.3: Activation",
                  subtitle = "pseudobulk") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                       vjust = 1,
                                                       angle = 30),
                   axis.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank()) +
    ggplot2::guides(fill = guide_colorbar(position = "bottom", 
                                          title.position = "top",
                                          title.hjust = 0.5)) +
    ggplot2::coord_flip())

(p <-
    (cl.DE.1 | 
        markerDE.DE.1));ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "DEA_DE_1.png"),
                        width = 5,
                        height = 3.75,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

# 2SD
.setting.meta.2 <-
  sim_data$sample_meta |>
  dplyr::filter(SD_Shift == "2SD")

cs.DE.2 <-
  flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(.setting.meta.2$fcs_path,
                            .setting.meta.2$Sample))

flowWorkspace::pData(cs.DE.2)$Sample <-
  flowWorkspace::sampleNames(cs.DE.2)

flowWorkspace::pData(cs.DE.2)$Treatment <-
  .setting.meta.2$Treatment[match(flowWorkspace::sampleNames(cs.DE.2),
                                    paste0(.setting.meta.2$Sample,
                                           ".fcs"))]

flowWorkspace::pData(cs.DE.2)$Batch <-
  .setting.meta.2$Batch[match(flowWorkspace::sampleNames(cs.DE.2),
                                paste0(.setting.meta.2$Sample,
                                       ".fcs"))]

set.seed(seed = 123)
lm.cells.DE.2 <-
  tinydenseR::RunTDR(
    x = cs.DE.2,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

# Unsupervised sample embedding
(smpl.pca.DE.2 <- 
    tinydenseR::plotSampleEmbedding(
      x = lm.cells.DE.2,
      .embedding = "pca",
      .color.by = "Treatment",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(2,1)],
      .panel.size = 1.5,
      .point.size = 3) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom"))

(p <- tinydenseR::plotSampleEmbedding(
  x = lm.cells.DE.2,
  .embedding = "pca",
  .color.by = "Batch",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(3,4)],
  .panel.size = 1.5,
  .point.size = 3
) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "left")); ggplot2::ggsave(plot = p,
                                                               filename = file.path(rd,
                                                                                    "sample_pca_DE_2_Batch.png"),
                                                               width = 3.5, 
                                                               height = 2.5,
                                                               dpi = 300,
                                                               bg = "white"); rm(p)

# supervised analysis
.meta.DE.2 <-
  lm.cells.DE.2@metadata |>
  dplyr::mutate(Treatment = factor(x = Treatment, levels = c("Baseline", "Activation")))

.design.2 <-
  model.matrix(object = ~ Treatment + Batch,
               data = .meta.DE.2) |> 
  (\(x)
   `colnames<-`(x = x,
                value = colnames(x = x) |>
                  gsub(pattern = "^Treatment|^Batch",
                       replacement = "",
                       fixed = FALSE))
  )()

lm.cells.DE.2 <-
  tinydenseR::get.lm(
    x = lm.cells.DE.2,
    .design = .design.2)

# Create reduced models to embed samples quantitatively along the variables of interest
noTreatment.design.2 <- 
  model.matrix(object = ~ Batch,
               data = lm.cells.DE.2@metadata)

noBatch.design.2 <- 
  model.matrix(object = ~ Treatment,
               data = lm.cells.DE.2@metadata)

lm.cells.DE.2 <- 
  tinydenseR::get.lm(
    x = lm.cells.DE.2,
    .design = noTreatment.design.2,
    .model.name = "noTreatment",
    .verbose = TRUE
  ) |>
  tinydenseR::get.lm(
    .design = noBatch.design.2,
    .model.name = "noBatch",
    .verbose = TRUE 
  )

# update stats results to get sample embedding
lm.cells.DE.2 <-
  tinydenseR::get.embedding(
    x = lm.cells.DE.2,
    .full.model = "default",
    .red.model = "noTreatment",
    .term.of.interest = "Treatment",
    .verbose = FALSE 
  ) |>
  tinydenseR::get.embedding(
    .full.model = "default",
    .red.model = "noBatch",
    .term.of.interest = "Batch",
    .verbose = FALSE 
  )

# Embed samples based on differences along the variable of interest
(smpl.pePC.DE.2 <- 
    tinydenseR::plotSampleEmbedding(
      x = lm.cells.DE.2,
      .embedding = "pePC",
      .sup.embed.slot = "Treatment",
      .color.by = "Treatment",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(2,1)],
      .panel.size = 1.5,
      .point.size = 3) +
    ggplot2::labs(title = "pePC") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom"))

(p <-
    ((smpl.pca.DE.2 | 
        smpl.pePC.DE.2) +
       patchwork::plot_layout(guides = "collect") &
       ggplot2::theme(legend.position = "bottom",
                      legend.justification = "center")));ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "sample_embed_DE_2.png"),
                        width = 4.5,
                        height = 3,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

(p <- tinydenseR::plotSampleEmbedding(
  x = lm.cells.DE.2,
  .embedding = "pePC",
  .sup.embed.slot = "Batch",
  .color.by = "Batch",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(3,4)],
  .panel.size = 1.5,
  .point.size = 2
) +
    ggplot2::labs(title = "pePC: Batch") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none")); ggplot2::ggsave(plot = p,
                                                               filename = file.path(rd,
                                                                                    "sample_pePC_DE_2.Batch.png"),
                                                               width = 2.5, 
                                                               height = 2.5,
                                                               dpi = 300,
                                                               bg = "white"); rm(p)

(p <- lapply(X = names(lm.cells.DE.2@cells),
             FUN = function(s) flowCore::exprs(cs.DE.2[[s]])) |>
    do.call(what = rbind) |>
    (\(x)
     (((Matrix::t(x = x[,lm.cells.DE.2@landmark.embed$pca$HVG]) - lm.cells.DE.2@landmark.embed$pca$center) /
         lm.cells.DE.2@landmark.embed$pca$scale) |>
         Matrix::t()) %*%
       lm.cells.DE.2@landmark.embed$pca$rotation
    )() |>
    as.matrix() |>
    as.data.frame() |>
    (\(x)
     dplyr::mutate(
       .data = x,
       Treatment = final_data_DE |>
         dplyr::filter(Sample %in% gsub(pattern = "\\.fcs$", replacement = "", names(x = lm.cells.DE.2@cells))) |>
         dplyr::pull(Treatment))
    )() |>
    (\(x)
     ggplot2::ggplot(data = x,
                     mapping = ggplot2::aes(x = PC1,
                                            y = PC2)) +
       ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
       ggplot2::theme_bw() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      legend.position = "none") +
       ggplot2::labs(title = "ground truth") +
       ggplot2::stat_bin_hex(bins = 128) + 
       ggplot2::scale_fill_viridis_c(trans = "log") +
       ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                                 units = "in"),
                               rows = grid::unit(x = 2,
                                                 units = "in"))
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd,
                                               "ground_truth_DE_2.png"),
                          width = 5, 
                          height = 3,
                          dpi = 300); rm(p)

# similar to above, but now a simple heatmap of ground truth mean marker expression in target vs other
(p <- 
  flowCore::fsApply(
  x = lm.cells.DE.2$config$source.env$cs,
  FUN = flowCore::exprs) |>
  as.data.frame() |>
  cbind(final_data_DE[(final_data_DE$SD_Shift == "2SD"),]) |>
  dplyr::filter(CellType == "target") |>
  dplyr::group_by(Sample, CellType) |>
  dplyr::summarize(
    dplyr::across(dplyr::starts_with(match = "Marker"), mean),
    Treatment = unique(x = Treatment),
    .groups = "drop") |>
  tidyr::pivot_longer(cols = dplyr::starts_with("Marker"),
                      names_to = "name",
                      values_to = "mean\nexpr") |>
  ggplot2::ggplot(mapping = ggplot2::aes(x = Treatment,
                                         y = name,
                                         fill = `mean\nexpr`)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        hjust = 0.5,
                                                        vjust = 1),
                   axis.title = ggplot2::element_blank(),
                   legend.position = "right",
                  axis.text.y = ggplot2::element_text()) +
    ggplot2::labs() +
    ggplot2::scale_y_discrete(position = "right", limits = rev) +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::geom_raster() +
    ggh4x::force_panelsizes(cols = grid::unit(x = 0.5,
                                              units = "in"),
                            rows = grid::unit(x = 1.5,
                                              units = "in"))); ggplot2::ggsave(
       filename = file.path(rd,
                            "ground_truth_markers_DE_2.png"),
       width = 2, 
       height = 2.5,
       dpi = 300); rm(p)

tinydenseR::plotPCA(x = lm.cells.DE.2,
                    .feature = lm.cells.DE.2@metadata$Treatment[lm.cells.DE.2@config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(DE.2.dens <- tinydenseR::plotPCA(x = lm.cells.DE.2,
                                    .feature = lm.cells.DE.2@results$lm$default$fit$coefficients[,"Activation"],
                                    .plot.title = "Activation vs Baseline",
                                    .color.label = "density\nlog2(+0.5)FC",
                                    .panel.size = 1.5,
                                    .point.size = 1,
                                    .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(DE.2.q <- tinydenseR::plotPCA(
  x = lm.cells.DE.2,
  .feature =
    ifelse(
      test = lm.cells.DE.2@results$lm$default$fit$coefficients[,"Activation"] < 0,
      yes = "lower density",
      no = "higher density") |>
    ifelse(
      test = lm.cells.DE.2@results$lm$default$fit$pca.weighted.q[,"Activation"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("lower density",
                      "higher density",
                      "not sig.")),
  .plot.title = "Activation vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2,6)],
  .point.size = 1,
  .panel.size = 1.5,
  .legend.position = "bottom")   +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(cl.DE.2 <- tinydenseR::plotPCA(
  x = lm.cells.DE.2,
  .feature = lm.cells.DE.2@landmark.annot$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 1.5,
  .legend.position = "bottom") |> 
    (\(x)
     x +
       ggplot2::guides(color = ggplot2::guide_legend(ncol = 2,
                                                     override.aes = list(size = I(x = 5)))) +
       ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
       ggrepel::geom_text_repel(data = x$data |>
                            dplyr::group_by(feature) |>
                            dplyr::summarize(PC1 = mean(x = PC1),
                                             PC2 = mean(x = PC2),
                                             .groups = "drop"),
                          mapping = ggplot2::aes(label = feature),
                          size = I(x = 4),
                          color = "black",
                          seed = 123)
    )())

(p <- tinydenseR::plotTradStats(
  x = lm.cells.DE.2,
  .model.name = "default")); ggplot2::ggsave(plot = p,
                                             filename = file.path(rd, 
                                                                  "TradStats_DE_2.png"),
                                             width = 7, 
                                             height = 5,
                                             dpi = 300); rm(p)

stat.test.percentages.DE.2 <-
  lm.cells.DE.2@density$composition$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.2@metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DE.2@results$lm$default$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = lm.cells.DE.2@results$lm$default$trad$clustering$fit$adj.p[name,"Activation"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "treatment")

(p <- tinydenseR::plotTradPerc(
  x = lm.cells.DE.2,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.2, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))); ggplot2::ggsave(plot = p,
                                                                                            filename = file.path(rd, 
                                                                                                                 "TradPerc_DE_2.png"),
                                                                                            width = 3.5, 
                                                                                            height = 3,
                                                                                            dpi = 300); rm(p)

(p <- tinydenseR::plotDensity(
  x = lm.cells.DE.2,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density")); ggplot2::ggsave(
      plot = p,
                                                                      filename = file.path(rd, 
                                                                                           "density_DE_2.png"),
                                                                      width = 3.5, 
                                                                      height = 3,
                                                                      dpi = 300); rm(p)

(DE.2.bees <- tinydenseR::plotBeeswarm(
  x = lm.cells.DE.2,
  .model.name = "default",
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE,
  .q.from = "pca.weighted.q",
  .legend.position = "bottom") +
    ggplot2::geom_vline(xintercept = -1,
                        color = "red",
                        linetype = "dashed"))

lm.cells.DE.2 <-
  tinydenseR::get.plsD(
    x = lm.cells.DE.2, 
    .coef.col = "Activation",
    .model.name = "default",
    .verbose = TRUE
  )

(p <-
    tinydenseR::plotPlsD(
      x = lm.cells.DE.2,
      .coef.col = "Activation",
      .point.size = 2
    )); rm(p)

(DE.2.plsD1 <-
    tinydenseR::plotPlsD(
      x = lm.cells.DE.2,
      .coef.col = "Activation",
      .plsD.dim = 1,
      .embed = "pca",
      .panel.size = 1.5
    )[[1]])

(p <-
    ((DE.2.dens +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                        title.hjust = 0.5)) +
        ggplot2::labs(title = "") +
        ggplot2::theme(legend.position = "bottom",
                       legend.margin = ggplot2::margin(t = -0.1, 
                                                       unit = "in"))) | 
       (DE.2.plsD1 +
          ggplot2::labs(title = "") +
          ggplot2::theme(legend.margin = ggplot2::margin(t = 0.1, 
                                                         unit = "in")))) +
    patchwork::plot_annotation(title = "Density contrast: Activation") &
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(hjust = 0.5,
                              margin = ggplot2::margin(t = -0.1, 
                                                       unit = "in")))); ggplot2::ggsave(
                                                         plot = p,
                                                         filename = file.path(rd,
                                                                              "q_DE_2.png"),
                                                         width = 5.5,
                                                         height = 3.5,
                                                         units = "in",
                                                         dpi = 300,
                                                         bg = "white");rm(p)

(p <-
    tinydenseR::plotPlsDHeatmap(
      x = lm.cells.DE.2,
      .coef.col = "Activation",
      .plsD.dim = 1,
      .order.by = "plsD.dim",
      .panel.height = 2,
      .feature.font.size = I(x = 8)
    )); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "sim.DE.2.plsD1.hm.png"),
                        width = 6.5, 
                        height = 4,
                        dpi = 300,
                        bg = "white"); rm(p) 

lm.cells.DE.2 <-
  tinydenseR::get.pbDE(
    x = lm.cells.DE.2,
    .design = .design.2,
    .id = "cluster.4"
  )

(markerDE.DE.2 <- 
    tinydenseR::plotPbDE(
      x = lm.cells.DE.2,
      .population.name = "cluster.4",
      .coefs = "Activation",
      .row.space.scaler = 0.065,
      .col.space.scaler = 0.2,
      .order.by = "none") +
    ggplot2::labs(title = "cluster.4: Activation",
                  subtitle = "pseudobulk") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                       vjust = 1,
                                                       angle = 30),
                   axis.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank()) +
    ggplot2::guides(fill = guide_colorbar(position = "bottom", 
                                          title.position = "top",
                                          title.hjust = 0.5)) +
    ggplot2::coord_flip())

(p <-
    (cl.DE.2 | 
        markerDE.DE.2));ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "DEA_DE_2.png"),
                        width = 5,
                        height = 3.75,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

# permutation tests
file.path(script.path |>
  dirname(),
"perm_utils.R") |>
  source()

# For each condition (0.5SD, 1SD, 2SD), run exact permutation test
conditions <- list(
  "0.5SD" = list(lm_obj = lm.cells.DE.0.5,
                 meta = .meta.DE.0.5),
  "1SD"   = list(lm_obj = lm.cells.DE.1, 
                 meta = .meta.DE.1),
  "2SD"   = list(lm_obj = lm.cells.DE.2, 
                 meta = .meta.DE.2)
)

results_list <- 
  names(x = conditions) |>
  lapply(FUN = function(nm) {
    
    print(nm)
    
    cond <- 
      conditions[[nm]]
    
    res  <-
      run_stratified_permutation_test(lm_obj = cond$lm_obj, 
                                      meta_df = cond$meta, 
                                      coef_name = "Activation")
    res$condition <-
      nm
    
    res
  })

# Row-bind
results <-
  dplyr::bind_rows(results_list)

# Sanity check MCC values
cat("\n=== MCC Validation ===\n")
for (cond_name in names(conditions)) {
  cond_results <- results[results$condition == cond_name, ]
  cat(sprintf("\n%s:\n", cond_name))
  cat("  Observed MCC:", unique(cond_results$mcc[cond_results$type == "observed"]), "\n")
  cat("  Complement MCC:", unique(cond_results$mcc[cond_results$type == "complement"]), "\n")
  cat("  MCC range for permuted:", 
      paste(range(cond_results$mcc[cond_results$type == "permuted"], na.rm = TRUE), collapse = " to "), "\n")
}

# Verify expected values across all conditions
stopifnot(all(results$mcc[results$type == "observed"] == 1.0, na.rm = TRUE))
stopifnot(all(abs(results$mcc[results$type == "complement"] - (-1.0)) < 1e-10, na.rm = TRUE))
cat("\n✓ MCC validation passed!\n\n")

# Quick summary
results |>
  group_by(condition, type) |>
  summarise(mean_sig = mean(n_sig),
            mean_mcc = mean(mcc, na.rm = TRUE),
            .groups = "drop")

# Plot results
# Plot results
(p <- 
  plot_n_sig(results = results)); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "perm_test_n_sig.png"),
                        width = 9, 
                        height = 4,
                        dpi = 300); rm(p)
(p <-
  plot_min_q(results = results)); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "perm_test_min_q.png"),
                        width = 9, 
                        height = 4,
                        dpi = 300); rm(p)

results[results$mcc == 0,] |>
  dplyr::group_by(condition, type) |>
  dplyr::summarize(mean_min_q = mean(min_q < 0.1),
                   .groups = "drop")

sessionInfo()
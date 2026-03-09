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
library(ggpubr)
library(rstatix)
library(patchwork)
library(ggh4x)
library(Matrix)
library(flowCore)
library(flowWorkspace)

# Simulate DA data and write FCS files
sim_data <- simulate_DA_data()
final_data_DA <- sim_data$cell_meta

DA_res_DA <-
  final_data_DA |>
  dplyr::group_by(Sample, Treatment, Batch, Setting, CellType) |>
  dplyr::summarize(n = dplyr::n()) |> 
  dplyr::mutate(prop = n /sum(n)) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |> 
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$Setting) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(setting) {
            lm(formula = log2(x = prop) ~ Treatment + Batch,
               data = x[x$Setting == setting,]) |> 
              (\(mod)
               c(coef = unname(obj = coef(object = mod)["TreatmentDepletion"]),
                 confint(object = mod)["TreatmentDepletion",])
              )()
          })
  )() |>
  dplyr::bind_rows(.id = "id") 

DA_res_DA |>
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
     ggplot2::labs(title = "DAA",
                   x = "log2FC in %",
                   y = "") +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 1,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 1,
                                                  units = "in"))
  )()

# 0.5%
.setting.meta.0.5 <-
  sim_data$sample_meta |>
  dplyr::filter(Setting == "0.5%")

cs.DA.0.5 <-
  flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(.setting.meta.0.5$fcs_path,
                            .setting.meta.0.5$Sample))

flowWorkspace::pData(cs.DA.0.5)$Sample <-
  flowWorkspace::sampleNames(cs.DA.0.5)
flowWorkspace::pData(cs.DA.0.5)$Treatment <-
  .setting.meta.0.5$Treatment[match(flowWorkspace::sampleNames(cs.DA.0.5),
                                    .setting.meta.0.5$Sample)]
flowWorkspace::pData(cs.DA.0.5)$Batch <-
  .setting.meta.0.5$Batch[match(flowWorkspace::sampleNames(cs.DA.0.5),
                                .setting.meta.0.5$Sample)]

set.seed(seed = 123)
lm.cells.DA.0.5 <-
  tinydenseR::RunTDR(
    cs.DA.0.5,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

.meta.DA.0.5 <-
  lm.cells.DA.0.5@metadata

.design.0.5 <-
  model.matrix(object = ~ Treatment + Batch,
               data = .meta.DA.0.5) |> 
  (\(x)
   `colnames<-`(x = x,
                value = colnames(x = x) |>
                  gsub(pattern = "^Treatment|^Batch",
                       replacement = "",
                       fixed = FALSE))
  )()

# New API: get.lm() returns updated .tdr.obj with results in $results$lm[[.model.name]]
lm.cells.DA.0.5 <-
  tinydenseR::get.lm(
    .tdr.obj = lm.cells.DA.0.5,
    .design = .design.0.5)

lapply(X = names(lm.cells.DA.0.5@cells),
       FUN = function(s) flowCore::exprs(cs.DA.0.5[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.0.5$landmark.embed$pca$HVG]) - lm.cells.DA.0.5$landmark.embed$pca$center) /
       lm.cells.DA.0.5$landmark.embed$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.0.5$landmark.embed$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.0.5@cells)) |>
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
  )()

final_data_DA |>
  dplyr::filter(Sample %in% names(x = lm.cells.DA.0.5@cells)) |>
  dplyr::mutate(CellType = factor(x = CellType,
                                  levels = c("target",
                                             "other"))) |>
  tidyr::pivot_longer(
    cols = c(Marker1, Marker2, Marker3, Marker4, Marker5),
    names_to = "Marker",
    values_to = "value"
  ) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = Marker,
                                          y = value,
                                          color = CellType,
                                          fill = CellType)) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "right",
                    axis.text.x = ggplot2::element_text(hjust = 1,
                                                        vjust = 1,
                                                        angle = 30)) +
     ggplot2::labs(title = "ground truth (no expression diff. with depletion)",
                   x = "",
                   y = "expression level") +
     ggplot2::scale_fill_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) + 
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                       alpha = 1)),
                     fill = "none") +
     ggplot2::scale_color_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) + 
     ggplot2::scale_y_log10(breaks = scales::trans_breaks(trans = "log10",
                                                          inv =  function(x) 10^x),
                            labels = scales::trans_format(trans = "log10", 
                                                          format = scales::math_format(10^.x))) + 
     ggplot2::annotation_logticks(sides = "l") +
     scattermore::geom_scattermore(position = ggplot2::position_jitterdodge(jitter.width = 0.5,
                                                                            jitter.height = 0,
                                                                            dodge.width = 0.5,
                                                                            seed = 123),
                                   alpha = 0.1) + 
     ggplot2::geom_violin(color = "black",
                          alpha = 0,
                          position = ggplot2::position_dodge(width = 0.5),
                          quantiles = 0.5, 
                          quantile.linetype = "solid") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = names(lm.cells.DA.0.5@cells),
       FUN = function(s) flowCore::exprs(cs.DA.0.5[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.0.5$landmark.embed$pca$HVG]) - lm.cells.DA.0.5$landmark.embed$pca$center) /
       lm.cells.DA.0.5$landmark.embed$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.0.5$landmark.embed$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.0.5@cells)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.0.5@cells)) |>
       dplyr::pull(Batch))
  )() |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2)) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "right") +
     ggplot2::labs(title = "ground truth",
                   color = "") +
     scattermore::geom_scattermore(mapping = ggplot2::aes(color = Batch),
                                   pointsize = I(x = 3)) + 
     ggplot2::scale_color_viridis_d() +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

tinydenseR::plotPCA(.tdr.obj = lm.cells.DA.0.5,
                    .feature = lm.cells.DA.0.5$metadata$Treatment[lm.cells.DA.0.5$config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.tdr.obj = lm.cells.DA.0.5,
                     .feature = lm.cells.DA.0.5$results$lm$default$fit$coefficients[,"Depletion"],
                     .plot.title = "Depletion vs Baseline",
                     .color.label = "density\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DA.0.5,
  .feature =
    ifelse(
      test = lm.cells.DA.0.5$results$lm$default$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = lm.cells.DA.0.5$results$lm$default$fit$pca.weighted.q[,"Depletion"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("less abundant",
                      "not sig.",
                      "more abundant")),
  .plot.title = "Depletion vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2)   +
    ggplot2::labs(subtitle = "hypothesis testing"))

tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DA.0.5,
  .feature = lm.cells.DA.0.5$landmark.annot$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 2) |> 
  (\(x)
   x +
     ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
     ggplot2::geom_text(data = x$data |>
                          dplyr::group_by(feature) |>
                          dplyr::summarize(PC1 = mean(x = PC1),
                                           PC2 = mean(x = PC2),
                                           .groups = "drop"),
                        mapping = ggplot2::aes(label = feature),
                        size = I(x = 3),
                        color = "black")
  )()

tinydenseR::plotTradStats(
  .tdr.obj = lm.cells.DA.0.5,
  .model.name = "default")

stat.test.percentages.DA.0.5 <-
  lm.cells.DA.0.5$density$composition$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DA.0.5$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DA.0.5$results$lm$default$trad$clustering$fit$adj.p[name,"Depletion"],
                p.adj = lm.cells.DA.0.5$results$lm$default$trad$clustering$fit$adj.p[name,"Depletion"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "treatment")

(tinydenseR::plotTradPerc(
  .tdr.obj = lm.cells.DA.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DA.0.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

(tinydenseR::plotDensity(
  .tdr.obj = lm.cells.DA.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density"))

plotBeeswarm(
  .tdr.obj = lm.cells.DA.0.5,
  .model.name = "default",
  .coefs = "Depletion",
  .swarm.title = "Depletion vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE) +
  ggplot2::geom_vline(xintercept = -1,
                      color = "red",
                      linetype = "dashed")

.dea.0.5 <-
  tinydenseR::get.pbDE(
    .tdr.obj = lm.cells.DA.0.5,
    .design = .design.0.5,
    .id = "cluster.4"
  )

tinydenseR::plotPbDE(
  .tdr.obj = lm.cells.DA.0.5,
  .dea.obj = .dea.0.5)

.dea.0.5$adj.p
.dea.0.5$coefficients

lm.cells.DA.0.5 <-
  tinydenseR::get.markerDE(
    .tdr.obj = lm.cells.DA.0.5,
    .id1 = "cluster.4",
    .comparison.name = "cluster4_vs_all"
  )

tinydenseR::plotMarkerDE(
  .tdr.obj = lm.cells.DA.0.5,
  .comparison.name = "cluster4_vs_all",
  .coefs = ".id1")

lm.cells.DA.0.5$results$marker$default$cluster4_vs_all$adj.p
lm.cells.DA.0.5$results$marker$default$cluster4_vs_all$coefficients

# 5%
.setting.meta.5 <-
  sim_data$sample_meta |>
  dplyr::filter(Setting == "5%")

cs.DA.5 <-
  flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(.setting.meta.5$fcs_path,
                            .setting.meta.5$Sample))

flowWorkspace::pData(cs.DA.5)$Sample <-
  flowWorkspace::sampleNames(cs.DA.5)
flowWorkspace::pData(cs.DA.5)$Treatment <-
  .setting.meta.5$Treatment[match(flowWorkspace::sampleNames(cs.DA.5),
                                  .setting.meta.5$Sample)]
flowWorkspace::pData(cs.DA.5)$Batch <-
  .setting.meta.5$Batch[match(flowWorkspace::sampleNames(cs.DA.5),
                              .setting.meta.5$Sample)]

set.seed(seed = 123)
lm.cells.DA.5 <-
  tinydenseR::RunTDR(
    cs.DA.5,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

.meta.DA.5 <-
  lm.cells.DA.5@metadata

.design.5 <-
  model.matrix(object = ~ Treatment + Batch,
               data = .meta.DA.5) |> 
  (\(x)
   `colnames<-`(x = x,
                value = colnames(x = x) |>
                  gsub(pattern = "^Treatment|^Batch",
                       replacement = "",
                       fixed = FALSE))
  )()

# New API: get.lm() returns updated .tdr.obj with results in $results$lm[[.model.name]]
lm.cells.DA.5 <-
  tinydenseR::get.lm(
    .tdr.obj = lm.cells.DA.5,
    .design = .design.5)

lapply(X = names(lm.cells.DA.5@cells),
       FUN = function(s) flowCore::exprs(cs.DA.5[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.5$landmark.embed$pca$HVG]) - lm.cells.DA.5$landmark.embed$pca$center) /
       lm.cells.DA.5$landmark.embed$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.5$landmark.embed$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.5@cells)) |>
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
  )()

final_data_DA |>
  dplyr::filter(Sample %in% names(x = lm.cells.DA.5@cells)) |>
  dplyr::mutate(CellType = factor(x = CellType,
                                  levels = c("target",
                                             "other"))) |>
  tidyr::pivot_longer(
    cols = c(Marker1, Marker2, Marker3, Marker4, Marker5),
    names_to = "Marker",
    values_to = "value"
  ) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = Marker,
                                          y = value,
                                          color = CellType,
                                          fill = CellType)) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "right",
                    axis.text.x = ggplot2::element_text(hjust = 1,
                                                        vjust = 1,
                                                        angle = 30)) +
     ggplot2::labs(title = "ground truth (no expression diff. with depletion)",
                   x = "",
                   y = "expression level") +
     ggplot2::scale_fill_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) + 
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                       alpha = 1)),
                     fill = "none") +
     ggplot2::scale_color_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) + 
     ggplot2::scale_y_log10(breaks = scales::trans_breaks(trans = "log10",
                                                          inv =  function(x) 10^x),
                            labels = scales::trans_format(trans = "log10", 
                                                          format = scales::math_format(10^.x))) + 
     ggplot2::annotation_logticks(sides = "l") +
     scattermore::geom_scattermore(position = ggplot2::position_jitterdodge(jitter.width = 0.5,
                                                                            jitter.height = 0,
                                                                            dodge.width = 0.5,
                                                                            seed = 123),
                                   alpha = 0.1) + 
     ggplot2::geom_violin(color = "black",
                          alpha = 0,
                          position = ggplot2::position_dodge(width = 0.5),
                          quantiles = 0.5,
                          quantile.linetype = "solid") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = names(lm.cells.DA.5@cells),
       FUN = function(s) flowCore::exprs(cs.DA.5[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.5$landmark.embed$pca$HVG]) - lm.cells.DA.5$landmark.embed$pca$center) /
       lm.cells.DA.5$landmark.embed$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.5$landmark.embed$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.5@cells)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.5@cells)) |>
       dplyr::pull(Batch))
  )() |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2)) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "right") +
     ggplot2::labs(title = "ground truth",
                   color = "") +
     scattermore::geom_scattermore(mapping = ggplot2::aes(color = Batch),
                                   pointsize = I(x = 3)) + 
     ggplot2::scale_color_viridis_d() +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

tinydenseR::plotPCA(.tdr.obj = lm.cells.DA.5,
                    .feature = lm.cells.DA.5$metadata$Treatment[lm.cells.DA.5$config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.tdr.obj = lm.cells.DA.5,
                     .feature = lm.cells.DA.5$results$lm$default$fit$coefficients[,"Depletion"],
                     .plot.title = "Depletion vs Baseline",
                     .color.label = "density\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DA.5,
  .feature =
    ifelse(
      test = lm.cells.DA.5$results$lm$default$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = lm.cells.DA.5$results$lm$default$fit$pca.weighted.q[,"Depletion"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("less abundant",
                      "not sig.",
                      "more abundant")),
  .plot.title = "Depletion vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2)   +
    ggplot2::labs(subtitle = "hypothesis testing"))

tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DA.5,
  .feature = lm.cells.DA.5$landmark.annot$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 2) |> 
  (\(x)
   x +
     ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
     ggplot2::geom_text(data = x$data |>
                          dplyr::group_by(feature) |>
                          dplyr::summarize(PC1 = mean(x = PC1),
                                           PC2 = mean(x = PC2),
                                           .groups = "drop"),
                        mapping = ggplot2::aes(label = feature),
                        size = I(x = 3),
                        color = "black")
  )()

tinydenseR::plotTradStats(
  .tdr.obj = lm.cells.DA.5,
  .model.name = "default")

stat.test.percentages.DA.5 <-
  lm.cells.DA.5$density$composition$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DA.5$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DA.5$results$lm$default$trad$clustering$fit$adj.p[name,"Depletion"],
                p.adj = lm.cells.DA.5$results$lm$default$trad$clustering$fit$adj.p[name,"Depletion"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "treatment")


(tinydenseR::plotTradPerc(
  .tdr.obj = lm.cells.DA.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DA.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

(tinydenseR::plotDensity(
  .tdr.obj = lm.cells.DA.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density"))

tinydenseR::plotBeeswarm(
  .tdr.obj = lm.cells.DA.5,
  .model.name = "default",
  .coefs = "Depletion",
  .swarm.title = "Depletion vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE) +
  ggplot2::geom_vline(xintercept = -1,
                      color = "red",
                      linetype = "dashed")


.dea.5 <-
  tinydenseR::get.pbDE(
    .tdr.obj = lm.cells.DA.5,
    .design = .design.5,
    .id = "cluster.3"
  )

tinydenseR::plotPbDE(
  .tdr.obj = lm.cells.DA.5,
  .dea.obj = .dea.5)

.dea.5$adj.p
.dea.5$coefficients

lm.cells.DA.5 <-
  tinydenseR::get.markerDE(
    .tdr.obj = lm.cells.DA.5,
    .id1 = "cluster.3",
    .comparison.name = "cluster3_vs_all"
  )

tinydenseR::plotMarkerDE(
  .tdr.obj = lm.cells.DA.5,
  .comparison.name = "cluster3_vs_all",
  .coefs = ".id1")

lm.cells.DA.5$results$marker$default$cluster3_vs_all$adj.p
lm.cells.DA.5$results$marker$default$cluster3_vs_all$coefficients

# 50%
.setting.meta.50 <-
  sim_data$sample_meta |>
  dplyr::filter(Setting == "50%")

cs.DA.50 <-
  flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(.setting.meta.50$fcs_path,
                            .setting.meta.50$Sample))

flowWorkspace::pData(cs.DA.50)$Sample <-
  flowWorkspace::sampleNames(cs.DA.50)
flowWorkspace::pData(cs.DA.50)$Treatment <-
  .setting.meta.50$Treatment[match(flowWorkspace::sampleNames(cs.DA.50),
                                   .setting.meta.50$Sample)]
flowWorkspace::pData(cs.DA.50)$Batch <-
  .setting.meta.50$Batch[match(flowWorkspace::sampleNames(cs.DA.50),
                               .setting.meta.50$Sample)]

set.seed(seed = 123)
lm.cells.DA.50 <-
  tinydenseR::RunTDR(
    cs.DA.50,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

.meta.DA.50 <-
  lm.cells.DA.50@metadata

.design.50 <-
  model.matrix(object = ~ Treatment + Batch,
               data = .meta.DA.50) |> 
  (\(x)
   `colnames<-`(x = x,
                value = colnames(x = x) |>
                  gsub(pattern = "^Treatment|^Batch",
                       replacement = "",
                       fixed = FALSE))
  )()

lm.cells.DA.50 <-
  tinydenseR::get.lm(
    .tdr.obj = lm.cells.DA.50,
    .design = .design.50)

lapply(X = names(lm.cells.DA.50@cells),
       FUN = function(s) flowCore::exprs(cs.DA.50[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.50$landmark.embed$pca$HVG]) - lm.cells.DA.50$landmark.embed$pca$center) /
       lm.cells.DA.50$landmark.embed$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.50$landmark.embed$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.50@cells)) |>
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
  )() 

final_data_DA |>
  dplyr::filter(Sample %in% names(x = lm.cells.DA.50@cells)) |>
  dplyr::mutate(CellType = factor(x = CellType,
                                  levels = c("target",
                                             "other"))) |>
  tidyr::pivot_longer(
    cols = c(Marker1, Marker2, Marker3, Marker4, Marker5),
    names_to = "Marker",
    values_to = "value"
  ) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = Marker,
                                          y = value,
                                          color = CellType,
                                          fill = CellType)) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "right",
                    axis.text.x = ggplot2::element_text(hjust = 1,
                                                        vjust = 1,
                                                        angle = 30)) +
     ggplot2::labs(title = "ground truth (no expression diff. with depletion)",
                   x = "",
                   y = "expression level") +
     ggplot2::scale_fill_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) + 
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                       alpha = 1)),
                     fill = "none") +
     ggplot2::scale_color_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) + 
     ggplot2::scale_y_log10(breaks = scales::trans_breaks(trans = "log10",
                                                          inv =  function(x) 10^x),
                            labels = scales::trans_format(trans = "log10", 
                                                          format = scales::math_format(10^.x))) + 
     ggplot2::annotation_logticks(sides = "l") +
     scattermore::geom_scattermore(position = ggplot2::position_jitterdodge(jitter.width = 0.5,
                                                                            jitter.height = 0,
                                                                            dodge.width = 0.5,
                                                                            seed = 123),
                                   alpha = 0.1) + 
     ggplot2::geom_violin(color = "black",
                          alpha = 0,
                          position = ggplot2::position_dodge(width = 0.5),
                          quantiles = 0.5, 
                          quantile.linetype = "solid") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )() 

lapply(X = names(lm.cells.DA.50@cells),
       FUN = function(s) flowCore::exprs(cs.DA.50[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.50$landmark.embed$pca$HVG]) - lm.cells.DA.50$landmark.embed$pca$center) /
       lm.cells.DA.50$landmark.embed$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.50$landmark.embed$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.50@cells)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DA |>
       dplyr::filter(Sample %in% names(x = lm.cells.DA.50@cells)) |>
       dplyr::pull(Batch))
  )() |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2)) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "right") +
     ggplot2::labs(title = "ground truth",
                   color = "") +
     scattermore::geom_scattermore(mapping = ggplot2::aes(color = Batch),
                                   pointsize = I(x = 3)) + 
     ggplot2::scale_color_viridis_d() +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )() 

tinydenseR::plotPCA(.tdr.obj = lm.cells.DA.50,
                    .feature = lm.cells.DA.50$metadata$Treatment[lm.cells.DA.50$config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.50,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.tdr.obj = lm.cells.DA.50,
                     .feature = lm.cells.DA.50$results$lm$default$fit$coefficients[,"Depletion"],
                     .plot.title = "Depletion vs Baseline",
                     .color.label = "density\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank())) 

(tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DA.50,
  .feature =
    ifelse(
      test = lm.cells.DA.50$results$lm$default$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = lm.cells.DA.50$results$lm$default$fit$pca.weighted.q[,"Depletion"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("less abundant",
                      "not sig.",
                      "more abundant")),
  .plot.title = "Depletion vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2)   +
    ggplot2::labs(subtitle = "hypothesis testing"))

tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DA.50,
  .feature = lm.cells.DA.50$landmark.annot$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 2) |> 
  (\(x)
   x +
     ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
     ggplot2::geom_text(data = x$data |>
                          dplyr::group_by(feature) |>
                          dplyr::summarize(PC1 = mean(x = PC1),
                                           PC2 = mean(x = PC2),
                                           .groups = "drop"),
                        mapping = ggplot2::aes(label = feature),
                        size = I(x = 3),
                        color = "black")
  )()

tinydenseR::plotTradStats(
  .tdr.obj = lm.cells.DA.50,
  .model.name = "default")

stat.test.percentages.DA.50 <-
  lm.cells.DA.50$density$composition$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DA.50$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DA.50$results$lm$default$trad$clustering$fit$adj.p[name,"Depletion"],
                p.adj = lm.cells.DA.50$results$lm$default$trad$clustering$fit$adj.p[name,"Depletion"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "treatment")


(tinydenseR::plotTradPerc(
  .tdr.obj = lm.cells.DA.50,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DA.50, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.050, 0.150))))

(tinydenseR::plotDensity(
  .tdr.obj = lm.cells.DA.50,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density"))

tinydenseR::plotBeeswarm(
  .tdr.obj = lm.cells.DA.50,
  .model.name = "default",
  .coefs = "Depletion",
  .swarm.title = "Depletion vs Baseline",
  .row.space.scaler = 0.50,
  .perc.plot = FALSE) +
  ggplot2::geom_vline(xintercept = -1,
                      color = "red",
                      linetype = "dashed")

.dea.50 <-
  tinydenseR::get.pbDE(
    .tdr.obj = lm.cells.DA.50,
    .design = .design.50,
    .id = "cluster.2"
  )

tinydenseR::plotPbDE(
  .tdr.obj = lm.cells.DA.50,
  .dea.obj = .dea.50) 

.dea.50$adj.p
.dea.50$coefficients

lm.cells.DA.50 <-
  tinydenseR::get.markerDE(
    .tdr.obj = lm.cells.DA.50,
    .id1 = "cluster.2",
    .comparison.name = "cluster2_vs_all"
  )

tinydenseR::plotMarkerDE(
  .tdr.obj = lm.cells.DA.50,
  .comparison.name = "cluster2_vs_all",
  .coefs = ".id1")

lm.cells.DA.50$results$marker$default$cluster2_vs_all$adj.p
lm.cells.DA.50$results$marker$default$cluster2_vs_all$coefficients

# permutation tests
source(file = "https://raw.githubusercontent.com/Novartis/tinydenseR/refs/heads/main/inst/scripts/perm_utils.R")

# For each condition (0.5%, 5%, 50%), run exact permutation test
conditions <- list(
  "0.5%" = list(lm_obj = lm.cells.DA.0.5, 
                meta = .meta.DA.0.5),
  "5%"   = list(lm_obj = lm.cells.DA.5, 
                meta = .meta.DA.5),
  "50%"   = list(lm_obj = lm.cells.DA.50, 
                 meta = .meta.DA.50)
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
                                      coef_name = "Depletion")
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
plot_n_sig(results = results)
plot_min_q(results = results)

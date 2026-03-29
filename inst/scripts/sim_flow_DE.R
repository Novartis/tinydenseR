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

# Simulate DE data and write FCS files
sim_data <- simulate_DE_data()
final_data_DE <- sim_data$cell_meta |>
  dplyr::mutate(Treatment = factor(x = Treatment,
                                   levels = c("Baseline", "Activation")))

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
          FUN = function(sd_shift) {
            lm(formula = log2(x = prop) ~ Treatment + Batch,
               data = x[x$SD_Shift == sd_shift,]) |> 
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

final_data_DE |>
  dplyr::group_by(Sample, Treatment, Batch, SD_Shift, CellType) |>
  dplyr::summarize(Marker1 = mean(log(Marker1))) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |>
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$SD_Shift) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(sd_shift) {
            lm(formula = Marker1 ~ Treatment + Batch,
               data = x[x$SD_Shift == sd_shift,]) |>
              summary()
          })
  )()

final_data_DE |>
  dplyr::group_by(Sample, Treatment, Batch, SD_Shift, CellType) |>
  dplyr::summarize(Marker2 = mean(log(Marker2))) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |>
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$SD_Shift) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(sd_shift) {
            lm(formula = Marker2 ~ Treatment + Batch,
               data = x[x$SD_Shift == sd_shift,]) |>
              summary()
          })
  )()

final_data_DE |>
  dplyr::group_by(Sample, Treatment, Batch, SD_Shift, CellType) |>
  dplyr::summarize(Marker3 = mean(log(Marker3))) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |>
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$SD_Shift) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(sd_shift) {
            lm(formula = Marker3 ~ Treatment + Batch,
               data = x[x$SD_Shift == sd_shift,]) |>
              summary()
          })
  )()

final_data_DE |>
  dplyr::group_by(Sample, Treatment, Batch, SD_Shift, CellType) |>
  dplyr::summarize(Marker4 = mean(log(Marker4))) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |>
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$SD_Shift) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(sd_shift) {
            lm(formula = Marker4 ~ Treatment + Batch,
               data = x[x$SD_Shift == sd_shift,]) |>
              summary()
          })
  )()

final_data_DE |>
  dplyr::group_by(Sample, Treatment, Batch, SD_Shift, CellType) |>
  dplyr::summarize(Marker5 = mean(log(Marker5))) |> 
  dplyr::ungroup() |>
  dplyr::filter(CellType == "target") |>
  as.data.frame() |>
  (\(x)
   lapply(X = unique(x = x$SD_Shift) |> 
            (\(y)
             setNames(object = y,
                      nm = y)
            )(),
          FUN = function(sd_shift) {
            lm(formula = Marker5 ~ Treatment + Batch,
               data = x[x$SD_Shift == sd_shift,]) |>
              summary()
          })
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
                                    .setting.meta.0.5$Sample)]
flowWorkspace::pData(cs.DE.0.5)$Batch <-
  .setting.meta.0.5$Batch[match(flowWorkspace::sampleNames(cs.DE.0.5),
                                .setting.meta.0.5$Sample)]

set.seed(seed = 123)
lm.cells.DE.0.5 <-
  tinydenseR::RunTDR(
    cs.DE.0.5,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

.meta.DE.0.5 <-
  lm.cells.DE.0.5@metadata

lm.cells.DE.0.5 <-
  tinydenseR::get.map(.tdr.obj = lm.cells.DE.0.5)

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

# New API: get.lm() returns updated .tdr.obj with results in $map$lm[[.model.name]]
lm.cells.DE.0.5 <-
  tinydenseR::get.lm(
    .tdr.obj = lm.cells.DE.0.5,
    .design = .design.0.5)

lapply(X = names(lm.cells.DE.0.5@cells),
       FUN = function(s) flowCore::exprs(cs.DE.0.5[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DE.0.5$pca$HVG]) - lm.cells.DE.0.5$pca$center) /
       lm.cells.DE.0.5$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DE.0.5$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.0.5@cells)) |>
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
     ggplot2::labs(title = "ground truth (no % diff. with activation)") +
     ggplot2::stat_bin_hex(bins = 128) + 
     ggplot2::scale_fill_viridis_c(trans = "log") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

final_data_DE |>
  dplyr::filter(Sample %in% names(x = lm.cells.DE.0.5@cells)) |>
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
     ggplot2::labs(title = "ground truth",
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

lapply(X = names(lm.cells.DE.0.5@cells),
       FUN = function(s) flowCore::exprs(cs.DE.0.5[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DE.0.5$pca$HVG]) - lm.cells.DE.0.5$pca$center) /
       lm.cells.DE.0.5$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DE.0.5$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.0.5@cells)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.0.5@cells)) |>
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

tinydenseR::plotPCA(.tdr.obj = lm.cells.DE.0.5,
                    .feature = lm.cells.DE.0.5$metada$Treatment[lm.cells.DE.0.5$config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.tdr.obj = lm.cells.DE.0.5,
                     .feature = lm.cells.DE.0.5$map$lm$default$fit$coefficients[,"Activation"],
                     .plot.title = "Activation vs Baseline",
                     .color.label = "density\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DE.0.5,
  .feature =
    ifelse(
      test = lm.cells.DE.0.5$map$lm$default$fit$coefficients[,"Activation"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = lm.cells.DE.0.5$map$lm$default$fit$pca.weighted.q[,"Activation"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("less abundant",
                      "not sig.",
                      "more abundant")),
  .plot.title = "Activation vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2)   +
    ggplot2::labs(subtitle = "hypothesis testing"))

tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DE.0.5,
  .feature = lm.cells.DE.0.5$landmark.annot$clustering$ids,
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
    .tdr.obj = lm.cells.DE.0.5,
    .model.name = "default")

stat.test.percentages.DE.0.5 <-
  lm.cells.DE.0.5$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.0.5$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DE.0.5$map$lm$default$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = lm.cells.DE.0.5$map$lm$default$trad$clustering$fit$adj.p[name,"Activation"]) |>
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
  .tdr.obj = lm.cells.DE.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.0.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))) |>
  gridExtra::grid.arrange()

(tinydenseR::plotDensity(
  .tdr.obj = lm.cells.DE.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density"))

tinydenseR::plotBeeswarm(
  .tdr.obj = lm.cells.DE.0.5,
  .model.name = "default",
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE)

.dea.0.5 <-
  tinydenseR::get.pbDE(
    .tdr.obj = lm.cells.DE.0.5,
    .design = .design.0.5,
    .id = "cluster.4"
  )

(tinydenseR::plotPbDE(
  .tdr.obj = lm.cells.DE.0.5,
  .dea.obj = .dea.0.5, 
  .coefs = c("Activation", "Batch2"),
  .col.space.scaler = 0.15) +
    ggplot2::coord_flip() +
    ggplot2::labs(subtitle = "for cluster.4") +
    ggplot2::scale_fill_gradientn(
      colours = c(unname(obj = tinydenseR::Color.Palette[1,6]),
                  unname(obj = tinydenseR::Color.Palette[1,2])),
      limits  = c(0, 0.5),
      breaks = c(0, 0.25, 0.5),
      values  = scales::rescale(x = c(0, 0.5),
                                from = c(0, 0.5))) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                     shape = 21))))

.dea.0.5$adj.p
.dea.0.5$coefficients

lm.cells.DE.0.5 <-
  tinydenseR::get.pbDE(
    x = lm.cells.DE.0.5,
    .mode = "marker",
    .id = "cluster.4",
    .result.name = "cluster4_vs_all"
  )

tinydenseR::plotMarkerDE(
  .tdr.obj = lm.cells.DE.0.5,
  .comparison.name = "cluster4_vs_all",
  .coefs = ".id1") |>
  gridExtra::grid.arrange()

lm.cells.DE.0.5$markerDE$default$cluster4_vs_all$adj.p
lm.cells.DE.0.5$markerDE$default$cluster4_vs_all$coefficients

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
                                    .setting.meta.1$Sample)]
flowWorkspace::pData(cs.DE.1)$Batch <-
  .setting.meta.1$Batch[match(flowWorkspace::sampleNames(cs.DE.1),
                                .setting.meta.1$Sample)]

set.seed(seed = 123)
lm.cells.DE.1 <-
  tinydenseR::RunTDR(
    cs.DE.1,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

.meta.DE.1 <-
  lm.cells.DE.1@metadata

lm.cells.DE.1 <-
  tinydenseR::get.map(.tdr.obj = lm.cells.DE.1)

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

# New API: get.lm() returns updated .tdr.obj with results in $map$lm[[.model.name]]
lm.cells.DE.1 <-
  tinydenseR::get.lm(
    .tdr.obj = lm.cells.DE.1,
    .design = .design.1)

lapply(X = names(lm.cells.DE.1@cells),
       FUN = function(s) flowCore::exprs(cs.DE.1[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DE.1$pca$HVG]) - lm.cells.DE.1$pca$center) /
       lm.cells.DE.1$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DE.1$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.1@cells)) |>
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
     ggplot2::labs(title = "ground truth (no % diff. with activation)") +
     ggplot2::stat_bin_hex(bins = 128) + 
     ggplot2::scale_fill_viridis_c(trans = "log") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

final_data_DE |>
  dplyr::filter(Sample %in% names(x = lm.cells.DE.1@cells)) |>
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
     ggplot2::labs(title = "ground truth",
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

lapply(X = names(lm.cells.DE.1@cells),
       FUN = function(s) flowCore::exprs(cs.DE.1[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DE.1$pca$HVG]) - lm.cells.DE.1$pca$center) /
       lm.cells.DE.1$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DE.1$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.1@cells)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.1@cells)) |>
       dplyr::pull(Batch))
  )() |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2)) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1),
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

tinydenseR::plotPCA(.tdr.obj = lm.cells.DE.1,
                    .feature = lm.cells.DE.1$metada$Treatment[lm.cells.DE.1$config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.tdr.obj = lm.cells.DE.1,
                     .feature = lm.cells.DE.1$map$lm$default$fit$coefficients[,"Activation"],
                     .plot.title = "Activation vs Baseline",
                     .color.label = "density\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DE.1,
  .feature =
    ifelse(
      test = lm.cells.DE.1$map$lm$default$fit$coefficients[,"Activation"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = lm.cells.DE.1$map$lm$default$fit$pca.weighted.q[,"Activation"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("less abundant",
                      "not sig.",
                      "more abundant")),
  .plot.title = "Activation vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2)   +
    ggplot2::labs(subtitle = "hypothesis testing"))

tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DE.1,
  .feature = lm.cells.DE.1$landmark.annot$clustering$ids,
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
    .tdr.obj = lm.cells.DE.1,
    .model.name = "default")

stat.test.percentages.DE.1 <-
  lm.cells.DE.1$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.1$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DE.1$map$lm$default$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = lm.cells.DE.1$map$lm$default$trad$clustering$fit$adj.p[name,"Activation"]) |>
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
  .tdr.obj = lm.cells.DE.1,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.1, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

(tinydenseR::plotDensity(
  .tdr.obj = lm.cells.DE.1,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density"))

tinydenseR::plotBeeswarm(
  .tdr.obj = lm.cells.DE.1,
  .model.name = "default",
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE)

.dea.1 <-
  tinydenseR::get.pbDE(
    .tdr.obj = lm.cells.DE.1,
    .design = .design.1,
    .id = "cluster.3"
  )

(tinydenseR::plotPbDE(
  .tdr.obj = lm.cells.DE.1,
  .dea.obj = .dea.1, 
  .coefs = c("Activation", "Batch2"),
  .col.space.scaler = 0.15) +
    ggplot2::coord_flip() +
    ggplot2::labs(subtitle = "for cluster.3") +
    ggplot2::scale_fill_gradientn(
      colours = c(unname(obj = tinydenseR::Color.Palette[1,6]),
                  unname(obj = tinydenseR::Color.Palette[1,2])),
      limits  = c(0, 1),
      breaks = c(0, 0.5, 1),
      values  = scales::rescale(x = c(0, 1),
                                from = c(0, 1))) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                     shape = 21))))

.dea.1$adj.p
.dea.1$coefficients

lm.cells.DE.1 <-
  tinydenseR::get.pbDE(
    x = lm.cells.DE.1,
    .mode = "marker",
    .id = "cluster.3",
    .result.name = "cluster3_vs_all"
  )

tinydenseR::plotMarkerDE(
  .tdr.obj = lm.cells.DE.1,
  .comparison.name = "cluster3_vs_all",
  .coefs = ".id1")

lm.cells.DE.1$markerDE$default$cluster3_vs_all$adj.p
lm.cells.DE.1$markerDE$default$cluster3_vs_all$coefficients

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
                                    .setting.meta.2$Sample)]
flowWorkspace::pData(cs.DE.2)$Batch <-
  .setting.meta.2$Batch[match(flowWorkspace::sampleNames(cs.DE.2),
                                .setting.meta.2$Sample)]

set.seed(seed = 123)
lm.cells.DE.2 <-
  tinydenseR::RunTDR(
    cs.DE.2,
    .sample.var = "Sample",
    .assay.type = "cyto",
    .markers = paste0("Marker", 1:5),
    .seed = 123,
    .verbose = TRUE,
    .cl.resolution.parameter = 0.5)

.meta.DE.2 <-
  lm.cells.DE.2@metadata

lm.cells.DE.2 <-
  tinydenseR::get.map(.tdr.obj = lm.cells.DE.2)

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

# New API: get.lm() returns updated .tdr.obj with results in $map$lm[[.model.name]]
lm.cells.DE.2 <-
  tinydenseR::get.lm(
    .tdr.obj = lm.cells.DE.2,
    .design = .design.2)

lapply(X = names(lm.cells.DE.2@cells),
       FUN = function(s) flowCore::exprs(cs.DE.2[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DE.2$pca$HVG]) - lm.cells.DE.2$pca$center) /
       lm.cells.DE.2$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DE.2$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.2@cells)) |>
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
     ggplot2::labs(title = "ground truth (no % diff. with activation)") +
     ggplot2::stat_bin_hex(bins = 128) + 
     ggplot2::scale_fill_viridis_c(trans = "log") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

final_data_DE |>
  dplyr::filter(Sample %in% names(x = lm.cells.DE.2@cells)) |>
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
     ggplot2::labs(title = "ground truth",
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

lapply(X = names(lm.cells.DE.2@cells),
       FUN = function(s) flowCore::exprs(cs.DE.2[[s]])) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DE.2$pca$HVG]) - lm.cells.DE.2$pca$center) /
       lm.cells.DE.2$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DE.2$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.2@cells)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DE |>
       dplyr::filter(Sample %in% names(x = lm.cells.DE.2@cells)) |>
       dplyr::pull(Batch))
  )() |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2)) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
     ggplot2::facet_grid(cols = ggplot2::vars(Treatment)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1),
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

tinydenseR::plotPCA(.tdr.obj = lm.cells.DE.2,
                    .feature = lm.cells.DE.2$metada$Treatment[lm.cells.DE.2$config$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.tdr.obj = lm.cells.DE.2,
                     .feature = lm.cells.DE.2$map$lm$default$fit$coefficients[,"Activation"],
                     .plot.title = "Activation vs Baseline",
                     .color.label = "density\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DE.2,
  .feature =
    ifelse(
      test = lm.cells.DE.2$map$lm$default$fit$coefficients[,"Activation"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = lm.cells.DE.2$map$lm$default$fit$pca.weighted.q[,"Activation"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("less abundant",
                      "not sig.",
                      "more abundant")),
  .plot.title = "Activation vs Baseline",
  .color.label = "q < 0.1",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2)   +
    ggplot2::labs(subtitle = "hypothesis testing"))

tinydenseR::plotPCA(
  .tdr.obj = lm.cells.DE.2,
  .feature = lm.cells.DE.2$landmark.annot$clustering$ids,
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
    .tdr.obj = lm.cells.DE.2,
    .model.name = "default")

stat.test.percentages.DE.2 <-
  lm.cells.DE.2$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.2$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = lm.cells.DE.2$map$lm$default$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = lm.cells.DE.2$map$lm$default$trad$clustering$fit$adj.p[name,"Activation"]) |>
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
  .tdr.obj = lm.cells.DE.2,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.2, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))))

(tinydenseR::plotDensity(
  .tdr.obj = lm.cells.DE.2,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster density"))

tinydenseR::plotBeeswarm(
    .tdr.obj = lm.cells.DE.2,
    .model.name = "default",
    .coefs = "Activation",
    .swarm.title = "Activation vs Baseline",
    .row.space.scaler = 0.5,
    .perc.plot = FALSE)

.dea.2 <-
  tinydenseR::get.pbDE(
    .tdr.obj = lm.cells.DE.2,
    .design = .design.2,
    .id = "cluster.4"
  )

(tinydenseR::plotPbDE(
  .tdr.obj = lm.cells.DE.2,
  .dea.obj = .dea.2, 
  .coefs = c("Activation", "Batch2"),
  .col.space.scaler = 0.15) +
    ggplot2::coord_flip() +
    ggplot2::labs(subtitle = "for cluster.4") +
    ggplot2::scale_fill_gradientn(
      colours = c(unname(obj = tinydenseR::Color.Palette[1,6]),
                  unname(obj = tinydenseR::Color.Palette[1,2])),
      limits  = c(0, 2),
      breaks = c(0, 1, 2),
      values  = scales::rescale(x = c(0, 2),
                                from = c(0, 2))) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                     shape = 21))))

.dea.2$adj.p
.dea.2$coefficients

lm.cells.DE.2 <-
  tinydenseR::get.pbDE(
    x = lm.cells.DE.2,
    .mode = "marker",
    .id = "cluster.4",
    .result.name = "cluster4_vs_all"
  )

tinydenseR::plotMarkerDE(
  .tdr.obj = lm.cells.DE.2,
  .comparison.name = "cluster4_vs_all",
  .coefs = ".id1")

lm.cells.DE.2$markerDE$default$cluster4_vs_all$adj.p
lm.cells.DE.2$markerDE$default$cluster4_vs_all$coefficients

# permutation tests
source(file = "https://raw.githubusercontent.com/Novartis/tinydenseR/inst/scripts/perm_utils.R")

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
plot_n_sig(results = results)
plot_min_q(results = results)

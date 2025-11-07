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

wd <- "path/to/your/working/directory/"
rd <- file.path(wd, "res")

setwd(dir = wd)

set.seed(42)

# Parameters
groups <- c("Baseline", "Depletion")
batches <- c("Batch1", "Batch2")
settings <- c(0.005, 0.05, 0.5)  # Proportions for Baseline
samples_per_group <- 6
mean_cells <- 50000
sd_cells <- 500  # Noise in total cell count

# Initialize storage
data_list_DA <- list()

# Simulate data
for (setting in settings) {
  for (group in groups) {
    for (sample_id in 1:samples_per_group) {
      batch <- 
        if(sample_id %% 2 == 0) {
          "Batch2"
        } else {
          "Batch1"
        }
      total_cells <- round(rnorm(1, mean = mean_cells, sd = sd_cells))
      total_cells <- max(total_cells, 1000)  # Ensure a minimum number of cells
      
      # Two-fold difference in proportions
      proportion <- if (group == "Baseline") setting else setting / 2
      
      num_interest <- round(total_cells * proportion)
      num_other <- total_cells - num_interest
      cell_types <- c(rep("target", num_interest), rep("other", num_other))
      cell_types <- sample(cell_types)
      
      # Simulate expression data
      marker1 <- numeric(total_cells)
      marker2 <- rlnorm(total_cells, meanlog = 0, sdlog = 1.5)
      marker3 <- rlnorm(total_cells, meanlog = 0, sdlog = 2.5)
      marker4 <- numeric(total_cells)  
      marker5 <- numeric(total_cells)  
      
      # Assign Marker1, Marker4 Marker5 based on cell type
      marker1[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 2)
      marker1[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + 5, sdlog = 2)  # Shift by 5 SD
      marker4[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.2)
      marker4[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + 3, sdlog = 1.2)  # Shift by 3 SD
      marker5[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.8)
      marker5[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + 7, sdlog = 1.8)  # Shift by 7 SD
      
      # Add batch effect
      if (batch == "Batch2") {
        marker1 <- marker1 * rlnorm(total_cells, meanlog = 0.1, sdlog = 0.5)
        marker2 <- marker2 * rlnorm(total_cells, meanlog = 0.2, sdlog = 0.3)
        marker3 <- marker3 * rlnorm(total_cells, meanlog = 0.3, sdlog = 0.4)
        marker4 <- marker4 * rlnorm(total_cells, meanlog = 0.4, sdlog = 0.3)
        marker5 <- marker5 * rlnorm(total_cells, meanlog = 0.5, sdlog = 0.5)
      }
      
      sample_name <- paste0(group, "_S", sample_id, "_Set", setting * 100)
      df <- data.frame(
        Sample = sample_name,
        Treatment = group,
        Batch = batch,
        Setting = paste0(setting * 100, "%"),
        CellType = cell_types,
        Marker1 = marker1,
        Marker2 = marker2,
        Marker3 = marker3,
        Marker4 = marker4,
        Marker5 = marker5
      )
      data_list_DA[[length(data_list_DA) + 1]] <- df
    }
  }
}

# Combine
final_data_DA <- do.call(rbind, data_list_DA)

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
.cells.DA.0.5 <-
  final_data_DA |>
  dplyr::filter(Setting == "0.5%") |>
  dplyr::pull(Sample) |>
  unique() |> 
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  lapply(FUN = function(sample.id){
    
    uri <- tempfile(fileext = ".RDS")
    
    saveRDS(object = final_data_DA[final_data_DA$Sample == sample.id,
                                   c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")] |>
              as.matrix() |>
              log(),
            file = uri,
            compress = FALSE)
    
    return(uri)
    
  })

.meta.DA.0.5 <-
  final_data_DA |>
  dplyr::filter(Sample %in% names(x = .cells.DA.0.5)) |>
  dplyr::select(Sample, Treatment, Batch) |>
  dplyr::distinct() |> 
  (\(x)
   `rownames<-`(x = x[match(x = names(x = .cells.DA.0.5),
                            table = x$Sample),
                      c("Treatment", "Batch")],
                value = names(x = .cells.DA.0.5))
  )()

set.seed(seed = 123)
lm.cells.DA.0.5 <-
  tinydenseR::setup.lm.obj(
    .cells = .cells.DA.0.5,
    .meta = .meta.DA.0.5,
    .assay.type = "cyto") |>
  tinydenseR::get.landmarks() |>
  tinydenseR::get.graph(.cl.resolution.parameter = 0.5)

lm.cells.DA.0.5 <-
  tinydenseR::get.map(.lm.obj = lm.cells.DA.0.5)

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

Treatment.stats.0.5 <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells.DA.0.5,
    .design = .design.0.5)

lapply(X = .cells.DA.0.5,
       FUN = readRDS) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.0.5$pca$HVG]) - lm.cells.DA.0.5$pca$center) /
       lm.cells.DA.0.5$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.0.5$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.0.5)) |>
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
  dplyr::filter(Sample %in% names(x = .cells.DA.0.5)) |>
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
                          draw_quantiles = 0.5) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = .cells.DA.0.5,
       FUN = readRDS) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.0.5$pca$HVG]) - lm.cells.DA.0.5$pca$center) /
       lm.cells.DA.0.5$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.0.5$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.0.5)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.0.5)) |>
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

tinydenseR::plotPCA(.lm.obj = lm.cells.DA.0.5,
                    .feature = lm.cells.DA.0.5$metada$Treatment[lm.cells.DA.0.5$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.lm.obj = lm.cells.DA.0.5,
                     .feature = Treatment.stats.0.5$fit$coefficients[,"Depletion"],
                     .plot.title = "Depletion vs Baseline",
                     .color.label = "abundance\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DA.0.5,
  .feature =
    ifelse(
      test = Treatment.stats.0.5$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.0.5$fit$density.weighted.bh.fdr[,"Depletion"] < 0.1,
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

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DA.0.5,
  .feature =
    ifelse(
      test = Treatment.stats.0.5$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.0.5$fit$pca.weighted.q[,"Depletion"] < 0.1,
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
  .lm.obj = lm.cells.DA.0.5,
  .feature = lm.cells.DA.0.5$graph$clustering$ids,
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
  .lm.obj = lm.cells.DA.0.5,
  .stats.obj = Treatment.stats.0.5)

stat.test.percentages.DA.0.5 <-
  lm.cells.DA.0.5$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DA.0.5$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = Treatment.stats.0.5$trad$clustering$fit$adj.p[name,"Depletion"],
                p.adj = Treatment.stats.0.5$trad$clustering$fit$adj.p[name,"Depletion"]) |>
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
  .lm.obj = lm.cells.DA.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DA.0.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells.DA.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster abundance"))

plotBeeswarm(
  .lm.obj = lm.cells.DA.0.5,
  .stats.obj = Treatment.stats.0.5,
  .coefs = "Depletion",
  .swarm.title = "Depletion vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE) +
  ggplot2::geom_vline(xintercept = -1,
                      color = "red",
                      linetype = "dashed")

.dea.0.5 <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells.DA.0.5,
    .design = .design.0.5,
    .id = "cluster.4"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DA.0.5,
  .dea.obj = .dea.0.5)

.dea.0.5$adj.p
.dea.0.5$coefficients

.subset.dea.0.5 <-
  tinydenseR::get.marker(
    .lm.obj = lm.cells.DA.0.5,
    .id1 = "cluster.4"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DA.0.5,
  .dea.obj = .subset.dea.0.5,
  .coefs = ".id1")

.subset.dea.0.5$adj.p
.subset.dea.0.5$coefficients

# 5%
.cells.DA.5 <-
  final_data_DA |>
  dplyr::filter(Setting == "5%") |>
  dplyr::pull(Sample) |>
  unique() |> 
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  lapply(FUN = function(sample.id){
    
    uri <- tempfile(fileext = ".RDS")
    
    saveRDS(object = final_data_DA[final_data_DA$Sample == sample.id,
                                   c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")] |>
              as.matrix() |>
              log(),
            file = uri,
            compress = FALSE)
    
    return(uri)
    
  })

.meta.DA.5 <-
  final_data_DA |>
  dplyr::filter(Sample %in% names(x = .cells.DA.5)) |>
  dplyr::select(Sample, Treatment, Batch) |>
  dplyr::distinct() |> 
  (\(x)
   `rownames<-`(x = x[match(x = names(x = .cells.DA.5),
                            table = x$Sample),
                      c("Treatment", "Batch")],
                value = names(x = .cells.DA.5))
  )()

set.seed(seed = 123)
lm.cells.DA.5 <-
  tinydenseR::setup.lm.obj(
    .cells = .cells.DA.5,
    .meta = .meta.DA.5,
    .assay.type = "cyto") |>
  tinydenseR::get.landmarks() |>
  tinydenseR::get.graph(.cl.resolution.parameter = 0.5)

lm.cells.DA.5 <-
  tinydenseR::get.map(.lm.obj = lm.cells.DA.5)

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

Treatment.stats.5 <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells.DA.5,
    .design = .design.5)

lapply(X = .cells.DA.5,
       FUN = readRDS) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.5$pca$HVG]) - lm.cells.DA.5$pca$center) /
       lm.cells.DA.5$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.5$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.5)) |>
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
  dplyr::filter(Sample %in% names(x = .cells.DA.5)) |>
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
                          draw_quantiles = 0.5) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = .cells.DA.5,
       FUN = readRDS) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.5$pca$HVG]) - lm.cells.DA.5$pca$center) /
       lm.cells.DA.5$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.5$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.5)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.5)) |>
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

tinydenseR::plotPCA(.lm.obj = lm.cells.DA.5,
                    .feature = lm.cells.DA.5$metada$Treatment[lm.cells.DA.5$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.lm.obj = lm.cells.DA.5,
                     .feature = Treatment.stats.5$fit$coefficients[,"Depletion"],
                     .plot.title = "Depletion vs Baseline",
                     .color.label = "abundance\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DA.5,
  .feature =
    ifelse(
      test = Treatment.stats.5$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.5$fit$density.weighted.bh.fdr[,"Depletion"] < 0.1,
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

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DA.5,
  .feature =
    ifelse(
      test = Treatment.stats.5$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.5$fit$pca.weighted.q[,"Depletion"] < 0.1,
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
  .lm.obj = lm.cells.DA.5,
  .feature = lm.cells.DA.5$graph$clustering$ids,
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
  .lm.obj = lm.cells.DA.5,
  .stats.obj = Treatment.stats.5)

stat.test.percentages.DA.5 <-
  lm.cells.DA.5$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DA.5$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = Treatment.stats.5$trad$clustering$fit$adj.p[name,"Depletion"],
                p.adj = Treatment.stats.5$trad$clustering$fit$adj.p[name,"Depletion"]) |>
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
  .lm.obj = lm.cells.DA.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DA.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells.DA.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster abundance"))

tinydenseR::plotBeeswarm(
  .lm.obj = lm.cells.DA.5,
  .stats.obj = Treatment.stats.5,
  .coefs = "Depletion",
  .swarm.title = "Depletion vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE) +
  ggplot2::geom_vline(xintercept = -1,
                      color = "red",
                      linetype = "dashed")


.dea.5 <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells.DA.5,
    .design = .design.5,
    .id = "cluster.3"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DA.5,
  .dea.obj = .dea.5)

.dea.5$adj.p
.dea.5$coefficients

.subset.dea.5 <-
  tinydenseR::get.marker(
    .lm.obj = lm.cells.DA.5,
    .id1 = "cluster.3"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DA.5,
  .dea.obj = .subset.dea.5,
  .coefs = ".id1")

.subset.dea.5$adj.p
.subset.dea.5$coefficients

# 50%
.cells.DA.50 <-
  final_data_DA |>
  dplyr::filter(Setting == "50%") |>
  dplyr::pull(Sample) |>
  unique() |> 
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  lapply(FUN = function(sample.id){
    
    uri <- tempfile(fileext = ".RDS")
    
    saveRDS(object = final_data_DA[final_data_DA$Sample == sample.id,
                                   c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")] |>
              as.matrix() |>
              log(),
            file = uri,
            compress = FALSE)
    
    return(uri)
    
  })

.meta.DA.50 <-
  final_data_DA |>
  dplyr::filter(Sample %in% names(x = .cells.DA.50)) |>
  dplyr::select(Sample, Treatment, Batch) |>
  dplyr::distinct() |> 
  (\(x)
   `rownames<-`(x = x[match(x = names(x = .cells.DA.50),
                            table = x$Sample),
                      c("Treatment", "Batch")],
                value = names(x = .cells.DA.50))
  )()

set.seed(seed = 123)
lm.cells.DA.50 <-
  tinydenseR::setup.lm.obj(
    .cells = .cells.DA.50,
    .meta = .meta.DA.50,
    .assay.type = "cyto") |>
  tinydenseR::get.landmarks() |>
  tinydenseR::get.graph(.cl.resolution.parameter = 0.5)

lm.cells.DA.50 <-
  tinydenseR::get.map(.lm.obj = lm.cells.DA.50)

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

Treatment.stats.50 <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells.DA.50,
    .design = .design.50)

lapply(X = .cells.DA.50,
       FUN = readRDS) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.50$pca$HVG]) - lm.cells.DA.50$pca$center) /
       lm.cells.DA.50$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.50$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.50)) |>
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
  dplyr::filter(Sample %in% names(x = .cells.DA.50)) |>
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
                          draw_quantiles = 0.5) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )() 

lapply(X = .cells.DA.50,
       FUN = readRDS) |>
  do.call(what = rbind) |>
  (\(x)
   (((Matrix::t(x = x[,lm.cells.DA.50$pca$HVG]) - lm.cells.DA.50$pca$center) /
       lm.cells.DA.50$pca$scale) |>
       Matrix::t()) %*%
     lm.cells.DA.50$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |>
  (\(x)
   dplyr::mutate(
     .data = x,
     Treatment = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.50)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DA |>
       dplyr::filter(Sample %in% names(x = .cells.DA.50)) |>
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

tinydenseR::plotPCA(.lm.obj = lm.cells.DA.50,
                    .feature = lm.cells.DA.50$metada$Treatment[lm.cells.DA.50$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.50,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.lm.obj = lm.cells.DA.50,
                     .feature = Treatment.stats.50$fit$coefficients[,"Depletion"],
                     .plot.title = "Depletion vs Baseline",
                     .color.label = "abundance\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank())) 

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DA.50,
  .feature =
    ifelse(
      test = Treatment.stats.50$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.50$fit$density.weighted.bh.fdr[,"Depletion"] < 0.1,
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

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DA.50,
  .feature =
    ifelse(
      test = Treatment.stats.50$fit$coefficients[,"Depletion"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.50$fit$pca.weighted.q[,"Depletion"] < 0.1,
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
  .lm.obj = lm.cells.DA.50,
  .feature = lm.cells.DA.50$graph$clustering$ids,
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
  .lm.obj = lm.cells.DA.50,
  .stats.obj = Treatment.stats.50)

stat.test.percentages.DA.50 <-
  lm.cells.DA.50$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DA.50$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = Treatment.stats.50$trad$clustering$fit$adj.p[name,"Depletion"],
                p.adj = Treatment.stats.50$trad$clustering$fit$adj.p[name,"Depletion"]) |>
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
  .lm.obj = lm.cells.DA.50,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DA.50, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.050, 0.150))))

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells.DA.50,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster abundance"))

tinydenseR::plotBeeswarm(
  .lm.obj = lm.cells.DA.50,
  .stats.obj = Treatment.stats.50,
  .coefs = "Depletion",
  .swarm.title = "Depletion vs Baseline",
  .row.space.scaler = 0.50,
  .perc.plot = FALSE) +
  ggplot2::geom_vline(xintercept = -1,
                      color = "red",
                      linetype = "dashed")

.dea.50 <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells.DA.50,
    .design = .design.50,
    .id = "cluster.2"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DA.50,
  .dea.obj = .dea.50) 

.dea.50$adj.p
.dea.50$coefficients

.subset.dea.50 <-
  tinydenseR::get.marker(
    .lm.obj = lm.cells.DA.50,
    .id1 = "cluster.2"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DA.50,
  .dea.obj = .subset.dea.50,
  .coefs = ".id1")

.subset.dea.50$adj.p
.subset.dea.50$coefficients


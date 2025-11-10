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

set.seed(42)

# Parameters
groups <- c("Baseline", "Activation")
batches <- c("Batch1", "Batch2")
sd_shifts <- c(0.5, 1, 2)  # SD differences in Marker2 for Baseline
samples_per_group <- 6
mean_cells <- 50000
sd_cells <- 500  # Noise in total cell count

# Initialize storage
data_list_DE <- list()

# Simulate data
for (sd_shift in sd_shifts) {
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
      
      # Fixed proportion of cell type of interest
      proportion <- 0.05
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
      
      # Marker2: differential expression between groups in activated only
      marker2[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.5)
      if (group == "Baseline") {
        marker2[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0, sdlog = 1.5)
      } else {
        marker2[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + sd_shift, sdlog = 1.5)
      }
      
      # Add batch effect
      if (batch == "Batch2") {
        marker1 <- marker1 * rlnorm(total_cells, meanlog = 0.1, sdlog = 0.5)
        marker2 <- marker2 * rlnorm(total_cells, meanlog = 0.2, sdlog = 0.3)
        marker3 <- marker3 * rlnorm(total_cells, meanlog = 0.3, sdlog = 0.4)
        marker4 <- marker4 * rlnorm(total_cells, meanlog = 0.4, sdlog = 0.3)
        marker5 <- marker5 * rlnorm(total_cells, meanlog = 0.5, sdlog = 0.5)
      }
      
      sample_name <- paste0(group, "_S", sample_id, "_Shift", sd_shift, "SD")
      df <- data.frame(
        Sample = sample_name,
        Treatment = group,
        Batch = batch,
        SD_Shift = paste0(sd_shift, "SD"),
        CellType = cell_types,
        Marker1 = marker1,
        Marker2 = marker2,
        Marker3 = marker3,
        Marker4 = marker4,
        Marker5 = marker5
      )
      data_list_DE[[length(data_list_DE) + 1]] <- df
    }
  }
}

# Combine
final_data_DE <- do.call(rbind, data_list_DE) |>
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
.cells.DE.0.5 <-
  final_data_DE |>
  dplyr::filter(SD_Shift == "0.5SD") |>
  dplyr::pull(Sample) |>
  unique() |> 
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  lapply(FUN = function(sample.id){
    
    uri <- tempfile(fileext = ".RDS")
    
    saveRDS(object = final_data_DE[final_data_DE$Sample == sample.id,
                                   c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")] |>
              as.matrix() |>
              log(),
            file = uri,
            compress = FALSE)
    
    return(uri)
    
  })

.meta.DE.0.5 <-
  final_data_DE |>
  dplyr::filter(Sample %in% names(x = .cells.DE.0.5)) |>
  dplyr::select(Sample, Treatment, Batch) |>
  dplyr::distinct() |> 
  (\(x)
   `rownames<-`(x = x[match(x = names(x = .cells.DE.0.5),
                            table = x$Sample),
                      c("Treatment", "Batch")],
                value = names(x = .cells.DE.0.5))
  )()

set.seed(seed = 123)
lm.cells.DE.0.5 <-
  tinydenseR::setup.lm.obj(
    .cells = .cells.DE.0.5,
    .meta = .meta.DE.0.5,
    .assay.type = "cyto") |>
  tinydenseR::get.landmarks() |>
  tinydenseR::get.graph(.cl.resolution.parameter = 0.5)

lm.cells.DE.0.5 <-
  tinydenseR::get.map(.lm.obj = lm.cells.DE.0.5)

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

Treatment.stats.0.5 <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells.DE.0.5,
    .design = .design.0.5)

lapply(X = .cells.DE.0.5,
       FUN = readRDS) |>
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
       dplyr::filter(Sample %in% names(x = .cells.DE.0.5)) |>
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
  dplyr::filter(Sample %in% names(x = .cells.DE.0.5)) |>
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
                          draw_quantiles = 0.5) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = .cells.DE.0.5,
       FUN = readRDS) |>
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
       dplyr::filter(Sample %in% names(x = .cells.DE.0.5)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DE |>
       dplyr::filter(Sample %in% names(x = .cells.DE.0.5)) |>
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

tinydenseR::plotPCA(.lm.obj = lm.cells.DE.0.5,
                    .feature = lm.cells.DE.0.5$metada$Treatment[lm.cells.DE.0.5$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.lm.obj = lm.cells.DE.0.5,
                     .feature = Treatment.stats.0.5$fit$coefficients[,"Activation"],
                     .plot.title = "Activation vs Baseline",
                     .color.label = "abundance\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DE.0.5,
  .feature =
    ifelse(
      test = Treatment.stats.0.5$fit$coefficients[,"Activation"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.0.5$fit$pca.weighted.q[,"Activation"] < 0.1,
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
  .lm.obj = lm.cells.DE.0.5,
  .feature = lm.cells.DE.0.5$graph$clustering$ids,
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
    .lm.obj = lm.cells.DE.0.5,
    .stats.obj = Treatment.stats.0.5)

stat.test.percentages.DE.0.5 <-
  lm.cells.DE.0.5$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.0.5$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = Treatment.stats.0.5$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = Treatment.stats.0.5$trad$clustering$fit$adj.p[name,"Activation"]) |>
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
  .lm.obj = lm.cells.DE.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.0.5, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))) |>
  gridExtra::grid.arrange()

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells.DE.0.5,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster abundance"))

tinydenseR::plotBeeswarm(
  .lm.obj = lm.cells.DE.0.5,
  .stats.obj = Treatment.stats.0.5,
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE)

.dea.0.5 <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells.DE.0.5,
    .design = .design.0.5,
    .id = "cluster.4"
  )

(tinydenseR::plotDEA(
  .lm.obj = lm.cells.DE.0.5,
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

.subset.dea.0.5 <-
  tinydenseR::get.marker(
    .lm.obj = lm.cells.DE.0.5,
    .id1 = "cluster.4"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DE.0.5,
  .dea.obj = .subset.dea.0.5, 
  .coefs = ".id1") |>
  gridExtra::grid.arrange()

.subset.dea.0.5$adj.p
.subset.dea.0.5$coefficients

# 1SD
.cells.DE.1 <-
  final_data_DE |>
  dplyr::filter(SD_Shift == "1SD") |>
  dplyr::pull(Sample) |>
  unique() |> 
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  lapply(FUN = function(sample.id){
    
    uri <- tempfile(fileext = ".RDS")
    
    saveRDS(object = final_data_DE[final_data_DE$Sample == sample.id,
                                   c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")] |>
              as.matrix() |>
              log(),
            file = uri,
            compress = FALSE)
    
    return(uri)
    
  })

.meta.DE.1 <-
  final_data_DE |>
  dplyr::filter(Sample %in% names(x = .cells.DE.1)) |>
  dplyr::select(Sample, Treatment, Batch) |>
  dplyr::distinct() |> 
  (\(x)
   `rownames<-`(x = x[match(x = names(x = .cells.DE.1),
                            table = x$Sample),
                      c("Treatment", "Batch")],
                value = names(x = .cells.DE.1))
  )()

set.seed(seed = 123)
lm.cells.DE.1 <-
  tinydenseR::setup.lm.obj(
    .cells = .cells.DE.1,
    .meta = .meta.DE.1,
    .assay.type = "cyto") |>
  tinydenseR::get.landmarks() |>
  tinydenseR::get.graph(.cl.resolution.parameter = 0.5)

lm.cells.DE.1 <-
  tinydenseR::get.map(.lm.obj = lm.cells.DE.1)

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

Treatment.stats.1 <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells.DE.1,
    .design = .design.1)

lapply(X = .cells.DE.1,
       FUN = readRDS) |>
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
       dplyr::filter(Sample %in% names(x = .cells.DE.1)) |>
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
  dplyr::filter(Sample %in% names(x = .cells.DE.1)) |>
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
                          draw_quantiles = 0.5) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = .cells.DE.1,
       FUN = readRDS) |>
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
       dplyr::filter(Sample %in% names(x = .cells.DE.1)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DE |>
       dplyr::filter(Sample %in% names(x = .cells.DE.1)) |>
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

tinydenseR::plotPCA(.lm.obj = lm.cells.DE.1,
                    .feature = lm.cells.DE.1$metada$Treatment[lm.cells.DE.1$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.lm.obj = lm.cells.DE.1,
                     .feature = Treatment.stats.1$fit$coefficients[,"Activation"],
                     .plot.title = "Activation vs Baseline",
                     .color.label = "abundance\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DE.1,
  .feature =
    ifelse(
      test = Treatment.stats.1$fit$coefficients[,"Activation"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.1$fit$pca.weighted.q[,"Activation"] < 0.1,
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
  .lm.obj = lm.cells.DE.1,
  .feature = lm.cells.DE.1$graph$clustering$ids,
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
    .lm.obj = lm.cells.DE.1,
    .stats.obj = Treatment.stats.1)

stat.test.percentages.DE.1 <-
  lm.cells.DE.1$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.1$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = Treatment.stats.1$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = Treatment.stats.1$trad$clustering$fit$adj.p[name,"Activation"]) |>
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
  .lm.obj = lm.cells.DE.1,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.1, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells.DE.1,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster abundance"))

tinydenseR::plotBeeswarm(
  .lm.obj = lm.cells.DE.1,
  .stats.obj = Treatment.stats.1,
  .coefs = "Activation",
  .swarm.title = "Activation vs Baseline",
  .row.space.scaler = 0.5,
  .perc.plot = FALSE)

.dea.1 <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells.DE.1,
    .design = .design.1,
    .id = "cluster.3"
  )

(tinydenseR::plotDEA(
  .lm.obj = lm.cells.DE.1,
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

.subset.dea.1 <-
  tinydenseR::get.marker(
    .lm.obj = lm.cells.DE.1,
    .id1 = "cluster.3"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DE.1,
  .dea.obj = .subset.dea.1, 
  .coefs = ".id1")

.subset.dea.1$adj.p
.subset.dea.1$coefficients

# 2SD
.cells.DE.2 <-
  final_data_DE |>
  dplyr::filter(SD_Shift == "2SD") |>
  dplyr::pull(Sample) |>
  unique() |> 
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  lapply(FUN = function(sample.id){
    
    uri <- tempfile(fileext = ".RDS")
    
    saveRDS(object = final_data_DE[final_data_DE$Sample == sample.id,
                                   c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")] |>
              as.matrix() |>
              log(),
            file = uri,
            compress = FALSE)
    
    return(uri)
    
  })

.meta.DE.2 <-
  final_data_DE |>
  dplyr::filter(Sample %in% names(x = .cells.DE.2)) |>
  dplyr::select(Sample, Treatment, Batch) |>
  dplyr::distinct() |> 
  (\(x)
   `rownames<-`(x = x[match(x = names(x = .cells.DE.2),
                            table = x$Sample),
                      c("Treatment", "Batch")],
                value = names(x = .cells.DE.2))
  )()

set.seed(seed = 123)
lm.cells.DE.2 <-
  tinydenseR::setup.lm.obj(
    .cells = .cells.DE.2,
    .meta = .meta.DE.2,
    .assay.type = "cyto") |>
  tinydenseR::get.landmarks() |>
  tinydenseR::get.graph(.cl.resolution.parameter = 0.5)

lm.cells.DE.2 <-
  tinydenseR::get.map(.lm.obj = lm.cells.DE.2)

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

Treatment.stats.2 <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells.DE.2,
    .design = .design.2)

lapply(X = .cells.DE.2,
       FUN = readRDS) |>
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
       dplyr::filter(Sample %in% names(x = .cells.DE.2)) |>
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
  dplyr::filter(Sample %in% names(x = .cells.DE.2)) |>
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
                          draw_quantiles = 0.5) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

lapply(X = .cells.DE.2,
       FUN = readRDS) |>
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
       dplyr::filter(Sample %in% names(x = .cells.DE.2)) |>
       dplyr::pull(Treatment),
     Batch = final_data_DE |>
       dplyr::filter(Sample %in% names(x = .cells.DE.2)) |>
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

tinydenseR::plotPCA(.lm.obj = lm.cells.DE.2,
                    .feature = lm.cells.DE.2$metada$Treatment[lm.cells.DE.2$key],
                    .cat.feature.color = tinydenseR::Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Treatment")

(tinydenseR::plotPCA(.lm.obj = lm.cells.DE.2,
                     .feature = Treatment.stats.2$fit$coefficients[,"Activation"],
                     .plot.title = "Activation vs Baseline",
                     .color.label = "abundance\nlog2(+0.5)FC",
                     .panel.size = 2,
                     .point.size = 1,
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()))

(tinydenseR::plotPCA(
  .lm.obj = lm.cells.DE.2,
  .feature =
    ifelse(
      test = Treatment.stats.2$fit$coefficients[,"Activation"] < 0,
      yes = "less abundant",
      no = "more abundant") |>
    ifelse(
      test = Treatment.stats.2$fit$pca.weighted.q[,"Activation"] < 0.1,
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
  .lm.obj = lm.cells.DE.2,
  .feature = lm.cells.DE.2$graph$clustering$ids,
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
    .lm.obj = lm.cells.DE.2,
    .stats.obj = Treatment.stats.2)

stat.test.percentages.DE.2 <-
  lm.cells.DE.2$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(treatment = lm.cells.DE.2$metadata$Treatment) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ treatment) |>
  dplyr::mutate(p = Treatment.stats.2$trad$clustering$fit$adj.p[name,"Activation"],
                p.adj = Treatment.stats.2$trad$clustering$fit$adj.p[name,"Activation"]) |>
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
  .lm.obj = lm.cells.DE.2,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages.DE.2, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))))

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells.DE.2,
  .x.split = "Treatment",
  .x.space.scaler = 0.3
) + 
    ggplot2::labs(title = "within-cluster abundance"))

tinydenseR::plotBeeswarm(
    .lm.obj = lm.cells.DE.2,
    .stats.obj = Treatment.stats.2,
    .coefs = "Activation",
    .swarm.title = "Activation vs Baseline",
    .row.space.scaler = 0.5,
    .perc.plot = FALSE)

.dea.2 <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells.DE.2,
    .design = .design.2,
    .id = "cluster.4"
  )

(tinydenseR::plotDEA(
  .lm.obj = lm.cells.DE.2,
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

.subset.dea.2 <-
  tinydenseR::get.marker(
    .lm.obj = lm.cells.DE.2,
    .id1 = "cluster.4"
  )

tinydenseR::plotDEA(
  .lm.obj = lm.cells.DE.2,
  .dea.obj = .subset.dea.2, 
  .coefs = ".id1")

.subset.dea.2$adj.p
.subset.dea.2$coefficients

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
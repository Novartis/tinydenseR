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

# Fetch trajectory data from miloR repository
trajectory_data <- fetch_trajectory_data()

# Extract components
sim_trajectory.meta <- trajectory_data$meta
sim_trajectory <- trajectory_data$SCE

## Process metadata for tinydenseR format
#sim_trajectory.meta <- sim_trajectory.meta[, c("Condition", "Replicate", "Sample")] |>
#  unique()
#rownames(sim_trajectory.meta) <- sim_trajectory.meta$Sample

# Create .meta object containing sample-level data
.meta <- get.meta(.obj = sim_trajectory,
                  .sample.var = "Sample",
                  .verbose = TRUE)

# Create .cells object using SCE method
.cells <- get.cells(.exprs = sim_trajectory,
                    .meta = .meta,
                    .sample.var = "Sample")[rownames(x = .meta)]

set.seed(seed = 123)
lm.cells <-
  tinydenseR::setup.lm.obj(
    .cells = .cells,
    .meta = .meta,
    .assay.type = "RNA",
    .prop.landmarks = 0.15) |>
  tinydenseR::get.landmarks(.nHVG = 500,
                            .nPC = 3) |>
  tinydenseR::get.graph(.cl.resolution.parameter = 2e2, 
                        .k = 10)

lm.cells <-
  tinydenseR::get.map(.lm.obj = lm.cells)

.design <-
  model.matrix(object = ~ Condition,
               data = .meta)

condition.stats <-
  tinydenseR::get.stats(
    .lm.obj = lm.cells,
    .design = .design)

# Extract count matrices from SCE object for plotting
count_matrices <- list()
for(sample in sort(x = .meta$Sample)) {
  sample_cells <- SummarizedExperiment::colData(sim_trajectory)$Sample == sample
  count_matrices[[sample]] <- SingleCellExperiment::counts(sim_trajectory)[, sample_cells]
}

count_matrices |>
  do.call(what = cbind) |>
  Matrix::t() |>
  (\(x)
   log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
  )() |> 
  (\(x)
   (((Matrix::t(x = x[,lm.cells$pca$HVG]) - lm.cells$pca$center) /
       lm.cells$pca$scale) |>
       Matrix::t()) %*%
     lm.cells$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |> 
  (\(x)
   dplyr::mutate(.data = x,
                 cell_id = rownames(x = x)) 
  )() |>
  dplyr::left_join(y = sim_trajectory.meta[,c("cell_id", "Condition", "Sample")],
                   by = "cell_id") |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2,
                                          color = Condition)) +
     ggplot2::facet_grid(cols = ggplot2::vars(Condition)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    legend.position = "none") +
     ggplot2::labs(title = "ground truth") +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 5)))) +
     ggplot2::scale_color_manual(values = unname(obj = Color.Palette[1,c(1,2)])) +
     ggplot2::xlim(c(-21,21)) +
     ggplot2::ylim(c(-7,20)) +
     ggplot2::geom_point(size = I(x = 1),
                         color = "black",
                         alpha = 0.1) +
     ggplot2::geom_density2d(adjust = 1/3,
                             contour_var = "ndensity") +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

count_matrices |>
  do.call(what = cbind) |>
  Matrix::t() |>
  (\(x)
   log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
  )() |> 
  (\(x)
   (((Matrix::t(x = x[,lm.cells$pca$HVG]) - lm.cells$pca$center) /
       lm.cells$pca$scale) |>
       Matrix::t()) %*%
     lm.cells$pca$rotation
  )() |>
  as.matrix() |>
  as.data.frame() |> 
  (\(x)
   dplyr::mutate(.data = x,
                 cell_id = rownames(x = x)) 
  )() |>
  dplyr::left_join(y = sim_trajectory.meta[,c("cell_id", "Condition", "Sample")],
                   by = "cell_id") |>
  dplyr::rename(Group = Condition) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = PC1,
                                          y = PC2,
                                          color = Group)) +
     ggplot2::facet_wrap(facets = ggplot2::vars(Sample), 
                         nrow = 2) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
     ggplot2::labs(title = "samples") +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 5)))) +
     ggplot2::scale_color_manual(values = unname(obj = Color.Palette[1,c(1,2)])) +
     ggplot2::xlim(c(-21,21)) +
     ggplot2::ylim(c(-7,20)) +
     ggplot2::geom_point(size = I(x = 1)) +
     ggh4x::force_panelsizes(cols = grid::unit(x = 1,
                                               units = "in"),
                             rows = grid::unit(x = 1,
                                               units = "in"))
  )()


(tinydenseR::plotPCA(.lm.obj = lm.cells,
                     .feature = "black",
                     .cat.feature.color = "black",
                     .plot.title = "landmarks",
                     .panel.size = 2,
                     .point.size = 1,
                     .color.label = "Group") +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::xlim(c(-21,21)) +
    ggplot2::ylim(c(-7,20)))

tinydenseR::plotPCA(.lm.obj = lm.cells,
                    .feature = lm.cells$metadata$Condition[lm.cells$key],
                    .cat.feature.color = Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Condition")

(tinydenseR::plotPCA(.lm.obj = lm.cells,
                     .feature = condition.stats$fit$coefficients[,"ConditionB"],
                     .panel.size = 2,
                     .point.size = 1,
                     .plot.title = "estimated abundance\nlog2(+0.5) fold change",
                     .color.label = "Group B vs A",
                     .midpoint = 0) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
    ggplot2::xlim(c(-21,21)) +
    ggplot2::ylim(c(-7,20)))

(tinydenseR::plotPCA(
  .lm.obj = lm.cells,
  .feature =
    ifelse(
      test = condition.stats$fit$coefficients[,"ConditionB"] < 0,
      yes = "abundance down",
      no = "abundance up") |>
    ifelse(
      test = condition.stats$fit$density.weighted.bh.fdr[,"ConditionB"] < 0.1,
      no = "not sig.")  |>
    factor(levels = c("abundance down",
                      "not sig.",
                      "abundance up")),
  .plot.title = "hypothesis testing",
  .color.label = "Group B vs A\nq < 0.1",
  .cat.feature.color = Color.Palette[1,c(1,6,2)],
  .point.size = 1,
  .panel.size = 2) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
    ggplot2::xlim(c(-21,21)) +
    ggplot2::ylim(c(-7,20)))

tinydenseR::plotPCA(
  .lm.obj = lm.cells,
  .feature = lm.cells$graph$clustering$ids,
  .plot.title = "clustering",
  .point.size = 1,
  .panel.size = 2) |> 
  (\(x)
   x +
     ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
     ggplot2::xlim(c(-21,21)) +
     ggplot2::ylim(c(-7,20)) +
     ggplot2::geom_text(data = x$data |>
                          dplyr::group_by(feature) |>
                          dplyr::summarize(PC1 = mean(x = PC1),
                                           PC2 = mean(x = PC2),
                                           .groups = "drop"),
                        mapping = ggplot2::aes(label = feature),
                        size = I(x = 3),
                        color = "black")
  )() 

stat.test.percentages <-
  lm.cells$map$clustering$cell.perc |>
  dplyr::as_tibble() |>
  dplyr::mutate(condition = lm.cells$metadata$Condition) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
  dplyr::group_by(name) |>
  rstatix::t_test(formula = value ~ condition) |>
  dplyr::mutate(p = condition.stats$trad$clustering$fit$adj.p[name,"ConditionB"],
                p.adj = condition.stats$trad$clustering$fit$adj.p[name,"ConditionB"]) |>
  rstatix::add_significance() |>
  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
                               yes = formatC(x = p.adj,
                                             digits = 0,
                                             format = "e"),
                               no = formatC(x = p.adj,
                                            digits = 2,
                                            format = "f"))) |>
  rstatix::add_xy_position(x = "condition")

(tinydenseR::plotTradPerc(
  .lm.obj = lm.cells,
  .x.split = "Condition",
  .x.space.scaler = 0.3
) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5,
                                                       angle = 0)) + 
    ggplot2::labs(title = "cluster %") + 
    ggpubr::stat_pvalue_manual(data = stat.test.percentages, 
                               label = "p.adj",
                               label.size = I(x = 3)) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))) 

(tinydenseR::plotAbundance(
  .lm.obj = lm.cells,
  .x.split = "Condition",
  .x.space.scaler = 0.3
) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5,
                                                       angle = 0)) + 
    ggplot2::labs(title = "within-cluster abundance"))

tinydenseR::plotBeeswarm(
  .lm.obj = lm.cells,
  .stats.obj = condition.stats,
  .coefs = "ConditionB",
  .row.space.scaler = 0.3,
  .perc.plot = FALSE,
  .point.size = 2)

.dea <-
  tinydenseR::get.dea(
    .lm.obj = lm.cells,
    .design = .design,
  )

sort(x = .dea$coefficients[,"ConditionB"],
     decreasing = TRUE)[1] |>
  names() |> 
  (\(x)
   tinydenseR::plotPCA(.lm.obj = lm.cells,
                       .feature = lm.cells$lm[,x],
                       .plot.title = "Top Gene Up in B vs A",
                       .panel.size = 2,
                       .point.size = 1,
                       .color.label = x) +
     ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  )()

sort(x = .dea$coefficients[,"ConditionB"],
     decreasing = FALSE)[1] |>
  names() |> 
  (\(x)
   tinydenseR::plotPCA(.lm.obj = lm.cells,
                       .feature = lm.cells$lm[,x],
                       .plot.title = "Top Gene Down in B vs A",
                       .panel.size = 2,
                       .point.size = 1,
                       .color.label = x) +
     ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  )()

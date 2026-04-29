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

library(princurve)
library(patchwork)
library(SingleCellExperiment)
library(dplyr)
library(tinydenseR)
library(tidyverse)
library(ggpubr)
library(rstatix)

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
rd <- file.path(dirname(script.path), "results", "sim_scRNAseq_trajectory_tdr")

if(!dir.exists(paths = rd)) dir.create(path = rd, recursive = TRUE)

## -------------------------------
## Reproducibility seed
## -------------------------------

set.seed(seed = 42)

## -------------------------------
## User-set parameters (interactive)
## -------------------------------

ngenes <- 500
ncells <- 5000
props <- 
  # proportion of Condition A per milestone
  c(seq(from = 0.5, to = 0.1, by = -0.1),
           seq(from = 0.8, to = 1, by = 0.05))
milestones <- length(x = props)

## -------------------------------
## Load pre-computed simulation (stored as package data)
## To regenerate, see data-raw/sim_trajectory_tdr.R
## -------------------------------

data(sim_trajectory_tdr, package = "tinydenseR")
sim <- sim_trajectory_tdr
rm(sim_trajectory_tdr)

sim_trajectory.meta <- sim$meta
sim_trajectory <- sim$sce

plot(sim_trajectory@int_colData$reducedDims$PCA,
     col = as.factor(sim_trajectory.meta$Condition),
     pch = 16)

# plot proportions along trajectory
(p <-
    ggplot2::ggplot(data = data.frame(A = props,
                                  B = 1 - props,
                                  milestones = 1:milestones) |>
                  tidyr::pivot_longer(cols = c("A","B")),
                mapping = ggplot2::aes(x = as.factor(x = milestones),
                                       y = value*100,
                                       fill = name)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                 plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::labs(
    title = "Cell % along linear trajectory milestones",
    subtitle = "Dev inhibition and accumulation in B at step 5",
    x = paste0("linear trajectory milestones\n(",
               ngenes,
               " genes x ",
               ncells,
               " cells in ",
               unique(x = sim$meta$Sample) |>
                 length(),
               " samples)"),
    y = "% of cells",
    fill = "Condition"
  ) +
  ggplot2::scale_fill_manual(values = unname(obj = tinydenseR::Color.Palette[1,1:2])) +
  ggplot2::geom_bar(position = "stack", 
                    stat = "identity") +
  ggh4x::force_panelsizes(cols = ggplot2::unit(x = 3,
                                               units = "in"),
                          rows = ggplot2::unit(x = 2,
                                               units = "in"))); ggplot2::ggsave(
                                                 plot = p,
                                                                                filename = file.path(rd,
                                                                                                     "ground_truth_traj_settings.png"),
                                                                                width = 4.5, 
                                                                                height = 3.25,
                                                                                dpi = 300); rm(p)


# Create .meta object containing sample-level data
#.meta <-
#  tinydenseR::get.meta(.obj = sim_trajectory,
#                       .sample.var = "Sample",
#                       .verbose = TRUE)
#
## Create .cells object using SCE method
#.cells <- 
#  tinydenseR::get.cells(.exprs = sim_trajectory,
#                        .meta = .meta,
#                        .sample.var = "Sample")[rownames(x = .meta)]
#set.seed(seed = 123)
#lm.cells <-
#  tinydenseR::setup.tdr.obj(
#    .cells = .cells,
#    .meta = .meta,
#    .assay.type = "RNA") |>
#  tinydenseR::get.landmarks(.nHVG = 500) |>
#  tinydenseR::get.graph(.cl.resolution.parameter = 10)
#
#lm.cells <-
#  tinydenseR::get.map(x = lm.cells,
#  .label.confidence = 0.5)

#set.seed(seed = 123)
sim_trajectory <-
  tinydenseR::RunTDR(
    x = sim_trajectory,
    .sample.var = "Sample",
    .assay.type = "RNA",
  .nHVG = 500,
  .cl.resolution.parameter = 0.1)

# NOTE: The explicit workflow above and RunTDR() are numerically equivalent
# for the same inputs, seed, and arguments.  Using RunTDR() followed by
# GetTDR(sim_trajectory) in place of lm.cells will produce identical results.

.design <-
  model.matrix(object = ~ Condition,
               data = tinydenseR::GetTDR(sim_trajectory)@metadata)

# get.lm() returns the updated x with results in $results$lm[[.model.name]]
sim_trajectory <-
  tinydenseR::get.lm(
    x = sim_trajectory,
    .design = .design)

# Fit trajectory for plotting
traj <-
  princurve::principal_curve(x = tinydenseR::GetTDR(sim_trajectory)$landmark.embed$pca$coord[,c("PC1","PC2")],
                             plot = FALSE)

traj.df <-
  data.frame(PC1 = traj$s[traj$ord,1],
             PC2 = traj$s[traj$ord,2])

# Extract count matrices from SCE object for plotting
count_matrices <- list()
for(sample in sort(x = tinydenseR::GetTDR(sim_trajectory)@metadata$Sample)) {
  sample_cells <- SummarizedExperiment::colData(sim_trajectory)$Sample == sample
  count_matrices[[sample]] <- SingleCellExperiment::counts(sim_trajectory)[, sample_cells]
}

(p <- 
    count_matrices |>
    do.call(what = cbind) |>
    Matrix::t() |>
    (\(x)
     log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
    )() |> 
    (\(x)
     (((Matrix::t(x = x[,tinydenseR::GetTDR(sim_trajectory)$pca$HVG]) - tinydenseR::GetTDR(sim_trajectory)$pca$center) /
         tinydenseR::GetTDR(sim_trajectory)$pca$scale) |>
         Matrix::t()) %*%
       tinydenseR::GetTDR(sim_trajectory)$pca$rotation
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
       ggplot2::facet_grid(cols = ggplot2::vars(Condition),
                           labeller = ggplot2::as_labeller(c("A" = "Condition A",
                                                             "B" = "Condition B"))) +
       ggplot2::theme_bw() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                         size = I(x = 14)),
                      legend.position = "none") +
       ggplot2::labs(title = "ground truth differential density") +
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 5)))) +
       ggplot2::scale_color_manual(values = unname(obj = Color.Palette[1,c(1,2)])) +
       ggplot2::geom_density2d(adjust = 1/2,
                               contour_var = "ndensity") +
       ggplot2::geom_path(data = traj.df,
                          color = "black") +
       ggh4x::force_panelsizes(cols = grid::unit(x = 1.5,
                                                 units = "in"),
                               rows = grid::unit(x = 1.5,
                                                 units = "in"))
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd,
                                               "ground_truth_traj.png"),
                          width = 4, 
                          height = 2.75,
                          dpi = 300); rm(p)

cor(x = SingleCellExperiment::reducedDim(x = sim$sce, type = "PCA")[,"PC1"],
    y = SummarizedExperiment::assay(x = sim$sce, slot = "counts") |>
      Matrix::t() |>
      (\(x)
       log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
      )() |>
      as.matrix(),
    method = "spearman")[1,] |>
  sort()

cor(x = SingleCellExperiment::reducedDim(x = sim$sce, type = "PCA")[,"PC2"],
    y = SummarizedExperiment::assay(x = sim$sce, slot = "counts") |>
      Matrix::t() |>
      (\(x)
       log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
      )() |>
      as.matrix(),
    method = "spearman")[1,] |>
  sort()

(p <- 
    #count_matrices |>
    #do.call(what = cbind) |>
    SummarizedExperiment::assay(x = sim$sce, slot = "counts") |>
    Matrix::t() |>
    (\(x)
     log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
    )() |> 
    (\(x)
     #((((Matrix::t(x = x[,tinydenseR::GetTDR(sim_trajectory)$pca$HVG]) - tinydenseR::GetTDR(sim_trajectory)$pca$center) /
     #    tinydenseR::GetTDR(sim_trajectory)$pca$scale) |>
     #    Matrix::t()) %*%
     #  tinydenseR::GetTDR(sim_trajectory)$pca$rotation) |>
       cbind(G2 = x[,"G2"],
             G500 = x[,"G500"],
             G244 = x[,"G244"],
             G302 = x[,"G302"],
             traj = sim$pt
       )
    )() |>
    as.matrix()  |>
    as.data.frame() |> 
    (\(x)
     tidyr::pivot_longer(data = x[,c("traj", "G2", "G500", "G244", "G302")],
                         cols = c("G2", "G500", "G244", "G302"))
    )() |>
    dplyr::mutate(name = factor(x = name,
                              levels = c("G2", "G244", "G500", "G302"))) |>
    (\(x)
              ggplot2::ggplot(data = x,
                              mapping = ggplot2::aes(x = traj,
                                                     y = value,
                                                     color = value)) +
                ggplot2::facet_wrap(facets = ggplot2::vars(name),
                                    ncol = 2,
                                    labeller = ggplot2::as_labeller(c("G2" = "Down",
                                                                      "G500" = "Up",
                                                                      "G244" = "Down and Up",
                                                                      "G302" = "Up and Down"))) +
                ggplot2::theme_bw() +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                                  size = I(x = 14))) +
                ggplot2::labs(title = "gene expression patterns",
                              y = "exprs",
                              color = "exprs") +
                #ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 5)))) +
                ggplot2::scale_color_viridis_c() +
                ggplot2::geom_point(size = I(x = 0.1)) +
                ggh4x::force_panelsizes(cols = grid::unit(x = 1,
                                                          units = "in"),
                                        rows = grid::unit(x = 1,
                                                          units = "in"))         
            
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd,
                                               "ground_truth_traj_gene_expression.png"),
                          width = 3.5, 
                          height = 3.5,
                          dpi = 300); rm(p)

(p <- 
    count_matrices |>
    do.call(what = cbind) |>
    Matrix::t() |>
    (\(x)
     log2(x = (x / (Matrix::rowSums(x = x) / mean(x = Matrix::rowSums(x = x)))) + 1)
    )() |> 
    (\(x)
     (((Matrix::t(x = x[,tinydenseR::GetTDR(sim_trajectory)$pca$HVG]) - tinydenseR::GetTDR(sim_trajectory)$pca$center) /
         tinydenseR::GetTDR(sim_trajectory)$pca$scale) |>
         Matrix::t()) %*%
       tinydenseR::GetTDR(sim_trajectory)$pca$rotation
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
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                         size = I(x = 20))) +
       ggplot2::labs(title = "samples",
                     color = "Condition") +
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 5)))) +
       ggplot2::scale_color_manual(values = unname(obj = Color.Palette[1,c(1,2)])) +
       ggplot2::geom_point(size = I(x = 0.1)) +
       ggplot2::geom_path(data = traj.df,
                          color = "black") +
       ggh4x::force_panelsizes(cols = grid::unit(x = 1,
                                                 units = "in"),
                               rows = grid::unit(x = 1,
                                                 units = "in"))
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd,
                                               "sim.scRNAseq.samples.png"),
                          width = 4.5,
                          height = 3.5,
                          dpi = 300); rm(p)

(p <- 
    tinydenseR::plotPCA(x = sim_trajectory,
                        .feature = "black",
                        .cat.feature.color = "black",
                        .plot.title = "landmarks",
                        .panel.size = 2,
                        .point.size = 2,
                        .color.label = "Group") +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank(),
                   legend.position = "none")); ggplot2::ggsave(plot = p,
                                               filename = file.path(rd,
                                                                    "sim.scRNAseq.landmarks.png"),
                                               width = 3,
                                               height = 3,
                                               dpi = 300,
                                               bg = "white"); rm(p)

tinydenseR::plotPCA(x = sim_trajectory,
                    .feature = tinydenseR::GetTDR(sim_trajectory)$metadata$Condition[tinydenseR::GetTDR(sim_trajectory)$config$key],
                    .cat.feature.color = Color.Palette[1,1:2],
                    .panel.size = 1.5,
                    .point.size = 1,
                    .color.label = "Condition")

(p <- 
    tinydenseR::plotPCA(x = sim_trajectory,
                        .feature = tinydenseR::GetTDR(sim_trajectory)$results$lm[["default"]]$fit$coefficients[,"ConditionB"],
                        .panel.size = 1.5,
                        .point.size = 1,
                        .plot.title = "landmarks",
                        .color.label = "density log2(+0.5)FC",
                        .midpoint = 0,
                        .legend.position = "bottom") +
    ggplot2::geom_path(data = traj.df,
                       color = "black") +
    ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                    title.hjust = 0.5)) +
    ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_blank(),
                   legend.margin = ggplot2::margin(t = -0.1, 
                                                   unit = "in"))); ggplot2::ggsave(plot = p,
                                               filename = file.path(rd,
                                                                    "ground_truth_traj_log2FC.png"),
                                               width = 2.75, 
                                               height = 3,
                                               dpi = 300); rm(p)


head(tinydenseR::GetTDR(sim_trajectory)$density$norm[,order(colnames(tinydenseR::GetTDR(sim_trajectory)$density$norm))],3) |>
  round(2)

tail(tinydenseR::GetTDR(sim_trajectory)$density$norm[,order(colnames(tinydenseR::GetTDR(sim_trajectory)$density$norm))],1) |>
  round(2)

(p <-
    tinydenseR::plotPCA(
      x = sim_trajectory,
      .feature =
        ifelse(
          test = tinydenseR::GetTDR(sim_trajectory)$results$lm[["default"]]$fit$coefficients[,"ConditionB"] < 0,
          yes = "up in A",
          no = "up in B") |>
        ifelse(
          test = tinydenseR::GetTDR(sim_trajectory)$results$lm[["default"]]$fit$pca.weighted.q[,"ConditionB"] < 0.1,
          no = "not sig.")  |>
        factor(levels = c("up in A",
                          "not sig.",
                          "up in B")),
      .plot.title = "landmarks",
      .cat.feature.color = Color.Palette[1,c(1,6,2)],
      .point.size = 1,
      .panel.size = 1.5, 
      .legend.position = "bottom") +
    ggplot2::geom_path(data = traj.df,
                       color = "black") +
    ggplot2::guides(color = ggplot2::guide_legend(title.position = "top",
                                                  title.hjust = 0.5,
                                                  override.aes = list(size = I(x = 5)))) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
    ggplot2::labs(color = "q < 0.1")); ggplot2::ggsave(plot = p,
                                               filename = file.path(rd,
                                                                    "ground_truth_traj_abundance.pca.weighted.q.png"),
                                               width = 2.75, 
                                               height = 3,
                                               dpi = 300, 
                                               bg = "white"); rm(p)

#(p <-
#    tinydenseR::plotPCA(
#      x = sim_trajectory,
#      .feature = tinydenseR::GetTDR(sim_trajectory)$landmark.annot$clustering$ids,
#      .plot.title = "clustering",
#      .point.size = 1,
#      .panel.size = 1.5) |> 
#    (\(x)
#     x +
#       ggplot2::geom_path(data = traj.df,
#                          color = "black") +
#       ggplot2::theme(plot.subtitle = ggplot2::element_blank()) +
#       ggplot2::geom_text(data = x$data |>
#                            dplyr::group_by(feature) |>
#                            dplyr::summarize(PC1 = mean(x = PC1),
#                                             PC2 = mean(x = PC2),
#                                             .groups = "drop"),
#                          mapping = ggplot2::aes(label = feature),
#                          size = I(x = 3),
#                          color = "black")
#    )()); ggplot2::ggsave(plot = p,
#                          filename = file.path(rd,
#                                               "ground_truth_traj_cl.png"),
#                          width = 4, 
#                          height = 2.5,
#                          dpi = 300); rm(p)
#
#stat.test.percentages <-
#  tinydenseR::GetTDR(sim_trajectory)$density$composition$clustering$cell.perc |>
#  dplyr::as_tibble() |>
#  dplyr::mutate(condition = tinydenseR::GetTDR(sim_trajectory)$metadata$Condition) |>
#  tidyr::pivot_longer(cols = dplyr::starts_with(match = "cluster.")) |>
#  dplyr::group_by(name) |>
#  rstatix::t_test(formula = value ~ condition) |>
#  dplyr::mutate(p = tinydenseR::GetTDR(sim_trajectory)$results$lm[["default"]]$trad$clustering$fit$adj.p[name,"ConditionB"],
#                p.adj = tinydenseR::GetTDR(sim_trajectory)$results$lm[["default"]]$trad$clustering$fit$adj.p[name,"ConditionB"]) |>
#  rstatix::add_significance() |>
#  dplyr::mutate(p.adj = ifelse(test = p.adj < 0.01,
#                               yes = formatC(x = p.adj,
#                                             digits = 0,
#                                             format = "e"),
#                               no = formatC(x = p.adj,
#                                            digits = 2,
#                                            format = "f"))) |>
#  rstatix::add_xy_position(x = "condition")
#
#(p <- 
#    tinydenseR::plotTradPerc(
#      x = sim_trajectory,
#      .x.split = "Condition",
#      .x.space.scaler = 0.3
#    ) + 
#    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5,
#                                                       angle = 0)) + 
#    ggplot2::labs(title = "cluster %") + 
#    ggpubr::stat_pvalue_manual(data = stat.test.percentages, 
#                               label = "p.adj",
#                               label.size = I(x = 3)) +
#    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))); ggplot2::ggsave(plot = p,
#                                                                                            filename = file.path(rd, 
#                                                                                                                 "TradPerc_traj.png"),
#                                                                                            width = 4, 
#                                                                                            height = 3,
#                                                                                            dpi = 300); rm(p)
#
#(p <-
#    tinydenseR::plotDensity(
#      x = sim_trajectory,
#      .x.split = "Condition",
#      .x.space.scaler = 0.3
#    ) + 
#    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5,
#                                                       angle = 0)) + 
#    ggplot2::labs(title = "within-cluster density")); ggplot2::ggsave(plot = p,
#                                                                      filename = file.path(rd, 
#                                                                                           "Abundance_traj.png"),
#                                                                      width = 4, 
#                                                                      height = 3,
#                                                                      dpi = 300); rm(p)
#
#(p <-
#    tinydenseR::plotBeeswarm(
#      x = sim_trajectory,
#      .model.name = "default",
#      .coefs = "ConditionB",
#      .row.space.scaler = 0.3,
#      .perc.plot = FALSE,
#      .point.size = 2,
#      .q.from = "pca.weighted.q") +
#    ggplot2::labs(title = "Condition B vs A")); ggplot2::ggsave(plot = p,
#                                                                filename = file.path(rd, 
#                                                                                     "Beeswarm_traj.pca.weighted.q.png"),
#                                                                width = 4, 
#                                                                height = 3,
#                                                                dpi = 300); rm(p)
#
sim_trajectory <-
  tinydenseR::get.pbDE(
    x = sim_trajectory,
    .design = .design,
  )

(p <-
    sort(x = tinydenseR::GetTDR(sim_trajectory)$pbDE$default$all$coefficients[,"ConditionB"],
         decreasing = TRUE)[1] |>
    names() |> 
    (\(x)
     tinydenseR::plotPCA(x = sim_trajectory,
                         .feature = tinydenseR::GetTDR(sim_trajectory)$landmarks[,x],
                         .plot.title = "Top Gene Up in B vs A",
                         .panel.size = 2,
                         .point.size = 1,
                         .color.label = x) +
       ggplot2::geom_path(data = traj.df,
                          color = "black") +
       ggplot2::theme(plot.subtitle = ggplot2::element_blank())
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd, 
                                               "sim.scRNAseq_traj_toUp.png"),
                          width = 3.5, 
                          height = 3,
                          dpi = 300); rm(p)

(p <-
    sort(x = tinydenseR::GetTDR(sim_trajectory)$pbDE$default$all$coefficients[,"ConditionB"],
         decreasing = FALSE)[1] |>
    names() |> 
    (\(x)
     tinydenseR::plotPCA(x = sim_trajectory,
                         .feature = tinydenseR::GetTDR(sim_trajectory)$landmarks[,x],
                         .plot.title = "Top Gene Down in B vs A",
                         .panel.size = 2,
                         .point.size = 1,
                         .color.label = x) +
       ggplot2::geom_path(data = traj.df,
                          color = "black") +
       ggplot2::theme(plot.subtitle = ggplot2::element_blank())
    )()); ggplot2::ggsave(plot = p,
                          filename = file.path(rd, 
                                               "sim.scRNAseq_traj_toDown.png"),
                          width = 3.5, 
                          height = 3,
                          dpi = 300); rm(p)

# Create reduced model to embed samples quantitatively along the Condition axis
noCondition.design <-
  model.matrix(object = ~ 1,
               data = tinydenseR::GetTDR(sim_trajectory)$metadata)

# Store reduced model stats with a different .model.name
sim_trajectory <- 
  tinydenseR::get.lm(
    x = sim_trajectory,
    .design = noCondition.design,
    .model.name = "noCondition",
    .verbose = FALSE 
  )

# Compute sample embeddings using nested model comparison
# Embeddings stored in tinydenseR::GetTDR(sim_trajectory)$results$embedding$pePC[[slot.name]]
sim_trajectory <-
  tinydenseR::get.embedding(
    x = sim_trajectory,
    .full.model = "default",
    .red.model = "noCondition",
    .term.of.interest = "Condition",
    .verbose = FALSE 
  )

# Embed samples based on major axis of variance
(smpl.pca.traj <- tinydenseR::plotSampleEmbedding(
  x = sim_trajectory,
  .embedding = "pca",
  .color.by = "Condition",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2)],
  .panel.size = 1.5,
  .point.size = 3
) +
    ggplot2::labs(title = "PCA") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom")); ggplot2::ggsave(plot = smpl.pca.traj,
                                                                 filename = file.path(rd,
                                                                                      "sample_embedding_pca.png"),
                                                                 width = 2.5, 
                                                                 height = 3,
                                                                 dpi = 300,
                                                                 bg = "white")

(p <-
    tinydenseR::plotSampleEmbedding(
      x = sim_trajectory,
      .embedding = "pca",
      .color.by = "Replicate",
      .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2)],
      .panel.size = 1.5,
      .point.size = 3
    ) +
    ggplot2::labs(title = "PCA",
                  color = "Batch") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom")); ggplot2::ggsave(plot = p,
                                                                 filename = file.path(rd,
                                                                                      "sample_embedding_pca_batch.png"),
                                                                 width = 2.5, 
                                                                 height = 3,
                                                                 dpi = 300,
                                                                 bg = "white"); rm(p)

#(p <-
#    tinydenseR::plotSampleEmbedding(
#      x = sim_trajectory,
#      .embedding = "traj",
#      .color.by = "Condition",
#      .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2)],
#      .panel.size = 1.5,
#      .point.size = 3
#    ) +
#    ggplot2::labs(title = "trajectory") +
#    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
#                   legend.position = "bottom")); ggplot2::ggsave(plot = p,
#                                                                 filename = file.path(rd,
#                                                                                      "sample_embedding_traj.png"),
#                                                                 width = 2.5, 
#                                                                 height = 3,
#                                                                 dpi = 300,
#                                                                 bg = "white"); rm(p)

#(p <-
#    tinydenseR::plotSampleEmbedding(
#      x = sim_trajectory,
#      .embedding = "traj",
#      .color.by = "Replicate",
#      .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2)],
#      .panel.size = 1.5,
#      .point.size = 3
#    ) +
#    ggplot2::labs(title = "traj",
#                  color = "Batch") +
#    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
#                   legend.position = "bottom")); ggplot2::ggsave(plot = p,
#                                                                 filename = file.path(rd,
#                                                                                      "sample_embedding_traj_batch.png"),
#                                                                 width = 2.75, 
#                                                                 height = 3,
#                                                                 dpi = 300,
#                                                                 bg = "white"); rm(p)

# Supervised embedding of samples
(smpl.pePC.traj <- tinydenseR::plotSampleEmbedding(
  x = sim_trajectory,
  .embedding = "pePC",
  .sup.embed.slot = "Condition",
  .color.by = "Condition",
  .cat.feature.color = tinydenseR::Color.Palette[1,c(1,2)],
  .panel.size = 1.5,
  .point.size = 3
) +
    ggplot2::labs(title = "partial-effect PC") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "bottom")); ggplot2::ggsave(plot = smpl.pePC.traj,
                                                                 filename = file.path(rd,
                                                                                      "sample_embedding_pePC.png"),
                                                                 width = 2.5, 
                                                                 height = 3,
                                                                 dpi = 300,
                                                                 bg = "white")

(p <-
    ((smpl.pca.traj | 
        smpl.pePC.traj) +
       patchwork::plot_layout(guides = "collect") &
       ggplot2::theme(legend.position = "bottom",
                      legend.justification = "center"))); ggplot2::ggsave(
                        plot = p,
                        filename = file.path(rd,
                                             "sample_embed_traj.png"),
                        width = 4.5,
                        height = 3,
                        units = "in",
                        dpi = 300,
                        bg = "white"); rm(p)

data.frame(coef = tinydenseR::GetTDR(sim_trajectory)$results$lm[["default"]]$fit$coefficients[,"ConditionB"],
           rotation = tinydenseR::GetTDR(sim_trajectory)$sample.embed$pepc$Condition$rotation[,"PC1"]) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = coef,
                                          y = rotation)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
     ggplot2::labs(title = "landmarks",
                   x = "abundance log2(+0.5)FC",
                   y = "partial-effect loading") +
     ggpubr::stat_cor(label.y = min(x$rotation),
                      size = 3) +
     ggplot2::geom_point(size = I(x = 3),
                         color = "black") +
     ggplot2::geom_smooth(method = "lm",
                          color = "red",
                          se = TRUE)  +
     ggh4x::force_panelsizes(cols = grid::unit(x = 2,
                                               units = "in"),
                             rows = grid::unit(x = 2,
                                               units = "in"))
  )()

## =============================================================================
## Run t-test on pePC1 to confirm it captures Condition variation (sanity check)
## =============================================================================

t.test(
  formula = tinydenseR::GetTDR(sim_trajectory)$sample.embed$pepc$Condition$coord[,1] ~ 
    tinydenseR::GetTDR(sim_trajectory)$metadata$Condition)$p.value

sim_trajectory <-
  tinydenseR::get.plsD(
    x = sim_trajectory, 
    .coef.col = "ConditionB",
    .model.name = "default",
    .verbose = TRUE
  )

(p <-
    plotPlsD(
      x = sim_trajectory,
      .coef.col = "ConditionB"
    )); rm(p)

plsD1.p <-
    plotPlsD(
      x = sim_trajectory,
      .coef.col = "ConditionB",
      .plsD.dim = 1,
      .embed = "pca",
      .panel.size = 1.5
    ) #&
    #ggplot2::theme(plot.title = ggplot2::element_blank(),
    #               plot.subtitle = ggplot2::element_blank())) +
    #patchwork::plot_annotation(title = "plsD1")

plsD1.p[[1]] <-
  plsD1.p[[1]] +
  ggplot2::labs(title = "",
                subtitle = "") +
  ggplot2::geom_path(data = traj.df,
                     color = "black")

plsD1.p[[2]] <-
  plsD1.p[[2]] +
  ggplot2::labs(title = "",
                subtitle = "")
  
plsD1.p <-
  plsD1.p + 
  patchwork::plot_annotation(
    title = "Graph-Diffused, Density Contrast-Aligned\nPLS Decomposition") &
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, 
                                                    margin = ggplot2::margin(t = -0.1,
                                                                             unit = "in"),
                                                    size = I(x = 12)),
                 plot.subtitle = ggplot2::element_blank())

plsD1.p; ggplot2::ggsave(
  plot = plsD1.p,
  filename = file.path(rd,
                       "sim.scRNAseq.plsD1.png"),
  width = 4, 
  height = 3.75,
  dpi = 300)

(p <-
    plotPlsD(
      x = sim_trajectory,
      .coef.col = "ConditionB",
      .plsD.dim = 2,
      .embed = "pca",
      .panel.size = 1.5
    )); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "sim.scRNAseq.plsD2.png"),
                        width = 4.6, 
                        height = 3.25,
                        dpi = 300); rm(p)


(p <-
    tinydenseR::plotPlsD(
      x = sim_trajectory,
      .coef.col = "ConditionB",
      .plsD.dim = 3,
      .embed = "pca",
      .panel.size = 1.5
    )); rm(p)

(p <-
    plotPlsDHeatmap(
  x = sim_trajectory,
  .coef.col = "ConditionB",
  .order.by = "plsD.dim",
  .add.annot = cbind(traj = -(traj$lambda - mean(x = traj$lambda))),
  .plsD.dim = 1,
  .panel.height = 3,
  .feature.font.size = 4
)); ggplot2::ggsave(
  plot = p,
  filename = file.path(rd,
                       "sim.scRNAseq.plsD1.hm.png"),
  width = 6, 
  height = 4.75,
  dpi = 300,
  bg = "white"); rm(p)


# permutation tests
file.path(script.path |>
  dirname(),
"perm_utils.R") |>
  source()

# Single condition structure: Condition B vs A
# Create permutation metadata with Treatment/Batch columns
# matching perm_utils.R conventions
.perm.meta <- data.frame(
  Treatment = tinydenseR::GetTDR(sim_trajectory)@metadata$Condition,
  Batch = 1,
  row.names = rownames(tinydenseR::GetTDR(sim_trajectory)@metadata)
)

results <-
  run_stratified_permutation_test(
    lm_obj = sim_trajectory,
    meta_df = .perm.meta,
    formula = ~ Treatment,
    coef_name = "B")

results$condition <- "ConditionB"

# Sanity check MCC values
cat("\n=== MCC Validation ===\n")
cat("  Observed MCC:", unique(results$mcc[results$type == "observed"]), "\n")
cat("  Complement MCC:", unique(results$mcc[results$type == "complement"]), "\n")
cat("  MCC range for permuted:", 
    paste(range(results$mcc[results$type == "permuted"], na.rm = TRUE), collapse = " to "), "\n")

# Verify expected values
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
(p <-
  plot_n_sig(results = results) &
  theme(plot.title = element_blank())); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "perm_test_n_sig.png"),
                        width = 3, 
                        height = 4,
                        dpi = 300); rm(p)

(p <-
  plot_min_q(results = results) &
  theme(plot.title = element_blank())); ggplot2::ggsave(plot = p,
                        filename = file.path(rd,
                                             "perm_test_min_q.png"),
                        width = 3, 
                        height = 4,
                        dpi = 300); rm(p)

## =============================================================================
## Permutation test: pePC sample embedding (t-test on pePC1)
## =============================================================================

results_pePC <-
  run_stratified_permutation_test_pePC(
    lm_obj = sim_trajectory,
    meta_df = .perm.meta,
    formula = ~ Treatment,
    coef_name = "B",
    term_of_interest = "Condition")

results_pePC$condition <- "ConditionB"

# Quick summary
cat("\n=== pePC permutation test results ===\n")
results_pePC |>
  dplyr::group_by(condition, type) |>
  dplyr::summarise(
    mean_sig    = mean(n_sig),
    mean_pePC1_p = mean(pePC1_p, na.rm = TRUE),
    mean_mcc    = mean(mcc, na.rm = TRUE),
    .groups = "drop") |>
  print()

sessionInfo()
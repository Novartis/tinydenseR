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

## ---------------------------------------------------------------
## Code to prepare `sim_trajectory_tdr` dataset
##
## This script reproduces the simulation from
##   inst/analysis/sim_scRNAseq_trajectory_tdr.R
## and saves the result as package data.
## ---------------------------------------------------------------

library(dyntoy)
library(SingleCellExperiment)
library(irlba)
library(dplyr)

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
## Simulation function
## -------------------------------

simulate_linear_trajectory <- function(n.milestones = 3, num.cells, num.features, p.vec) {
  
  dataset <- generate_dataset(
    model = model_linear(num_milestones = n.milestones),
    num_cells = num.cells,
    num_features = num.features,
    sample_mean_count = function()
      runif(n = 1, min = 10, max = 100),
    differentially_expressed_rate = 0.8
  )
  
  # reorder genes along trajectory to facilitate identification later
  # first get princurve
  init.pca <-
    prcomp_irlba(x = dataset$expression,
                 center = TRUE,
                 n = 2)$x
  
  pt <-
    princurve::principal_curve(x = init.pca,
                               plot = FALSE)$lambda
  
  # get order
  reord.de.g <-
    cor(x = pt,
        y = as.matrix(x = dataset$expression),
        method = "spearman")[1,] |>
    order()
  
  # re-order
  dataset$expression <-
    dataset$expression[,reord.de.g] |> 
    (\(x)
     `colnames<-`(x = x,
                  value = paste0("G",
                                 1:500))
    )()
  
  dataset$counts <-
    dataset$counts[,reord.de.g] |> 
    (\(x)
     `colnames<-`(x = x,
                  value = paste0("G",
                                 1:500))
    )()
  
  ## Build metadata
  cnts <- 
    Matrix::t(x = dataset$counts)
  
  branches <-
    dataset$prior_information$groups_id
  
  coldata_df <-
    data.frame(cell_id = colnames(cnts)) |>
    dplyr::left_join(y = branches)
  
  ## Assign condition (same logic as original)
  n_groups <- 
    unique(x = coldata_df$group_id) |>
    length()
  
  a.cells <- 
    vector(mode = "integer",
           length = nrow(x = coldata_df))
  
  for (i in seq_len(n_groups)) {
    g <- paste0("M", i)
    p <- p.vec[i]
    m.A <- sample(
      coldata_df$cell_id[coldata_df$group_id == g],
      size = floor(sum(coldata_df$group_id == g) * p)
    )
    a.cells <- c(a.cells, m.A)
  }
  
  coldata_df <- 
    coldata_df |>
    mutate(Condition = ifelse(cell_id %in% a.cells, "A", "B"))
  
  ## -------------------------------
  ## Replicates
  ## -------------------------------
  
  coldata_df <- coldata_df |>
    group_by(group_id) |>
    mutate(
      Replicate = c(
        rep("R1", floor(n() * 0.3)),
        rep("R2", floor(n() * 0.3)),
        rep("R3", n() - 2 * floor(n() * 0.3))
      )
    )
  
  coldata_df$Sample <- paste(coldata_df$Condition, coldata_df$Replicate, sep = "_")
  
  ## -------------------------------
  ## Dimensionality reduction & SCE
  ## -------------------------------
  
  pca <- prcomp_irlba(
    dataset$expression,
    n = 50,
    scale. = TRUE,
    center = TRUE
  )
  
  sce <- SingleCellExperiment(
    assays = list(
      counts = t(dataset$counts),
      logcounts = t(dataset$expression)
    ),
    colData = as.data.frame(coldata_df),
    reducedDims = list(PCA = pca$x)
  )
  
  list(sce = sce,
       meta = as.data.frame(x = coldata_df),
       pt = pt)
}

## -------------------------------
## Run simulation
## -------------------------------

sim_trajectory_tdr <- simulate_linear_trajectory(
  n.milestones = milestones,
  num.cells = ncells,
  num.features = ngenes,
  p.vec = props
)

## -------------------------------
## Save as package data
## -------------------------------

usethis::use_data(sim_trajectory_tdr, overwrite = TRUE)

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

# Compare Seurat reference mapping vs tinydenseR on COVID PBMC datasets
wd <- 
  "path/to/your/working/directory/" # change to your working directory
rd <-
  file.path(wd,
            "res")

setwd(dir = wd)

if(!(file.path(wd,
               "Seurat.data") |>
     dir.exists())){
  file.path(wd,
            "Seurat.data") |>
    dir.create(recursive = TRUE,
               showWarnings = TRUE)
}

# Follow instructions to create Seurat object from h5ad files here
# https://github.com/satijalab/seurat/blob/e35abd442520808a20025e589f861620ddc315af/vignettes/seurat5_bpcells_interaction_vignette.Rmd
# and just to line 83 here: 
# https://github.com/satijalab/seurat/blob/30f82df52159ac5f0feb80b149698abbd876b779/vignettes/COVID_SCTMapping.Rmd

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
  )()

# Build symphony reference to use with tinydenseR
reference.symphony <-
  symphony::buildReference(
    exp_ref = reference[["SCT"]]$counts, # ONLY SCT FOUND IN THIS DATASET!
    metadata_ref = reference@meta.data,
    verbose = TRUE,
    save_uwot_path = file.path(wd ,"Seurat.data/symphony_ref_uwot_model"), # file path to save uwot UMAP model
    seed = 123
  )

# get intersect of genes
.genes.names <-
  names(x = object@assays$RNA@layers) |>
  grep(pattern = "counts",
       value = TRUE,
       fixed = TRUE) |>
  lapply(FUN = function(dataset){
    
    SeuratObject::LayerData(object = object, 
                            assay = "RNA", 
                            layer = dataset) |>
      rownames()
  }) |>
  Reduce(f = base::intersect)

# get the .cells object
if(!(file.path(wd,
               "Seurat.data/RDS") |>
     dir.exists())){
  file.path(wd,
            "Seurat.data/RDS") |>
    dir.create(recursive = TRUE,
               showWarnings = TRUE)
}

.cells <-
  # create a named character vector with each sample
  unique(x = object$donor_id_disease) |>
  (\(x)
   setNames(object = x,
            nm = x)
  )() |>
  # loop over each sample
  lapply(FUN = function(smpl){
    
    # define the on-disk location to save the matrix
    uri <-
      file.path(wd,
                "Seurat.data/RDS",
                paste0(smpl,
                       ".RDS"))
    
    # loop over each dataset and combine cells from the same sample
    names(x = object@assays$RNA@layers) |>
      grep(pattern = "counts",
           value = TRUE,
           fixed = TRUE) |>
      lapply(FUN = function(dataset){
        
        .cells.in.smpl <-
          object$donor_id_disease[
            object$publication == gsub(pattern = "counts.",
                                       replacement = "",
                                       x = dataset,
                                       fixed = TRUE)
          ] == smpl
        
        if(sum(.cells.in.smpl) == 0){
          return(NULL)
        }
        
        # extract the count matrix for the sample
        SeuratObject::LayerData(object = object, 
                                assay = "RNA", 
                                layer = dataset)[.genes.names,
                                                 .cells.in.smpl] |>
          as(Class = "dgCMatrix") |>
          (\(x)
           # save each matrix to the on-disk location
           saveRDS(
             object = x,
             file = uri,
             compress = FALSE)
          )()
        
        # progress
        (which(x = unique(x = object$donor_id_disease)  == smpl) * 100 / 
            (unique(x = object$donor_id_disease) |>
               length())) |>
          round(digits = 2) |>
          print()
        
        return(uri)
        
      }) |> 
      (\(x)
       x[lengths(x = x) > 0]
      )()
    
  }) |>
  unlist() |>
  as.list()

#.cells <-
#  list.files(path = file.path(wd,"Seurat.data/RDS"),
#             pattern = "\\.RDS$",
#             full.names = TRUE) |>
#  (\(x)
#   setNames(object = x,
#            nm = basename(path = x) |>
#              gsub(pattern = "\\.RDS$",
#                   replacement = "")) 
#  )() |>
#  as.list()

# get .meta object
.meta <-
  object@meta.data[,c("publication","sex","donor_id", "disease","donor_id_disease")] |>
  dplyr::distinct() |>
  (\(x)
   `rownames<-`(x = x, 
                value = x$donor_id_disease)
  )() |> 
  (\(x)
   x[match(x = names(x = .cells),
           table = rownames(x = x)),]
  )() |>
  droplevels()

covid.lm.cells <-
  tinydenseR::setup.lm.obj(
    .meta = .meta[.meta$donor_id_disease %in%
                    (table(object$donor_id_disease) |> 
                       (\(x)
                        names(x = x[x > 1000])
                       )()),],
    .cells = .cells[.meta$donor_id_disease %in% 
                      (table(object$donor_id_disease) |> 
                         (\(x)
                          names(x = x[x > 1000])
                         )())],
    .harmony.var = "publication",
    .assay.type = "RNA",
    .verbose = TRUE) |>
  tinydenseR::get.landmarks(
    .verbose = TRUE)

covid.lm.cells <-
  tinydenseR::get.graph(
    .lm.obj = covid.lm.cells,
    .cl.resolution.parameter = 2,
    .verbose = TRUE)

tinydenseR::plotUMAP(
  .lm.obj = covid.lm.cells,
  .panel.size = 2
)

tinydenseR::plotPCA(
  .lm.obj = covid.lm.cells,
  .feature = covid.lm.cells$metadata$publication[covid.lm.cells$key],
  .panel.size = 2
)

tinydenseR::plotUMAP(
  .lm.obj = covid.lm.cells,
  .feature = covid.lm.cells$metadata$publication[covid.lm.cells$key],
  .panel.size = 2
)

covid.lm.cells <-
  tinydenseR::get.map(
    .lm.obj = covid.lm.cells,
    .ref.obj = reference.symphony,
    .integrate.vars = "LibraryId",
    .celltype.col.name = "celltype.l2",
    .verbose = TRUE)

tinydenseR::plotUMAP(
  .lm.obj = covid.lm.cells,
  .feature = covid.lm.cells$graph$celltyping$ids,
  .panel.size = 2
)

tinydenseR::plotSamplePCA(
  .lm.obj = covid.lm.cells,
  .labels.from = "disease",
  .point.size = 1
)

tinydenseR::plotSamplePCA(
  .lm.obj = covid.lm.cells,
  .labels.from = "publication",
  .point.size = 1
)

tinydenseR::plotSamplePCA(
  .lm.obj = covid.lm.cells,
  .labels.from = "sex",
  .point.size = 1
)

tinydenseR::plotSamplePCA(
  .lm.obj = covid.lm.cells,
  .labels.from = "log10.n.cells",
  .point.size = 1
)

COVID.tdr.stat.percentages <-
  as.data.frame(x = covid.lm.cells$map$celltyping$cell.perc) |>
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
  as.data.frame(x = covid.lm.cells$map$celltyping$cell.perc) |>
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
  )()

disease.covid.design <-
  stats::model.matrix(object = ~ 0 + disease + publication + sex,
                      data = covid.lm.cells$metadata)

# check group names
colnames(x = disease.covid.design) <-
  colnames(x = disease.covid.design) |>
  gsub(pattern = "^disease|^publication|^sex",
       replacement = "",
       fixed = FALSE) |>
  make.names()

# set contrasts
disease.covid.contrasts <- 
  limma::makeContrasts(
    COVID.19 - normal,
    levels = colnames(x = disease.covid.design)
  )

# get stats
disease.covid.stats.res <-
  tinydenseR::get.stats(
    .lm.obj = covid.lm.cells,
    .design = disease.covid.design,
    .contrasts = disease.covid.contrasts
  )

tinydenseR::plotBeeswarm(
    .lm.obj = covid.lm.cells,
    .stats.obj = disease.covid.stats.res,
    .coefs = "COVID.19 - normal",
    .split.by = "clustering",
    .swarm.title = "PBMC",
    .legend.position = "left",
    .FDR = 0.1,
    .perc.plot = FALSE) +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::labs(subtitle = "COVID.19 vs normal")


table(covid.lm.cells$graph$clustering$ids,
      covid.lm.cells$graph$celltyping$ids) |>
  as.data.frame.matrix() |> 
  (\(x)
   dplyr::mutate(.data = x,
                 cluster = rownames(x = x))
  )() |>
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "celltype",
                      values_to = "n_cells") |>
  dplyr::group_by(celltype, cluster) |>
  dplyr::summarize(n_cells = sum(n_cells),
                   .groups = "drop") |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(perc_cluster = n_cells / sum(n_cells) * 100) |>
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = celltype,
                                          y = cluster,
                                          label = round(x = perc_cluster,
                                                        digits = 0))) +
     ggplot2::theme_bw(base_size = 10) +
     ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                        vjust = 1,
                                                        angle = 30)) +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    axis.title.x = ggplot2::element_blank(),
                    axis.title.y = ggplot2::element_blank()) +
     ggplot2::labs(title = "cell typing via ref. mapping",
                   fill = "% of cluster\n(row sums = 100)") +
     ggplot2::scale_y_discrete(limits = rev) +
     ggplot2::scale_fill_viridis_c() +
     ggplot2::geom_tile(mapping = ggplot2::aes(fill = perc_cluster)) +
     ggplot2::geom_text(size = I(x = 3))+
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = (unique(x = x$celltype) |>
                                                         length()) * 0.2,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 3,
                                                  units = "in"))
   
  )()

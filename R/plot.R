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

# Helper function to check if ggiraph is available
.is_ggiraph_available <- function() {
  requireNamespace("ggiraph", quietly = TRUE)
}

#' Plot PCA
#'
#' Visualizes landmarks in PCA space with flexible coloring by features, clusters, or statistical 
#' results. Supports both continuous (e.g., gene expression, fold changes) and categorical 
#' (e.g., clusters, cell types) overlays. Optional interactive hover shows landmark feature signatures.
#'
#' @param .tdr.obj A tinydenseR object with PCA computed via \code{get.landmarks()}.
#' @param .PC.x Integer specifying x-axis principal component (default 1).
#' @param .PC.y Integer specifying y-axis principal component (default 2).
#' @param .feature Numeric vector or factor of length \code{nrow(.tdr.obj$landmarks)} to color points by. 
#'   Defaults to cluster IDs from \code{$graph$clustering$ids}. Can be:
#'   \itemize{
#'     \item Numeric: gene expression, fold changes, statistics (uses diverging color scale)
#'     \item Factor: clusters, cell types, conditions (uses discrete color palette)
#'   }
#' @param .cat.feature.color Character vector of colors for categorical features. Default uses 
#'   \code{Color.Palette[1,1:5]} (5-color palette, interpolated to match factor levels).
#' @param .panel.size Numeric panel width/height in inches. Default 2 for numeric, 3 for categorical 
#'   (to accommodate legend).
#' @param .midpoint Numeric midpoint for diverging color scale (continuous features only). 
#'   Defaults to median of \code{.feature}.
#' @param .plot.title Character plot title (default "").
#' @param .color.label Character legend title (default "").
#' @param .legend.position Character legend position: "right" (default), "left", "top", "bottom", or "none".
#' @param .point.size Numeric point size (default 0.1).
#' @param .seed Integer random seed for plot point ordering (default 123).
#' @param .hover.stats Character specifying hover information: "none" (default) or "marker" (shows 
#'   landmark feature signatures from \code{get.lm.features.stats()}). Requires \pkg{ggiraph}.
#'   
#' @return Plot object (class depends on interactivity):
#'   \describe{
#'     \item{Static ggplot}{\code{.hover.stats = "none"} returns a \code{ggplot} object}
#'     \item{Interactive girafe}{\code{.hover.stats != "none"} returns a \code{girafe} object (if \pkg{ggiraph} installed), otherwise falls back to static \code{ggplot} with warning}
#'   }
#'   
#' @note Interactive hover features require the \pkg{ggiraph} package. Install with 
#'   \code{install.packages("ggiraph")}.
#'   
#' @seealso \code{\link{plotUMAP}} for UMAP visualization, \code{\link{get.lm.features.stats}} 
#'   for hover feature computation
#' 
#' @examples
#' \dontrun{
#' # From README: Basic PCA visualization
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500, .nPC = 3) |>
#'   get.graph()
#' 
#' plotPCA(lm.cells, .point.size = 1, .panel.size = 1.5)
#' 
#' # Color by fold change from get.lm()
#' condition.stats <- get.lm(lm.cells, .design = design)
#' plotPCA(lm.cells, 
#'         .feature = condition.stats$fit$coefficients[,"ConditionB"],
#'         .color.label = "log2 FC", 
#'         .midpoint = 0)
#' 
#' # Interactive with feature hover
#' lm.cells <- get.lm.features.stats(lm.cells)
#' plotPCA(lm.cells, .hover.stats = "marker")
#' }
#' 
#' @export
plotPCA <-
  function(.tdr.obj,
           .PC.x = 1,
           .PC.y = 2,
           .feature = .tdr.obj$graph$clustering$ids,
           .cat.feature.color = Color.Palette[1,1:5],
           .panel.size = if(is.numeric(x = .feature)) 2 else 3,
           .midpoint = NULL,
           .plot.title = "",
           .color.label = "",
           .legend.position = "right",
           .point.size = 0.1,
           .seed = 123,
           .hover.stats = "none"){
    
    # R CMD check appeasement
    topFeatTab <- feature <- NULL
    
    # Validate PC selection
    if(any(!c(.PC.x, .PC.y) %in%
           1:ncol(x = .tdr.obj$pca$embed))){
      
      stop(".PC.x and .PC.y must be valid PC indices.\n",
           "Available PCs: ", paste(x = 1:ncol(x = .tdr.obj$pca$embed), collapse = ", "),
           "\nProvided: .PC.x = ", .PC.x, ", .PC.y = ", .PC.y)
    }
    
    .hover.stats <-
      match.arg(arg = .hover.stats,
                choices = c(
                  "none",
                  "marker"
                ))
    
    # Convert PC indices to column names
    .PC.x <-
      colnames(x = .tdr.obj$pca$embed)[.PC.x]
    
    .PC.y <-
      colnames(x = .tdr.obj$pca$embed)[.PC.y]
    
    # Set midpoint to median for continuous features
    if(is.null(x = .midpoint) &
       is.numeric(x = .feature)){
      .midpoint <- stats::median(x = .feature)
    }
    
    # Build plotting data frame
    set.seed(seed = .seed)
    dat.df <-
      as.data.frame(x = .tdr.obj$pca$embed[,c(.PC.x,.PC.y)]) |>
      cbind(feature = .feature,
            point.size = .point.size)
    
    # Attach hover tooltips if requested
    if(.hover.stats == "marker"){
      
      dat.df$topFeatTab <-
        .tdr.obj$interact.plot$lm.features$html
      
    }
    
    # Randomize row order to avoid plotting bias (one class on top)
    dat.df <-
      dat.df  |>
      (\(x)
       # See: https://stackoverflow.com/a/29325361
       x[nrow(x = x) |>
           seq_len() |>
           sample(),]
      )()
    
    # Create base plot
    p <-
      dat.df  |>
      (\(x)
       ggplot2::ggplot(data = x,
                       mapping = ggplot2::aes(x = !!rlang::sym(x = .PC.x),
                                              y = !!rlang::sym(x = .PC.y),
                                              color = feature,
                                              tooltip = topFeatTab)) +
         ggplot2::theme_bw() +
         ggplot2::theme(legend.position = .legend.position,
                        plot.title = ggplot2::element_text(hjust = 0.5),
                        plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
         ggplot2::labs(title = .plot.title,
                       subtitle = if(.hover.stats != "none"){
                         paste0(.hover.stats,
                                " stats (hover)")
                       } else {
                         ""
                       },
                       color = .color.label)
      )()
    
    # Apply color scale based on feature type
    if(is.numeric(x = dat.df$feature)){
      # Diverging scale for continuous features (e.g., fold changes)
      p <-
        p +
        ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
                                       mid = unname(obj = Color.Palette[1,6]),
                                       high = unname(obj = Color.Palette[1,2]),
                                       midpoint = .midpoint)
    } else {
      # Discrete palette for categorical features (e.g., clusters)
      p <-
        p  +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
        ggplot2::scale_color_manual(
          values = grDevices::colorRampPalette(
            colors = unname(
              obj = .cat.feature.color)
          )(length(x = unique(x = dat.df$feature))))
    }
    
    # Fix panel dimensions
    p <-
      p +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
    # Return static or interactive plot
    if(.hover.stats == "none"){
      
      p +
        ggplot2::geom_point(size = I(x = dat.df$point.size))
      
    } else {
      
      # Add interactivity if ggiraph available
      if (.is_ggiraph_available()) {
        ggiraph::girafe(ggobj = p +
                          ggiraph::geom_point_interactive(size = I(x = dat.df$point.size))) |>
          ggiraph::girafe_options(
            ggiraph::opts_zoom(max = 10,
                               min = 0.25)
          )
      } else {
        warning("ggiraph package not available. Install with install.packages('ggiraph') for interactive hover features.\n",
                "Returning static plot.", 
                call. = FALSE)
        p +
          ggplot2::geom_point(size = I(x = dat.df$point.size))
      }
      
    }
    
    
  }

#' Plot UMAP
#' Plot UMAP
#'
#' Visualizes landmarks in UMAP space with flexible coloring by features, clusters, or statistical 
#' results. Supports both continuous (e.g., gene expression, fold changes) and categorical 
#' (e.g., clusters, cell types) overlays. Optional interactive hover shows landmark feature signatures.
#' Identical interface to \code{plotPCA()} for easy comparison.
#'
#' @param .tdr.obj A tinydenseR object with UMAP computed via \code{get.graph()}.
#' @param .feature Numeric vector or factor of length \code{nrow(.tdr.obj$landmarks)} to color points by. 
#'   Defaults to cluster IDs from \code{$graph$clustering$ids}. Can be:
#'   \itemize{
#'     \item Numeric: gene expression, fold changes, statistics (uses diverging color scale)
#'     \item Factor: clusters, cell types, conditions (uses discrete color palette)
#'   }
#' @param .cat.feature.color Character vector of colors for categorical features. Default uses 
#'   \code{Color.Palette[1,1:5]} (5-color palette, interpolated to match factor levels).
#' @param .panel.size Numeric panel width/height in inches (default 2).
#' @param .midpoint Numeric midpoint for diverging color scale (continuous features only). 
#'   Defaults to median of \code{.feature}.
#' @param .plot.title Character plot title (default "").
#' @param .color.label Character legend title (default "").
#' @param .legend.position Character legend position: "right" (default), "left", "top", "bottom", or "none".
#' @param .point.size Numeric point size (default 0.1).
#' @param .seed Integer random seed for plot point ordering (default 123).
#' @param .hover.stats Character specifying hover information: "none" (default) or "marker" (shows 
#'   landmark feature signatures from \code{get.lm.features.stats()}). Requires \pkg{ggiraph}.
#'   
#' @return Plot object (class depends on interactivity):
#'   \describe{
#'     \item{Static ggplot}{\code{.hover.stats = "none"} returns a \code{ggplot} object}
#'     \item{Interactive girafe}{\code{.hover.stats != "none"} returns a \code{girafe} object (if \pkg{ggiraph} installed), otherwise falls back to static \code{ggplot} with warning}
#'   }
#'   
#' @note Requires \code{get.graph()} to have been run. Interactive hover features require the 
#'   \pkg{ggiraph} package. Install with \code{install.packages("ggiraph")}.
#'   
#' @seealso \code{\link{plotPCA}} for PCA visualization, \code{\link{get.graph}}, 
#'   \code{\link{get.lm.features.stats}} for hover feature computation
#' 
#' @examples
#' \dontrun{
#' # After graph construction
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#' 
#' # Basic UMAP
#' plotUMAP(lm.cells, .point.size = 1)
#' 
#' # Color by gene expression
#' plotUMAP(lm.cells, .feature = lm.cells$landmarks[,"CD4"], .color.label = "CD4")
#' 
#' # Color by fold change from get.lm()
#' condition.stats <- get.lm(lm.cells, .design = design)
#' plotUMAP(lm.cells, 
#'          .feature = condition.stats$fit$coefficients[,"ConditionB"],
#'          .color.label = "log2 FC", 
#'          .midpoint = 0)
#' 
#' # Interactive with feature hover
#' lm.cells <- get.lm.features.stats(lm.cells)
#' plotUMAP(lm.cells, .hover.stats = "marker")
#' }
#' 
#' @export
plotUMAP <-
  function(.tdr.obj,
           .feature = .tdr.obj$graph$clustering$ids,
           .cat.feature.color = Color.Palette[1,1:5],
           .panel.size = 2,
           .midpoint = NULL,
           .plot.title = "",
           .color.label = "",
           .legend.position = "right",
           .point.size = 0.1,
           .seed = 123,
           .hover.stats = "none"){
    
    # R CMD check appeasement
    umap.1 <- umap.2 <- topFeatTab <- feature <- NULL
    
    if(is.null(x = .tdr.obj$graph)){
      stop("Graph component missing. Run get.graph() before plotting UMAP.")
    }
    
    .hover.stats <-
      match.arg(arg = .hover.stats,
                choices = c(
                  "none",
                  "marker"
                ))
    
    # Set midpoint to median for continuous features
    if(is.null(x = .midpoint) &
       is.numeric(x = .feature)){
      .midpoint <- stats::median(x = .feature)
    }
    
    # Build plotting data frame
    set.seed(seed = .seed)
    dat.df <-
      as.data.frame(x = .tdr.obj$graph$uwot$embedding) |>
      cbind(feature = .feature,
            point.size = .point.size)
    
    # Attach hover tooltips if requested
    if(.hover.stats == "marker"){
      
      dat.df$topFeatTab <-
        .tdr.obj$interact.plot$lm.features$html
      
    }
    
    # Randomize row order to avoid plotting bias
    old.seed2 <- .Random.seed
    on.exit(expr = assign(x = ".Random.seed",
                          value = old.seed2,
                          envir = .GlobalEnv),
            add = TRUE)
    set.seed(seed = 123)
    dat.df <-
      dat.df  |>
      (\(x)
       x[nrow(x = x) |>
           seq_len() |>
           sample(),]
      )()
    
    p <-
      dat.df |>
      (\(x)
       ggplot2::ggplot(data = x,
                       mapping = ggplot2::aes(x = umap.1,
                                              y = umap.2,
                                              color = feature,
                                              tooltip = topFeatTab)) +
         ggplot2::theme_bw() +
         ggplot2::theme(legend.position = .legend.position,
                        plot.title = ggplot2::element_text(hjust = 0.5),
                        plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
         ggplot2::labs(title = .plot.title,
                       subtitle = if(.hover.stats != "none"){
                         paste0(.hover.stats,
                                " stats (hover)")
                       } else {
                         ""
                       },
                       color = .color.label)
      )()
    
    if(is.numeric(x = .feature)){
      p <-
        p +
        ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
                                       mid = unname(obj = Color.Palette[1,6]),
                                       high = unname(obj = Color.Palette[1,2]),
                                       midpoint = .midpoint)
    } else {
      p <-
        p  +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
        ggplot2::scale_color_manual(
          values = grDevices::colorRampPalette(
            colors = unname(
              obj = .cat.feature.color)
          )(length(x = unique(x = .feature))))
    }
    
    p <-
      p +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
    if(.hover.stats == "none"){
      
      p +
        ggplot2::geom_point(size = I(x = dat.df$point.size))
      
    } else {
      
      # Check if ggiraph is available for interactive features
      if (.is_ggiraph_available()) {
        ggiraph::girafe(ggobj = p +
                          ggiraph::geom_point_interactive(size = I(x = dat.df$point.size))) |>
          ggiraph::girafe_options(
            ggiraph::opts_zoom(max = 10,
                               min = 0.25)
          )
      } else {
        warning("ggiraph package not available. Interactive hover stats not available. Returning static plot.", 
                call. = FALSE)
        p +
          ggplot2::geom_point(size = I(x = dat.df$point.size))
      }
      
    }
    
  }

#' Bee Swarm Plot of Density Estimate Change
#'
#' Creates a beeswarm plot showing effect sizes (log fold changes) from differential density 
#' testing, with points colored by significance. Optionally splits results by cluster or cell type 
#' and displays mean cell percentages alongside for biological context.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()} and \code{get.lm()}.
#'   Statistical results should be stored in \code{.tdr.obj$map$lm}.
#' @param .model.name Character string naming which model fit to use from \code{.tdr.obj$map$lm}
#'   (default "default"). Must match a name used in \code{get.lm(.model.name = ...)}.
#' @param .coefs Character vector of coefficient names from the design matrix to plot. Must match 
#'   column names in the model fit's coefficients. Can plot single or multiple coefficients.
#' @param .q Numeric q-value threshold for significance coloring (default 0.1). Points with 
#'   q < threshold are colored by direction (red/blue), otherwise gray.
#' @param .q.from Character specifying q-value source: "pca.weighted.q" (default, PCA-variance weighted) 
#'   or "density.weighted.bh.fdr" (density weighted). See \code{get.lm()} for details.
#' @param .split.by Character controlling plot faceting: "none" (default if multiple coefficients), 
#'   "clustering" (split by clusters), or "celltyping" (split by cell types). Default "clustering" 
#'   for single coefficient.
#' @param .swarm.title Character plot title. If NULL and no facets, uses coefficient name.
#' @param .label.substr.rm Character substring to remove from axis labels (default "").
#' @param .point.size Numeric point size (default 0.1).
#' @param .facet.scales Character facet scales: "fixed" (default), "free", "free_x", "free_y".
#' @param .row.space.scaler Numeric scaling factor for row height when splitting (default 0.2 inches per row).
#' @param .col.space.scaler Numeric scaling factor for column width when splitting (default 0.1).
#' @param .panel.width Numeric panel width in inches when plotting multiple coefficients (default 1.5).
#' @param .legend.position Character: "right" (default), "bottom", "left", "top", "none".
#' @param .perc.plot Logical whether to show cell percentage plot alongside beeswarm (default TRUE). 
#'   Only applies when \code{.split.by != "none"}.
#' @param .order.ids Logical whether to order IDs based on dendrogram order (default FALSE).
#'   
#' @return A \code{patchwork} object combining plots:
#'   \describe{
#'     \item{With percentages}{\code{.perc.plot = TRUE} and \code{.split.by != "none"}: percentage plot + beeswarm plot side-by-side}
#'     \item{Beeswarm only}{Otherwise: beeswarm plot alone (still wrapped in patchwork)}
#'   }
#'   
#' @details
#' The beeswarm layout prevents overlapping points, making it easier to assess the distribution 
#' of effect sizes across landmarks. Red indicates significant increases (log FC > 0), blue indicates 
#' significant decreases, and gray indicates non-significant changes.
#' 
#' When split by clustering/celltyping, the percentage plot shows mean cell type abundances across 
#' samples, helping interpret whether changes occur in rare or common populations.
#' 
#' @seealso \code{\link{get.lm}}, \code{\link{plotDensity}}
#' 
#' @examples
#' \dontrun{
#' # After statistical testing
#' lm.obj <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.map()
#' 
#' design <- model.matrix(~ Condition + Replicate, data = .meta)
#' lm.obj <- get.lm(lm.obj, .design = design)
#' 
#' # Basic beeswarm split by clusters
#' plotBeeswarm(lm.obj, .coefs = "ConditionB")
#' 
#' # Multiple coefficients without splitting
#' plotBeeswarm(lm.obj, .coefs = c("ConditionB", "ConditionC"), .split.by = "none")
#' 
#' # Using a different model fit
#' plotBeeswarm(lm.obj, .model.name = "full", .coefs = "ConditionB")
#' }
#' 
#' @export
plotBeeswarm <-
  function(
    .tdr.obj,
    .model.name = "default",
    .coefs,
    .q = 0.1,
    .q.from = "pca.weighted.q",
    .split.by = if(length(x = .coefs) > 1) "none" else "clustering",
    .swarm.title = NULL,
    .label.substr.rm = "",
    .point.size = 0.1,
    .facet.scales = "fixed",
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.1,
    .panel.width = 1.5,
    .legend.position = "right",
    .perc.plot = TRUE,
    .order.ids = FALSE) {
    
    sig <- adj.p <- facets <- pos.t <- neg.t <- q.bars <- dat.df <- value <- split.by <- .coef <- cell.perc <- pop <- NULL
    
    # Validate model exists
    if(is.null(x = .tdr.obj$map$lm[[.model.name]])){
      stop(paste0("Model '", .model.name, "' not found in .tdr.obj$map$lm. ",
                  "Available models: ", paste(names(.tdr.obj$map$lm), collapse = ", ")))
    }
    .stats.obj <- .tdr.obj$map$lm[[.model.name]]
    
    .split.by <-
      match.arg(arg = .split.by,
                choices = c("none",
                            "clustering",
                            "celltyping"))
    
    .legend.position <-
      match.arg(arg = .legend.position,
                choices = c("bottom",
                            "left",
                            "right",
                            "top",
                            "none"))
    
    .q.from <-
      match.arg(arg = .q.from,
                choices = c("pca.weighted.q",
                            "density.weighted.bh.fdr"))
    
    if(.split.by != "none"){
      
      perc.plot <-
        data.frame(x = as.factor(x = 1),
                   pop = colnames(x = .tdr.obj$map[[.split.by]]$cell.perc),
                   cell.perc = Matrix::colMeans(x = .tdr.obj$map[[.split.by]]$cell.perc)) |>
        dplyr::mutate(pop = factor(x = pop,
                                   levels = if(isTRUE(x = .order.ids)){
                                     
                                     .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$labels[
                                       .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$order
                                     ]
                                     
                                   } else {
                                     
                                     pop
                                     
                                   })) |>
        (\(x)
         ggplot2::ggplot(data = x,
                         mapping = ggplot2::aes(x = x,
                                                y = pop,
                                                size = cell.perc)) +
           ggplot2::theme_minimal() +
           ggplot2::theme(legend.margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
                          plot.title = ggplot2::element_text(hjust = 0.5),
                          legend.position = "bottom") +
           ggplot2::labs(title = "percentages",
                         x = "",
                         y = "",
                         size = "%") +
           ggplot2::scale_x_discrete(labels = "mean % across\nall samples") +
           ggplot2::scale_y_discrete(limits = rev) +
           ggplot2::scale_size_continuous(range = log2(x = x$cell.perc + 1) |>
                                            range() |>
                                            I()) +
           ggh4x::force_panelsizes(cols = ggplot2::unit(x = 1.5,
                                                        units = "in"),
                                   rows = ggplot2::unit(x = length(x = x$pop) * .row.space.scaler,
                                                        units = "in")) +
           ggplot2::geom_point()
        )()
      
    }
    
    .facets <-
      ((length(x = .coefs) > 1) &
         .split.by %in% c("clustering",
                          "celltyping"))
    
    if(!(any(.coefs %in%
             colnames(x = .stats.obj$fit$coefficients)))){
      stop(paste0("given your model, .coefs must be one of the following: ",
                  paste(x = colnames(x = .stats.obj$fit$coefficients),
                        collapse = ", ")))
    }
    
    if((length(x = .coefs) > 1)){
      
      dat.df <-
        as.data.frame(x = .stats.obj$fit$coefficients[!(.tdr.obj$graph[[.split.by]]$ids %in% 
                                                          .tdr.obj$map$cl.ct.to.ign),.coefs]) |>
        tidyr::pivot_longer(cols = tidyselect::everything(),
                            names_to = "split.by",
                            values_to = "value",
                            cols_vary = "slowest") |>
        dplyr::bind_cols(sig = as.data.frame(x = (.stats.obj$fit[[.q.from]][!(.tdr.obj$graph[[.split.by]]$ids %in% 
                                                                                .tdr.obj$map$cl.ct.to.ign),.coefs] < .q)) |>
                           tidyr::pivot_longer(cols = tidyselect::everything(),
                                               names_to = "split.by",
                                               values_to = "adj.p",
                                               cols_vary = "slowest") |>
                           dplyr::select(adj.p) |>
                           unlist(use.names = FALSE)) |>
        dplyr::mutate(split.by = factor(x = split.by,
                                        levels = colnames(x = .stats.obj$fit$coefficients)[
                                          match(x = .coefs,
                                                table = colnames(x = .stats.obj$fit$coefficients))] |>
                                          rev()),
                      sig = ifelse(test = sig,
                                   yes = ifelse(test = value > 0,
                                                yes = "more abundant",
                                                no = "less abundant"),
                                   no = "not sig.")) |>
        dplyr::mutate(split.by =
                        gsub(pattern = .label.substr.rm,
                             replacement = "",
                             x = split.by,
                             fixed = FALSE) |>
                        (\(x)
                         factor(x = x,
                                levels = unique(x = x))
                        )())
      
      if(isTRUE(x = .facets)){
        
        dat.df <-
          dplyr::mutate(.data = dat.df,
                        facets = split.by,
                        split.by = rep(x = .tdr.obj$graph[[.split.by]]$ids[!(.tdr.obj$graph[[.split.by]]$ids %in% 
                                                                               .tdr.obj$map$cl.ct.to.ign)] |>
                                         as.character(),
                                       times = length(x = .coefs))  |> 
                          factor(levels = if(isTRUE(x = .order.ids)){
                            
                            .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$labels[
                              .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$order
                            ]
                            
                          } else {
                            
                            levels(x = .tdr.obj$graph[[.split.by]]$ids)
                            
                          }) |>
                          droplevels())
        
      }
      
    } else {
      
      if(.split.by == "celltyping"){
        
        if(is.null(x = .tdr.obj$graph$celltyping$ids)){
          
          stop(".tdr.obj$graph$celltyping$ids could not be found")
          
        }
        
      }
      
      dat.df <-
        data.frame(
          value =
            .stats.obj$fit$coefficients[,.coefs],
          sig =
            ifelse(test = .stats.obj$fit[[.q.from]][,.coefs] < .q,
                   yes = ifelse(test = .stats.obj$fit$coefficients[,.coefs] > 0,
                                yes = "more abundant",
                                no = "less abundant"),
                   no = "not sig."),
          split.by =
            as.character(x = .tdr.obj$graph[[.split.by]]$ids) |> 
            factor(levels = if(isTRUE(x = .order.ids)){
              
              .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$labels[
                .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$order
              ]
              
            } else {
              
              levels(x = .tdr.obj$graph[[.split.by]]$ids)
              
            })
        ) |> 
        (\(x)
         x[!(.tdr.obj$graph[[.split.by]]$ids %in% 
               .tdr.obj$map$cl.ct.to.ign),]
        )() |>
        droplevels()
      
    }
    
    old.seed <- .Random.seed
    on.exit(expr = assign(x = ".Random.seed",
                          value = old.seed,
                          envir = .GlobalEnv),
            add = TRUE)
    set.seed(seed = 123)
    other.plot <-
      dat.df |>
      (\(x)
       ggplot2::ggplot(data = x,
                       mapping = ggplot2::aes(x = value,
                                              y = split.by,
                                              color = sig)) +
         ggplot2::theme_bw() +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.position = .legend.position,
                        axis.title.y = ggplot2::element_blank()) +
         {if(isTRUE(x = .perc.plot)){
           ggplot2::theme(axis.text.y = ggplot2::element_blank())
           
         }} +
         ggplot2::labs(
           title = if(!is.null(x = .swarm.title)) .swarm.title else if(isTRUE(x = .facets)) "" else {
             gsub(pattern = .label.substr.rm,
                  replacement = "",
                  x = .coefs,
                  fixed = FALSE)  
           },
           x = "density\nlog2(+0.5)FC",
           color = paste0("q < ", .q)) +
         {
           if(isTRUE(x = .facets)) {
             ggplot2::facet_grid(cols = ggplot2::vars(facets), 
                                 scales = .facet.scales)
           }
         } +
         ggplot2::scale_y_discrete(limits = rev) +
         ggplot2::geom_point(position = ggplot2::position_jitter(width = 0,
                                                                 height = 0.25),
                             size = I(x = .point.size)) +
         ggplot2::geom_violin(color = "black",
                              alpha = 0,
                              scale = "width",
                              quantiles = 0.5,
                              quantile.linetype = "solid") +
         ggplot2::geom_vline(xintercept = 0,
                             linetype = "dotted",
                             color = "black") +
         ggplot2::scale_color_manual(
           values = stats::setNames(object = Color.Palette[1,c(1,6,2)],
                                    nm = c("less abundant", "not sig.", "more abundant"))) +
         #ggplot2::scale_x_continuous(labels = ~ ifelse(test = .x >= 0,
         #                                              yes = round(x = 2^.x,
         #                                                          digits = 1),
         #                                              no = -round(x = 1 / 2^.x,
         #                                                          digits = 1))) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
         ggh4x::force_panelsizes(
           col = grid::unit(
             x = if(isTRUE(x = .facets)) {
               (as.character(x = dat.df$facets) |>
                  unique() |>
                  nchar() |>
                  max()) * .col.space.scaler
             } else {
               .panel.width
             },
             units = "in"),
           rows = grid::unit(x = (unique(x = dat.df$split.by) |>
                                    length()) * .row.space.scaler,
                             units = "in"))
       
      )()
    
    if(isTRUE(x = .perc.plot)){
      patchwork::wrap_plots(perc.plot + other.plot, 
                            ncol = 2)  
    } else {
      patchwork::wrap_plots(other.plot, 
                            ncol = 1)  
    }
    
  }

#' Bidimensional Hexbin Plot for Marker Expression
#'
#' Creates hexagonal binning plots to visualize the joint distribution of two markers across 
#' landmarks. Useful for exploring marker coexpression patterns and identifying cell populations. 
#' Optionally overlays reference density from all landmarks when plotting a specific cluster/celltype.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.landmarks()} and \code{get.graph()}.
#' @param .id Optional character: cluster or celltype ID to highlight. If NULL (default), plots all landmarks.
#' @param .id.from Character: "clustering" (default) or "celltyping". Source of \code{.id}.
#' @param .x.feature Character: column name from \code{.tdr.obj$landmarks} for x-axis (default "CD3").
#' @param .y.feature Character: column name from \code{.tdr.obj$landmarks} for y-axis (default "CD20").
#' @param .bins Integer: number of hexagonal bins for main plot (default 128). Higher values = finer resolution.
#' @param .legend.position Character: "right" (default), "left", "top", "bottom", or "none".
#' @param .plot.title Character: plot title (default "").
#' @param .panel.size Numeric: panel width/height in inches (default 1.5).
#' @param .reference Logical: if TRUE (default) and \code{.id} is specified, overlay reference density 
#'   contours from all landmarks for comparison.
#' @param .density.bins Integer: number of bins for reference density contours (default 32).
#' @param .sd.range Numeric vector: range of standard deviations for outlier exclusion (default c(-3, 6)). 
#'   Currently not implemented in the function.
#'   
#' @return A \code{ggplot} object with hexagonal binning showing marker coexpression.
#'   
#' @details
#' The function creates a hexagonal heatmap where:
#' \itemize{
#'   \item Each hexagon represents a bin in 2D marker expression space
#'   \item Color intensity (red gradient) indicates cell density (log2 scale)
#'   \item For cytometry data, expression values are divided by 50 for scaling
#'   \item When \code{.id} is specified and \code{.reference = TRUE}, gray density contours 
#'     show the distribution of all landmarks for context
#' }
#' 
#' This visualization helps identify:
#' \itemize{
#'   \item Marker coexpression patterns (e.g., CD3+CD4+ vs CD3+CD8+ populations)
#'   \item Whether a cluster is truly distinct from the background
#'   \item Outlier populations in marker expression space
#' }
#' 
#' @seealso \code{\link{plotPCA}}, \code{\link{plotUMAP}} for dimensionality reduction visualizations
#' 
#' @examples
#' \dontrun{
#' # After landmark identification and graph construction
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#' 
#' # All landmarks
#' plot2Markers(lm.cells, .x.feature = "CD3", .y.feature = "CD4")
#' 
#' # Specific cluster with reference density
#' plot2Markers(lm.cells, 
#'              .id = "1", 
#'              .x.feature = "CD3", 
#'              .y.feature = "CD8",
#'              .reference = TRUE)
#' }
#' 
#' @export
plot2Markers <-
  function(.tdr.obj,
           #           .color.from = NULL,
           .id = NULL,
           .id.from = "clustering",
           .x.feature = "CD3",
           .y.feature = "CD20",
           .bins = 128,
           .legend.position = "right",
           .plot.title = "",
           .panel.size = 1.5,
           .reference = TRUE,
           .density.bins = 32,
           .sd.range = c(-3, 6)){
    
    level <- NULL
    
    if(!is.null(x = .id)){
      
      .id.from <-
        match.arg(arg = .id.from,
                  choices = c("clustering",
                              "celltyping"))
      
      if(.id.from == "clustering"){
        
        .id <- .tdr.obj$graph$clustering$ids == .id
        
      } else if(.id.from == "celltyping"){
        
        if(is.null(x = .tdr.obj$graph$celltyping$ids)){
          stop(".tdr.obj$graph$celltyping$ids could not be found")
        }
        .id <- .tdr.obj$graph$celltyping$ids == .id
      }
    }
    
    if(.tdr.obj$config$assay.type == "RNA"){
      
      dat.df <-
        (.tdr.obj$landmarks[if(!is.null(x = .id)) .id else 1:nrow(x = .tdr.obj$landmarks),
                            colnames(x = .tdr.obj$landmarks) %in% c(.x.feature,.y.feature)]) |>
        as.data.frame()
      
      
    } else {
      
      dat.df <-
        (.tdr.obj$landmarks[if(!is.null(x = .id)) .id else 1:nrow(x = .tdr.obj$landmarks),
                            colnames(x = .tdr.obj$landmarks) %in% c(.x.feature,.y.feature)]) |>
        as.data.frame()
      
    }
    
    p <-
      ggplot2::ggplot(data = dat.df,
                      mapping = ggplot2::aes(x = !!rlang::sym(x = .x.feature),
                                             y = !!rlang::sym(x = .y.feature))) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = .legend.position,
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = .plot.title)
    
    if(.tdr.obj$config$assay.type != "RNA"){
      
      p <-
        p +
        ggplot2::scale_x_continuous(limits = (.tdr.obj$landmarks[,colnames(x = .tdr.obj$landmarks) %in%
                                                                   .x.feature, drop = TRUE]) |>
                                      (\(x)
                                       range(x[!((x > (mean(x = x) + .sd.range[2]*stats::sd(x = x))) |
                                                   (x < (mean(x = x) + .sd.range[1]*stats::sd(x = x))))])
                                      )()) +
        ggplot2::scale_y_continuous(limits = (.tdr.obj$landmarks[,colnames(x = .tdr.obj$landmarks) %in%
                                                                   .y.feature, drop = TRUE]) |>
                                      (\(x)
                                       range(x[!((x > (mean(x = x) + .sd.range[2]*stats::sd(x = x))) |
                                                   (x < (mean(x = x) + .sd.range[1]*stats::sd(x = x))))])
                                      )())
      
    }
    
    if(!is.null(x = .id) &
       isTRUE(x = .reference)) {
      
      p <-
        p  +
        ggplot2::stat_density_2d(data = as.data.frame(x = .tdr.obj$landmarks),
                                 mapping = ggplot2::aes(fill = ggplot2::after_stat(x = log2(x = level))),
                                 geom = "polygon",
                                 bins = .density.bins) +
        ggplot2::scale_fill_gradient(low = unname(obj = Color.Palette[2,7]),
                                     high = unname(obj = Color.Palette[6,7]),
                                     guide = "none")
      
    }
    
    p <-
      p +
      ggnewscale::new_scale(new_aes = "fill") +
      ggplot2::geom_hex(bins = .bins) +
      ggplot2::scale_fill_gradient(low = unname(obj = Color.Palette[1,1]),
                                   high = unname(obj = Color.Palette[1,2]),
                                   trans = "log2")
    
    
    p +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
  }

#' Plot Sample PCA (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' 
#' This function is deprecated. Please use \code{\link{plotSampleEmbedding}} with 
#' \code{.embedding = "pca"} instead.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()} and \code{get.embedding()}.
#' @param .labels.from Character specifying metadata column for coloring points.
#' @param .cat.feature.color Character vector of colors for categorical labels.
#' @param .point.size Numeric point size (default 1).
#' @param .panel.size Numeric panel width/height in inches (default 2).
#' @param .midpoint Numeric midpoint for diverging color scale (continuous labels only).
#'   
#' @return A \code{ggplot} object showing PC1 vs PC2 of sample-level density profiles.
#'
#' @seealso \code{\link{plotSampleEmbedding}} for the recommended replacement
#' 
#' @export
#'
plotSamplePCA <-
  function(
    .tdr.obj,
    .labels.from = colnames(x = .tdr.obj$metadata)[1],
    .cat.feature.color = Color.Palette[1,1:5],
    .point.size = 1,
    .panel.size = 2,
    .midpoint = if(is.numeric(x = .tdr.obj$metadata[[.labels.from]])) stats::median(x = .tdr.obj$metadata[[.labels.from]]) else NA
  ){
    
    .Deprecated(new = "plotSampleEmbedding", 
                msg = "plotSamplePCA is deprecated. Use plotSampleEmbedding(.embedding = 'pca') instead.")
    
    plotSampleEmbedding(
      .tdr.obj = .tdr.obj,
      .embedding = "pca",
      .color.by = .labels.from,
      .cat.feature.color = .cat.feature.color,
      .point.size = .point.size,
      .panel.size = .panel.size,
      .midpoint = .midpoint
    )
    
  }


#' Plot Sample Embedding from get.embedding
#'
#' Visualizes sample-level embeddings computed by \code{\link{get.embedding}}. Supports
#' three embedding types: supervised partial-effect PCA (pePC), unsupervised PCA on
#' landmark densities (pca), and diffusion map trajectory (traj).
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()} and 
#'   \code{get.embedding()}. All embeddings are stored in \code{.tdr.obj$map$embedding}:
#'   unsupervised (\code{pca}, \code{traj}) and supervised (\code{pePC$<name>}).
#' @param .embedding Character string specifying which embedding type to plot. One of:
#'   \describe{
#'     \item{"pePC"}{(default) Supervised partial-effect PCA. Requires \code{.sup.embed.slot}.}
#'     \item{"pca"}{Unsupervised PCA on log-transformed landmark densities. Uses \code{.tdr.obj$map$embedding$pca}.}
#'     \item{"traj"}{Diffusion map trajectory embedding. Uses \code{.tdr.obj$map$embedding$traj}.}
#'   }
#' @param .sup.embed.slot Character string specifying which supervised embedding to plot.
#'   Only used when \code{.embedding = "pePC"}. Must match a name in 
#'   \code{.tdr.obj$map$embedding$pePC} (contrast name for FWL method, term name for nested model).
#'   Ignored (with warning) when \code{.embedding} is "pca" or "traj".
#' @param .color.by Character string specifying which metadata column to use for coloring
#'   points. Defaults to \code{.sup.embed.slot} if it exists in metadata (for pePC),
#'   otherwise uses the first metadata column.
#' @param .x.by Character string specifying which metadata column to use for the x-axis
#'   in single-PC embeddings. Defaults to \code{NULL}, which uses \code{.sup.embed.slot}
#'   if it exists in metadata. Required when \code{.sup.embed.slot} is not a metadata
#'   column (e.g., FWL contrast names).
#' @param .pc.x Integer specifying which PC/DC to plot on x-axis (default 1). Only used
#'   when embedding has more than 1 dimension.
#' @param .pc.y Integer specifying which PC/DC to plot on y-axis (default 2). Only used
#'   when embedding has more than 1 dimension.
#' @param .cat.feature.color Character vector of colors for categorical features.
#'   Defaults to first 5 colors from \code{Color.Palette} row 1.
#' @param .point.size Numeric point size (default 1).
#' @param .panel.size Numeric panel size in inches (default 2).
#' @param .midpoint Numeric midpoint for continuous color scale.
#'   Defaults to median of \code{.color.by} column.
#'
#' @return A \code{ggplot} object showing the sample embedding colored by the
#'   metadata variable specified by \code{.color.by}.
#'
#' @details
#' This function plots the sample coordinates from embeddings computed by
#' \code{\link{get.embedding}}:
#'
#' \describe{
#'   \item{pePC (supervised)}{Partial-effect PCA isolating variation attributable to a
#'     specific contrast or set of covariates. Requires running \code{get.embedding()} with
#'     \code{.contrast.of.interest} or \code{.red.model}. Stored in \code{.tdr.obj$map$embedding$pePC$<name>}.}
#'   \item{pca (unsupervised)}{Standard PCA on log-transformed landmark densities. Good for
#'     exploratory visualization and QC. Stored in \code{.tdr.obj$map$embedding$pca}.}
#'   \item{traj (unsupervised)}{Diffusion map trajectory capturing continuous variation.
#'     Useful for developmental or temporal trajectories. Stored in \code{.tdr.obj$map$embedding$traj}.}
#' }
#'
#' For pePC embeddings with rank > 1 (e.g., from nested model comparison with multiple
#' covariates), you can plot different PC combinations using \code{.pc.x} and \code{.pc.y}.
#'
#' When pePC embedding has only 1 PC (e.g., from a single contrast), the function
#' creates either a scatter plot (for continuous x-axis variable) or a boxplot 
#' (for categorical x-axis variable), with the x-axis determined by \code{.x.by}
#' or \code{.sup.embed.slot}.
#'
#' @seealso \code{\link{get.embedding}} for computing embeddings,
#'   \code{\link{plotSamplePCA}} for legacy PCA visualization
#'
#' @examples
#' \dontrun{
#' # Compute unsupervised embeddings
#' lm.obj <- get.embedding(.tdr.obj = lm.obj)
#'
#' # Plot unsupervised PCA
#' plotSampleEmbedding(lm.obj, .embedding = "pca", .color.by = "Condition")
#'
#' # Plot diffusion map trajectory
#' plotSampleEmbedding(lm.obj, .embedding = "traj", .color.by = "Timepoint")
#'
#' # Compute supervised embedding
#' lm.obj <- get.embedding(lm.obj, .contrast.of.interest = "TrtVsCtrl")
#'
#' # Plot supervised pePC embedding
#' plotSampleEmbedding(lm.obj, .embedding = "pePC",
#'                     .sup.embed.slot = "TrtVsCtrl", .color.by = "Group")
#' }
#'
#' @export
#'
plotSampleEmbedding <-
  function(
    .tdr.obj,
    .embedding = "pePC",
    .sup.embed.slot = NULL,
    .color.by = NULL,
    .x.by = NULL,
    .pc.x = 1,
    .pc.y = 2,
    .cat.feature.color = Color.Palette[1,1:5],
    .point.size = 1,
    .panel.size = 2,
    .midpoint = NULL
  ){
    
    PC_x <- PC_y <- PC1 <- x.var <- color.var <- labels.from <- NULL
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    .embedding <- 
      match.arg(arg = .embedding,
                choices = c("pePC", "pca", "traj"))
    
    # For pePC, require .sup.embed.slot to exist in .tdr.obj$map$embedding$pePC
    if(.embedding == "pePC"){
      
      if(is.null(x = .sup.embed.slot) || 
         !is.character(x = .sup.embed.slot) || 
         length(x = .sup.embed.slot) != 1){
        stop("'.sup.embed.slot' must be a single character string when .embedding = 'pePC'")
      }
      
      if(is.null(x = .tdr.obj$map$embedding$pePC[[.sup.embed.slot]])){
        available <- if(!is.null(.tdr.obj$map$embedding$pePC)) names(.tdr.obj$map$embedding$pePC) else "none"
        stop(paste0("'", .sup.embed.slot, "' not found in .tdr.obj$map$embedding$pePC. ",
                    "Run get.embedding() with .contrast.of.interest or .red.model first. ",
                    "Available: ", paste(available, collapse = ", ")))
      }
      
    } else {
      
      # For pca/traj, .sup.embed.slot should not be provided
      if(!is.null(x = .sup.embed.slot)){
        warning("'.sup.embed.slot' is ignored when .embedding = '", .embedding, 
                "'. It is only used for supervised (pePC) embeddings.",
                call. = FALSE)
      }
      
      # Validate the corresponding slot exists in .tdr.obj$map$embedding
      if(is.null(x = .tdr.obj$map$embedding[[.embedding]])){
        stop(paste0("'", .embedding, "' embedding not found in .tdr.obj$map$embedding. Run get.embedding() first."))
      }
      
    }
    
    # Determine color.by: use provided value, fall back to .sup.embed.slot if in metadata (for pePC), else first column
    if(is.null(x = .color.by)){
      if(.embedding == "pePC" && 
         !is.null(x = .sup.embed.slot) && 
         .sup.embed.slot %in% colnames(x = .tdr.obj$metadata)){
        .color.by <- .sup.embed.slot
      } else {
        .color.by <- colnames(x = .tdr.obj$metadata)[1]
        if(.embedding == "pePC"){
          message("'.sup.embed.slot' not in metadata. Coloring by '", .color.by, "'")
        }
      }
    }
    
    if(!(.color.by %in% colnames(x = .tdr.obj$metadata))){
      stop(paste0("'.color.by' = '", .color.by, "' not found in metadata columns: ",
                  paste(colnames(.tdr.obj$metadata), collapse = ", ")))
    }
    
    # Validate .x.by if provided
    if(!is.null(x = .x.by)){
      if(!(.x.by %in% colnames(x = .tdr.obj$metadata))){
        stop(paste0("'.x.by' = '", .x.by, "' not found in metadata columns: ",
                    paste(colnames(.tdr.obj$metadata), collapse = ", ")))
      }
    }
    
    # Set midpoint default
    if(is.null(x = .midpoint) && is.numeric(x = .tdr.obj$metadata[[.color.by]])){
      .midpoint <- stats::median(x = .tdr.obj$metadata[[.color.by]])
    }
    
    # -------------------------------------------------------------------------
    # Extract embedding coordinates based on .embedding type
    # -------------------------------------------------------------------------
    
    var.explained <- NULL
    
    if(.embedding == "pePC"){
      # Supervised embeddings are now in .tdr.obj$map$embedding$pePC[[.sup.embed.slot]]
      embed <- .tdr.obj$map$embedding$pePC[[.sup.embed.slot]]
      axis.prefix <- "pePC"
      # pePC stores perc.tot.var.exp directly
      if(!is.null(x = embed$perc.tot.var.exp)){
        var.explained <- embed$perc.tot.var.exp
      }
    } else if(.embedding == "pca"){
      # Unsupervised embeddings are in .tdr.obj$map$embedding
      embed <- .tdr.obj$map$embedding$pca
      axis.prefix <- "PC"
      # PCA stores sdev, compute variance explained as proportion of TOTAL variance
      # (not sum of truncated sdev^2, since we use truncated PCA)
      if(!is.null(x = embed$sdev)){
        total.var <- matrixStats::rowVars(x = .tdr.obj$map$Y) |> sum()
        var.explained <- 100 * embed$sdev^2 / total.var
        names(x = var.explained) <- paste0("PC", seq_along(along.with = var.explained))
      }
    } else { # traj
      embed <- .tdr.obj$map$embedding$traj
      axis.prefix <- "DC"
      # Diffusion map eigenvalues don't represent variance explained in same way
    }
    
    if(is.null(x = embed$coord)){
      stop("Embedding does not contain coordinates. Check get.embedding() output.")
    }
    
    n.pcs <- ncol(x = embed$coord)
    
    # -------------------------------------------------------------------------
    # Handle single PC case: scatter for continuous x, boxplot for categorical x
    # (only applies to pePC embeddings)
    # -------------------------------------------------------------------------
    
    if(n.pcs == 1 && .embedding == "pePC"){
      
      # Determine x-axis variable
      if(!is.null(x = .x.by)){
        x.by.col <- .x.by
      } else if(!is.null(x = .sup.embed.slot) && 
                .sup.embed.slot %in% colnames(x = .tdr.obj$metadata)){
        x.by.col <- .sup.embed.slot
      } else {
        stop(paste0("'.sup.embed.slot' = '", .sup.embed.slot, 
                    "' is not a metadata column. Please specify '.x.by' as one of: ",
                    paste(colnames(.tdr.obj$metadata), collapse = ", ")))
      }
      
      plot.data <-
        data.frame(
          PC1 = embed$coord[, 1],
          x.var = .tdr.obj$metadata[[x.by.col]],
          color.var = .tdr.obj$metadata[[.color.by]]
        )
      
      # Build y-axis label with variance explained if available
      y.label <- "pePC1"
      if(!is.null(x = var.explained) && length(x = var.explained) >= 1){
        y.label <- paste0("pePC1 (", round(x = var.explained[1], digits = 1), "%)")
      }
      
      # Set midpoint for color scale if x.by differs from color.by
      if(is.null(x = .midpoint) && is.numeric(x = plot.data$color.var)){
        .midpoint <- stats::median(x = plot.data$color.var)
      }
      
      if(is.numeric(x = plot.data$x.var)){
        
        # ----- Scatter plot for continuous x-axis -----
        p <-
          ggplot2::ggplot(data = plot.data,
                          mapping = ggplot2::aes(x = x.var,
                                                 y = PC1,
                                                 color = color.var)) +
          ggplot2::theme_bw() +
          ggplot2::geom_point(size = I(x = .point.size)) +
          ggplot2::labs(
            x = x.by.col,
            y = y.label,
            color = .color.by
          ) +
          ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                    units = "in"),
                                  col = grid::unit(x = .panel.size,
                                                   units = "in"))
        
      } else {
        
        # ----- Boxplot for categorical x-axis -----
        p <-
          ggplot2::ggplot(data = plot.data,
                          mapping = ggplot2::aes(x = x.var,
                                                 y = PC1)) +
          ggplot2::theme_bw() +
          ggplot2::geom_boxplot(outliers = FALSE) +
          ggplot2::geom_point(mapping = ggplot2::aes(color = color.var),
                              size = I(x = .point.size),
                              position = ggplot2::position_jitter(height = 0,
                                                                  width = 0.2,
                                                                  seed = 123)) +
          ggplot2::labs(
            x = x.by.col,
            y = y.label,
            color = .color.by
          ) +
          ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                    units = "in"),
                                  col = grid::unit(x = .panel.size,
                                                   units = "in"))
        
      }
      
      # Apply color scale for points
      if(is.numeric(x = plot.data$color.var)){
        p <-
          p +
          ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
                                         mid = unname(obj = Color.Palette[1,6]),
                                         high = unname(obj = Color.Palette[1,2]),
                                         midpoint = .midpoint)
      } else {
        p <-
          p +
          ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
          ggplot2::scale_color_manual(
            values = grDevices::colorRampPalette(
              colors = unname(
                obj = .cat.feature.color)
            )(length(x = unique(x = plot.data$color.var))))
      }
      
      return(p)
      
    }
    
    # -------------------------------------------------------------------------
    # Multi-PC/DC case: validate requested components
    # -------------------------------------------------------------------------
    
    if(.pc.x > n.pcs || .pc.y > n.pcs){
      stop(paste0("Requested component (", max(.pc.x, .pc.y), ") exceeds available (", n.pcs, ")"))
    }
    
    # -------------------------------------------------------------------------
    # Build axis labels with variance explained if available
    # -------------------------------------------------------------------------
    
    if(!is.null(x = var.explained) && 
       length(x = var.explained) >= max(.pc.x, .pc.y)){
      x.label <- paste0(axis.prefix, .pc.x, " (", round(x = var.explained[.pc.x], digits = 1), "%)")
      y.label <- paste0(axis.prefix, .pc.y, " (", round(x = var.explained[.pc.y], digits = 1), "%)")
    } else {
      x.label <- paste0(axis.prefix, .pc.x)
      y.label <- paste0(axis.prefix, .pc.y)
    }
    
    # -------------------------------------------------------------------------
    # Build plot data
    # -------------------------------------------------------------------------
    
    plot.data <-
      data.frame(
        PC_x = embed$coord[, .pc.x],
        PC_y = embed$coord[, .pc.y],
        labels.from = .tdr.obj$metadata[[.color.by]]
      )
    
    # -------------------------------------------------------------------------
    # Create plot
    # -------------------------------------------------------------------------
    
    p <-
      ggplot2::ggplot(data = plot.data,
                      mapping = ggplot2::aes(x = PC_x,
                                             y = PC_y,
                                             color = labels.from)) +
      ggplot2::theme_bw() +
      ggplot2::geom_point(size = I(x = .point.size)) +
      ggplot2::labs(
        x = x.label,
        y = y.label,
        color = .color.by
      ) +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
    # -------------------------------------------------------------------------
    # Apply color scale
    # -------------------------------------------------------------------------
    
    if(is.numeric(x = plot.data$labels.from)){
      p <-
        p +
        ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
                                       mid = unname(obj = Color.Palette[1,6]),
                                       high = unname(obj = Color.Palette[1,2]),
                                       midpoint = .midpoint)
    } else {
      p <-
        p +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
        ggplot2::scale_color_manual(
          values = grDevices::colorRampPalette(
            colors = unname(
              obj = .cat.feature.color)
          )(length(x = unique(x = plot.data$labels.from))))
    }
    
    return(p)
    
  }


#' Plot Traditional Statistics
#'
#' Visualizes results from traditional cluster/cell type-level differential abundance testing. 
#' Shows effect sizes as heatmap with significance markers, providing a complementary view to 
#' landmark-based analysis for easier interpretation at the population level.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()} and \code{get.lm()}.
#'   Statistical results should be stored in \code{.tdr.obj$map$lm}.
#' @param .model.name Character string naming which model fit to use from \code{.tdr.obj$map$lm}
#'   (default "default"). Must match a name used in \code{get.lm(.model.name = ...)}.
#' @param .split.by Character: "clustering" (default) or "celltyping" - which population grouping to use.
#' @param .coefs Character vector of coefficient names to plot. Defaults to all coefficients from 
#'   traditional model.
#' @param .q Numeric q-value threshold for significance stars (default 0.1).
#' @param .row.space.scaler Numeric scaling for row height (default 0.2 inches per population).
#' @param .col.space.scaler Numeric scaling for column width (default 0.5 inches per coefficient).
#' @param .label.substr.rm Character substring to remove from labels (default "").
#' @param .order.ids Logical whether to order IDs based on dendrogram order (default FALSE).
#'   
#' @return A \code{patchwork} object combining two plots: (1) a dot plot showing mean cell 
#'   percentages per population, and (2) a heatmap showing log fold changes colored by magnitude 
#'   with significance indicated by asterisks.
#'   
#' @details
#' Traditional analysis tests for DA at the cluster/celltype level by aggregating cell counts 
#' per sample, then using standard models (typically edgeR). This provides:
#' \itemize{
#'   \item Easier biological interpretation (known populations vs abstract landmarks)
#'   \item Comparison to landmark-based results for validation
#'   \item Detection of broad shifts affecting entire populations
#' }
#' 
#' Note: Traditional analysis is less sensitive to subtle within-cluster variation that 
#' landmark-based methods can detect.
#' 
#' @seealso \code{\link{get.lm}}, \code{\link{plotTradPerc}}, \code{\link{plotBeeswarm}}
#' 
#' @examples
#' \dontrun{
#' # After get.lm with traditional analysis
#' lm.obj <- get.lm(lm.obj, .design = design)
#' 
#' # Heatmap of cluster-level changes
#' plotTradStats(lm.obj, .split.by = "clustering")
#' 
#' # Cell type-level changes (using a different model)
#' plotTradStats(lm.obj, .model.name = "full", .split.by = "celltyping")
#' }
#' 
#' @export
#'
plotTradStats <-
  function(
    .tdr.obj,
    .model.name = "default",
    .split.by = "clustering",
    .coefs = NULL,
    .q = 0.1,
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.07,
    .label.substr.rm = "",
    .order.ids = FALSE
  ){
    
    sig.shape <- sig <- adj.p <- level <- pop <- term <- coef <- cell.perc <- NULL
    
    # Validate model exists
    if(is.null(x = .tdr.obj$map$lm[[.model.name]])){
      stop(paste0("Model '", .model.name, "' not found in .tdr.obj$map$lm. ",
                  "Available models: ", paste(names(.tdr.obj$map$lm), collapse = ", ")))
    }
    .stats.obj <- .tdr.obj$map$lm[[.model.name]]
    
    .split.by <-
      match.arg(arg = .split.by,
                choices = c("clustering",
                            "celltyping"))
    
    # Default .coefs to all coefficients if not provided
    if(is.null(x = .coefs)){
      .coefs <- colnames(x = .stats.obj$trad[[.split.by]]$fit$coefficients)
    }
    
    perc.plot <-
      data.frame(x = as.factor(x = 1),
                 pop = colnames(x = .tdr.obj$map[[.split.by]]$cell.perc),
                 cell.perc = Matrix::colMeans(x = .tdr.obj$map[[.split.by]]$cell.perc)) |>
      dplyr::mutate(pop = factor(x = pop,
                                 levels = if(isTRUE(x = .order.ids)){
                                   
                                   .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$labels[
                                     .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$order
                                   ]
                                   
                                 } else {
                                   
                                   pop
                                   
                                 })) |>
      droplevels() |>
      (\(x)
       ggplot2::ggplot(data = x,
                       mapping = ggplot2::aes(x = x,
                                              y = pop,
                                              size = cell.perc)) +
         ggplot2::theme_minimal() +
         ggplot2::theme(legend.margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
                        plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.position = "bottom") +
         ggplot2::labs(title = "percentages",
                       x = "",
                       y = "",
                       size = "%") +
         ggplot2::scale_x_discrete(labels = "mean % across\nall samples") +
         ggplot2::scale_y_discrete(limits = rev) +
         ggplot2::scale_size_continuous(range = log2(x = x$cell.perc + 1) |>
                                          range()) +
         ggplot2::geom_point() +
         ggh4x::force_panelsizes(cols = ggplot2::unit(x = 1.5,
                                                      units = "in"),
                                 rows = ggplot2::unit(x = pmax(length(x = x$pop) * .row.space.scaler,
                                                               3.5),
                                                      units = "in")) 
      )()
    
    coef.df <-
      cbind(.stats.obj$trad[[.split.by]]$fit$coefficients[,.coefs,drop = FALSE]) |>
      dplyr::as_tibble(rownames = "pop") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "coef")
    
    adj.p.df <-
      cbind(.stats.obj$trad[[.split.by]]$fit$adj.p[,.coefs,drop = FALSE])  |>
      dplyr::as_tibble(rownames = "pop") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "adj.p")
    
    dat.df <-
      dplyr::full_join(x = coef.df,
                       y = adj.p.df,
                       by = c("pop","term")) |>
      dplyr::mutate(sig =
                      ifelse(
                        test = coef < 0,
                        yes = "fewer cells",
                        no = "more cells") |>
                      ifelse(
                        test = adj.p < .q,
                        no = "not sig.")  |>
                      factor(levels = c("fewer cells",
                                        "not sig.",
                                        "more cells"))) |>
      dplyr::mutate(sig.shape = ifelse(
        test = sig == "not sig.",
        yes = NA,
        no = 1))
    
    dat.df$term <-
      gsub(pattern = .label.substr.rm,
           replacement = "",
           x = dat.df$term,
           fixed = FALSE) |>
      (\(x)
       factor(x = x,
              levels = unique(x = x))
      )()
    
    dat.df <-
      dat.df |>
      dplyr::mutate(pop = factor(x = pop,
                                 levels = if(isTRUE(x = .order.ids)){
                                   
                                   .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$labels[
                                     .tdr.obj$graph[[.split.by]]$pheatmap$tree_row$order
                                   ]
                                   
                                 } else {
                                   
                                   colnames(x = .tdr.obj$map[[.split.by]]$cell.perc)
                                   
                                 })) |>
      droplevels()
    
    other.plot <-
      ggplot2::ggplot(data = dat.df,
                      mapping = ggplot2::aes(x = term,
                                             y = pop,
                                             fill = coef,
                                             stroke = sig.shape
                      )) +
      ggplot2::guides(stroke = ggplot2::guide_legend(override.aes = list(size = I(x = 5)),
                                                     order = 1),
                      size = ggplot2::guide_legend(override.aes = list(fill = unname(obj = Color.Palette[1,6]),
                                                                       stroke = NA),
                                                   order = 2),
                      fill = ggplot2::guide_colorbar(order = 3)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),
                     axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5),
                     legend.position = "right",
                     axis.text.y = ggplot2::element_blank()) +
      ggplot2::scale_fill_gradient2(low = unname(obj = Color.Palette[1,1]),
                                    mid = unname(obj = Color.Palette[1,6]),
                                    high = unname(obj = Color.Palette[1,2]),
                                    midpoint = 0
      ) +
      ggplot2::scale_size_continuous(labels = ~ formatC(x = 10^(-(.x)),
                                                        format = "g",
                                                        digits = 1),
                                     range = I(x = c(3,6))) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::geom_point(shape = 21,
                          mapping = ggplot2::aes(size = -log10(x = adj.p)),
                          color = "black",
                          show.legend = TRUE) +
      {if(!(is.na(x = dat.df$sig.shape) |> all())){
        ggplot2::continuous_scale(aesthetics = "stroke",
                                  name =  paste0("q < ",
                                                 .q),
                                  palette = function(x){scales::rescale(x = x, c(0, 1))},
                                  labels = "TRUE")
      }} +
      {if(!(is.na(x = dat.df$sig.shape) |> all())){
        ggplot2::geom_point(shape = 8,
                            mapping = ggplot2::aes(size = -log10(x = adj.p)),
                            color = "black",
                            show.legend = TRUE)
      }} +
      ggplot2::labs(title = "difference in cell %",
                    x = "",
                    y = "",
                    fill = "log2(+0.5)FC",
                    size = if((is.na(x = dat.df$sig.shape) |> all())) {
                      paste0("q\n(",
                             "all > ",
                             .q,
                             ")")} else "q") +
      ggh4x::force_panelsizes(
        col = grid::unit(x = pmax(
          (as.character(x = dat.df$term) |> unique() |> nchar() |> max(na.rm = TRUE)) *
            (unique(x = dat.df$term) |> length()) *
            .col.space.scaler,
          0.5, 
          na.rm = TRUE),
          units = "in"),
        rows = grid::unit(x = pmax((unique(x = dat.df$pop) |> length()) * .row.space.scaler,
                                   3.5, 
                                   na.rm = TRUE),
                          units = "in"))
    
    patchwork::wrap_plots(perc.plot + other.plot, 
                          ncol = 2)
    
  }

#' Plot Traditional Percentages
#'
#' Creates dot/line plots showing cell percentages per sample for specific populations. Useful 
#' for visually inspecting distribution of cell abundances across conditions and identifying 
#' paired/longitudinal patterns. Complements statistical test results from \code{plotTradStats()}.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .x.split Character specifying metadata column for x-axis grouping. Defaults to first 
#'   column.
#' @param .x.split.subset Optional character vector to subset \code{.x.split} categories. Default NULL.
#' @param .pop Character vector of population names to plot. If NULL, plots all populations from 
#'   \code{.pop.from}.
#' @param .pop.from Character: "clustering" (default) or "celltyping" - which grouping to plot.
#' @param .order.pop Logical whether to order populations based on dendrogram order (default FALSE).
#' @param .line.by Character metadata column for connecting paired samples with lines (e.g., 
#'   "Subject" for longitudinal data). Default NULL (no lines).
#' @param .dodge.by Character metadata column for coloring/dodging points. Default NULL (all black).
#' @param .x.space.scaler Numeric scaling factor for x-axis panel width (default 0.25 inches per group).
#' @param .height Numeric plot height in inches (default 1.5).
#' @param .cat.feature.color Character vector of colors for \code{.dodge.by} categories (default 
#'   \code{Color.Palette[1,1:5]}).
#' @param .seed Integer random seed for x-axis jitter (default 123).
#' @param .orientation Character: "wide" (default, all populations in one row) or "square" (facet grid).
#' @param .log2.y Logical whether to log2-transform y-axis percentages (default FALSE).
#'   
#' @return A \code{ggplot} object showing cell percentages with optional paired connections.
#'   
#' @details
#' This function visualizes the raw cell percentages used in traditional DA testing. Points show 
#' individual samples, and lines (if \code{.line.by} specified) connect repeated measures from 
#' the same subject/mouse. Useful for:
#' \itemize{
#'   \item Inspecting data distribution before statistical testing
#'   \item Identifying outliers or batch effects
#'   \item Visualizing paired/longitudinal designs
#'   \item Confirming significant results have biological meaning
#' }
#' 
#' @seealso \code{\link{plotTradStats}}, \code{\link{get.lm}}
#' 
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.map()
#' 
#' # Plot CD4 T cells across conditions
#' plotTradPerc(lm.cells, .pop = "CD4.T.cells", .x.split = "Condition")
#' 
#' # Paired design with subject lines
#' plotTradPerc(lm.cells, 
#'              .pop = "CD4.T.cells",
#'              .line.by = "Subject",
#'              .dodge.by = "Timepoint")
#' }
#' 
#' @export
#'
plotTradPerc <-
  function(
    .tdr.obj,
    .x.split = colnames(x = .tdr.obj$metadata)[1],
    .x.split.subset = NULL,
    .pop = NULL,
    .pop.from = "clustering",
    .order.pop = FALSE,
    .line.by = NULL,
    .dodge.by = NULL,
    .x.space.scaler = 0.25,
    .height = 1.5,
    .cat.feature.color = Color.Palette[1,1:5],
    .seed = 123,
    .orientation = "wide",
    .log2.y = FALSE
  ){
    
    dodge <- value <- name <- color <- group <- x <- y <- NULL
    
    if(length(x = .x.split) != 1){
      stop(".x.split must be length 1")
    }
    
    if(!(.x.split %in% colnames(x = .tdr.obj$metadata))){
      stop(paste0(".x.split must be one of the following: ",
                  paste(x = colnames(x = .tdr.obj$metadata),
                        collapse = ", ")))
    }
    
    if(!is.null(x = .x.split.subset)){
      if(!all(.x.split.subset %in% .tdr.obj$metadata[[.x.split]])){
        stop(paste0(".x.split.subset must within the following: ",
                    paste(x = unique(x = .tdr.obj$metadata[[.x.split]]),
                          collapse = ", ")))
      }
    }
    
    if(!is.null(x = .dodge.by)){
      if(length(x = .dodge.by) != 1){
        stop(".dodge.by must be length 1")
      }
      
      if(!(.dodge.by %in% colnames(x = .tdr.obj$metadata))){
        stop(paste0(".dodge.by must be one of the following: ",
                    paste(x = colnames(x = .tdr.obj$metadata),
                          collapse = ", ")))
      }
      
    }
    
    if(!is.null(x = .line.by)){
      if(length(x = .line.by) != 1){
        stop(".line.by must be length 1")
      }
      
      if(!(.line.by %in% colnames(x = .tdr.obj$metadata))){
        stop(paste0(".line.by must be one of the following: ",
                    paste(x = colnames(x = .tdr.obj$metadata),
                          collapse = ", ")))
      }
      
    }
    
    if((!is.null(x = .dodge.by)) &
       (!is.null(x = .line.by))){
      stop("only one of .dodge.by or .line.by can be provided")
    }
    
    .orientation <-
      match.arg(arg = .orientation,
                choices = c("wide",
                            "square"))
    
    dat.df <-
      data.frame(
        sample = rownames(x = .tdr.obj$metadata),
        x = .tdr.obj$metadata[[.x.split]]
      )
    
    if(!is.null(x = .dodge.by)){
      dat.df$dodge <-
        .tdr.obj$metadata[[.dodge.by]]
    }
    
    if(!is.null(x = .line.by)){
      dat.df$group <-
        .tdr.obj$metadata[[.line.by]]
    }
    
    if(is.null(x = .pop)){
      
      dat.df <-
        cbind(dat.df,
              .tdr.obj$map[[.pop.from]]$cell.perc) |>
        tidyr::pivot_longer(
          cols = colnames(x = .tdr.obj$map[[.pop.from]]$cell.perc),
          cols_vary = "slowest"
        ) |>
        as.data.frame()
      
    } else {
      
      dat.df$value <-
        .tdr.obj$map[[.pop.from]]$cell.perc[,.pop]
      
    }
    
    if(!is.null(x = .x.split.subset)){
      
      dat.df <- 
        droplevels(x = dat.df[dat.df$x %in% .x.split.subset,])
      
    }
    
    if(isTRUE(x = .order.pop)){
      
      dat.df$name <-
        as.character(x = dat.df$name) |>
        factor(levels = .tdr.obj$graph[[.pop.from]]$pheatmap$tree_row$labels[
          .tdr.obj$graph[[.pop.from]]$pheatmap$tree_row$order
        ])
      
    }
    
    p <-
      ggplot2::ggplot(data = dat.df,
                      mapping = ggplot2::aes(x = x,
                                             y = value)) +
      {
        if(is.null(x = .pop)){
          if(.orientation == "wide"){
            ggplot2::facet_grid(cols = ggplot2::vars(name))
          } else {
            ggplot2::facet_wrap(facets = ggplot2::vars(name), 
                                ncol = unique(x = dat.df$name) |>
                                  length() |>
                                  sqrt() |>
                                  ceiling(),
                                scales = "free_y")    
          }
        } else {
          ggplot2::ggtitle(label = .pop)
        }
      } +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                         vjust = 1,
                                                         angle = 30),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(x = .x.split,
                    y = "% of cells") +
      ggh4x::force_panelsizes(cols = ggplot2::unit(x = (unique(x = dat.df$x) |>
                                                          length()) * 
                                                     .x.space.scaler,
                                                   units = "in"),
                              rows = ggplot2::unit(x = .height,
                                                   units = "in"))
    
    if(!is.null(x = .dodge.by)){
      
      p <-
        p +
        ggplot2::geom_boxplot(mapping = ggplot2::aes(fill = dodge),
                              outliers = FALSE) +
        ggplot2::geom_point(mapping = ggplot2::aes(fill = dodge),
                            size = I(x = 1),
                            position = ggplot2::position_jitterdodge(jitter.width =  0.1,
                                                                     jitter.height = 0,
                                                                     seed = .seed)) +
        ggplot2::scale_fill_manual(
          values = grDevices::colorRampPalette(
            colors = unname(
              obj = .cat.feature.color)
          )(length(x = unique(x = dat.df$dodge)))) 
      
    } else if(!is.null(x = .line.by)){
      
      p <-
        p +
        ggplot2::geom_boxplot(outliers = FALSE,
                              color = "grey50") +
        ggplot2::geom_point(size = I(x = 1)) +
        ggplot2::geom_line(mapping = ggplot2::aes(group = group),
                           linetype = "dashed",
                           alpha = 1/4)
      
    } else {
      p <-
        p +
        ggplot2::geom_boxplot(outliers = FALSE,
                              color = "grey50") +
        ggplot2::geom_point(size = I(x = 1),
                            position = ggplot2::position_jitter(width = 0.2,
                                                                height = 0,
                                                                seed = .seed))
    }
    
    if(isTRUE(x = .log2.y)){
      p <-
        p +
        ggplot2::scale_y_continuous(transform = "log2")
    }
    
    p
    
  }

#' Plot Density
#'
#' Creates dot/line plots showing log2-transformed landmark fuzzy densities per sample. 
#' Unlike \code{plotTradPerc()} which plots cell percentages at the cluster/celltype level, this 
#' function plots landmark-level densities Useful for inspecting distributions and 
#' paired/longitudinal patterns at the landmark resolution.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.map()}.
#' @param .x.split Character specifying metadata column for x-axis grouping (default first column).
#' @param .pop Character vector of population names to plot. If NULL, plots all landmarks. If 
#'   specified, plots only landmarks belonging to that population (from \code{.pop.from}).
#' @param .pop.from Character: "clustering" (default) or "celltyping" - which grouping to use for 
#'   filtering landmarks when \code{.pop} is specified.
#' @param .subject.id Character metadata column for connecting paired samples with lines (e.g., 
#'   "Subject"). Default NULL.
#' @param .color.by Character metadata column for coloring points. Default NULL (all black).
#' @param .x.space.scaler Numeric x-axis width scaling (default 0.25 inches per group).
#' @param .height Numeric plot height in inches (default 1.5).
#' @param .cat.feature.color Character vector of colors (default \code{Color.Palette[1,1:5]}).
#' @param .seed Integer random seed for jitter (default 123).
#' @param .orientation Character: "wide" (default) or "square" faceting.
#'   
#' @return A \code{ggplot} object showing log2-transformed landmark densities.
#'   
#' @seealso \code{\link{plotTradPerc}}, \code{\link{plotTradStats}}
#' 
#' @examples
#' \dontrun{
#' # Basic density plot
#' plotDensity(lm.cells, .pop = "B.cells")
#' 
#' # With paired subject lines
#' plotDensity(lm.cells, 
#'               .pop = "B.cells",
#'               .subject.id = "MouseID",
#'               .color.by = "Treatment")
#' }
#' 
#' @export
#'
plotDensity <-
  function(
    .tdr.obj,
    .x.split = colnames(x = .tdr.obj$metadata)[1],
    .pop = NULL,
    .pop.from = "clustering",
    .subject.id = NULL,
    .color.by = NULL,
    .x.space.scaler = 0.25,
    .height = 1.5,
    .cat.feature.color = Color.Palette[1,1:5],
    .seed = 123,
    .orientation = "wide"
  ){
    
    value <- name <- color <- group <- x <- y <- NULL
    
    if(length(x = .x.split) != 1){
      stop(".x.split must be length 1")
    }
    
    if(!(.x.split %in% colnames(x = .tdr.obj$metadata))){
      stop(paste0(".x.split must be one of the following: ",
                  paste(x = colnames(x = .tdr.obj$metadata),
                        collapse = ", ")))
    }
    
    if(!is.null(x = .color.by)){
      if(length(x = .color.by) != 1){
        stop(".color.by must be length 1")
      }
      
      if(!(.color.by %in% colnames(x = .tdr.obj$metadata))){
        stop(paste0(".color.by must be one of the following: ",
                    paste(x = colnames(x = .tdr.obj$metadata),
                          collapse = ", ")))
      }
      
    }
    
    .orientation <-
      match.arg(arg = .orientation,
                choices = c("wide",
                            "square"))
    
    dat.df <-
      data.frame(
        sample = rownames(x = .tdr.obj$metadata),
        x = .tdr.obj$metadata[[.x.split]],
        color = if(is.null(x = .color.by)) "black" else .tdr.obj$metadata[[.color.by]]
      ) |> 
      (\(x)
       `rownames<-`(x = x,
                    value = rownames(x = .tdr.obj$metadata))
      )()
    
    if(is.null(x = .pop)){
      
      dat.df <-
        cbind(dat.df,
              Matrix::t(x = log2(x = .tdr.obj$map$fdens + 0.5))) |>
        tidyr::pivot_longer(
          cols = rownames(x = .tdr.obj$map$fdens),
          cols_vary = "slowest"
        ) |>
        dplyr::mutate(name = as.character(x = .tdr.obj$graph[[.pop.from]]$ids) |>
                        stats::setNames(nm = rownames(x = .tdr.obj$map$fdens)) |>
                        (\(x)
                         x[name]
                        )()) |>
        as.data.frame()
      
    } else {
      
      dat.df <-
        cbind(dat.df,
              Matrix::t(x = log2(x = .tdr.obj$map$fdens[.tdr.obj$graph[[.pop.from]]$ids == .pop,] + 0.5))) |>
        tidyr::pivot_longer(
          cols = rownames(.tdr.obj$map$fdens),
          cols_vary = "slowest"
        ) |>
        as.data.frame()
      
    }
    
    if(!is.null(x = .subject.id)){
      dat.df$group <-
        .tdr.obj$metadata[[.subject.id]][match(x = dat.df$sample,
                                               table = rownames(x = .tdr.obj$metadata))]  
    }
    
    p <-
      ggplot2::ggplot(data = dat.df,
                      mapping = ggplot2::aes(x = x,
                                             y = value)) +
      {
        if(is.null(x = .pop)){
          if(.orientation == "wide"){
            ggplot2::facet_grid(cols = ggplot2::vars(name))
          } else {
            ggplot2::facet_wrap(facets = ggplot2::vars(name), 
                                ncol = unique(x = dat.df$name) |>
                                  length() |>
                                  sqrt() |>
                                  ceiling(),
                                scales = "free_y")    
          }
        } else {
          ggplot2::ggtitle(label = .pop)
        }
      } +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1,
                                                         vjust = 1,
                                                         angle = 45),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(x = .x.split,
                    y = "density\nlog2(+0.5)",
                    color = .color.by) +
      ggplot2::geom_boxplot(outliers = FALSE,
                            color = "grey50") +
      {
        if(!is.null(x = .subject.id)){
          
          ggplot2::geom_line(mapping = ggplot2::aes(group = group))}
        
      } +
      ggh4x::force_panelsizes(cols = ggplot2::unit(x = (unique(x = dat.df$x) |>
                                                          length()) * 
                                                     .x.space.scaler,
                                                   units = "in"),
                              rows = ggplot2::unit(x = .height,
                                                   units = "in"))
    
    if(!is.null(x = .color.by)){
      p <-
        p +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
        scattermore::geom_scattermore(mapping = ggplot2::aes(color = color),
                                      pointsize = 0,
                                      position = ggplot2::position_jitter(width = 0.2,
                                                                          height = 0,
                                                                          seed = .seed)) +
        ggplot2::scale_color_manual(
          values = grDevices::colorRampPalette(
            colors = unname(
              obj = Color.Palette)
          )(length(x = unique(x = dat.df$color)))) 
      
    } else {
      p <-
        p +
        scattermore::geom_scattermore(pointsize = 0,
                                      position = ggplot2::position_jitter(width = 0.2,
                                                                          height = 0,
                                                                          seed = .seed))
    }
    
    p
    
  }

#' Plot Pseudobulk Differential Expression Results
#'
#' Visualizes pseudobulk differential expression results as a heatmap showing log fold changes 
#' for genes/markers across coefficients. Rows ordered by cluster/celltype expression patterns, 
#' with significance indicated by asterisks. Helps identify which features drive population-level changes.
#'
#' @param .tdr.obj A tinydenseR object processed through \code{get.pbDE()}.
#' @param .de.obj Optional: differential expression results object. If NULL (default), retrieves 
#'   results from \code{.tdr.obj$pbDE[[.model.name]][[.population.name]]}.
#' @param .model.name Character: model name to retrieve from \code{.tdr.obj$pbDE} (default "default").
#'   Ignored if \code{.de.obj} is provided.
#' @param .population.name Character: population name to retrieve (default "all"). 
#'   Ignored if \code{.de.obj} is provided.
#' @param .coefs Character vector of coefficient names to plot. Defaults to all coefficients.
#' @param .order.by Character: "clustering" (default) or "celltyping" - order rows by mean expression 
#'   in these groups.
#' @param .markers Character vector of feature names (genes/proteins) to plot. Defaults to features 
#'   shown in cluster/celltype heatmap (top PC contributors for RNA, all markers for cytometry).
#' @param .q Numeric adjusted p-value threshold for significance marking (default 0.1).
#' @param .row.space.scaler Numeric row height scaling (default 0.2 inches per feature).
#' @param .col.space.scaler Numeric column width scaling (default 0.065 inches per coefficient).
#' @param .label.substr.rm Character substring to remove from labels (default "").
#'   
#' @return A \code{ggplot} heatmap showing effect sizes (log2 fold changes for RNA, estimated 
#'   differences for cytometry) with point size indicating adjusted p-values and asterisks marking 
#'   features meeting the significance threshold.
#'   
#' @details
#' This function shows which genes/markers are differentially expressed between conditions, 
#' organized by their expression patterns across clusters/cell types. The ordering helps identify:
#' \itemize{
#'   \item Marker genes defining specific populations
#'   \item Broadly vs. specifically regulated features
#'   \item Cell type-specific transcriptional responses
#' }
#' 
#' For RNA data, typically shows top PC-loading genes. For cytometry, shows all markers.
#' 
#' @seealso \code{\link{get.pbDE}}, \code{\link{plotBeeswarm}}
#' 
#' @examples
#' \dontrun{
#' # After pbDE analysis
#' design <- model.matrix(~ Condition, data = .meta)
#' lm.cells <- get.pbDE(lm.cells, .design = design)
#' 
#' # Heatmap of DE genes (uses .tdr.obj$pbDE$default$all)
#' plotPbDE(lm.cells, .coefs = "ConditionB", .order.by = "clustering")
#' 
#' # Plot results from specific population
#' lm.cells <- get.pbDE(lm.cells, .design = design, .id = "1", 
#'                      .id.from = "clustering", .population.name = "cluster1")
#' plotPbDE(lm.cells, .population.name = "cluster1", .coefs = "ConditionB")
#' 
#' # Focus on specific markers
#' plotPbDE(lm.cells, 
#'          .coefs = "ConditionB",
#'          .markers = c("CD4", "CD8A", "CD3D"))
#' }
#' 
#' @export
#'
plotPbDE <-
  function(
    .tdr.obj,
    .de.obj = NULL,
    .model.name = "default",
    .population.name = "all",
    .coefs = NULL,
    .order.by = "none",
    .markers = .tdr.obj$config$markers,
    .q = 0.1,
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.065,
    .label.substr.rm = ""
  ){
    
    # R CMD check appeasement
    sig.shape <- sig <- adj.p <- level <- marker <- term <- coef <- cell.perc <- NULL
    
    # Get DE results from .tdr.obj$pbDE or from provided .de.obj
    if(is.null(x = .de.obj)){
      if(is.null(x = .tdr.obj$pbDE[[.model.name]][[.population.name]])){
        stop(sprintf("No results found at .tdr.obj$pbDE$%s$%s. Run get.pbDE() first.", 
                     .model.name, .population.name))
      }
      .de.obj <- .tdr.obj$pbDE[[.model.name]][[.population.name]]
    }
    
    # Default coefficients
    if(is.null(x = .coefs)){
      .coefs <- colnames(x = .de.obj$coefficients)
    }
    
    # check .order.by
    .order.by <-
      match.arg(arg = .order.by,
                choices = c("none",
                            "clustering",
                            "celltyping"))
    
    # Extract log fold changes for selected markers and coefficients
    coef.df <-
      cbind(.de.obj$coefficients[.markers,.coefs,drop = FALSE]) |>
      dplyr::as_tibble(rownames = "marker") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "coef")
    
    # Extract adjusted p-values
    adj.p.df <-
      cbind(.de.obj$adj.p[.markers,.coefs,drop = FALSE])  |>
      dplyr::as_tibble(rownames = "marker") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "adj.p")
    
    # Combine coefficients and p-values, determine significance
    dat.df <-
      dplyr::full_join(x = coef.df,
                       y = adj.p.df,
                       by = c("marker","term")) |>
      # Order rows by cluster/celltype mean expression patterns
      dplyr::mutate(marker =
                      factor(
                        x = marker,
                        levels = if(.order.by == "none") .tdr.obj$config$markers else {
                          rev(x = .tdr.obj$graph[[.order.by]]$pheatmap$tree_col$labels[
                          .tdr.obj$graph[[.order.by]]$pheatmap$tree_col$order]) |>
                          (\(x)
                           x[x %in% marker]
                          )()
                          }),
                    # Classify as higher/lower/not significant
                    sig =
                      ifelse(
                        test = coef < 0,
                        yes = "lower",
                        no = "higher") |>
                      ifelse(
                        test = adj.p < .q,
                        no = "not sig.")  |>
                      factor(levels = c("lower",
                                        "not sig.",
                                        "higher"))) |>
      # Create indicator for asterisk overlay
      dplyr::mutate(sig.shape = ifelse(
        test = sig == "not sig.",
        yes = NA,
        no = 1))
    
    # Clean coefficient labels
    dat.df$term <-
      gsub(pattern = .label.substr.rm,
           replacement = "",
           x = dat.df$term,
           fixed = FALSE) |>
      (\(x)
       factor(x = x,
              levels = unique(x = x))
      )()
    
    # Create heatmap with fold changes colored and point size by -log10(p)
    other.plot <-
      ggplot2::ggplot(data = dat.df,
                      mapping = ggplot2::aes(x = term,
                                             y = marker,
                                             fill = coef,
                                             stroke = sig.shape
                      )) +
      # Order legend: significance, adj.p, then fold change
      ggplot2::guides(stroke = ggplot2::guide_legend(override.aes = list(size = I(x = 5)),
                                                     order = 1),
                      size = ggplot2::guide_legend(override.aes = list(fill = unname(obj = Color.Palette[1,6]),
                                                                       stroke = NA),
                                                   order = 2),
                      fill = ggplot2::guide_colorbar(order = 3)) +
      # Plot circles sized by -log10(adj.p)
      ggplot2::geom_point(shape = 21,
                          mapping = ggplot2::aes(size = -log10(x = adj.p)),
                          color = "black",
                          show.legend = TRUE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),
                     axis.text.x = ggplot2::element_text(angle = 0,
                                                         hjust = 0.5),
                     legend.position = "right") +
      # Diverging color scale: blue (lower) to white (no change) to red (higher)
      ggplot2::scale_fill_gradient2(low = unname(obj = Color.Palette[1,1]),
                                    mid = unname(obj = Color.Palette[1,6]),
                                    high = unname(obj = Color.Palette[1,2]),
                                    midpoint = 0#,
                                    #labels = ~ ifelse(test = .x >= 0,
                                    #                  yes = round(x = 2^.x,
                                    #                              digits = 1),
                                    #                  no = -round(x = 1 / 2^.x,
                                    #                              digits = 1))
      ) +
      # Size legend shows original p-values
      ggplot2::scale_size_continuous(labels = ~ formatC(x = 10^(-(.x)),
                                                        format = "g",
                                                        digits = 1),
                                     range = I(x = c(3,6))) +
      # Add significance stroke aesthetic if any features meet q-value threshold
      {if(!(is.na(x = dat.df$sig.shape) |> all())){
        ggplot2::continuous_scale(aesthetics = "stroke",
                                  name =  paste0("q < ",
                                                 .q),
                                  palette = function(x){scales::rescale(x = x, c(0, 1))},
                                  labels = "TRUE") 
      }} +
      # Overlay asterisks for significant features
      {if(!(is.na(x = dat.df$sig.shape) |> all())){
        ggplot2::geom_point(shape = 8,
                            mapping = ggplot2::aes(size = -log10(x = adj.p)),
                            color = "black",
                            show.legend = TRUE)
      }} +
      ggplot2::labs(title = "difference in expression",
                    x = "",
                    y = "",
                    fill = if(.tdr.obj$config$assay.type == "RNA"){"log2(+0.5)FC"} else {"estimated\ndifference"},
                    size = if((is.na(x = dat.df$sig.shape) |> all())) {
                      paste0("adj.p\n(",
                             "all > ",
                             .q,
                             ")")} else "adj.p") +
      ggh4x::force_panelsizes(
        col = grid::unit(x = pmax(
          (as.character(x = dat.df$term) |> unique() |> nchar() |> max(na.rm = TRUE)) *
            (unique(x = dat.df$term) |> length()) *
            .col.space.scaler,
          0.5, 
          na.rm = TRUE),
          units = "in"),
        rows = grid::unit(x = pmax((unique(x = dat.df$marker) |> length()) * .row.space.scaler,
                                   3.5, 
                                   na.rm = TRUE),
                          units = "in"))
    
    return(other.plot)
    
  }

#' Plot Differential Expression Analysis Results (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' \code{plotDEA()} has been renamed to \code{\link{plotPbDE}()} for clarity. This function
#' is provided for backward compatibility and will be removed in a future version.
#'
#' @param .tdr.obj A tinydenseR object.
#' @param .dea.obj Differential expression results from \code{get.dea()} or \code{get.pbDE()}.
#' @param .coefs Character vector of coefficient names to plot.
#' @param .order.by Character: "clustering" or "celltyping".
#' @param .markers Character vector of feature names.
#' @param .q Numeric adjusted p-value threshold.
#' @param .row.space.scaler Numeric row height scaling.
#' @param .col.space.scaler Numeric column width scaling.
#' @param .label.substr.rm Character substring to remove from labels.
#'
#' @return A \code{ggplot} heatmap.
#' @seealso \code{\link{plotPbDE}}
#' @keywords internal
#' @export
plotDEA <- function(
    .tdr.obj,
    .dea.obj,
    .coefs = colnames(x = .dea.obj$coefficients),
    .order.by = "clustering",
    .markers = colnames(x = .tdr.obj$graph[[.order.by]]$median.exprs),
    .q = 0.1,
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.065,
    .label.substr.rm = ""
) {
  .Deprecated("plotPbDE",
              msg = "plotDEA() is deprecated. Use plotPbDE() instead.")
  
  plotPbDE(
    .tdr.obj = .tdr.obj,
    .de.obj = .dea.obj,
    .coefs = .coefs,
    .order.by = .order.by,
    .markers = .markers,
    .q = .q,
    .row.space.scaler = .row.space.scaler,
    .col.space.scaler = .col.space.scaler,
    .label.substr.rm = .label.substr.rm
  )
}

#' Plot Marker DE Results
#'
#' Visualizes marker identification results from \code{get.markerDE} as a colored heatmap showing 
#' log fold changes with significance overlays. Uses the new storage pattern where results are 
#' stored in \code{.tdr.obj$markerDE[[.model.name]][[.comparison.name]]}.
#'
#' @param .tdr.obj A tinydenseR object with marker results from \code{get.markerDE}.
#' @param .de.obj Optional: Direct DE results object for backward compatibility. If NULL, fetches 
#'   from \code{.tdr.obj$markerDE[[.model.name]][[.comparison.name]]}.
#' @param .model.name Character identifying the model (default "default").
#' @param .comparison.name Character identifying the comparison (e.g., "cluster.1_vs_all"). Required 
#'   if \code{.de.obj} is NULL.
#' @param .coefs Character vector of coefficient names to plot. Default includes "(Intercept)" and 
#'   ".id1" to show overall expression and group-specific enrichment.
#' @param .order.by Source for marker ordering: "clustering" or "celltyping".
#' @param .markers Features to include (default: all markers from \code{.order.by}).
#' @param .q Numeric FDR threshold for significance markers (default 0.1).
#' @param .row.space.scaler Numeric row height scaling (default 0.2).
#' @param .col.space.scaler Numeric column width scaling (default 0.065).
#' @param .label.substr.rm Character pattern to remove from coefficient labels.
#'
#' @return A \code{ggplot} heatmap showing log fold changes per feature and coefficient.
#'
#' @details The plot shows:
#' \itemize{
#'   \item Fill color: Log fold change (blue = lower in .id1, red = higher in .id1)
#'   \item Point size: -log10(adjusted p-value)
#'   \item Asterisk overlay: Features passing the q-value threshold
#' }
#'
#' For marker analysis, the key coefficient is typically ".id1" which represents the difference 
#' between group 1 and group 2 (or all other landmarks).
#'
#' @seealso \code{\link{get.markerDE}}, \code{\link{plotPbDE}}
#'
#' @examples
#' \dontrun{
#' # After running get.markerDE
#' lm.cells <- get.markerDE(lm.cells, .id1 = "cluster.3", 
#'                          .comparison.name = "cluster3_markers")
#'
#' # Visualize marker results
#' plotMarkerDE(lm.cells, .comparison.name = "cluster3_markers",
#'              .coefs = c("(Intercept)", ".id1"))
#' }
#'
#' @export
plotMarkerDE <- function(
    .tdr.obj,
    .de.obj = NULL,
    .model.name = "default",
    .comparison.name = NULL,
    .coefs = NULL,
    .order.by = "none",
    .markers = .tdr.obj$config$markers,
    .q = 0.1,
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.065,
    .label.substr.rm = ""
) {
  
  # Get DE results from .tdr.obj$markerDE or from provided .de.obj
  if(is.null(x = .de.obj)){
    if(is.null(x = .comparison.name)){
      stop("Please provide .comparison.name to fetch results from .tdr.obj$markerDE")
    }
    if(is.null(x = .tdr.obj$markerDE[[.model.name]][[.comparison.name]])){
      stop(sprintf("No results found at .tdr.obj$markerDE$%s$%s. Run get.markerDE() first.", 
                   .model.name, .comparison.name))
    }
    .de.obj <- .tdr.obj$markerDE[[.model.name]][[.comparison.name]]
  }
  
  # Use plotPbDE with the extracted DE object
  plotPbDE(
    .tdr.obj = .tdr.obj,
    .de.obj = .de.obj,
    .coefs = .coefs,
    .order.by = .order.by,
    .markers = .markers,
    .q = .q,
    .row.space.scaler = .row.space.scaler,
    .col.space.scaler = .col.space.scaler,
    .label.substr.rm = .label.substr.rm
  )
}

#' Scatter Plot with Feature Coloring
#'
#' Creates a 2D scatter plot with flexible x/y features and optional coloring by a third feature. 
#' Useful for exploring relationships between any features in the tinydenseR object (scaled expression, 
#' PCA coordinates, graph embeddings, cluster IDs, metadata, etc.).
#'
#' @param .x.feature Numeric vector for x-axis values (e.g., \code{.tdr.obj$scaled.landmarks[,"CD3"]} or 
#'   \code{.tdr.obj$pca$embed[,"PC1"]}).
#' @param .y.feature Numeric vector for y-axis values.
#' @param .color.feature Optional vector for point colors. Can be numeric (continuous coloring with 
#'   diverging blue-white-red scale) or categorical (discrete colors). If \code{NULL}, all points 
#'   receive default ggplot2 coloring.
#' @param .x.label Character x-axis label (default "").
#' @param .y.label Character y-axis label (default "").
#' @param .color.label Character color legend label (default "").
#' @param .cat.feature.color Character vector of colors for categorical features 
#'   (default \code{Color.Palette[1,1:5]}). Automatically interpolated if more categories than 
#'   colors exist.
#' @param .panel.size Numeric vector \code{c(width, height)} in inches (default \code{c(2,2)}).
#' @param .midpoint Numeric midpoint for continuous color gradients. Defaults to median of 
#'   \code{.color.feature}. Useful for centering diverging scales.
#' @param .plot.title Character plot title (default "").
#' @param .legend.position Character: "right", "left", "top", "bottom", or "none" (default "right").
#' @param .point.size Numeric point size (default 0.1). Increase for smaller datasets.
#' @param .seed Integer random seed for point order randomization (default 123). Prevents systematic 
#'   overplotting of one group by another.
#' 
#' @return A \code{ggplot} object.
#' 
#' @details
#' This flexible plotting function enables custom visualizations beyond the standard \code{plotPCA} 
#' and \code{plotUMAP} interfaces. Use cases include:
#' \itemize{
#'   \item Plotting Laplacian Eigenmap coordinates from \code{.tdr.obj$graph$LE$embed}
#'   \item Exploring relationships between markers (e.g., CD3 vs CD4)
#'   \item Overlaying metadata on any 2D embedding
#'   \item Creating custom QC plots (e.g., library size vs mitochondrial %)
#' }
#' 
#' Points are plotted in randomized order (controlled by \code{.seed}) to avoid systematic 
#' overplotting bias when groups overlap. 
#' 
#' For numeric \code{.color.feature}, applies a diverging color scale (blue-white-red) centered 
#' at \code{.midpoint} (defaults to median). For categorical features, generates discrete colors 
#' by interpolating \code{.cat.feature.color}.
#' 
#' @seealso \code{\link{plotPCA}}, \code{\link{plotUMAP}}
#' 
#' @examples
#' \dontrun{
#' # After processing
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph()
#' 
#' # Plot PC1 vs PC2 colored by cluster
#' scatterPlot(.x.feature = lm.cells$pca$embed[,"PC1"],
#'             .y.feature = lm.cells$pca$embed[,"PC2"],
#'             .color.feature = lm.cells$graph$clustering$ids,
#'             .x.label = "PC1", .y.label = "PC2",
#'             .color.label = "Cluster")
#' 
#' # Marker expression relationship
#' scatterPlot(.x.feature = lm.cells$scaled.landmarks[,"CD4"],
#'             .y.feature = lm.cells$scaled.landmarks[,"CD8A"],
#'             .color.feature = .meta$Condition,
#'             .x.label = "CD4", .y.label = "CD8A",
#'             .color.label = "Condition")
#' 
#' # Laplacian Eigenmap coordinates
#' scatterPlot(.x.feature = lm.cells$graph$LE$embed[,1],
#'             .y.feature = lm.cells$graph$LE$embed[,2],
#'             .color.feature = lm.cells$graph$clustering$ids,
#'             .x.label = "LE1", .y.label = "LE2")
#' }
#' 
#' @export
#'
scatterPlot <-
  function(
    .x.feature,
    .y.feature,
    .color.feature = NULL,
    .x.label = "",
    .y.label = "",
    .color.label = "",
    .cat.feature.color = Color.Palette[1,1:5],
    .panel.size = c(2,2),
    .midpoint = NULL,
    .plot.title = "",
    .legend.position = "right",
    .point.size = 0.1,
    .seed = 123) {
    
    # Set midpoint to median for continuous features if not specified
    if(is.null(x = .midpoint) &
       is.numeric(x = .color.feature)){
      .midpoint <- stats::median(x = .color.feature)
    }
    
    # Build data frame from input vectors
    dat.df <-
      data.frame(.x.feature = .x.feature,
                 .y.feature = .y.feature)
    
    # Add color feature if provided
    if(!is.null(x = .color.feature)){
      dat.df <-
        cbind(dat.df,
              .color.feature = .color.feature)
    } else {
      dat.df$.color.feature <- 1
    }
    
    old.seed <- .Random.seed
    on.exit(expr = assign(x = ".Random.seed",
                          value = old.seed,
                          envir = .GlobalEnv),
            add = TRUE)
    set.seed(seed = .seed)
    p <-
      dat.df  |>
      (\(x)
       # randomize order to avoid plotting one class in front of the other
       # see: https://stackoverflow.com/a/29325361
       x[nrow(x = x) |>
           seq_len() |>
           sample(),]
      )() |>
      (\(x)
       ggplot2::ggplot(data = x,
                       mapping = ggplot2::aes(x = .x.feature,
                                              y = .y.feature,
                                              color = .color.feature)) +
         ggplot2::theme_bw() +
         ggplot2::theme(legend.position = .legend.position,
                        plot.title = ggplot2::element_text(hjust = 0.5),
                        plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
         ggplot2::labs(title = .plot.title,
                       x = .x.label,
                       y = .y.label,
                       color = .color.label) +
         ggplot2::geom_point(size = I(x = .point.size))
      )()
    
    # Apply color scale based on feature type
    if(is.numeric(x = .color.feature)){
      # Diverging gradient for continuous features
      p <-
        p +
        ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
                                       mid = unname(obj = Color.Palette[1,6]),
                                       high = unname(obj = Color.Palette[1,2]),
                                       midpoint = .midpoint)
    } else {
      # Discrete colors for categorical features
      p <-
        p  +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
        ggplot2::scale_color_manual(
          values = grDevices::colorRampPalette(
            colors = unname(
              obj = .cat.feature.color)
          )(length(x = unique(x = .color.feature))))
    }
    
    # Force panel dimensions
    p +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in")[2],
                              col = grid::unit(x = .panel.size,
                                               units = "in")[1])
  }

#' Plot Mean Expression Heatmap
#' 
#' Displays a heatmap of mean marker/gene expression across clusters or cell types. Shows the expression 
#' patterns that define each population, helping to validate cluster annotations and identify marker genes. 
#' The heatmap is automatically generated during \code{get.graph()} and stored for quick display.
#' 
#' @param .tdr.obj A tinydenseR object processed through \code{get.graph()}.
#' @param .id.from Character: "clustering" or "celltyping" (default "clustering"). Determines which 
#'   population definitions to show. Use "clustering" after \code{get.graph()}, or "celltyping" after 
#'   manual annotation with \code{celltyping()}.
#'   
#' @return A \code{gtable}/\code{grob} object (from \code{pheatmap}) rendered to the graphics device. 
#'   The heatmap includes hierarchical clustering of both features and populations.
#' 
#' @details
#' The heatmap content depends on assay type:
#' \itemize{
#'   \item \strong{Cytometry}: All markers from \code{.tdr.obj$marker}
#'   \item \strong{RNA-seq}: Top 3 positive and 3 negative genes per PC (from highly variable genes)
#' }
#' 
#' Rows (features) are ordered by hierarchical clustering to group co-expressed markers/genes. 
#' Columns (populations) also clustered to reveal relationships between cell types.
#' 
#' The heatmap is generated once during \code{get.graph()} and cached in 
#' \code{.tdr.obj$graph[[.id.from]]$pheatmap}. This function simply displays the cached version.
#' 
#' @seealso \code{\link{get.graph}}, \code{\link{celltyping}}
#' 
#' @examples
#' \dontrun{
#' # After clustering
#' lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |> get.graph()
#' 
#' # Show cluster marker heatmap
#' plotHeatmap(lm.cells, .id.from = "clustering")
#' 
#' # After manual annotation
#' lm.cells <- celltyping(lm.cells, 
#'                        .celltyping.map = list("CD4_T" = c("cluster.1", "cluster.3"),
#'                                               "CD8_T" = c("cluster.2")))
#' plotHeatmap(lm.cells, .id.from = "celltyping")
#' }
#' 
#' @export
#'
plotHeatmap <-
  function(
    .tdr.obj,
    .id.from = "clustering"
  ){
    
    # Verify graph component exists
    if(is.null(x = .tdr.obj$graph[[.id.from]])){
      stop(paste0("Please run get.graph first."))
    }
    
    # Validate id.from parameter
    .id.from <- 
      match.arg(
        arg = .id.from,
        choices = c("clustering",
                    "celltyping")
      )
    
    # Display cached pheatmap generated during get.graph()
    return(gridExtra::grid.arrange(.tdr.obj$graph[[.id.from]]$pheatmap$gtable))
    
  }

#' Plot specDE Components
#'
#' Visualizes spectral differential expression (specDE) results. When a component
#' is specified, shows UMAP colored by scores alongside Y-vs-specDE scatter.
#' When component is NULL, displays a diagnostic overview of all components.
#'
#' @param .tdr.obj A tinydenseR object after \code{get.specDE()}.
#' @param .coef.col Character: coefficient name matching a slot in \code{.tdr.obj$specDE}.
#' @param .specDE.dim Integer or NULL: component to visualize (1-indexed). If NULL,
#'   plots Ak vs Sk diagnostic scatter for all components.
#' @param .embed Character: either \code{"umap"} or \code{"pca"}. Default "umap".
#' @param .point.size Numeric: point size for scatter plots. Default 1.
#' @param .label.size Numeric: label size for diagnostic scatter plots. Default 3. Applies only if .specDE.dim is NULL.
#' @param .panel.size Numeric: panel size in inches. Default 2.
#' @param .seed Integer: random seed for point ordering. Default 123.
#'
#' @return A ggplot2 object (if \code{.specDE.dim = NULL}) or a patchwork composition
#'   (if \code{.specDE.dim} is integer).
#'
#' @details
#' \strong{Component view} (\code{.specDE.dim = integer}):
#' \itemize{
#'   \item Left panel: UMAP embedding colored by specDE scores
#'   \item Right panel: Scatter of Y vs specDE scores (essential for interpretation)
#' }
#'
#' \strong{Diagnostic view} (\code{.specDE.dim = NULL}):
#' \itemize{
#'   \item X-axis: Smoothness (Sk); higher = large-scale graph-smooth structure
#'   \item Y-axis: Y-alignment (Ak); higher = stronger density coupling
#'   \item Color: Variance explained (Vk)
#'   \item Labels: Component names (specDE1, specDE2, ...)
#' }
#'
#' \strong{Interpreting the diagnostic plot:}
#' \itemize{
#'   \item DA-dominated (typical): upper-right, high Ak + high Sk
#'   \item DE-dominated (functional): lower region, moderate/low Ak
#'   \item Structural/constraint: low Ak with high Vk (bright color), regardless of Sk
#' }
#' No single metric determines DA vs DE---use Y vs score scatter to confirm.
#'
#' @seealso \code{\link{get.specDE}} for computing specDE
#'
#' @examples
#' \dontrun{
#' # After running specDE
#' lm.obj <- get.specDE(lm.obj, .coef.col = "Infection")
#'
#' # Diagnostic overview
#' plotSpecDE(lm.obj, .coef.col = "Infection", .specDE.dim = NULL)
#'
#' # Visualize first component
#' plotSpecDE(lm.obj, .coef.col = "Infection", .specDE.dim = 1)
#' }
#'
#' @export
#'
plotSpecDE <-
  function(
    .tdr.obj,
    .coef.col,
    .specDE.dim = NULL,
    .embed = "umap",
    .point.size = 1,
    .label.size = 3,
    .panel.size = 2,
    .seed = 123
  ) {
    
    # R CMD check appeasement
    Sk <- Ak <- Vk <- component <- PC1 <- PC2 <- umap.1 <- umap.2 <- Y <- score <- NULL
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    .embed <-
      match.arg(arg = .embed,
                choices = c(
                  "umap",
                  "pca"
                ))
    
    if (is.null(x = .tdr.obj$specDE[[.coef.col]])) {
      avail <-
        if (is.null(x = .tdr.obj$specDE)) {
          "none (run get.specDE() first)"
        } else {
          paste(names(x = .tdr.obj$specDE), collapse = ", ")
        }
      stop("specDE results for '", .coef.col, "' not found.\n",
           "Available: ", avail)
    }
    
    specDE.res <-
      .tdr.obj$specDE[[.coef.col]]
    
    # -------------------------------------------------------------------------
    # Diagnostic view: Sk x Ak scatter
    # -------------------------------------------------------------------------
    
    if (is.null(x = .specDE.dim)) {
      
      diag.df <-
        data.frame(
          component = gsub(pattern = "specDE",
                           replacement = "",
                           x = names(x = specDE.res$smoothness)),
          Sk = specDE.res$smoothness,
          Ak = specDE.res$Y.alignment,
          Vk = specDE.res$var.explained
        )
      
      p <-
        ggplot2::ggplot(data = diag.df,
                        mapping = ggplot2::aes(x = Sk,
                                               y = Ak,
                                               color = Vk*100,
                                               label = component)) +
        ggplot2::geom_point(size = I(x = .point.size)) +
        ggrepel::geom_text_repel(size = I(x = .label.size),
                                 max.overlaps = 20, 
                                 seed = .seed) +
        ggplot2::scale_color_gradient(low = unname(obj = Color.Palette[1,1]),
                                      high = unname(obj = Color.Palette[1,2]),
                                      name = "% of var") +
        ggplot2::labs(
          title = paste0("specDE diagnostics: ", .coef.col),
          x = "Smoothness (Sk)",
          y = "Y-alignment (Ak)"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                  units = "in"),
                                cols = grid::unit(x = .panel.size,
                                                  units = "in"))
      
      return(p)
      
    }
    
    # -------------------------------------------------------------------------
    # Component view: UMAP + Y-vs-score scatter
    # -------------------------------------------------------------------------
    
    if (!is.numeric(x = .specDE.dim) ||
        length(x = .specDE.dim) != 1 ||
        .specDE.dim < 1 ||
        .specDE.dim > ncol(x = specDE.res$scores)) {
      stop(".specDE.dim must be an integer between 1 and ", ncol(x = specDE.res$scores))
    }
    
    .specDE.dim <-
      as.integer(x = .specDE.dim)
    
    comp.name <-
      colnames(x = specDE.res$scores)[.specDE.dim]
    
    score.vec <-
      specDE.res$scores[, .specDE.dim]
    
    # Check UMAP exists
    if (is.null(x = .tdr.obj$graph$uwot$embedding)) {
      stop("UMAP embedding not found. Run get.graph() first.")
    }
    
    # Build data frames
    old.seed <- .Random.seed
    on.exit(expr = assign(x = ".Random.seed",
                          value = old.seed,
                          envir = .GlobalEnv),
            add = TRUE)
    set.seed(seed = .seed)
    
    embed.df <-
      (if(.embed == "umap"){
        as.data.frame(x = .tdr.obj$graph$uwot$embedding)
      } else {
        as.data.frame(x = .tdr.obj$pca$embed)
      }) |>
      cbind(score = score.vec) |>
      (\(x)
       x[sample(x = nrow(x = x)), ]
      )()
    
    # Get Y (centered)
    Y.vec <- specDE.res$Y

    scatter.df <-
      data.frame(
        Y = Y.vec,
        score = score.vec,
        col = .tdr.obj$map$lm[[specDE.res$params$model.name]]$fit$coefficients[, specDE.res$params$coef.col]
      )
    
    # Panel 1: UMAP colored by scores
    p.embed <-
      ggplot2::ggplot(data = embed.df,
                      mapping = ggplot2::aes(x = if(.embed == "umap") umap.1 else PC1,
                                             y = if(.embed == "umap") umap.2 else PC2,
                                             color = score)) +
      ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                      title.hjust = 0.5)) +
      ggplot2::geom_point(size = I(x = .point.size)) +
      ggplot2::scale_color_gradient2(
        low = unname(obj = Color.Palette[1, 1]),
        mid = unname(obj = Color.Palette[1, 6]),
        high = unname(obj = Color.Palette[1, 2]),
        midpoint = 0,
        name = comp.name
      ) +
      ggplot2::labs(
        title = comp.name,
        x = if(.embed == "umap") "umap.1" else "PC1",
        y = if(.embed == "umap") "umap.2" else "PC2"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.margin = ggplot2::margin(t = -0.1, 
                                        unit = "in")
      ) +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              cols = grid::unit(x = .panel.size,
                                                units = "in"))
    
    # Panel 2: Y vs score scatter
    Y.label <- 
      "Density contrast"
    
    p.scatter <-
      ggplot2::ggplot(data = scatter.df,
                      mapping = ggplot2::aes(x = Y,
                                             y = score,
                                             color = col)) +
      ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top",
                                                      title.hjust = 0.5)) +
      ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
                                     mid = unname(obj = Color.Palette[1,6]),
                                     high = unname(obj = Color.Palette[1,2]),
                                     midpoint = 0)+
      ggplot2::geom_vline(xintercept = 0,
                          color = "black",
                          linetype = "dashed",
                          linewidth = I(x = 1/2)) +
      ggplot2::geom_hline(yintercept = 0,
                          color = "black",
                          linetype = "dashed",
                          linewidth = I(x = 1/2)) +
      ggplot2::geom_point(size = I(x = .point.size)) +
      ggplot2::labs(
        title = paste0(.coef.col, " vs ", comp.name),
        subtitle = sprintf("Ak = %.2f, Sk = %.2f, Vk = %.0f%%",
                           specDE.res$Y.alignment[.specDE.dim],
                           specDE.res$smoothness[.specDE.dim],
                           specDE.res$var.explained[.specDE.dim] * 100),
        x = paste0("Density contrast", "\n(centered ", .coef.col, ")"),
        y = comp.name,
        color = paste0(.coef.col,
                       " (raw)")
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.margin = ggplot2::margin(t = -0.1, 
                                        unit = "in")
      ) +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              cols = grid::unit(x = .panel.size,
                                                units = "in"))
    
    # Combine with patchwork
    p.embed + p.scatter +
      patchwork::plot_layout(ncol = 2)
    
  }

#' Plot specDE Expression Heatmap
#'
#' Creates an efficient expression heatmap with landmarks as columns and top-loaded
#' features as rows. Landmarks are ordered by specDE scores (or custom ordering),
#' and features are ranked by absolute loading. For RNA data, expression is properly
#' normalized, log-transformed, and centered.
#'
#' @param .tdr.obj A tinydenseR object after \code{get.specDE()}.
#' @param .coef.col Character: coefficient name matching a slot in \code{.tdr.obj$specDE}.
#' @param .specDE.dim Integer or integer vector: specDE component(s) to use for:
#'   (1) ranking features by absolute loading, and (2) ordering landmarks (if
#'   \code{.order.by = NULL}). If vector (e.g., \code{c(1, 2)}), the first dimension
#'   is used for feature ranking, and landmarks are sorted by all dimensions sequentially.
#'   Default 1.
#' @param .model.name Character: name of the fitted model (default "default").
#' @param .n.features Integer: maximum number of features to display. For RNA, features
#'   are ranked by absolute loading and top \code{.n.features} are shown. For cytometry,
#'   all markers are shown and this parameter is ignored. Default 50.
#' @param .order.by Numeric matrix with column names: custom ordering for landmarks.
#'   If NULL (default), uses centered density contrast from \code{.coef.col}. If provided,
#'   landmarks are sorted by columns sequentially (first column primary, etc.) and
#'   each column is displayed as an annotation strip with the column name as legend title.
#' @param .add.annot Numeric matrix with column names: displayed as annotation strips.
#'   These are added as additional annotation strips after the density contrast and before 
#'   specDE score strips.
#' @param .order.decreasing Logical: sort landmarks in decreasing order? Default FALSE.
#' @param .viridis.options.annot Character vector: viridis color options for annotation
#'   strips (density contrast + specDE dimensions + .order.by columns). Cycled if fewer
#'   options than strips. Default \code{c("viridis", "plasma", "cividis")}.
#' @param .annot.panel.width Numeric: width of annotation strips in inches. Default 4.
#' @param .annot.panel.height Numeric: height of each annotation strip in inches. Default 0.15.
#' @param .panel.width Numeric: width of expression heatmap panel in inches. Default 4.
#' @param .panel.height Numeric: height of expression heatmap panel in inches. Default 3.
#' @param .feature.font.size Numeric: font size for feature labels. Default 7.
#' @param .show.landmark.labels Logical: show landmark IDs on x-axis? Default FALSE
#'   (usually too many to display).
#'
#' @return A ggplot2 object (or patchwork composition with annotation strips).
#'
#' @details
#' \strong{Feature selection and ordering}:
#' Features are \emph{selected} by absolute loading value for the first element of \code{.specDE.dim}
#' (for RNA, the top \code{.n.features}; for cytometry, all markers). However, the selected features 
#' are then \emph{ordered} by signed loading (highest positive at top, most negative at bottom), 
#' making it easy to see which genes are up vs down-regulated along the specDE component.
#'
#' \strong{Expression normalization} (RNA only):
#' To avoid processing the full expression matrix, the function:
#' \enumerate{
#'   \item Computes row sums (per-landmark library size) from the entire \code{raw.landmarks}
#'   \item Subsets to top features
#'   \item Applies size factor normalization using pre-computed row sums
#'   \item Log2-transforms: \code{log2(x + 1)}
#'   \item Centers each feature (row) to mean 0
#' }
#'
#' \strong{Landmark ordering}:
#' When \code{.order.by} is multi-dimensional (matrix with >1 column), landmarks are
#' sorted by the first column, then by the second within ties, etc. (equivalent to
#' \code{dplyr::arrange(df, col1, col2, ...)}). If \code{.order.by = NULL}, landmarks
#' are ordered by specDE scores from \code{.specDE.dim}.
#'
#' \strong{Annotation strips}:
#' The following are shown as annotation strips at the top:
#' \enumerate{
#'   \item Density contrast (centered Y from \code{get.lm})
#'   \item specDE scores for all dimensions in \code{.specDE.dim}
#'   \item If \code{.order.by} is provided: each column of \code{.order.by}, with
#'     column names used as legend titles
#' }
#' This ensures the relationship between expression patterns and the specDE structure
#' is always visible, plus any custom ordering variables.
#'
#' \strong{Visualization}:
#' Uses \code{ggplot2::geom_raster()} for efficient rendering of large heatmaps.
#' Expression uses a diverging color scale centered at 0. Annotation strips use
#' sequential viridis scales.
#'
#' @seealso \code{\link{get.specDE}} for computing specDE, \code{\link{plotSpecDE}}
#'   for score visualization
#'
#' @examples
#' \dontrun{
#' # After running specDE
#' lm.obj <- get.specDE(lm.obj, .coef.col = "Infection")
#'
#' # Basic heatmap using specDE1
#' plotSpecDEHeatmap(lm.obj, .coef.col = "Infection", .specDE.dim = 1)
#'
#' # Order by two dimensions (e.g., to see 2D structure)
#' plotSpecDEHeatmap(lm.obj, .coef.col = "Infection", .specDE.dim = c(1, 2))
#'
#' # Custom ordering by density contrast Y
#' plotSpecDEHeatmap(lm.obj, .coef.col = "Infection", .specDE.dim = 1,
#'                   .order.by = lm.obj$specDE$Infection$Y)
#'
#' # More features for cytometry
#' plotSpecDEHeatmap(lm.obj, .coef.col = "Treatment", .specDE.dim = 1,
#'                   .n.features = 100)
#' }
#'
#' @export
#'
plotSpecDEHeatmap <-
  function(
    .tdr.obj,
    .coef.col,
    .specDE.dim = 1,
    .model.name = "default",
    .n.features = 50,
    .order.by = NULL,
    .add.annot = NULL,
    .order.decreasing = FALSE,
    .viridis.options.annot = c("cividis", "rocket", "inferno", "mako", "magma"),
    .annot.panel.width = 4,
    .annot.panel.height = 0.15,
    .panel.width = 4,
    .panel.height = 3,
    .feature.font.size = 7,
    .show.landmark.labels = FALSE
  ) {
    
    # R CMD check appeasement
    loading <- x <- landmark <- feature <- expr <- annot_dim <- annot_value <- NULL
    
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    
    if (is.null(x = .tdr.obj$specDE[[.coef.col]])) {
      avail <-
        if (is.null(x = .tdr.obj$specDE)) {
          "none (run get.specDE() first)"
        } else {
          paste(names(x = .tdr.obj$specDE), collapse = ", ")
        }
      stop("specDE results for '", .coef.col, "' not found.\n",
           "Available: ", avail)
    }
    
    specDE.res <-
      .tdr.obj$specDE[[.coef.col]]
    
    # Validate .specDE.dim
    if (!is.numeric(x = .specDE.dim) ||
        any(.specDE.dim < 1) ||
        any(.specDE.dim > ncol(x = specDE.res$scores))) {
      stop(".specDE.dim must be integer(s) between 1 and ",
           ncol(x = specDE.res$scores))
    }
    
    .specDE.dim <-
      as.integer(x = .specDE.dim)
    
    if (is.null(x = .tdr.obj$raw.landmarks)) {
      stop("Raw landmarks not found. Run get.landmarks() first.")
    }
    
    # -------------------------------------------------------------------------
    # Feature selection: rank by signed loading
    # -------------------------------------------------------------------------
    
    # Use the first .specDE.dim for feature ranking
    loading.dim <-
      .specDE.dim[1]
    
    loadings <-
      specDE.res$loadings[, loading.dim]
    
    # Sort by signed loading (most positive first)
    feat.signed.order <-
      order(loadings, decreasing = TRUE)
    
    if (.tdr.obj$config$assay.type == "RNA") {
      
      # RNA: select half from top positive, half from top negative loadings
      n.select <-
        min(.n.features, length(x = loadings))
      
      n.positive <- ceiling(x = n.select / 2)
      n.negative <- floor(x = n.select / 2)
      
      top.features <- names(x = loadings)[c(
        head(x = feat.signed.order, n = n.positive),
        tail(x = feat.signed.order, n = n.negative)
      )]
      
    } else {
      
      # Cytometry: use all markers (sorted by signed loading)
      top.features <-
        names(x = loadings)[feat.signed.order]
      
    }
    
    # Order selected features by signed loading (high to low)
    top.loadings <-
      loadings[top.features]
    
    top.features <-
      top.features[order(top.loadings, decreasing = TRUE)]
    
    # -------------------------------------------------------------------------
    # Expression matrix: normalize and subset
    # -------------------------------------------------------------------------
    
    if (.tdr.obj$config$assay.type == "RNA") {
      
      # Step 1: Get row (landmark) sums from full matrix BEFORE subsetting
      lm.libsize <-
        Matrix::rowSums(x = .tdr.obj$raw.landmarks)
      
      # Step 2: Subset to top features (raw.landmarks is landmarks x features)
      X.sub <-
        .tdr.obj$raw.landmarks[, top.features, drop = FALSE]
      
      # Step 3: Size factor normalization using pre-computed library sizes
      # Size factor = libsize / mean(libsize)
      size.factors <-
        lm.libsize / mean(x = lm.libsize)
      
      X.norm <-
        X.sub / size.factors
      
      # Step 4: Log2 transform
      X.log <-
        log2(x = as.matrix(x = X.norm) + 1)
      
      # Step 5: Center each feature (column) to mean 0
      X.centered <-
        scale(x = X.log, center = TRUE, scale = FALSE)
      
    } else {
      
      # Cytometry: use raw landmarks directly, just center
      X.sub <-
        .tdr.obj$raw.landmarks[, top.features, drop = FALSE]
      
      X.centered <-
        scale(x = as.matrix(x = X.sub), center = TRUE, scale = FALSE)
      
    }
    
    # Transpose for heatmap: features x landmarks -> we want features as rows
    # X.centered is landmarks x features, so we transpose
    expr.mat <-
      t(x = X.centered)
    
    # -------------------------------------------------------------------------
    # Landmark ordering
    # -------------------------------------------------------------------------
    
    .user.order.by <- !is.null(x = .order.by)
    .user.add.annot <- !is.null(x = .add.annot)
    
    if (is.null(x = .order.by)) {
      
      # Default: use Y
      .order.by <-
        matrix(data = specDE.res$Y,
               ncol = 1,
               dimnames = list(names(x = specDE.res$Y),
                               "density\ncontrast\n(centered)"))
      
    } else {
      
      # Ensure matrix format with column names
      if (!inherits(x = .order.by, what = "matrix")) {
        stop(".order.by must be a matrix")
      }
      
      if (is.null(x = colnames(x = .order.by))) {
        stop(".order.by must be a matrix with column names.")
      }
      
      if (nrow(x = .order.by) != ncol(x = expr.mat)) {
        stop(".order.by must have ", ncol(x = expr.mat), " rows (one per landmark), ",
             "but has ", nrow(x = .order.by))
      }
      
    }
    
    # Validate .add.annot (does NOT affect ordering)
    if (isTRUE(x = .user.add.annot)) {
      
      if (!inherits(x = .add.annot, what = "matrix")) {
        stop(".add.annot must be a matrix")
      }
      
      if (is.null(x = colnames(x = .add.annot))) {
        stop(".add.annot must be a matrix with column names.")
      }
      
      if (nrow(x = .add.annot) != ncol(x = expr.mat)) {
        stop(".add.annot must have ", ncol(x = expr.mat), " rows (one per landmark), ",
             "but has ", nrow(x = .add.annot))
      }
      
    }
    
    # Sort landmarks by .order.by columns sequentially
    # Create a data frame for dplyr::arrange-style ordering
    order.df <-
      as.data.frame(x = .order.by)
    
    order.df$.idx <-
      seq_len(length.out = nrow(x = order.df))
    
    # Sort by all columns (sequential tie-breaking)
    if (.order.decreasing) {
      # Decreasing: multiply by -1 for numeric sorting
      for (col in colnames(x = .order.by)) {
        order.df[[col]] <- -order.df[[col]]
      }
    }
    
    order.df <-
      order.df[do.call(what = order,
                       args = order.df[, colnames(x = .order.by), drop = FALSE]), ]
    
    lm.order <-
      order.df$.idx
    
    # Reorder expression matrix columns
    expr.mat <-
      expr.mat[, lm.order, drop = FALSE]
    
    # Reorder .order.by for annotation
    .order.by <-
      .order.by[lm.order, , drop = FALSE]
    
    # Reorder .add.annot for annotation (does NOT affect ordering)
    if (isTRUE(x = .user.add.annot)) {
      .add.annot <- .add.annot[lm.order, , drop = FALSE]
    }
    
    # -------------------------------------------------------------------------
    # Build heatmap data frame for geom_raster
    # -------------------------------------------------------------------------
    
    n.landmarks <-
      ncol(x = expr.mat)
    
    n.features.plot <-
      nrow(x = expr.mat)
    
    # Create unique string labels for landmarks to avoid issues with numeric names
    # Use positional indices prefixed with "lm_" to ensure proper factor ordering
    lm.labels <-
      paste0("lm_", sprintf(fmt = "%05d", seq_len(length.out = n.landmarks)))
    
    # Build data frame manually to guarantee correct alignment
    # Features as rows, landmarks as columns in expr.mat
    heat.df <-
      data.frame(
        landmark = rep(x = lm.labels, times = n.features.plot),
        feature = rep(x = rownames(x = expr.mat), each = n.landmarks),
        expr = as.vector(x = t(x = expr.mat)),
        stringsAsFactors = FALSE
      )
    
    # Set factor levels to control plot ordering
    heat.df$landmark <-
      factor(x = heat.df$landmark,
             levels = lm.labels)
    
    # Features: reverse so highest loading is at top of heatmap (y-axis goes bottom to top)
    heat.df$feature <-
      factor(x = heat.df$feature,
             levels = rev(x = rownames(x = expr.mat)))
    
    # -------------------------------------------------------------------------
    # Build annotation strips: add.annot, order.by, specDE scores, density contrast
    # -------------------------------------------------------------------------
    
    # Get Y (density contrast) - always centered
    Y.ordered <- specDE.res$Y[lm.order]
    
    # Get specDE scores for selected dimensions, reordered
    specDE.scores.ordered <-
      specDE.res$scores[lm.order, .specDE.dim, drop = FALSE]
    
    # Build annotation data: .add.annot (top), then .order.by, then specDE dims, then density
    if (isTRUE(x = .user.add.annot) && isTRUE(x = .user.order.by)) {
      
      annot.names <- c(
        colnames(x = .order.by),
        "density\ncontrast\n(centered)",
        colnames(x = .add.annot),
        colnames(x = specDE.scores.ordered)
      )
      annot.mat <- cbind(.order.by, Y.ordered, .add.annot, specDE.scores.ordered)
      
    } else if (isTRUE(x = .user.add.annot) && !isTRUE(x = .user.order.by)) {
      
      annot.names <- c(
        "density\ncontrast\n(centered)",
        colnames(x = .add.annot),
        colnames(x = specDE.scores.ordered)
      )
      annot.mat <- cbind(Y.ordered, .add.annot, specDE.scores.ordered)
      
    } else if (!isTRUE(x = .user.add.annot) && isTRUE(x = .user.order.by)) {
      
      annot.names <- c(
        colnames(x = .order.by),
        "density\ncontrast\n(centered)",
        colnames(x = specDE.scores.ordered)
      )
      annot.mat <- cbind(.order.by, Y.ordered, specDE.scores.ordered)
      
    } else {
      
      annot.names <- c("density\ncontrast\n(centered)", colnames(x = specDE.scores.ordered))
      annot.mat <- cbind(Y.ordered , specDE.scores.ordered)
      
    }
    
    colnames(x = annot.mat) <- annot.names
    n.annot <- ncol(x = annot.mat)
    
    # Cycle viridis options for annotations
    .viridis.options.annot <-
      rep_len(x = .viridis.options.annot,
              length.out = n.annot)
    
    annot.df.list <-
      lapply(X = seq_len(length.out = n.annot),
             FUN = function(d) {
               data.frame(
                 landmark = lm.labels,
                 annot_dim = annot.names[d],
                 annot_value = annot.mat[, d],
                 stringsAsFactors = FALSE
               )
             })
    
    annot.df <-
      do.call(what = rbind,
              args = annot.df.list)
    
    # Set factor levels to ensure same ordering as expression heatmap
    annot.df$landmark <-
      factor(x = annot.df$landmark,
             levels = lm.labels)
    
    annot.df$annot_dim <-
      factor(x = annot.df$annot_dim,
             levels = annot.names)
    
    # row annotation data (signed loadings)
    row.annot.df <- 
      data.frame(
        feature = rownames(x = expr.mat) |> 
          (\(x)
           factor(x = x,
                  levels = rev(x = x))
          )(),
        loading = top.loadings[rownames(expr.mat)],
        x = 1
      )
    
    # -------------------------------------------------------------------------
    # Build ggplot: expression heatmap
    # -------------------------------------------------------------------------
    
    p.heat <-
      ggplot2::ggplot(data = heat.df,
                      mapping = ggplot2::aes(x = landmark,
                                             y = feature,
                                             fill = expr)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradient2(low = unname(obj = Color.Palette[1, 1]),
                                    mid = unname(obj = Color.Palette[1, 6]),
                                    high = unname(obj = Color.Palette[1, 2]),
                                    midpoint = 0,
                                    name = "Expression\n(centered)",
                                    guide = ggplot2::guide_colorbar(
                                      direction = "vertical",
                                      barwidth = grid::unit(x = 0.15, units = "in"),
                                      barheight = grid::unit(x = 0.5, units = "in"),
                                      title.position = "top",
                                      title.hjust = 0
                                    )) +
      ggplot2::labs(
        x = "Landmarks",
        y = "Features"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = I(x = .feature.font.size)),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 7),
        legend.text = ggplot2::element_text(size = 6)
      )
    
    
    # row annotation plot (right side)
    p.row.annot <-
      ggplot2::ggplot(row.annot.df,
                      ggplot2::aes(x = x, 
                                   y = feature, 
                                   fill = loading)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_viridis_c(
        option = "turbo",
        name = "loadings",
        guide = ggplot2::guide_colorbar(
          direction = "vertical",
          barwidth = grid::unit(x = 0.15, units = "in"),
          barheight = grid::unit(x = 0.5, units = "in"),
          title.position = "top",
          title.hjust = 0
        )) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right",
                     legend.title = ggplot2::element_text(size = 7),
                     legend.text = ggplot2::element_text(size = 6))
    
    # Hide landmark labels if requested (usually too many)
    if (!.show.landmark.labels) {
      p.heat <-
        p.heat +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    } else {
      p.heat <-
        p.heat +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.5,
                                              size = I(x = 6))
        )
    }
    
    # -------------------------------------------------------------------------
    # Build annotation strips (one per annotation)
    # -------------------------------------------------------------------------
    
    annot.plots <-
      seq_len(length.out = n.annot) |>
      lapply(FUN = function(d) {
        
        dim.name <-
          annot.names[d]
        
        dim.df <-
          annot.df[annot.df$annot_dim == dim.name, ]
        
        .option.this <-
          if (identical(x = dim.name, y = "density\ncontrast\n(centered)")) {
            "viridis"
          } else if (dim.name %in% colnames(x = specDE.scores.ordered)) {
            "plasma"
          } else {
            .viridis.options.annot[d]
          }
        
        ggplot2::ggplot(data = dim.df,
                        mapping = ggplot2::aes(x = landmark,
                                               y = 1,
                                               fill = annot_value)) +
          ggplot2::geom_raster() +
          ggplot2::scale_fill_viridis_c(option = .option.this,
                                        name = dim.name,
                                        guide = ggplot2::guide_colorbar(
                                          direction = "vertical",
                                          barwidth = grid::unit(x = 0.15, units = "in"),
                                          barheight = grid::unit(x = 0.5, units = "in"),
                                          title.position = "top",
                                          title.hjust = 0
                                        )) +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::theme_void() +
          ggplot2::theme(
            legend.position = "right",
            legend.title = ggplot2::element_text(size = 7),
            legend.text = ggplot2::element_text(size = 6)
          )
      })
    
    # -------------------------------------------------------------------------
    # Assemble via gtable for deterministic sizing
    # -------------------------------------------------------------------------
    
    # --- Helper: force panel size in a gtable ---
    .force.panel <- function(gt, w, h) {
      panel.pos <- gt$layout[gt$layout$name == "panel", ]
      gt$widths[panel.pos$l]  <- grid::unit(x = w, units = "in")
      gt$heights[panel.pos$t] <- grid::unit(x = h, units = "in")
      gt
    }
    
    # --- Helper: extract legend from gtable, return list(gt.no.leg, legend) ---
    .extract.legend <- function(gt) {
      leg.idx <- grep(pattern = "guide-box", x = gt$layout$name)
      if (length(x = leg.idx) == 0) {
        return(list(body = gt, legend = grid::nullGrob()))
      }
      legend <- gt$grobs[[leg.idx[1]]]
      gt$grobs[[leg.idx[1]]] <- grid::nullGrob()
      list(body = gt, legend = legend)
    }
    
    # --- Convert all plots to gtables and force panel sizes ---
    
    # 1. Heatmap
    g.heat <- ggplot2::ggplotGrob(x = p.heat)
    g.heat <- .force.panel(gt = g.heat,
                           w = .panel.width,
                           h = .panel.height)
    heat.parts <- .extract.legend(gt = g.heat)
    g.heat.body <- heat.parts$body
    leg.heat <- heat.parts$legend
    
    # 2. Row annotation
    g.row <- ggplot2::ggplotGrob(x = p.row.annot)
    g.row <- .force.panel(gt = g.row,
                          w = .annot.panel.height,
                          h = .panel.height)
    row.parts <- .extract.legend(gt = g.row)
    g.row.body <- row.parts$body
    leg.row <- row.parts$legend
    
    # 3. Column annotation strips
    g.annots <- lapply(X = annot.plots, FUN = function(p) {
      gt <- ggplot2::ggplotGrob(x = p)
      gt <- .force.panel(gt = gt,
                         w = .panel.width,
                         h = .annot.panel.height)
      .extract.legend(gt = gt)
    })
    g.annot.bodies <- lapply(X = g.annots, FUN = `[[`, "body")
    leg.annots <- lapply(X = g.annots, FUN = `[[`, "legend")
    
    # --- Align widths for vertically-stacked panels (annotation strips + heatmap) ---
    all.col.aligned <- c(g.annot.bodies, list(g.heat.body))
    max.widths <- do.call(what = grid::unit.pmax,
                          args = lapply(X = all.col.aligned,
                                        FUN = function(g) g$widths))
    g.annot.bodies <- lapply(X = g.annot.bodies, FUN = function(g) {
      g$widths <- max.widths
      g
    })
    g.heat.body$widths <- max.widths
    
    # --- Align heights for horizontally-adjacent panels (row annot + heatmap) ---
    max.heights <- grid::unit.pmax(g.row.body$heights, g.heat.body$heights)
    g.row.body$heights <- max.heights
    g.heat.body$heights <- max.heights
    
    # --- Zero out inner margins so heatmap + row annotation sit flush ---
    # Remove right margin of heatmap (columns to the right of panel)
    heat.panel.col <- g.heat.body$layout[g.heat.body$layout$name == "panel", "l"]
    if (heat.panel.col < ncol(x = g.heat.body)) {
      for (col in (heat.panel.col + 1):ncol(x = g.heat.body)) {
        g.heat.body$widths[col] <- grid::unit(x = 0, units = "pt")
      }
    }
    # Remove left margin of row annotation (columns to the left of panel)
    row.panel.col <- g.row.body$layout[g.row.body$layout$name == "panel", "l"]
    if (row.panel.col > 1) {
      for (col in seq_len(length.out = row.panel.col - 1)) {
        g.row.body$widths[col] <- grid::unit(x = 0, units = "pt")
      }
    }
    # Remove right margin of row annotation (to keep legends close)
    if (row.panel.col < ncol(x = g.row.body)) {
      for (col in (row.panel.col + 1):ncol(x = g.row.body)) {
        g.row.body$widths[col] <- grid::unit(x = 0, units = "pt")
      }
    }
    
    # Also zero out right-side columns of annotation strips to match
    for (j in seq_along(along.with = g.annot.bodies)) {
      g <- g.annot.bodies[[j]]
      annot.panel.col <- g$layout[g$layout$name == "panel", "l"]
      if (annot.panel.col < ncol(x = g)) {
        for (col in (annot.panel.col + 1):ncol(x = g)) {
          g$widths[col] <- grid::unit(x = 0, units = "pt")
        }
      }
      g.annot.bodies[[j]] <- g
    }
    
    # --- Build bottom row: heatmap | row annotation ---
    g.bottom <- cbind(g.heat.body, g.row.body, size = "first")
    
    # --- Stack: annotation strips on top, bottom row below ---
    # Pad annotation strip widths to match bottom row
    g.annot.bodies <- lapply(X = g.annot.bodies, FUN = function(g) {
      # Add empty columns on the right to account for g.row.body width
      n.row.cols <- ncol(x = g.row.body)
      for (i in seq_len(length.out = n.row.cols)) {
        g <- gtable::gtable_add_cols(x = g,
                                     widths = g.row.body$widths[i])
      }
      g
    })
    
    # Align all widths across annotation strips + bottom row for final rbind
    all.to.stack <- c(g.annot.bodies, list(g.bottom))
    final.max.widths <- do.call(what = grid::unit.pmax,
                                args = lapply(X = all.to.stack,
                                              FUN = function(g) g$widths))
    all.to.stack <- lapply(X = all.to.stack, FUN = function(g) {
      g$widths <- final.max.widths
      g
    })
    
    # Stack vertically
    g.final <- all.to.stack[[1]]
    if (length(x = all.to.stack) > 1) {
      for (i in 2:length(x = all.to.stack)) {
        g.final <- rbind(g.final, all.to.stack[[i]], size = "max")
      }
    }
    
    # --- Assemble legends column ---
    all.legends <- c(leg.annots, list(leg.row), list(leg.heat))
    # Filter out nullGrobs
    keep <- vapply(X = all.legends,
                   FUN = function(lg) !inherits(x = lg, what = "nullGrob"),
                   FUN.VALUE = logical(length = 1))
    all.legends <- all.legends[keep]
    
    if (length(x = all.legends) > 0) {
      # Build heights: legend height + spacer between each pair
      spacer.h <- grid::unit(x = 0.15, units = "in")
      height.list <- vector(mode = "list")
      for (i in seq_along(along.with = all.legends)) {
        if (i > 1) {
          height.list <- c(height.list, list(spacer.h))
        }
        h <- if (inherits(x = all.legends[[i]], what = "gtable")) {
          sum(all.legends[[i]]$heights)
        } else {
          grid::unit(x = 0.5, units = "in")
        }
        height.list <- c(height.list, list(h))
      }
      
      leg.column <- gtable::gtable(
        widths = grid::unit(x = 1, units = "in"),
        heights = do.call(what = grid::unit.c, args = height.list)
      )
      
      # Place legends in odd rows (1, 3, 5, ...), spacers in even rows
      row.idx <- seq(from = 1, by = 2, length.out = length(x = all.legends))
      for (i in seq_along(along.with = all.legends)) {
        leg.column <- gtable::gtable_add_grob(x = leg.column,
                                              grobs = all.legends[[i]],
                                              t = row.idx[i], l = 1)
      }
      
      # Combine: body | legends
      # Add spacer
      g.final <- gtable::gtable_add_cols(x = g.final,
                                         widths = grid::unit(x = 0.1, units = "in"))
      g.final <- gtable::gtable_add_cols(x = g.final,
                                         widths = grid::unit(x = 1, units = "in"))
      g.final <- gtable::gtable_add_grob(x = g.final,
                                         grobs = leg.column,
                                         t = 1,
                                         l = ncol(x = g.final),
                                         b = nrow(x = g.final))
    }
    
    # Return g.final as ggplot object
    return(ggplot2::ggplot() +
             ggplot2::annotation_custom(grob = g.final) +
             ggplot2::theme_void())
    
  }

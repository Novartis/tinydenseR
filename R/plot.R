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

#' Plot landmark feature discovery summary
#'
#' Creates an interactive scatter plot showing landmarks in reduced dimensions (UMAP, PCA, or 
#' custom features) with hover tooltips displaying the top characteristic features for each 
#' landmark. Requires \code{get.lm.features.stats()} to have been run first.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.lm.features.stats()}.
#' @param .n.feat Integer number of top features to display in hover tooltip (default 10). 
#'   Currently not implemented - displays all features from \code{$interact.plot$lm.features}.
#' @param .plot.dims Character specifying plot dimensions. Options:
#'   \itemize{
#'     \item \code{"umap"} - UMAP coordinates (default)
#'     \item \code{"pca"} - PCA coordinates (specify PCs with \code{.PC.x}, \code{.PC.y})
#'     \item Vector of 2 feature names from \code{colnames(.lm.obj$lm)} for custom biplot
#'   }
#' @param .PC.x Integer specifying x-axis PC when \code{.plot.dims = "pca"} (default 1).
#' @param .PC.y Integer specifying y-axis PC when \code{.plot.dims = "pca"} (default 2).
#' @param .panel.size Numeric panel width/height in inches (default 1.5).
#' @param .point.size Numeric point size (default 2).
#'   
#' @return If \pkg{ggiraph} is available, returns an interactive \code{girafe} object with hover 
#'   tooltips. Otherwise returns a static \code{ggplot} object with a warning.
#'   
#' @note Interactive features require the \pkg{ggiraph} package. Install with 
#'   \code{install.packages("ggiraph")} for full functionality.
#'   
#' @seealso \code{\link{get.lm.features.stats}} for computing feature signatures
#' 
#' @examples
#' \dontrun{
#' # Typical workflow for interactive feature exploration
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500) |>
#'   get.graph() |>
#'   get.lm.features.stats()
#' 
#' # Interactive UMAP with feature tooltips
#' interactFeatPlot(lm.cells)
#' 
#' # PCA biplot with PC2 vs PC3
#' interactFeatPlot(lm.cells, .plot.dims = "pca", .PC.x = 2, .PC.y = 3)
#' 
#' # Custom feature biplot (cytometry example)
#' interactFeatPlot(lm.cells, .plot.dims = c("CD3", "CD4"))
#' }
#' 
#' @export
interactFeatPlot <-
  function(
    .lm.obj,
    .n.feat = 10,
    .plot.dims = "umap",
    .PC.x = 1,
    .PC.y = 2,
    .panel.size = 1.5,
    .point.size = 2
  ){
    
    # R CMD check appeasement
    .plot.x <- .plot.y <- topFeatTab <- NULL
    
    if(length(x = .plot.dims) > 2){
      stop(".plot.dims must be 'umap', 'pca', or a vector of 2 feature names.\n",
           "Current length: ", length(x = .plot.dims))
    }
    
    # Validate plot dimensions
    .plot.dims <-
      match.arg(arg = .plot.dims,
                choices = c(
                  "umap",
                  "pca",
                  colnames(x = .lm.obj$lm)
                ),
                several.ok = TRUE)
    
    # Extract x-axis coordinates based on plot type
    .plot.x <-
      if(.plot.dims == "umap"){
        .lm.obj$graph$uwot$embedding[,1]
      } else if(.plot.dims == "pca"){
        .lm.obj$pca$embed[,.PC.x]
      } else {
        .lm.obj$lm[,.plot.dims[1]]
      }
    
    # Extract y-axis coordinates based on plot type
    .plot.y <-
      if(.plot.dims == "umap"){
        .lm.obj$graph$uwot$embedding[,2]
      } else if(.plot.dims == "pca"){
        .lm.obj$pca$embed[,.PC.y]
      } else {
        .lm.obj$lm[,.plot.dims[2]]
      }
    
    # Build data frame for plotting
    dat.df <-
      cbind(.plot.x,
            .plot.y) |>
      as.data.frame()
    
    # Set column names to match plot type
    colnames(x = dat.df) <-
      if(.plot.dims == "umap"){
        paste0("umap.",1:2)
      } else if(.plot.dims == "pca"){
        paste0("PC",c(.PC.x,.PC.y))
      } else {
        .plot.dims
      }
    
    # Attach HTML tooltips with feature signatures
    dat.df$topFeatTab <-
      .lm.obj$interact.plot$lm.features$html
    
    # Create base plot
    p <-
      ggplot2::ggplot(data = dat.df,
                      mapping = ggplot2::aes(x = !!rlang::sym(colnames(x = dat.df)[1]),
                                             y = !!rlang::sym(colnames(x = dat.df)[2]),
                                             tooltip = topFeatTab)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "top marker",
                    subtitle = "stats (hover)") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
    # Add interactivity if ggiraph is available
    if (.is_ggiraph_available()) {
      p <- p + ggiraph::geom_point_interactive(size = I(x = .point.size))
      
      # Return interactive plot with zoom
      ggiraph::girafe(ggobj = p) |>
        ggiraph::girafe_options(
          ggiraph::opts_zoom(max = 10,
                             min = 1)
        )
    } else {
      # Fallback to static plot if ggiraph not installed
      warning("ggiraph package not available. Install with install.packages('ggiraph') for interactive features.\n",
              "Returning static plot instead.", 
              call. = FALSE)
      p + ggplot2::geom_point(size = I(x = .point.size))
    }
  }

#' Plot PCA
#'
#' Visualizes landmarks in PCA space with flexible coloring by features, clusters, or statistical 
#' results. Supports both continuous (e.g., gene expression, fold changes) and categorical 
#' (e.g., clusters, cell types) overlays. Optional interactive hover shows landmark feature signatures.
#'
#' @param .lm.obj A tinydenseR object with PCA computed via \code{get.landmarks()}.
#' @param .PC.x Integer specifying x-axis principal component (default 1).
#' @param .PC.y Integer specifying y-axis principal component (default 2).
#' @param .feature Numeric vector or factor of length \code{nrow(.lm.obj$lm)} to color points by. 
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
#' @param .hover.n.features Integer number of features per direction in hover (default 10). 
#'   Currently not implemented.
#'   
#' @return A \code{ggplot} object (static) or \code{girafe} object (interactive, if \pkg{ggiraph} 
#'   available and \code{.hover.stats != "none"}).
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
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks(.nHVG = 500, .nPC = 3) |>
#'   get.graph()
#' 
#' plotPCA(lm.cells, .point.size = 1, .panel.size = 1.5)
#' 
#' # Color by fold change from get.stats()
#' condition.stats <- get.stats(lm.cells, .design = design)
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
  function(.lm.obj,
           .PC.x = 1,
           .PC.y = 2,
           .feature = .lm.obj$graph$clustering$ids,
           .cat.feature.color = Color.Palette[1,1:5],
           .panel.size = if(is.numeric(x = .feature)) 2 else 3,
           .midpoint = NULL,
           .plot.title = "",
           .color.label = "",
           .legend.position = "right",
           .point.size = 0.1,
           .seed = 123,
           .hover.stats = "none",
           .hover.n.features = 10){
    
    # R CMD check appeasement
    topFeatTab <- feature <- NULL
    
    # Validate PC selection
    if(any(!c(.PC.x, .PC.y) %in%
           1:ncol(x = .lm.obj$pca$embed))){
      
      stop(".PC.x and .PC.y must be valid PC indices.\n",
           "Available PCs: ", paste(x = 1:ncol(x = .lm.obj$pca$embed), collapse = ", "),
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
      colnames(x = .lm.obj$pca$embed)[.PC.x]
    
    .PC.y <-
      colnames(x = .lm.obj$pca$embed)[.PC.y]
    
    # Set midpoint to median for continuous features
    if(is.null(x = .midpoint) &
       is.numeric(x = .feature)){
      .midpoint <- stats::median(x = .feature)
    }
    
    # Build plotting data frame
    set.seed(seed = .seed)
    dat.df <-
      as.data.frame(x = .lm.obj$pca$embed[,c(.PC.x,.PC.y)]) |>
      cbind(feature = .feature,
            point.size = .point.size)
    
    # Attach hover tooltips if requested
    if(.hover.stats == "marker"){
      
      dat.df$topFeatTab <-
        .lm.obj$interact.plot$lm.features$html
      
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
#' @param .lm.obj A tinydenseR object with UMAP computed via \code{get.graph()}.
#' @param .feature Numeric vector or factor of length \code{nrow(.lm.obj$lm)} to color points by. 
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
#' @param .hover.n.features Integer number of features per direction in hover (default 10). 
#'   Currently not implemented.
#'   
#' @return A \code{ggplot} object (static) or \code{girafe} object (interactive, if \pkg{ggiraph} 
#'   available and \code{.hover.stats != "none"}).
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
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph()
#' 
#' # Basic UMAP
#' plotUMAP(lm.cells, .point.size = 1)
#' 
#' # Color by gene expression
#' plotUMAP(lm.cells, .feature = lm.cells$lm[,"CD4"], .color.label = "CD4")
#' 
#' # Color by fold change from get.stats()
#' condition.stats <- get.stats(lm.cells, .design = design)
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
  function(.lm.obj,
           .feature = .lm.obj$graph$clustering$ids,
           .cat.feature.color = Color.Palette[1,1:5],
           .panel.size = 2,
           .midpoint = NULL,
           .plot.title = "",
           .color.label = "",
           .legend.position = "right",
           .point.size = 0.1,
           .seed = 123,
           .hover.stats = "none",
           .hover.n.features = 10){
    
    # R CMD check appeasement
    umap.1 <- umap.2 <- topFeatTab <- feature <- NULL
    
    if(is.null(x = .lm.obj$graph)){
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
      as.data.frame(x = .lm.obj$graph$uwot$embedding) |>
      cbind(feature = .feature,
            point.size = .point.size)
    
    # Attach hover tooltips if requested
    if(.hover.stats == "marker"){
      
      dat.df$topFeatTab <-
        .lm.obj$interact.plot$lm.features$html
      
    }
    
    # Randomize row order to avoid plotting bias
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

#' Bee swarm Plot of Density Estimate Change
#'
#' Creates a beeswarm plot showing effect sizes (log fold changes) from differential abundance 
#' testing, with points colored by significance. Optionally splits results by cluster or cell type 
#' and displays mean cell percentages alongside for biological context.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .stats.obj Statistical results from \code{get.stats()}.
#' @param .coefs Character vector of coefficient names from the design matrix to plot. Must match 
#'   column names in \code{.stats.obj$fit$coefficients}. Can plot single or multiple coefficients.
#' @param .q Numeric q-value threshold for significance coloring (default 0.1). Points with 
#'   q < threshold are colored by direction (red/blue), otherwise gray.
#' @param .q.from Character specifying q-value source: "pca.weighted.q" (default, PCA-variance weighted) 
#'   or "density.weighted.bh.fdr" (density weighted). See \code{get.stats()} for details.
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
#'   
#' @return A \code{ggplot} object (single plot) or combined plot using \code{patchwork} (when 
#'   \code{.perc.plot = TRUE} and \code{.split.by != "none"}).
#'   
#' @details
#' The beeswarm layout prevents overlapping points, making it easier to assess the distribution 
#' of effect sizes across landmarks. Red indicates significant increases (log FC > 0), blue indicates 
#' significant decreases, and gray indicates non-significant changes.
#' 
#' When split by clustering/celltyping, the percentage plot shows mean cell type abundances across 
#' samples, helping interpret whether changes occur in rare or common populations.
#' 
#' @seealso \code{\link{get.stats}}, \code{\link{plotAbundance}}
#' 
#' @examples
#' \dontrun{
#' # After statistical testing
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.map()
#' 
#' design <- model.matrix(~ Condition + Replicate, data = .meta)
#' stats <- get.stats(lm.cells, .design = design)
#' 
#' # Basic beeswarm split by clusters
#' plotBeeswarm(lm.cells, stats, .coefs = "ConditionB")
#' 
#' # Multiple coefficients without splitting
#' plotBeeswarm(lm.cells, stats, 
#'              .coefs = c("ConditionB", "ConditionC"),
#'              .split.by = "none")
#' }
#' 
#' @export
plotBeeswarm <-
  function(
    .lm.obj,
    .stats.obj,
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
    .perc.plot = TRUE) {
    
    sig <- adj.p <- facets <- pos.t <- neg.t <- q.bars <- dat.df <- value <- split.by <- .coef <- cell.perc <- pop <- NULL
    
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
                   pop = colnames(x = .lm.obj$map[[.split.by]]$cell.perc),
                   cell.perc = Matrix::colMeans(x = .lm.obj$map[[.split.by]]$cell.perc)) |>
        (\(x)
         ggplot2::ggplot(data = x,
                         mapping = ggplot2::aes(x = x,
                                                y = pop,
                                                size = cell.perc)) +
           ggplot2::scale_y_discrete(limits = rev) +
           ggplot2::theme_minimal() +
           ggplot2::theme(legend.margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
                          plot.title = ggplot2::element_text(hjust = 0.5),
                          legend.position = "bottom") +
           ggplot2::labs(title = "percentages",
                         x = "",
                         y = "",
                         size = "%") +
           ggplot2::scale_x_discrete(labels = "mean % across\nall samples") +
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
        as.data.frame(x = .stats.obj$fit$coefficients[!(.lm.obj$graph[[.split.by]]$ids %in% 
                                                          .lm.obj$map$cl.ct.to.ign),.coefs]) |>
        tidyr::pivot_longer(cols = tidyselect::everything(),
                            names_to = "split.by",
                            values_to = "value",
                            cols_vary = "slowest") |>
        dplyr::bind_cols(sig = as.data.frame(x = (.stats.obj$fit[[.q.from]][!(.lm.obj$graph[[.split.by]]$ids %in% 
                                                                                .lm.obj$map$cl.ct.to.ign),.coefs] < .q)) |>
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
                        split.by = rep(x = .lm.obj$graph[[.split.by]]$ids[!(.lm.obj$graph[[.split.by]]$ids %in% 
                                                                              .lm.obj$map$cl.ct.to.ign)] |>
                                         as.character(),
                                       times = length(x = .coefs)))
        
      }
      
    } else if(.split.by == "clustering"){
      
      dat.df <-
        data.frame(
          value =
            .stats.obj$fit$coefficients[!(.lm.obj$graph$clustering$ids %in% 
                                            .lm.obj$map$cl.ct.to.ign),
                                        .coefs],
          sig =
            ifelse(test = .stats.obj$fit[[.q.from]][,.coefs] < .q,
                   yes = ifelse(test = .stats.obj$fit$coefficients[,.coefs] > 0,
                                yes = "more abundant",
                                no = "less abundant"),
                   no = "not sig.")[!(.lm.obj$graph$clustering$ids %in% 
                                        .lm.obj$map$cl.ct.to.ign)],
          split.by =
            .lm.obj$graph$clustering$ids[!(.lm.obj$graph$clustering$ids %in% 
                                             .lm.obj$map$cl.ct.to.ign)] |>
            as.character()) |>
        droplevels()
      
    } else if(.split.by == "celltyping"){
      
      if(is.null(x = .lm.obj$graph$celltyping$ids)){
        
        stop(".lm.obj$graph$celltyping$ids could not be found")
        
      }
      
      dat.df <-
        data.frame(
          value =
            .stats.obj$fit$coefficients[!(.lm.obj$graph$celltyping$ids %in% 
                                            .lm.obj$map$cl.ct.to.ign),.coefs],
          sig =
            ifelse(test = .stats.obj$fit[[.q.from]][,.coefs] < .q,
                   yes = ifelse(test = .stats.obj$fit$coefficients[,.coefs] > 0,
                                yes = "more abundant",
                                no = "less abundant"),
                   no = "not sig.")[!(.lm.obj$graph$celltyping$ids %in% 
                                        .lm.obj$map$cl.ct.to.ign)],
          split.by =
            .lm.obj$graph$celltyping$ids[!(.lm.obj$graph$celltyping$ids %in% 
                                             .lm.obj$map$cl.ct.to.ign)] |>
            as.character()) |>
        droplevels()
      
    }
    
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
           x = "abundance\nlog2(+0.5)FC",
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
                              draw_quantiles = 0.5) +
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

#' Bidimensional plots for markers expression visualization
#'
#' This function plots the expression of two markers in a bidimensional plot.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks and get.graph.
#' @param .id The id of the cluster or celltype to plot. If NULL, all landmarks are plotted.
#' @param .id.from The source of the id. Can be set to "clustering" or "celltyping". Defaults to "clustering".
#' @param .x.feature The x-axis marker. Defaults to "CD3".
#' @param .y.feature The y-axis marker. Defaults to "CD20".
#' @param .bins The number of bins for the hexagonal plot. Defaults to 128.
#' @param .legend.position The position of the legend. Defaults to "right" and can be set to "left", "top", "bottom" or "none".
#' @param .plot.title The title of the plot.
#' @param .panel.size The size of the panel.
#' @param .reference If TRUE, the plot will show the density of the landmarks. Defaults to TRUE.
#' @param .density.bins The number of bins for the density plot. Defaults to 32.
#' @param .sd.range The range of the standard deviation outside which cells are considered outliers and excluded from plot. Defaults to c(-3, 6).
#'
#' @export
plot2Markers <-
  function(.lm.obj,
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
        
        .id <- .lm.obj$graph$clustering$ids == .id
        
      } else if(.id.from == "celltyping"){
        
        if(is.null(x = .lm.obj$graph$celltyping$ids)){
          stop(".lm.obj$graph$celltyping$ids could not be found")
        }
        .id <- .lm.obj$graph$celltyping$ids == .id
      }
    }
    
    if(.lm.obj$assay.type == "RNA"){
      
      dat.df <-
        (.lm.obj$lm[if(!is.null(x = .id)) .id else 1:nrow(x = .lm.obj$lm),
                    colnames(x = .lm.obj$lm) %in% c(.x.feature,.y.feature)]) |>
        as.data.frame()
      
      
    } else {
      
      dat.df <-
        (.lm.obj$lm[if(!is.null(x = .id)) .id else 1:nrow(x = .lm.obj$lm),
                    colnames(x = .lm.obj$lm) %in% c(.x.feature,.y.feature)]/50) |>
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
    
    if(.lm.obj$assay.type != "RNA"){
      
      p <-
        p +
        ggplot2::scale_x_continuous(labels = scales::label_math(expr = 10^.x),
                                    limits = (.lm.obj$lm[,colnames(x = .lm.obj$lm) %in%
                                                           .x.feature, drop = TRUE]/50) |>
                                      (\(x)
                                       range(x[!((x > (mean(x = x) + .sd.range[2]*stats::sd(x = x))) |
                                                   (x < (mean(x = x) + .sd.range[1]*stats::sd(x = x))))])
                                      )()) +
        ggplot2::scale_y_continuous(labels = scales::label_math(expr = 10^.x),
                                    limits = (.lm.obj$lm[,colnames(x = .lm.obj$lm) %in%
                                                           .y.feature, drop = TRUE]/50) |>
                                      (\(x)
                                       range(x[!((x > (mean(x = x) + .sd.range[2]*stats::sd(x = x))) |
                                                   (x < (mean(x = x) + .sd.range[1]*stats::sd(x = x))))])
                                      )()) +
        ggplot2::annotation_logticks(base = 10)
      
    }
    
    if(!is.null(x = .id) &
       isTRUE(x = .reference)) {
      
      p <-
        p  +
        ggplot2::stat_density_2d(data = as.data.frame(x = .lm.obj$lm / (if(.lm.obj$assay.type == "RNA") 1 else 50)),
                                 mapping = ggplot2::aes(fill = ggplot2::after_stat(x = log2(x = level))),
                                 geom = "polygon",
                                 bins = .density.bins) +
        ggplot2::scale_fill_gradient(low = unname(obj = Color.Palette[2,7]),
                                     high = unname(obj = Color.Palette[6,7]),
                                     guide = "none")
      
    }
    
    #if(is.null(x = .color.from)){
    
    p <-
      p +
      ggnewscale::new_scale(new_aes = "fill") +
      ggplot2::geom_hex(bins = .bins) +
      ggplot2::scale_fill_gradient(low = unname(obj = Color.Palette[1,1]),
                                   high = unname(obj = Color.Palette[1,2]),
                                   trans = "log2")
    
    #} else {
    #
    #p <-
    #  p +
    #  ggplot2::geom_point(mapping = ggplot2::aes(color = t),
    #                      size = I(x = 0.1)) +
    #  ggplot2::scale_color_gradient2(low = unname(obj = Color.Palette[1,1]),
    #                                 mid = unname(obj = Color.Palette[1,6]),
    #                                 high = unname(obj = Color.Palette[1,2]),
    #                                 midpoint = 0)
    #
    #}
    
    p +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
  }

#' Plot Sample PCA
#'
#' Visualizes samples in PCA space based on their landmark density profiles. This "reverse PCA" 
#' treats samples as observations and landmark densities as features, revealing sample-level 
#' patterns and clustering. Useful for quality control and identifying batch effects.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .labels.from Character specifying metadata column for coloring points. Defaults to 
#'   first column of \code{.lm.obj$metadata} (often "Condition" or "Sample").
#' @param .cat.feature.color Character vector of colors for categorical labels (default 
#'   \code{Color.Palette[1,1:5]}).
#' @param .point.size Numeric point size (default 1).
#' @param .panel.size Numeric panel width/height in inches (default 2).
#' @param .midpoint Numeric midpoint for diverging color scale (continuous labels only). 
#'   Defaults to median of \code{.labels.from} column.
#' @param .seed Integer random seed (default 123).
#'   
#' @return A \code{ggplot} object showing PC1 vs PC2 of sample-level density profiles.
#'   
#' @details
#' This function performs PCA on log2-transformed landmark densities (samples × landmarks matrix). 
#' Samples with similar cellular composition will cluster together in PC space. The analysis:
#' \itemize{
#'   \item Log-transforms densities: log2(density + 0.5) for variance stabilization
#'   \item Filters zero-variance landmarks
#'   \item Centers and scales features before PCA
#'   \item Computes top 2 PCs using truncated SVD
#' }
#' 
#' Useful for identifying:
#' \itemize{
#'   \item Batch effects (samples clustering by technical rather than biological factors)
#'   \item Outlier samples with aberrant cellular profiles
#'   \item Main axes of variation in cellular composition
#' }
#' 
#' @seealso \code{\link{get.map}}, \code{\link{plotPCA}}
#' 
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
#'   get.landmarks() |>
#'   get.graph() |>
#'   get.map()
#' 
#' # Color by condition
#' plotSamplePCA(lm.cells, .labels.from = "Condition")
#' 
#' # Color by replicate to check for batch effects
#' plotSamplePCA(lm.cells, .labels.from = "Replicate")
#' }
#' 
#' @export
#'
plotSamplePCA <-
  function(
    .lm.obj,
    .labels.from = colnames(x = .lm.obj$metadata)[1],
    .cat.feature.color = Color.Palette[1,1:5],
    .point.size = 1,
    .panel.size = 2,
    .midpoint = if(is.numeric(x = .lm.obj$metadata[[.labels.from]])) stats::median(x = .lm.obj$metadata[[.labels.from]]) else NA,
    .seed = 123){
    
    harmony_1 <- harmony_2 <- PC1 <- PC2 <- labels.from <- NULL
    
    if(length(x = .labels.from) != 1){
      stop(".labels.from must be length 1")
    }
    
    if(!(.labels.from %in% colnames(x = .lm.obj$metadata))){
      stop(paste0(".labels.from must be one of the following: ",
                  paste(x = colnames(x = .lm.obj$metadata),
                        collapse = ", ")))
    }
    
    if(is.null(x = .lm.obj$map$fdens)){
      stop("run get.map first.")
    }
    
    set.seed(seed = .seed)
    y.pca <-
      log2(x = .lm.obj$map$fdens + 0.5) |> 
      (\(x)
       x[matrixStats::rowVars(x = x) > 0,]
      )() |>
      Matrix::t() |>
      (\(x)
       irlba::prcomp_irlba(x = x,
                           center = Matrix::colMeans(x = x),
                           scale. = matrixStats::colSds(x = x),
                           rank. = 2))() |>
      (\(x)
       as.data.frame(x = x$x)
      )() |>
      cbind(labels.from = .lm.obj$metadata[[.labels.from]])
    
    p <-
      ggplot2::ggplot(data = y.pca,
                      mapping = ggplot2::aes(x = PC1,
                                             y = PC2,
                                             color = labels.from)) +
      ggplot2::theme_bw() +
      ggplot2::geom_point(size = I(x = .point.size)) +
      ggplot2::labs(color = .labels.from) +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in"),
                              col = grid::unit(x = .panel.size,
                                               units = "in"))
    
    if(is.numeric(x = y.pca$labels.from)){
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
          )(length(x = unique(x = y.pca$labels.from))))
    }
    
    return(p)
    
  }


#' Plot Traditional Statistics
#'
#' Visualizes results from traditional cluster/cell type-level differential abundance testing. 
#' Shows effect sizes as heatmap with significance markers, providing a complementary view to 
#' landmark-based analysis for easier interpretation at the population level.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .stats.obj Statistical results from \code{get.stats()} (must include traditional analysis).
#' @param .coefs Character vector of coefficient names to plot. Defaults to all coefficients from 
#'   traditional model.
#' @param .split.by Character: "clustering" (default) or "celltyping" - which population grouping to use.
#' @param .q Numeric q-value threshold for significance stars (default 0.1).
#' @param .row.space.scaler Numeric scaling for row height (default 0.2 inches per population).
#' @param .col.space.scaler Numeric scaling for column width (default 0.5 inches per coefficient).
#' @param .label.substr.rm Character substring to remove from labels (default "").
#'   
#' @return A \code{ggplot} heatmap with log fold changes colored by magnitude and significance 
#'   indicated by asterisks.
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
#' @seealso \code{\link{get.stats}}, \code{\link{plotTradPerc}}, \code{\link{plotBeeswarm}}
#' 
#' @examples
#' \dontrun{
#' # After get.stats with traditional analysis
#' stats <- get.stats(lm.cells, .design = design)
#' 
#' # Heatmap of cluster-level changes
#' plotTradStats(lm.cells, stats, .split.by = "clustering")
#' 
#' # Cell type-level changes
#' plotTradStats(lm.cells, stats, .split.by = "celltyping")
#' }
#' 
#' @export
#'
plotTradStats <-
  function(
    .lm.obj,
    .stats.obj,
    .split.by = "clustering",
    .coefs = colnames(x = .stats.obj$trad[[.split.by]]$fit$coefficients),
    .q = 0.1,
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.07,
    .label.substr.rm = ""
  ){
    
    sig.shape <- sig <- adj.p <- level <- pop <- term <- coef <- cell.perc <- NULL
    
    .split.by <-
      match.arg(arg = .split.by,
                choices = c("clustering",
                            "celltyping"))
    
    perc.plot <-
      data.frame(x = as.factor(x = 1),
                 pop = colnames(x = .lm.obj$map[[.split.by]]$cell.perc),
                 cell.perc = Matrix::colMeans(x = .lm.obj$map[[.split.by]]$cell.perc)) |>
      (\(x)
       ggplot2::ggplot(data = x,
                       mapping = ggplot2::aes(x = x,
                                              y = pop,
                                              size = cell.perc)) +
         ggplot2::scale_y_discrete(limits = rev) +
         ggplot2::theme_minimal() +
         ggplot2::theme(legend.margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
                        plot.title = ggplot2::element_text(hjust = 0.5),
                        legend.position = "bottom") +
         ggplot2::labs(title = "percentages",
                       x = "",
                       y = "",
                       size = "%") +
         ggplot2::scale_x_discrete(labels = "mean % across\nall samples") +
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
      ggplot2::scale_y_discrete(limits = rev) +
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
                                    midpoint = 0#,
                                    #labels = ~ ifelse(test = .x >= 0,
                                    #                  yes = round(x = 2^.x,
                                    #                              digits = 1),
                                    #                  no = -round(x = 1 / 2^.x,
                                    #                              digits = 1))
      ) +
      ggplot2::scale_size_continuous(labels = ~ formatC(x = 10^(-(.x)),
                                                        format = "g",
                                                        digits = 1),
                                     range = I(x = c(3,6))) +
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
      ggh4x::force_panelsizes(col = grid::unit(x = (as.character(x = dat.df$term) |> unique() |> nchar() |> max()) *
                                                 (unique(x = dat.df$term) |> length()) *
                                                 .col.space.scaler,
                                               units = "in"),
                              rows = grid::unit(x = pmax((unique(x = dat.df$pop) |> length()) * .row.space.scaler,
                                                         3.5),
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
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .x.split Character specifying metadata column for x-axis grouping. Defaults to first 
#'   column (often "Condition").
#' @param .pop Character vector of population names to plot. If NULL, plots all populations from 
#'   \code{.pop.from}.
#' @param .pop.from Character: "clustering" (default) or "celltyping" - which grouping to plot.
#' @param .line.by Character metadata column for connecting paired samples with lines (e.g., 
#'   "Subject" for longitudinal data). Default NULL (no lines).
#' @param .dodge.by Character metadata column for coloring/dodging points. Default NULL (all black).
#' @param .x.space.scaler Numeric scaling factor for x-axis panel width (default 0.25 inches per group).
#' @param .height Numeric plot height in inches (default 1.5).
#' @param .cat.feature.color Character vector of colors for \code{.dodge.by} categories (default 
#'   \code{Color.Palette[1,1:5]}).
#' @param .seed Integer random seed for x-axis jitter (default 123).
#' @param .orientation Character: "wide" (default, all populations in one row) or "square" (facet grid).
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
#' @seealso \code{\link{plotTradStats}}, \code{\link{get.stats}}
#' 
#' @examples
#' \dontrun{
#' # After mapping
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
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
    .lm.obj,
    .x.split = colnames(x = .lm.obj$metadata)[1],
    .pop = NULL,
    .pop.from = "clustering",
    .line.by = NULL,
    .dodge.by = NULL,
    .x.space.scaler = 0.25,
    .height = 1.5,
    .cat.feature.color = Color.Palette[1,1:5],
    .seed = 123,
    .orientation = "wide"
  ){
    
    dodge <- value <- name <- color <- group <- x <- y <- NULL
    
    if(length(x = .x.split) != 1){
      stop(".x.split must be length 1")
    }
    
    if(!(.x.split %in% colnames(x = .lm.obj$metadata))){
      stop(paste0(".x.split must be one of the following: ",
                  paste(x = colnames(x = .lm.obj$metadata),
                        collapse = ", ")))
    }
    
    if(!is.null(x = .dodge.by)){
      if(length(x = .dodge.by) != 1){
        stop(".dodge.by must be length 1")
      }
      
      if(!(.dodge.by %in% colnames(x = .lm.obj$metadata))){
        stop(paste0(".dodge.by must be one of the following: ",
                    paste(x = colnames(x = .lm.obj$metadata),
                          collapse = ", ")))
      }
      
    }
    
    if(!is.null(x = .line.by)){
      if(length(x = .line.by) != 1){
        stop(".line.by must be length 1")
      }
      
      if(!(.line.by %in% colnames(x = .lm.obj$metadata))){
        stop(paste0(".line.by must be one of the following: ",
                    paste(x = colnames(x = .lm.obj$metadata),
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
        sample = rownames(x = .lm.obj$metadata),
        x = .lm.obj$metadata[[.x.split]]
      )
    
    if(!is.null(x = .dodge.by)){
      dat.df$dodge <-
        .lm.obj$metadata[[.dodge.by]]
    }
    
    if(!is.null(x = .line.by)){
      dat.df$group <-
        .lm.obj$metadata[[.line.by]]
    }
    
    if(is.null(x = .pop)){
      
      dat.df <-
        cbind(dat.df,
              .lm.obj$map[[.pop.from]]$cell.perc) |>
        tidyr::pivot_longer(
          cols = colnames(x = .lm.obj$map[[.pop.from]]$cell.perc),
          cols_vary = "slowest"
        ) |>
        as.data.frame()
      
    } else {
      
      dat.df$value <-
        .lm.obj$map[[.pop.from]]$cell.perc[,.pop]
      
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
    
    p
    
  }

#' Plot Abundance
#'
#' Creates dot/line plots showing cell percentages per sample. Similar to \code{plotTradPerc()} 
#' but with slightly different interface and styling options. Visualizes raw abundance data to 
#' inspect distributions and paired/longitudinal patterns.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .x.split Character specifying metadata column for x-axis grouping (default first column).
#' @param .pop Character vector of population names to plot. If NULL, plots all.
#' @param .pop.from Character: "clustering" (default) or "celltyping".
#' @param .subject.id Character metadata column for connecting paired samples with lines (e.g., 
#'   "Subject"). Default NULL.
#' @param .color.by Character metadata column for coloring points. Default NULL (all black).
#' @param .x.space.scaler Numeric x-axis width scaling (default 0.25 inches per group).
#' @param .height Numeric plot height in inches (default 1.5).
#' @param .cat.feature.color Character vector of colors (default \code{Color.Palette[1,1:5]}).
#' @param .seed Integer random seed for jitter (default 123).
#' @param .orientation Character: "wide" (default) or "square" faceting.
#'   
#' @return A \code{ggplot} object showing cell abundances.
#'   
#' @seealso \code{\link{plotTradPerc}}, \code{\link{plotTradStats}}
#' 
#' @examples
#' \dontrun{
#' # Basic abundance plot
#' plotAbundance(lm.cells, .pop = "B.cells")
#' 
#' # With paired subject lines
#' plotAbundance(lm.cells, 
#'               .pop = "B.cells",
#'               .subject.id = "MouseID",
#'               .color.by = "Treatment")
#' }
#' 
#' @export
#'
plotAbundance <-
  function(
    .lm.obj,
    .x.split = colnames(x = .lm.obj$metadata)[1],
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
    
    if(!(.x.split %in% colnames(x = .lm.obj$metadata))){
      stop(paste0(".x.split must be one of the following: ",
                  paste(x = colnames(x = .lm.obj$metadata),
                        collapse = ", ")))
    }
    
    if(!is.null(x = .color.by)){
      if(length(x = .color.by) != 1){
        stop(".color.by must be length 1")
      }
      
      if(!(.color.by %in% colnames(x = .lm.obj$metadata))){
        stop(paste0(".color.by must be one of the following: ",
                    paste(x = colnames(x = .lm.obj$metadata),
                          collapse = ", ")))
      }
      
    }
    
    .orientation <-
      match.arg(arg = .orientation,
                choices = c("wide",
                            "square"))
    
    dat.df <-
      data.frame(
        sample = rownames(x = .lm.obj$metadata),
        x = .lm.obj$metadata[[.x.split]],
        color = if(is.null(x = .color.by)) "black" else .lm.obj$metadata[[.color.by]]
      ) |> 
      (\(x)
       `rownames<-`(x = x,
                    value = rownames(x = .lm.obj$metadata))
      )()
    
    if(is.null(x = .pop)){
      
      dat.df <-
        cbind(dat.df,
              Matrix::t(x = log2(x = .lm.obj$map$fdens + 0.5))) |>
        tidyr::pivot_longer(
          cols = rownames(x = .lm.obj$map$fdens),
          cols_vary = "slowest"
        ) |>
        dplyr::mutate(name = as.character(x = .lm.obj$graph[[.pop.from]]$ids) |>
                        stats::setNames(nm = rownames(x = .lm.obj$map$fdens)) |>
                        (\(x)
                         x[name]
                        )()) |>
        as.data.frame()
      
    } else {
      
      dat.df <-
        cbind(dat.df,
              Matrix::t(x = log2(x = .lm.obj$map$fdens[.lm.obj$graph[[.pop.from]]$ids == .pop,] + 0.5))) |>
        tidyr::pivot_longer(
          cols = rownames(.lm.obj$map$fdens),
          cols_vary = "slowest"
        ) |>
        as.data.frame()
      
    }
    
    if(!is.null(x = .subject.id)){
      dat.df$group <-
        .lm.obj$metadata[[.subject.id]][match(x = dat.df$sample,
                                              table = rownames(x = .lm.obj$metadata))]  
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
                    y = "abundance\nlog2(+0.5)",
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

#' Plot Differential Expression Analysis Results
#'
#' Visualizes differential expression results as a heatmap showing log fold changes for genes/markers 
#' across coefficients. Rows ordered by cluster/celltype expression patterns, with significance 
#' indicated by asterisks. Helps identify which features drive population-level changes.
#'
#' @param .lm.obj A tinydenseR object processed through \code{get.map()}.
#' @param .dea.obj Differential expression results from \code{get.dea()}.
#' @param .coefs Character vector of coefficient names to plot.
#' @param .order.by Character: "clustering" (default) or "celltyping" - order rows by mean expression 
#'   in these groups.
#' @param .markers Character vector of feature names (genes/proteins) to plot. Defaults to features 
#'   shown in cluster/celltype heatmap (top PC contributors for RNA, all markers for cytometry).
#' @param .q Numeric q-value threshold for significance stars (default 0.1).
#' @param .row.space.scaler Numeric row height scaling (default 0.2 inches per feature).
#' @param .col.space.scaler Numeric column width scaling (default 0.5 inches per coefficient).
#' @param .label.substr.rm Character substring to remove from labels (default "").
#'   
#' @return A \code{ggplot} heatmap with log2 fold changes colored and significance marked.
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
#' @seealso \code{\link{get.dea}}, \code{\link{plotBeeswarm}}
#' 
#' @examples
#' \dontrun{
#' # After DEA
#' design <- model.matrix(~ Condition, data = .meta)
#' dea <- get.dea(lm.cells, .design = design)
#' 
#' # Heatmap of DE genes
#' plotDEA(lm.cells, dea, .coefs = "ConditionB", .order.by = "clustering")
#' 
#' # Focus on specific markers
#' plotDEA(lm.cells, dea, 
#'         .coefs = "ConditionB",
#'         .markers = c("CD4", "CD8A", "CD3D"))
#' }
#' 
#' @export
#'
plotDEA <-
  function(
    .lm.obj,
    .dea.obj,
    .coefs,
    .order.by = "clustering",
    .markers = colnames(x = .lm.obj$graph[[.order.by]]$mean.exprs),
    .q = 0.1,
    .row.space.scaler = 0.2,
    .col.space.scaler = 0.065,
    .label.substr.rm = ""
  ){
    
    # R CMD check appeasement
    sig.shape <- sig <- adj.p <- level <- marker <- term <- coef <- cell.perc <- NULL
    
    # Extract log fold changes for selected markers and coefficients
    coef.df <-
      cbind(.dea.obj$coefficients[.markers,.coefs,drop = FALSE]) |>
      dplyr::as_tibble(rownames = "marker") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "coef")
    
    # Extract adjusted p-values
    adj.p.df <-
      cbind(.dea.obj$adj.p[.markers,.coefs,drop = FALSE])  |>
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
                        levels = rev(x = .lm.obj$graph[[.order.by]]$pheatmap$tree_col$labels[
                          .lm.obj$graph[[.order.by]]$pheatmap$tree_col$order]) |>
                          (\(x)
                           x[x %in% marker]
                          )()),
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
                    fill = if(.lm.obj$assay.type == "RNA"){"log2(+0.5)FC"} else {"estimated\ndifference"},
                    size = if((is.na(x = dat.df$sig.shape) |> all())) {
                      paste0("adj.p\n(",
                             "all > ",
                             .q,
                             ")")} else "adj.p") +
      ggh4x::force_panelsizes(col = grid::unit(x = (as.character(x = dat.df$term) |> unique() |> nchar() |> max()) *
                                                 (unique(x = dat.df$term) |> length()) *
                                                 .col.space.scaler,
                                               units = "in"),
                              rows = grid::unit(x = pmax((unique(x = dat.df$marker) |> length()) * .row.space.scaler,
                                                         3.5),
                                                units = "in"))
    
    return(other.plot)
    
  }

#' Scatter Plot with Feature Coloring
#'
#' Creates a 2D scatter plot with flexible x/y features and optional coloring by a third feature. 
#' Useful for exploring relationships between any features in the tinydenseR object (scaled expression, 
#' PCA coordinates, graph embeddings, cluster IDs, metadata, etc.).
#'
#' @param .x.feature Numeric vector for x-axis values (e.g., \code{.lm.obj$scaled.lm[,"CD3"]} or 
#'   \code{.lm.obj$pca$embed[,"PC1"]}).
#' @param .y.feature Numeric vector for y-axis values.
#' @param .color.feature Optional vector for point colors. Can be numeric (continuous coloring) or 
#'   categorical (discrete colors). If \code{NULL}, all points colored identically.
#' @param .x.label Character x-axis label (default "").
#' @param .y.label Character y-axis label (default "").
#' @param .color.label Character color legend label (default "").
#' @param .cat.feature.color Color palette for categorical features 
#'   (default \code{Color.Palette[1,1:5]}). Automatically interpolated if more categories exist.
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
#'   \item Plotting Laplacian Eigenmap coordinates from \code{.lm.obj$graph$LE$embed}
#'   \item Exploring relationships between markers (e.g., CD3 vs CD4)
#'   \item Overlaying metadata on any 2D embedding
#'   \item Creating custom QC plots (e.g., library size vs mitochondrial %)
#' }
#' 
#' Points are plotted in randomized order to avoid bias when groups overlap. For continuous color 
#' features, uses diverging blue-white-red scale. For categorical, generates discrete colors.
#' 
#' @seealso \code{\link{plotPCA}}, \code{\link{plotUMAP}}
#' 
#' @examples
#' \dontrun{
#' # After processing
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
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
#' scatterPlot(.x.feature = lm.cells$scaled.lm[,"CD4"],
#'             .y.feature = lm.cells$scaled.lm[,"CD8A"],
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
#' @param .lm.obj A tinydenseR object processed through \code{get.graph()}.
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
#'   \item \strong{Cytometry}: All markers from \code{.lm.obj$marker}
#'   \item \strong{RNA-seq}: Top 3 positive and 3 negative genes per PC (from highly variable genes)
#' }
#' 
#' Rows (features) are ordered by hierarchical clustering to group co-expressed markers/genes. 
#' Columns (populations) also clustered to reveal relationships between cell types.
#' 
#' The heatmap is generated once during \code{get.graph()} and cached in 
#' \code{.lm.obj$graph[[.id.from]]$pheatmap}. This function simply displays the cached version.
#' 
#' @seealso \code{\link{get.graph}}, \code{\link{celltyping}}
#' 
#' @examples
#' \dontrun{
#' # After clustering
#' lm.cells <- setup.lm.obj(.cells = .cells, .meta = .meta) |>
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
    .lm.obj,
    .id.from = "clustering"
  ){
    
    # Verify graph component exists
    if(is.null(x = .lm.obj$graph[[.id.from]])){
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
    return(gridExtra::grid.arrange(.lm.obj$graph[[.id.from]]$pheatmap$gtable))
    
  }

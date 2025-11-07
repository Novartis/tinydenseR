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
#' This function plots the summary of the landmark feature discovery.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.lm.features.
#' @param .n.feat The number of top features to plot. Defaults to 10.
#' @param .plot.dims The dimensions to plot. Can be set to "umap", "pca", or a combination of two column names from .lm.obj$lm. Defaults to "umap".
#' @param .PC.x The principal component to plot on the x-axis. Defaults to 1 (PC1).
#' @param .PC.y The principal component to plot on the y-axis. Defaults to 2 (PC2).
#' @param .panel.size The size of the panel.
#' @param .point.size The size of the points in the plot.
#' @return A ggiraph object if the ggiraph package is available, otherwise a ggplot object.
#' @note Interactive features require the ggiraph package. If ggiraph is not available, a static plot will be returned with a warning.
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
    
    .plot.x <- .plot.y <- topFeatTab <- NULL
    
    if(length(x = .plot.dims) > 2){
      stop(".plot.dims must be with umap, pca, or a combination of two column names from .lm.obj$lm")
    }
    
    .plot.dims <-
      match.arg(arg = .plot.dims,
                choices = c(
                  "umap",
                  "pca",
                  colnames(x = .lm.obj$lm)
                ),
                several.ok = TRUE)
    
    .plot.x <-
      if(.plot.dims == "umap"){
        .lm.obj$graph$uwot$embedding[,1]
      } else if(.plot.dims == "pca"){
        .lm.obj$pca$embed[,.PC.x]
      } else {
        .lm.obj$lm[,.plot.dims[1]]
      }
    
    .plot.y <-
      if(.plot.dims == "umap"){
        .lm.obj$graph$uwot$embedding[,2]
      } else if(.plot.dims == "pca"){
        .lm.obj$pca$embed[,.PC.y]
      } else {
        .lm.obj$lm[,.plot.dims[2]]
      }
    
    dat.df <-
      cbind(.plot.x,
            .plot.y) |>
      as.data.frame()
    
    colnames(x = dat.df) <-
      if(.plot.dims == "umap"){
        paste0("umap.",1:2)
      } else if(.plot.dims == "pca"){
        paste0("PC",c(.PC.x,.PC.y))
      } else {
        .plot.dims
      }
    
    dat.df$topFeatTab <-
      .lm.obj$interact.plot$lm.features$html
    
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
    
    # Check if ggiraph is available for interactive features
    if (.is_ggiraph_available()) {
      p <- p + ggiraph::geom_point_interactive(size = I(x = .point.size))
      
      ggiraph::girafe(ggobj = p) |>
        ggiraph::girafe_options(
          ggiraph::opts_zoom(max = 10,
                             min = 1)
        )
    } else {
      # Fallback to static plot
      warning("ggiraph package not available. Returning static plot instead of interactive plot.", 
              call. = FALSE)
      p + ggplot2::geom_point(size = I(x = .point.size))
    }
  }

#' Plot PCA
#'
#' This function plots the dimensionality reduction of the landmarks (PCA).
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks and get.graph.
#' @param .PC.x The principal component to plot on the x-axis. Defaults to 1 (PC1).
#' @param .PC.y The principal component to plot on the y-axis. Defaults to 2 (PC2).
#' @param .feature The .feature to plot. It can be a numeric vector or a factor of length equal to the number of landmarks (nrow(.lm.obj$lm)). If numeric, the plot will be colored in relation to the median value of the .feature, or by .midpoint. If a factor, the plot will be colored by the levels of the factor.
#' @param .cat.feature.color The color palette to use in case .feature is categorical. Defaults to  colors (row 1, cols 1:5).
#' @param .panel.size The size of the panel.
#' @param .midpoint The midpoint for the color scale. Only used if .feature is numeric. Defaults to the median of the .feature.
#' @param .plot.title The title of the plot.
#' @param .color.label The label of the color scale.
#' @param .legend.position The position of the legend. Defaults to "right" and can be set to "left", "top", "bottom" or "none".
#' @param .point.size The size of the points in the plot.
#' @param .seed The .seed for the random number generator.
#' @param .hover.stats The statistics to show in the hover. Defaults to "none" and can be set to "marker" to show top features associated with each neighborhood compared to all others.
#' @param .hover.n.features The number of top features to show in the hover in each fold change direction. Defaults to 10 (i.e., 20 total features, 10 in each direction).
#' @note Interactive hover features (when .hover.stats != "none") require the ggiraph package. If ggiraph is not available, a static plot will be returned with a warning.
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
    
    topFeatTab <- feature <- NULL
    
    if(any(!c(.PC.x, .PC.y) %in%
           1:ncol(x = .lm.obj$pca$embed))){
      
      stop(paste0(".PC.x and .PC.y must be one of the following: ",
                  paste(x = 1:ncol(x = .lm.obj$pca$embed),
                        collapse = ", ")))
    }
    
    .hover.stats <-
      match.arg(arg = .hover.stats,
                choices = c(
                  "none",
                  "marker"
                ))
    
    .PC.x <-
      colnames(x = .lm.obj$pca$embed)[.PC.x]
    
    .PC.y <-
      colnames(x = .lm.obj$pca$embed)[.PC.y]
    
    if(is.null(x = .midpoint) &
       is.numeric(x = .feature)){
      .midpoint <- stats::median(x = .feature)
    }
    
    set.seed(seed = .seed)
    dat.df <-
      as.data.frame(x = .lm.obj$pca$embed[,c(.PC.x,.PC.y)]) |>
      cbind(feature = .feature,
            point.size = .point.size)
    
    if(.hover.stats == "marker"){
      
      dat.df$topFeatTab <-
        .lm.obj$interact.plot$lm.features$html
      
    }
    
    dat.df <-
      dat.df  |>
      (\(x)
       # randomize order to avoid plotting one class in front of the other
       # see: https://stackoverflow.com/a/29325361
       x[nrow(x = x) |>
           seq_len() |>
           sample(),]
      )()
    
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
    
    if(is.numeric(x = dat.df$feature)){
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
          )(length(x = unique(x = dat.df$feature))))
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

#' Plot UMAP
#'
#' This function plots the dimensionality reduction of the landmarks (UMAP).
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks and get.graph.
#' @param .feature The .feature to plot. It can be a numeric vector or a factor of length equal to the number of landmarks (nrow(.lm.obj$lm)). If numeric, the plot will be colored in relation to the median value of the .feature, or by .midpoint. If a factor, the plot will be colored by the levels of the factor.
#' @param .cat.feature.color The color palette to use in case .feature is categorical. Defaults to  colors (row 1, cols 1:5).
#' @param .panel.size The size of the panel.
#' @param .midpoint The midpoint for the color scale. Only used if .feature is numeric. Defaults to the median of the .feature.
#' @param .plot.title The title of the plot.
#' @param .color.label The label of the color scale.
#' @param .legend.position The position of the legend. Defaults to "right" and can be set to "left", "top", "bottom" or "none".
#' @param .point.size The size of the points in the plot.
#' @param .seed The .seed for the random number generator.
#' @param .hover.stats The statistics to show in the hover. Defaults to "none" and can be set to "marker" to show top features associated with each neighborhood compared to all others.
#' @param .hover.n.features The number of top features to show in the hover in each fold change direction. Defaults to 10 (i.e., 20 total features, 10 in each direction).
#' @note Interactive hover features (when .hover.stats != "none") require the ggiraph package. If ggiraph is not available, a static plot will be returned with a warning.
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
    
    umap.1 <- umap.2 <- topFeatTab <- feature <- NULL
    
    if(is.null(x = .lm.obj$graph)){
      stop("First run get.graph")
    }
    
    .hover.stats <-
      match.arg(arg = .hover.stats,
                choices = c(
                  "none",
                  "marker"
                ))
    
    if(is.null(x = .midpoint) &
       is.numeric(x = .feature)){
      .midpoint <- stats::median(x = .feature)
    }
    
    set.seed(seed = .seed)
    dat.df <-
      as.data.frame(x = .lm.obj$graph$uwot$embedding) |>
      cbind(feature = .feature,
            point.size = .point.size)
    
    if(.hover.stats == "marker"){
      
      dat.df$topFeatTab <-
        .lm.obj$interact.plot$lm.features$html
      
    }
    
    set.seed(seed = 123)
    dat.df <-
      dat.df  |>
      (\(x)
       # randomize order to avoid plotting one class in front of the other
       # see: https://stackoverflow.com/a/29325361
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
#' This function plots the statistics of modeling outcomes as a function of the sample wise sum of cell-landmark neighboring probabilities
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.map.
#' @param .stats.obj  The list object that is the result of running get.stats.
#' @param .coefs The coefficient level to plot.
#' @param .q The q-value threshold to use for coloring the plot. Defaults to 0.1.
#' @param .q.from The source of the q values. Can be set to "density.weighted.bh.fdr", or "pca.weighted.q". Defaults to "pca.weighted.q".
#' @param .split.by The variable to split the plot by. Can be set to "none", "clustering", or "celltyping". If either of the two latter, .coefs must be provided and, if not, defaults to showing the first coefficient in formula. Defaults to "none".
#' @param .swarm.title The title of the swarm plot. If NULL and if there are no facets (facets are only used if length of `.coefs` is greater than 1 and pops are split by clustering or celltyping), the title will be set to the value of .coefs.
#' @param .label.substr.rm The substring to remove from the labels.
#' @param .point.size The size of the points in the plot. Defaults to 0.1.
#' @param .facet.scales The scales to use. Defaults to "fixed".
#' @param .row.space.scaler The scaler for the row space. Defaults to 0.2.
#' @param .col.space.scaler The scaler for the column space. Defaults to 0.1.
#' @param .panel.width The width of the panel. Defaults to 1.5 and is used only if the length of `.coefs` is greater than 1.
#' @param .legend.position The position of the legend. Defaults to "bottom" and can be set to "left", "top", "bottom", "right" or "none".
#' @param .perc.plot If TRUE, the percentages of cells in each cluster or celltype are plotted. Defaults to TRUE. If FALSE, only the beeswarm plot is shown.
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
#' This function plots the PCA of samples in abundance estimate space.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.map.
#' @param .labels.from The source of the labels. Must be a .lm.obj$metadata column name. Defaults to the first column name of .lm.obj$metadata.
#' @param .cat.feature.color The color palette for the categories. Defaults to the Color.Palette.
#' @param .point.size The size of the points in the plot.
#' @param .panel.size The size of the panel.
#' @param .midpoint The midpoint for the color gradient. Defaults to the median of the `.labels.from` column.
#' @param .seed The .seed for the random number generator.
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
#' This function plots the statistics of modeling outcomes as a function of the abundance of the cells in each cluster or celltype.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.map.
#' @param .stats.obj  The list object that is the result of running get.stats.
#' @param .coefs The coefficient level to plot.
#' @param .split.by The variable to split the plot by. Can be set to "clustering" or "celltyping". Defaults to "clustering".
#' @param .q The q-value threshold to use for coloring the plot.
#' @param .row.space.scaler The scaler for the row space. Defaults to 0.2.
#' @param .col.space.scaler The scaler for the column space. Defaults to 0.5.
#' @param .label.substr.rm The substring to remove from the labels.
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
#' This function plots the actual percentages from traditional analysis to facilitate visual inspection of value distributions.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.map.
#' @param .x.split The column of metadata to use in the x axis. Defaults to the first column of .lm.obj$metadata.
#' @param .pop The name of the pop to plot in the y axis.
#' @param .pop.from Either `"clustering"` or `"celltyping"`. Defaults to `"clustering"`.
#' @param .line.by If provided, draws a line connecting dots from the same subject.
#' @param .dodge.by The column of metadata to use for coloring the points. Defaults to NULL, which makes all points `"black"`.
#' @param .x.space.scaler The scaler for the x axis Defaults to `0.5`.
#' @param .height The height of the plot. Defaults to `1.5`.
#' @param .cat.feature.color The color palette for the categories. Defaults to the Color.Palette.
#' @param .seed The seed for the random number generator that control x-axis jitter. Defaults to `123`.
#' @param .orientation The orientation of the plot. Can be set to `"wide"` or `"square"`. Defaults to `"wide"`.
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
#' This function plots the actual percentages from traditional analysis to facilitate visual inspection of value distributions.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks, get.graph and get.map.
#' @param .x.split The column of metadata to use in the x axis. Defaults to the first column of .lm.obj$metadata.
#' @param .pop The name of the pop to plot in the y axis.
#' @param .pop.from Either `"clustering"` or `"celltyping"`. Defaults to `"clustering"`.
#' @param .subject.id If provided, draws a line connecting dots from the same subject.
#' @param .color.by The column of metadata to use for coloring the points. Defaults to NULL, which makes all points `"black"`.
#' @param .x.space.scaler The scaler for the x axis Defaults to `0.5`.
#' @param .height The height of the plot. Defaults to `1.5`.
#' @param .cat.feature.color The color palette for the categories. Defaults to the Color.Palette.
#' @param .seed The seed for the random number generator that control x-axis jitter. Defaults to `123`.
#' @param .orientation The orientation of the plot. Can be set to `"wide"` or `"square"`. Defaults to `"wide"`.
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
#' This function plots the statistics of modeling outcomes as a function of feature (e.g., genes or proteins) expression.
#'
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with `get.landmarks`, `get.graph` and `get.map`.
#' @param .dea.obj  The list object that is the result of running `get.dea`.
#' @param .coefs The coefficient level to plot.
#' @param .order.by The variable to order the row (proteins) by. Can be set to `"clustering"` or `"celltyping"`. Defaults to `"clustering"`.
#' @param .markers A character vector with the markers to plot. Defaults to the column names of the `"clustering"` or `"celltyping"` mean expression matrix.
#' @param .q The q-value threshold to use for coloring the plot.
#' @param .row.space.scaler The scaler for the row space. Defaults to `0.2`.
#' @param .col.space.scaler The scaler for the column space. Defaults to `0.5`.
#' @param .label.substr.rm The substring to remove from the labels.
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
    
    sig.shape <- sig <- adj.p <- level <- marker <- term <- coef <- cell.perc <- NULL
    
    coef.df <-
      cbind(.dea.obj$coefficients[.markers,.coefs,drop = FALSE]) |>
      dplyr::as_tibble(rownames = "marker") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "coef")
    
    adj.p.df <-
      cbind(.dea.obj$adj.p[.markers,.coefs,drop = FALSE])  |>
      dplyr::as_tibble(rownames = "marker") |>
      tidyr::pivot_longer(cols = dplyr::where(fn = is.numeric),
                          names_to = "term",
                          values_to = "adj.p")
    
    dat.df <-
      dplyr::full_join(x = coef.df,
                       y = adj.p.df,
                       by = c("marker","term")) |>
      dplyr::mutate(marker =
                      factor(
                        x = marker,
                        levels = rev(x = .lm.obj$graph[[.order.by]]$pheatmap$tree_col$labels[
                          .lm.obj$graph[[.order.by]]$pheatmap$tree_col$order]) |>
                          (\(x)
                           x[x %in% marker]
                          )()),
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
                                             y = marker,
                                             fill = coef,
                                             stroke = sig.shape
                      )) +
      ggplot2::guides(stroke = ggplot2::guide_legend(override.aes = list(size = I(x = 5)),
                                                     order = 1),
                      size = ggplot2::guide_legend(override.aes = list(fill = unname(obj = Color.Palette[1,6]),
                                                                       stroke = NA),
                                                   order = 2),
                      fill = ggplot2::guide_colorbar(order = 3)) +
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

#' Scatter plot
#'
#' This function plots a scatter plot of two features colored by a third feature.
#'
#' @param .x.feature The x-axis feature. Defaults to `.lm.obj$scaled.lm[,"CD3"]`.
#' @param .y.feature The y-axis feature. Defaults to `.lm.obj$pca$embed[,"PC1"]`.
#' @param .color.feature The feature to color the plot by. Defaults to `.lm.obj$graph$clustering$ids`.
#' @param .x.label The x-axis label. Defaults to "x".
#' @param .y.label The y-axis label. Defaults to "y".
#' @param .color.label The color label. Defaults to "cluster".
#' @param .cat.feature.color The color palette to use in case .feature is categorical. Defaults to `Color.Palette[1,1:5]`.
#' @param .panel.size The size of the panel.
#' @param .midpoint The midpoint of the color gradient. Defaults to the median of .color.feature. Only used if .color.feature is numeric.
#' @param .plot.title The title of the plot.
#' @param .legend.position The position of the legend. Defaults to "right" and can be set to "left", "top", "bottom" or "none".
#' @param .point.size The size of the points.
#' @param .seed The seed for the random number generator.
#' @export
#' @return A ggplot2 object.
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
    
    if(is.null(x = .midpoint) &
       is.numeric(x = .color.feature)){
      .midpoint <- stats::median(x = .color.feature)
    }
    
    dat.df <-
      data.frame(.x.feature = .x.feature,
                 .y.feature = .y.feature)
    
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
    
    if(is.numeric(x = .color.feature)){
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
          )(length(x = unique(x = .color.feature))))
    }
    
    p +
      ggh4x::force_panelsizes(rows = grid::unit(x = .panel.size,
                                                units = "in")[2],
                              col = grid::unit(x = .panel.size,
                                               units = "in")[1])
  }

#' Plot Heatmap
#' 
#' This function plots a heatmap of mean expression values across clusters or cell types. For cytometry experiments, all markers in `.lm.obj$marker` are returned,  while the top 3 positive and negative genes contributing to each principal component are returned for scRNA-seq experiments.
#' 
#' @param .lm.obj A list object initialized with setup.lm.obj and processed with get.landmarks and get.graph.
#' @param .id.from Either `"clustering"` or `"celltyping"`. Defaults to `"clustering"`.
#' @return Returns a `gtable`, `gTree`, `grob`, `gDesc` object created with `pheatmap` and draws on the current device.
#' @export
plotHeatmap <-
  function(
    .lm.obj,
    .id.from = "clustering"
  ){
    
    if(is.null(x = .lm.obj$graph[[.id.from]])){
      stop(paste0("Please run get.graph first."))
    }
    
    .id.from <- 
      match.arg(
        arg = .id.from,
        choices = c("clustering",
                    "celltyping")
      )
    
    return(gridExtra::grid.arrange(.lm.obj$graph[[.id.from]]$pheatmap$gtable))
    
  }

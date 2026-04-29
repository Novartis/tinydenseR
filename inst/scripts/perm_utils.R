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

# Requires: dplyr, ggplot2, patchwork, tinydenseR

# Function: Count discoveries for a given design
count_discoveries <-
  function(lm_obj, 
           design_mat, 
           coef_name,
           alpha = 0.10) {
    
    lm_obj <-
      tinydenseR::get.lm(
        x = lm_obj, 
        .design = design_mat, 
        .force.recalc = TRUE,
        .verbose = FALSE)
    
    q_vals <- 
      if(inherits(x = lm_obj,
      what = "SingleCellExperiment")) {
        tinydenseR::GetTDR(x = lm_obj)@results$lm[["default"]]$fit$pca.weighted.q[, coef_name]
      } else {
        lm_obj@results$lm[["default"]]$fit$pca.weighted.q[, coef_name]
      }
      
    
    dplyr::tibble(
      m_tests = length(x = q_vals),
      n_sig   = sum(q_vals <= alpha, 
                    na.rm = TRUE),
      min_q   = min(q_vals, 
                    na.rm = TRUE)
    )
  }

# Function: Build complement labels (swap levels)
make_complement_labels <- 
  function(labels) {
    
    lev <-
      levels(x = labels)
    
    ifelse(test = labels == lev[1],
           yes = lev[2], 
           no = lev[1]) |>
      factor(levels = lev)
    
  }

# Function: Encode labels to {0,1} for MCC calculation
encode01 <- 
  function(x, lev = levels(x)) {
    as.integer(x == lev[2])
  }

# Function: Calculate Matthews Correlation Coefficient
mcc_binary <- 
  function(y_true01, y_pred01) {
    tp <- sum(y_true01 == 1 & y_pred01 == 1)
    tn <- sum(y_true01 == 0 & y_pred01 == 0)
    fp <- sum(y_true01 == 0 & y_pred01 == 1)
    fn <- sum(y_true01 == 1 & y_pred01 == 0)
    denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if (denom == 0) return(NA_real_)  # undefined when all samples in one class
    (tp * tn - fp * fn) / denom
  }

# Function: Run exact stratified permutation test
run_stratified_permutation_test <- 
  function(lm_obj, 
           meta_df, 
           formula = ~ Treatment + Batch,
           coef_name, 
           alpha = 0.10) {
    
    # Observed design and stats
    obs_design <-
      model.matrix(object = formula, 
                   data = meta_df)
    
    colnames(x = obs_design) <- 
      gsub(pattern = "^Treatment|^Batch",
           replacement =  "",
           x =  colnames(x = obs_design),
           fixed = FALSE)
    
    obs <-
      count_discoveries(
        lm_obj = lm_obj, 
        design_mat = obs_design, 
        coef_name = coef_name, 
        alpha = alpha)
    
    obs$type <-
      "observed"
    
    obs$run <- 
      "0"
    
    obs$mcc <-
      1.0  # identical labels by definition
    
    # Prepare labels
    obs_labels <-
      factor(x = meta_df$Treatment) |>
      droplevels()
    
    comp_labels <-
      make_complement_labels(labels = obs_labels)
    
    lev <-
      levels(x = obs_labels)
    
    # Encode observed labels for MCC calculation
    obs01 <-
      encode01(obs_labels)
    
    # Generate all stratified permutations
    batch_levels <-
      unique(x = meta_df$Batch)
    
    batch_combos <-
      list()
    
    batch_indices <-
      list()
    
    for(bl in batch_levels) {
      
      idx <-
        which(x = meta_df$Batch == bl)
      
      batch_indices[[bl]] <-
        idx  # Store absolute indices
      
      n_depl <-
        sum(meta_df$Treatment[idx] == coef_name)
      
      batch_combos[[bl]] <-
        length(x = idx) |>
        combn(m = n_depl, 
              simplify = FALSE)
    }
    
    # Cartesian product of batch combinations
    all_combos <- 
      lapply(X = batch_combos, 
             FUN = seq_along) |>
      expand.grid(stringsAsFactors = FALSE)
    
    names(x = all_combos) <- 
      paste0("batch_",
             seq_along(along.with = batch_levels))
    
    # Remove original combination
    orig_idx <-
      which(x = meta_df$Treatment == coef_name)
    
    all_combos <-
      all_combos[
        !(apply(X = all_combos, 
                MARGIN = 1,
                FUN = function(row) {
                  
                  chosen <- 
                    Map(f = function(i, combos, idx) idx[combos[[i]]],
                        row,
                        batch_combos, 
                        batch_indices) |>
                    unlist()
                  
                  setequal(x = chosen, 
                           y = orig_idx)
                  
                })), , drop = FALSE]
    
    # Evaluate each permutation
    # Iterate row-wise (each row provides one combination)
    perm_results_list <-
      nrow(x = all_combos) |>
      seq_len() |>
      lapply(FUN = function(i) {
        
        # Extract the i-th combination as a list (preserve order)
        combo_idx <-
          as.list(x = all_combos[i, , drop = FALSE])
        
        # chosen: map relative positions to absolute indices
        chosen <- 
          Map(function(i, combos, idx) idx[combos[[i]]],
              combo_idx, 
              batch_combos, 
              batch_indices) |>
          unlist()
        
        # Build permuted metadata
        perm_meta <-
          meta_df
        
        perm_meta$Treatment <- 
          rep(x = lev[1],
              times = nrow(meta_df)) |>
          factor(levels = lev)
        
        perm_meta$Treatment[chosen] <- 
          lev[2]
        
        # Detect complement
        perm_type <-
          if(identical(x = perm_meta$Treatment, 
                       y = comp_labels)){ 
            "complement"
          } else {
            "permuted"
          }
        
        # Calculate MCC
        perm01 <-
          encode01(perm_meta$Treatment, lev = levels(obs_labels))
        
        perm_mcc <-
          mcc_binary(obs01, perm01)
        
        # Design matrix with cleaned column names
        perm_design <-
          model.matrix(object = formula, 
                       data = perm_meta)
        
        colnames(x = perm_design) <- 
          gsub(pattern = "^Treatment|^Batch",
               replacement = "", 
               x = colnames(x = perm_design),
               fixed = FALSE)
        
        # Compute result and annotate
        res <- 
          count_discoveries(
            lm_obj = lm_obj, 
            design_mat = perm_design, 
            coef_name = coef_name, 
            alpha = alpha)
        
        res$type <-
          perm_type
        
        res$run <-
          unlist(x = combo_idx) |>
          paste(collapse = "-")
        
        res$mcc <-
          perm_mcc
        
        res
      })
    
    # Bind rows at the end
    perm_results <-
      dplyr::bind_rows(perm_results_list)
    
    dplyr::bind_rows(obs, 
                     perm_results)
  }

# Function: Plot permutation test results - n_sig
plot_n_sig <-
  function(results) {
    
    dplyr::transmute(.data = results,
                     condition,
                     m_tests,
                     mcc,
                     n_sig) |>
      (\(x) {
        
        # Split by condition
        conditions <- unique(x$condition)
        
        # Create a plot for each condition
        plot_list <- lapply(conditions, function(cond) {
          
          x_cond <- x[x$condition == cond, ]
          
          # Top scatter plot (n_sig on x-axis, mcc on y-axis)
          p_top <- 
            ggplot2::ggplot(data = x_cond,
                            mapping = ggplot2::aes(x = n_sig,
                                                   y = mcc)) +
            ggplot2::theme_bw() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                           axis.title.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank()) +
            ggplot2::geom_point(size = 1) +
            ggplot2::labs(title = cond,
                          y = "MCC")
          
          # Extract the actual x-axis limits from scatter
          p_top_built <- ggplot2::ggplot_build(p_top)
          x_range <- p_top_built$layout$panel_params[[1]]$x.range
          
          # Manually calculate histogram bins
          breaks <- seq(x_range[1], x_range[2], length.out = 51)
          hist_data <- hist(x_cond$n_sig, breaks = breaks, plot = FALSE)
          
          # Create data frame for geom_rect
          hist_df <- data.frame(
            xmin = hist_data$breaks[-length(hist_data$breaks)],
            xmax = hist_data$breaks[-1],
            count = hist_data$counts
          )
          
          # Bottom histogram using geom_rect (y-axis will scale freely)
          p_bottom <- 
            ggplot2::ggplot(data = hist_df,
                            mapping = ggplot2::aes(xmin = xmin, 
                                                   xmax = xmax,
                                                   ymin = 0, 
                                                   ymax = count)) +
            ggplot2::geom_rect(fill = "grey95", color = "grey30") +
            ggplot2::scale_x_continuous(limits = x_range, 
                                        expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
            ggplot2::theme_bw() +
            ggplot2::labs(x = paste0("# q < 0.1"),
                          y = paste0("Counts"))
          
          # Combine for this condition (top scatter + bottom histogram)
          p_top / p_bottom + 
            patchwork::plot_layout(heights = c(1, 2))
        })
        
        # Combine all conditions
        patchwork::wrap_plots(plot_list, ncol = length(conditions)) +
          patchwork::plot_annotation(
            theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
          )
        
      })()
    
  }
# Function: Compute pePC t-test p-value for a given design
count_discoveries_pePC <-
  function(lm_obj,
           design_mat,
           coef_name,
           red_design_mat,
           term_of_interest,
           alpha = 0.10) {
    
    # Full model
    lm_obj <-
      tinydenseR::get.lm(
        x = lm_obj,
        .design = design_mat,
        .force.recalc = TRUE,
        .verbose = FALSE)
    
    q_vals <-
      if(inherits(x = lm_obj,
      what = "SingleCellExperiment")) {
        tinydenseR::GetTDR(x = lm_obj)@results$lm[["default"]]$fit$pca.weighted.q[, coef_name]
      } else {
        lm_obj@results$lm[["default"]]$fit$pca.weighted.q[, coef_name]
      }
      
    
    # Reduced model (intercept-only or otherwise)
    lm_obj <-
      tinydenseR::get.lm(
        x = lm_obj,
        .design = red_design_mat,
        .model.name = "noCondition",
        .force.recalc = TRUE,
        .verbose = FALSE)
    
    # Compute pePC embedding
    lm_obj <-
      tinydenseR::get.embedding(
        x = lm_obj,
        .full.model = "default",
        .red.model = "noCondition",
        .term.of.interest = term_of_interest,
        .verbose = FALSE)
    
    # t-test on pePC1 vs the grouping variable
    pePC1 <-
      lm_obj$sample.embed$pepc[[term_of_interest]]$coord[, 1]
    
    cond_labels <-
      if(inherits(x = lm_obj,
      what = "SingleCellExperiment")) {
        tinydenseR::GetTDR(x = lm_obj)@metadata[[term_of_interest]]
      } else {
        lm_obj@metadata[[term_of_interest]]
      }
      
    
    pePC_p <-
      tryCatch(
        t.test(formula = pePC1 ~ cond_labels)$p.value,
        error = function(e) NA_real_)
    
    dplyr::tibble(
      m_tests   = length(x = q_vals),
      n_sig     = sum(q_vals <= alpha, na.rm = TRUE),
      min_q     = min(q_vals, na.rm = TRUE),
      pePC1_p   = pePC_p
    )
  }

# Function: Run exact stratified permutation test with pePC embedding
run_stratified_permutation_test_pePC <-
  function(lm_obj,
           meta_df,
           formula = ~ Treatment + Batch,
           coef_name,
           term_of_interest,
           alpha = 0.10) {
    
    # Observed design and stats
    obs_design <-
      model.matrix(object = formula,
                   data = meta_df)
    
    colnames(x = obs_design) <-
      gsub(pattern = "^Treatment|^Batch",
           replacement =  "",
           x =  colnames(x = obs_design),
           fixed = FALSE)
    
    # Reduced model design (intercept-only)
    red_design <-
      model.matrix(object = ~ 1,
                   data = meta_df)
    
    obs <-
      count_discoveries_pePC(
        lm_obj = lm_obj,
        design_mat = obs_design,
        coef_name = coef_name,
        red_design_mat = red_design,
        term_of_interest = term_of_interest,
        alpha = alpha)
    
    obs$type <-
      "observed"
    
    obs$run <-
      "0"
    
    obs$mcc <-
      1.0
    
    # Prepare labels
    obs_labels <-
      factor(x = meta_df$Treatment) |>
      droplevels()
    
    comp_labels <-
      make_complement_labels(labels = obs_labels)
    
    lev <-
      levels(x = obs_labels)
    
    obs01 <-
      encode01(obs_labels)
    
    # Generate all stratified permutations
    batch_levels <-
      unique(x = meta_df$Batch)
    
    batch_combos <-
      list()
    
    batch_indices <-
      list()
    
    for(bl in batch_levels) {
      
      idx <-
        which(x = meta_df$Batch == bl)
      
      batch_indices[[bl]] <-
        idx
      
      n_depl <-
        sum(meta_df$Treatment[idx] == coef_name)
      
      batch_combos[[bl]] <-
        length(x = idx) |>
        combn(m = n_depl,
              simplify = FALSE)
    }
    
    all_combos <-
      lapply(X = batch_combos,
             FUN = seq_along) |>
      expand.grid(stringsAsFactors = FALSE)
    
    names(x = all_combos) <-
      paste0("batch_",
             seq_along(along.with = batch_levels))
    
    orig_idx <-
      which(x = meta_df$Treatment == coef_name)
    
    all_combos <-
      all_combos[
        !(apply(X = all_combos,
                MARGIN = 1,
                FUN = function(row) {
                  
                  chosen <-
                    Map(f = function(i, combos, idx) idx[combos[[i]]],
                        row,
                        batch_combos,
                        batch_indices) |>
                    unlist()
                  
                  setequal(x = chosen,
                           y = orig_idx)
                  
                })), , drop = FALSE]
    
    perm_results_list <-
      nrow(x = all_combos) |>
      seq_len() |>
      lapply(FUN = function(i) {
        
        combo_idx <-
          as.list(x = all_combos[i, , drop = FALSE])
        
        chosen <-
          Map(function(i, combos, idx) idx[combos[[i]]],
              combo_idx,
              batch_combos,
              batch_indices) |>
          unlist()
        
        perm_meta <-
          meta_df
        
        perm_meta$Treatment <-
          rep(x = lev[1],
              times = nrow(meta_df)) |>
          factor(levels = lev)
        
        perm_meta$Treatment[chosen] <-
          lev[2]
        
        perm_type <-
          if(identical(x = perm_meta$Treatment,
                       y = comp_labels)){
            "complement"
          } else {
            "permuted"
          }
        
        perm01 <-
          encode01(perm_meta$Treatment, lev = levels(obs_labels))
        
        perm_mcc <-
          mcc_binary(obs01, perm01)
        
        perm_design <-
          model.matrix(object = formula,
                       data = perm_meta)
        
        colnames(x = perm_design) <-
          gsub(pattern = "^Treatment|^Batch",
               replacement = "",
               x = colnames(x = perm_design),
               fixed = FALSE)
        
        # Must also update metadata on the object so that
        # get.embedding can find the term of interest labels
        perm_lm_obj <- lm_obj
        if(inherits(
          x = perm_lm_obj,
          what = "SingleCellExperiment")) {
          
          tmp.tdr <-
            tinydenseR::GetTDR(x = perm_lm_obj)
          
          tmp.tdr@metadata[[term_of_interest]] <-
            as.character(x = perm_meta$Treatment)
          
          perm_lm_obj <-
            tinydenseR::SetTDR(x = perm_lm_obj, tdr = tmp.tdr)

          rm(tmp.tdr)
          
        } else {
            perm_lm_obj@metadata[[term_of_interest]] <-
              as.character(x = perm_meta$Treatment)
        }
        
        res <-
          count_discoveries_pePC(
            lm_obj = perm_lm_obj,
            design_mat = perm_design,
            coef_name = coef_name,
            red_design_mat = red_design,
            term_of_interest = term_of_interest,
            alpha = alpha)
        
        res$type <-
          perm_type
        
        res$run <-
          unlist(x = combo_idx) |>
          paste(collapse = "-")
        
        res$mcc <-
          perm_mcc
        
        res
      })
    
    perm_results <-
      dplyr::bind_rows(perm_results_list)
    
    dplyr::bind_rows(obs,
                     perm_results)
  }


# Function: Plot permutation test results - min_q
plot_min_q <-
  function(results){
    
    dplyr::transmute(.data = results,
                     condition,
                     m_tests,
                     mcc,
                     min_q) |>
      dplyr::mutate(min_q = -log10(x = min_q)) |>
      (\(x) {
        
        # Split by condition
        conditions <- unique(x$condition)
        
        # Create a plot for each condition
        plot_list <- lapply(conditions, function(cond) {
          
          x_cond <- x[x$condition == cond, ]
          
          # Top scatter plot (min_q on x-axis, mcc on y-axis)
          p_top <- 
            ggplot2::ggplot(data = x_cond,
                            mapping = ggplot2::aes(x = min_q,
                                                   y = mcc)) +
            ggplot2::theme_bw() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                           axis.title.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank()) +
            ggplot2::geom_point(size = 1) +
            ggplot2::labs(title = cond,
                          y = "MCC")
          
          # Extract the actual x-axis limits from scatter
          p_top_built <- ggplot2::ggplot_build(p_top)
          x_range <- p_top_built$layout$panel_params[[1]]$x.range
          
          # Manually calculate histogram bins
          breaks <- seq(x_range[1], x_range[2], length.out = 51)
          hist_data <- hist(x_cond$min_q, breaks = breaks, plot = FALSE)
          
          # Create data frame for geom_rect
          hist_df <- data.frame(
            xmin = hist_data$breaks[-length(hist_data$breaks)],
            xmax = hist_data$breaks[-1],
            count = hist_data$counts
          )
          
          # Bottom histogram using geom_rect (y-axis will scale freely)
          p_bottom <- 
            ggplot2::ggplot(data = hist_df,
                            mapping = ggplot2::aes(xmin = xmin, 
                                                   xmax = xmax,
                                                   ymin = 0, 
                                                   ymax = count)) +
            ggplot2::geom_vline(xintercept = 1,
                                color = "red",
                                linetype = "dashed") +
            ggplot2::geom_rect(fill = "grey95", color = "grey30") +
            ggplot2::scale_x_continuous(limits = x_range, 
                                        expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
            ggplot2::theme_bw() +
            ggplot2::labs(x = paste0("min q\n(neg. log10)"),
                          y = paste0("Counts"))
          
          # Combine for this condition (top scatter + bottom histogram)
          p_top / p_bottom + 
            patchwork::plot_layout(heights = c(1, 2))
        })
        
        # Combine all conditions
        patchwork::wrap_plots(plot_list, ncol = length(conditions)) +
          patchwork::plot_annotation(
            theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
          )
        
      })()
    
    
  }


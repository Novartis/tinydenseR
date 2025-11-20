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

# Helper functions for testing

#' Create a minimal .lm.obj for testing
#' @param n_cells Number of cells per sample
#' @param n_markers Number of markers
#' @param n_samples Number of samples
create_test_lm_obj <- function(n_cells = 10, n_markers = 3, n_samples = 2) {
  .cells <- lapply(X = 1:n_samples, FUN = function(i) {
    mat <- matrix(
      data = runif(n = n_cells * n_markers),
      nrow = n_cells,
      ncol = n_markers,
      dimnames = list(
        paste0("sample", i, "_cell_", 1:n_cells),
        paste0("marker_", 1:n_markers)
      )
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(object = mat, file = uri, compress = FALSE)
    return(uri)
  })
  names(.cells) <- paste0("sample", 1:n_samples)
  
  .meta <- data.frame(
    row.names = names(.cells),
    group = sample(x = LETTERS[1:2], size = n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  return(list(cells = .cells, meta = .meta))
}

#' Create a mock .lm.obj with graph data for plotting tests
create_mock_graph_obj <- function(n_points = 20, n_clusters = 3) {
  list(
    graph = list(
      uwot = list(embedding = matrix(runif(n_points * 2), ncol = 2)),
      clustering = list(ids = factor(sample(1:n_clusters, n_points, replace = TRUE)))
    ),
    pca = list(embed = matrix(runif(n_points * 2), ncol = 2, dimnames = list(NULL, c("PC1", "PC2")))),
    lm = matrix(runif(n_points * 2), ncol = 2)
  )
}

#' Clean up temporary files created by create_test_lm_obj
cleanup_test_files <- function(test_obj) {
  if ("cells" %in% names(test_obj)) {
    lapply(X = test_obj$cells, FUN = unlink)
  }
}

#' Create test design matrix
create_test_design <- function(n_samples, include_intercept = TRUE) {
  if (include_intercept) {
    cbind(
      intercept = 1,
      group = sample(c(0, 1), n_samples, replace = TRUE)
    )
  } else {
    matrix(sample(c(0, 1), n_samples, replace = TRUE), ncol = 1)
  }
}

#' Run full pipeline for stability testing
#' @param test_data Output from create_test_lm_obj
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages
#' @param nPC Number of principal components (default NULL for auto)
#' @param k Number of neighbors for kNN graph (default NULL for auto)
run_full_pipeline <- function(test_data, seed = 123, verbose = FALSE, nPC = NULL, k = NULL) {
  # Determine appropriate nPC if not specified
  if (is.null(nPC)) {
    first_mat <- readRDS(test_data$cells[[1]])
    n_cells <- nrow(first_mat) * length(test_data$cells)
    n_markers <- ncol(first_mat)
    # Use min of: 30 (default), smallest sample size, or n_markers
    nPC <- min(30, n_cells, n_markers)
  }
  
  # Determine appropriate k if not specified
  if (is.null(k)) {
    first_mat <- readRDS(test_data$cells[[1]])
    n_landmarks <- ceiling(nrow(first_mat) * length(test_data$cells) * 0.15)  # default prop
    # k must be < n_landmarks; use min(5, n_landmarks-1)
    k <- min(5, max(2, n_landmarks - 1))
  }
  
  result <- setup.lm.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = names(test_data$cells[[1]] |> readRDS() |> colnames()),
    .assay.type = "cyto",
    .verbose = verbose
  ) |>
    get.landmarks(.verbose = verbose, .seed = seed, .nPC = nPC) |>
    get.graph(.k = k, .verbose = verbose, .seed = seed) |>
    get.map(.verbose = verbose, .seed = seed)
  
  return(result)
}

#' Compute Jaccard similarity between two character vectors
#' @param set1 First set
#' @param set2 Second set
#' @return Jaccard index (0 to 1)
jaccard_similarity <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

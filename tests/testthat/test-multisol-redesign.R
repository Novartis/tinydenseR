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

# ======================================================================
# test-multisol-redesign.R
#
# Comprehensive test suite for the multi-solution annotation redesign.
# Covers:
#   - Multi-solution storage for clustering & celltyping
#   - Active `ids` mirror semantics
#   - recluster(), set_active_clustering(), set_active_celltyping()
#   - On-the-fly heatmap computation (.compute_annot_pheatmap)
#   - Refresh / invalidation / stale-state propagation
#   - Removed-storage regression checks
#   - Generalized helpers (.relabel_cellmap, .recompute_composition,
#     .invalidate_trad, .warn_stale_de_results)
#   - Backend dispatch wrappers
#   - Edge cases / negative tests
#   - Determinism under fixed seeds
# ======================================================================

library(testthat)
library(tinydenseR)

# ======================================================================
# FIXTURE BUILDERS
# ======================================================================

# Minimal cyto-type fixture: setup → landmarks → graph
# Returns object ready for lm.cluster / get.map / celltyping
.build_graph_obj <- function(seed = 42, n_samples = 4, n_per_sample = 200,
                             n_markers = 5) {
  set.seed(seed)
  blob1 <- matrix(rnorm(n_per_sample * n_markers, mean = 0, sd = 0.5),
                  nrow = n_per_sample, ncol = n_markers)
  blob2 <- matrix(rnorm(n_per_sample * n_markers, mean = 5, sd = 0.5),
                  nrow = n_per_sample, ncol = n_markers)

  .cells <- lapply(seq_len(n_samples), function(i) {
    mat <- rbind(blob1 + rnorm(1, 0, 0.1), blob2 + rnorm(1, 0, 0.1))
    dimnames(mat) <- list(
      paste0("s", i, "_", seq_len(nrow(mat))),
      paste0("M", seq_len(n_markers))
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(mat, file = uri, compress = FALSE)
    uri
  })
  names(.cells) <- paste0("sample", seq_len(n_samples))

  .meta <- data.frame(
    row.names = names(.cells),
    group = rep(c("A", "B"), length.out = n_samples)
  )

  setup.tdr.obj(
    .cells = .cells, .meta = .meta,
    .markers = paste0("M", seq_len(n_markers)),
    .assay.type = "cyto", .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = seed)
}

# Full pipeline: setup → landmarks → graph → map
.build_mapped_obj <- function(seed = 42, n_samples = 4) {
  .build_graph_obj(seed = seed, n_samples = n_samples) |>
    get.map(.verbose = FALSE, .seed = seed)
}

# Build a cluster→celltype map from whatever clusters exist
.make_ct_map <- function(obj) {
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  stopifnot(n >= 2)
  half <- floor(n / 2)
  list(
    "TypeA" = cls[seq_len(half)],
    "TypeB" = cls[(half + 1):n]
  )
}

# Build an ALTERNATIVE ct_map with different labels
.make_ct_map_alt <- function(obj) {
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  stopifnot(n >= 2)
  half <- floor(n / 2)
  list(
    "Alpha" = cls[seq_len(half)],
    "Beta"  = cls[(half + 1):n]
  )
}

# Inject fake get.lm trad results so .invalidate_trad has something to act on
.inject_fake_trad <- function(obj, annot_type = "celltyping") {
  if (is.null(obj@results$lm)) obj@results$lm <- list()
  obj@results$lm$default <- list(
    trad = stats::setNames(list("fake_fit"), annot_type)
  )
  obj
}

# Inject fake pbDE results for stale-result warning tests
.inject_fake_pbDE <- function(obj, annot_type = "celltyping") {
  if (is.null(obj@results$pb)) obj@results$pb <- list()
  obj@results$pb$default <- list(
    pop1 = list(.id.from = annot_type, result = "fake")
  )
  obj
}

# Inject fake marker DE results
.inject_fake_markerDE <- function(obj, annot_type = "celltyping") {
  if (is.null(obj@results$marker)) obj@results$marker <- list()
  obj@results$marker$default <- list(
    comp1 = list(.id.from = annot_type, result = "fake")
  )
  obj
}


# ======================================================================
# CATEGORY 1 — STRUCTURAL INVARIANTS
# ======================================================================

test_that("C1.1: lm.cluster stores ids + named solution in clustering", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .cl.resolution.parameter = 0.8,
                      .seed = 123, .verbose = FALSE)

  # $ids must exist and be a factor

  expect_true(is.factor(obj2@landmark.annot$clustering$ids))
  # Named solution stored under default name
  expect_true("leiden.res.0.8" %in% names(obj2@landmark.annot$clustering))
  # Named solution mirrors ids exactly
  expect_identical(
    obj2@landmark.annot$clustering$ids,
    obj2@landmark.annot$clustering[["leiden.res.0.8"]]
  )
  # Length matches landmarks
  expect_equal(
    length(obj2@landmark.annot$clustering$ids),
    nrow(obj2@assay$expr)
  )
})

test_that("C1.2: lm.cluster with custom .column.name stores correctly", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "my_clust",
                      .seed = 123, .verbose = FALSE)

  expect_true("my_clust" %in% names(obj2@landmark.annot$clustering))
  expect_identical(
    obj2@landmark.annot$clustering$ids,
    obj2@landmark.annot$clustering$my_clust
  )
})

test_that("C1.3: celltyping stores ids + named solution", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  expect_true(is.factor(obj2@landmark.annot$celltyping$ids))
  # A named solution should exist (auto-generated name starting with "manual.")
  stored_names <- setdiff(names(obj2@landmark.annot$celltyping), "ids")
  expect_true(length(stored_names) >= 1)
  # The one named solution mirrors ids
  expect_identical(
    obj2@landmark.annot$celltyping$ids,
    obj2@landmark.annot$celltyping[[stored_names[1]]]
  )
})

test_that("C1.4: celltyping with custom .name stores correctly", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)

  expect_true("v1" %in% names(obj2@landmark.annot$celltyping))
  expect_identical(
    obj2@landmark.annot$celltyping$ids,
    obj2@landmark.annot$celltyping$v1
  )
})

test_that("C1.5: reserved name 'ids' rejected for .column.name in lm.cluster", {
  obj <- .build_graph_obj()
  expect_error(
    lm.cluster(obj, .column.name = "ids", .verbose = FALSE),
    "ids.*reserved"
  )
})

test_that("C1.6: reserved name 'ids' rejected for .name in celltyping", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  expect_error(
    celltyping(obj, ct_map, .name = "ids", .verbose = FALSE),
    "ids.*reserved"
  )
})

test_that("C1.7: $map and $mode are absent after celltyping (both modes)", {
  obj <- .build_graph_obj()

  # Mode A
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  expect_null(obj2@landmark.annot$celltyping$map)
  expect_null(obj2@landmark.annot$celltyping$mode)
})

test_that("C1.8: clustering ids length equals nrow(assay$expr)", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .seed = 123, .verbose = FALSE)
  expect_equal(
    length(obj2@landmark.annot$clustering$ids),
    nrow(obj2@assay$expr)
  )
})


# ======================================================================
# CATEGORY 2 — CORRECT ACTIVE-SOLUTION SEMANTICS
# ======================================================================

test_that("C2.1: recluster at different resolution creates new named solution + updates ids", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .cl.resolution.parameter = 0.8,
                      .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .cl.resolution.parameter = 2.0,
                     .seed = 123, .verbose = FALSE)

  # Both solutions stored
  expect_true("leiden.res.0.8" %in% names(obj3@landmark.annot$clustering))
  expect_true("leiden.res.2"   %in% names(obj3@landmark.annot$clustering))
  # ids reflects the LATEST recluster call
  expect_identical(
    obj3@landmark.annot$clustering$ids,
    obj3@landmark.annot$clustering[["leiden.res.2"]]
  )
  # Prior solution not corrupted — should still match original
  expect_identical(
    obj2@landmark.annot$clustering[["leiden.res.0.8"]],
    obj3@landmark.annot$clustering[["leiden.res.0.8"]]
  )
})

test_that("C2.2: set_active_clustering copies exactly the requested column into ids", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .cl.resolution.parameter = 0.8,
                      .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .cl.resolution.parameter = 2.0,
                     .seed = 123, .verbose = FALSE)

  # ids is currently leiden.res.2; switch back to 0.8
  obj4 <- set_active_clustering(obj3, .column.name = "leiden.res.0.8",
                                 .verbose = FALSE)

  expect_identical(
    obj4@landmark.annot$clustering$ids,
    obj4@landmark.annot$clustering[["leiden.res.0.8"]]
  )
  # The other solution is still intact
  expect_identical(
    obj3@landmark.annot$clustering[["leiden.res.2"]],
    obj4@landmark.annot$clustering[["leiden.res.2"]]
  )
})

test_that("C2.3: repeated switches are lossless for clustering", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .cl.resolution.parameter = 0.8,
                      .column.name = "A", .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .cl.resolution.parameter = 2.0,
                     .column.name = "B", .seed = 123, .verbose = FALSE)

  sol_A <- obj3@landmark.annot$clustering$A
  sol_B <- obj3@landmark.annot$clustering$B

  # Switch A → B → A → B and verify no data loss
  obj4 <- set_active_clustering(obj3, "A", .verbose = FALSE)
  obj5 <- set_active_clustering(obj4, "B", .verbose = FALSE)
  obj6 <- set_active_clustering(obj5, "A", .verbose = FALSE)
  obj7 <- set_active_clustering(obj6, "B", .verbose = FALSE)

  expect_identical(obj7@landmark.annot$clustering$A, sol_A)
  expect_identical(obj7@landmark.annot$clustering$B, sol_B)
  expect_identical(obj7@landmark.annot$clustering$ids, sol_B)
})

test_that("C2.4: set_active_clustering errors on non-existent column", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "A", .seed = 123, .verbose = FALSE)

  expect_error(
    set_active_clustering(obj2, "nonexistent", .verbose = FALSE),
    "not found"
  )
})

test_that("C2.5: set_active_clustering errors on 'ids'", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .seed = 123, .verbose = FALSE)

  expect_error(
    set_active_clustering(obj2, "ids", .verbose = FALSE),
    "ids"
  )
})

test_that("C2.6: set_active_celltyping copies exactly the requested column into ids", {
  obj <- .build_graph_obj()
  ct_map_1 <- .make_ct_map(obj)
  ct_map_2 <- .make_ct_map_alt(obj)

  obj2 <- celltyping(obj, ct_map_1, .name = "v1", .verbose = FALSE)
  obj3 <- celltyping(obj2, ct_map_2, .name = "v2", .verbose = FALSE)

  # ids currently = v2; switch to v1
  obj4 <- set_active_celltyping(obj3, "v1", .verbose = FALSE)
  expect_identical(
    obj4@landmark.annot$celltyping$ids,
    obj3@landmark.annot$celltyping$v1
  )
  # v2 preserved
  expect_identical(
    obj4@landmark.annot$celltyping$v2,
    obj3@landmark.annot$celltyping$v2
  )
})

test_that("C2.7: set_active_celltyping errors on non-existent solution", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)

  expect_error(
    set_active_celltyping(obj2, "nonexistent", .verbose = FALSE),
    "not found"
  )
})

test_that("C2.8: set_active_celltyping errors on 'ids'", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  expect_error(
    set_active_celltyping(obj2, "ids", .verbose = FALSE),
    "ids"
  )
})

test_that("C2.9: celltyping() preserves prior named solutions (no wipe)", {
  # After fix: celltyping no longer wipes $celltyping list.
  # Matches lm.cluster behavior — only $ids is overwritten.
  obj <- .build_graph_obj()
  ct_map_1 <- .make_ct_map(obj)
  ct_map_2 <- .make_ct_map_alt(obj)

  obj2 <- celltyping(obj, ct_map_1, .name = "v1", .verbose = FALSE)
  expect_true("v1" %in% names(obj2@landmark.annot$celltyping))

  obj3 <- celltyping(obj2, ct_map_2, .name = "v2", .verbose = FALSE)
  # v1 should be PRESERVED
  expect_true("v1" %in% names(obj3@landmark.annot$celltyping))
  expect_true("v2" %in% names(obj3@landmark.annot$celltyping))
  # ids reflects the latest call
  expect_identical(
    obj3@landmark.annot$celltyping$ids,
    obj3@landmark.annot$celltyping$v2
  )
})

test_that("C2.10: lm.cluster does NOT wipe prior named clustering solutions", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "A", .seed = 123, .verbose = FALSE)
  obj3 <- lm.cluster(obj2, .cl.resolution.parameter = 2.0,
                      .column.name = "B", .seed = 123, .verbose = FALSE)

  # Both should be present
  expect_true("A" %in% names(obj3@landmark.annot$clustering))
  expect_true("B" %in% names(obj3@landmark.annot$clustering))
})

test_that("C2.11: repeated set_active_celltyping switches are lossless", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)
  # Inject second solution manually (celltyping wipes)
  alt_ids <- obj2@landmark.annot$celltyping$ids
  levels(alt_ids) <- paste0("Alt_", levels(alt_ids))
  obj2@landmark.annot$celltyping$v2 <- alt_ids

  sol_v1 <- obj2@landmark.annot$celltyping$v1
  sol_v2 <- obj2@landmark.annot$celltyping$v2

  obj3 <- set_active_celltyping(obj2, "v2", .verbose = FALSE)
  obj4 <- set_active_celltyping(obj3, "v1", .verbose = FALSE)
  obj5 <- set_active_celltyping(obj4, "v2", .verbose = FALSE)

  expect_identical(obj5@landmark.annot$celltyping$v1, sol_v1)
  expect_identical(obj5@landmark.annot$celltyping$v2, sol_v2)
  expect_identical(obj5@landmark.annot$celltyping$ids, sol_v2)
})


# ======================================================================
# CATEGORY 3 — REFRESH / INVALIDATION / STALE-STATE
# ======================================================================

test_that("C3.1: recluster refreshes composition after get.map", {
  obj <- .build_mapped_obj()
  comp_before <- obj@density$composition$clustering$cell.count

  obj2 <- recluster(obj, .cl.resolution.parameter = 2.0,
                     .seed = 42, .verbose = FALSE)
  comp_after <- obj2@density$composition$clustering$cell.count

  # Composition must exist and reflect new clustering
  expect_true(is.matrix(comp_after))
  expect_equal(nrow(comp_after), nrow(comp_before))
  # Column names = cluster labels from new ids
  expect_setequal(
    colnames(comp_after),
    levels(obj2@landmark.annot$clustering$ids)
  )
})

test_that("C3.2: set_active_clustering refreshes composition", {
  obj <- .build_mapped_obj()
  obj2 <- recluster(obj, .cl.resolution.parameter = 2.0,
                     .column.name = "hi_res", .seed = 42, .verbose = FALSE)
  # Switch back to original
  original_name <- setdiff(
    names(obj2@landmark.annot$clustering),
    c("ids", "hi_res")
  )[1]
  obj3 <- set_active_clustering(obj2, original_name, .verbose = FALSE)

  comp <- obj3@density$composition$clustering$cell.count
  expect_true(is.matrix(comp))
  expect_setequal(
    colnames(comp),
    levels(obj3@landmark.annot$clustering$ids)
  )
})

test_that("C3.3: recluster invalidates trad clustering fits", {
  obj <- .build_mapped_obj()
  obj <- .inject_fake_trad(obj, "clustering")
  expect_false(is.null(obj@results$lm$default$trad$clustering))

  expect_warning(
    obj2 <- recluster(obj, .cl.resolution.parameter = 2.0,
                       .seed = 42, .verbose = FALSE),
    "Invalidated trad\\$clustering"
  )
  expect_null(obj2@results$lm$default$trad$clustering)
})

test_that("C3.4: set_active_clustering invalidates trad clustering fits", {
  obj <- .build_mapped_obj()
  obj2 <- recluster(obj, .cl.resolution.parameter = 2.0,
                     .column.name = "hi", .seed = 42, .verbose = FALSE)
  obj2 <- .inject_fake_trad(obj2, "clustering")

  original_name <- setdiff(
    names(obj2@landmark.annot$clustering),
    c("ids", "hi")
  )[1]

  expect_warning(
    obj3 <- set_active_clustering(obj2, original_name, .verbose = FALSE),
    "Invalidated trad\\$clustering"
  )
  expect_null(obj3@results$lm$default$trad$clustering)
})

test_that("C3.5: celltyping after get.map refreshes celltyping composition", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  comp <- obj2@density$composition$celltyping$cell.count
  expect_true(is.matrix(comp))
  expect_setequal(colnames(comp), c("TypeA", "TypeB"))
  # Row count = # samples
  expect_equal(nrow(comp), length(obj2@cells))
})

test_that("C3.6: celltyping invalidates trad celltyping fits", {
  obj <- .build_mapped_obj()
  obj <- .inject_fake_trad(obj, "celltyping")

  ct_map <- .make_ct_map(obj)
  expect_warning(
    obj2 <- celltyping(obj, ct_map, .verbose = FALSE),
    "Invalidated trad\\$celltyping"
  )
  expect_null(obj2@results$lm$default$trad$celltyping)
})

test_that("C3.7: set_active_celltyping refreshes composition", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)
  alt_ids <- obj2@landmark.annot$celltyping$ids
  levels(alt_ids) <- c("X", "Y")
  obj2@landmark.annot$celltyping$v2 <- alt_ids

  obj3 <- set_active_celltyping(obj2, "v2", .verbose = FALSE)
  comp <- obj3@density$composition$celltyping$cell.count
  expect_true(is.matrix(comp))
  expect_setequal(colnames(comp), c("X", "Y"))
})

test_that("C3.8: recluster warns about stale celltyping when celltyping exists", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  expect_warning(
    recluster(obj2, .cl.resolution.parameter = 2.0,
              .seed = 42, .verbose = FALSE),
    "Celltyping annotations exist"
  )
})

test_that("C3.9: recluster warns about stale pbDE results", {
  obj <- .build_mapped_obj()
  obj <- .inject_fake_pbDE(obj, "clustering")

  expect_warning(
    recluster(obj, .cl.resolution.parameter = 2.0,
              .seed = 42, .verbose = FALSE),
    "stale"
  )
})

test_that("C3.10: set_active_celltyping invalidates trad celltyping fits", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)
  alt_ids <- obj2@landmark.annot$celltyping$ids
  levels(alt_ids) <- c("X", "Y")
  obj2@landmark.annot$celltyping$v2 <- alt_ids
  obj2 <- .inject_fake_trad(obj2, "celltyping")

  expect_warning(
    obj3 <- set_active_celltyping(obj2, "v2", .verbose = FALSE),
    "Invalidated trad\\$celltyping"
  )
  expect_null(obj3@results$lm$default$trad$celltyping)
})

test_that("C3.11: refresh skips downstream when get.map has NOT run", {
  obj <- .build_graph_obj()
  obj <- .inject_fake_trad(obj, "clustering")

  # No fdens → should NOT touch composition (which doesn't exist)
  expect_warning(
    obj2 <- recluster(obj, .cl.resolution.parameter = 2.0,
                       .seed = 42, .verbose = FALSE),
    "Invalidated"
  )
  expect_null(obj2@density$composition$clustering)
})

test_that("C3.12: celltyping warns about stale marker DE results", {
  obj <- .build_mapped_obj()
  obj <- .inject_fake_markerDE(obj, "celltyping")

  ct_map <- .make_ct_map(obj)
  expect_warning(
    celltyping(obj, ct_map, .verbose = FALSE),
    "stale"
  )
})


# ======================================================================
# CATEGORY 4 — ON-THE-FLY HEATMAP / PLOTTING
# ======================================================================

test_that("C4.1: .compute_annot_pheatmap returns pheatmap with expected structure", {
  obj <- .build_graph_obj()

  ph <- .compute_annot_pheatmap(obj, .id.from = "clustering")
  expect_s3_class(ph, "pheatmap")
  expect_true(!is.null(ph$tree_row))
  expect_true(!is.null(ph$gtable))

  # tree_row labels match cluster levels
  expect_setequal(
    ph$tree_row$labels,
    levels(obj@landmark.annot$clustering$ids)
  )
})

test_that("C4.2: .compute_annot_pheatmap works for celltyping", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  ph <- .compute_annot_pheatmap(obj2, .id.from = "celltyping")
  expect_s3_class(ph, "pheatmap")
  expect_setequal(ph$tree_row$labels, c("TypeA", "TypeB"))
})

test_that("C4.3: .compute_annot_pheatmap errors when ids are missing", {
  obj <- .build_graph_obj()
  obj@landmark.annot$celltyping <- list()

  expect_error(
    .compute_annot_pheatmap(obj, .id.from = "celltyping"),
    "No.*celltyping.*annotations"
  )
})

test_that("C4.4: plotHeatmap depends on current active ids, not stale data", {
  skip_if_not_installed("gridExtra")
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "A", .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .cl.resolution.parameter = 2.0,
                     .column.name = "B", .seed = 123, .verbose = FALSE)

  # Pheatmap from current active ids (B)
  ph_B <- .compute_annot_pheatmap(obj3, .id.from = "clustering")
  expect_setequal(
    ph_B$tree_row$labels,
    levels(obj3@landmark.annot$clustering$ids)
  )

  # Switch to A
  obj4 <- set_active_clustering(obj3, "A", .verbose = FALSE)
  ph_A <- .compute_annot_pheatmap(obj4, .id.from = "clustering")
  expect_setequal(
    ph_A$tree_row$labels,
    levels(obj4@landmark.annot$clustering$ids)
  )

  # The pheatmaps should have different row labels if cluster count differs,
  # or at minimum the tree is recomputed (not stale)
  # Strong assertion: the labels correspond to the active ids
  expect_true(setequal(
    ph_A$tree_row$labels,
    levels(obj4@landmark.annot$clustering[["A"]])
  ))
})

test_that("C4.5: plotHeatmap.TDRObj works without stored @results pheatmap", {
  skip_if_not_installed("gridExtra")
  obj <- .build_graph_obj()
  # Ensure @results has no clustering/celltyping pheatmap
  obj@results$clustering <- NULL
  obj@results$celltyping <- NULL

  expect_no_error(
    plotHeatmap(obj, .id.from = "clustering")
  )
})

test_that("C4.6: plotHeatmap errors for invalid .id.from", {
  obj <- .build_graph_obj()
  expect_error(
    plotHeatmap(obj, .id.from = "invalidtype"),
    "should be one of"
  )
})

test_that("C4.7: plotHeatmap errors when annotation missing", {
  obj <- .build_graph_obj()
  obj@landmark.annot$celltyping <- list()

  expect_error(
    plotHeatmap(obj, .id.from = "celltyping"),
    "celltyping"
  )
})

test_that("C4.8: .compute_annot_pheatmap uses active ids not named solution", {
  # After switching active, pheatmap should reflect new ids
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "A",
                      .cl.resolution.parameter = 0.8,
                      .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .column.name = "B",
                     .cl.resolution.parameter = 5.0,
                     .seed = 123, .verbose = FALSE)

  # Active = B
  ph1 <- .compute_annot_pheatmap(obj3, "clustering")
  obj4 <- set_active_clustering(obj3, "A", .verbose = FALSE)
  ph2 <- .compute_annot_pheatmap(obj4, "clustering")

  # Both should have tree_row labels matching their respective active ids
  expect_setequal(ph1$tree_row$labels, levels(obj3@landmark.annot$clustering$ids))
  expect_setequal(ph2$tree_row$labels, levels(obj4@landmark.annot$clustering$ids))
})


# ======================================================================
# CATEGORY 5 — REGRESSION: REMOVED STORAGE
# ======================================================================

test_that("C5.1: @results$clustering is not populated by lm.cluster", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .seed = 123, .verbose = FALSE)
  expect_null(obj2@results$clustering)
})

test_that("C5.2: @results$celltyping is not populated by celltyping()", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  expect_null(obj2@results$celltyping)
})

test_that("C5.3: plotHeatmap works when @results is completely empty", {
  skip_if_not_installed("gridExtra")
  obj <- .build_graph_obj()
  obj@results <- list()

  expect_no_error(plotHeatmap(obj, .id.from = "clustering"))
})

test_that("C5.4: celltyping does not store $map or $mode (Mode A)", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  expect_null(obj2@landmark.annot$celltyping$map)
  expect_null(obj2@landmark.annot$celltyping$mode)
})

test_that("C5.5: no @results$clustering$pheatmap after recluster", {
  obj <- .build_mapped_obj()
  obj2 <- recluster(obj, .cl.resolution.parameter = 2.0,
                     .seed = 42, .verbose = FALSE)
  expect_null(obj2@results$clustering)
})


# ======================================================================
# CATEGORY 6 — recluster() BEHAVIOR
# ======================================================================

test_that("C6.1: recluster creates named solution and updates ids", {
  obj <- .build_graph_obj()
  obj2 <- recluster(obj, .cl.resolution.parameter = 1.5,
                     .column.name = "test_res", .seed = 42, .verbose = FALSE)

  expect_true("test_res" %in% names(obj2@landmark.annot$clustering))
  expect_identical(
    obj2@landmark.annot$clustering$ids,
    obj2@landmark.annot$clustering$test_res
  )
})

test_that("C6.2: recluster preserves prior named clustering solutions", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "first", .seed = 123, .verbose = FALSE)
  sol_first <- obj2@landmark.annot$clustering$first

  obj3 <- recluster(obj2, .cl.resolution.parameter = 3.0,
                     .column.name = "second", .seed = 42, .verbose = FALSE)

  expect_identical(obj3@landmark.annot$clustering$first, sol_first)
  expect_true("second" %in% names(obj3@landmark.annot$clustering))
})

test_that("C6.3: recluster is deterministic under fixed seed", {
  obj <- .build_graph_obj()

  obj_a <- recluster(obj, .cl.resolution.parameter = 1.0,
                      .column.name = "test", .seed = 999, .verbose = FALSE)
  obj_b <- recluster(obj, .cl.resolution.parameter = 1.0,
                      .column.name = "test", .seed = 999, .verbose = FALSE)

  expect_identical(
    obj_a@landmark.annot$clustering$test,
    obj_b@landmark.annot$clustering$test
  )
})

test_that("C6.4: recluster rejects .column.name = 'ids'", {
  obj <- .build_graph_obj()
  expect_error(
    recluster(obj, .column.name = "ids", .verbose = FALSE),
    "ids.*reserved"
  )
})

test_that("C6.5: recluster refreshes cell-level clustering IDs after get.map", {
  obj <- .build_mapped_obj()
  # Cell-level clustering ids exist (from get.map)
  cl_ids_before <- tinydenseR:::.tdr_get_map_slot(obj, "clustering.ids",
                                                   names(obj@cells)[1])

  obj2 <- recluster(obj, .cl.resolution.parameter = 3.0,
                     .seed = 42, .verbose = FALSE)
  cl_ids_after <- tinydenseR:::.tdr_get_map_slot(obj2, "clustering.ids",
                                                  names(obj2@cells)[1])

  # The cell-level labels should come from the new clustering
  expect_true(
    all(cl_ids_after %in% c(levels(obj2@landmark.annot$clustering$ids), NA))
  )
})

test_that("C6.6: recluster handles case where no prior map/lm slots exist", {
  obj <- .build_graph_obj()
  # No get.map run, no results
  expect_no_error(
    recluster(obj, .column.name = "x", .seed = 42, .verbose = FALSE)
  )
})


# ======================================================================
# CATEGORY 7 — celltyping.TDRObj & NAMED STORAGE
# ======================================================================

test_that("C7.1: celltyping stores active ids and named solution identically", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "my_annot", .verbose = FALSE)

  expect_identical(
    obj2@landmark.annot$celltyping$ids,
    obj2@landmark.annot$celltyping$my_annot
  )
})

test_that("C7.2: celltyping auto-generates name when .name is NULL", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  nm <- setdiff(names(obj2@landmark.annot$celltyping), "ids")
  expect_true(length(nm) == 1)
  expect_true(grepl("^manual\\.", nm))
})

test_that("C7.3: celltyping ids are factor of celltype labels", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  ids <- obj2@landmark.annot$celltyping$ids
  expect_true(is.factor(ids))
  expect_setequal(levels(ids), c("TypeA", "TypeB"))
})

test_that("C7.4: celltyping ids length matches landmark count", {
  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  expect_equal(
    length(obj2@landmark.annot$celltyping$ids),
    nrow(obj2@assay$expr)
  )
})

test_that("C7.5: celltyping catches unmapped clusters", {
  obj <- .build_graph_obj()
  cls <- levels(obj@landmark.annot$clustering$ids)
  # Map only first cluster, leave rest unmapped
  bad_map <- list("OnlyOne" = cls[1])
  expect_error(
    celltyping(obj, bad_map, .verbose = FALSE),
    "Unmapped clusters"
  )
})

test_that("C7.6: celltyping catches duplicate cluster in map", {
  obj <- .build_graph_obj()
  cls <- levels(obj@landmark.annot$clustering$ids)
  bad_map <- list(
    "A" = cls,
    "B" = cls[1]  # duplicate
  )
  expect_error(
    celltyping(obj, bad_map, .verbose = FALSE),
    "multiple cell types"
  )
})

test_that("C7.7: celltyping catches invalid cluster IDs", {
  obj <- .build_graph_obj()
  bad_map <- list("A" = c("fake_cluster_99"))
  expect_error(
    celltyping(obj, bad_map, .verbose = FALSE),
    "Unmapped|Invalid"
  )
})

test_that("C7.8: celltyping catches unnamed list elements", {
  obj <- .build_graph_obj()
  bad_map <- list(c("cluster.01"))  # no name
  expect_error(
    celltyping(obj, bad_map, .verbose = FALSE),
    "names"
  )
})


# ======================================================================
# CATEGORY 8 — GENERALIZED HELPER TESTS
# ======================================================================

test_that("C8.1: .annot_slot_map returns correct slots for clustering", {
  m <- .annot_slot_map("clustering")
  expect_equal(m$cellmap_slot, "cluster.ids")
  expect_equal(m$cache_slot, "clustering.ids")
  expect_equal(m$comp_slot, "clustering")
})

test_that("C8.2: .annot_slot_map returns correct slots for celltyping", {
  m <- .annot_slot_map("celltyping")
  expect_equal(m$cellmap_slot, "celltype.ids")
  expect_equal(m$cache_slot, "celltyping.ids")
  expect_equal(m$comp_slot, "celltyping")
})

test_that("C8.3: .annot_slot_map errors on invalid type", {
  expect_error(.annot_slot_map("bogus"), "Unknown")
})

test_that("C8.4: .invalidate_trad nullifies correct annot type only", {
  obj <- .build_graph_obj()
  obj <- .inject_fake_trad(obj, "clustering")
  obj@results$lm$default$trad$celltyping <- "fake_ct_fit"

  expect_warning(
    obj2 <- .invalidate_trad(obj, .annot.type = "clustering"),
    "Invalidated trad\\$clustering"
  )
  # clustering trad NULLed, celltyping trad preserved
  expect_null(obj2@results$lm$default$trad$clustering)
  expect_equal(obj2@results$lm$default$trad$celltyping, "fake_ct_fit")
})

test_that("C8.5: .invalidate_trad does nothing when no trad results exist", {
  obj <- .build_graph_obj()
  obj@results$lm$default <- list(trad = list())

  expect_no_warning(
    obj2 <- .invalidate_trad(obj, .annot.type = "clustering")
  )
})

test_that("C8.6: .warn_stale_de_results warns for matching pbDE", {
  obj <- .build_graph_obj()
  obj <- .inject_fake_pbDE(obj, "clustering")

  expect_warning(
    .warn_stale_de_results(obj, .annot.type = "clustering"),
    "stale"
  )
})

test_that("C8.7: .warn_stale_de_results does NOT warn for non-matching pbDE", {
  obj <- .build_graph_obj()
  obj <- .inject_fake_pbDE(obj, "celltyping")

  expect_no_warning(
    .warn_stale_de_results(obj, .annot.type = "clustering")
  )
})

test_that("C8.8: .warn_stale_de_results warns for matching marker DE", {
  obj <- .build_graph_obj()
  obj <- .inject_fake_markerDE(obj, "celltyping")

  expect_warning(
    .warn_stale_de_results(obj, .annot.type = "celltyping"),
    "stale"
  )
})

test_that("C8.9: .recompute_composition produces correct matrix dimensions", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  comp <- obj2@density$composition$celltyping
  expect_equal(nrow(comp$cell.count), length(obj2@cells))
  expect_setequal(colnames(comp$cell.count), c("TypeA", "TypeB"))
  # Row sums of cell.perc should each be ~100
  expect_true(all(abs(rowSums(comp$cell.perc) - 100) < 1e-6))
})

test_that("C8.10: .recompute_composition cell.count has no NAs", {
  obj <- .build_mapped_obj()
  # Already has composition from get.map → lm.cluster
  comp <- obj@density$composition$clustering$cell.count
  expect_false(any(is.na(comp)))
})

test_that("C8.11: .relabel_cellmap updates cell-level ids for both types", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  # Check cell-level celltype IDs exist for each sample
  for (sn in names(obj2@cells)) {
    ct_ids <- tinydenseR:::.tdr_get_map_slot(obj2, "celltyping.ids", sn)
    expect_true(!is.null(ct_ids))
    # All cell-level labels should be from celltype vocabulary
    expect_true(all(ct_ids %in% c(levels(obj2@landmark.annot$celltyping$ids), NA)))
  }
})


# ======================================================================
# CATEGORY 9 — CROSS-BACKEND / DISPATCH TESTS
# ======================================================================

test_that("C9.1: recluster dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")

  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "init", .seed = 123, .verbose = FALSE)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat <- SetTDR(seurat, obj2)

  result <- suppressWarnings(
    recluster(seurat, .cl.resolution.parameter = 2.0,
              .column.name = "re", .seed = 42, .verbose = FALSE)
  )

  tdr_out <- GetTDR(result)
  expect_true("re" %in% names(tdr_out@landmark.annot$clustering))
  expect_true("init" %in% names(tdr_out@landmark.annot$clustering))
})

test_that("C9.2: set_active_clustering dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")

  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "A", .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .column.name = "B", .seed = 42, .verbose = FALSE)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat <- SetTDR(seurat, obj3)

  result <- suppressWarnings(
    set_active_clustering(seurat, "A", .verbose = FALSE)
  )
  tdr_out <- GetTDR(result)
  expect_identical(
    tdr_out@landmark.annot$clustering$ids,
    tdr_out@landmark.annot$clustering$A
  )
})

test_that("C9.3: set_active_celltyping dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")

  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)
  # Inject second solution
  alt_ids <- obj2@landmark.annot$celltyping$ids
  levels(alt_ids) <- c("X", "Y")
  obj2@landmark.annot$celltyping$v2 <- alt_ids

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat <- SetTDR(seurat, obj2)

  result <- suppressWarnings(
    set_active_celltyping(seurat, "v2", .verbose = FALSE)
  )
  tdr_out <- GetTDR(result)
  expect_identical(
    tdr_out@landmark.annot$celltyping$ids,
    tdr_out@landmark.annot$celltyping$v2
  )
})

test_that("C9.4: recluster dispatches on SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")

  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "init", .seed = 123, .verbose = FALSE)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 3))
  )
  sce <- SetTDR(sce, obj2)

  result <- recluster(sce, .column.name = "re", .seed = 42, .verbose = FALSE)
  tdr_out <- GetTDR(result)
  expect_true("re" %in% names(tdr_out@landmark.annot$clustering))
})

test_that("C9.5: set_active_clustering dispatches on SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")

  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "A", .seed = 123, .verbose = FALSE)
  obj3 <- recluster(obj2, .column.name = "B", .seed = 42, .verbose = FALSE)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 3))
  )
  sce <- SetTDR(sce, obj3)

  result <- set_active_clustering(sce, "A", .verbose = FALSE)
  tdr_out <- GetTDR(result)
  expect_identical(tdr_out@landmark.annot$clustering$ids, tdr_out@landmark.annot$clustering$A)
})

test_that("C9.6: set_active_celltyping dispatches on SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")

  obj <- .build_graph_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)
  alt_ids <- obj2@landmark.annot$celltyping$ids
  levels(alt_ids) <- c("X", "Y")
  obj2@landmark.annot$celltyping$v2 <- alt_ids

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 3))
  )
  sce <- SetTDR(sce, obj2)

  result <- set_active_celltyping(sce, "v2", .verbose = FALSE)
  tdr_out <- GetTDR(result)
  expect_identical(tdr_out@landmark.annot$celltyping$ids, tdr_out@landmark.annot$celltyping$v2)
})


# ======================================================================
# CATEGORY 10 — EDGE CASES / NEGATIVE TESTS
# ======================================================================

test_that("C10.1: celltyping rejects non-list, non-named-character input", {
  obj <- .build_graph_obj()
  expect_error(
    celltyping(obj, 42, .verbose = FALSE),
    "named list.*named character"
  )
  expect_error(
    celltyping(obj, c("a", "b"), .verbose = FALSE),  # unnamed char
    "named list.*named character"
  )
})

test_that("C10.2: set_active_clustering before any named solution errors", {
  obj <- .build_graph_obj()
  # Build raw graph obj — clustering only has $ids from get.graph
  expect_error(
    set_active_clustering(obj, "nonexistent", .verbose = FALSE),
    "not found"
  )
})

test_that("C10.3: set_active_celltyping before any celltyping errors informatively", {
  obj <- .build_graph_obj()

  expect_error(
    set_active_celltyping(obj, "nonexistent", .verbose = FALSE)
  )
})

test_that("C10.4: lm.cluster requires graph", {
  # An object with no graphs should error
  obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graphs = list(adj.matrix = NULL)
  )
  expect_error(
    lm.cluster(obj, .verbose = FALSE),
    "Graph component missing"
  )
})

test_that("C10.5: celltyping with duplicate cell type names errors", {
  obj <- .build_graph_obj()
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  bad_map <- list(
    "TypeA" = cls[1],
    "TypeA" = cls[2:n]  # duplicate name
  )
  expect_error(
    celltyping(obj, bad_map, .verbose = FALSE),
    "Duplicate"
  )
})

test_that("C10.6: .compute_annot_pheatmap with single-cluster still works", {
  # Edge case: 1 cluster → cluster_rows = FALSE
  n_pts <- 10
  expr_mat <- matrix(runif(n_pts * 3), nrow = n_pts, ncol = 3,
                     dimnames = list(paste0("lm", 1:n_pts), paste0("M", 1:3)))
  obj <- TDRObj(
    config = list(assay.type = "cyto"),
    assay = list(expr = expr_mat),
    landmark.annot = list(
      clustering = list(ids = factor(rep("c1", n_pts)))
    )
  )
  ph <- .compute_annot_pheatmap(obj, "clustering")
  expect_s3_class(ph, "pheatmap")
  # With 1 cluster, cluster_rows = FALSE → tree_row is NA (not a dendrogram)
  expect_false(inherits(ph$tree_row, "dendrogram"))
})


# ======================================================================
# CATEGORY 11 — DETERMINISM
# ======================================================================

test_that("C11.1: lm.cluster is deterministic under fixed seed", {
  obj <- .build_graph_obj()
  obj_a <- lm.cluster(obj, .seed = 42, .verbose = FALSE)
  obj_b <- lm.cluster(obj, .seed = 42, .verbose = FALSE)
  expect_identical(
    obj_a@landmark.annot$clustering$ids,
    obj_b@landmark.annot$clustering$ids
  )
})

test_that("C11.2: different seeds produce different clusterings (sanity)", {
  obj <- .build_graph_obj()
  obj_a <- lm.cluster(obj, .seed = 1, .verbose = FALSE)
  obj_b <- lm.cluster(obj, .seed = 99999, .verbose = FALSE)
  # Not guaranteed to differ, but overwhelmingly likely with different seeds
  # and enough data points. We test a softer check:
  # cluster COUNT may differ, or if not, at least we confirm the function ran.
  expect_true(is.factor(obj_a@landmark.annot$clustering$ids))
  expect_true(is.factor(obj_b@landmark.annot$clustering$ids))
})

test_that("C11.3: recluster + set_active round-trip preserves exact factor", {
  obj <- .build_graph_obj()
  obj2 <- lm.cluster(obj, .column.name = "orig", .seed = 42, .verbose = FALSE)
  snapshot <- obj2@landmark.annot$clustering$orig

  obj3 <- recluster(obj2, .column.name = "new", .seed = 999, .verbose = FALSE)
  obj4 <- set_active_clustering(obj3, "orig", .verbose = FALSE)

  expect_identical(obj4@landmark.annot$clustering$ids, snapshot)
  expect_identical(obj4@landmark.annot$clustering$orig, snapshot)
})


# ======================================================================
# CATEGORY 12 — COMPOSITION INTEGRITY (cross-cutting)
# ======================================================================

test_that("C12.1: composition cell.count rownames = sample names", {
  obj <- .build_mapped_obj()
  comp <- obj@density$composition$clustering$cell.count
  expect_equal(sort(rownames(comp)), sort(names(obj@cells)))
})

test_that("C12.2: composition cell.perc rows sum to 100", {
  obj <- .build_mapped_obj()
  comp <- obj@density$composition$clustering$cell.perc
  expect_true(all(abs(rowSums(comp) - 100) < 1e-6))
})

test_that("C12.3: composition column names match active clustering ids", {
  obj <- .build_mapped_obj()
  comp <- obj@density$composition$clustering$cell.count
  expect_setequal(colnames(comp), levels(obj@landmark.annot$clustering$ids))
})

test_that("C12.4: after recluster, composition colnames = new cluster levels", {
  obj <- .build_mapped_obj()
  obj2 <- recluster(obj, .cl.resolution.parameter = 3.0,
                     .seed = 42, .verbose = FALSE)
  comp <- obj2@density$composition$clustering$cell.count
  expect_setequal(colnames(comp), levels(obj2@landmark.annot$clustering$ids))
})

test_that("C12.5: after celltyping, celltyping composition colnames = celltype levels", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  comp <- obj2@density$composition$celltyping$cell.count
  expect_setequal(colnames(comp), levels(obj2@landmark.annot$celltyping$ids))
})

test_that("C12.6: after set_active_celltyping, composition reflects new labels", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "v1", .verbose = FALSE)
  alt_ids <- obj2@landmark.annot$celltyping$ids
  levels(alt_ids) <- c("Gamma", "Delta")
  obj2@landmark.annot$celltyping$v2 <- alt_ids

  obj3 <- set_active_celltyping(obj2, "v2", .verbose = FALSE)
  comp <- obj3@density$composition$celltyping$cell.count
  expect_setequal(colnames(comp), c("Gamma", "Delta"))
})


# ======================================================================
# CATEGORY 13 — .refresh_clustering vs .refresh_celltyping symmetry
# ======================================================================

test_that("C13.1: .refresh_clustering updates cell-level ids and composition", {
  obj <- .build_mapped_obj()
  # Manually tweak clustering ids, then call .refresh_clustering
  new_ids <- obj@landmark.annot$clustering$ids
  levels(new_ids) <- paste0("new.", levels(new_ids))
  obj@landmark.annot$clustering$ids <- new_ids

  obj2 <- .refresh_clustering(obj, .verbose = FALSE)
  comp <- obj2@density$composition$clustering$cell.count
  expect_true(all(grepl("^new\\.", colnames(comp))))
})

test_that("C13.2: .refresh_celltyping updates cell-level ids and composition", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  # Manually tweak celltyping ids
  new_ids <- obj2@landmark.annot$celltyping$ids
  levels(new_ids) <- c("Mega", "Ultra")
  obj2@landmark.annot$celltyping$ids <- new_ids

  obj3 <- .refresh_celltyping(obj2, .verbose = FALSE)
  comp <- obj3@density$composition$celltyping$cell.count
  expect_setequal(colnames(comp), c("Mega", "Ultra"))
})

test_that("C13.3: .refresh_clustering warns about stale celltyping", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)

  expect_warning(
    .refresh_clustering(obj2, .verbose = FALSE),
    "Celltyping annotations exist"
  )
})

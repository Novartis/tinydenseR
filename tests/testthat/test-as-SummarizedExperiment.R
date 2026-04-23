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

library(testthat)
library(tinydenseR)

# ======================================================================
# Helper: build a minimal TDRObj with get.map() completed
# ======================================================================
create_mapped_tdr <- function(n_lm = 10, n_samples = 4, n_markers = 5,
                              n_clusters = 3, with_fuzzy_graphs = TRUE,
                              with_celltyping = TRUE) {
  set.seed(42)
  lm_names   <- paste0("lm", seq_len(n_lm))
  smpl_names <- paste0("sample", seq_len(n_samples))

  # Expression matrix (landmarks x markers)
  expr_mat <- matrix(runif(n_lm * n_markers), nrow = n_lm, ncol = n_markers,
                     dimnames = list(lm_names, paste0("marker", seq_len(n_markers))))

  # fdens: L x N
  fdens <- matrix(runif(n_lm * n_samples, min = 0.1, max = 5),
                  nrow = n_lm, ncol = n_samples,
                  dimnames = list(lm_names, smpl_names))

  # Y: log2(fdens + 0.5)
  Y <- log2(fdens + 0.5)

  # metadata
  meta <- data.frame(
    group = rep(c("A", "B"), length.out = n_samples),
    row.names = smpl_names,
    stringsAsFactors = FALSE
  )

  # clustering
  cl_ids <- factor(sample(seq_len(n_clusters), n_lm, replace = TRUE))
  names(cl_ids) <- lm_names
  clustering <- list(ids = cl_ids)

  # additional named clustering solution
  cl_alt <- factor(sample(c("x", "y"), n_lm, replace = TRUE))
  names(cl_alt) <- lm_names
  clustering[["alt_cl"]] <- cl_alt

  # celltyping
  celltyping <- list()
  if (with_celltyping) {
    ct_ids <- factor(sample(c("Tcell", "Bcell", "Mono"), n_lm, replace = TRUE))
    names(ct_ids) <- lm_names
    celltyping$ids <- ct_ids
    ct_alt <- factor(sample(c("CD4", "CD8"), n_lm, replace = TRUE))
    names(ct_alt) <- lm_names
    celltyping[["manual"]] <- ct_alt
  }

  # fuzzy graphs (in-memory, cells x landmarks or landmarks x cells)
  fg_list <- NULL
  if (with_fuzzy_graphs) {
    fg_list <- stats::setNames(
      lapply(seq_len(n_samples), function(i) {
        n_cells <- 20L
        # landmarks in rows, cells in columns
        fg <- Matrix::sparseMatrix(
          i = sample(n_lm, n_cells * 3, replace = TRUE),
          j = rep(seq_len(n_cells), each = 3),
          x = runif(n_cells * 3),
          dims = c(n_lm, n_cells),
          dimnames = list(lm_names,
                          paste0(smpl_names[i], "_cell", seq_len(n_cells)))
        )
        fg
      }),
      smpl_names
    )
  }

  # PCA / UMAP embeddings (for plot functions)
  pca_coord <- matrix(runif(n_lm * 2), ncol = 2,
                       dimnames = list(lm_names, c("PC1", "PC2")))
  umap_coord <- matrix(runif(n_lm * 2), ncol = 2,
                        dimnames = list(lm_names, c("UMAP1", "UMAP2")))

  TDRObj(
    cells          = stats::setNames(
      as.list(paste0("/fake/path/", smpl_names, ".rds")), smpl_names),
    metadata       = meta,
    config         = list(assay.type = "cyto", markers = colnames(expr_mat)),
    assay          = list(expr = expr_mat),
    landmark.embed = list(
      pca  = list(coord = pca_coord),
      umap = list(coord = umap_coord),
      le   = list()
    ),
    landmark.annot = list(clustering = clustering, celltyping = celltyping),
    graphs         = list(adj.matrix = NULL, snn = NULL, fgraph = NULL),
    density        = list(fdens = fdens, Y = Y),
    cellmap        = list(fuzzy.graphs = fg_list),
    sample.embed   = list(),
    results        = list()
  )
}


# ======================================================================
# 1. Structural validity
# ======================================================================
test_that("as.SummarizedExperiment produces valid SE with correct dimensions", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)

  expect_true(methods::validObject(se))
  expect_equal(nrow(se), nrow(tdr@density$fdens))
  expect_equal(ncol(se), ncol(tdr@density$fdens))
  expect_true("normcounts" %in% SummarizedExperiment::assayNames(se))
  expect_true("logcounts"  %in% SummarizedExperiment::assayNames(se))
})

test_that("counts assay present when fuzzy graphs available", {
  tdr <- create_mapped_tdr(with_fuzzy_graphs = TRUE)
  se  <- as.SummarizedExperiment(tdr)
  expect_true("counts" %in% SummarizedExperiment::assayNames(se))
})


# ======================================================================
# 2. Assay fidelity
# ======================================================================
test_that("normcounts matches fdens exactly", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_identical(
    SummarizedExperiment::assay(se, "normcounts"),
    as.matrix(tdr@density$fdens)
  )
})

test_that("logcounts matches Y exactly", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_identical(
    SummarizedExperiment::assay(se, "logcounts"),
    as.matrix(tdr@density$Y)
  )
})

test_that("assay dimnames match density dimnames", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_identical(rownames(SummarizedExperiment::assay(se, "normcounts")),
                   rownames(tdr@density$fdens))
  expect_identical(colnames(SummarizedExperiment::assay(se, "normcounts")),
                   colnames(tdr@density$fdens))
})


# ======================================================================
# 3. Counts assay correctness
# ======================================================================
test_that("counts equals per-sample rowSums of fuzzy graphs", {
  tdr <- create_mapped_tdr(with_fuzzy_graphs = TRUE)
  se  <- as.SummarizedExperiment(tdr)

  counts <- SummarizedExperiment::assay(se, "counts")
  lm_names <- rownames(tdr@density$fdens)

  for (sn in names(tdr@cellmap$fuzzy.graphs)) {
    fg <- tdr@cellmap$fuzzy.graphs[[sn]]
    if (identical(rownames(fg), lm_names)) {
      expected <- Matrix::rowSums(fg)
    } else {
      expected <- Matrix::colSums(fg)
    }
    expect_equal(unname(counts[, sn]), unname(as.numeric(expected)),
                 tolerance = 1e-10)
  }
})


# ======================================================================
# 4. rowData correctness
# ======================================================================
test_that("rowData contains clustering and celltyping columns", {
  tdr <- create_mapped_tdr(with_celltyping = TRUE)
  se  <- as.SummarizedExperiment(tdr)
  rd  <- SummarizedExperiment::rowData(se)

  expect_equal(nrow(rd), nrow(tdr@density$fdens))

  # Clustering
  cl_names <- setdiff(names(tdr@landmark.annot$clustering), "ids")
  for (nm in cl_names) {
    col_nm <- paste0("clustering_", nm)
    expect_true(col_nm %in% colnames(rd), info = paste("Missing column:", col_nm))
    expect_true(is.factor(rd[[col_nm]]))
  }
  expect_true("clustering_active" %in% colnames(rd))

  # Celltyping
  ct_names <- setdiff(names(tdr@landmark.annot$celltyping), "ids")
  for (nm in ct_names) {
    col_nm <- paste0("celltyping_", nm)
    expect_true(col_nm %in% colnames(rd), info = paste("Missing column:", col_nm))
    expect_true(is.factor(rd[[col_nm]]))
  }
  expect_true("celltyping_active" %in% colnames(rd))
})

test_that("factor levels preserved in rowData", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  rd  <- SummarizedExperiment::rowData(se)
  expect_identical(
    levels(rd[["clustering_active"]]),
    levels(tdr@landmark.annot$clustering$ids)
  )
})


# ======================================================================
# 5. colData correctness
# ======================================================================
test_that("colData matches metadata", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  expect_equal(cd, tdr@metadata)
  expect_identical(rownames(SummarizedExperiment::colData(se)),
                   colnames(tdr@density$fdens))
})


# ======================================================================
# 6. Metadata round-trip
# ======================================================================
test_that("GetTDR round-trips the TDRObj", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)

  tdr2 <- GetTDR(se)
  expect_true(is.TDRObj(tdr2))
  expect_identical(tdr2@density$fdens, tdr@density$fdens)
  expect_identical(tdr2@density$Y, tdr@density$Y)
})

test_that("metadata contains version and date", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_false(is.null(S4Vectors::metadata(se)$tinydenseR_version))
  expect_false(is.null(S4Vectors::metadata(se)$conversion_date))
})


# ======================================================================
# 7. Dispatch wrappers
# ======================================================================
test_that("plotPCA dispatches on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  p   <- plotPCA(se)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotUMAP dispatches on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  p   <- plotUMAP(se)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("get.landmarks errors on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_error(get.landmarks(se), "requires cell-level data")
})

test_that("get.map errors on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_error(get.map(se), "requires cell-level data")
})

test_that("get.pbDE errors on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_error(get.pbDE(se), "requires cell-level data")
})

test_that("get.markerDE errors on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_error(get.markerDE(se), "requires cell-level data")
})

test_that("goi.summary errors on SummarizedExperiment", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)
  expect_error(goi.summary(se), "requires cell-level data")
})


# ======================================================================
# 8. Fail-fast validation
# ======================================================================
test_that("conversion fails when fdens is NULL", {
  tdr <- create_mapped_tdr()
  tdr@density$fdens <- NULL
  expect_error(as.SummarizedExperiment(tdr),
               "get.map\\(\\) must be run")
})

test_that("conversion fails when Y is NULL", {
  tdr <- create_mapped_tdr()
  tdr@density$Y <- NULL
  expect_error(as.SummarizedExperiment(tdr),
               "get.map\\(\\) must be run")
})

test_that("conversion fails when clustering$ids is NULL", {
  tdr <- create_mapped_tdr()
  tdr@landmark.annot$clustering$ids <- NULL
  expect_error(as.SummarizedExperiment(tdr),
               "get.map\\(\\) must be run")
})

test_that("conversion fails on metadata row mismatch", {
  tdr <- create_mapped_tdr(n_samples = 4)
  tdr@metadata <- tdr@metadata[1:2, , drop = FALSE]
  expect_error(as.SummarizedExperiment(tdr), "Mismatch")
})


# ======================================================================
# 9. Graceful degradation — no fuzzy graphs
# ======================================================================
test_that("SE created without counts when fuzzy graphs are NULL", {
  tdr <- create_mapped_tdr(with_fuzzy_graphs = FALSE)
  expect_message(
    se <- as.SummarizedExperiment(tdr),
    NA  # no message expected; counts simply omitted
  )
  expect_false("counts" %in% SummarizedExperiment::assayNames(se))
  expect_true("normcounts" %in% SummarizedExperiment::assayNames(se))
})

test_that("SE created without counts when fuzzy graph cache expired", {
  tdr <- create_mapped_tdr(with_fuzzy_graphs = FALSE)
  # Simulate cache paths that don't exist
  tdr@cellmap$fuzzy.graphs <- stats::setNames(
    as.list(paste0("/nonexistent/path/", colnames(tdr@density$fdens), ".rds")),
    colnames(tdr@density$fdens)
  )
  expect_message(
    se <- as.SummarizedExperiment(tdr),
    "Fuzzy graph cache unavailable"
  )
  expect_false("counts" %in% SummarizedExperiment::assayNames(se))
})


# ======================================================================
# 10. SetTDR.SummarizedExperiment
# ======================================================================
test_that("SetTDR stores and retrieves TDRObj in SE", {
  tdr <- create_mapped_tdr()
  se  <- as.SummarizedExperiment(tdr)

  # Modify the TDRObj and store back
  tdr2 <- GetTDR(se)
  tdr2@config$modified <- TRUE
  se2  <- SetTDR(se, tdr2)

  expect_true(GetTDR(se2)@config$modified)
})

test_that("GetTDR errors on empty SE", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(x = matrix(1, 2, 2))
  )
  expect_error(GetTDR(se), "No TDRObj found")
})

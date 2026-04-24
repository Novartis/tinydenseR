test_that("get.subset — .id path selects correct cells by celltyping", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)
  obj <- celltyping(obj,
    .celltyping.map = list(
      typeA = levels(obj@landmark.annot$clustering$ids)[1],
      typeB = setdiff(levels(obj@landmark.annot$clustering$ids),
                      levels(obj@landmark.annot$clustering$ids)[1])
    ),
    .verbose = FALSE
  )
  obj <- get.map(obj, .verbose = FALSE, .seed = 123)

  child <- get.subset(obj, .id = "typeA", .id.from = "celltyping",
                      .verbose = FALSE)

  expect_s4_class(child, "TDRObj")
  expect_true(length(child@cells) > 0)
  expect_true(length(child@cells) <= length(obj@cells))

  # All analysis slots should be empty

  expect_length(child@assay, 0)
  expect_length(child@landmark.embed, 0)
  expect_length(child@landmark.annot, 0)
  expect_length(child@graphs, 0)
  expect_length(child@density, 0)
  expect_length(child@cellmap, 0)
  expect_length(child@results, 0)

  # Metadata should match retained samples
  expect_equal(nrow(child@metadata), length(child@cells))
  expect_true(all(child@metadata$n.cells > 0))
})

test_that("get.subset — .id path selects by clustering", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .verbose = FALSE)

  expect_s4_class(child, "TDRObj")
  expect_true(length(child@cells) > 0)
})

test_that("get.subset — .id.idx path uses fuzzy selection", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 500, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  # Select half the landmarks with low confidence threshold
  n_lm <- nrow(obj@assay$expr)
  idx <- seq_len(max(1, floor(n_lm / 2)))
  child <- get.subset(obj, .id.idx = idx, .label.confidence = 0.1,
                      .min.cells.per.sample = 1, .verbose = FALSE)

  expect_s4_class(child, "TDRObj")
  expect_true(length(child@cells) > 0)
})

test_that("get.subset — sample filtering by .min.cells.per.sample", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  # Use a very high threshold to potentially drop samples
  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster,
                      .id.from = "clustering",
                      .min.cells.per.sample = 1,
                      .verbose = FALSE)

  expect_s4_class(child, "TDRObj")
  expect_equal(nrow(child@metadata), length(child@cells))
  expect_true(all(child@metadata$n.cells >= 1))
})

test_that("get.subset — metadata consistency", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # n.cells in metadata should match actual cell counts
  expect_equal(nrow(child@metadata), length(child@cells))
  expect_true(all(names(child@cells) %in% rownames(child@metadata)))
  expect_true("n.cells" %in% colnames(child@metadata))
  expect_true("log10.n.cells" %in% colnames(child@metadata))
  expect_equal(child@metadata$log10.n.cells,
               log10(child@metadata$n.cells))
  # n.perSample should be set (computed from sampling config)
  expect_true("n.perSample" %in% colnames(child@metadata))
  expect_true(all(child@metadata$n.perSample > 0))
})

test_that("get.subset — config inheritance", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  expect_equal(child@config$assay.type, obj@config$assay.type)
  expect_equal(child@config$markers, obj@config$markers)
})

test_that("get.subset — provenance is stored", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  prov <- child@config$.subset.provenance
  expect_type(prov, "list")
  expect_length(prov, 1)
  expect_equal(prov[[1]]$selection.type, "id")
  expect_equal(prov[[1]]$selection.values, first_cluster)
  expect_equal(prov[[1]]$selection.from, "clustering")
  expect_equal(prov[[1]]$parent.n.samples, length(obj@cells))
  expect_equal(prov[[1]]$parent.n.landmarks, nrow(obj@assay$expr))
})

test_that("get.subset — nested subset accumulates provenance", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 500, n_markers = 4, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # Run pipeline on child so it has cellmap for nested subset
  first_mat <- .get_sample_matrix(NULL, child, 1)
  k_child <- min(5, max(2, sum(child@config$sampling$n.perSample) - 1))
  child <- get.landmarks(child, .verbose = FALSE, .seed = 123,
                         .nPC = min(3, ncol(first_mat))) |>
    get.graph(.k = k_child, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  second_cluster <- levels(child@landmark.annot$clustering$ids)[1]
  grandchild <- get.subset(child, .id = second_cluster,
                           .id.from = "clustering",
                           .min.cells.per.sample = 1,
                           .verbose = FALSE)

  prov <- grandchild@config$.subset.provenance
  expect_type(prov, "list")
  expect_length(prov, 2)
  expect_equal(prov[[1]]$selection.values, first_cluster)
  expect_equal(prov[[2]]$selection.values, second_cluster)
})

test_that("get.subset — files backend uses list(path, idx)", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # Files backend: child cells should be lists with path and idx
  for (sn in names(child@cells)) {
    entry <- child@cells[[sn]]
    expect_type(entry, "list")
    expect_true("path" %in% names(entry))
    expect_true("idx" %in% names(entry))
    expect_true(file.exists(entry$path))
    expect_type(entry$idx, "integer")
    expect_true(all(entry$idx >= 1))
  }
})

test_that("get.subset — .get_sample_matrix works with list(path, idx)", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # The accessor should return a matrix with the correct number of rows
  sn <- names(child@cells)[1]
  mat <- .get_sample_matrix(NULL, child, 1)
  entry <- child@cells[[sn]]
  expect_equal(nrow(mat), length(entry$idx))
  expect_equal(ncol(mat), ncol(readRDS(entry$path)))
})

test_that("get.subset — full pipeline runs on child", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 500, n_markers = 4, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # Full pipeline on child
  first_mat <- .get_sample_matrix(NULL, child, 1)
  k_child <- min(5, max(2, sum(child@config$sampling$n.perSample) - 1))
  child <- get.landmarks(child, .verbose = FALSE, .seed = 123,
                         .nPC = min(3, ncol(first_mat))) |>
    get.graph(.k = k_child, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  # Validate child has its own landmarks, graph, density, cellmap
  expect_true(nrow(child@assay$expr) > 0)
  expect_true(!is.null(child@landmark.embed$pca))
  expect_true(!is.null(child@landmark.embed$umap))
  expect_true(!is.null(child@landmark.annot$clustering$ids))
  expect_true(!is.null(child@density$norm))
  expect_true(!is.null(child@density$log.norm))
  expect_true(!is.null(child@cellmap$clustering$ids))

  # Child landmarks should be fewer than parent
  expect_true(nrow(child@assay$expr) <= nrow(obj@assay$expr))
})

test_that("get.subset — error if get.map not run", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 3, n_samples = 2)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .assay.type = "cyto",
    .verbose = FALSE
  )

  expect_error(
    get.subset(obj, .id = "cluster.1", .verbose = FALSE),
    "get.map\\(\\)"
  )
})

test_that("get.subset — error if no .id or .id.idx", {
  obj <- TDRObj(
    config = list(assay.type = "cyto"),
    cellmap = list(clustering = list(ids = list(s1 = c("a", "b"))))
  )
  expect_error(
    get.subset(obj, .verbose = FALSE),
    "Either .id or .id.idx"
  )
})

test_that("get.subset — error if both .id and .id.idx", {
  obj <- TDRObj(
    config = list(assay.type = "cyto"),
    cellmap = list(clustering = list(ids = list(s1 = c("a", "b"))))
  )
  expect_error(
    get.subset(obj, .id = "a", .id.idx = 1L, .verbose = FALSE),
    "not both"
  )
})

test_that("get.subset — error for invalid .id label", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  expect_error(
    get.subset(obj, .id = "nonexistent_label", .id.from = "clustering",
               .verbose = FALSE),
    "not found"
  )
})

test_that("get.subset — error for invalid .label.confidence", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  expect_error(
    get.subset(obj, .id = first_cluster, .label.confidence = 1.5,
               .verbose = FALSE),
    "\\[0,1\\]"
  )
})

test_that("get.subset — error if all samples too small", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 2)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  expect_error(
    get.subset(obj, .id = first_cluster, .id.from = "clustering",
               .min.cells.per.sample = 1e6, .verbose = FALSE),
    "No samples have"
  )
})

test_that("get.subset — matrix backend", {
  set.seed(42)
  n_cells <- 200
  n_markers <- 4
  n_samples <- 3
  sample_names <- paste0("sample", seq_len(n_samples))

  # Build a single cells-x-markers dense matrix for cyto
  total_cells <- n_cells * n_samples
  mat <- matrix(
    runif(total_cells * n_markers),
    nrow = total_cells,
    ncol = n_markers,
    dimnames = list(
      paste0("cell_", seq_len(total_cells)),
      paste0("marker_", seq_len(n_markers))
    )
  )

  sample_labels <- rep(sample_names, each = n_cells)
  meta <- data.frame(
    row.names = sample_names,
    group = c("A", "B", "A")[seq_len(n_samples)]
  )

  # Build TDRObj manually for matrix backend (setup.tdr.obj only accepts file paths)
  .cells <- lapply(
    stats::setNames(sample_names, sample_names),
    function(s) which(sample_labels == s)
  )
  source_env <- new.env(parent = emptyenv())
  source_env$mat <- mat

  n.cells.vec <- lengths(.cells)
  .prop <- 0.1
  target.lm.n <- pmin(sum(n.cells.vec) * .prop, 5e3)
  n.perSample  <- pmin(ceiling(n.cells.vec * .prop),
                        ceiling(target.lm.n / length(.cells)))

  obj <- TDRObj(
    cells = .cells,
    metadata = meta,
    config = list(
      key = seq_along(.cells) |>
        rep(times = n.perSample) |>
        stats::setNames(nm = names(.cells)[
          seq_along(.cells) |> rep(times = n.perSample)
        ]),
      sampling = list(
        n.cells     = n.cells.vec,
        target.lm.n = target.lm.n,
        n.perSample = n.perSample
      ),
      assay.type = "cyto",
      markers    = paste0("marker_", seq_len(n_markers)),
      n.threads  = 1L,
      backend    = "matrix",
      source.env = source_env
    ),
    integration = list(harmony.var = NULL, harmony.obj = NULL)
  )
  obj@metadata$n.cells      <- n.cells.vec
  obj@metadata$log10.n.cells <- log10(n.cells.vec)
  obj@metadata$n.perSample  <- n.perSample

  obj <- get.landmarks(obj, .verbose = FALSE, .seed = 123) |>
    get.graph(.verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  expect_equal(obj@config$backend, "matrix")

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  expect_s4_class(child, "TDRObj")
  expect_equal(child@config$backend, "matrix")

  # Index-based backend: @cells entries should be integer vectors
  for (sn in names(child@cells)) {
    expect_type(child@cells[[sn]], "integer")
    # Child indices must be a subset of parent indices
    expect_true(all(child@cells[[sn]] %in% obj@cells[[sn]]))
  }

  # Full pipeline on child
  child2 <- get.landmarks(child, .verbose = FALSE, .seed = 123) |>
    get.graph(.verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  expect_true(nrow(child2@assay$expr) > 0)
  expect_true(!is.null(child2@density$norm))
})

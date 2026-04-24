test_that("no-copy: files backend — child reuses same RDS paths", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]

  # Count RDS files before
  rds_before <- list.files(tempdir(), pattern = "\\.RDS$", full.names = TRUE)

  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # Count RDS files after — should be the same (no new files created)
  rds_after <- list.files(tempdir(), pattern = "\\.RDS$", full.names = TRUE)
  new_rds <- setdiff(rds_after, rds_before)
  expect_length(new_rds, 0)

  # Child paths should be identical to parent paths
  for (sn in names(child@cells)) {
    child_entry <- child@cells[[sn]]
    parent_path <- obj@cells[[sn]]
    expect_identical(child_entry$path, parent_path)
  }
})

test_that("no-copy: child object size is much smaller than parent", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  parent_size <- as.numeric(object.size(obj))
  child_size <- as.numeric(object.size(child))

  # Child should be much smaller (no landmarks, graphs, density, results)
  expect_true(child_size < parent_size * 0.5)
})

test_that("no-copy: matrix backend — child shares source.env", {
  set.seed(42)
  n_cells <- 200
  n_markers <- 4
  n_samples <- 3
  sample_names <- paste0("sample", seq_len(n_samples))

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

  # Build TDRObj manually for matrix backend
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

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # Source environment should be shared (same pointer, not a copy)
  expect_identical(child@config$source.env, obj@config$source.env)
  expect_identical(child@config$source.env$mat, obj@config$source.env$mat)
})

test_that("no-copy: .get_sample_matrix returns correct subset for files", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 200, n_markers = 4, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  for (i in seq_along(child@cells)) {
    sn <- names(child@cells)[i]
    child_mat <- .get_sample_matrix(NULL, child, i)
    entry <- child@cells[[sn]]

    # Number of rows should match the subset index count
    expect_equal(nrow(child_mat), length(entry$idx))

    # The subset should be actual rows from the parent matrix
    parent_mat <- readRDS(entry$path)
    expected_mat <- parent_mat[entry$idx, , drop = FALSE]
    expect_equal(child_mat, expected_mat)
  }
})

test_that("no-copy: nested subset still references original files", {
  set.seed(42)
  test_data <- create_test_lm_obj(n_cells = 500, n_markers = 4, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  obj <- run_full_pipeline(test_data, verbose = FALSE)
  original_paths <- vapply(obj@cells, identity, character(1))

  first_cluster <- levels(obj@landmark.annot$clustering$ids)[1]
  child <- get.subset(obj, .id = first_cluster, .id.from = "clustering",
                      .min.cells.per.sample = 1, .verbose = FALSE)

  # Run pipeline on child
  child <- get.landmarks(child, .verbose = FALSE, .seed = 123) |>
    get.graph(.verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  second_cluster <- levels(child@landmark.annot$clustering$ids)[1]
  grandchild <- get.subset(child, .id = second_cluster,
                           .id.from = "clustering",
                           .min.cells.per.sample = 1,
                           .verbose = FALSE)

  # Grandchild should still reference the original RDS files
  for (sn in names(grandchild@cells)) {
    expect_true(grandchild@cells[[sn]]$path %in% original_paths)
  }

  # Grandchild indices should be a subset of child indices
  for (sn in names(grandchild@cells)) {
    gc_entry <- grandchild@cells[[sn]]
    ch_entry <- child@cells[[sn]]
    # The grandchild's indices (into the original file) should be within
    # the child's indices (also into the original file)
    expect_true(all(gc_entry$idx %in% ch_entry$idx))
  }
})

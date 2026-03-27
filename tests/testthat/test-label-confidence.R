test_that(".tdr_validate_label_confidence enforces [0,1] scalar finite numeric", {
  expect_no_error(tinydenseR:::.tdr_validate_label_confidence(0))
  expect_no_error(tinydenseR:::.tdr_validate_label_confidence(0.5))
  expect_no_error(tinydenseR:::.tdr_validate_label_confidence(1))

  expect_error(tinydenseR:::.tdr_validate_label_confidence(NA_real_), "\\[0,1\\]")
  expect_error(tinydenseR:::.tdr_validate_label_confidence(Inf), "\\[0,1\\]")
  expect_error(tinydenseR:::.tdr_validate_label_confidence(-1e-6), "\\[0,1\\]")
  expect_error(tinydenseR:::.tdr_validate_label_confidence(1 + 1e-6), "\\[0,1\\]")
  expect_error(tinydenseR:::.tdr_validate_label_confidence(c(0.2, 0.3)), "single finite numeric")
  expect_error(tinydenseR:::.tdr_validate_label_confidence("0.5"), "single finite numeric")
})

test_that(".tdr_transfer_labels fuzzy mode computes expected thresholded assignments", {
  fgraph <- Matrix::Matrix(
    data = c(
      0.80, 0.20, 0.00,
      0.40, 0.60, 0.00,
      0.49, 0.51, 0.00
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )

  labels <- c("A", "B", "B")

  pred_06 <- tinydenseR:::.tdr_transfer_labels(
    .method = "fuzzy",
    .label.confidence = 0.6,
    .n.cells = 3,
    .cell.names = c("c1", "c2", "c3"),
    .fgraph = fgraph,
    .landmark.labels = labels
  )

  expect_identical(unname(pred_06), c("A", "B", "..low.confidence.."))

  pred_00 <- tinydenseR:::.tdr_transfer_labels(
    .method = "fuzzy",
    .label.confidence = 0,
    .n.cells = 3,
    .cell.names = c("c1", "c2", "c3"),
    .fgraph = fgraph,
    .landmark.labels = labels
  )

  expect_true(all(pred_00 %in% c("A", "B")))

  pred_10 <- tinydenseR:::.tdr_transfer_labels(
    .method = "fuzzy",
    .label.confidence = 1,
    .n.cells = 3,
    .cell.names = c("c1", "c2", "c3"),
    .fgraph = fgraph,
    .landmark.labels = labels
  )

  expect_identical(unname(pred_10), rep("..low.confidence..", 3))
})

test_that(".tdr_transfer_labels knn_vote mode computes vote-frequency confidence", {
  nn.idx <- matrix(
    c(
      1, 1, 1, 1, 1, 2, 2, 2, 3, 4,
      2, 2, 2, 2, 2, 2, 1, 3, 4, 5,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )

  ref.labels <- c("A", "B", "A", "B", "B")

  pred_06 <- tinydenseR:::.tdr_transfer_labels(
    .method = "knn_vote",
    .label.confidence = 0.6,
    .n.cells = 3,
    .cell.names = c("k1", "k2", "k3"),
    .nn.idx = nn.idx,
    .ref.labels = ref.labels
  )

  expect_identical(unname(pred_06), c("A", "B", "A"))

  pred_09 <- tinydenseR:::.tdr_transfer_labels(
    .method = "knn_vote",
    .label.confidence = 0.9,
    .n.cells = 3,
    .cell.names = c("k1", "k2", "k3"),
    .nn.idx = nn.idx,
    .ref.labels = ref.labels
  )

  expect_identical(unname(pred_09), c("..low.confidence..", "..low.confidence..", "A"))
})

test_that("deprecated wrappers validate .label.confidence before dispatch", {
  expect_error(
    tinydenseR::get.dea(.tdr.obj = NULL, .design = matrix(1, 1, 1), .label.confidence = -0.1),
    "\\[0,1\\]"
  )

  expect_error(
    tinydenseR::get.marker(.tdr.obj = NULL, .label.confidence = 1.1),
    "\\[0,1\\]"
  )
})

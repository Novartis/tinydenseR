#' @noRd
#' @keywords internal
.tdr_validate_label_confidence <- function(.label.confidence,
                                           .arg = ".label.confidence") {
  if (!is.numeric(x = .label.confidence) ||
      length(x = .label.confidence) != 1 ||
      is.na(x = .label.confidence) ||
      !is.finite(x = .label.confidence) ||
      .label.confidence < 0 ||
      .label.confidence > 1) {
    stop(
      .arg,
      " must be a single finite numeric value in [0,1].",
      call. = FALSE
    )
  }

  invisible(.label.confidence)
}

#' @noRd
#' @keywords internal
.tdr_transfer_labels <- function(
  .method = c("fuzzy", "knn_vote"),
  .label.confidence,
  .n.cells,
  .cell.names,
  .default.label = "..low.confidence..",
  .fgraph = NULL,
  .landmark.labels = NULL,
  .nn.idx = NULL,
  .ref.labels = NULL
) {
  i <- j <- cell <- landmark <- label <- x <- confidence <- ref.idx <- N <- NULL

  .method <- match.arg(arg = .method)
  .tdr_validate_label_confidence(.label.confidence)

  if (!is.numeric(x = .n.cells) || length(x = .n.cells) != 1 || .n.cells < 1) {
    stop(".n.cells must be a positive scalar.", call. = FALSE)
  }

  if (length(x = .cell.names) != .n.cells) {
    stop("length(.cell.names) must equal .n.cells.", call. = FALSE)
  }

  if (.method == "fuzzy") {
    if (is.null(x = .fgraph) || is.null(x = .landmark.labels)) {
      stop("For method = 'fuzzy', .fgraph and .landmark.labels are required.", call. = FALSE)
    }

    if (length(x = .landmark.labels) != ncol(x = .fgraph)) {
      stop("length(.landmark.labels) must equal ncol(.fgraph).", call. = FALSE)
    }

    votes <-
      Matrix::summary(object = .fgraph) |>
      dplyr::rename(cell = i,
                    landmark = j) |>
      dplyr::mutate(label = as.character(x = .landmark.labels[landmark])) |>
      dplyr::select(-landmark) |>
      collapse::fgroup_by(cell,
                          label,
                          sort = FALSE) |>
      collapse::fsum() |>
      collapse::fgroup_by(cell,
                          sort = FALSE) |>
      collapse::fmutate(confidence = x / collapse::fsum(x)) |>
      collapse::fungroup() |>
      collapse::fsubset(confidence >= .label.confidence) |>
      collapse::roworderv(cols = "confidence",
                          decreasing = TRUE) |>
      collapse::fgroup_by(cell,
                          sort = FALSE) |>
      collapse::ffirst()

  } else {
    if (is.null(x = .nn.idx) || is.null(x = .ref.labels)) {
      stop("For method = 'knn_vote', .nn.idx and .ref.labels are required.", call. = FALSE)
    }

    if (nrow(x = .nn.idx) != .n.cells) {
      stop("nrow(.nn.idx) must equal .n.cells.", call. = FALSE)
    }

    votes <-
      data.frame(cell = rep(x = seq_len(length.out = .n.cells),
                            times = ncol(x = .nn.idx)),
                 ref.idx = as.vector(x = .nn.idx)) |>
      dplyr::mutate(label = as.character(x = .ref.labels[ref.idx])) |>
      dplyr::select(-ref.idx) |>
      collapse::fgroup_by(cell,
                          label,
                          sort = FALSE) |>
      collapse::fcount() |>
      collapse::fgroup_by(cell,
                          label,
                          sort = FALSE) |>
      collapse::fmutate(confidence = N / ncol(x = .nn.idx)) |>
      collapse::fungroup() |>
      collapse::fsubset(confidence >= .label.confidence) |>
      collapse::roworderv(cols = "confidence",
                          decreasing = TRUE) |>
      collapse::fgroup_by(cell,
                          sort = FALSE) |>
      collapse::ffirst()
  }

  out <- rep(x = .default.label,
             times = .n.cells)

  if (nrow(x = votes) > 0) {
    out[votes$cell] <- votes$label
  }

  stats::setNames(object = out,
                  nm = .cell.names)
}

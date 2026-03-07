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

# ──────────────────────────────────────────────────────────────────────
# HDF5-level utilities for subsetting .h5ad files without full
# materialization.
# ──────────────────────────────────────────────────────────────────────

#' Detect the encoding-type attribute on an HDF5 group/dataset
#'
#' @param file Character scalar. Path to the HDF5 file.
#' @param path Character scalar. HDF5 path within the file.
#' @return Character scalar (e.g., \code{"csr_matrix"}) or
#'   \code{NA_character_} if the attribute is missing.
#' @keywords internal
#' @noRd
.h5ad_detect_encoding <- function(file, path) {
  attrs <- rhdf5::h5readAttributes(file, path)
  enc <- attrs[["encoding-type"]]
  if (is.null(enc)) NA_character_ else as.character(enc)
}

#' Copy all HDF5 attributes from one path to another
#'
#' @param src_file Character scalar. Source HDF5 file path.
#' @param dest_file Character scalar. Destination HDF5 file path.
#' @param path Character scalar. HDF5 path whose attributes to copy.
#' @keywords internal
#' @noRd
.h5ad_copy_attrs <- function(src_file, dest_file, path) {
  attrs <- rhdf5::h5readAttributes(src_file, path)
  for (nm in names(attrs)) {
    val <- attrs[[nm]]
    # rhdf5 cannot write zero-length character vectors via h5writeAttribute;
    # fall back to the low-level API in that case.
    if (is.character(val) && length(val) == 0L) {
      .h5ad_write_empty_string_attr(dest_file, path, nm)
      next
    }
    is_scalar <- is.atomic(val) && length(val) == 1L
    rhdf5::h5writeAttribute(val, dest_file, nm,
                            h5loc = path, asScalar = is_scalar)
  }
}

#' Write an empty character-vector attribute via the low-level HDF5 API
#'
#' \code{rhdf5::h5writeAttribute(character(0), ...)} errors out;
#' this helper creates a zero-length fixed-string attribute directly.
#'
#' @keywords internal
#' @noRd
.h5ad_write_empty_string_attr <- function(file, h5path, attr_name) {
  fid <- rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(fid), add = TRUE)

  loc <- if (h5path == "/") {
    fid
  } else {
    obj <- rhdf5::H5Oopen(fid, h5path)
    on.exit(rhdf5::H5Oclose(obj), add = TRUE)
    obj
  }

  tid <- rhdf5::H5Tcopy("H5T_C_S1")
  rhdf5::H5Tset_size(tid, 1L)
  sid <- rhdf5::H5Screate_simple(0)
  aid <- rhdf5::H5Acreate(loc, attr_name, tid, sid)
  rhdf5::H5Aclose(aid)
  rhdf5::H5Sclose(sid)
}

#' Copy a dataframe-encoded HDF5 group (/obs or /var), subsetting rows
#'
#' Handles regular array columns, categorical columns (codes + categories),
#' and the index column. Copies all HDF5 attributes.
#'
#' @param src_file Character scalar. Source file path.
#' @param dest_file Character scalar. Destination file path.
#' @param src_path Character scalar. HDF5 path (e.g., \code{"/obs"}).
#' @param idx Integer vector. 1-based row indices to keep.
#' @param compression Integer 0-9. gzip level.
#' @keywords internal
#' @noRd
.h5ad_copy_dataframe_group <- function(src_file, dest_file, src_path,
                                       idx, compression) {
  rhdf5::h5createGroup(dest_file, src_path)

  # List children of this group
  ls_df <- rhdf5::h5ls(src_file, recursive = FALSE)
  # Filter to children of src_path
  ls_all <- rhdf5::h5ls(src_file, recursive = TRUE)
  # Direct children: group == src_path (or "/" if src_path == "/")
  parent <- if (src_path == "/") "/" else src_path
  children <- ls_all[ls_all$group == parent, , drop = FALSE]

  for (i in seq_len(nrow(children))) {
    child_name <- children$name[i]
    child_type <- children$otype[i]
    child_path <- paste0(src_path, "/", child_name)

    if (child_type == "H5I_GROUP") {
      # Categorical column: has codes + categories
      enc <- .h5ad_detect_encoding(src_file, child_path)
      if (!is.na(enc) && enc == "categorical") {
        rhdf5::h5createGroup(dest_file, child_path)

        # Subset codes by idx
        codes <- rhdf5::h5read(src_file, paste0(child_path, "/codes"))
        if (is.factor(codes)) codes <- as.integer(codes) - 1L
        rhdf5::h5write(codes[idx], dest_file, paste0(child_path, "/codes"))

        # Copy categories as-is
        cats <- rhdf5::h5read(src_file, paste0(child_path, "/categories"))
        if (is.factor(cats)) cats <- as.character(cats)
        rhdf5::h5write(cats, dest_file, paste0(child_path, "/categories"))

        # Copy attributes on the categorical group
        .h5ad_copy_attrs(src_file, dest_file, child_path)
      } else {
        # Unknown sub-group inside obs/var — deep copy
        .h5ad_copy_group_recursive(src_file, dest_file, child_path)
      }
    } else {
      # H5I_DATASET — regular array column or the index column
      data <- rhdf5::h5read(src_file, child_path)
      # rhdf5 can return factors for enum-typed datasets; convert back
      if (is.factor(data)) data <- as.character(data)
      rhdf5::h5write(data[idx], dest_file, child_path)
      # Copy any attributes on this dataset
      .h5ad_copy_attrs(src_file, dest_file, child_path)
    }
  }

  # Copy attributes on the dataframe group itself
  .h5ad_copy_attrs(src_file, dest_file, src_path)
}

#' Recursively copy an HDF5 group without subsetting
#'
#' Used for groups like /uns that are copied as-is.
#'
#' @param src_file Character scalar.
#' @param dest_file Character scalar.
#' @param path Character scalar. HDF5 path of the group.
#' @keywords internal
#' @noRd
.h5ad_copy_group_recursive <- function(src_file, dest_file, path) {
  fid <- rhdf5::H5Fopen(dest_file)
  exists <- rhdf5::H5Lexists(fid, path)
  rhdf5::H5Fclose(fid)
  if (!exists) rhdf5::h5createGroup(dest_file, path)

  ls_all <- rhdf5::h5ls(src_file, recursive = TRUE)
  children <- ls_all[ls_all$group == path, , drop = FALSE]

  for (i in seq_len(nrow(children))) {
    child_name <- children$name[i]
    child_type <- children$otype[i]
    child_path <- paste0(path, "/", child_name)

    if (child_type == "H5I_GROUP") {
      .h5ad_copy_group_recursive(src_file, dest_file, child_path)
    } else {
      data <- rhdf5::h5read(src_file, child_path)
      if (is.factor(data)) data <- as.character(data)
      rhdf5::h5write(data, dest_file, child_path)
      .h5ad_copy_attrs(src_file, dest_file, child_path)
    }
  }

  .h5ad_copy_attrs(src_file, dest_file, path)
}

#' Copy an embedding group (/obsm or /varm), subsetting rows
#'
#' Each child of obsm/varm is a dense 2-D array. Rows are subset by idx.
#'
#' @param src_file Character scalar.
#' @param dest_file Character scalar.
#' @param group_path Character scalar (e.g., \code{"/obsm"}).
#' @param idx Integer vector. 1-based row indices.
#' @param compression Integer 0-9.
#' @keywords internal
#' @noRd
.h5ad_copy_embedding_group <- function(src_file, dest_file, group_path,
                                       idx, compression) {
  fid <- rhdf5::H5Fopen(src_file)
  exists <- rhdf5::H5Lexists(fid, group_path)
  rhdf5::H5Fclose(fid)
  if (!exists) return(invisible(NULL))

  rhdf5::h5createGroup(dest_file, group_path)

  ls_all <- rhdf5::h5ls(src_file, recursive = TRUE)
  children <- ls_all[ls_all$group == group_path, , drop = FALSE]

  for (i in seq_len(nrow(children))) {
    child_name <- children$name[i]
    child_type <- children$otype[i]
    child_path <- paste0(group_path, "/", child_name)

    if (child_type == "H5I_DATASET") {
      mat <- rhdf5::h5read(src_file, child_path)
      if (is.matrix(mat)) {
        rhdf5::h5write(mat[idx, , drop = FALSE], dest_file, child_path)
      } else {
        # 1-D or other; subset directly
        rhdf5::h5write(mat[idx], dest_file, child_path)
      }
      .h5ad_copy_attrs(src_file, dest_file, child_path)
    } else {
      # Sparse embedding — handle like a matrix group
      .h5ad_copy_group_recursive(src_file, dest_file, child_path)
    }
  }

  .h5ad_copy_attrs(src_file, dest_file, group_path)
}

#' Copy a pairwise group (/obsp or /varp), subsetting rows and columns
#'
#' @param src_file Character scalar.
#' @param dest_file Character scalar.
#' @param group_path Character scalar (e.g., \code{"/obsp"}).
#' @param idx Integer vector. 1-based indices for both dimensions.
#' @param compression Integer 0-9.
#' @keywords internal
#' @noRd
.h5ad_copy_pairwise_group <- function(src_file, dest_file, group_path,
                                      idx, compression) {
  fid <- rhdf5::H5Fopen(src_file)
  exists <- rhdf5::H5Lexists(fid, group_path)
  rhdf5::H5Fclose(fid)
  if (!exists) return(invisible(NULL))

  rhdf5::h5createGroup(dest_file, group_path)

  ls_all <- rhdf5::h5ls(src_file, recursive = TRUE)
  children <- ls_all[ls_all$group == group_path, , drop = FALSE]

  for (i in seq_len(nrow(children))) {
    child_name <- children$name[i]
    child_type <- children$otype[i]
    child_path <- paste0(group_path, "/", child_name)

    enc <- .h5ad_detect_encoding(src_file, child_path)
    if (!is.na(enc) && enc %in% c("csr_matrix", "csc_matrix")) {
      # Sparse pairwise matrix: subset both dimensions
      .h5ad_copy_sparse_pairwise(src_file, dest_file, child_path, idx,
                                 compression)
    } else if (child_type == "H5I_DATASET") {
      mat <- rhdf5::h5read(src_file, child_path)
      if (is.matrix(mat)) {
        rhdf5::h5write(mat[idx, idx, drop = FALSE], dest_file, child_path)
      } else {
        rhdf5::h5write(mat, dest_file, child_path)
      }
      .h5ad_copy_attrs(src_file, dest_file, child_path)
    } else {
      .h5ad_copy_group_recursive(src_file, dest_file, child_path)
    }
  }

  .h5ad_copy_attrs(src_file, dest_file, group_path)
}

#' Subset a sparse pairwise matrix (obsp/varp child)
#'
#' Reads the sparse matrix, subsets both rows and columns, writes back.
#' These matrices are typically small enough to materialize.
#'
#' @keywords internal
#' @noRd
.h5ad_copy_sparse_pairwise <- function(src_file, dest_file, path, idx,
                                       compression) {
  enc <- .h5ad_detect_encoding(src_file, path)
  attrs <- rhdf5::h5readAttributes(src_file, path)
  shape <- as.integer(attrs[["shape"]])

  data    <- rhdf5::h5read(src_file, paste0(path, "/data"))
  indices <- rhdf5::h5read(src_file, paste0(path, "/indices"))
  indptr  <- rhdf5::h5read(src_file, paste0(path, "/indptr"))

  # Reconstruct as dgCMatrix/dgRMatrix, subset, write back
  n <- shape[1]; m <- shape[2]
  if (enc == "csr_matrix") {
    # Convert to triplet form
    rows <- rep(seq_len(n), diff(indptr)) # 1-based
    cols <- indices + 1L                   # 0-based to 1-based
    sp <- Matrix::sparseMatrix(i = rows, j = cols, x = as.numeric(data),
                               dims = c(n, m))
    sub <- sp[idx, idx, drop = FALSE]
  } else {
    cols <- rep(seq_len(m), diff(indptr))
    rows <- indices + 1L
    sp <- Matrix::sparseMatrix(i = rows, j = cols, x = as.numeric(data),
                               dims = c(n, m))
    sub <- sp[idx, idx, drop = FALSE]
  }

  # Write back as CSR
  sub_csr <- as(sub, "RsparseMatrix")
  rhdf5::h5createGroup(dest_file, path)
  rhdf5::h5write(sub_csr@x, dest_file, paste0(path, "/data"))
  rhdf5::h5write(sub_csr@j, dest_file, paste0(path, "/indices"))
  rhdf5::h5write(sub_csr@p, dest_file, paste0(path, "/indptr"))

  n_out <- length(idx)
  rhdf5::h5writeAttribute(c(n_out, n_out), dest_file, "shape",
                          h5loc = path)
  rhdf5::h5writeAttribute(enc, dest_file, "encoding-type",
                          h5loc = path, asScalar = TRUE)
  enc_ver <- attrs[["encoding-version"]]
  if (!is.null(enc_ver)) {
    rhdf5::h5writeAttribute(enc_ver, dest_file, "encoding-version",
                            h5loc = path, asScalar = TRUE)
  }
}

#' Merge contiguous or overlapping 0-based [start, end) ranges
#'
#' @param starts Integer vector of 0-based start positions.
#' @param ends Integer vector of 0-based end positions (exclusive).
#' @return List of 2-element integer vectors \code{c(start, end)}.
#' @keywords internal
#' @noRd
.merge_h5_ranges <- function(starts, ends) {
  ord <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]
  merged <- list()
  cs <- starts[1L]; ce <- ends[1L]
  for (i in seq_along(starts)[-1L]) {
    if (starts[i] <= ce) {
      ce <- max(ce, ends[i])
    } else {
      merged[[length(merged) + 1L]] <- c(cs, ce)
      cs <- starts[i]; ce <- ends[i]
    }
  }
  merged[[length(merged) + 1L]] <- c(cs, ce)
  merged
}

#' Chunked copy of a matrix group, dispatching on encoding-type
#'
#' @param src_file Character scalar.
#' @param dest_file Character scalar.
#' @param src_path Character scalar. HDF5 path of the source matrix.
#' @param dest_path Character scalar. HDF5 path for the output matrix.
#' @param cell_idx Integer vector. 1-based row indices.
#' @param gene_idx Integer vector or NULL. 1-based column indices.
#' @param n_obs_src Integer. Number of source rows (cells).
#' @param n_var_src Integer. Number of source columns (genes).
#' @param chunk_size Integer.
#' @param compression Integer.
#' @param verbose Logical.
#' @keywords internal
#' @noRd
.h5ad_copy_matrix_chunked <- function(src_file, dest_file,
                                      src_path, dest_path,
                                      cell_idx, gene_idx,
                                      n_obs_src, n_var_src,
                                      chunk_size, compression,
                                      verbose) {
  enc <- .h5ad_detect_encoding(src_file, src_path)
  if (is.na(enc)) enc <- "array"

  n_obs_out <- length(cell_idx)
  n_var_out <- if (is.null(gene_idx)) n_var_src else length(gene_idx)

  if (enc == "array") {
    .h5ad_copy_dense_chunked(src_file, dest_file, src_path, dest_path,
                             cell_idx, gene_idx,
                             n_obs_out, n_var_out,
                             chunk_size, compression, verbose)
  } else if (enc == "csr_matrix") {
    .h5ad_copy_csr_chunked(src_file, dest_file, src_path, dest_path,
                           cell_idx, gene_idx,
                           n_obs_src, n_var_src,
                           n_obs_out, n_var_out,
                           chunk_size, compression, verbose)
  } else if (enc == "csc_matrix") {
    .h5ad_copy_csc_chunked(src_file, dest_file, src_path, dest_path,
                           cell_idx, gene_idx,
                           n_obs_src, n_var_src,
                           n_obs_out, n_var_out,
                           chunk_size, compression, verbose)
  } else {
    stop("Unsupported matrix encoding: ", enc,
         " (at ", src_path, ")")
  }
}

#' Chunked copy of a dense matrix
#' @keywords internal
#' @noRd
.h5ad_copy_dense_chunked <- function(src_file, dest_file,
                                     src_path, dest_path,
                                     cell_idx, gene_idx,
                                     n_obs_out, n_var_out,
                                     chunk_size, compression,
                                     verbose) {
  # Pre-create the output dataset
  rhdf5::h5createDataset(dest_file, dest_path,
                         dims = c(n_obs_out, n_var_out),
                         storage.mode = "double",
                         chunk = c(min(chunk_size, n_obs_out), n_var_out),
                         level = compression)

  col_sel <- if (is.null(gene_idx)) seq_len(n_var_out) else gene_idx

  for (chunk_start in seq(1L, n_obs_out, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1L, n_obs_out)
    src_rows <- cell_idx[chunk_start:chunk_end]

    chunk <- rhdf5::h5read(src_file, src_path,
                           index = list(src_rows, col_sel))
    rhdf5::h5write(chunk, dest_file, dest_path,
                   index = list(chunk_start:chunk_end, seq_len(n_var_out)))

    if (verbose) {
      cat(sprintf("  dense: rows %d-%d / %d\n",
                  chunk_start, chunk_end, n_obs_out))
    }
  }

  # Write encoding attributes
  rhdf5::h5writeAttribute("array", dest_file, "encoding-type",
                          h5loc = dest_path, asScalar = TRUE)
  src_attrs <- rhdf5::h5readAttributes(src_file, src_path)
  enc_ver <- src_attrs[["encoding-version"]]
  if (!is.null(enc_ver)) {
    rhdf5::h5writeAttribute(enc_ver, dest_file, "encoding-version",
                            h5loc = dest_path, asScalar = TRUE)
  }
}

#' Chunked copy of a CSR sparse matrix
#' @keywords internal
#' @noRd
.h5ad_copy_csr_chunked <- function(src_file, dest_file,
                                   src_path, dest_path,
                                   cell_idx, gene_idx,
                                   n_obs_src, n_var_src,
                                   n_obs_out, n_var_out,
                                   chunk_size, compression,
                                   verbose) {
  # Read full indptr (small: n_obs_src + 1 integers)
  indptr_src <- rhdf5::h5read(src_file, paste0(src_path, "/indptr"))
  # indptr_src values are 0-based; indptr_src[k] (R 1-based k) gives the
  # 0-based start position for source row k.

  # Gene remapping: 0-based source col -> 0-based dest col
  gene_remap <- NULL
  if (!is.null(gene_idx)) {
    gene_remap <- rep(NA_integer_, n_var_src)
    gene_remap[gene_idx] <- seq_along(gene_idx) - 1L
  }

  out_data    <- list()
  out_indices <- list()
  out_indptr  <- integer(n_obs_out + 1L)
  out_indptr[1L] <- 0L
  cumulative_nnz <- 0L

  for (chunk_start in seq(1L, n_obs_out, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1L, n_obs_out)
    src_rows  <- cell_idx[chunk_start:chunk_end]

    # 0-based start/end for each source row from indptr
    # R-indexed: source row r has indptr_src[r] .. indptr_src[r+1]-1
    row_starts <- indptr_src[src_rows]
    row_ends   <- indptr_src[src_rows + 1L]
    row_nnz    <- row_ends - row_starts
    nonempty   <- which(row_nnz > 0L)

    if (length(nonempty) > 0L) {
      # Merge read ranges for efficiency
      read_ranges <- .merge_h5_ranges(row_starts[nonempty], row_ends[nonempty])

      # Read all needed data and indices
      all_data    <- vector("list", length(read_ranges))
      all_indices <- vector("list", length(read_ranges))
      # Build a lookup: for each merged range, what 0-based positions it covers
      range_starts <- vapply(read_ranges, `[`, integer(1), 1L)
      range_ends   <- vapply(read_ranges, `[`, integer(1), 2L)

      for (r in seq_along(read_ranges)) {
        rng_start_0 <- read_ranges[[r]][1]
        rng_end_0   <- read_ranges[[r]][2]
        rng_1based  <- (rng_start_0 + 1L):rng_end_0
        all_data[[r]]    <- rhdf5::h5read(src_file,
                                          paste0(src_path, "/data"),
                                          index = list(rng_1based))
        all_indices[[r]] <- rhdf5::h5read(src_file,
                                          paste0(src_path, "/indices"),
                                          index = list(rng_1based))
      }

      # Concatenate all read data/indices with a global offset map
      # For each merged range, compute the cumulative offset in
      # the concatenated buffer
      buf_data    <- unlist(all_data)
      buf_indices <- unlist(all_indices)
      buf_offsets <- c(0L, cumsum(vapply(all_data, length, integer(1))))
      # Position p (0-based in source file) maps to buffer position:
      #   find which merged range contains it, then offset within that range

      # Process each row in this chunk
      for (k in seq_along(src_rows)) {
        if (row_nnz[k] == 0L) {
          out_indptr[chunk_start + k] <- cumulative_nnz
          next
        }
        rs <- row_starts[k]  # 0-based start in source
        re <- row_ends[k]    # 0-based end (exclusive)

        # Find which merged range contains this row's data
        # rs >= range_starts & rs < range_ends
        ri <- which(range_starts <= rs & range_ends >= re)
        if (length(ri) != 1L) ri <- ri[1L]  # should always be exactly 1

        buf_start <- buf_offsets[ri] + (rs - range_starts[ri]) + 1L
        buf_end   <- buf_start + (re - rs) - 1L

        row_data    <- buf_data[buf_start:buf_end]
        row_indices <- buf_indices[buf_start:buf_end]

        if (!is.null(gene_remap)) {
          keep <- !is.na(gene_remap[row_indices + 1L])
          row_data    <- row_data[keep]
          row_indices <- gene_remap[row_indices[keep] + 1L]
        }

        if (length(row_data) > 0L) {
          out_data[[length(out_data) + 1L]]       <- row_data
          out_indices[[length(out_indices) + 1L]]  <- row_indices
        }
        cumulative_nnz <- cumulative_nnz + length(row_data)
        out_indptr[chunk_start + k] <- cumulative_nnz
      }
    } else {
      # All rows in chunk are empty
      seq_range <- (chunk_start + 1L):(chunk_end + 1L)
      out_indptr[seq_range] <- cumulative_nnz
    }

    if (verbose) {
      cat(sprintf("  csr: rows %d-%d / %d (nnz so far: %d)\n",
                  chunk_start, chunk_end, n_obs_out, cumulative_nnz))
    }
  }

  # Concatenate accumulators
  out_data_vec    <- if (length(out_data) > 0L) unlist(out_data) else numeric(0)
  out_indices_vec <- if (length(out_indices) > 0L) unlist(out_indices) else integer(0)

  # Ensure correct types
  out_indices_vec <- as.integer(out_indices_vec)
  out_indptr      <- as.integer(out_indptr)

  # Write output CSR
  rhdf5::h5createGroup(dest_file, dest_path)
  rhdf5::h5write(out_data_vec,    dest_file, paste0(dest_path, "/data"))
  rhdf5::h5write(out_indices_vec, dest_file, paste0(dest_path, "/indices"))
  rhdf5::h5write(out_indptr,      dest_file, paste0(dest_path, "/indptr"))

  # Attributes
  src_attrs <- rhdf5::h5readAttributes(src_file, src_path)
  rhdf5::h5writeAttribute(c(n_obs_out, n_var_out), dest_file, "shape",
                          h5loc = dest_path)
  rhdf5::h5writeAttribute("csr_matrix", dest_file, "encoding-type",
                          h5loc = dest_path, asScalar = TRUE)
  enc_ver <- src_attrs[["encoding-version"]]
  if (!is.null(enc_ver)) {
    rhdf5::h5writeAttribute(enc_ver, dest_file, "encoding-version",
                            h5loc = dest_path, asScalar = TRUE)
  }
}

#' Chunked copy of a CSC sparse matrix
#' @keywords internal
#' @noRd
.h5ad_copy_csc_chunked <- function(src_file, dest_file,
                                   src_path, dest_path,
                                   cell_idx, gene_idx,
                                   n_obs_src, n_var_src,
                                   n_obs_out, n_var_out,
                                   chunk_size, compression,
                                   verbose) {
  # Read full indptr (small: n_var_src + 1 integers)
  indptr_src <- rhdf5::h5read(src_file, paste0(src_path, "/indptr"))

  # Cell remapping: 0-based source row -> 0-based dest row
  cell_remap <- rep(NA_integer_, n_obs_src)
  cell_remap[cell_idx] <- seq_along(cell_idx) - 1L

  # Columns to iterate over
  col_order <- if (is.null(gene_idx)) seq_len(n_var_src) else gene_idx

  out_data    <- list()
  out_indices <- list()
  out_indptr  <- integer(n_var_out + 1L)
  out_indptr[1L] <- 0L
  cumulative_nnz <- 0L

  for (chunk_start in seq(1L, n_var_out, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1L, n_var_out)
    src_cols  <- col_order[chunk_start:chunk_end]

    # 0-based start/end for each source column
    col_starts <- indptr_src[src_cols]
    col_ends   <- indptr_src[src_cols + 1L]
    col_nnz    <- col_ends - col_starts
    nonempty   <- which(col_nnz > 0L)

    if (length(nonempty) > 0L) {
      read_ranges <- .merge_h5_ranges(col_starts[nonempty], col_ends[nonempty])
      range_starts <- vapply(read_ranges, `[`, integer(1), 1L)
      range_ends   <- vapply(read_ranges, `[`, integer(1), 2L)

      all_data    <- vector("list", length(read_ranges))
      all_indices <- vector("list", length(read_ranges))
      for (r in seq_along(read_ranges)) {
        rng_1based <- (read_ranges[[r]][1] + 1L):read_ranges[[r]][2]
        all_data[[r]]    <- rhdf5::h5read(src_file,
                                          paste0(src_path, "/data"),
                                          index = list(rng_1based))
        all_indices[[r]] <- rhdf5::h5read(src_file,
                                          paste0(src_path, "/indices"),
                                          index = list(rng_1based))
      }
      buf_data    <- unlist(all_data)
      buf_indices <- unlist(all_indices)
      buf_offsets <- c(0L, cumsum(vapply(all_data, length, integer(1))))

      for (k in seq_along(src_cols)) {
        if (col_nnz[k] == 0L) {
          out_indptr[chunk_start + k] <- cumulative_nnz
          next
        }
        cs <- col_starts[k]
        ce <- col_ends[k]
        ri <- which(range_starts <= cs & range_ends >= ce)
        if (length(ri) != 1L) ri <- ri[1L]

        buf_start <- buf_offsets[ri] + (cs - range_starts[ri]) + 1L
        buf_end   <- buf_start + (ce - cs) - 1L

        col_data    <- buf_data[buf_start:buf_end]
        col_indices <- buf_indices[buf_start:buf_end]

        # Filter to only rows in cell_idx, remap row indices
        keep <- !is.na(cell_remap[col_indices + 1L])
        col_data    <- col_data[keep]
        col_indices <- cell_remap[col_indices[keep] + 1L]

        if (length(col_data) > 0L) {
          # Sort by new row index for valid CSC
          ord <- order(col_indices)
          out_data[[length(out_data) + 1L]]      <- col_data[ord]
          out_indices[[length(out_indices) + 1L]] <- col_indices[ord]
        }
        cumulative_nnz <- cumulative_nnz + length(col_data)
        out_indptr[chunk_start + k] <- cumulative_nnz
      }
    } else {
      seq_range <- (chunk_start + 1L):(chunk_end + 1L)
      out_indptr[seq_range] <- cumulative_nnz
    }

    if (verbose) {
      cat(sprintf("  csc: cols %d-%d / %d (nnz so far: %d)\n",
                  chunk_start, chunk_end, n_var_out, cumulative_nnz))
    }
  }

  out_data_vec    <- if (length(out_data) > 0L) unlist(out_data) else numeric(0)
  out_indices_vec <- if (length(out_indices) > 0L) unlist(out_indices) else integer(0)
  out_indices_vec <- as.integer(out_indices_vec)
  out_indptr      <- as.integer(out_indptr)

  rhdf5::h5createGroup(dest_file, dest_path)
  rhdf5::h5write(out_data_vec,    dest_file, paste0(dest_path, "/data"))
  rhdf5::h5write(out_indices_vec, dest_file, paste0(dest_path, "/indices"))
  rhdf5::h5write(out_indptr,      dest_file, paste0(dest_path, "/indptr"))

  src_attrs <- rhdf5::h5readAttributes(src_file, src_path)
  rhdf5::h5writeAttribute(c(n_obs_out, n_var_out), dest_file, "shape",
                          h5loc = dest_path)
  rhdf5::h5writeAttribute("csc_matrix", dest_file, "encoding-type",
                          h5loc = dest_path, asScalar = TRUE)
  enc_ver <- src_attrs[["encoding-version"]]
  if (!is.null(enc_ver)) {
    rhdf5::h5writeAttribute(enc_ver, dest_file, "encoding-version",
                            h5loc = dest_path, asScalar = TRUE)
  }
}

#' Subset and reorder an on-disk h5ad file without full materialization
#'
#' Performs a chunked, streaming copy of an HDF5-backed .h5ad file,
#' subsetting and reordering cells (rows of X) and optionally genes
#' (columns of X) without ever loading the full matrix into memory.
#'
#' Requires the \pkg{rhdf5} Bioconductor package.
#'
#' @param .source.path Character scalar. Path to the source .h5ad file.
#' @param .dest.path Character scalar. Path for the output .h5ad file.
#'   Must differ from \code{.source.path}.
#' @param .cell.idx Integer vector. 1-based row (cell) indices to keep,
#'   in the desired output order. Must be valid indices into the source obs.
#' @param .gene.idx Integer vector or \code{NULL}. 1-based column (gene)
#'   indices to keep, in the desired output order. \code{NULL} keeps all genes
#'   in original order.
#' @param .layers Character vector or \code{NULL}. Names of additional layers
#'   (under \code{/layers/}) to copy. \code{NULL} copies all layers found.
#'   Use \code{character(0)} to skip layers entirely.
#' @param .copy.obsm Logical. Whether to copy \code{/obsm} (cell embeddings).
#'   Default \code{TRUE}.
#' @param .copy.obsp Logical. Whether to copy \code{/obsp} (cell-cell graphs).
#'   Default \code{FALSE} (these are large and rarely needed after subsetting).
#' @param .copy.varm Logical. Whether to copy \code{/varm} (gene loadings).
#'   Default \code{TRUE}.
#' @param .copy.varp Logical. Whether to copy \code{/varp}. Default \code{FALSE}.
#' @param .copy.uns  Logical. Whether to copy \code{/uns} (unstructured metadata).
#'   Default \code{TRUE}.
#' @param .chunk.size Integer. Number of cells (rows) per chunk when copying
#'   matrices. Controls peak memory. Default \code{5000L}.
#' @param .compression Integer 0-9. gzip compression level for output datasets.
#'   \code{0L} = no compression (fastest). Default \code{0L}.
#' @param .verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisible character scalar: the \code{.dest.path}.
#'
#' @details
#' \strong{Supported matrix encodings}: \code{"csr_matrix"}, \code{"csc_matrix"},
#' and \code{"array"} (dense). The encoding is detected from the
#' \code{encoding-type} HDF5 attribute on each matrix group/dataset.
#'
#' \strong{Peak memory}: approximately \code{.chunk.size * n_genes_out} numeric
#' values, regardless of total dataset size.
#'
#' \strong{Categorical columns}: obs/var columns stored as h5ad categoricals
#' (codes + categories) are correctly handled.
#'
#' @seealso \code{\link[anndataR]{read_h5ad}} to open the result.
#'
#' @examples
#' \dontrun{
#' subset_h5ad(
#'   .source.path = "full_data.h5ad",
#'   .dest.path   = "subset.h5ad",
#'   .cell.idx    = which(metadata$keep),
#'   .gene.idx    = which(!duplicated(gene_symbols)),
#'   .chunk.size  = 5000L
#' )
#' h5ad <- anndataR::read_h5ad("subset.h5ad", as = "HDF5AnnData")
#' }
#'
#' @export
subset_h5ad <- function(.source.path,
                        .dest.path,
                        .cell.idx,
                        .gene.idx    = NULL,
                        .layers      = NULL,
                        .copy.obsm   = TRUE,
                        .copy.obsp   = FALSE,
                        .copy.varm   = TRUE,
                        .copy.varp   = FALSE,
                        .copy.uns    = TRUE,
                        .chunk.size  = 5000L,
                        .compression = 0L,
                        .verbose     = TRUE) {

  # ── dependency check ──
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required for subset_h5ad(). ",
         "Install it from Bioconductor: ",
         "BiocManager::install(\"rhdf5\")",
         call. = FALSE)
  }

  # ── input validation ──
  stopifnot(
    is.character(.source.path), length(.source.path) == 1L,
    is.character(.dest.path),   length(.dest.path) == 1L,
    file.exists(.source.path),
    normalizePath(.source.path, mustWork = TRUE) !=
      normalizePath(.dest.path, mustWork = FALSE),
    is.integer(.cell.idx) || is.numeric(.cell.idx),
    length(.cell.idx) > 0L,
    is.null(.gene.idx) || is.integer(.gene.idx) || is.numeric(.gene.idx),
    is.logical(.copy.obsm),  is.logical(.copy.obsp),
    is.logical(.copy.varm),  is.logical(.copy.varp),
    is.logical(.copy.uns),   is.logical(.verbose)
  )

  .cell.idx  <- as.integer(.cell.idx)
  if (!is.null(.gene.idx)) .gene.idx <- as.integer(.gene.idx)
  .chunk.size  <- as.integer(.chunk.size)
  .compression <- as.integer(.compression)

  # ── determine source dimensions from /obs and /var ──
  src_file <- .source.path

  obs_attrs <- rhdf5::h5readAttributes(src_file, "/obs")
  var_attrs <- rhdf5::h5readAttributes(src_file, "/var")
  obs_index_name <- obs_attrs[["_index"]]
  var_index_name <- var_attrs[["_index"]]

  obs_index <- rhdf5::h5read(src_file, paste0("/obs/", obs_index_name))
  var_index <- rhdf5::h5read(src_file, paste0("/var/", var_index_name))
  n_obs_src <- length(obs_index)
  n_var_src <- length(var_index)

  stopifnot(
    all(.cell.idx >= 1L & .cell.idx <= n_obs_src),
    is.null(.gene.idx) || all(.gene.idx >= 1L & .gene.idx <= n_var_src)
  )

  n_obs_out <- length(.cell.idx)
  n_var_out <- if (is.null(.gene.idx)) n_var_src else length(.gene.idx)

  # ── create destination file ──
  if (file.exists(.dest.path)) file.remove(.dest.path)
  rhdf5::h5createFile(.dest.path)

  # Write root attributes
  rhdf5::h5writeAttribute("anndata", .dest.path, "encoding-type",
                          h5loc = "/", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.1.0", .dest.path, "encoding-version",
                          h5loc = "/", asScalar = TRUE)

  # ── 1. Copy /obs ──
  if (.verbose) cat("Copying /obs ...\n")
  .h5ad_copy_dataframe_group(src_file, .dest.path, "/obs",
                             .cell.idx, .compression)

  # ── 2. Copy /var ──
  if (.verbose) cat("Copying /var ...\n")
  var_idx <- if (is.null(.gene.idx)) seq_len(n_var_src) else .gene.idx
  .h5ad_copy_dataframe_group(src_file, .dest.path, "/var",
                             var_idx, .compression)

  # ── 3. Copy /X ──
  if (.verbose) cat("Copying /X ...\n")
  .h5ad_copy_matrix_chunked(src_file, .dest.path,
                            "/X", "/X",
                            .cell.idx, .gene.idx,
                            n_obs_src, n_var_src,
                            .chunk.size, .compression, .verbose)

  # ── 4. Copy /layers ──
  fid <- rhdf5::H5Fopen(src_file)
  has_layers <- rhdf5::H5Lexists(fid, "/layers")
  rhdf5::H5Fclose(fid)

  if (has_layers && !identical(.layers, character(0))) {
    ls_all <- rhdf5::h5ls(src_file, recursive = TRUE)
    layer_children <- ls_all[ls_all$group == "/layers", , drop = FALSE]
    layer_names <- layer_children$name

    if (!is.null(.layers)) {
      layer_names <- intersect(layer_names, .layers)
    }

    if (length(layer_names) > 0L) {
      rhdf5::h5createGroup(.dest.path, "/layers")
      for (ln in layer_names) {
        if (.verbose) cat(sprintf("Copying /layers/%s ...\n", ln))
        .h5ad_copy_matrix_chunked(src_file, .dest.path,
                                  paste0("/layers/", ln),
                                  paste0("/layers/", ln),
                                  .cell.idx, .gene.idx,
                                  n_obs_src, n_var_src,
                                  .chunk.size, .compression, .verbose)
      }
    }
  }

  # ── 5. Copy /obsm ──
  if (.copy.obsm) {
    if (.verbose) cat("Copying /obsm ...\n")
    .h5ad_copy_embedding_group(src_file, .dest.path, "/obsm",
                               .cell.idx, .compression)
  }

  # ── 6. Copy /varm ──
  if (.copy.varm) {
    if (.verbose) cat("Copying /varm ...\n")
    .h5ad_copy_embedding_group(src_file, .dest.path, "/varm",
                               var_idx, .compression)
  }

  # ── 7. Copy /obsp ──
  if (.copy.obsp) {
    if (.verbose) cat("Copying /obsp ...\n")
    .h5ad_copy_pairwise_group(src_file, .dest.path, "/obsp",
                              .cell.idx, .compression)
  }

  # ── 8. Copy /varp ──
  if (.copy.varp) {
    if (.verbose) cat("Copying /varp ...\n")
    .h5ad_copy_pairwise_group(src_file, .dest.path, "/varp",
                              var_idx, .compression)
  }

  # ── 9. Copy /uns ──
  if (.copy.uns) {
    fid <- rhdf5::H5Fopen(src_file)
    has_uns <- rhdf5::H5Lexists(fid, "/uns")
    rhdf5::H5Fclose(fid)
    if (has_uns) {
      if (.verbose) cat("Copying /uns ...\n")
      .h5ad_copy_group_recursive(src_file, .dest.path, "/uns")
    }
  }

  if (.verbose) cat("Done. Output: ", .dest.path, "\n")
  invisible(.dest.path)
}

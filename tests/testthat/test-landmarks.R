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

# Test for setup.tdr.obj

test_that("setup.tdr.obj returns a list with correct names", {
  .cells <- list(sample1 = matrix(data = runif(n = 30), 
                                  nrow = 10, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample1",
                                                         formatC(x = 1:10,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8"))),
                 sample2 = matrix(data = runif(n = 30), 
                                  nrow = 10, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample2",
                                                         formatC(x = 1:10,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8")))) |>
    lapply(FUN = function(x){
      uri <- tempfile(fileext = ".RDS")
      
      saveRDS(object = x,
              file = uri,
              compress = FALSE)
      
      return(uri)
    })
  .meta <- data.frame(row.names = c("sample1", "sample2"),
                      group = c("A", "B"))
  result <- setup.tdr.obj(.cells = .cells,
                         .meta = .meta,
                         .markers = c("CD3", "CD4", "CD8"),
                         .assay.type = "cyto",
                         .verbose = FALSE)
  expect_true(is.TDRObj(result))
  expect_true(all(c("cells",
                    "metadata",
                    "config",
                    "integration",
                    "assay",
                    "landmark.embed",
                    "landmark.annot",
                    "graphs",
                    "density",
                    "sample.embed",
                    "cellmap",
                    "results") %in% names(x = result)))
  # Check nested config structure
  expect_true(all(c("key", "sampling", "assay.type", "markers", "n.threads") %in% names(x = result$config)))
  # Check nested integration structure
  expect_true(all(c("harmony.var", "harmony.obj") %in% names(x = result$integration)))
})

# Test error conditions for setup.tdr.obj

test_that("setup.tdr.obj throws error when .cells names don't match .meta rownames", {
  .cells <- list(wrong_name = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.tdr.obj(.cells = .cells, .meta = .meta),
               "Sample names mismatch between .cells and .meta")
})

test_that("setup.tdr.obj throws error when .harmony.var not in metadata", {
  .cells <- list(sample1 = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.tdr.obj(.cells = .cells, .meta = .meta, .harmony.var = "missing_var"),
               "Variables not found in metadata")
})

test_that("setup.tdr.obj throws error when .markers insufficient for cyto", {
  .cells <- list(sample1 = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.tdr.obj(.cells = .cells, .meta = .meta, .markers = c("CD3"), .assay.type = "cyto"),
               ".markers must contain at least 3 markers for meaningful dimensionality reduction")
})

test_that("setup.tdr.obj throws error when .markers used with RNA", {
  .cells <- list(sample1 = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.tdr.obj(.cells = .cells, .meta = .meta, .markers = c("CD3"), .assay.type = "RNA"),
               ".markers argument only applies to cytometry data")
})

# Test for get.landmarks

test_that("get.landmarks validates input properly", {
  # Test that get.landmarks expects a proper .tdr.obj structure
  expect_error(get.landmarks(list()),
               "no applicable method")  # S3 dispatch error for non-TDRObj
  
  # Test with incomplete object
  incomplete_obj <- list(
    config = list(assay.type = "RNA"),
    cells = list()
  )
  expect_error(get.landmarks(incomplete_obj),
               "no applicable method")  # S3 dispatch error for non-TDRObj
})

# Test for HVG exclusion pattern

test_that("HVG exclusion pattern excludes VDJ variable-region genes", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"

  # TCR alpha/beta V/D/J segments
  tcr_ab <- c("TRAV1-1", "TRAV12-1", "TRAV38-2DV8", "TRBV2", "TRBV28",
              "TRAJ1", "TRAJ61", "TRBD1", "TRBD2", "TRBJ1-1", "TRBJ2-7")
  expect_true(all(grepl(pattern, tcr_ab, ignore.case = TRUE)))


  # TCR gamma/delta V/D/J segments
  tcr_gd <- c("TRDV1", "TRDV2", "TRDV3", "TRDD1", "TRDD2", "TRDJ1",
              "TRGV1", "TRGV2", "TRGV9", "TRGV10", "TRGJ1", "TRGJ2")
  expect_true(all(grepl(pattern, tcr_gd, ignore.case = TRUE)))

  # Ig V/D/J segments
  ig_vdj <- c("IGHV1-2", "IGHV3-23", "IGHV4-34", "IGHD1-1", "IGHD3-10",
              "IGHJ1", "IGHJ6", "IGKV1-5", "IGKV3-20", "IGKJ1", "IGKJ5",
              "IGLV1-40", "IGLV2-14", "IGLJ1", "IGLJ7")
  expect_true(all(grepl(pattern, ig_vdj, ignore.case = TRUE)))
})

test_that("HVG exclusion pattern retains constant-region genes", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"

  # TCR constant regions
  tcr_const <- c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2")
  expect_false(any(grepl(pattern, tcr_const, ignore.case = TRUE)))

  # Ig constant regions and isotype genes
  ig_const <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",
                "IGHM", "IGHD", "IGHE", "IGKC",
                "IGLC1", "IGLC2", "IGLC3", "IGLC7",
                "IGLL1", "IGLL5", "JCHAIN")
  expect_false(any(grepl(pattern, ig_const, ignore.case = TRUE)))
})

test_that("HVG exclusion pattern excludes mitochondrial genes", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"
  mito <- c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND4L",
            "MT-ATP6", "MT-ATP8", "MT-CYB", "MT-RNR1", "MT-RNR2", "MT-TA")
  expect_true(all(grepl(pattern, mito, ignore.case = TRUE)))

  # Nuclear-encoded MT pseudogenes should NOT match
  expect_false(any(grepl(pattern, c("MTRNR2L1", "MTRNR2L8", "MTRNR2L12"), ignore.case = TRUE)))
})

test_that("HVG exclusion pattern excludes ribosomal protein genes", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"

  rps <- c("RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS6",
           "RPS15A", "RPS27A", "RPS28", "RPS29")
  expect_true(all(grepl(pattern, rps, ignore.case = TRUE)))

  rpl <- c("RPL3", "RPL7A", "RPL7L1", "RPL10A", "RPL22L1", "RPL26L1",
           "RPL35A", "RPL36A", "RPL36AL", "RPL39L", "RPL41")
  expect_true(all(grepl(pattern, rpl, ignore.case = TRUE)))

  rplp_rpsa <- c("RPLP0", "RPLP1", "RPLP2", "RPSA")
  expect_true(all(grepl(pattern, rplp_rpsa, ignore.case = TRUE)))
})

test_that("HVG exclusion pattern does not match non-ribosomal RP-prefixed genes", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"

  non_ribo <- c("RPA1", "RPA2", "RPA3", "RPE", "RPE65",
                "RPGR", "RPGRIP1", "RPH3A", "RPAIN", "RPF1",
                "RPP14", "RPP25", "RPP30", "RPP40",
                "RPRM", "RPRD1A", "RPRD2", "RPTOR",
                "RPUSD1", "RPUSD2", "RPN1", "RPN2",
                "RPS6KA1", "RPS6KB1", "RPS6KC1")
  expect_false(any(grepl(pattern, non_ribo, ignore.case = TRUE)))
})

test_that("HVG exclusion pattern does not match housekeeping/marker genes", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"
  safe <- c("GAPDH", "ACTB", "B2M", "CD3D", "CD4", "CD8A", "FOXP3",
            "HLA-A", "HLA-B", "NKG7", "GZMB", "PTPRC",
            "TRAP1", "TRAF1", "TRAM1", "IGF1", "IGFBP1", "IGFBP3")
  expect_false(any(grepl(pattern, safe, ignore.case = TRUE)))
})

test_that("HVG exclusion pattern rejects known false positives", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"

  # RPS6KL1 (kinase, not ribosomal protein), ribosomal pseudogenes
  expect_false(any(grepl(pattern, c("RPS6KL1", "RPL3P4", "RPS10P2"), ignore.case = TRUE)))
})

test_that("HVG exclusion pattern covers HGNC edge-case ribosomal symbols", {
  pattern <- "^TR[ABDG][VDJ]\\d|^IG[KHL][VDJ]\\d|^MT-|^RP[SL]\\d{1,2}[AXYL]?L?\\d?$|^RPLP[012]$|^RPSA$"

  # "like" genes, sex-linked variants, and other HGNC-approved edge cases
  hgnc_edge <- c("RPL3L", "RPL10L", "RPL36AL", "RPS4Y2", "RPS27L")
  expect_true(all(grepl(pattern, hgnc_edge, ignore.case = TRUE)))
})

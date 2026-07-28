# Per-library doublet detection for the mouse-brain snRNA-seq workflow.
# This file defines helper functions and is sourced by
# snRNAseq_mouseBrain_integration.R before the two libraries are merged.

suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
})

scdblfinder_randomized_pca <- function(e, dims = 20, ...) {
  sce <- SingleCellExperiment(
    assays = list(counts = e)
  )
  sce <- scuttle::logNormCounts(sce)
  pca <- scater::calculatePCA(
    SummarizedExperiment::assay(sce, "logcounts"),
    ncomponents = dims,
    subset_row = seq_len(nrow(e)),
    ntop = nrow(e),
    BSPARAM = BiocSingular::RandomParam(deferred = TRUE)
  )
  if (is.list(pca)) {
    pca <- pca$x
  }
  rownames(pca) <- colnames(e)
  pca
}

run_scdblfinder <- function(counts, library_id, seed = 20260728) {
  # Expected 10x multiplet rate: approximately 1% per 1,000 recovered nuclei.
  expected_rate <- 0.01 * ncol(counts) / 1000

  sce <- SingleCellExperiment(
    assays = list(counts = counts)
  )
  set.seed(seed)
  sce <- scDblFinder(
    sce,
    clusters = NULL,
    dbr = expected_rate,
    processing = scdblfinder_randomized_pca,
    verbose = TRUE
  )

  calls <- data.frame(
    library = library_id,
    barcode = colnames(sce),
    scDblFinder_score = sce$scDblFinder.score,
    scDblFinder_class = as.character(sce$scDblFinder.class),
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    library = library_id,
    input_nuclei = ncol(sce),
    expected_rate = expected_rate,
    expected_doublets = round(ncol(sce) * expected_rate),
    flagged_doublets = sum(calls$scDblFinder_class == "doublet"),
    flagged_rate = mean(calls$scDblFinder_class == "doublet"),
    removed_doublets = sum(calls$scDblFinder_class == "doublet"),
    retained_singlets = sum(calls$scDblFinder_class == "singlet"),
    stringsAsFactors = FALSE
  )

  list(
    singlet_barcodes =
      calls$barcode[calls$scDblFinder_class == "singlet"],
    calls = calls,
    summary = summary
  )
}

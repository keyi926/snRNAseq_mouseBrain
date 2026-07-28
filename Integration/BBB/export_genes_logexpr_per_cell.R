# Export per-cell log-normalized expression + metadata for the genes shown in
# VlnPlot_endothelial_key_BBB_split.svg and VlnPlot_pericyte_key_BBB_split.svg
# (same gene sets / cell subsets / assay+slot as snRNAseq_mouseBrain_BBB_integrity.R section 5)

suppressMessages({ library(Seurat); library(dplyr); library(tidyr) })

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(file.path(getwd(), "BBB"))
}
root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(file.path(root, "Integration"))
outdir <- "BBB"

mb <- readRDS("mousebrain_integrated_final.rds")
DefaultAssay(mb) <- "RNA"
mb$sample <- factor(mb$sample, levels = c("Sample1", "Sample2"))

meta_cols <- c("sample", "celltype", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt")

extract_long <- function(obj, genes, panel_name) {
  obj <- subset(obj, features = genes)
  dat <- GetAssayData(obj, assay = "RNA", layer = "data")  # log-normalized expression
  dat <- as.data.frame(t(as.matrix(dat)))
  dat$cell <- rownames(dat)
  meta <- obj@meta.data[, meta_cols]
  meta$cell <- rownames(meta)
  wide <- left_join(meta, dat, by = "cell")
  long <- wide %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "log_expr") %>%
    mutate(panel = panel_name) %>%
    select(panel, cell, gene, log_expr, everything())
  long
}

endo <- subset(mb, subset = celltype == "Endothelial")
key_endo <- intersect(c("Cldn5","Ocln","Tjp1","Mfsd2a","Slc2a1","Plvap","Vcam1","Icam1"), rownames(endo))
endo_long <- extract_long(endo, key_endo, "Endothelial_key_BBB")

peri <- subset(mb, subset = celltype == "Pericytes")
key_peri <- intersect(c("Pdgfrb","Rgs5","Kcnj8","Anpep","Vtn","Notch3"), rownames(peri))
peri_long <- extract_long(peri, key_peri, "Pericyte_key_BBB")

out <- bind_rows(endo_long, peri_long) %>%
  arrange(panel, gene, sample, cell)

write.csv(out, file.path(outdir, "genes_logexpr_per_cell_with_meta.csv"), row.names = FALSE)

message(sprintf("Wrote %d rows (%d endothelial-panel cells x %d genes + %d pericyte-panel cells x %d genes) to %s",
                 nrow(out), ncol(endo), length(key_endo), ncol(peri), length(key_peri),
                 file.path(outdir, "genes_logexpr_per_cell_with_meta.csv")))

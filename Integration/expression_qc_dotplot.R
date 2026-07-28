## ============================================================
## Per-cell expression matrix + QC + DotPlot (bubble plot)
## Integrated dataset (Sample1 + Sample2)
## Output format matches reference: genes_logexpr_per_cell_wide
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(svglite)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
setwd(root)
dir.create("Gene_expression", showWarnings = FALSE)
dir.create("DotPlot",         showWarnings = FALSE)
dir.create("QC",              showWarnings = FALSE)

cat("Loading integrated object...\n")
mb <- readRDS(file.path(root, "mousebrain_integrated_final.rds"))
DefaultAssay(mb) <- "RNA"

## ---- gene_list (same as per-sample script) ----
gene_list <- c(
  "App","Apoe","Bace1",
  "Trem2","Tyrobp","Aif1","Cd68","Cx3cr1","Itgam",
  "Ccl3","Ccl4","Ccl5",
  "Gfap","Vim","S100b",
  "Casp1","Ptgs2","Stat3",
  "Tgfb1","Arg1","Mrc1",
  "Syn1","Rbfox3","Ptprc","Tubb3","Map2","Dcx","Syp",
  "Dlg4","Camk2a","Atf3","Fos","Jun","Bax","Egr1",
  "Casp3","P2ry12","Il1b","Tnf","Ccl2","Ifit3","Isg15",
  "Tmem119","Aldh1l1","C3","Serping1","Slc1a2",
  "Slc1a3","Gja1","Sox9","Aqp4",
  "Lcn2","Clu"
)
gene_list <- unique(gene_list)

all_genes <- rownames(GetAssayData(mb, assay = "RNA", layer = "data"))
present_genes <- intersect(gene_list, all_genes)
missing_genes <- setdiff(gene_list, all_genes)
cat("Present:", length(present_genes), " Missing:", length(missing_genes), "\n")
if (length(missing_genes)) {
  cat("  Missing:", paste(missing_genes, collapse = ", "), "\n")
  writeLines(missing_genes, "Gene_expression/missing_genes.txt")
}

## =========================================
## 1) Per-cell wide table (cell × gene)
## =========================================
cat("Building per-cell wide table...\n")
mat <- GetAssayData(mb, assay = "RNA", layer = "data")[present_genes, , drop = FALSE]
mat_t <- as.matrix(Matrix::t(mat))  # cells x genes

meta <- mb@meta.data
cells <- rownames(mat_t)

# cluster = integrated cluster id (numeric, matching reference style)
cluster_id <- as.integer(as.character(meta[cells, "clusters_integrated"]))

# Sheet 1: reference format (cell, cluster, gene1..gene49) -- exact match
wide_ref <- data.frame(cell = cells, cluster = cluster_id, mat_t,
                       check.names = FALSE, stringsAsFactors = FALSE)

# Sheet 2: with meta (sample, celltype, QC) inserted for context
wide_meta <- data.frame(
  cell        = cells,
  sample      = meta[cells, "sample"],
  cluster     = cluster_id,
  celltype    = as.character(meta[cells, "celltype"]),
  nCount_RNA  = meta[cells, "nCount_RNA"],
  nFeature_RNA= meta[cells, "nFeature_RNA"],
  percent.mt  = meta[cells, "percent.mt"],
  mat_t,
  check.names = FALSE, stringsAsFactors = FALSE
)

# Also a long/tidy version
cat("Building long table...\n")
long_df <- wide_meta |>
  dplyr::select(cell, sample, cluster, celltype, all_of(present_genes)) |>
  tidyr::pivot_longer(cols = all_of(present_genes),
                      names_to = "gene", values_to = "logexpr")

## Write outputs
cat("Writing xlsx...\n")
wb <- createWorkbook()
addWorksheet(wb, "genes_logexpr_per_cell_wide")
writeData(wb, "genes_logexpr_per_cell_wide", wide_ref)
freezePane(wb, "genes_logexpr_per_cell_wide", firstActiveRow = 2, firstActiveCol = 3)

addWorksheet(wb, "with_meta")
writeData(wb, "with_meta", wide_meta)
freezePane(wb, "with_meta", firstActiveRow = 2, firstActiveCol = 5)

saveWorkbook(wb, "Gene_expression/genes_logexpr_per_cell_wide.xlsx", overwrite = TRUE)

# CSV backups
write.csv(wide_ref,  "Gene_expression/genes_logexpr_per_cell_wide.csv",      row.names = FALSE)
write.csv(wide_meta, "Gene_expression/genes_logexpr_per_cell_with_meta.csv", row.names = FALSE)
write.csv(long_df,   "Gene_expression/genes_logexpr_per_cell_long.csv",      row.names = FALSE)

cat("Wide table: ", nrow(wide_ref), " cells x ", ncol(wide_ref), " cols (", length(present_genes), " genes)\n", sep = "")

## Also per-cluster and per-celltype averages (handy summaries)
cat("Building averages per cluster / per celltype...\n")
Idents(mb) <- "clusters_integrated"
avg_clu <- AverageExpression(mb, assays = "RNA", features = present_genes, layer = "data")$RNA
colnames(avg_clu) <- sub("^g", "", colnames(avg_clu))
write.csv(round(avg_clu, 4), "Gene_expression/genes_mean_logexpr_per_cluster_wide.csv", row.names = TRUE)

Idents(mb) <- "celltype"
avg_ct <- AverageExpression(mb, assays = "RNA", features = present_genes, layer = "data")$RNA
write.csv(round(avg_ct, 4), "Gene_expression/genes_mean_logexpr_per_celltype_wide.csv", row.names = TRUE)

Idents(mb) <- "celltype"
avg_ct_sample <- AverageExpression(mb, assays = "RNA", features = present_genes,
                                   layer = "data", group.by = c("celltype", "sample"))$RNA
write.csv(round(avg_ct_sample, 4),
          "Gene_expression/genes_mean_logexpr_per_celltype_x_sample.csv", row.names = TRUE)

## =========================================
## 2) QC plots
## =========================================
cat("QC plots...\n")
qc_theme <- theme_classic(base_size = 11)

# QC1: violin by sample
p_qc1 <- VlnPlot(mb,
                 features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                 group.by = "sample", pt.size = 0, ncol = 3) &
  theme(legend.position = "none")
ggsave("QC/QC1_VlnPlot_by_sample.svg", p_qc1,
       width = 10, height = 4, device = svglite::svglite)

# QC2: violin by celltype, faceted by sample
p_qc2 <- VlnPlot(mb,
                 features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                 group.by = "celltype", split.by = "sample",
                 pt.size = 0, ncol = 1, log = FALSE) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("QC/QC2_VlnPlot_by_celltype.svg", p_qc2,
       width = 11, height = 10, device = svglite::svglite)

# QC3: scatter nCount vs nFeature, colored by percent.mt
p_qc3 <- FeatureScatter(mb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "sample") + qc_theme
p_qc4 <- FeatureScatter(mb, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        group.by = "sample") + qc_theme
ggsave("QC/QC3_FeatureScatter.svg", p_qc3 | p_qc4,
       width = 11, height = 4.5, device = svglite::svglite)

# QC4: tSNE colored by QC metrics
p_qc5 <- FeaturePlot(mb,
                     features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                     reduction = "tsne.integrated", ncol = 3, order = TRUE) &
  scale_color_viridis_c() &
  theme_void(base_size = 10)
ggsave("QC/QC4_tSNE_QCmetrics.svg", p_qc5,
       width = 12, height = 4, device = svglite::svglite)

# QC summary table by sample x celltype
qc_summary <- meta |>
  dplyr::group_by(sample, celltype) |>
  dplyr::summarise(
    n_cells     = dplyr::n(),
    median_nCount   = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    median_pct_mt   = round(median(percent.mt), 3),
    .groups = "drop"
  ) |>
  dplyr::arrange(celltype, sample)
write.csv(qc_summary, "QC/QC_summary_sample_x_celltype.csv", row.names = FALSE)

## =========================================
## 3) Bubble / Dot plots
## =========================================
cat("DotPlots...\n")

# DotPlot 1: gene x celltype (combined samples)
Idents(mb) <- "celltype"
p_dot1 <- DotPlot(mb, features = present_genes, cols = c("lightgrey", "#b2182b"),
                  dot.scale = 6) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10)) +
  labs(title = "Gene expression by cell type (all cells)",
       x = NULL, y = "Cell type")
ggsave("DotPlot/DotPlot_gene_x_celltype.svg", p_dot1,
       width = 14, height = 5.5, device = svglite::svglite)

# DotPlot 2: gene x celltype, split by sample (side-by-side)
p_dot2 <- DotPlot(mb, features = present_genes, cols = c("#3b8bba", "#b2182b"),
                  dot.scale = 5, split.by = "sample") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  labs(title = "Gene expression by cell type, split by sample",
       x = NULL, y = "Cell type / Sample")
ggsave("DotPlot/DotPlot_gene_x_celltype_split_sample.svg", p_dot2,
       width = 14, height = 9, device = svglite::svglite)

# DotPlot 3: gene x cluster (integrated clusters)
Idents(mb) <- "clusters_integrated"
mb$clusters_integrated <- factor(mb$clusters_integrated,
                                 levels = sort(unique(as.integer(as.character(mb$clusters_integrated)))))
Idents(mb) <- "clusters_integrated"
p_dot3 <- DotPlot(mb, features = present_genes, cols = c("lightgrey", "#08519c"),
                  dot.scale = 5) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  labs(title = "Gene expression by integrated cluster",
       x = NULL, y = "Cluster")
ggsave("DotPlot/DotPlot_gene_x_cluster.svg", p_dot3,
       width = 14, height = 7, device = svglite::svglite)

# DotPlot 4: 4 focal cell types only (Excitatory, Inhibitory, Microglia, Astrocytes)
focal_ct <- c("Excitatory neurons","Inhibitory neurons","Microglia","Astrocytes")
mb_focal <- subset(mb, subset = celltype %in% focal_ct)
Idents(mb_focal) <- "celltype"
p_dot4 <- DotPlot(mb_focal, features = present_genes,
                  cols = c("#3b8bba","#b2182b"), dot.scale = 6,
                  split.by = "sample") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  labs(title = "Focal cell types × sample",
       x = NULL, y = "Cell type / Sample")
ggsave("DotPlot/DotPlot_focal_celltypes_split_sample.svg", p_dot4,
       width = 14, height = 6, device = svglite::svglite)

cat("\nALL DONE.\n")
cat("Outputs:\n")
cat("  Gene_expression/genes_logexpr_per_cell_wide.xlsx  (matches reference)\n")
cat("  Gene_expression/genes_logexpr_per_cell_*.csv\n")
cat("  Gene_expression/genes_mean_logexpr_per_{cluster,celltype}*.csv\n")
cat("  QC/QC{1..4}*.svg, QC/QC_summary_sample_x_celltype.csv\n")
cat("  DotPlot/DotPlot_*.svg\n")

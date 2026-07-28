# =============================================================
# snRNA-seq Integration: Sample1 + Sample2 (Seurat v5)
# Reproduces the Fig.3-style panels from:
#   Nature 2025, "Lithium deficiency and the onset of Alzheimer's disease"
#   https://www.nature.com/articles/s41586-025-09335-x/figures/3
# Panels:
#   a. Joint t-SNE coloured by cell type + t-SNE split by sample
#   b. DEG counts per cell type (up / down) bar chart
#   c. Heatmap of top DEGs per cell type (sample-wise log2FC)
#   d. Microglia focus: Cx3cr1 (homeostatic) vs Apoe (DAM/AD-like)
# Visualization: t-SNE (matches the per-sample pipeline). Clustering is
# graph-based (FindNeighbors/FindClusters on the integrated reduction)
# and does not depend on the 2-D embedding.
# =============================================================

library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

# -------------------- paths --------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
root <- normalizePath(file.path(script_dir, ".."))
setwd(script_dir)
source(file.path(script_dir, "scDblFinder_doublet_detection.R"))

h5_s1 <- file.path(root, "Sample1", "filtered_feature_bc_matrix.h5")
h5_s2 <- file.path(root, "Sample2", "filtered_feature_bc_matrix.h5")

dir.create("QC",        showWarnings = FALSE, recursive = TRUE)
dir.create("TSNE",      showWarnings = FALSE, recursive = TRUE)
dir.create("Celltypes", showWarnings = FALSE, recursive = TRUE)
dir.create("DEG",       showWarnings = FALSE, recursive = TRUE)
dir.create("Microglia", showWarnings = FALSE, recursive = TRUE)

palette_DotPlot <- c("#4589C6", "white", "#CE544D")

# ============================================================
# 1) Load each sample and remove predicted doublets independently
# ============================================================
s1.data <- Read10X_h5(h5_s1)
s2.data <- Read10X_h5(h5_s2)

# Doublets are called separately for each physical 10x library before merging.
# The raw count matrices are retained on disk; only predicted singlets enter
# normalization, integration, clustering, and downstream analyses.
dbl_s1 <- run_scdblfinder(s1.data, library_id = "Sample1")
dbl_s2 <- run_scdblfinder(s2.data, library_id = "Sample2")

write.csv(
  rbind(dbl_s1$summary, dbl_s2$summary),
  "QC/scDblFinder_summary.csv",
  row.names = FALSE
)
write.csv(
  rbind(dbl_s1$calls, dbl_s2$calls),
  "QC/scDblFinder_calls.csv",
  row.names = FALSE
)

s1.data <- s1.data[, dbl_s1$singlet_barcodes, drop = FALSE]
s2.data <- s2.data[, dbl_s2$singlet_barcodes, drop = FALSE]

s1 <- CreateSeuratObject(counts = s1.data, project = "Sample1",
                         min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = s2.data, project = "Sample2",
                         min.cells = 3, min.features = 200)
s1$sample <- "Sample1"
s2$sample <- "Sample2"

# Merge the retained singlets and continue with the standard Seurat workflow.
mb <- merge(s1, y = s2, add.cell.ids = c("S1", "S2"),
            project = "mousebrain_integrated")

# QC metric
mb[["percent.mt"]] <- PercentageFeatureSet(mb, pattern = "^mt-")

VlnPlot(mb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample", ncol = 3)
ggsave("QC/VlnPlot_by_sample.svg", width = 9, height = 6)

# ============================================================
# 2) Seurat v5: merge() above already produced per-sample layers
#    (counts.Sample1 / counts.Sample2). Normalize / HVG / Scale / PCA
#    on the unintegrated object, then IntegrateLayers().
# ============================================================
mb <- NormalizeData(mb)
mb <- FindVariableFeatures(mb, selection.method = "vst", nfeatures = 2000)
mb <- ScaleData(mb)
mb <- RunPCA(mb, npcs = 30, verbose = FALSE)

# Unintegrated baseline (for the "before" t-SNE, useful to visualize batch)
mb <- FindNeighbors(mb, reduction = "pca", dims = 1:30)
mb <- FindClusters(mb, resolution = 0.5, cluster.name = "clusters_unintegrated")
mb <- RunTSNE(mb, reduction = "pca", dims = 1:30,
              reduction.name = "tsne.unintegrated",
              reduction.key  = "tsneUnint_")

DimPlot(mb, reduction = "tsne.unintegrated", group.by = "sample") +
  ggtitle("Before integration")
ggsave("TSNE/tSNE_unintegrated_by_sample.svg", width = 6, height = 6)

# CCA integration (Seurat v5). Switch to RPCAIntegration / HarmonyIntegration
# if memory is tight or batch is severe.
mb <- IntegrateLayers(
  object         = mb,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  verbose        = FALSE
)

# Re-join layers so that downstream DE / averaging works on a single matrix.
mb[["RNA"]] <- JoinLayers(mb[["RNA"]])

mb <- FindNeighbors(mb, reduction = "integrated.cca", dims = 1:30)
mb <- FindClusters(mb, resolution = 0.5, cluster.name = "clusters_integrated")
mb <- RunTSNE(mb, reduction = "integrated.cca", dims = 1:30,
              reduction.name = "tsne.integrated",
              reduction.key  = "tsneInt_")

# Set the integrated clusters as active identity for the rest of the script.
Idents(mb) <- "clusters_integrated"

# ============================================================
# 3) Cell-type annotation (reuse the marker panel from per-sample scripts)
#    Score each integrated cluster against grouped cell-type marker sets.
# ============================================================
celltypes_grouped <- list(
  "Excitatory neurons" = c("Slc17a7","Camk2a","Tbr1","Reln","Satb2"),
  "Inhibitory neurons" = c("Gad1","Gad2","Slc6a1"),
  "Astrocytes"         = c("S100b","Gfap","Aqp4","Aldoc","Slc1a3"),
  "Oligodendrocytes"   = c("Mog","Mag","Mbp","Plp1","Cnp"),
  "OPC"                = c("Pdgfra","Cspg4"),
  "Microglia"          = c("Itgam","Cx3cr1","P2ry12","Aif1","Tmem119"),
  "Endothelial"        = c("Pecam1","Kdr","Cldn5"),
  "Pericytes"          = c("Rgs5","Kcnj8","Pdgfrb"),
  "Ependymal"          = c("Foxj1","Dynlrb2"),
  "Choroid plexus"     = c("Ttr","Klk8")
)

avg_expr <- AverageExpression(mb, assays = "RNA", slot = "data")$RNA
all_genes_avg <- rownames(avg_expr)
# Seurat 5's AverageExpression prepends "g" to numeric cluster names ("0" → "g0")
# to make valid variable names. Strip it so the keys line up with Idents(mb).
cluster_ids   <- sub("^g", "", colnames(avg_expr))
colnames(avg_expr) <- cluster_ids

score_mat <- sapply(celltypes_grouped, function(genes) {
  g <- intersect(genes, all_genes_avg)
  if (length(g) == 0) return(rep(0, length(cluster_ids)))
  colMeans(avg_expr[g, , drop = FALSE])
})

best_celltype <- apply(score_mat, 1, function(x) names(which.max(x)))
second_best   <- apply(score_mat, 1, function(x) names(sort(x, decreasing = TRUE)[2]))
top_score     <- apply(score_mat, 1, max)

celltype_annotation <- data.frame(
  cluster       = cluster_ids,
  celltype      = best_celltype,
  second_choice = second_best,
  top_score     = round(top_score, 4),
  row.names     = NULL
)
write.csv(celltype_annotation,
          "Celltypes/integrated_celltype_scores.csv", row.names = FALSE)
write.csv(round(score_mat, 4),
          "Celltypes/integrated_celltype_scorematrix.csv")

# Attach cell type label to each cell
names(best_celltype) <- cluster_ids
mb$celltype <- unname(best_celltype[as.character(Idents(mb))])

# Save the integrated + annotated object so downstream can resume here.
saveRDS(mb, file = "mousebrain_integrated.rds")

# ============================================================
# Panel a — Joint t-SNE by cell type + split by sample
# ============================================================
p_celltype <- DimPlot(mb, reduction = "tsne.integrated",
                      group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Integrated — cell type")
p_sample <- DimPlot(mb, reduction = "tsne.integrated",
                    group.by = "sample") +
  ggtitle("Integrated — sample")

(p_celltype | p_sample)
ggsave("TSNE/tSNE_integrated_celltype_and_sample.svg",
       width = 14, height = 6)

DimPlot(mb, reduction = "tsne.integrated",
        group.by = "celltype", split.by = "sample",
        label = TRUE, repel = TRUE) +
  ggtitle("Cell type, split by sample")
ggsave("TSNE/tSNE_integrated_split_by_sample.svg",
       width = 14, height = 6)

# ============================================================
# Panel b — DEG counts per cell type (Sample2 vs Sample1), up / down
# ============================================================
# Use FindMarkers' built-in group.by + subset.ident: this avoids fragile
# subset() + Idents<- juggling and works correctly after JoinLayers.
Idents(mb) <- "celltype"
celltypes_present <- sort(unique(mb$celltype))

deg_per_celltype <- lapply(celltypes_present, function(ct) {
  # require both samples represented in this cell type with >= 3 cells each
  tab <- table(mb$celltype == ct, mb$sample)
  if (!("TRUE" %in% rownames(tab))) return(NULL)
  n_per_sample <- tab["TRUE", ]
  if (length(n_per_sample) < 2 || min(n_per_sample) < 3) {
    message(sprintf("Skip '%s': counts per sample = %s",
                    ct, paste(n_per_sample, collapse = "/")))
    return(NULL)
  }
  res <- tryCatch(
    FindMarkers(mb,
                ident.1       = "Sample2",
                ident.2       = "Sample1",
                group.by      = "sample",
                subset.ident  = ct,
                logfc.threshold = 0.25,
                min.pct         = 0.1),
    error = function(e) { message("FindMarkers failed for '", ct, "': ",
                                  conditionMessage(e)); NULL }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$gene     <- rownames(res)
  res$celltype <- ct
  res
})
names(deg_per_celltype) <- celltypes_present
deg_all <- bind_rows(deg_per_celltype)

write.csv(deg_all, "DEG/DEG_sample2_vs_sample1_by_celltype.csv",
          row.names = FALSE)

if (nrow(deg_all) == 0) {
  warning("No DEGs returned for any cell type — skipping Panels b/c/c2.")
  saveRDS(mb, "mousebrain_integrated_final.rds")
  message("Run halted early. Check the log for 'FindMarkers failed' lines.")
  quit(save = "no", status = 0)
}

# Significant DEGs: adj p < 0.05 & |log2FC| > 0.25
sig <- deg_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))

deg_counts <- sig %>%
  count(celltype, direction) %>%
  mutate(n_signed = ifelse(direction == "Down", -n, n))

ggplot(deg_counts, aes(x = reorder(celltype, abs(n_signed)),
                       y = n_signed, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#CE544D", "Down" = "#4589C6")) +
  labs(x = NULL, y = "# DEGs (Sample2 vs Sample1)",
       title = "DEGs per cell type") +
  theme_classic(base_size = 12)
ggsave("DEG/DEG_counts_per_celltype.svg", width = 8, height = 5)
write.csv(deg_counts, "DEG/DEG_counts_per_celltype.csv", row.names = FALSE)

# ============================================================
# Panel c — Heatmap of top DEGs across cell types
# ============================================================
top_deg <- sig %>%
  group_by(celltype) %>%
  slice_max(abs(avg_log2FC), n = 8, with_ties = FALSE) %>%
  ungroup()

# log2FC matrix: rows = genes, cols = cell types
fc_mat <- deg_all %>%
  filter(gene %in% top_deg$gene) %>%
  select(gene, celltype, avg_log2FC) %>%
  pivot_wider(names_from = celltype, values_from = avg_log2FC,
              values_fill = 0) %>%
  as.data.frame()
rownames(fc_mat) <- fc_mat$gene
fc_mat$gene <- NULL

fc_long <- fc_mat %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "celltype", values_to = "log2FC")

ggplot(fc_long, aes(x = celltype, y = gene, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4589C6", mid = "white", high = "#CE544D",
                       midpoint = 0) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Top DEGs (Sample2 vs Sample1) per cell type")
ggsave("DEG/Heatmap_top_DEGs_by_celltype.svg",
       width = 8, height = max(6, 0.18 * nrow(fc_mat)))

write.csv(fc_mat, "DEG/Heatmap_top_DEGs_matrix.csv", row.names = TRUE)

# ============================================================
# Panel c2 — GO enrichment of DEGs (BP) per cell type, up vs down
# Requires: clusterProfiler, org.Mm.eg.db, AnnotationDbi
#   install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler","org.Mm.eg.db","AnnotationDbi"))
# ============================================================
has_go <- requireNamespace("clusterProfiler", quietly = TRUE) &&
          requireNamespace("org.Mm.eg.db",    quietly = TRUE) &&
          requireNamespace("AnnotationDbi",   quietly = TRUE)

if (!has_go) {
  message("Skipping GO enrichment: install clusterProfiler + org.Mm.eg.db to enable.")
} else {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
  })
  dir.create("GO", showWarnings = FALSE, recursive = TRUE)

  universe_sym <- rownames(mb)
  universe_eid <- AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = universe_sym, keytype = "SYMBOL",
    column = "ENTREZID", multiVals = "first"
  )
  universe_eid <- unique(na.omit(unname(universe_eid)))

  go_results <- list()
  for (ct in unique(sig$celltype)) {
    for (dir_ in c("Up", "Down")) {
      genes_sym <- sig %>%
        filter(celltype == ct, direction == dir_) %>%
        pull(gene) %>% unique()
      if (length(genes_sym) < 10) next
      genes_eid <- AnnotationDbi::mapIds(
        org.Mm.eg.db, keys = genes_sym, keytype = "SYMBOL",
        column = "ENTREZID", multiVals = "first"
      )
      genes_eid <- unique(na.omit(unname(genes_eid)))
      if (length(genes_eid) < 10) next

      ego <- tryCatch(
        enrichGO(gene          = genes_eid,
                 universe      = universe_eid,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = "ENTREZID",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE),
        error = function(e) NULL
      )
      if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next

      tab <- as.data.frame(ego)
      tab$celltype  <- ct
      tab$direction <- dir_
      go_results[[paste(ct, dir_, sep = "_")]] <- tab
    }
  }

  if (length(go_results) > 0) {
    go_all <- bind_rows(go_results)
    write.csv(go_all, "GO/GO_BP_per_celltype_direction.csv", row.names = FALSE)

    top_go <- go_all %>%
      group_by(celltype, direction) %>%
      slice_min(p.adjust, n = 5, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        group   = paste(celltype, direction, sep = " | "),
        log10p  = -log10(p.adjust),
        Term    = factor(Description, levels = rev(unique(Description)))
      )

    ggplot(top_go,
           aes(x = log10p, y = Term, fill = direction)) +
      geom_col() +
      facet_wrap(~ celltype, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = c("Up" = "#CE544D", "Down" = "#4589C6")) +
      labs(x = expression(-log[10]~adj.~italic(P)), y = NULL,
           title = "Top GO:BP terms (Sample2 vs Sample1) per cell type") +
      theme_classic(base_size = 10) +
      theme(strip.text = element_text(face = "bold"))
    ggsave("GO/GO_BP_top_terms_per_celltype.svg",
           width = 14, height = max(6, 1.6 * length(unique(top_go$celltype))))
  } else {
    message("GO enrichment ran but no terms passed thresholds.")
  }
}

# ============================================================
# Panel d — Microglia focus (Cx3cr1 vs Apoe), echoing the
#           homeostatic → DAM/AD-like shift in the paper.
# ============================================================
microglia_markers <- c("Cx3cr1", "P2ry12", "Tmem119",   # homeostatic
                       "Apoe", "Trem2", "Tyrobp",       # DAM / AD-like
                       "Itgax", "Cst7", "Spp1")         # broader DAM panel

present <- intersect(microglia_markers, rownames(mb))

FeaturePlot(mb, features = c("Cx3cr1", "Apoe"),
            reduction = "tsne.integrated", split.by = "sample",
            order = TRUE) &
  theme(legend.position = "right")
ggsave("Microglia/FeaturePlot_Cx3cr1_Apoe_split_by_sample.svg",
       width = 12, height = 8)

VlnPlot(mb, features = present, group.by = "celltype", split.by = "sample",
        pt.size = 0, ncol = 3) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Microglia/VlnPlot_microglia_markers_split.svg",
       width = 14, height = 10)

# Microglia subset re-clustering: are there sample-specific microglia states?
mg <- subset(mb, subset = celltype == "Microglia")
if (ncol(mg) >= 50) {
  mg <- FindVariableFeatures(mg, nfeatures = 1500)
  mg <- ScaleData(mg)
  mg <- RunPCA(mg, npcs = 20, verbose = FALSE)
  mg <- FindNeighbors(mg, reduction = "pca", dims = 1:15)
  mg <- FindClusters(mg, resolution = 0.4)
  mg <- RunTSNE(mg, reduction = "pca", dims = 1:15,
                reduction.name = "tsne", reduction.key = "tSNE_")

  DimPlot(mg, reduction = "tsne",
          group.by = "seurat_clusters", label = TRUE) +
    DimPlot(mg, reduction = "tsne", group.by = "sample")
  ggsave("Microglia/tSNE_microglia_subclusters.svg",
         width = 12, height = 5)

  # microglia subcluster composition by sample
  mg_comp <- mg@meta.data %>%
    count(sample, seurat_clusters) %>%
    group_by(sample) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

  ggplot(mg_comp, aes(x = seurat_clusters, y = frac, fill = sample)) +
    geom_col(position = "dodge") +
    labs(x = "Microglia subcluster", y = "Fraction of microglia",
         title = "Microglia subcluster composition by sample") +
    theme_classic()
  ggsave("Microglia/Microglia_subcluster_composition.svg",
         width = 8, height = 5)

  write.csv(mg_comp, "Microglia/Microglia_subcluster_composition.csv",
            row.names = FALSE)
  saveRDS(mg, "Microglia/microglia_subset.rds")
}

# DotPlot of the AD / inflammation panel from the per-sample script,
# split by sample so reader can compare side by side.
interests <- c(
  "App","Apoe","Bace1",
  "Trem2","Tyrobp","Cd68","Cx3cr1","Itgam","P2ry12","Tmem119",
  "Ccl3","Ccl4","Ccl5",
  "Gfap","Vim","S100b","Aqp4","Slc1a3",
  "Casp1","Ptgs2","Stat3",
  "Tgfb1","Arg1","Mrc1"
)
interests <- intersect(interests, rownames(mb))

DotPlot(mb, features = interests, group.by = "celltype",
        split.by = "sample", cols = c("#4589C6", "#CE544D")) +
  RotatedAxis()
ggsave("DEG/DotPlot_AD_panel_split_by_sample.svg",
       width = 14, height = 7)

# ============================================================
# Overall cell-type composition by sample
# ============================================================
comp <- mb@meta.data %>%
  count(sample, celltype) %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

ggplot(comp, aes(x = celltype, y = frac, fill = sample)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells",
       title = "Cell-type composition by sample")
ggsave("TSNE/Celltype_composition_by_sample.svg", width = 9, height = 5)
write.csv(comp, "TSNE/Celltype_composition_by_sample.csv", row.names = FALSE)

# ============================================================
# Final save
# ============================================================
saveRDS(mb, "mousebrain_integrated_final.rds")
message("Integration complete. Outputs under: ", getwd())

# =============================================================
# Downstream-only: resume from mousebrain_integrated.rds and run
# Panels b / c / c2 / d + composition. Use this when the
# integration step has already completed but the DEG / GO /
# microglia panels still need to be generated.
# =============================================================
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
root <- normalizePath(file.path(script_dir, ".."))
setwd(script_dir)

dir.create("TSNE",      showWarnings = FALSE, recursive = TRUE)
dir.create("Celltypes", showWarnings = FALSE, recursive = TRUE)
dir.create("DEG",       showWarnings = FALSE, recursive = TRUE)
dir.create("Microglia", showWarnings = FALSE, recursive = TRUE)

palette_DotPlot <- c("#4589C6", "white", "#CE544D")

# ---------- Load + (re)annotate -----------------------------
mb <- readRDS("mousebrain_integrated.rds")

# the saved rds (from the failed run) does not carry $celltype yet.
# Rebuild from Celltypes/integrated_celltype_scores.csv produced earlier.
if (!"celltype" %in% colnames(mb@meta.data)) {
  ct_tab <- read.csv("Celltypes/integrated_celltype_scores.csv",
                     stringsAsFactors = FALSE)
  # Seurat 5's AverageExpression prepends "g" to numeric cluster names —
  # strip it so the keys line up with Idents(mb).
  ct_tab$cluster <- sub("^g", "", as.character(ct_tab$cluster))
  ct_map <- setNames(ct_tab$celltype, ct_tab$cluster)
  Idents(mb) <- "clusters_integrated"
  mb$celltype <- unname(ct_map[as.character(Idents(mb))])
  stopifnot("celltype mapping failed (all NA)" =
            !all(is.na(mb$celltype)))
  message("Re-attached celltype from CSV: ",
          length(unique(na.omit(mb$celltype))), " unique labels across ",
          sum(!is.na(mb$celltype)), " cells; ",
          sum(is.na(mb$celltype)), " cells unmapped.")
}
Idents(mb) <- "celltype"

# ---------- Compute t-SNE (overwrite any cached UMAP) -------
# Clustering already done on the graph from `integrated.cca`; here we
# only add the 2-D embedding for visualization. The unintegrated t-SNE
# on `pca` is for the "before integration" panel.
mb <- RunTSNE(mb, reduction = "pca", dims = 1:30,
              reduction.name = "tsne.unintegrated",
              reduction.key  = "tsneUnint_")
mb <- RunTSNE(mb, reduction = "integrated.cca", dims = 1:30,
              reduction.name = "tsne.integrated",
              reduction.key  = "tsneInt_")

# Regenerate the t-SNE panels that the upstream script would have made.
DimPlot(mb, reduction = "tsne.unintegrated", group.by = "sample") +
  ggtitle("Before integration")
ggsave("TSNE/tSNE_unintegrated_by_sample.svg", width = 6, height = 6)

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
# Panel b — DEG counts per cell type (Sample2 vs Sample1)
# ============================================================
celltypes_present <- sort(unique(mb$celltype))

deg_per_celltype <- lapply(celltypes_present, function(ct) {
  tab <- table(mb$celltype == ct, mb$sample)
  if (!("TRUE" %in% rownames(tab))) return(NULL)
  n_per_sample <- tab["TRUE", ]
  if (length(n_per_sample) < 2 || min(n_per_sample) < 3) {
    message(sprintf("Skip '%s': counts per sample = %s",
                    ct, paste(n_per_sample, collapse = "/")))
    return(NULL)
  }
  message("DE for: ", ct, "  (Sample1=", n_per_sample["Sample1"],
          ", Sample2=", n_per_sample["Sample2"], ")")
  res <- tryCatch(
    FindMarkers(mb,
                ident.1       = "Sample2",
                ident.2       = "Sample1",
                group.by      = "sample",
                subset.ident  = ct,
                logfc.threshold = 0.25,
                min.pct         = 0.1),
    error = function(e) { message("  failed: ", conditionMessage(e)); NULL }
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
message("DEG table rows: ", nrow(deg_all),
        "; cell types with DEGs: ", length(unique(deg_all$celltype)))

stopifnot(nrow(deg_all) > 0)

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
# Panel c2 — GO BP enrichment per cell type × direction
# ============================================================
has_go <- requireNamespace("clusterProfiler", quietly = TRUE) &&
          requireNamespace("org.Mm.eg.db",    quietly = TRUE) &&
          requireNamespace("AnnotationDbi",   quietly = TRUE)

if (!has_go) {
  message("Skip GO: install clusterProfiler + org.Mm.eg.db to enable.")
} else {
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Mm.eg.db)
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
        enrichGO(gene = genes_eid, universe = universe_eid,
                 OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
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
      mutate(log10p = -log10(p.adjust),
             Term   = factor(Description, levels = rev(unique(Description))))
    ggplot(top_go, aes(x = log10p, y = Term, fill = direction)) +
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
    message("GO ran but no terms passed thresholds.")
  }
}

# ============================================================
# Panel d — Microglia focus
# ============================================================
microglia_markers <- c("Cx3cr1","P2ry12","Tmem119",
                       "Apoe","Trem2","Tyrobp",
                       "Itgax","Cst7","Spp1")
present <- intersect(microglia_markers, rownames(mb))

FeaturePlot(mb, features = c("Cx3cr1","Apoe"),
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
  ggsave("Microglia/tSNE_microglia_subclusters.svg", width = 12, height = 5)

  mg_comp <- mg@meta.data %>%
    count(sample, seurat_clusters) %>%
    group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
  ggplot(mg_comp, aes(x = seurat_clusters, y = frac, fill = sample)) +
    geom_col(position = "dodge") +
    labs(x = "Microglia subcluster", y = "Fraction of microglia",
         title = "Microglia subcluster composition by sample") +
    theme_classic()
  ggsave("Microglia/Microglia_subcluster_composition.svg", width = 8, height = 5)
  write.csv(mg_comp, "Microglia/Microglia_subcluster_composition.csv",
            row.names = FALSE)
  saveRDS(mg, "Microglia/microglia_subset.rds")
}

# AD panel DotPlot split by sample
interests <- c("App","Apoe","Bace1",
               "Trem2","Tyrobp","Cd68","Cx3cr1","Itgam","P2ry12","Tmem119",
               "Ccl3","Ccl4","Ccl5",
               "Gfap","Vim","S100b","Aqp4","Slc1a3",
               "Casp1","Ptgs2","Stat3",
               "Tgfb1","Arg1","Mrc1")
interests <- intersect(interests, rownames(mb))

DotPlot(mb, features = interests, group.by = "celltype",
        split.by = "sample", cols = c("#4589C6", "#CE544D")) + RotatedAxis()
ggsave("DEG/DotPlot_AD_panel_split_by_sample.svg", width = 14, height = 7)

# Cell-type composition by sample
comp <- mb@meta.data %>%
  count(sample, celltype) %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()

ggplot(comp, aes(x = celltype, y = frac, fill = sample)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells",
       title = "Cell-type composition by sample")
ggsave("TSNE/Celltype_composition_by_sample.svg", width = 9, height = 5)
write.csv(comp, "TSNE/Celltype_composition_by_sample.csv", row.names = FALSE)

saveRDS(mb, "mousebrain_integrated_final.rds")
message("Downstream complete. Outputs under: ", getwd())

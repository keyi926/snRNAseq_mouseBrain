library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

#### Setup the Seurat Object####
# Load the mousebrain dataset
setwd("/Users/keyi/Documents/1_Data/single_cell_Sequencing/Single_cell_seq/R/Sample2")
mousebrain.data <- Read10X_h5("filtered_feature_bc_matrix.h5")

###############################################################
# Initialize the Seurat object with the raw (non-normalized data).
mousebrain <- CreateSeuratObject(counts = mousebrain.data, project = "mousebrain",
 min.cells = 3, min.features = 200)
mousebrain

# QC
dir.create("QC", showWarnings = FALSE, recursive = TRUE)

mousebrain[["percent.mt"]] <- PercentageFeatureSet(mousebrain, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot1 <- VlnPlot(mousebrain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot1
ggsave("QC/VlnPlot.svg", width = 9, height = 6)

plot1 <- FeatureScatter(mousebrain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mousebrain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ScatterPlot1 <- plot1 + plot2
ScatterPlot1
ggsave("QC/ScatterPlot.svg", width = 9, height = 6)

# Normalizing the data
mousebrain <- NormalizeData(mousebrain, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection
mousebrain <- FindVariableFeatures(mousebrain, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mousebrain), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mousebrain)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
VFPlot1 <- plot1 + plot2
VFPlot1
ggsave("QC/VFPlot.svg", width = 9, height = 6)

#Scalling the data
all.genes <- rownames(mousebrain)
mousebrain <- ScaleData(mousebrain, features = all.genes)

#Perform linear dimensional reduction (PCA)
mousebrain <- RunPCA(mousebrain, features = VariableFeatures(object = mousebrain))
print(mousebrain[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoading1 <- VizDimLoadings(mousebrain, dims = 1:2, reduction = "pca")
VizDimLoading1
ggsave("QC/PCA_VizDimLoading.svg", width = 9, height = 6)
DimPlot1 <- DimPlot(mousebrain, reduction = "pca") + NoLegend()
DimPlot1
ggsave("QC/PCA_DimPlot.svg", width = 6, height = 6)
DimHeatmap1 <- DimHeatmap(mousebrain, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap1
ggsave("QC/DimHeatmap.jpg", width = 9, height = 6)

#Determine the ‘dimensionality’ of the dataset
ElbowPlot1 <- ElbowPlot(mousebrain)
ElbowPlot1
ggsave("QC/ElbowPlot.svg", width = 9, height = 6)

#### Clustering ####
mousebrain <- FindNeighbors(mousebrain, dims = 1:10)
mousebrain <- FindClusters(mousebrain, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(mousebrain), 5)

#t-SNE
mousebrain <- RunTSNE(mousebrain, dims = 1:10)
DimPlot_tsne <- DimPlot(mousebrain, reduction = "tsne")
DimPlot_tsne
ggsave("DimPlot_tsne.svg", width = 6, height = 6)
saveRDS(mousebrain, file = "mousebrain_tsne.rds")

#### Cluster biomarkers#### 
#### Time-consuming step
mousebrain.markers <- FindAllMarkers(mousebrain, only.pos = TRUE)
mousebrain.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
head(mousebrain.markers)

write.csv(x = mousebrain.markers, file = "mousebrain_markers.csv")
saveRDS(mousebrain.markers, file = "mousebrain_markers.rds")

###############################################################
################### start from .rds############################
###############################################################
mousebrain <- readRDS("/Users/keyi/Documents/1_Data/single_cell_Sequencing/Single_cell_seq/R/Sample2/mousebrain_tsne.rds")
mousebrain.markers <- readRDS("/Users/keyi/Documents/1_Data/single_cell_Sequencing/Single_cell_seq/R/Sample2/mousebrain_markers.rds")

###############################################################
###############################################################
###############################################################

#################### by Celltypes ####################
dir.create("Celltypes", showWarnings = FALSE, recursive = TRUE)

celltypes <- c("Slc17a7","Camk2a","Tbr1","Reln","Satb2", #Excitatory neurons
              "Gad1","Gad2","Slc6a1", #Inhibitory neurons
              #"Map2", #Neurons
              "S100b","Gfap", "Aqp4","Aldoc","Slc1a3", #Astrocytes,
              "Mog","Mag","Mbp","Plp1","Cnp", #Oligodendrocytes
              "Pdgfra","Cspg4", #OPC
              "Itgam", "Cx3cr1","P2ry12","Aif1","Tmem119", #Microglia
              "Pecam1","Kdr","Cldn5", #Endothelial
              "Rgs5","Kcnj8","Pdgfrb", #Pericytes
              "Foxj1","Dynlrb2", #Ependymal
              "Ttr","Klk8" #Choroid plexus
)

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

############## Plots
palette_DotPlot <- c("#4589C6", "white" ,"#CE544D")

DotPlot_celltypes <- DotPlot(mousebrain, features = celltypes) + RotatedAxis() +
  scale_color_gradientn(colors = palette_DotPlot)
DotPlot_celltypes
ggsave("Celltypes/DotPlot_celltypes.svg", width = 12, height = 6)
DotPlot_celltypes_grouped <- DotPlot(mousebrain, features = celltypes_grouped) + RotatedAxis() +
  scale_color_gradientn(colors = palette_DotPlot)
DotPlot_celltypes_grouped
ggsave("Celltypes/DotPlot_celltypes_grouped.svg", width = 20, height = 6)

DoHeatmap_celltypes <- DoHeatmap(subset(mousebrain, downsample = 100), features = celltypes, size = 3)
DoHeatmap_celltypes
ggsave("Celltypes/DoHeatmap_celltypes.svg", width = 9, height = 6)

#RidgePlot_celltypes <- RidgePlot(mousebrain, features = celltypes, ncol = 1)
#RidgePlot_celltypes
#VlnPlot_celltypes <- VlnPlot(mousebrain, features = celltypes, slot = "counts", log = TRUE)
#VlnPlot_celltypes
#FeaturePlot_celltypes <- FeaturePlot(mousebrain, features = celltypes, slot = "counts")
#FeaturePlot_celltypes

############### Tables

avg_expr <- AverageExpression(mousebrain, assays = "RNA", slot = "data")$RNA
all_genes_avg <- rownames(avg_expr)
cluster_ids <- colnames(avg_expr)

score_mat <- sapply(celltypes_grouped, function(genes) {
  g <- intersect(genes, all_genes_avg)
  if (length(g) == 0) return(rep(0, length(cluster_ids)))
  colMeans(avg_expr[g, , drop = FALSE])
})

best_celltype <- apply(score_mat, 1, function(x) names(which.max(x)))
second_best   <- apply(score_mat, 1, function(x) names(sort(x, decreasing=TRUE)[2]))
top_score     <- apply(score_mat, 1, max)

celltype_annotation <- data.frame(
  cluster       = cluster_ids,
  celltype      = best_celltype,
  second_choice = second_best,
  top_score     = round(top_score, 4),
  row.names     = NULL
)

write.csv(
  data.frame(x = best_celltype, row.names = cluster_ids),
  file = "Celltypes/mousebrain_celltype.csv"
)

write.csv(celltype_annotation, file = "Celltypes/mousebrain_celltype_scores.csv", row.names = FALSE)
write.csv(round(score_mat, 4),  file = "Celltypes/mousebrain_celltype_scorematrix.csv")

message("✓ Exported: mousebrain_celltype.csv / _scores.csv / _scorematrix.csv")

# ==================== by Regions ============================
dir.create("Regions", showWarnings = FALSE, recursive = TRUE)

regions <- c(
  # Cortex 
  "Satb2", "Cux1", "Rorb", 
  # Hippocampus 
  "Prox1",  # DG granule cells
  "Calb1",  # CA1 pyramidal
  "Grin2a", "Grin2b", # NMDA receptor subunits
  # Striatum 
  "Drd1", "Drd2", "Tac1", "Penk",  # medium spiny neurons
  # Thalamus 
  "Tcf7l2", "Pou2f2", "Slc17a6",  # excitatory relay neurons
  # Hypothalamus
  "Avp", "Oxt",
  # Cerebellum 
  "Pcp2",  # Purkinje cells
  "Lhx1", # granule cells / inhibitory neurons
  # Amygdala 
  "Nr2f2", "Foxp2",
  # Olfactory bulb
  "Omp", "Gng8", # olfactory sensory neurons
  # Brainstem
  "Dbh", "Tph2"# monoaminergic neurons
)

regions_grouped <- list(
  "Cortex"          = c("Satb2", "Cux1", "Rorb"),
  "Hippocampus"     = c("Prox1", "Calb1", "Grin2a", "Grin2b"),
  "Striatum"        = c("Drd1", "Drd2", "Tac1", "Penk"),
  "Thalamus"        = c("Tcf7l2", "Pou2f2", "Slc17a6"),
  "Hypothalamus"    = c("Avp", "Oxt"),
  "Cerebellum"      = c("Pcp2", "Lhx1"),
  "Amygdala"        = c("Nr2f2", "Foxp2"),
  "Olfactory bulb"  = c("Omp", "Gng8"),
  "Brainstem"       = c("Dbh", "Tph2")
)

DotPlot_regions <- DotPlot(mousebrain, features = regions) + RotatedAxis() +
  scale_color_gradientn(colors = palette_DotPlot)
DotPlot_regions
ggsave("Regions/DotPlot_regions.svg", width = 12, height = 6)
DotPlot_regions_grouped <- DotPlot(mousebrain, features = regions_grouped) + RotatedAxis() +
  scale_color_gradientn(colors = palette_DotPlot)
DotPlot_regions_grouped
ggsave("Regions/DotPlot_regions_grouped.svg", width = 20, height = 6)

DoHeatmap_regions <- DoHeatmap(subset(mousebrain, downsample = 100), features = regions, size = 3)
DoHeatmap_regions
ggsave("Regions/DoHeatmap_regions.svg", width = 9, height = 6)

#FeaturePlot_regions <- FeaturePlot(mousebrain, features = regions)
#FeaturePlot_regions
#ggsave("FeaturePlot_regions.jpg", width = 6, height = 3)

avg_expr <- AverageExpression(mousebrain, assays = "RNA", slot = "data")$RNA
all_genes_avg <- rownames(avg_expr)
cluster_ids <- colnames(avg_expr)

region_score_mat <- sapply(regions_grouped, function(genes) {
  g <- intersect(genes, all_genes_avg)
  if (length(g) == 0) return(rep(0, length(cluster_ids)))
  colMeans(avg_expr[g, , drop = FALSE])
})

best_region  <- apply(region_score_mat, 1, function(x) names(which.max(x)))
second_region <- apply(region_score_mat, 1, function(x) names(sort(x, decreasing = TRUE)[2]))
top_region_score <- apply(region_score_mat, 1, max)

region_annotation <- data.frame(
  cluster       = cluster_ids,
  region        = best_region,
  second_choice = second_region,
  top_score     = round(top_region_score, 4),
  row.names     = NULL
)

write.csv(
  data.frame(x = best_region, row.names = cluster_ids),
  file = "Regions/mousebrain_region.csv"
)

write.csv(region_annotation, file = "Regions/mousebrain_region_scores.csv", row.names = FALSE)
write.csv(round(region_score_mat, 4), file = "Regions/mousebrain_region_scorematrix.csv")

message("✓ Exported: mousebrain_region.csv / _scores.csv / _scorematrix.csv")

# ========================= by Interests =======================
dir.create("Gene_of_Interest", showWarnings = FALSE, recursive = TRUE)

interests <- c(
  "App",  # Alzheimer’s disease related gene
  "Trem2","Tyrobp",  # complement & DAM (disease-associated microglia)
  "Cd68", "Cx3cr1", "Itgam", # classic microglia markers
  # Pro-inflammatory cytokines
  "Ccl3", "Ccl4", "Ccl5",
  # Astrocytes (reactive)
  "Gfap", "Vim", "S100b",
  # Inflammasome / signaling
  "Casp1", "Ptgs2", "Stat3",
  # Anti-inflammatory / resolution
  "Tgfb1", "Arg1", "Mrc1"
)
DotPlot_interests <- DotPlot(mousebrain, features = interests) + RotatedAxis() +
  scale_color_gradientn(colors = palette_DotPlot)
DotPlot_interests
ggsave("Gene_of_Interest/DotPlot_interests.svg", width = 9, height = 6)
#FeaturePlot_interests <- FeaturePlot(mousebrain, features = interests, ncol = 3)
#ggsave("FeaturePlot_interests.jpg", width = 6, height = 3)
#RidgePlot(mousebrain, features = interests, ncol = 1)
DoHeatmap_interests <- DoHeatmap(subset(mousebrain, downsample = 100), features = interests, size = 3)
DoHeatmap_interests
ggsave("Gene_of_Interest/DoHeatmap_interests.svg", width = 9, height = 6)
#VlnPlot(mousebrain, features = interests)
#VlnPlot(mousebrain, features = interests, slot = "counts", log = TRUE)

# ============================================================

mousebrain.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap_top10_in_each_cluster <- DoHeatmap(mousebrain, features = top10$gene) + NoLegend()
DoHeatmap_top10_in_each_cluster
ggsave("DoHeatmap_top10_in_each_cluster.jpg", width = 12, height = 18)


cluster0.markers <- FindMarkers(mousebrain, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(x = cluster0.markers, file = "cluster0.markers.csv")


# ============== Gene expression (.csv) ==========================
dir.create("Gene_expression", showWarnings = FALSE, recursive = TRUE)

# ---- Input ----
# gene_list: genes of interest
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
  "Tmem119","Aldh1l1","C3","Serping1","Slc1a2"
)

# 
assay_use <- "RNA"
slot_use  <- "data"   # log-normalized counts

# ---- Check genes present in object ----
all_genes <- rownames(GetAssayData(mousebrain, assay = "RNA", layer = "data"))
gene_list <- unique(gene_list)
present_genes <- intersect(gene_list, all_genes)
missing_genes <- setdiff(gene_list, all_genes)

if (length(missing_genes) > 0) {
  message("Missing genes (not found in object): ", paste(missing_genes, collapse = ", "))
}

stopifnot(length(present_genes) > 0)

# =========================
# 1) Per-cluster averages 
# =========================

avg_list <- AverageExpression(
  mousebrain,
  assays   = assay_use,
  features = present_genes,
  slot     = slot_use
)

avg_mat <- avg_list[[assay_use]]   # rows=genes, cols=clusters (log-normalized means)

# wide format
head(avg_mat)
write.csv(avg_mat, "Gene_expression/genes_mean_logexpr_per_cluster_wide.csv", row.names = TRUE)

# long/tidy format
avg_long <- as.data.frame(avg_mat)
avg_long$gene <- rownames(avg_long)
avg_long <- avg_long[, c(ncol(avg_long), 1:(ncol(avg_long)-1))]

suppressPackageStartupMessages({
  library(tidyr); library(dplyr)
})

avg_long <- avg_long |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "cluster",
    values_to = "mean_logexpr"
  )

head(avg_long)
write.csv(avg_long, "Gene_expression/genes_mean_logexpr_per_cluster_long.csv", row.names = FALSE)

# =========================
# 2) Per-cell values
# =========================
mat <- GetAssayData(mousebrain, assay = "RNA", layer = "data")[present_genes, , drop = FALSE]
clu <- as.character(Idents(mousebrain))  # factor -> character
names(clu) <- colnames(mousebrain)

# wide format
mat_wide <- as.data.frame(Matrix::t(mat), check.names = FALSE)  # cells x genes
mat_wide$cell    <- rownames(mat_wide)
mat_wide$cluster <- unname(clu[mat_wide$cell])

mat_wide <- mat_wide[, c("cell", "cluster", setdiff(colnames(mat_wide), c("cell","cluster")))]
gene_cols <- setdiff(colnames(mat_wide), c("cell", "cluster"))
write.csv(mat_wide, "Gene_expression/genes_logexpr_per_cell_wide.csv", row.names = FALSE)

# long format
long_df <- mat_wide |>
  tidyr::pivot_longer(
    cols = all_of(gene_cols),
    names_to = "gene",
    values_to = "logexpr"
  )

write.csv(long_df, "Gene_expression/genes_logexpr_per_cell_long.csv", row.names = FALSE)

# =========================
# 3) per-gene files
# =========================

invisible(lapply(gene_cols, function(g) {
  out <- long_df |> dplyr::filter(gene == g)
  write.csv(out,
            file.path("Gene_expression", paste0(g, "_logexpr_per_cell.csv")),
            row.names = FALSE)
}))

if (length(missing_genes) > 0) {
  writeLines(missing_genes,
             file.path("Gene_expression", "missing_genes.txt"))
  message("⚠ Missing Genes Record: ", paste(missing_genes, collapse = ", "))
}

message("✓ Output: the expression tables of ", length(gene_cols), " gene's have exported to Gene_expression/")


# ==========Assigning cell type identity to clusters ===========

# Cluster Annotation (Celltypes)
celltype_list <- read.csv("Celltypes/mousebrain_celltype.csv", row.names = 1)
new.cluster.ids.celltype <- celltype_list$x
new.cluster.ids.celltype

names(new.cluster.ids.celltype) <- levels(mousebrain)
mousebrain_celltype <- RenameIdents(mousebrain, new.cluster.ids.celltype)

DimPlot_celltype_tsne <- DimPlot(mousebrain_celltype, reduction = "tsne")
DimPlot_celltype_tsne
ggsave("Celltypes/DimPlot_celltype_tsne.svg", width = 6, height = 6)

DimPlot_celltype_tsne_NoLegend <- DimPlot(mousebrain_celltype, reduction = "tsne") + NoLegend ()
DimPlot_celltype_tsne_NoLegend
ggsave("Celltypes/DimPlot_celltype_tsne_NoLegend.svg", width = 6, height = 6)

saveRDS(mousebrain_celltype, file = "Celltypes/mousebrain_celltype.rds")


# Cluster Annotation (Regions)
region_list <- read.csv("Regions/mousebrain_region.csv", row.names = 1)
new.cluster.ids.region <- region_list$x
new.cluster.ids.region

names(new.cluster.ids.region) <- levels(mousebrain)
mousebrain_region <- RenameIdents(mousebrain, new.cluster.ids.region)

DimPlot_region_tsne <- DimPlot(mousebrain_region, reduction = "tsne")
DimPlot_region_tsne
ggsave("Regions/DimPlot_region_tsne.svg", width = 6, height = 6)

DimPlot_region_tsne_NoLegend <- DimPlot_region_tsne + NoLegend()
DimPlot_region_tsne_NoLegend
ggsave("Regions/DimPlot_region_tsne_NoLegend.svg", width = 6, height = 6)

saveRDS(mousebrain_region, file = "Regions/mousebrain_region.rds")

# ========================================================
saveRDS(mousebrain, file = "mousebrain_final.rds")

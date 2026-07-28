# =============================================================
# BBB integrity re-analysis (snRNA-seq only)
#   Sample1 = AD model (untreated)   Sample2 = AD model + ALNP (treated)
#   Question: does repeated ALNP dosing compromise BBB integrity?
#   -> Read barrier-related transcriptional programs in the vascular
#      compartment (endothelial + pericytes) and astrocyte endfeet.
#
# Method precedents (snRNA/scRNA as a BBB readout):
#   Yang et al. Nature 2022 (human brain vascular atlas, AD)
#   Vanlandewijck et al. Nature 2018 (brain vasculature atlas / mural markers)
#   Ben-Zvi 2014 Nature; Andreone 2017 Neuron (Mfsd2a / transcytosis; Plvap = leak)
#   Nitta 2003 JCB (Cldn5 sets the BBB sieve)
#   Nation 2019 Nat Med; Montagne 2020 Nature (pericyte loss = early BBB breakdown)
#   Tirosh 2016 Science (module scoring, AddModuleScore)
#   Squair 2021 Nat Commun (pseudoreplication caveat for cell-level DE)
#
# CAVEATS (n=1 per group; snRNA under-captures endothelium; transcriptome != barrier function).
# =============================================================

suppressMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(svglite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
root <- normalizePath(file.path(script_dir, ".."))
setwd(script_dir)
outdir <- "BBB"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

palette_DotPlot <- c("#4589C6", "white", "#CE544D")

mb <- readRDS("mousebrain_integrated_final.rds")
DefaultAssay(mb) <- "RNA"
mb$sample <- factor(mb$sample, levels = c("Sample1", "Sample2"))  # S1=AD, S2=AD+ALNP

# -------------------- BBB gene panels --------------------
panels <- list(
  TightJunction   = c("Cldn5","Ocln","Tjp1","Jam2","Cgn","Esam"),          # maintained = good
  Transcytosis    = c("Mfsd2a","Plvap","Cav1"),                            # Mfsd2a high / Plvap low = intact
  Transporter     = c("Slc2a1","Abcb1a","Abcg2","Slco1a4","Mfsd2a"),        # maintained = functional barrier
  EndoActivation  = c("Vcam1","Icam1","Vwf","Angpt2"),                      # up = activated/leaky = bad
  Pericyte        = c("Pdgfrb","Rgs5","Kcnj8","Anpep","Vtn","Notch3"),      # loss = early breakdown
  AstroEndfoot    = c("Aqp4","Gja1")                                        # loss/depol = dysfunction
)
panels <- lapply(panels, function(g) intersect(g, rownames(mb)))
all_bbb <- unique(unlist(panels))
writeLines(setdiff(unique(unlist(list(
  c("Cldn5","Ocln","Tjp1","Jam2","Cgn","Esam","Mfsd2a","Plvap","Cav1","Slc2a1",
    "Abcb1a","Abcg2","Slco1a4","Vcam1","Icam1","Vwf","Angpt2","Pdgfrb","Rgs5",
    "Kcnj8","Anpep","Vtn","Notch3","Aqp4","Gja1")))), rownames(mb)),
  file.path(outdir, "missing_genes.txt"))

# ============================================================
# 1) Composition: vascular abundance, pericyte:EC ratio, immune infiltration proxy
# ============================================================
comp <- as.data.frame(table(mb$sample, mb$celltype))
colnames(comp) <- c("sample","celltype","n")
comp <- comp %>% group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(comp, file.path(outdir, "Composition_celltype_by_sample.csv"), row.names = FALSE)

vasc_summary <- comp %>%
  filter(celltype %in% c("Endothelial","Pericytes")) %>%
  select(sample, celltype, n, frac) %>%
  pivot_wider(names_from = celltype, values_from = c(n, frac)) %>%
  mutate(pericyte_to_EC_ratio = n_Pericytes / n_Endothelial,
         vascular_frac = frac_Endothelial + frac_Pericytes)

# Ptprc+ (CD45) fraction per sample = peripheral-immune / infiltration proxy
cnt <- GetAssayData(mb, assay = "RNA", layer = "counts")
ptprc_pos <- if ("Ptprc" %in% rownames(cnt)) cnt["Ptprc", ] > 0 else rep(NA, ncol(mb))
immune <- data.frame(sample = mb$sample, ptprc_pos = as.integer(ptprc_pos)) %>%
  group_by(sample) %>% summarise(n_cells = n(), n_Ptprc_pos = sum(ptprc_pos),
                                 frac_Ptprc_pos = mean(ptprc_pos), .groups = "drop")
vasc_summary <- left_join(vasc_summary, immune, by = "sample")
write.csv(vasc_summary, file.path(outdir, "Vascular_composition_summary.csv"), row.names = FALSE)
message("[1] composition written")

# ============================================================
# 2) Per-gene mean log-expression + % expressing, by celltype x sample
#    (effect direction is the primary readout; see caveats)
# ============================================================
dat <- GetAssayData(mb, assay = "RNA", layer = "data")
meta <- mb@meta.data
grp <- interaction(meta$celltype, meta$sample, sep = "__", drop = TRUE)
genes_use <- intersect(all_bbb, rownames(dat))

mean_expr <- t(sapply(genes_use, function(g) tapply(dat[g, ], grp, mean)))
pct_expr  <- t(sapply(genes_use, function(g) tapply(cnt[g, ], grp, function(x) mean(x > 0))))

long_from <- function(m, value) {
  as.data.frame(as.table(as.matrix(m))) %>%
    setNames(c("gene","group",value)) %>%
    separate(group, into = c("celltype","sample"), sep = "__")
}
expr_tbl <- left_join(long_from(round(mean_expr,4), "mean_logexpr"),
                       long_from(round(pct_expr,4),  "pct_expr"),
                       by = c("gene","celltype","sample"))
write.csv(expr_tbl, file.path(outdir, "BBB_meanlogexpr_pct_by_celltype_x_sample.csv"), row.names = FALSE)

# focused wide table for the vascular compartment (delta = S2 - S1)
vasc_wide <- expr_tbl %>%
  filter(celltype %in% c("Endothelial","Pericytes")) %>%
  select(gene, celltype, sample, mean_logexpr) %>%
  pivot_wider(names_from = sample, values_from = mean_logexpr) %>%
  mutate(delta_S2_minus_S1 = round(Sample2 - Sample1, 4)) %>%
  arrange(celltype, desc(abs(delta_S2_minus_S1)))
write.csv(vasc_wide, file.path(outdir, "BBB_vascular_delta_S2_vs_S1.csv"), row.names = FALSE)
message("[2] per-gene expression tables written")

# ============================================================
# 3) Module scores in endothelial cells (AddModuleScore; Tirosh 2016)
#    higher TightJunction/Transporter = better barrier
#    higher Transcytosis(leak: Plvap/Cav1)/EndoActivation = worse
# ============================================================
endo <- subset(mb, subset = celltype == "Endothelial")
leak_set <- intersect(c("Plvap","Cav1"), rownames(endo))          # Mfsd2a scored separately
mod_lists <- list(
  TightJunction  = panels$TightJunction,
  Transporter    = intersect(c("Slc2a1","Abcb1a","Abcg2","Slco1a4","Mfsd2a"), rownames(endo)),
  Leak           = leak_set,
  EndoActivation = panels$EndoActivation
)
mod_lists <- mod_lists[sapply(mod_lists, length) > 0]
endo <- AddModuleScore(endo, features = mod_lists, name = "MOD_", seed = 42)
score_cols <- paste0("MOD_", seq_along(mod_lists))
names(score_cols) <- names(mod_lists)

scores_df <- endo@meta.data[, c("sample", score_cols)]
colnames(scores_df) <- c("sample", names(mod_lists))
write.csv(cbind(cell = rownames(scores_df), scores_df),
          file.path(outdir, "Endothelial_module_scores_percell.csv"), row.names = FALSE)

score_summary <- lapply(names(mod_lists), function(m) {
  s1 <- scores_df[[m]][scores_df$sample == "Sample1"]
  s2 <- scores_df[[m]][scores_df$sample == "Sample2"]
  wt <- tryCatch(wilcox.test(s2, s1)$p.value, error = function(e) NA)
  data.frame(module = m, genes = paste(mod_lists[[m]], collapse = ","),
             mean_S1 = round(mean(s1),4), mean_S2 = round(mean(s2),4),
             delta_S2_minus_S1 = round(mean(s2) - mean(s1),4),
             wilcox_p_celllevel = signif(wt, 3))
}) %>% bind_rows()
write.csv(score_summary, file.path(outdir, "Endothelial_module_scores_summary.csv"), row.names = FALSE)

# violin of module scores
scores_long <- scores_df %>% pivot_longer(-sample, names_to = "module", values_to = "score")
p_mod <- ggplot(scores_long, aes(sample, score, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.6) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Sample1" = "#4589C6", "Sample2" = "#CE544D")) +
  labs(x = NULL, y = "Module score", title = "Endothelial BBB module scores (S1=AD, S2=AD+ALNP)") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(outdir, "Endothelial_module_scores.svg"), p_mod, width = 12, height = 4)
message("[3] endothelial module scores written")

# ============================================================
# 4) DotPlot of BBB panels on the vascular subset, split by sample
# ============================================================
vasc <- subset(mb, subset = celltype %in% c("Endothelial","Pericytes"))
vasc$celltype <- droplevels(factor(vasc$celltype))

dp <- tryCatch(
  DotPlot(vasc, features = panels, group.by = "celltype", split.by = "sample",
          cols = c("#4589C6", "#CE544D")) + RotatedAxis(),
  error = function(e) { message("DotPlot split failed: ", conditionMessage(e)); NULL })
if (!is.null(dp)) ggsave(file.path(outdir, "DotPlot_BBB_vascular_split_by_sample.svg"), dp, width = 20, height = 5)

dp2 <- DotPlot(vasc, features = all_bbb, group.by = "celltype", split.by = "sample",
               cols = c("#4589C6", "#CE544D")) + RotatedAxis()
ggsave(file.path(outdir, "DotPlot_BBB_allgenes_vascular_split.svg"), dp2, width = 16, height = 5)

# ============================================================
# 5) Key-gene violins (endothelial + pericyte), split by sample
# ============================================================
key_endo <- intersect(c("Cldn5","Ocln","Tjp1","Mfsd2a","Slc2a1","Plvap","Vcam1","Icam1"), rownames(endo))
p_vln_e <- VlnPlot(endo, features = key_endo, group.by = "sample", pt.size = 0,
                   cols = c("#4589C6","#CE544D"), ncol = 4)
ggsave(file.path(outdir, "VlnPlot_endothelial_key_BBB_split.svg"), p_vln_e, width = 14, height = 6)

peri <- subset(mb, subset = celltype == "Pericytes")
key_peri <- intersect(c("Pdgfrb","Rgs5","Kcnj8","Anpep","Vtn","Notch3"), rownames(peri))
p_vln_p <- VlnPlot(peri, features = key_peri, group.by = "sample", pt.size = 0,
                   cols = c("#4589C6","#CE544D"), ncol = 3)
ggsave(file.path(outdir, "VlnPlot_pericyte_key_BBB_split.svg"), p_vln_p, width = 12, height = 6)
message("[4-5] dotplots + violins written")

# ============================================================
# 6) Endothelial DE (S2 vs S1) restricted to BBB panel
#    cell-level Wilcoxon: report log2FC/pct primarily; p is pseudoreplicated (Squair 2021)
# ============================================================
Idents(endo) <- "sample"
deg <- tryCatch(
  FindMarkers(endo, ident.1 = "Sample2", ident.2 = "Sample1",
              features = all_bbb, logfc.threshold = 0, min.pct = 0),
  error = function(e) { message("FindMarkers failed: ", conditionMessage(e)); NULL })
if (!is.null(deg)) {
  deg$gene <- rownames(deg)
  deg <- deg %>% arrange(desc(avg_log2FC)) %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)
  write.csv(deg, file.path(outdir, "Endothelial_DE_S2_vs_S1_BBBpanel.csv"), row.names = FALSE)
}

message("\n==== DONE. Outputs in Integration/BBB/ ====")
print(vasc_summary)
print(score_summary)

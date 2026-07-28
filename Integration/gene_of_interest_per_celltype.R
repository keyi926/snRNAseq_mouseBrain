# =============================================================
# Per-gene expression: Sample1 vs Sample2, WITHIN each cell type
#   + two-sided Wilcoxon rank-sum test (Mann-Whitney U)
# Cell types covered: Excitatory neurons, Inhibitory neurons,
#                     Microglia, Astrocytes.
# Outputs go to Gene_of_Interest_byCellType/<celltype>/...
# =============================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(Matrix)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
root <- normalizePath(file.path(script_dir, ".."))
setwd(script_dir)

out_root <- "Gene_of_Interest_byCellType"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

mb <- readRDS("mousebrain_integrated_final.rds")

gene_list <- unique(c(
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
))

all_genes <- rownames(mb)
present   <- intersect(gene_list, all_genes)
missing   <- setdiff(gene_list, all_genes)
if (length(missing) > 0) {
  writeLines(missing, file.path(out_root, "missing_genes.txt"))
  message("Missing genes (not in object): ", paste(missing, collapse = ", "))
}

cell_types <- c("Excitatory neurons", "Inhibitory neurons",
                "Microglia", "Astrocytes")

# Pull the expression matrix + metadata once, reuse across cell types
mat  <- GetAssayData(mb, assay = "RNA", layer = "data")[present, , drop = FALSE]
meta <- data.frame(
  cell     = colnames(mb),
  sample   = as.character(mb$sample),
  celltype = as.character(mb$celltype),
  stringsAsFactors = FALSE,
  row.names = colnames(mb)
)

sig_palette <- c("Sample1" = "#4589C6", "Sample2" = "#CE544D")

# ---- helper: run one cell-type ----
analyze_celltype <- function(ct) {
  ct_safe <- gsub("[ /]", "_", ct)
  ct_dir  <- file.path(out_root, ct_safe)
  dir.create(file.path(ct_dir, "per_gene"), showWarnings = FALSE, recursive = TRUE)

  cells <- meta$cell[meta$celltype == ct & !is.na(meta$celltype)]
  sub_meta <- meta[cells, ]
  sub_mat  <- mat[, cells, drop = FALSE]
  n_s1 <- sum(sub_meta$sample == "Sample1")
  n_s2 <- sum(sub_meta$sample == "Sample2")
  if (min(n_s1, n_s2) < 10) {
    message("Skip ", ct, ": cell counts S1=", n_s1, " S2=", n_s2)
    return(NULL)
  }
  message("Analyzing ", ct, ": S1=", n_s1, "  S2=", n_s2)

  rows <- lapply(present, function(g) {
    x  <- as.numeric(sub_mat[g, ])
    s1 <- x[sub_meta$sample == "Sample1"]
    s2 <- x[sub_meta$sample == "Sample2"]

    wt <- tryCatch(
      wilcox.test(s2, s1, alternative = "two.sided", exact = FALSE),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        wilcox.test(s2, s1, alternative = "two.sided", exact = FALSE)
      )
    )
    data.frame(
      gene     = g,
      celltype = ct,
      n_S1     = length(s1),
      n_S2     = length(s2),
      mean_S1  = mean(s1),
      mean_S2  = mean(s2),
      sd_S1    = sd(s1),
      sd_S2    = sd(s2),
      sem_S1   = sd(s1) / sqrt(length(s1)),
      sem_S2   = sd(s2) / sqrt(length(s2)),
      delta    = mean(s2) - mean(s1),
      W_stat   = if (is.null(wt)) NA_real_ else unname(wt$statistic),
      p_value  = if (is.null(wt)) NA_real_ else wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  stats <- bind_rows(rows) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           sig   = case_when(
             is.na(p_adj)  ~ "",
             p_adj < 0.001 ~ "***",
             p_adj < 0.01  ~ "**",
             p_adj < 0.05  ~ "*",
             TRUE          ~ "ns"
           )) |>
    arrange(desc(abs(delta)))

  write.csv(stats, file.path(ct_dir, "wilcoxon_results.csv"),
            row.names = FALSE)
  message("  ", ct, ": significant (adj.P<0.05) = ",
          sum(stats$sig != "ns" & stats$sig != ""))

  # ---- plotting tables ----
  plot_df <- bind_rows(
    stats |> transmute(gene, sample = "Sample1", mean_expr = mean_S1, sem = sem_S1),
    stats |> transmute(gene, sample = "Sample2", mean_expr = mean_S2, sem = sem_S2)
  )
  gene_order <- stats$gene
  plot_df$gene <- factor(plot_df$gene, levels = gene_order)
  stats$gene   <- factor(stats$gene,   levels = gene_order)

  sig_layer <- stats |>
    mutate(label = ifelse(sig %in% c("*","**","***"), sig, ""),
           y_top = pmax(mean_S1 + sem_S1, mean_S2 + sem_S2))

  # ---- Top 20 first (so it lands even if per-gene loop hits trouble) ----
  top_df <- stats[seq_len(min(20, nrow(stats))), ]
  ggplot(top_df, aes(x = reorder(gene, delta), y = delta,
                     fill = sig != "ns")) +
    geom_col() + coord_flip() +
    geom_text(aes(label = sig), hjust = -0.2, size = 3.2) +
    scale_fill_manual(values = c(`TRUE` = "#CE544D", `FALSE` = "#bbbbbb"),
                      guide = "none") +
    labs(x = NULL, y = "Δ mean log-expression (Sample2 − Sample1)",
         title = sprintf("%s — top 20 genes by |Δ|", ct)) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(face = "italic"))
  ggsave(file.path(ct_dir, "Top20_delta.svg"),
         width = 7, height = 7, device = "svg")

  # ---- Multi-panel overview ----
  p_all <- ggplot(plot_df, aes(x = sample, y = mean_expr, fill = sample)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = pmax(0, mean_expr - sem),
                      ymax = mean_expr + sem),
                  width = 0.25, linewidth = 0.3) +
    geom_text(data = sig_layer,
              aes(x = 1.5, y = y_top * 1.15, label = label),
              inherit.aes = FALSE, size = 4) +
    facet_wrap(~ gene, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = sig_palette) +
    labs(x = NULL, y = "log-normalized expression",
         title = sprintf("Gene-of-interest in %s (n=S1:%d, S2:%d)",
                         ct, n_s1, n_s2),
         subtitle = "mean ± SEM, two-sided Wilcoxon rank-sum on per-cell expression") +
    theme_classic(base_size = 9) +
    theme(legend.position = "top",
          strip.text = element_text(face = "italic"),
          panel.spacing = unit(0.3, "lines"))
  n_rows <- ceiling(length(present) / 6)
  ggsave(file.path(ct_dir, "Barplot_all_genes.svg"),
         p_all, width = 14, height = max(6, n_rows * 2.2),
         device = "svg")

  # ---- Per-gene SVGs ----
  fails <- character()
  for (g in present) {
    ok <- tryCatch({
      pd <- plot_df |> filter(gene == g)
      s  <- stats   |> filter(gene == g)
      sub_txt <- sprintf("Wilcoxon W = %.3g,  p = %.3g  (adj. p = %.3g,  %s)",
                         s$W_stat, s$p_value, s$p_adj, s$sig)
      p <- ggplot(pd, aes(x = sample, y = mean_expr, fill = sample)) +
        geom_col(width = 0.5, colour = "black", linewidth = 0.3) +
        geom_errorbar(aes(ymin = pmax(0, mean_expr - sem),
                          ymax = mean_expr + sem),
                      width = 0.18, linewidth = 0.4) +
        scale_fill_manual(values = sig_palette) +
        labs(title = g, subtitle = sub_txt,
             x = NULL, y = "log-normalized expression") +
        theme_classic(base_size = 12) +
        theme(legend.position = "none",
              plot.title    = element_text(face = "italic", size = 14),
              plot.subtitle = element_text(size = 9, colour = "#555555"))
      ggsave(file.path(ct_dir, "per_gene", paste0(g, ".svg")),
             plot = p, width = 4, height = 4, device = "svg")
      TRUE
    }, error = function(e) {
      message("    plot failed for ", g, " in ", ct, ": ",
              conditionMessage(e)); FALSE
    })
    if (!isTRUE(ok)) fails <- c(fails, g)
  }
  if (length(fails) > 0) {
    writeLines(fails, file.path(ct_dir, "per_gene", "_failed.txt"))
  }

  stats
}

# ---- Run all four cell types ----
all_stats <- lapply(cell_types, analyze_celltype)
names(all_stats) <- cell_types
all_stats <- all_stats[!sapply(all_stats, is.null)]

# ---- Combined long table ----
combined <- bind_rows(all_stats)
write.csv(combined,
          file.path(out_root, "wilcoxon_all_celltypes_combined.csv"),
          row.names = FALSE)
message("Combined table: ", nrow(combined), " rows  (",
        length(unique(combined$gene)), " genes × ",
        length(unique(combined$celltype)), " cell types)")

# ---- Heatmap: gene × celltype, Δ with significance overlay ----
combined$celltype <- factor(combined$celltype, levels = cell_types)

# Order genes by max |delta| across the four cell types
gene_max_delta <- combined |>
  group_by(gene) |>
  summarise(max_abs = max(abs(delta), na.rm = TRUE)) |>
  arrange(desc(max_abs))
combined$gene <- factor(combined$gene, levels = rev(gene_max_delta$gene))

ggplot(combined, aes(x = celltype, y = gene, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(sig %in% c("*","**","***"), sig, "")),
            size = 3) +
  scale_fill_gradient2(low = "#4589C6", mid = "white", high = "#CE544D",
                       midpoint = 0, name = "Δ log-expr\n(S2 − S1)") +
  labs(x = NULL, y = NULL,
       title = "Δ mean expression by cell type (Sample2 − Sample1)",
       subtitle = "* / ** / *** = adj.P < 0.05 / 0.01 / 0.001  (Wilcoxon, BH)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "italic"),
        panel.grid  = element_blank())

ggsave(file.path(out_root, "Heatmap_delta_by_celltype.svg"),
       width = 6.5, height = max(8, 0.22 * length(unique(combined$gene))),
       device = "svg")

# ---- Significance heatmap: p_adj direction count per gene ----
sig_summary <- combined |>
  mutate(direction = case_when(sig == "ns" | sig == "" ~ "ns",
                               delta > 0  ~ "Up",
                               TRUE       ~ "Down")) |>
  count(gene, direction) |>
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) |>
  arrange(desc(Up + Down))
write.csv(sig_summary,
          file.path(out_root, "Significance_summary_per_gene.csv"),
          row.names = FALSE)

message("Done. Outputs under: ", out_root)

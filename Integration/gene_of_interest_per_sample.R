# =============================================================
# Per-gene expression bar charts (Sample1 vs Sample2)
#   + two-sided Wilcoxon rank-sum test (Mann-Whitney U)
# Genes from the per-sample script's gene_list.
# Loads the integrated rds and uses the log-normalized RNA layer.
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
dir.create("Gene_of_Interest",         showWarnings = FALSE, recursive = TRUE)
dir.create("Gene_of_Interest/per_gene", showWarnings = FALSE, recursive = TRUE)

mb <- readRDS("mousebrain_integrated_final.rds")

# ---- the gene_list from the per-sample scripts ----
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
  writeLines(missing, "Gene_of_Interest/missing_genes.txt")
  message("Missing genes (not in object): ", paste(missing, collapse = ", "))
}
message("Genes present: ", length(present), " / ", length(gene_list))

# ---- Long-format expression table ----
mat <- GetAssayData(mb, assay = "RNA", layer = "data")[present, , drop = FALSE]
expr_wide <- as.data.frame(as.matrix(Matrix::t(mat)), check.names = FALSE)
expr_wide$cell   <- rownames(expr_wide)
expr_wide$sample <- mb$sample[expr_wide$cell]

long <- expr_wide |>
  pivot_longer(cols = all_of(present),
               names_to = "gene", values_to = "logexpr")

# ---- Wilcoxon rank-sum test per gene (two-sided, Sample2 vs Sample1) ----
# exact = FALSE forces the normal approximation; with ~4 300 cells per
# sample the exact test is impractical and the approximation is reliable.
stats <- long |>
  group_by(gene) |>
  summarise(
    n_S1      = sum(sample == "Sample1"),
    n_S2      = sum(sample == "Sample2"),
    mean_S1   = mean(logexpr[sample == "Sample1"]),
    mean_S2   = mean(logexpr[sample == "Sample2"]),
    median_S1 = median(logexpr[sample == "Sample1"]),
    median_S2 = median(logexpr[sample == "Sample2"]),
    sd_S1     = sd(logexpr[sample == "Sample1"]),
    sd_S2     = sd(logexpr[sample == "Sample2"]),
    sem_S1    = sd_S1 / sqrt(n_S1),
    sem_S2    = sd_S2 / sqrt(n_S2),
    delta     = mean_S2 - mean_S1,
    wt        = list(tryCatch(
      wilcox.test(logexpr ~ sample,
                  alternative = "two.sided",
                  exact = FALSE),
      error = function(e) NULL)),
    .groups = "drop"
  ) |>
  mutate(
    W_stat  = sapply(wt, function(x) if (is.null(x)) NA else unname(x$statistic)),
    p_value = sapply(wt, function(x) if (is.null(x)) NA else x$p.value)
  ) |>
  select(-wt) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      is.na(p_adj)       ~ "",
      p_adj < 0.001      ~ "***",
      p_adj < 0.01       ~ "**",
      p_adj < 0.05       ~ "*",
      TRUE               ~ "ns"
    )
  ) |>
  arrange(desc(abs(delta)))

write.csv(stats, "Gene_of_Interest/gene_of_interest_wilcoxon_by_sample.csv",
          row.names = FALSE)
message("Wilcoxon table written: ", nrow(stats),
        " genes; significant (adj.P<0.05): ",
        sum(stats$sig != "ns" & stats$sig != ""))

# ---- Plot data: mean ± SEM ----
plot_df <- bind_rows(
  stats |> transmute(gene, sample = "Sample1", mean_expr = mean_S1, sem = sem_S1),
  stats |> transmute(gene, sample = "Sample2", mean_expr = mean_S2, sem = sem_S2)
)

# Order genes by |delta|, descending
gene_order <- stats$gene
plot_df$gene <- factor(plot_df$gene, levels = gene_order)
stats$gene   <- factor(stats$gene,   levels = gene_order)

sig_palette <- c("Sample1" = "#4589C6", "Sample2" = "#CE544D")

# ---- 1) Multi-panel overview (all genes, faceted) ----
sig_layer <- stats |>
  mutate(label = ifelse(sig %in% c("*","**","***"), sig, ""),
         y_top = pmax(mean_S1 + sem_S1, mean_S2 + sem_S2))

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
       title = "Gene-of-interest expression: Sample1 vs Sample2",
       subtitle = "mean ± SEM, two-sided Wilcoxon rank-sum test on per-cell expression; * adj.P<0.05, ** <0.01, *** <0.001") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top",
        strip.text = element_text(face = "italic"),
        panel.spacing = unit(0.3, "lines"))

n_rows <- ceiling(length(present) / 6)
ggsave("Gene_of_Interest/Barplot_all_genes_by_sample.svg", p_all,
       width = 14, height = max(6, n_rows * 2.2))

# ---- 2) Top-divergent genes summary (top 20 by |delta|), runs first
#         so it gets produced even if the per-gene loop hits issues.
top_n <- min(20, nrow(stats))
top_df <- stats[seq_len(top_n), ]

ggplot(top_df, aes(x = reorder(gene, delta), y = delta,
                   fill = sig != "ns")) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = sig), hjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = "#CE544D", `FALSE` = "#bbbbbb"),
                    guide = "none") +
  labs(x = NULL, y = "Δ mean log-expression (Sample2 − Sample1)",
       title = "Top 20 genes by absolute Δ-expression between samples") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic"))
ggsave("Gene_of_Interest/Top20_delta_genes.svg", width = 7, height = 7)

# ---- 3) Per-gene SVGs (each wrapped so one bad gene doesn't kill the loop)
fails <- character()
for (g in present) {
  ok <- tryCatch({
    pd <- plot_df  |> filter(gene == g)
    s  <- stats    |> filter(gene == g)
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
    ggsave(sprintf("Gene_of_Interest/per_gene/%s.svg", g),
           plot = p, width = 4, height = 4, device = "svg")
    TRUE
  }, error = function(e) {
    message("Plot failed for ", g, ": ", conditionMessage(e))
    FALSE
  })
  if (!isTRUE(ok)) fails <- c(fails, g)
}
if (length(fails) > 0) {
  writeLines(fails, "Gene_of_Interest/per_gene/_failed.txt")
  message("Per-gene plots failed: ", length(fails), " → see _failed.txt")
}

message("Done. Outputs under: Gene_of_Interest/")

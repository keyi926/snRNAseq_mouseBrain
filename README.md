# Mouse Brain snRNA-seq Analysis

A Seurat v5 workflow for mouse-brain single-nucleus RNA sequencing. The project includes per-sample quality control and annotation, two-sample CCA integration, cell-type-specific differential expression, GO enrichment, genes-of-interest analyses, and a transcriptional assessment of blood-brain barrier (BBB)-related programs.

> **Important limitation:** the project contains one 10x sample per condition (`n = 1`). All between-sample results are descriptive and hypothesis-generating. They do not support group-level statistical inference or establish treatment safety.

## 1. Study design

| Sample | Biological condition | Comparison role |
| --- | --- | --- |
| `Sample1` | Untreated AD model mouse | Reference |
| `Sample2` | AD model mouse after repeated ALNP treatment | Treated sample |

The primary comparison is:

```text
Sample2 (AD + ALNP) - Sample1 (untreated AD)
```

## 2. Workflow

```text
10x filtered_feature_bc_matrix.h5
        |
        +-- Per-sample analysis (Sample1 / Sample2)
        |   QC -> LogNormalize -> HVG -> ScaleData -> PCA
        |   -> clustering -> t-SNE -> markers -> cell-type/region annotation
        |
        +-- Two-sample integration (Integration)
            merge -> Seurat v5 CCAIntegration -> joint clustering/t-SNE
            -> cell-type DEG -> GO -> microglia -> genes of interest
            -> expression/QC export -> BBB transcriptional assessment
```

## 3. Repository layout

```text
snRNAseq_mouseBrain/
├── README.md
├── .gitignore
├── Sample1/
│   └── snRNAseq_mouseBrain_sample1_tsne.R
├── Sample2/
│   └── snRNAseq_mouseBrain_sample 2_tsne.R

└── Integration/
    ├── snRNAseq_mouseBrain_integration.R
    ├── snRNAseq_mouseBrain_integration_downstream.R
    ├── expression_qc_dotplot.R
    ├── gene_of_interest_per_sample.R
    ├── gene_of_interest_per_celltype.R
    └── snRNAseq_mouseBrain_BBB_integrity.R
```

Raw inputs and all generated analysis outputs are local artifacts and are excluded from Git. The repository contains source code, repository configuration, and this README only.

## 4. Main scripts

| Script | Purpose | Main input | Main output |
| --- | --- | --- | --- |
| `Sample1/snRNAseq_mouseBrain_sample1_tsne.R` | Sample1 QC, clustering, annotation, and expression export | Sample1 10x H5 | Sample1 result directories |
| `Sample2/snRNAseq_mouseBrain_sample 2_tsne.R` | Equivalent Sample2 workflow | Sample2 10x H5 | Sample2 result directories |
| `Integration/snRNAseq_mouseBrain_integration.R` | CCA integration and principal analyses from the two H5 files | Two 10x H5 files | Integrated objects, t-SNE, DEG, GO, and microglia results |
| `Integration/snRNAseq_mouseBrain_integration_downstream.R` | Resume downstream analyses from an integrated object | `mousebrain_integrated.rds` | DEG, GO, microglia, and composition |
| `Integration/expression_qc_dotplot.R` | Expression matrices, QC summaries, and dot plots | Final integrated object | CSV/XLSX tables and figures |
| `Integration/gene_of_interest_per_sample.R` | Sample-level genes-of-interest comparison | Final integrated object | Wilcoxon results and SVG figures |
| `Integration/gene_of_interest_per_celltype.R` | Cell-type-specific genes-of-interest comparison | Final integrated object | Cell-type statistics and SVG figures |
| `Integration/snRNAseq_mouseBrain_BBB_integrity.R` | BBB-related composition, module, and gene assessment | Final integrated object | `Integration/BBB/` |
| `Integration/BBB/export_genes_logexpr_per_cell.R` | Export per-cell source data for BBB violin plots | Final integrated object | Metadata-bearing CSV |

## 5. Software environment

The project has been run with:

- R 4.3.3
- Seurat 5.2.1

Core R packages:

```r
Seurat
dplyr
tidyr
tibble
ggplot2
patchwork
svglite
Matrix
openxlsx
```

GO enrichment additionally requires:

```r
clusterProfiler
org.Mm.eg.db
```

The repository does not currently contain an `renv.lock` file, so the complete package environment is not locked.

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

A positive `avg_log2FC` or an `Up` label therefore means higher expression in `Sample2`. A negative value or a `Down` label means lower expression in `Sample2`.

The samples differ in `Xist` expression, suggesting that sex may be a potential confounder. The Fig.3-style presentation was inspired by the 2025 *Nature* paper “Lithium deficiency and the onset of Alzheimer's disease,” but this dataset and study design must be interpreted independently.

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

Principal per-sample settings:

- `CreateSeuratObject(min.cells = 3, min.features = 200)`
- `LogNormalize(scale.factor = 10000)`
- `FindVariableFeatures(selection.method = "vst", nfeatures = 2000)`
- PCA and neighbor graph: `dims = 1:10`
- `FindClusters(resolution = 0.5)`
- t-SNE for two-dimensional visualization
- Cluster-level marker-panel scores for cell-type and brain-region annotation

The integration workflow uses Seurat v5 `IntegrateLayers(method = CCAIntegration)` and clusters cells in the `integrated.cca` reduction.

## 3. Repository layout

```text
snRNAseq_mouseBrain/
├── README.md
├── .gitignore
├── Sample1/
│   ├── filtered_feature_bc_matrix.h5
│   ├── snRNAseq_mouseBrain_sample1_tsne.R
│   ├── QC/
│   ├── Celltypes/
│   ├── Regions/
│   ├── Gene_of_Interest/
│   └── Gene_expression/
├── Sample2/
│   ├── filtered_feature_bc_matrix.h5
│   ├── snRNAseq_mouseBrain_sample 2_tsne.R
│   ├── QC/
│   ├── Celltypes/
│   ├── Regions/
│   ├── Gene_of_Interest/
│   └── Gene_expression/
└── Integration/
    ├── snRNAseq_mouseBrain_integration.R
    ├── snRNAseq_mouseBrain_integration_downstream.R
    ├── expression_qc_dotplot.R
    ├── gene_of_interest_per_sample.R
    ├── gene_of_interest_per_celltype.R
    ├── snRNAseq_mouseBrain_BBB_integrity.R
    ├── BBB/
    ├── Celltypes/
    ├── DEG/
    ├── DotPlot/
    ├── GO/
    ├── Gene_expression/
    ├── Gene_of_Interest/
    ├── Gene_of_Interest_byCellType/
    ├── Microglia/
    ├── QC/
    └── TSNE/
```

Serialized R objects (`.rds`) and generated HTML reports are local artifacts and are excluded from Git.

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


## 6. Integrated dataset summary

The integrated object contains 8,614 nuclei:

| Cell type | Sample1 | Sample2 |
| --- | ---: | ---: |
| Excitatory neurons | 1,786 | 2,414 |
| Inhibitory neurons | 1,542 | 773 |
| Astrocytes | 282 | 305 |
| Oligodendrocytes | 361 | 408 |
| OPC | 86 | 116 |
| Microglia | 84 | 113 |
| Endothelial | 57 | 93 |
| Pericytes | 46 | 45 |
| Choroid plexus | 44 | 27 |
| Ependymal | 21 | 11 |
| **Total** | **4,309** | **4,305** |

The neuronal subtype proportions differ substantially between the samples. With one sample per condition, this difference may reflect treatment, individual variation, sampling, or technical composition.

## 7. BBB integrity analysis

`Integration/snRNAseq_mouseBrain_BBB_integrity.R` evaluates endothelial, pericyte, and astrocyte-endfoot transcriptional features:

- Tight-junction and transporter modules are directionally higher in `Sample2`.
- Endothelial activation and leak/transcytosis modules are directionally lower in `Sample2`.
- Pericyte and `Ptprc`-positive fractions do not show an evident increase.

These observations provide descriptive evidence that no strong BBB transcriptional injury signal was detected. They do not establish physical barrier integrity and cannot replace tracer leakage, IgG/fibrinogen extravasation, TEER, or independent animal replication. Detailed values, references, and limitations are documented in [`Integration/BBB/BBB_summary.md`](Integration/BBB/BBB_summary.md).


## 8. Git data policy

The repository retains:

- `.R`R analysis scripts;
- selected `.csv` and `.xlsx` outputs;

The repository excludes:
- the 10x `filtered_feature_bc_matrix.h5` inputs of two samples;
- `.svg`, `.jpg` images
- `.rds` / `.RDS` serialized objects;
- `.RData`, `.Rhistory`, logs, temporary files, and local backups;

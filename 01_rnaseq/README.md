# Bulk RNA-seq Differential Expression & Pathway Pipeline

End-to-end R workflow for bulk RNA-seq QC, outlier diagnostics, `DESeq2` differential expression, pre-ranked `fgsea`, and `clusterProfiler` Gene Ontology enrichment.

## Methodological Highlights

1. **Pre-fit exploratory QC**: Variance-stabilising transformation (`vst`) followed by `prcomp` sample/loading biplots (PC1–PC3), scree plots, top-100 variable gene heatmaps, and sample-to-sample Pearson correlation hierarchical clustering to validate outliers before model fitting.
2. **All-tested-gene pre-ranked GSEA**: Ranks the complete tested transcriptome by `DESeq2` Wald statistic (`stat`, combining effect direction and variance) for `fgseaMultilevel`, avoiding the statistical distortion caused by ranking only pre-filtered DEGs.
3. **Species-configurable, redundancy-reduced GO BP enrichment**: Runs `clusterProfiler::enrichGO` with any Bioconductor `OrgDb` (`--orgdb org.Hs.eg.db` or `org.Mm.eg.db`) against the tested gene universe (not the entire genome), applies consistent $|\log_2\mathrm{FC}|$ thresholds across combined, upregulated, and downregulated gene sets, and collapses semantically redundant terms with `simplify(cutoff = 0.7)`.

## Usage

```bash
Rscript scripts/run_deseq2_pathway_pipeline.R \
  --counts path/to/counts.csv \
  --metadata path/to/metadata.csv \
  --target TREATED \
  --reference CONTROL \
  --orgdb org.Hs.eg.db \
  --outliers "SAMPLE_OUTLIER_1" \
  --pathways path/to/msigdb_hallmark.gmt \
  --outdir outputs/rnaseq
```

### Arguments

- `--counts`: Raw integer count matrix (CSV/TSV with gene symbols as row names, or `.rdata` containing `countdata`).
- `--metadata`: Sample metadata CSV with sample IDs in column 1 and experimental group column (`--group-col`, default `group`).
- `--target` / `--reference`: Contrast levels (`log2FC = target - reference`).
- `--orgdb`: Bioconductor `OrgDb` annotation package name (`org.Mm.eg.db` by default; e.g. `org.Hs.eg.db` for human).
- `--outliers`: Optional comma-separated sample IDs to validate via correlation clustering and drop prior to `DESeq()`.
- `--pathways`: Optional `.rdata` or `.gmt` gene-set file for `fgseaMultilevel`.
- `--extra-contrasts`: Optional semicolon-separated `target:reference` pairs (e.g. `MUT_A:WT;MUT_B:WT`) for multi-contrast Venn overlap analysis.
- `--lfc-thresh`: Absolute $\log_2\mathrm{FC}$ threshold for volcano labelling and GO enrichment (default `1.0`).
- `--padj-thresh`: Adjusted $p$-value threshold for significance (default `0.05`).

## Testing

```bash
Rscript tests/test_rnaseq_cli.R
```

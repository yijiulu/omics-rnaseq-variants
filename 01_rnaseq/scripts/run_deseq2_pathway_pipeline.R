#!/usr/bin/env Rscript
# End-to-end bulk RNA-seq pipeline:
#   VST -> PCA / scree / top-var heatmap -> correlation outlier QC ->
#   outlier-filtered DESeq2 -> volcano -> all-ranked Wald-stat fgsea ->
#   GO BP enrichment -> top-DEG heatmap -> optional multi-contrast Venn overlap.

parse_cli <- function(args) {
  opts <- list(
    counts          = NULL,
    metadata        = NULL,
    target          = NULL,
    reference       = NULL,
    group_col       = "group",
    orgdb           = "org.Mm.eg.db",
    outliers        = "",
    pathways        = NULL,
    extra_contrasts = "",
    outdir          = "outputs/rnaseq",
    lfc_thresh      = 1.0,
    padj_thresh     = 0.05,
    n_volcano_label = 12L,
    min_row_sum     = 10L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript run_deseq2_pathway_pipeline.R --counts <csv|rdata> --metadata <csv>",
        "--target <group> --reference <group> [--group-col group] [--orgdb org.Mm.eg.db|org.Hs.eg.db]",
        "[--outliers id1,id2] [--pathways <rdata|gmt>] [--extra-contrasts 'T2:R2;T3:R3'] [--outdir dir]\n"
      )
      quit(status = 0)
    }
    if (i == length(args)) stop(sprintf("Missing value for %s", key), call. = FALSE)
    val <- args[i + 1L]
    ckey <- gsub("-", "_", sub("^--", "", key))
    if (!ckey %in% names(opts)) stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    opts[[ckey]] <- val
    i <- i + 2L
  }
  for (req in c("counts", "metadata", "target", "reference")) {
    if (is.null(opts[[req]]) || !nzchar(opts[[req]])) {
      stop(sprintf("Required argument --%s is missing. Pass --help for usage.", req), call. = FALSE)
    }
  }
  opts$lfc_thresh      <- as.numeric(opts$lfc_thresh)
  opts$padj_thresh     <- as.numeric(opts$padj_thresh)
  opts$n_volcano_label <- as.integer(opts$n_volcano_label)
  opts$min_row_sum     <- as.integer(opts$min_row_sum)
  opts
}

load_count_matrix <- function(path) {
  if (!file.exists(path)) stop(sprintf("Count file not found: %s", path), call. = FALSE)
  if (grepl("\\.rdata$", path, ignore.case = TRUE) || !grepl("\\.(csv|tsv|txt)$", path, ignore.case = TRUE)) {
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    if ("countdata" %in% loaded) return(as.matrix(env$countdata))
    return(as.matrix(env[[loaded[1]]]))
  }
  sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
  df <- read.table(path, header = TRUE, sep = sep, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  as.matrix(df)
}

plot_heatmap_pdf <- function(mat, filename, col_side_colors, group_levels, group_cols, heat_cols) {
  mat <- as.matrix(mat)
  pdf(filename, width = 11.5, height = 11)
  on.exit(dev.off(), add = TRUE)
  par(oma = c(0, 0, 0.2, 0))
  gplots::heatmap.2(
    mat,
    col = heat_cols,
    trace = "none",
    srtCol = 45,
    margins = c(12, 9),
    cexCol = 0.85,
    cexRow = if (nrow(mat) > 60) 0.42 else 0.85,
    density.info = "none",
    ColSideColors = col_side_colors,
    key = TRUE,
    key.title = "Colour key",
    key.xlab = "Row Z-score",
    key.par = list(mar = c(3.2, 3.0, 2.0, 0.5), mgp = c(1.8, 0.6, 0), cex.axis = 0.85, cex.lab = 0.90),
    lmat = rbind(c(5, 4, 6), c(0, 1, 6), c(3, 2, 6)),
    lhei = c(1.7, 0.22, 7.4),
    lwid = c(1.65, 7.6, 2.35),
    main = "",
    scale = "row",
    extrafun = function() {
      par(mar = c(0, 0, 0, 0))
      plot.new()
      legend(
        x = 0.02, y = 0.50,
        legend = group_levels,
        fill = group_cols,
        border = NA, bty = "n", cex = 0.78, title = "Group",
        y.intersp = 1.02, xjust = 0, yjust = 0.5
      )
    }
  )
}

plot_pca_pair <- function(pca_obj, eig, point_cols, group_levels, group_cols, pcs, fname, title_suffix = "") {
  pdf(fname, width = 16, height = 7.5)
  on.exit(dev.off(), add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4, 2.4, 4))

  par(mar = c(5.5, 5.6, 3.8, 0.8), mgp = c(3.2, 0.9, 0), las = 1, cex.axis = 1.18, cex.lab = 1.30, cex.main = 1.30)
  plot(
    pca_obj$x[, pcs],
    type = "n",
    main = sprintf("Sample plot (PC%d-PC%d%s)", pcs[1], pcs[2], title_suffix),
    xlab = sprintf("PC%d (%.1f%% variance)", pcs[1], eig[pcs[1]]),
    ylab = sprintf("PC%d (%.1f%% variance)", pcs[2], eig[pcs[2]])
  )
  grid(col = "grey92", lty = 1, lwd = 0.6)
  abline(v = 0, h = 0, col = "grey75", lty = 2)
  points(pca_obj$x[, pcs], pch = 21, bg = point_cols, col = "grey25", cex = 2.35, lwd = 0.65)

  par(mar = c(4.5, 0, 3, 0))
  plot.new()
  legend(
    "center",
    legend = group_levels, pt.bg = group_cols, col = "grey25",
    pch = 21, pt.cex = 2.25, cex = 1.22, bty = "n", title = "Group", y.intersp = 1.18
  )

  par(mar = c(5.5, 5.6, 3.8, 0.8), mgp = c(3.2, 0.9, 0), las = 1, cex.axis = 1.18, cex.lab = 1.30, cex.main = 1.30)
  plot(
    pca_obj$rotation[, pcs],
    type = "n",
    main = sprintf("Gene loadings (PC%d-PC%d%s)", pcs[1], pcs[2], title_suffix),
    xlab = sprintf("PC%d", pcs[1]), ylab = sprintf("PC%d", pcs[2])
  )
  grid(col = "grey92", lty = 1, lwd = 0.6)
  abline(v = 0, h = 0, col = "grey75", lty = 2)
  points(pca_obj$rotation[, pcs], pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.25))
}

run_pca_suite <- function(vsd, group_col, group_levels, group_cols, heat_cols, outdir, tag = "") {
  suffix <- if (nzchar(tag)) paste0("_", tag) else ""
  title_suffix <- if (nzchar(tag)) paste0(", ", tag) else ""

  pca_theme <- ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 18, face = "bold"),
      axis.title      = ggplot2::element_text(size = 17),
      axis.text       = ggplot2::element_text(size = 14),
      legend.title    = ggplot2::element_text(size = 15),
      legend.text     = ggplot2::element_text(size = 13),
      legend.key.size = grid::unit(0.75, "cm")
    )

  p_pca <- DESeq2::plotPCA(vsd, intgroup = group_col, returnData = FALSE) +
    pca_theme +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))
  ggplot2::ggsave(file.path(outdir, sprintf("pca_plot%s.pdf", suffix)), p_pca, dpi = 300, width = 8, height = 6.5)

  expr_t <- t(SummarizedExperiment::assay(vsd))
  pca_obj <- stats::prcomp(expr_t, scale. = FALSE)
  eig <- (pca_obj$sdev)^2
  names(eig) <- seq_along(eig)
  eig <- 100 * eig / sum(eig)

  p_cols <- group_cols[as.integer(SummarizedExperiment::colData(vsd)[[group_col]])]
  for (pair in list(c(1, 2), c(1, 3), c(2, 3))) {
    if (ncol(pca_obj$x) >= max(pair)) {
      plot_pca_pair(
        pca_obj, eig, p_cols, group_levels, group_cols, pair,
        file.path(outdir, sprintf("PCA_%d%d%s.pdf", pair[1], pair[2], suffix)),
        title_suffix = title_suffix
      )
    }
  }

  pdf(file.path(outdir, sprintf("scree%s.pdf", suffix)), width = 9, height = 6)
  par(mar = c(5.0, 5.0, 3.2, 1), mgp = c(2.9, 0.85, 0), las = 1, cex.axis = 1.18, cex.lab = 1.25, cex.main = 1.25)
  barplot(eig, main = sprintf("Scree plot%s", title_suffix), xlab = "PC", ylab = "% explained variance", col = "steelblue", border = NA)
  dev.off()

  expr_mat <- SummarizedExperiment::assay(vsd)
  top_n <- min(100L, nrow(expr_mat))
  top_genes <- names(sort(apply(expr_mat, 1, stats::var), decreasing = TRUE))[seq_len(top_n)]
  plot_heatmap_pdf(
    expr_mat[top_genes, , drop = FALSE],
    filename = file.path(outdir, sprintf("top100_var_heatmap%s.pdf", suffix)),
    col_side_colors = p_cols,
    group_levels = group_levels,
    group_cols = group_cols,
    heat_cols = heat_cols
  )
}

extract_deg_symbols <- function(dds, group_col, target, reference, padj_thresh = 0.05) {
  res <- DESeq2::results(dds, contrast = c(group_col, target, reference))
  df <- as.data.frame(res)
  keep <- !is.na(df$log2FoldChange) & !is.na(df$padj) & df$padj < padj_thresh
  rownames(df)[keep]
}

load_orgdb <- function(pkg_name) {
  if (!is.character(pkg_name) || length(pkg_name) != 1L || !nzchar(pkg_name)) {
    stop("--orgdb must be a non-empty Bioconductor OrgDb package name (e.g. org.Hs.eg.db or org.Mm.eg.db).", call. = FALSE)
  }
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(sprintf("Annotation package '%s' is not installed. Install via BiocManager::install('%s').", pkg_name, pkg_name), call. = FALSE)
  }
  getExportedValue(pkg_name, pkg_name)
}

main <- function() {
  opts <- parse_cli(commandArgs(trailingOnly = TRUE))

  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(DESeq2)
    library(RColorBrewer)
    library(ggplot2)
    library(ggrepel)
    library(gplots)
    library(fgsea)
    library(VennDiagram)
    library(clusterProfiler)
    library(enrichplot)
    library(stringr)
  })
  org_db <- load_orgdb(opts$orgdb)

  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  counts_all <- load_count_matrix(opts$counts)
  metadata <- read.table(opts$metadata, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
  sample_ids <- as.character(metadata[[1]])
  counts <- counts_all[, sample_ids, drop = FALSE]
  counts <- counts[rowSums(counts) > opts$min_row_sum, , drop = FALSE]
  metadata[[opts$group_col]] <- factor(metadata[[opts$group_col]])

  design_formula <- stats::as.formula(paste("~", opts$group_col))
  dds_full <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design_formula)
  vsd_full <- vst(dds_full, blind = FALSE)

  group_levels <- levels(colData(vsd_full)[[opts$group_col]])
  n_groups <- length(group_levels)
  group_cols <- grDevices::colorRampPalette(brewer.pal(min(12, max(3, n_groups)), "Paired"))(n_groups)
  heat_cols <- grDevices::colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

  run_pca_suite(vsd_full, opts$group_col, group_levels, group_cols, heat_cols, opts$outdir, tag = "full")

  # Outlier QC via sample-to-sample Pearson correlation clustering
  outliers <- if (nzchar(opts$outliers)) trimws(strsplit(opts$outliers, ",", fixed = TRUE)[[1]]) else character(0)
  outliers <- outliers[nzchar(outliers)]

  sample_cor <- stats::cor(assay(vsd_full), method = "pearson")
  sample_hc <- stats::hclust(stats::as.dist(1 - sample_cor), method = "complete")
  outlier_side <- ifelse(colnames(sample_cor) %in% outliers, "firebrick3", "grey85")
  pdf(file.path(opts$outdir, "sample_correlation_outlier_QC.pdf"), width = 11.5, height = 11)
  heatmap.2(
    sample_cor,
    Rowv = stats::as.dendrogram(sample_hc),
    Colv = stats::as.dendrogram(sample_hc),
    dendrogram = "both", symm = TRUE, col = heat_cols,
    trace = "none", density.info = "none",
    ColSideColors = outlier_side, srtCol = 45, margins = c(12, 9),
    cexCol = 0.85, cexRow = 0.85, lhei = c(2, 8), keysize = 0.9, main = ""
  )
  dev.off()

  if (length(outliers) > 0) {
    missing_out <- setdiff(outliers, colnames(counts))
    if (length(missing_out) > 0) {
      stop("Outlier sample(s) not found in count matrix: ", paste(missing_out, collapse = ", "), call. = FALSE)
    }
    keep_samples <- setdiff(colnames(counts), outliers)
    meta_clean <- metadata[metadata[[1]] %in% keep_samples, , drop = FALSE]
    counts_clean <- counts[, keep_samples, drop = FALSE]
    dds_clean <- DESeqDataSetFromMatrix(countData = counts_clean, colData = meta_clean, design = design_formula)
    vsd_clean <- vst(dds_clean, blind = FALSE)
    run_pca_suite(vsd_clean, opts$group_col, group_levels, group_cols, heat_cols, opts$outdir, tag = "clean")
  } else {
    dds_clean <- dds_full
    vsd_clean <- vsd_full
  }

  # Fit DESeq2 model and extract primary contrast
  dds <- DESeq(dds_clean)
  res <- results(dds, contrast = c(opts$group_col, opts$target, opts$reference))
  res_all <- as.data.frame(res)
  res_all$symbol <- rownames(res_all)
  write.csv(res_all, file.path(opts$outdir, sprintf("%s_minus_%s_deseq2_all.csv", opts$target, opts$reference)), row.names = FALSE)

  deg <- res_all[!is.na(res_all$log2FoldChange) & !is.na(res_all$padj), , drop = FALSE]
  deg <- deg[order(deg$padj), , drop = FALSE]
  deg_sig <- deg[deg$padj < opts$padj_thresh, , drop = FALSE]
  write.csv(deg_sig, file.path(opts$outdir, sprintf("%s_minus_%s_DEGs.csv", opts$target, opts$reference)), row.names = FALSE)

  # Volcano plot with dynamic top-N up/down labelling
  deg$state <- ifelse(
    deg$padj < opts$padj_thresh & deg$log2FoldChange > opts$lfc_thresh, "Upregulated",
    ifelse(deg$padj < opts$padj_thresh & deg$log2FoldChange < -opts$lfc_thresh, "Downregulated", "Not significant")
  )
  up_idx <- head(which(deg$state == "Upregulated"), opts$n_volcano_label)
  dn_idx <- head(which(deg$state == "Downregulated"), opts$n_volcano_label)
  label_idx <- c(up_idx, dn_idx)
  deg$point_label <- ""
  deg$point_label[label_idx] <- deg$symbol[label_idx]

  volcano_plot <- ggplot(deg, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = state), alpha = 0.7, size = 1.6) +
    geom_text_repel(aes(label = point_label), max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4, size = 3, segment.alpha = 0.4) +
    scale_color_manual(values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4", "Not significant" = "grey75")) +
    labs(x = "log2(fold change)", y = "-log10(adjusted p-value)", color = NULL, title = sprintf("%s minus %s", opts$target, opts$reference)) +
    theme_minimal() +
    theme(legend.position = "bottom", panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.5))

  ggsave(file.path(opts$outdir, sprintf("%s_minus_%s_volcano.pdf", opts$target, opts$reference)), volcano_plot, dpi = 300, width = 10, height = 8)

  # Top volcano-labelled genes clustered heatmap
  if (length(label_idx) >= 2) {
    top_symbols <- deg$symbol[label_idx]
    writeLines(top_symbols, file.path(opts$outdir, "top_volcano_genes.txt"))
    expr_mat <- assay(vsd_clean)
    avail <- intersect(top_symbols, rownames(expr_mat))
    if (length(avail) >= 2) {
      col_group <- group_cols[as.integer(colData(vsd_clean)[[opts$group_col]])]
      plot_heatmap_pdf(
        expr_mat[avail, , drop = FALSE],
        filename = file.path(opts$outdir, sprintf("%s_minus_%s_heatmap_top_genes.pdf", opts$target, opts$reference)),
        col_side_colors = col_group,
        group_levels = group_levels,
        group_cols = group_cols,
        heat_cols = heat_cols
      )
    }
  }

  # All-tested-gene Wald-statistic pre-ranked fgsea
  if (!is.null(opts$pathways) && nzchar(opts$pathways) && file.exists(opts$pathways)) {
    if (grepl("\\.gmt$", opts$pathways, ignore.case = TRUE)) {
      pathways <- fgsea::gmtPathways(opts$pathways)
    } else {
      penv <- new.env(parent = emptyenv())
      pname <- load(opts$pathways, envir = penv)
      pathways <- penv[[pname[1]]]
    }

    res_all$entrez <- mapIds(org_db, keys = res_all$symbol, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    fg_df <- res_all[!is.na(res_all$entrez) & !is.na(res_all$stat), , drop = FALSE]
    fg_df <- fg_df[order(abs(fg_df$stat), decreasing = TRUE), , drop = FALSE]
    fg_df <- fg_df[!duplicated(fg_df$entrez), , drop = FALSE]

    ranks <- sort(stats::setNames(fg_df$stat, fg_df$entrez), decreasing = TRUE)
    set.seed(42)
    fgsea_res <- as.data.frame(fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500))
    fgsea_res <- fgsea_res[order(fgsea_res$padj), , drop = FALSE]

    fg_out <- fgsea_res
    if ("leadingEdge" %in% names(fg_out)) {
      fg_out$leadingEdge <- vapply(fg_out$leadingEdge, paste, collapse = "/", FUN.VALUE = character(1))
    }
    write.csv(fg_out, file.path(opts$outdir, sprintf("%s_minus_%s_fgsea_all_ranked.csv", opts$target, opts$reference)), row.names = FALSE)

    fg_plot <- fgsea_res[abs(fgsea_res$NES) > 0.5 & fgsea_res$padj < opts$padj_thresh, , drop = FALSE]
    if (nrow(fg_plot) > 0) {
      fg_plot <- head(fg_plot[order(abs(fg_plot$NES), decreasing = TRUE), , drop = FALSE], 30)
      fp <- ggplot(fg_plot, aes(x = stats::reorder(pathway, NES), y = NES, fill = NES > 0)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"), guide = "none") +
        labs(x = "Pathway", y = "Normalised Enrichment Score (NES)", title = sprintf("%s vs %s", opts$target, opts$reference)) +
        theme_minimal()
      ggsave(file.path(opts$outdir, sprintf("%s_minus_%s_fgsea.pdf", opts$target, opts$reference)), fp, dpi = 300, width = 12, height = 8)
    }
  }

  # GO BP enrichment (hypergeometric test against tested gene universe)
  gene_universe <- res_all$symbol
  run_go <- function(genes, suffix, label) {
    if (length(genes) < 5) return(invisible(NULL))
    ego <- enrichGO(
      gene = genes, universe = gene_universe, OrgDb = org_db,
      keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(invisible(NULL))
    ego_simp <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    dp <- dotplot(ego_simp, showCategory = 15) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
      ggtitle(sprintf("GO BP enrichment: %s (%s vs %s)", label, opts$target, opts$reference)) +
      theme(axis.text.y = element_text(size = 9, lineheight = 0.85), plot.title = element_text(size = 11))
    ggsave(file.path(opts$outdir, sprintf("%s_minus_%s_GO_%s.pdf", opts$target, opts$reference, suffix)), dp, dpi = 300, width = 10, height = 10)
    invisible(ego_simp)
  }

  run_go(deg_sig$symbol[abs(deg_sig$log2FoldChange) > opts$lfc_thresh], "all", sprintf("all DEGs (|log2FC| > %g)", opts$lfc_thresh))
  run_go(deg_sig$symbol[deg_sig$log2FoldChange > opts$lfc_thresh], "up", paste("up in", opts$target))
  run_go(deg_sig$symbol[deg_sig$log2FoldChange < -opts$lfc_thresh], "down", paste("down in", opts$target))

  # Optional multi-contrast Venn diagram
  if (nzchar(opts$extra_contrasts)) {
    venn_items <- list()
    venn_items[[sprintf("%s_vs_%s", opts$target, opts$reference)]] <- deg_sig$symbol
    for (spec in strsplit(opts$extra_contrasts, ";", fixed = TRUE)[[1]]) {
      parts <- trimws(strsplit(spec, ":", fixed = TRUE)[[1]])
      if (length(parts) == 2) {
        venn_items[[sprintf("%s_vs_%s", parts[1], parts[2])]] <- extract_deg_symbols(dds, opts$group_col, parts[1], parts[2], opts$padj_thresh)
      }
    }
    if (length(venn_items) >= 2 && length(venn_items) <= 5) {
      palette <- c("dodgerblue", "darkorange", "darkblue", "forestgreen", "purple")[seq_along(venn_items)]
      vd <- venn.diagram(x = venn_items, category.names = names(venn_items), filename = NULL, output = TRUE, lty = "blank", fill = palette, cex = 0.8, fontface = "bold", cat.cex = 0.8)
      pdf(file.path(opts$outdir, "multi_contrast_venn.pdf"))
      grid::grid.draw(vd)
      dev.off()
    }
  }

  cat(sprintf("[OK] RNA-seq pipeline outputs written to %s\n", opts$outdir))
}

if (sys.nframe() == 0) {
  main()
}

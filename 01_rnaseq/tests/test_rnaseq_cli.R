#!/usr/bin/env Rscript
# Lightweight base-R unit tests for RNA-seq CLI parsing, count-matrix loading, and OrgDb validation.

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
test_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
script_path <- normalizePath(file.path(test_dir, "..", "scripts", "run_deseq2_pathway_pipeline.R"))

env <- new.env(parent = globalenv())
source(script_path, local = env)

# 1. Default CLI parsing
opts <- env$parse_cli(c(
  "--counts", "counts.csv",
  "--metadata", "meta.csv",
  "--target", "TREATED",
  "--reference", "CONTROL"
))
stopifnot(
  identical(opts$counts, "counts.csv"),
  identical(opts$orgdb, "org.Mm.eg.db"),
  identical(opts$group_col, "group"),
  opts$lfc_thresh == 1.0,
  opts$padj_thresh == 0.05,
  opts$n_volcano_label == 12L
)

# 2. Custom OrgDb, outliers, thresholds, and extra contrasts
opts2 <- env$parse_cli(c(
  "--counts", "counts.tsv",
  "--metadata", "meta.csv",
  "--target", "MUT",
  "--reference", "WT",
  "--orgdb", "org.Hs.eg.db",
  "--outliers", "S1,S2",
  "--lfc-thresh", "1.5",
  "--padj-thresh", "0.01",
  "--extra-contrasts", "KO:WT;HET:WT"
))
stopifnot(
  identical(opts2$orgdb, "org.Hs.eg.db"),
  identical(opts2$outliers, "S1,S2"),
  opts2$lfc_thresh == 1.5,
  opts2$padj_thresh == 0.01,
  identical(opts2$extra_contrasts, "KO:WT;HET:WT")
)

# 3. Missing required arguments and unknown flags raise errors
err_missing <- tryCatch({
  env$parse_cli(c("--counts", "counts.csv", "--target", "A", "--reference", "B"))
  FALSE
}, error = function(e) grepl("--metadata", conditionMessage(e), fixed = TRUE))
stopifnot(isTRUE(err_missing))

err_unknown <- tryCatch({
  env$parse_cli(c("--counts", "a", "--metadata", "b", "--target", "c", "--reference", "d", "--bogus", "1"))
  FALSE
}, error = function(e) grepl("Unknown argument", conditionMessage(e), fixed = TRUE))
stopifnot(isTRUE(err_unknown))

# 4. CSV and TSV count matrix loading
tmp_csv <- tempfile(fileext = ".csv")
tmp_tsv <- tempfile(fileext = ".tsv")
on.exit(unlink(c(tmp_csv, tmp_tsv)), add = TRUE)

writeLines(c("gene,S1,S2", "TP53,120,45", "MDM2,300,80"), tmp_csv)
writeLines(c("gene\tS1\tS2", "TP53\t120\t45", "MDM2\t300\t80"), tmp_tsv)

mat_csv <- env$load_count_matrix(tmp_csv)
mat_tsv <- env$load_count_matrix(tmp_tsv)
stopifnot(
  is.matrix(mat_csv),
  identical(rownames(mat_csv), c("TP53", "MDM2")),
  identical(colnames(mat_csv), c("S1", "S2")),
  mat_csv["TP53", "S1"] == 120,
  identical(mat_csv, mat_tsv)
)

# 5. Missing OrgDb package produces informative error
err_orgdb <- tryCatch({
  env$load_orgdb("org.Nonexistent.eg.db")
  FALSE
}, error = function(e) grepl("BiocManager::install", conditionMessage(e), fixed = TRUE))
stopifnot(isTRUE(err_orgdb))

cat("[OK] All 5 RNA-seq CLI and loader unit tests passed.\n")

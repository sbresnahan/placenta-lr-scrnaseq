#!/usr/bin/env Rscript
# process_pilot.R — counts -> Seurat QC -> SingleR placental composition -> isoform table.
#
# Input : a count matrix directory. Accepts either
#           (a) 10x-style MEX (matrix.mtx[.gz], barcodes.tsv[.gz], features/genes.tsv[.gz]), or
#           (b) a plain genes x cells counts file (.rds / .csv / .tsv).
#         scnanoseq/oarfish counts are converted to MEX by bin/oarfish_to_mex.py.
# Output: an analysis directory with a Seurat object, QC table, SingleR labels,
#         cell-composition table, and an isoform-level table (if transcript counts given).
#
# Example:
#   Rscript R/process_pilot.R \
#     --matrix_dir results/gene_mex \
#     --tx_matrix  results/tx_mex \
#     --reference  references/placenta_atlas.rds \
#     --gtf        references/HPLRv2.1.annotated.PANTRY.gtf \
#     --outdir     results/r_analysis

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(Seurat)
  library(SingleR)
  library(SingleCellExperiment)
  library(data.table)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--matrix_dir", type = "character", help = "MEX dir or counts file (gene level)"),
  make_option("--tx_matrix",  type = "character", default = NULL, help = "Optional transcript-level MEX dir / counts file"),
  make_option("--reference",  type = "character", default = NULL, help = "placenta_atlas.rds (SingleCellExperiment) for SingleR"),
  make_option("--ref_labels", type = "character", default = NULL, help = "colData column in reference with cell-type labels (auto-detect if NULL)"),
  make_option("--gtf",        type = "character", default = NULL, help = "HPLRv2.1 GTF (for isoform biotype/flag annotation)"),
  make_option("--min_cells",  type = "integer", default = 3),
  # Long-read scRNA-seq detects far fewer genes/UMIs per cell than short-read
  # (each molecule = 1 read), so the short-read default of 200 would drop every cell.
  make_option("--min_features", type = "integer", default = 20,
              help = "min genes/cell (long-read: ~20; short-read: ~200)"),
  make_option("--max_mt",     type = "double", default = 20.0, help = "max percent.mt (HPLRv2.1 chrM)"),
  make_option("--outdir",     type = "character", default = "results/r_analysis")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## ---- helpers ---------------------------------------------------------------
read_mex <- function(dir) {
  # version-robust 10x MEX reader (Read10X works across Seurat versions)
  mat <- Read10X(data.dir = dir)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]  # multi-assay MEX
  mat
}
read_counts <- function(path) {
  if (dir.exists(path)) return(read_mex(path))   # MEX directory
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  m <- fread(path)
  rn <- m[[1]]; m <- as.matrix(m[, -1]); rownames(m) <- rn; as(m, "dgCMatrix")
}
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

## ---- load gene counts ------------------------------------------------------
message("Loading gene counts from: ", opt$matrix_dir)
counts <- read_counts(opt$matrix_dir)
if (any(grepl("^ENSG", rownames(counts)))) rownames(counts) <- strip_version(rownames(counts))

# guard: min.features is applied at object creation; if it would remove every cell
# (e.g. a short-read threshold on long-read data), warn and relax it.
n_pass <- sum(Matrix::colSums(counts > 0) >= opt$min_features)
if (n_pass == 0) {
  warning("min_features=", opt$min_features, " removes ALL cells (long-read data has few ",
          "genes/cell). Relaxing to min_features=20.")
  opt$min_features <- 20
}
seu <- CreateSeuratObject(counts = counts, min.cells = opt$min_cells,
                          min.features = opt$min_features, project = "pilot")

## ---- QC (percent.mt works because HPLRv2.1 adds chrM / MT- genes) ----------
mt_pattern <- "^MT-"
if (!any(grepl(mt_pattern, rownames(seu)))) {
  message("NOTE: no MT- genes found — is the matrix annotated with HPLRv2.1? percent.mt set to 0.")
}
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)

qc_before <- data.table(metric = c("cells", "median_genes", "median_UMI", "median_pct_mt"),
                        value = c(ncol(seu), median(seu$nFeature_RNA), median(seu$nCount_RNA),
                                  round(median(seu$percent.mt), 2)))
fwrite(qc_before, file.path(opt$outdir, "qc_summary_prefilter.csv"))

seu <- subset(seu, subset = nFeature_RNA >= opt$min_features & percent.mt <= opt$max_mt)

## ---- normalize / cluster / UMAP -------------------------------------------
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

## ---- SingleR against the placental atlas ----------------------------------
if (!is.null(opt$reference) && file.exists(opt$reference)) {
  message("Running SingleR against: ", opt$reference)
  ref <- readRDS(opt$reference)
  ref_label_col <- opt$ref_labels
  if (is.null(ref_label_col)) {
    # Deng-lab term atlas: annotate_general covers all cells (7 placental types).
    cand <- c("annotate_general", "cell_type", "celltype", "annotation_merged",
              "annotation", "cluster", "label", "CellType")
    cand <- cand[cand %in% colnames(colData(ref))]
    # prefer the candidate with the fewest NAs
    na_count <- sapply(cand, function(cc) sum(is.na(colData(ref)[[cc]])))
    ref_label_col <- cand[which.min(na_count)]
    if (length(ref_label_col) == 0) stop("Could not auto-detect reference label column; pass --ref_labels")
  }
  message("  reference label column: ", ref_label_col)
  # align gene identifiers (strip versions on reference too)
  if (any(grepl("^ENSG", rownames(ref)))) rownames(ref) <- strip_version(rownames(ref))
  # SingleR needs log-normalized values in the reference
  if (!"logcounts" %in% assayNames(ref)) {
    message("  computing logcounts for reference (scran/scater)")
    suppressPackageStartupMessages({library(scran); library(scater)})
    ref <- computeSumFactors(ref)
    ref <- logNormCounts(ref)
  }
  sce <- as.SingleCellExperiment(seu)
  sce <- logNormCounts(sce)   # ensure test data also has logcounts
  pred <- SingleR(test = sce, ref = ref, labels = colData(ref)[[ref_label_col]])
  seu$singleR_label <- pred$labels[match(colnames(seu), rownames(pred))]
  seu$singleR_pruned <- pred$pruned.labels[match(colnames(seu), rownames(pred))]
  fwrite(as.data.table(pred, keep.rownames = "cell"),
         file.path(opt$outdir, "singler_per_cell.csv"))
} else {
  message("No reference provided — skipping SingleR (marker-based fallback in report).")
  seu$singleR_label <- seu$seurat_clusters
}

## ---- composition table -----------------------------------------------------
comp <- as.data.table(seu@meta.data)[, .N, by = .(singleR_label)][order(-N)]
comp[, fraction := round(N / sum(N), 4)]
fwrite(comp, file.path(opt$outdir, "cell_composition.csv"))

## ---- isoform-level table (optional) ---------------------------------------
if (!is.null(opt$tx_matrix) && file.exists(opt$tx_matrix)) {
  message("Processing transcript-level matrix")
  tx <- read_counts(opt$tx_matrix)
  tx_dt <- data.table(transcript_id = rownames(tx), total_counts = Matrix::rowSums(tx),
                      cells_detected = Matrix::rowSums(tx > 0))
  tx_dt[, transcript_id_stripped := strip_version(transcript_id)]
  tx_dt[, is_novel_HPLRT := grepl("^HPLRT_", transcript_id)]
  # annotate from GTF (biotype + population flags) if provided
  if (!is.null(opt$gtf) && file.exists(opt$gtf)) {
    gtf <- fread(cmd = paste("grep -v '^##'", shQuote(opt$gtf)), sep = "\t", header = FALSE, quote = "")
    tx_lines <- gtf[V3 == "transcript"]
    tx_anno <- tx_lines[, .(
      transcript_id = sub('.*transcript_id "([^"]+)".*', "\\1", V9),
      gene_name     = sub('.*gene_name "([^"]+)".*', "\\1", V9),
      tx_biotype    = sub('.*transcript_biotype "([^"]+)".*', "\\1", V9),
      EAS_term      = sub('.*EAS_term "([^"]+)".*', "\\1", V9),
      EUR_preterm   = sub('.*EUR_preterm "([^"]+)".*', "\\1", V9),
      EUR_term      = sub('.*EUR_term "([^"]+)".*', "\\1", V9)
    )]
    tx_dt <- merge(tx_dt, tx_anno, by = "transcript_id", all.x = TRUE)
  }
  setorder(tx_dt, -total_counts)
  fwrite(tx_dt, file.path(opt$outdir, "isoform_table.csv"))
  message("  placental HPLRT isoforms quantified: ", sum(tx_dt$is_novel_HPLRT, na.rm = TRUE))
}

## ---- save ------------------------------------------------------------------
saveRDS(seu, file.path(opt$outdir, "seurat_object.rds"))
qc_after <- data.table(metric = c("cells_retained", "median_genes", "median_UMI", "median_pct_mt"),
                       value = c(ncol(seu), median(seu$nFeature_RNA), median(seu$nCount_RNA),
                                 round(median(seu$percent.mt), 2)))
fwrite(qc_after, file.path(opt$outdir, "qc_summary_postfilter.csv"))
message("Done. Outputs in: ", opt$outdir)

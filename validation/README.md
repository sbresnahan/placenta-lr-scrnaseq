# Validation

The pipeline was validated two ways. No public **placental long-read** scRNA-seq
dataset exists yet, so the long-read arm is validated on a real PromethION 10x
single-cell dataset (LongBench, lung cancer cell lines) and the R/QC/cell-composition
arm is validated on a real **term placenta** 10x dataset (GSE247036).

## Arm A — long-read processing (BLAZE → minimap2 → oarfish)

**Data**: LongBench `SC_ONT.fastq.gz` (GSE303762; full PromethION flow cell, 10x 3' v3).
Subsampled the first **500,000 reads** via S3 range request (`longbench_subset.sh`).

**Steps run** (the same chain scnanoseq orchestrates on mhgcp):

| Step | Tool | Result |
|---|---|---|
| Barcode detection | BLAZE2 (`--10x-kit-version 3v3 --expect-cells 5000`) | **5,030 cells**; **79.0%** of reads assigned to cells |
| Genome alignment | minimap2 `-ax splice:hq -uf` vs GRCh38 | **97.99%** primary mapping rate |
| Transcriptome build | gffread (HPLRv2.1 GTF + hg38) | 51,543 transcripts (19,739 HPLRT_) |
| Quantification | oarfish 0.10.3 `--single-cell` (transcriptome-aligned, CB-tagged, CB-sorted BAM) | **5,027 cells × 51,543 transcripts** → **14,034 genes** |

Notes:
- LongBench is 10x **3' v3** chemistry → validated with the `3v3` whitelist. The pilot
  uses GEM-X **3' v4** (`--barcode_format 10X_3v4`); identical read structure, different
  whitelist. The mechanics validated here transfer directly.
- oarfish single-cell requires a **transcriptome-aligned** BAM with a `CB:Z` tag on every
  record, **sorted/collated by CB** (`samtools sort -t CB`). See `docs/WORKFLOW.md`.
- The HPLRv2.1 chrM addition works end-to-end: mitochondrial genes (MT-CO2, MT-RNR2) are
  the top-expressed genes, so MT% QC is computable.
- On mhgcp these steps run inside scnanoseq's Singularity containers; here they were run
  directly (no container runtime in the validation sandbox).
- The LongBench data are lung-cancer cell lines, so SingleR-vs-placenta composition is
  not biologically meaningful — Arm A validates the **mechanics** (barcode detection,
  alignment, quantification, MEX conversion, R ingestion). Arm B validates composition.

## Arm B — R / QC / cell composition (Seurat + SingleR)

**Data**: GSE247036 `GSM7882491` term placenta villi snRNA-seq (10x 3' v3.1),
`filtered_feature_bc_matrix` (MEX) — matches the term gestation of the pilot samples.

**Run**:
```bash
Rscript R/process_pilot.R \
  --matrix_dir validation/term_late1_mex \
  --reference references/placenta_atlas.rds \
  --ref_labels annotate_general \
  --gtf references/HPLRv2.1.annotated.PANTRY.gtf \
  --outdir validation/r_analysis_term
```

**Results** (`r_analysis_term/`):
- 4,439 cells; median 1,675 genes & 3,935 UMIs/cell; median 0.27% mitochondrial.
- **SingleR composition** (vs Deng-lab term atlas, `annotate_general` labels):

| Cell type | Fraction |
|---|---|
| Syncytiotrophoblast (STB) | 0.904 |
| Endothelial cells | 0.070 |
| Stroma | 0.010 |
| Cytotrophoblast (CTB) | 0.006 |
| Perivascular | 0.006 |
| Hofbauer cells | 0.003 |
| EVT | 0.0007 |

STB dominance is expected for **term villi snRNA-seq** (syncytiotrophoblast nuclei are
the most abundant in term villi). All major placental lineages are recovered.

**Rendered report**: `pilot_report_validation.html` — the knitted Rmd showing the exact
output format (QC table + violins, UMAP by SingleR label, composition chart, isoform
section). Figures also in `figures/`.

## Reference

`references/placenta_atlas.rds` = Deng-lab integrated **term** placenta atlas
(Zenodo doi:10.5281/zenodo.11097690, `term_final.h5ad`, 23,378 cells) converted to
SingleCellExperiment. Label column: `annotate_general` (7 placental types, full coverage).

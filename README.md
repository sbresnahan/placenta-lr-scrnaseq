# placenta-lr-scrnaseq

Pipeline to process, QC, and explore placental **long-read single-cell RNA-seq** data
(10x GEM-X 3' v4 → Oxford Nanopore PromethION), with cell composition assessed against
an integrated placental cell reference. Built for the **Baylor College of Medicine
`mhgcp` SLURM cluster** and the Baylor PeriBank pilot.

Scope: **quantification, QC, and cell-type composition** — no novel transcript discovery.

- **Heavy arm** (FASTQ → counts): [nf-core/scnanoseq](https://nf-co.re/scnanoseq) via
  Nextflow + Singularity (BLAZE barcode detection → minimap2 → UMI-tools → **oarfish**
  transcript-level quantification against the HPLRv2.1 transcriptome).
- **Analysis arm** (counts → report): R script (`R/process_pilot.R`) + Rmd report
  (`reports/pilot_report.Rmd`) using Seurat + SingleR against the Deng-lab integrated
  term placenta atlas.
- **Reference**: `references/HPLRv2.1.annotated.PANTRY.gtf.zip` — the HPLRv2.0 placental
  long-read transcriptome **plus GENCODE v45 chrM** (for mitochondrial QC). See
  `references/README.md`.

## Quickstart (BCM mhgcp)

```bash
# 0. One-time setup ----------------------------------------------------------
git clone <this-repo> && cd placenta-lr-scrnaseq

# Reference genome (hg38) — skip if you already have it
bash references/fetch_hg38.sh            # downloads + indexes GRCh38 primary assembly

# 10x GEM-X 3' v4 whitelist (for barcode correction)
bash bin/fetch_whitelist.sh              # downloads 3M-3pgex-may-2023 whitelist

# 1. Long-read arm (SLURM) ---------------------------------------------------
# Edit samplesheet.csv:  sample,fastq,cell_count
sbatch run_scnanoseq.slurm               # submits Nextflow to the mhgcp partition
# (builds the HPLRv2.1 transcriptome FASTA for oarfish automatically)

# 2. Analysis arm (R) --------------------------------------------------------
# Create the conda env (conda-forge + bioconda ONLY — no Anaconda channels)
conda env create -f envs/r-report.yml
conda activate placenta-lr-report

# Convert oarfish transcript counts -> gene + transcript MEX
python bin/oarfish_to_mex.py \
  --oarfish_prefix results/oarfish/OUT \
  --gtf references/HPLRv2.1.annotated.PANTRY.gtf \
  --out_gene results/gene_mex \
  --out_tx   results/tx_mex

# Process counts -> Seurat QC -> SingleR composition -> isoform table
Rscript R/process_pilot.R \
  --matrix_dir results/gene_mex \
  --tx_matrix  results/tx_mex \
  --reference references/placenta_atlas.rds \
  --outdir results/r_analysis

# 3. Render the report -------------------------------------------------------
Rscript -e "rmarkdown::render('reports/pilot_report.Rmd', \
  params=list(analysis_dir='results/r_analysis'))"
```

## Layout

| Path | Purpose |
|---|---|
| `run_scnanoseq.slurm` | sbatch wrapper for the long-read arm |
| `conf/mhgcp.config` | SLURM executor settings for mhgcp |
| `conf/mhgcp_resources.config` | PromethION-scale process resources |
| `envs/r-report.yml` | conda env (conda-forge/bioconda only) |
| `references/` | HPLRv2.1 GTF zip + transcriptome builder + hg38/whitelist fetchers |
| `R/process_pilot.R` | counts → QC → SingleR composition → isoform table |
| `reports/pilot_report.Rmd` | parameterized report template |
| `bin/oarfish_to_mex.py` | oarfish transcript counts → gene/transcript MEX |
| `bin/ocm_demux.py` | OPTIONAL OCM 4-plex demux (post-pilot) |
| `validation/` | how the pipeline was validated + results |
| `docs/WORKFLOW.md` | end-to-end core-handoff → report |
| `docs/OCM_4PLEX.md` | extending to 4-plex demultiplexing |

## Important notes

- **conda licensing**: all conda specs here use `conda-forge` + `bioconda` only
  (`nodefaults`). Do **not** add the Anaconda `defaults` channel.
- **Never activate conda on the mhgcp head node** — use an interactive/SLURM job.
- The long-read arm runs in **Singularity containers**; no conda needed for it.
- The pilot is a **single sample** (no OCM demultiplexing). For a future 4-plex run see
  `docs/OCM_4PLEX.md` and `bin/ocm_demux.py`.

## Validation

The pipeline was validated two ways (see `validation/README.md`):
1. **Long-read arm** — a 500k-read subset of the LongBench PromethION 10x single-cell
   dataset processed end-to-end through BLAZE → minimap2 → oarfish.
2. **Placental composition arm** — a term placenta villi snRNA-seq dataset (GSE247036)
   processed through the R arm, recovering the expected syncytiotrophoblast-dominant
   composition against the term placenta atlas.

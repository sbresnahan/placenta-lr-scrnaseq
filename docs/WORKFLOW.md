# End-to-end workflow: core handoff → report

This documents the full path of the Baylor PeriBank placental pilot from the cores to
a rendered report.

## 1. Wet-lab (cores)

| Step | Owner |
|---|---|
| Cell prep + suspension | BCM Single Cell Genomics (SCG) Core |
| 10x GEM-X 3' v4 library prep (up to ~20k cells) | SCG Core |
| Full-length barcoded cDNA → ONT library (Ligation Kit V14) | GARP Core |
| PromethION sequencing (R10.4.1, super-accurate basecalling, min 200 bp) | GARP Core |

**Deliverable from GARP**: one basecalled FASTQ (`.fastq.gz`) for the pilot sample.

## 2. Long-read processing (mhgcp)

```bash
# place FASTQ in data/, edit samplesheet.csv (sample,fastq,cell_count)
sbatch run_scnanoseq.slurm
```

scnanoseq runs: BLAZE (barcode detection) → minimap2 (spliced alignment to hg38) →
barcode tagging/correction → UMI-tools dedup → **oarfish** (transcript-level
quantification against the HPLRv2.1 transcriptome) → preliminary Seurat QC → MultiQC.

Scope note: this pipeline is for **quantification, QC, and cell-type composition only** —
no novel transcript discovery. oarfish quantifies known transcripts from the HPLRv2.1
annotation (which already includes the placental HPLRT_ isoforms).

Key outputs under `results/`:
- `oarfish/` — per-cell transcript count matrix (MEX-style)
- `multiqc/` — run-level QC
- aligned, barcode-tagged BAM

## 3. Convert counts for R

oarfish emits a per-cell **transcript** count matrix (`<prefix>.count.mtx` +
`.barcodes.txt` + `.features.txt`). Aggregate transcripts to genes and write a 10x MEX
directory for Seurat:

```bash
python bin/oarfish_to_mex.py \
  --oarfish_prefix results/oarfish/OUT \
  --gtf references/HPLRv2.1.annotated.PANTRY.gtf \
  --out_gene results/gene_mex \
  --out_tx   results/tx_mex
```

If you run the components manually, note two oarfish single-cell requirements:

1. **Transcriptome alignment**: oarfish quantifies against a transcriptome FASTA, not
   the genome. Build it once with `references/build_transcriptome.sh` (gffread), then
   align with `minimap2 -ax map-ont -N 100`.
2. **CB tag + CB-sorted BAM**: every record needs a `CB:Z` tag, and the BAM must be
   sorted/collated by CB (`samtools sort -t CB`) so all records for a cell are adjacent.
   BLAZE's demultiplexed FASTQ encodes the barcode/UMI in the read name
   (`@<CB>_<UMI>#<uuid>_<strand>`); re-add tags with `bin/tag_bam_cb_ub.py` first.

## 4. R analysis

```bash
conda activate placenta-lr-report
Rscript R/process_pilot.R \
  --matrix_dir results/gene_mex \
  --tx_matrix  results/tx_mex \
  --reference  references/placenta_atlas.rds \
  --gtf        references/HPLRv2.1.annotated.PANTRY.gtf \
  --outdir     results/r_analysis
```

## 5. Report

```bash
Rscript -e "rmarkdown::render('reports/pilot_report.Rmd', \
  params=list(analysis_dir='results/r_analysis'))"
# -> reports/pilot_report.html
```

## Notes

- Run everything from a SLURM allocation / interactive node — never on the head node.
- The long-read arm uses Singularity containers; only the R arm needs the conda env.
- For a future 4-plex OCM run, see `docs/OCM_4PLEX.md`.

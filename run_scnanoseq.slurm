#!/bin/bash
#SBATCH --job-name=scnanoseq
#SBATCH --partition=mhgcp        # Nextflow driver node (long-lived)
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=scnanoseq_%j.log

# run_scnanoseq.slurm — submit the nf-core/scnanoseq long-read arm to mhgcp.
# Nextflow (this job) is the SLURM *driver*; it submits each process as its own
# SLURM job per conf/mhgcp.config. Edit the params below, then:  sbatch run_scnanoseq.slurm

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

# --- user-editable parameters ------------------------------------------------
SAMPLESHEET="samplesheet.csv"                       # sample,fastq,cell_count
OUTDIR="results"
GENOME_FA="references/hg38/GRCh38.primary_assembly.genome.fa"
GTF="references/HPLRv2.1.annotated.PANTRY.gtf"      # unzip the shipped .zip first
# oarfish quantifies against a TRANSCRIPTOME FASTA derived from the GTF.
TRANSCRIPTOME_FA="references/HPLRv2.1.transcriptome.fa"  # built by references/build_transcriptome.sh
BARCODE_FORMAT="10X_3v4"                            # GEM-X 3' v4 (pilot). Use 10X_3v3 for v3.1.
QUANTIFIER="oarfish"                                # transcript-level quantification (no novel discovery)
# -----------------------------------------------------------------------------

# Nextflow itself (Java) — load via module or a conda-forge nextflow install.
# Prefer a site module if present; otherwise:  conda create -c conda-forge -n nf nextflow
if ! command -v nextflow >/dev/null 2>&1; then
  echo "nextflow not on PATH; load a module or activate a conda-forge nextflow env" >&2
  exit 1
fi

# Unzip the reference GTF if only the zip is present
if [ ! -f "$GTF" ] && [ -f "${GTF}.zip" ]; then unzip -o "${GTF}.zip" -d references/; fi

# Build the transcriptome FASTA (oarfish quantifies against transcripts)
if [ ! -f "$TRANSCRIPTOME_FA" ]; then bash references/build_transcriptome.sh "$GTF" "$GENOME_FA" "$TRANSCRIPTOME_FA"; fi

nextflow run nf-core/scnanoseq \
  -r 1.2.0 \
  -profile singularity \
  -c conf/mhgcp.config \
  -c conf/mhgcp_resources.config \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  --genome_fasta "$GENOME_FA" \
  --gtf "$GTF" \
  --transcript_fasta "$TRANSCRIPTOME_FA" \
  --barcode_format "$BARCODE_FORMAT" \
  --dedup_tool umitools \
  --quantifier "$QUANTIFIER" \
  --split_amount 500000 \
  -resume

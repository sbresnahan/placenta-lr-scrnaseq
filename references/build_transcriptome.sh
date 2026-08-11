#!/usr/bin/env bash
# build_transcriptome.sh — build the HPLRv2.1 transcriptome FASTA for oarfish.
# oarfish quantifies long reads against a transcriptome (not the genome), so we
# extract transcript sequences from the GTF + genome with gffread.
#
# Usage: bash build_transcriptome.sh <gtf> <genome_fa> <out_fa>
set -euo pipefail
GTF="${1:?gtf}"; GENOME="${2:?genome fa}"; OUT="${3:?out fasta}"

if ! command -v gffread >/dev/null 2>&1; then
  echo "gffread not found; install with: conda install -c bioconda gffread" >&2
  exit 1
fi
echo "Building transcriptome: $OUT"
gffread -w "$OUT" -g "$GENOME" "$GTF"
echo "Transcripts: $(grep -c '^>' "$OUT")"

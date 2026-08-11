#!/usr/bin/env bash
# fetch_hg38.sh — download + index the GRCh38 primary assembly for minimap2/oarfish.
# Skip if you already have an hg38 FASTA on mhgcp; then point run_scnanoseq.slurm at it.
set -euo pipefail
OUT="${1:-references/hg38}"
mkdir -p "$OUT"
FA="$OUT/GRCh38.primary_assembly.genome.fa"

if [ ! -f "$FA" ]; then
  echo "Downloading GRCh38 primary assembly (GENCODE)"
  curl -sSL -o "$FA.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
  echo "Decompressing"
  gunzip "$FA.gz"
fi

echo "Indexing with samtools (if available)"
if command -v samtools >/dev/null 2>&1; then samtools faidx "$FA"; fi

echo "Done -> $FA"

#!/usr/bin/env bash
# gse247036_fetch.sh — fetch a term placenta villi snRNA-seq sample (10x 3' v3.1)
# from GSE247036 to validate the R/QC/cell-composition arm on real placental data.
#
# The term samples are GSM7882491-96 (late1-6). Each is a processed
# filtered_feature_bc_matrix (MEX) tarball inside the series RAW.tar.
#
# NOTE: GSE145036 was considered first but its SRA runs are single-end 96-bp
# bulk-like reads, not 10x paired scRNA-seq, so GSE247036 is used instead.
set -euo pipefail
OUT="${1:-validation/term_late1_mex}"
SAMPLE="${2:-GSM7882491}"   # term "late1" villi
mkdir -p "$OUT"

RAW_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247036/suppl/GSE247036_RAW.tar"
echo "Fetching $SAMPLE MEX from GSE247036 RAW.tar"
# stream the series tar and extract only the one sample's matrix tarball, then unpack it
curl -s "$RAW_URL" | tar -xO --wildcards "*${SAMPLE}*filtered_feature_bc_matrix.tar.gz" \
  | tar -xzf - -C "$OUT"
echo "Done -> $OUT"
ls "$OUT"

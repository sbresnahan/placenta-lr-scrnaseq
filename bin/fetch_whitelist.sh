#!/usr/bin/env bash
# fetch_whitelist.sh — download 10x cell-barcode whitelists for barcode correction.
# GEM-X 3' v4 (pilot): 3M-3pgex-may-2023   |   Next GEM 3' v3.1: 3M-february-2018
set -euo pipefail
OUT="${1:-references/whitelists}"
mkdir -p "$OUT"
BASE="https://gcp-public-data--broad-references.storage.googleapis.com/optimus_whitelists"

for wl in 3M-3pgex-may-2023 3M-february-2018; do
  echo "Downloading $wl"
  curl -sSL -o "$OUT/${wl}.txt" "${BASE}/${wl}.txt"
  echo "  $(wc -l < "$OUT/${wl}.txt") barcodes"
done
echo "Done -> $OUT"

#!/usr/bin/env python3
"""ocm_demux.py — OPTIONAL OCM 4-plex demultiplexing (not needed for the 1-sample pilot).

10x GEM-X Universal Multiplex (OCM) pools up to 4 samples per GEM well, each with its
own gel-bead barcode list (OB1-OB4). After scnanoseq/BLAZE has *corrected* cell
barcodes, assign each corrected barcode to an OB list and split reads per sample.

Officially OCM demux is done by `cellranger multi` (v9.0+) with `ocm_barcode_ids`.
This script is a lightweight alternative that works on the corrected-barcode table.

Usage:
    python ocm_demux.py \
        --barcodes corrected_barcodes.tsv \   # one corrected barcode per line (or BAM CB tags)
        --ob_lists OB1.txt OB2.txt OB3.txt OB4.txt \
        --out_prefix demux

Outputs one file per OB list with the barcodes assigned to it, plus an "unassigned"
file. Feed the per-OB barcode lists back into read-splitting (e.g. `samtools view -D CB:file`
or a FASTQ filter) to get per-sample reads.
"""
import argparse
import sys


def load_barcodes(path):
    with open(path) as fh:
        return {line.strip().split("\t")[0] for line in fh if line.strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--barcodes", required=True, help="Corrected barcodes (one per line)")
    ap.add_argument("--ob_lists", nargs="+", required=True, help="OB1..OB4 barcode list files")
    ap.add_argument("--out_prefix", default="demux")
    args = ap.parse_args()

    corrected = load_barcodes(args.barcodes)
    ob_sets = []
    for i, p in enumerate(args.ob_lists, 1):
        ob_sets.append((f"OB{i}", load_barcodes(p)))

    assigned = {name: set() for name, _ in ob_sets}
    unassigned = set()
    for bc in corrected:
        hits = [name for name, s in ob_sets if bc in s]
        if len(hits) == 1:
            assigned[hits[0]].add(bc)
        else:
            unassigned.add(bc)  # 0 or >1 (barcode collision across OB lists)

    for name, s in assigned.items():
        with open(f"{args.out_prefix}_{name}.txt", "w") as fh:
            fh.write("\n".join(sorted(s)) + "\n")
        print(f"{name}: {len(s)} barcodes", file=sys.stderr)
    with open(f"{args.out_prefix}_unassigned.txt", "w") as fh:
        fh.write("\n".join(sorted(unassigned)) + "\n")
    print(f"unassigned: {len(unassigned)} barcodes", file=sys.stderr)


if __name__ == "__main__":
    main()

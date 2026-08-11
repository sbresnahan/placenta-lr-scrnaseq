# OCM 4-plex demultiplexing (post-pilot)

The pilot is a **single sample**, so no demultiplexing is needed. When you move to the
4-plex OCM design (4 samples pooled per GEM well), use one of the routes below.

## Background

10x GEM-X Universal Multiplex (OCM) gives each of the 4 pooled samples its own
gel-bead barcode list — **OB1–OB4** (3' kit: Blue=OB1, Red=OB2, Yellow=OB3, Green=OB4),
each ~921,600 barcodes (~3.6M total). A read's *corrected* cell barcode identifies
which sample it came from.

## Route A — cellranger multi (official)

10x's supported path is `cellranger multi` (v9.0+) with `ocm_barcode_ids`. This is
short-read oriented, so for long-read data Route B is usually simpler.

## Route B — barcode-list assignment (this repo)

After scnanoseq/BLAZE has produced **corrected** cell barcodes:

```bash
# 1. obtain the OB1-OB4 barcode lists (from 10x) into references/ocm/
# 2. assign corrected barcodes to OB lists
python bin/ocm_demux.py \
  --barcodes results/blaze/corrected_barcodes.tsv \
  --ob_lists references/ocm/OB1.txt references/ocm/OB2.txt \
             references/ocm/OB3.txt references/ocm/OB4.txt \
  --out_prefix results/demux/demux

# 3. split the barcode-tagged BAM per sample by CB tag
samtools view -D CB:results/demux/demux_OB1.txt -b results/tagged.bam > results/OB1.bam
# ... repeat for OB2-OB4, then re-run quantification per sample
```

Each per-sample BAM is then quantified (oarfish) and run through the R arm separately.

## Caveats

- Barcode collisions across OB lists (a corrected barcode matching >1 list) are sent to
  `_unassigned.txt`; with ~3.6M distinct barcodes these are rare.
- Confirm the exact OB list files with 10x / the SCG core for your kit lot.

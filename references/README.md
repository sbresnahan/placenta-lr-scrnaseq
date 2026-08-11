# references/

## HPLRv2.1.annotated.PANTRY.gtf.zip  (shipped in this repo)

The HPLRv2.0 placental long-read transcriptome (PANTRY-normalized) **plus the
GENCODE v45 mitochondrial chromosome**. HPLRv2.0 contains chr1–22/X/Y only — no
chrM — so mitochondrial-percent QC was impossible. v2.1 appends the 37 GENCODE v45
chrM genes (13 MT protein-coding, 22 Mt_tRNA, 2 Mt_rRNA) with their transcripts,
exons, and CDS, converted to the HPLR attribute style. **chrM is for QC only**; it
does not change the placental isoform annotation.

- Genes: 13,918 (13,881 HPLRv2.0 + 37 chrM)
- Unzip before use:  `unzip HPLRv2.1.annotated.PANTRY.gtf.zip`

## hg38 genome FASTA

minimap2 (genome) + oarfish (transcriptome) need the genome FASTA. If you already have hg38 on mhgcp, point
`run_scnanoseq.slurm` at it. Otherwise:

```bash
bash references/fetch_hg38.sh        # -> references/hg38/GRCh38.primary_assembly.genome.fa
```

## 10x whitelist

GEM-X 3' v4 (pilot) uses the `3M-3pgex-may-2023` whitelist; v3.1 uses
`3M-february-2018`. scnanoseq/BLAZE bundles these, but to fetch directly:

```bash
bash bin/fetch_whitelist.sh          # -> references/whitelists/
```

## Placental cell reference (for SingleR)

`R/process_pilot.R` expects `references/placenta_atlas.rds` — the Deng-lab integrated
placenta atlas (first-trimester + term) converted to a SingleCellExperiment. Build it
from the Zenodo h5ad (CC-BY, doi:10.5281/zenodo.11097690):

```r
library(zellkonverter); library(SingleCellExperiment)
sce <- readH5AD("first_trimester_final.h5ad")   # and/or term_final.h5ad
saveRDS(sce, "references/placenta_atlas.rds")
```

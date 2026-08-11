#!/usr/bin/env python3
"""oarfish_to_mex.py — convert scnanoseq/oarfish single-cell output to 10x MEX.

oarfish single-cell mode emits a per-cell TRANSCRIPT count matrix:
    <prefix>.count.mtx        (cells x transcripts, MatrixMarket)
    <prefix>.barcodes.txt     (cell barcodes, one per row)
    <prefix>.features.txt     (transcript ids, one per column)

This script:
  1. reads the oarfish matrix,
  2. maps transcripts -> genes using the GTF (so HPLRT_ novel placental
     isoforms roll up to their parent gene),
  3. writes a gene-level MEX directory (for Seurat/SingleR) and, optionally,
     a transcript-level MEX directory (for isoform exploration).

Counts are EM estimates (non-integer); they are rounded for Seurat.

Usage:
    python oarfish_to_mex.py \
        --oarfish_prefix results/oarfish/OUT \
        --gtf references/HPLRv2.1.annotated.PANTRY.gtf \
        --out_gene results/gene_mex \
        --out_tx   results/tx_mex
"""
import argparse
import gzip
import os
import re

import numpy as np
import pandas as pd
from scipy import io, sparse


def read_tx2gene(gtf):
    """Parse transcript_id -> (gene_id, gene_name) from a GTF (HPLR attribute style)."""
    tx2gene = {}
    tx_re = re.compile(r'transcript_id "([^"]+)"')
    gid_re = re.compile(r'gene_id "([^"]+)"')
    gname_re = re.compile(r'gene_name "([^"]+)"')
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m_t = tx_re.search(line)
            if not m_t:
                continue
            m_gid = gid_re.search(line)
            m_gname = gname_re.search(line)
            gid = m_gid.group(1) if m_gid else None
            gname = m_gname.group(1) if m_gname else None
            tx2gene[m_t.group(1)] = (gid, gname)
    return tx2gene


def write_mex(mat, barcodes, features, out, feature_type="Gene Expression"):
    """mat: sparse (features x cells)."""
    os.makedirs(out, exist_ok=True)
    with gzip.open(os.path.join(out, "matrix.mtx.gz"), "wb") as fh:
        io.mmwrite(fh, sparse.csr_matrix(mat))
    with gzip.open(os.path.join(out, "barcodes.tsv.gz"), "wt") as fh:
        fh.write("\n".join(map(str, barcodes)) + "\n")
    with gzip.open(os.path.join(out, "features.tsv.gz"), "wt") as fh:
        for f in features:
            fh.write(f"{f}\t{f}\t{feature_type}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oarfish_prefix", required=True,
                    help="oarfish output prefix (files <prefix>.count.mtx/.barcodes.txt/.features.txt)")
    ap.add_argument("--gtf", required=True, help="HPLRv2.1 GTF for transcript->gene mapping")
    ap.add_argument("--out_gene", required=True, help="output gene-level MEX dir")
    ap.add_argument("--out_tx", default=None, help="optional transcript-level MEX dir")
    args = ap.parse_args()

    p = args.oarfish_prefix
    barcodes = pd.read_csv(f"{p}.barcodes.txt", header=None)[0].astype(str).tolist()
    features = pd.read_csv(f"{p}.features.txt", header=None)[0].astype(str).tolist()
    mat = io.mmread(f"{p}.count.mtx")  # cells x transcripts
    mat = sparse.csr_matrix(mat)
    print(f"oarfish matrix: {mat.shape[0]} cells x {mat.shape[1]} transcripts")

    # round EM counts to integers for Seurat
    mat.data = np.rint(mat.data)

    # transcript-level MEX (features x cells)
    if args.out_tx:
        write_mex(mat.T.tocsr(), barcodes, features, args.out_tx)
        print(f"Wrote transcript MEX -> {args.out_tx}")

    # gene-level aggregation. Use gene NAME (falls back to gene_id, then tx id) so the
    # matrix is directly compatible with Seurat MT% QC (^MT-) and SingleR symbol matching.
    tx2gene = read_tx2gene(args.gtf)
    gene_labels = []
    for tx in features:
        gid, gname = tx2gene.get(tx, (None, None))
        gene_labels.append((gname or gid or tx, gid))
    # A few symbols map to >1 gene_id (e.g. PAR genes on X/Y). Seurat requires unique
    # rownames and disallows some punctuation, so disambiguate duplicated symbols with a
    # safe "-<gene_id>" suffix (dash + alphanumeric only).
    from collections import Counter
    name_counts = Counter(lbl for lbl, _ in gene_labels)

    def make_key(lbl, gid):
        if name_counts[lbl] > 1 and gid:
            return f"{lbl}-{gid.split('.')[0]}"
        return lbl

    uniq_genes = []
    seen = set()
    for lbl, gid in gene_labels:
        key = make_key(lbl, gid)
        if key not in seen:
            seen.add(key)
            uniq_genes.append(key)
    g_index = {g: i for i, g in enumerate(uniq_genes)}
    col_map = np.array([g_index[make_key(lbl, gid)] for lbl, gid in gene_labels])

    # aggregate columns (transcripts) into genes: gene_mat = mat @ indicator
    ind = sparse.csr_matrix(
        (np.ones(len(features)), (np.arange(len(features)), col_map)),
        shape=(len(features), len(uniq_genes)),
    )
    gene_mat = (mat @ ind).tocsr()  # cells x genes
    write_mex(gene_mat.T.tocsr(), barcodes, uniq_genes, args.out_gene)
    print(f"Wrote gene MEX -> {args.out_gene}  ({len(uniq_genes)} genes x {len(barcodes)} cells)")


if __name__ == "__main__":
    main()

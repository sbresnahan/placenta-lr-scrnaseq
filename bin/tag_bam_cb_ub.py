#!/usr/bin/env python3
"""tag_bam_cb_ub.py — add CB/UB BAM tags parsed from BLAZE-style read headers.

BLAZE's demultiplexed FASTQ encodes the cell barcode and UMI in the read name:
    @<CB>_<UMI>#<read_uuid>_<strand>
and minimap2 preserves that full name as the BAM QNAME but drops the CB/UB tags.
oarfish single-cell mode needs a CB:Z tag on every record (and the BAM must then be
sorted/collated by CB). This script streams a BAM, parses QNAME, and appends
CB:Z: and UB:Z: tags.

Usage:
    python tag_bam_cb_ub.py in.bam out.bam
"""
import re
import sys

import pysam

# QNAME: <CB>_<UMI>#<uuid>_<strand>
PAT = re.compile(r"^([ACGTN]+)_([ACGTN]+)#")


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    n_tagged = n_untagged = 0
    with pysam.AlignmentFile(inp, "rb") as ibam, \
         pysam.AlignmentFile(outp, "wb", header=ibam.header) as obam:
        for read in ibam:
            m = PAT.match(read.query_name)
            if m:
                read.set_tag("CB", m.group(1), value_type="Z")
                read.set_tag("UB", m.group(2), value_type="Z")
                n_tagged += 1
            else:
                n_untagged += 1
            obam.write(read)
    print(f"tagged={n_tagged} untagged={n_untagged}", file=sys.stderr)


if __name__ == "__main__":
    main()

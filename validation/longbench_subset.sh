#!/usr/bin/env bash
# longbench_subset.sh — subsample the LongBench PromethION 10x single-cell FASTQ.
# LongBench SC_ONT is a full 88 GB flow cell; we stream a byte-range and take the
# first ~N reads for pipeline validation. S3 is public (no-sign-request).
set -euo pipefail
N_READS="${1:-2000000}"          # default ~2M reads
OUT="${2:-validation/longbench_subset.fastq.gz}"
URL="https://longbench-data.s3.ap-southeast-2.amazonaws.com/raw/fastq/SC_ONT.fastq.gz"

# Each ONT read ~1-2 kb -> 4 lines; grab a generous byte range then truncate to N reads.
BYTES=$(( N_READS * 1200 ))
echo "Streaming first ~${BYTES} bytes of SC_ONT.fastq.gz"
curl -sSL -r 0-${BYTES} "$URL" | gzip -dc 2>/dev/null | head -n $(( N_READS * 4 )) | gzip -c > "$OUT"
echo "Reads written: $(( $(gzip -dc "$OUT" | wc -l) / 4 )) -> $OUT"

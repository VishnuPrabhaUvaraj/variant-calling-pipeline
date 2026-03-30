#!/bin/bash
set -euo pipefail

ALIGN_DIR="results/alignment"
SAMPLE="NA12878"
THREADS=4

echo "[BAM] Marking duplicates..."
picard MarkDuplicates \
  INPUT="$ALIGN_DIR/${SAMPLE}.sorted.bam" \
  OUTPUT="$ALIGN_DIR/${SAMPLE}.markdup.bam" \
  METRICS_FILE="$ALIGN_DIR/${SAMPLE}.dup_metrics.txt" \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true \
  2> logs/picard_markdup.log

echo "[BAM] Calculating coverage depth..."
samtools depth \
  -a "$ALIGN_DIR/${SAMPLE}.markdup.bam" \
  > "$ALIGN_DIR/${SAMPLE}.coverage.txt"

awk '{sum+=$3; count++} END {printf "Mean coverage depth: %.2fx\n", sum/count}' \
  "$ALIGN_DIR/${SAMPLE}.coverage.txt" \
  | tee "$ALIGN_DIR/${SAMPLE}.coverage_summary.txt"

echo "[BAM] Collecting insert size metrics..."
picard CollectInsertSizeMetrics \
  INPUT="$ALIGN_DIR/${SAMPLE}.markdup.bam" \
  OUTPUT="$ALIGN_DIR/${SAMPLE}.insert_metrics.txt" \
  HISTOGRAM_FILE="$ALIGN_DIR/${SAMPLE}.insert_histogram.pdf" \
  2> logs/picard_insertsize.log

echo "[BAM] Done."
cat "$ALIGN_DIR/${SAMPLE}.dup_metrics.txt" | grep -A 2 "ESTIMATED_LIBRARY_SIZE"

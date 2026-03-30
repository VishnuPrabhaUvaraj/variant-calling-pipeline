#!/bin/bash
set -euo pipefail

REF="data/reference/chr20.fa"
R1="data/trimmed/NA12878_R1_trimmed.fastq.gz"
R2="data/trimmed/NA12878_R2_trimmed.fastq.gz"
ALIGN_DIR="results/alignment"
THREADS=4
SAMPLE="NA12878"

mkdir -p "$ALIGN_DIR"

echo "[ALIGN] Running BWA-MEM..."
bwa mem \
  -t "$THREADS" \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  "$REF" "$R1" "$R2" \
  2> logs/bwa.log \
| samtools view -bS - \
| samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${SAMPLE}.sorted.bam"

echo "[ALIGN] Indexing BAM..."
samtools index "$ALIGN_DIR/${SAMPLE}.sorted.bam"

echo "[ALIGN] Generating alignment stats..."
samtools flagstat \
  "$ALIGN_DIR/${SAMPLE}.sorted.bam" \
  > "$ALIGN_DIR/${SAMPLE}.flagstat.txt"

cat "$ALIGN_DIR/${SAMPLE}.flagstat.txt"
echo "[ALIGN] Done."

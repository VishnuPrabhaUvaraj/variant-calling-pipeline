#!/bin/bash
set -euo pipefail

RAW_R1="data/raw/NA12878_R1.fastq.gz"
RAW_R2="data/raw/NA12878_R2.fastq.gz"
TRIMMED_DIR="data/trimmed"
QC_DIR="results/qc"
THREADS=4
ADAPTER="ADAPTER_PATH_HERE"

mkdir -p "$TRIMMED_DIR" "$QC_DIR/raw" "$QC_DIR/trimmed"

echo "[QC] Running FastQC on raw reads..."
fastqc "$RAW_R1" "$RAW_R2" --outdir "$QC_DIR/raw" --threads "$THREADS"

echo "[QC] Running Trimmomatic..."
trimmomatic PE -threads "$THREADS" \
  "$RAW_R1" "$RAW_R2" \
  "$TRIMMED_DIR/NA12878_R1_trimmed.fastq.gz" \
  "$TRIMMED_DIR/NA12878_R1_unpaired.fastq.gz" \
  "$TRIMMED_DIR/NA12878_R2_trimmed.fastq.gz" \
  "$TRIMMED_DIR/NA12878_R2_unpaired.fastq.gz" \
  ILLUMINACLIP:"$ADAPTER":2:30:10:2:keepBothReads \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
  2> logs/trimmomatic.log

echo "[QC] Running FastQC on trimmed reads..."
fastqc "$TRIMMED_DIR/NA12878_R1_trimmed.fastq.gz" \
       "$TRIMMED_DIR/NA12878_R2_trimmed.fastq.gz" \
       --outdir "$QC_DIR/trimmed" --threads "$THREADS"

echo "[QC] Generating MultiQC report..."
multiqc "$QC_DIR" --outdir "$QC_DIR" --filename multiqc_report

echo "[QC] Done!"

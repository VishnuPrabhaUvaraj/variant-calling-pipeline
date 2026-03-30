#!/bin/bash
set -euo pipefail

VAR_DIR="results/variants"
ANNOT_DIR="results/annotation"
SAMPLE="NA12878"

mkdir -p "$ANNOT_DIR"

# Remove empty output file from previous failed run
rm -f "$ANNOT_DIR/${SAMPLE}.annotated.vcf"

echo "[ANNOT] Running SnpEff annotation..."
snpEff ann \
  -Xmx16g \
  -v \
  -stats "$ANNOT_DIR/${SAMPLE}.snpeff_stats.html" \
  GRCh38.99 \
  "$VAR_DIR/${SAMPLE}.final.vcf.gz" \
  > "$ANNOT_DIR/${SAMPLE}.annotated.vcf"

echo "[ANNOT] Compressing and indexing..."
bgzip "$ANNOT_DIR/${SAMPLE}.annotated.vcf"
tabix -p vcf "$ANNOT_DIR/${SAMPLE}.annotated.vcf.gz"

echo "[ANNOT] Extracting variant table..."
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/ANN\n' \
  "$ANNOT_DIR/${SAMPLE}.annotated.vcf.gz" \
  > "$ANNOT_DIR/${SAMPLE}.variant_table.tsv"

echo "[ANNOT] Total annotated variants:"
bcftools view -H "$ANNOT_DIR/${SAMPLE}.annotated.vcf.gz" | wc -l

echo "[ANNOT] Done!"

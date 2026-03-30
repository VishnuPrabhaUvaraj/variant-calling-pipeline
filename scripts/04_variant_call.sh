#!/bin/bash
set -euo pipefail

REF="data/reference/chr20.fa"
ALIGN_DIR="results/alignment"
VAR_DIR="results/variants"
SAMPLE="NA12878"
THREADS=4

mkdir -p "$VAR_DIR"

# ── Step 1: HaplotypeCaller ─────────────────────────
echo "[GATK] Running HaplotypeCaller..."
gatk HaplotypeCaller \
  -I "$ALIGN_DIR/${SAMPLE}.markdup.bam" \
  -R "$REF" \
  -O "$VAR_DIR/${SAMPLE}.raw.vcf.gz" \
  -L chr20 \
  --native-pair-hmm-threads "$THREADS" \
  2> logs/haplotypecaller.log

echo "[GATK] HaplotypeCaller done."

# ── Step 2: Split SNPs and Indels ──────────────────
echo "[GATK] Separating SNPs and indels..."
gatk SelectVariants \
  -R "$REF" \
  -V "$VAR_DIR/${SAMPLE}.raw.vcf.gz" \
  --select-type-to-include SNP \
  -O "$VAR_DIR/${SAMPLE}.raw_snps.vcf.gz"

gatk SelectVariants \
  -R "$REF" \
  -V "$VAR_DIR/${SAMPLE}.raw.vcf.gz" \
  --select-type-to-include INDEL \
  -O "$VAR_DIR/${SAMPLE}.raw_indels.vcf.gz"

# ── Step 3: Hard filter SNPs ───────────────────────
echo "[GATK] Hard filtering SNPs..."
gatk VariantFiltration \
  -R "$REF" \
  -V "$VAR_DIR/${SAMPLE}.raw_snps.vcf.gz" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  -O "$VAR_DIR/${SAMPLE}.filtered_snps.vcf.gz" \
  2> logs/filter_snps.log

# ── Step 4: Hard filter Indels ─────────────────────
echo "[GATK] Hard filtering indels..."
gatk VariantFiltration \
  -R "$REF" \
  -V "$VAR_DIR/${SAMPLE}.raw_indels.vcf.gz" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  -O "$VAR_DIR/${SAMPLE}.filtered_indels.vcf.gz" \
  2> logs/filter_indels.log

# ── Step 5: Merge filtered VCFs ────────────────────
echo "[GATK] Merging filtered VCFs..."
gatk MergeVcfs \
  -I "$VAR_DIR/${SAMPLE}.filtered_snps.vcf.gz" \
  -I "$VAR_DIR/${SAMPLE}.filtered_indels.vcf.gz" \
  -O "$VAR_DIR/${SAMPLE}.final.vcf.gz"

# ── Step 6: Count PASS variants ────────────────────
echo "[GATK] PASS variant counts:"
zcat "$VAR_DIR/${SAMPLE}.final.vcf.gz" \
  | grep -v "^#" \
  | awk '$7=="PASS"' \
  | wc -l
echo "variants passed filters"

echo "[GATK] Done!"

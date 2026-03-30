configfile: "config/config.yaml"

SAMPLE  = config["sample"]
REF     = config["reference"]
THREADS = config["threads"]
ADAPTER = config["adapter"]

rule all:
    input:
        f"results/qc/raw/{SAMPLE}_R1_fastqc.html",
        f"results/alignment/{SAMPLE}.sorted.bam",
        f"results/alignment/{SAMPLE}.markdup.bam",
        f"results/variants/{SAMPLE}.final.vcf.gz",
        f"results/annotation/{SAMPLE}.annotated.vcf.gz",
        f"figures/variant_summary.png"

rule fastqc_raw:
    input:
        r1=f"data/raw/{SAMPLE}_R1.fastq.gz",
        r2=f"data/raw/{SAMPLE}_R2.fastq.gz"
    output:
        f"results/qc/raw/{SAMPLE}_R1_fastqc.html"
    threads: THREADS
    shell:
        "fastqc {input.r1} {input.r2} --outdir results/qc/raw -t {threads}"

rule trimmomatic:
    input:
        r1=f"data/raw/{SAMPLE}_R1.fastq.gz",
        r2=f"data/raw/{SAMPLE}_R2.fastq.gz"
    output:
        r1=f"data/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        r2=f"data/trimmed/{SAMPLE}_R2_trimmed.fastq.gz"
    threads: THREADS
    log: "logs/trimmomatic.log"
    shell:
        """
        trimmomatic PE -threads {threads} \
          {input.r1} {input.r2} \
          {output.r1} data/trimmed/{SAMPLE}_R1_unpaired.fastq.gz \
          {output.r2} data/trimmed/{SAMPLE}_R2_unpaired.fastq.gz \
          ILLUMINACLIP:{ADAPTER}:2:30:10:2:keepBothReads \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
          2> {log}
        """

rule bwa_align:
    input:
        r1=f"data/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        r2=f"data/trimmed/{SAMPLE}_R2_trimmed.fastq.gz",
        ref=REF
    output:
        bam=f"results/alignment/{SAMPLE}.sorted.bam",
        bai=f"results/alignment/{SAMPLE}.sorted.bam.bai"
    threads: THREADS
    log: "logs/bwa.log"
    shell:
        """
        bwa mem -t {threads} \
          -R "@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1" \
          {input.ref} {input.r1} {input.r2} 2>{log} \
          | samtools view -bS - \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule markduplicates:
    input:
        f"results/alignment/{SAMPLE}.sorted.bam"
    output:
        bam=f"results/alignment/{SAMPLE}.markdup.bam",
        metrics=f"results/alignment/{SAMPLE}.dup_metrics.txt"
    shell:
        """
        picard MarkDuplicates \
          INPUT={input} \
          OUTPUT={output.bam} \
          METRICS_FILE={output.metrics} \
          CREATE_INDEX=true \
          VALIDATION_STRINGENCY=LENIENT
        """

rule haplotypecaller:
    input:
        bam=f"results/alignment/{SAMPLE}.markdup.bam",
        ref=REF
    output:
        f"results/variants/{SAMPLE}.raw.vcf.gz"
    threads: THREADS
    shell:
        """
        gatk HaplotypeCaller \
          -I {input.bam} \
          -R {input.ref} \
          -O {output} \
          -L chr20 \
          --native-pair-hmm-threads {threads}
        """

rule filter_variants:
    input:
        vcf=f"results/variants/{SAMPLE}.raw.vcf.gz",
        ref=REF
    output:
        f"results/variants/{SAMPLE}.final.vcf.gz"
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} \
          --select-type-to-include SNP \
          -O results/variants/{SAMPLE}.raw_snps.vcf.gz

        gatk SelectVariants -R {input.ref} -V {input.vcf} \
          --select-type-to-include INDEL \
          -O results/variants/{SAMPLE}.raw_indels.vcf.gz

        gatk VariantFiltration -R {input.ref} \
          -V results/variants/{SAMPLE}.raw_snps.vcf.gz \
          --filter-expression "QD < 2.0" --filter-name "QD2" \
          --filter-expression "FS > 60.0" --filter-name "FS60" \
          --filter-expression "MQ < 40.0" --filter-name "MQ40" \
          -O results/variants/{SAMPLE}.filtered_snps.vcf.gz

        gatk VariantFiltration -R {input.ref} \
          -V results/variants/{SAMPLE}.raw_indels.vcf.gz \
          --filter-expression "QD < 2.0" --filter-name "QD2" \
          --filter-expression "FS > 200.0" --filter-name "FS200" \
          -O results/variants/{SAMPLE}.filtered_indels.vcf.gz

        gatk MergeVcfs \
          -I results/variants/{SAMPLE}.filtered_snps.vcf.gz \
          -I results/variants/{SAMPLE}.filtered_indels.vcf.gz \
          -O {output}
        """

rule annotate:
    input:
        f"results/variants/{SAMPLE}.final.vcf.gz"
    output:
        f"results/annotation/{SAMPLE}.annotated.vcf.gz"
    shell:
        """
        conda run -n snpeff-env snpEff ann \
          -Xmx4g GRCh38.99 {input} \
          | bgzip > {output}
        tabix -p vcf {output}
        """

rule python_analysis:
    input:
        f"results/variants/{SAMPLE}.final.vcf.gz"
    output:
        "figures/variant_summary.png"
    shell:
        "python notebooks/variant_analysis.py"

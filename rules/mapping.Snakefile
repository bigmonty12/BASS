rule trim_galore_pe:
    input:
        unpack(get_fastq)
    output:
        temp("trimmed/{sample}.{unit}.1_val_1.fq.gz"),
        temp("trimmed/{sample}.{unit}.1.fastq.gz_trimming_report.txt"),
        temp("trimmed/{sample}.{unit}.2_val_2.fq.gz"),
        temp("trimmed/{sample}.{unit}.2.fastq.gz_trimming_report.txt")
    threads: config["params"]["trim-galore"]["cores"]
    params:
        extra=config["params"]["trim-galore"]["options"]
    log:
        "logs/trim_galore/{sample}.{unit}.log"
    wrapper:
        "0.36.0/bio/trim_galore/pe"

rule map_reads:
    input:
        sample=get_trimmed_reads
    output:
        temp("mapped/{sample}.{unit}.bam")
    log:
        "logs/bowtie2/{sample}.{unit}.log"
    threads: config["params"]["bowtie2"]["cores"]
    params:
        index=config["ref"]["genome"],
        extra=config["params"]["bowtie2"]["options"]
    wrapper:
        "0.36.0/bio/bowtie2/align"

rule first_index:
    input:
        "mapped/{sample}.{unit}.bam"
    output:
        temp("mapped/{sample}.{unit}.bam.bai")
    wrapper:
        "0.36.0/bio/samtools/index"

rule samtools_view:
    input:
        "mapped/{sample}.{unit}.bam"
    output:
        temp("mapped/{sample}.{unit}.filt.bam")
    threads: config["params"]["samtools-view"]["cores"]
    params:
        config["params"]["samtools-view"]["options"]
    wrapper:
        "0.36.0/bio/samtools/view"

rule mark_duplicates:
    input:
        "mapped/{sample}.{unit}.filt.bam"
    output:
        bam=protected("dedup/{sample}.{unit}.bam"),
        metrics="qc/dedup/{sample}.{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.36.0/bio/picard/markduplicates"

rule samtools_index:
    input:
        "dedup/{sample}.{unit}.bam"
    output:
        "dedup/{sample}.{unit}.bam.bai"
    wrapper:
        "0.36.0/bio/samtools/index"

rule samtools_merge:
    input:
        sample=get_sample_bams
    output:
        protected("dedup/{sample}-merged.bam")
    params:
        ""
    threads:
        8
    wrapper:
        "0.36.0/bio/samtools/merge"

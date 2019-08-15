rule trim_galore_pe:
    input:
        unpack(get_fastq)
    output:
        temp("trimmed/{sample}.{unit}.1_val_1.fq.gz"),
        temp("trimmed/{sample}.{unit}.1.fastq.gz_trimming_report.txt"),
        temp("trimmed/{sample}.{unit}.2_val_2.fq.gz"),
        temp("trimmed/{sample}.{unit}.2.fastq.gz_trimming_report.txt")
    threads: 6
    params:
        extra="--cores 6"
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
    params:
        index=config["ref"]["genome"],
        extra="--very-sensitive-local -X 2000"
    threads: 8
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
    threads: 4
    params:
        "-b -@ 4 -h -q 30 chr2L chr2R chr3L chr3R chrX"
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

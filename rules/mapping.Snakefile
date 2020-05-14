rule trim_galore_pe:
    input:
        get_fastq
    output:
        "trimmed/{sample}.{unit,\d+}.1_val_1.fq.gz",
        "trimmed/{sample}.{unit,\d+}.1.fastq.gz_trimming_report.txt",
        "trimmed/{sample}.{unit,\d+}.2_val_2.fq.gz",
        "trimmed/{sample}.{unit,\d+}.2.fastq.gz_trimming_report.txt"
    threads: config["params"]["trim-galore"]["cores"]
    params:
        extra=config["params"]["trim-galore"]["options"]
    log:
        "logs/trim_galore/{sample}.{unit,\d+}.log"
    conda: "../envs/trim_galore.yaml"
    shell:
        """
        trim_galore \
        {params.extra} \
        --paired \
        --output_dir trimmed \
        {input} \
        2> {log}

        """

rule map_reads:
    input:
        sample=get_trimmed_reads
    output:
        temp("mapped/{sample}.{unit,\d+}.bam")
    log:
        "logs/bowtie2/{sample}.{unit,\d+}.log"
    threads: config["params"]["bowtie2"]["cores"]
    params:
        index=config["ref"]["genome"],
        extra=config["params"]["bowtie2"]["options"]
    wrapper:
        "0.36.0/bio/bowtie2/align"

rule sort_reads:
    input:
        "mapped/{sample}.{unit,\d+}.bam"
    output:
        temp("mapped/{sample}.{unit,\d+}.sorted.bam")
    params:
        config["params"]["samtools-sort"]["options"]
    threads: config["params"]["samtools-sort"]["cores"]
    wrapper:
        "0.36.0/bio/samtools/sort"

rule first_index:
    input:
        "mapped/{sample}.{unit,\d+}.sorted.bam"
    output:
        temp("mapped/{sample}.{unit,\d}.sorted.bam.bai")
    wrapper:
        "0.36.0/bio/samtools/index"

rule samtools_view:
    input:
        bam="mapped/{sample}.{unit,\d+}.sorted.bam",
        bai="mapped/{sample}.{unit,\d+}.sorted.bam.bai"
    output:
        temp("mapped/{sample}.{unit,\d+}.filt.bam")
    threads: config["params"]["samtools-view"]["cores"]
    params:
        config["params"]["samtools-view"]["options"]
    conda: "../envs/samtools.yaml"
    shell:
        """
        samtools view {input.bam} {params} > {output}
        """

rule mark_duplicates:
    input:
        "mapped/{sample}.{unit,\d+}.filt.bam"
    output:
        bam=protected("dedup/{sample}.{unit,\d+}.bam"),
        metrics="qc/dedup/{sample}.{unit,\d+}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.{unit,\d+}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.36.0/bio/picard/markduplicates"

rule insert_size:
    input:
        "dedup/{sample}.{unit,\d+}.bam"
    output:
        txt="qc/dedup/{sample}.{unit,\d+}.isize.txt",
        pdf="qc/dedup/{sample}.{unit,\d+}.isize.pdf"
    log:
        "logs/picard/insert_size/{sample}.{unit,\d+}.log"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    wrapper:
        "0.36.0/bio/picard/collectinsertsizemetrics"

rule samtools_index:
    input:
        "dedup/{sample}.{unit,\d+}.bam"
    output:
        "dedup/{sample}.{unit,\d+}.bam.bai"
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

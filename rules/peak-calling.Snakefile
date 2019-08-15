rule peak_calling:
    input:
        bam="dedup/{sample}.{unit}.bam"
    output:
        peaks="macs2/{sample}.{unit}_peaks.narrowPeak",
        excel="mac2/{sample}.{unit}_peaks.xls",
        bed="macs2/{sample}.{unit}_summits.bed"
    conda: "envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.bam} -f BAM -g dm --name {wildcards.sample} --outdir peaks --nomodel --shift -37 --extsize 73 --call-summits --keep-dup all"

rule make_bigwig:
    input:
        "dedup/{sample}.{unit}.bam"
    output:
        "bw/{sample}.{unit}_coverage.bw"
    log:
        "logs/deeptools/{sample}.{unit}.log"
    conda: "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input} -o {output} --extendReads"

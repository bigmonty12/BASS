rule ind_peak_calling:
    input:
        bam="dedup/{sample}.{unit}.bam"
    output:
        peaks="macs2/{sample}.{unit}_peaks.narrowPeak",
        excel="mac2/{sample}.{unit}_peaks.xls",
        bed="macs2/{sample}.{unit}_summits.bed"
    conda: "envs/macs2.yaml"
    params:
        options=config["params"]["macs2"]["options"]
    shell:
        "macs2 callpeak -t {input.bam} --name {wildcards.sample} --outdir macs2 {params.options}"

rule make_bigwig:
    input:
        "dedup/{sample}.{unit}.bam"
    output:
        "bw/{sample}.{unit}_coverage.bw"
    log:
        "logs/deeptools/{sample}.{unit}.log"
    conda: "envs/deeptools.yaml"
    params:
        options=config["params"]["bam-coverage"]["options"]
    shell:
        "bamCoverage -b {input} -o {output} {params.options}"

rule merged_peak_calling:
    input:
        bam="dedup/{sample}-merged.bam"
    output:
        peaks="macs2/{sample}-merged_peaks.narrowPeak",
        excel="macs2/{sample}-merged_peaks.xls",
        bed="macs2/{sample}-merged_summits.bed"
    conda: "envs/macs2.yaml"
    params:
        options=config["params"]["macs2"]["options"]
    shell:
        "macs2 callpeak -t {input.bam} --name {wildcards.sample} --outdir macs2 {params.options}"

rule intersect:
    input:
        a=get1_merged_peaks,
        b=get_rest_merged_peaks
    output:
        "macs2/consensus_peaks.narrowPeak"
    conda: "envs/bedtools.yaml"
    params:
        ""
    shell:
        "bedtools intersect -a {input.a} -b {input.b} > {output}"

rule create_saf:
    input:
        "macs2/consensus_peaks.narrowPeak"
    output:
        "macs2/consensus_peaks.saf"
    shell:
        """
        awk  '{{OFS="\t";print $1"."$2+1"."$3, $1, $2+1, $3, "."}}' {input} > {output}"
        """

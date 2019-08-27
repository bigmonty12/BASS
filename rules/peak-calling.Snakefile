rule ind_peak_calling:
    input:
        bam="dedup/{sample}.{unit}.bam"
    output:
        peaks="macs2/{sample}.{unit}_peaks.narrowPeak",
        excel="macs2/{sample}.{unit}_peaks.xls",
        bed="macs2/{sample}.{unit}_summits.bed"
    log:
        "logs/macs2/{sample}.{unit}.log"
    conda: "../envs/macs2.yaml"
    params:
        options=config["params"]["macs2"]["options"]
    shell:
        """
        macs2 callpeak -t {input.bam} --name {wildcards.sample}.{wildcards.unit} --outdir macs2 {params.options} 2> {log}
        """
rule homer:
    input:
        peaks="macs2/{sample}.{unit}_peaks.narrowPeak"
    output:
        anno="macs2/homer/annotate.{sample}.{unit}.diffexp.txt",
        stats="macs2/homer/{sample}.{unit}.AnnotationStats.txt"
    params:
        genome=config["ref"]["name"]
    conda: "../envs/homer.yaml"
    shell:
        """
        perl "$CONDA_PREFIX"/share/homer-4.9.1-6/configureHomer.pl -install {params.genome}
        annotatePeaks.pl {input.peaks} {params.genome} -annStats {output.stats} > {output.anno}
        """
rule make_bigwig:
    input:
        bam="dedup/{sample}.{unit}.bam",
        bai="dedup/{sample}.{unit}.bam.bai"
    output:
        "bw/{sample}.{unit}_coverage.bw"
    log:
        "logs/deeptools/{sample}.{unit}.log"
    conda: "../envs/deeptools.yaml"
    params:
        options=config["params"]["bam-coverage"]["options"]
    shell:
        """
        bamCoverage -b {input.bam} -o {output} {params.options} 2> {log}
        """

rule merged_peak_calling:
    input:
        bam="dedup/{sample}-merged.bam"
    output:
        peaks="macs2/{sample}-merged_peaks.narrowPeak",
        excel="macs2/{sample}-merged_peaks.xls",
        bed="macs2/{sample}-merged_summits.bed"
    log:
        "logs/macs2/{sample}-merged.log"
    conda: "../envs/macs2.yaml"
    params:
        options=config["params"]["macs2"]["options"]
    shell:
        """
        macs2 callpeak -t {input.bam} --name {wildcards.sample}-merged --outdir macs2 {params.options} 2> {log}
        """

rule intersect:
    input:
        a=get1_merged_peaks,
        b=get_rest_merged_peaks
    output:
        "macs2/consensus_peaks.narrowPeak"
    conda: "../envs/bedtools.yaml"
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
        awk  '{{OFS="\t";print $1"."$2+1"."$3, $1, $2+1, $3, "."}}' {input} > {output}
        """

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}.{unit,\d+}.html",
        zip="qc/fastqc/{sample}.{unit,\d+}_fastqc.zip"
    wrapper:
        "0.36.0/bio/fastqc"

rule samtools_stats:
    input:
        "dedup/{sample}.{unit,\d+}.bam"
    output:
        "qc/samtools-stats/{sample}.{unit,\d+}.txt"
    log:
        "logs/samtools-stats/{sample}.{unit,\d+}.log"
    wrapper:
        "0.36.0/bio/samtools/stats"

rule peaks_qc:
    input:
        bam="dedup/{sample}.{unit,\d+}.bam",
        peaks="macs2/{sample}-merged_peaks.narrowPeak",
        peak="macs2/{sample}.{unit,\d+}_peaks.narrowPeak"
    output:
        peaks="qc/peaks_qc/{sample}.{unit,\d+}.peaks.counts_mqc.tsv",
        frip="qc/peaks_qc/{sample}.{unit,\d+}.peaks.FRiP_mqc.tsv"
    conda: "../envs/samtools.yaml"
    params:
        peaks=qc_peaks,
        frip=qc_frip
    threads: 4
    shell:
        """
        name={wildcards.sample}.{wildcards.unit}

        peak_counts=`cat {input.peak} | wc -l`

        mapped_reads=`samtools view -@ {threads} -c {input.bam}`

        reads_in_peaks=`samtools view -@ {threads} -c -L {input.peaks} {input.bam}`

        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`

        printf "%s\t%s\n" $name $peak_counts | cat {params.peaks} - > {output.peaks}
        printf "%s\t%8.3f\n" $name $frip | cat {params.frip} - > {output.frip}

        """
rule homer_qc:
    input:
        "macs2/homer/{sample}.{unit,\d+}.AnnotationStats.txt"
    output:
        "qc/homer/{sample}.{unit,\d+}.features_mqc.tsv"
    params:
        homer=qc_homer
    shell:
        """
        cut -f1,2 {input} | head -10 | sed '1d' | sort -n -k2 | cat {params.homer} - > {output}
        """

rule df_homer_qc:
    input:
        "results/diffexp/homer/{contrast}.AnnotationStats.txt"
    output:
        "qc/homer/{contrast}.df_features_mqc.tsv"
    params:
        homer=qc_homer
    shell:
        """
        cut -f1,2 {input} | head -10 | sed '1d' | sort -n -k2 | cat {params.homer} - > {output}
        """

rule multiqc:
    input:
        expand(["qc/samtools-stats/{sample_unit}.txt",
                "qc/fastqc/{sample_unit}_fastqc.zip",
                "qc/dedup/{sample_unit}.metrics.txt",
                "qc/dedup/{sample_unit}.isize.pdf",
                "qc/peaks_qc/{sample_unit}.peaks.FRiP_mqc.tsv",
                "qc/peaks_qc/{sample_unit}.peaks.counts_mqc.tsv",
                "qc/homer/{contrast}.df_features_mqc.tsv",
                "qc/homer/{sample_unit}.features_mqc.tsv"],
               sample_unit=SAMPLE_UNITS,
               contrast=config["deseq2"]["contrasts"]),
        "qc/deseq2/featureCounts.summary",
        "qc/deseq2/pca.vals_mqc.tsv",
        "qc/deseq2/heatplot.vars_mqc.tsv"
    output:
        report("qc/multiqc.html",
               caption="../report/multiqc.rst",
               category="Quality Control")
    conda: "../envs/multiqc.yaml"
    params:
        config = qc_config
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc --config {params.config} \
        -m custom_content \
        -m samtools \
        -m picard \
        -m featureCounts \
        -m fastqc \
        --force \
        -o qc \
        -n multiqc.html \
        qc/peaks_qc qc/fastqc qc/samtools-stats qc/deseq2 qc/dedup qc/homer > logs/multiqc/multiqc.log 2>&1"

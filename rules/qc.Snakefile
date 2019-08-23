rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}.{unit}.html",
        zip="qc/fastqc/{sample}.{unit}_fastqc.zip"
    wrapper:
        "0.36.0/bio/fastqc"

rule samtools_stats:
    input:
        "dedup/{sample}.{unit}.bam"
    output:
        "qc/samtools-stats/{sample}.{unit}.txt"
    log:
        "logs/samtools-stats/{sample}.{unit}.log"
    wrapper:
        "0.36.0/bio/samtools/stats"

rule peaks_qc:
    input:
        bam="dedup/{sample}.{unit}.bam",
        peaks="macs2/{sample}-merged_peaks.narrowPeak",
        peak="macs2/{sample}.{unit}_peaks.narrowPeak"
    output:
        peaks="qc/peaks_qc/{sample}.{unit}.peaks.counts_mqc.tsv",
        frip="qc/peaks_qc/{sample}.{unit}.peaks.FRiP_mqc.tsv"
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

rule multiqc:
    input:
        expand(["qc/samtools-stats/{sample}.{unit}.txt",
                "qc/fastqc/{sample}.{unit}_fastqc.zip",
                "qc/dedup/{sample}.{unit}.metrics.txt",
                "qc/dedup/{sample}.{unit}.isize.pdf",
                "qc/peaks_qc/{sample}.{unit}.peaks.FRiP_mqc.tsv",
                "qc/peaks_qc/{sample}.{unit}.peaks.counts_mqc.tsv"],
               sample=units.index.get_level_values('sample').unique().values,
               unit=units.index.get_level_values('unit').unique().values),
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
        qc/peaks_qc qc/fastqc qc/samtools-stats qc/deseq2 qc/dedup > logs/multiqc/multiqc.log 2>&1"

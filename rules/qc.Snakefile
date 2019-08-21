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
        bam = "dedup/{sample}.{unit}.bam",
        peaks = "macs2/{sample}-merged_peaks.narrowPeak",
        peak = "macs2/{sample}.{unit}_peaks.narrowPeak"
    output:
        "qc/peaks_qc/{sample}.{unit}.peaks.qc.txt"
    conda: "../envs/samtools.yaml"
    threads: 4
    shell:
        """
        name={wildcards.sample}.{wildcards.unit}

        peak_counts=`cat {input.peak} | wc -l`

        mapped_reads=`samtools view -@ {threads} -c {input.bam}`

        reads_in_peaks=`samtools view -@ {threads} -c -L {input.peaks} {input.bam}`

        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`

        printf "%s\t%s\t%8.3f\n" $name $peak_counts $frip > {output}

        """

rule merge_peaks_qc:
    input:
        expand("qc/peaks_qc/{sample}.{unit}.peaks.qc.txt",
               sample=units.index.get_level_values('sample').unique().values,
               unit=units.index.get_level_values('unit').unique().values)
    output:
        report("qc/peaks_qc/peaks.qc.txt",
               caption="../report/peaks.rst",
               category="Quality Control")
    shell:
        """
        echo -e "sample\tpeaks\tFRiP" | cat - {input} > {output}
        """

rule multiqc:
    input:
        expand(["qc/samtools-stats/{sample}.{unit}.txt",
                "qc/fastqc/{sample}.{unit}_fastqc.zip",
                "qc/dedup/{sample}.{unit}.metrics.txt"],
               sample=units.index.get_level_values('sample').unique().values,
               unit=units.index.get_level_values('unit').unique().values)
    output:
        report("qc/multiqc.html",
               caption="../report/multiqc.rst",
               category="Quality Control")
    log:
        "logs/multiqc.log"
    wrapper:
        "0.36.0/bio/multiqc"

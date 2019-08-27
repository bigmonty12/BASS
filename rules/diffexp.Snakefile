def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        samples=config["deseq2"]["samples"],
        bams="dedup",
        saf="macs2/consensus_peaks.saf"
    output:
        dds="deseq2/all.rds",
        fc="qc/deseq2/featureCounts.summary"
    conda: "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca_heatmap:
    input:
        "deseq2/all.rds"
    output:
        pca="results/pca.vals.txt",
        heatplot="results/heatplot.vals.txt"
    params:
        pca_labels=config["pca"]["labels"]
    conda: "../envs/deseq2.yaml"
    log:
        "logs/deseq2/pca.log"
    script:
        "../scripts/plot-pca.R"

rule deseq_multiqc:
    input:
        pca="results/pca.vals.txt",
        heatmap="results/heatplot.vals.txt"
    output:
        pca="qc/deseq2/pca.vals_mqc.tsv",
        heatmap="qc/deseq2/heatplot.vars_mqc.tsv"
    params:
        pca=qc_pca,
        heatmap=qc_heatmap
    shell:
        """
        cat {params.pca} {input.pca} > {output.pca}
        cat {params.heatmap} {input.heatmap} > {output.heatmap}
        """


def get_contrast(wildcards):
    return config["deseq2"]["contrasts"][wildcards.contrast]


rule deseq2_report:
    input:
        "deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv",
                     caption="../report/de_table.rst",
                     category="DESeq2"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg",
                       caption="../report/ma_plot.rst",
                       category="DESeq2")
    params:
        contrast=get_contrast
    conda: "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule annotate:
    input:
        "results/diffexp/{contrast}.diffexp.tsv"
    output:
        "results/diffexp/annotate.{contrast}.diffexp.txt"
    params:
        genome=config["ref"]["name"]
    conda: "../envs/homer.yaml"
    shell:
        """
        perl "$CONDA_PREFIX"/share/homer-4.9.1-6/configureHomer.pl -install {params.genome}
        annotatePeaks.pl {input} {params.genome} > {output}
        """

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def get_design(wildcards):
    return config["deseq2"]["designs"][wildcards.design]

def get_contrast(wildcards):
   x = str(wildcards.design) + "." + str( wildcards.contrast)
   return config["deseq2"]["contrasts"][x]

rule deseq2_init:
    input:
        samples=config["deseq2"]["samples"],
        # bams="dedup",
        saf="macs2/consensus_peaks.saf"
    output:
        dds="deseq2/{design}.rds",
        fc="qc/deseq2/{design}.featureCounts.summary"
    params:
        design=get_design
    conda: "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{design}.init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca_heatmap:
    input:
        "deseq2/condition.rds"
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


rule deseq2_report:
    input:
        rds="deseq2/{design}.rds"
    output:
        table=report("results/diffexp/{design}.{contrast}.diffexp.tsv",
                     caption="../report/de_table.rst",
                     category="DESeq2"),
        ma_plot=report("results/diffexp/{design}.{contrast}.ma-plot.svg",
                       caption="../report/ma_plot.rst",
                       category="DESeq2")
    params:
        contrast_name=get_contrast,
        pval=config["deseq2"]["pval"]
    conda: "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{design}.{contrast}.diffexp.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"

rule annotate_diffexp:
    input:
        "results/diffexp/{design}.{contrast}.diffexp.tsv"
    output:
        anno="results/diffexp/homer/annotate.{design}.{contrast}.diffexp.txt",
        stats="results/diffexp/homer/{design}.{contrast}.AnnotationStats.txt",
        go=directory("results/diffexp/homer/{design}.{contrast}.go")
    params:
        genome=config["ref"]["name"]
    conda: "../envs/homer.yaml"
    shell:
        """
        perl "$CONDA_PREFIX"/share/homer-4.10-0/configureHomer.pl -install {params.genome}
        annotatePeaks.pl {input} {params.genome} -annStats {output.stats} -go {output.go} > {output.anno}
        """

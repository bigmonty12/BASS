def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        samples=config["deseq2"]["samples"],
        bams=expand(("dedup/{sample}.{unit}.bam"),
                    sample=units.index.get_level_values('sample').unique().values,
                    unit=units.index.get_level_values('unit').unique().values),
        saf="macs2/consensus_peaks.saf"
    output:
        "deseq2/all.rds"
    conda: "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/all.rds"
    output:
        report("results/pca.svg", caption="../report/pca.rst", category="DESeq2")
    params:
        pca_labels=config["pca"]["labels"]
    conda: "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


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

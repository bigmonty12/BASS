# Run with --directory option (and path of directory containing Fastq folder)

include: srcdir("rules/common.Snakefile")

# Target rules #

rule all:
    input:
        "qc/multiqc.html",
        expand("bw/{sample}.{unit}_coverage.bw",
               sample=units.index.get_level_values('sample').unique().values,
               unit=units.index.get_level_values('unit').unique().values),
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["deseq2"]["contrasts"])


# Modules #

include: srcdir("rules/mapping.Snakefile")
include: srcdir("rules/peak-calling.Snakefile")
include: srcdir("rules/diffexp.Snakefile")
include: srcdir("rules/qc.Snakefile")

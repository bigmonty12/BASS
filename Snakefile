# Run with --directory option (and path of directory containing Fastq folder)

include: srcdir("rules/common.Snakefile")

# Target rules #

SAMPLE_UNITS = create_combos(sample_units)

rule all:
    input:
        "qc/multiqc.html",
        expand(["bw/{sample_unit}_coverage.bw",
               "macs2/homer/annotate.{sample_unit}.diffexp.txt"],
               sample_unit=SAMPLE_UNITS),
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["deseq2"]["contrasts"]),



# Modules #

include: srcdir("rules/mapping.Snakefile")
include: srcdir("rules/peak-calling.Snakefile")
include: srcdir("rules/diffexp.Snakefile")
include: srcdir("rules/qc.Snakefile")

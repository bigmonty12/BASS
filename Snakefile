include: "rules/common.Snakefile"

# Target rules #

rule all:
    input:
        "qc/multiqc.html",
        expand("macs2/{sample}.{unit}_peaks.narrowPeak",
               sample=units.index.get_level_values('sample').unique().values,
               unit=units.index.get_level_values('unit').unique().values),
        expand("bw/{sample}.{unit}_coverage.bw",
               sample=units.index.get_level_values('sample').unique().values,
               unit=units.index.get_level_values('unit').unique().values)

# Modules #

include: "rules/mapping.Snakefile"
include: "rules/peak-calling.Snakefile"
include: "rules/qc.Snakefile"

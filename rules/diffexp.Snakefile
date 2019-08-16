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
    conda: "envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "scripts/deseq2-init.R"

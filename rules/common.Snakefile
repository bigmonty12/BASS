import pandas as pd
report: srcdir("../report/workflow.rst")

# Config file and sample sheets #
configfile: "config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)


units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"],
                                                            drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index


# create dict of samples with respective replicates
sample_units = pd.read_table(config["units"], dtype=str).iloc[:,0:2].groupby('sample')['unit'].apply(lambda x: x.tolist()).to_dict()

# Multiqc files
qc_config = srcdir("../multiqc/multiqc_config.yaml")
qc_peaks = srcdir("../multiqc/peak_count_header.txt")
qc_frip = srcdir("../multiqc/frip_score_header.txt")
qc_pca = srcdir("../multiqc/deseq2_pca_header.txt")
qc_heatmap = srcdir("../multiqc/deseq2_clustering_header.txt")
qc_homer = srcdir("../multiqc/homer_anno_header.txt")


# Helper Functions #
def create_combos(d):
    x = list()
    for key,value in d.items():
        for v in value:
            x.append(key + "." + str(v))
    return x


def get_fastq(wildcards):
    """ Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit),
                       ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return [fastqs.fq1, fastqs.fq2]
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """ Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """ Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """ Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        fastqs = expand("trimmed/{sample}.{unit}.{group}_val_{group}.fq.gz",
                        group=[1, 2], **wildcards)
        return fastqs
    # single end sample
    return "trimmed/{sample}.{unit}.{group}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """ Get all aligned reads of given sample."""
    bams = expand("dedup/{sample}.{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)
    return bams


def get1_merged_peaks(wildcards):
    peaks = expand("macs2/{sample}-merged_peaks.narrowPeak",
                   sample=units.index.get_level_values('sample').unique().values)
    return peaks[0]


def get_rest_merged_peaks(wildcards):
    peaks = expand("macs2/{sample}-merged_peaks.narrowPeak",
                   sample=units.index.get_level_values('sample').unique().values)
    return peaks[1:]

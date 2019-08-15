import pandas as pd

samples = pd.read_table("samples.tsv").set_index(["sample"], drop=False)

units = pd.read_table("units.tsv").set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

def get_fastq():
    fastqs = units.loc[(sample, unit), ["fq1", "fq2"]].dropna()
    print(fastqs)
    return

# get_fastq()

fastqs = units.loc[("EtOH", "2"), ["fq1", "fq2"]].dropna()

d = {"r1": fastqs.fq1, "r2": fastqs.fq2}
print(units)

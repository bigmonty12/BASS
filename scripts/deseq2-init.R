log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("dplyr")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

threads <- snakemake@threads

# make sure colData and countData have the same sample order
samples <- read.table(
    snakemake@input[["samples"]],
    header = T)
print(samples)
counts <- readRDS(snakemake@input[["counts"]])
print(counts)
dds_design = snakemake@params[["design"]]
print(dds_design)
if (dds_design == 'condition') {
    dds <- DESeqDataSetFromMatrix(
        counts,
        samples,
        design = ~condition) 
} else if (dds_design == 'Factor') {
    dds <- DESeqDataSetFromMatrix(
        counts,
        samples,
        design = ~Factor)
} else if (dds_design == 'Tissue') {
    dds <- DESeqDataSetFromMatrix(
        counts,
        samples,
        design = ~Tissue)
} else if (dds_design == 'genotype') {
    dds <- DESeqDataSetFromMatrix(
        counts,
        samples,
        design = ~genotype)
}


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])

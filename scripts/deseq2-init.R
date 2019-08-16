log <- file(snakemake@log[[1]]), open="wt")
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
    threads <- snakemake@threads
}


# make sure colData and countData have the same sample order
samples <- reads.table(
    snakemake@input[["samples"]],
    header = T)

filenames <- snakemake@input[["bams"]]
filenames <- strsplit(filenames, " ")

saf <- snakemake@input[["saf"]]

f_counts <- featureCounts(
    files = filenames,
    annot.ext = saf,
    isGTFAnnotationFile = FALSE,
    isPairedEnd = TRUE,
    nthreads = threads)

colnames(f_counts$counts) <- samples$SampleID

counts <- as.data.frame(
    f_counts$counts)

dds <- DESeqDataSetFromMatrix(
    counts,
    samples,
    design = ~condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])


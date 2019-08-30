log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("dplyr")
library("Rsubread")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

threads <- snakemake@threads
print(threads)

# make sure colData and countData have the same sample order
samples <- read.table(
    snakemake@input[["samples"]],
    header = T)

dir <- snakemake@input[["bams"]]

filenames <- file.path(
    dir,
    paste0(samples$SampleID, ".bam"))

print(samples$SampleID)
print(filenames)

saf <- snakemake@input[["saf"]]

f_counts <- featureCounts(
    files = filenames,
    annot.ext = saf,
    isGTFAnnotationFile = FALSE,
    isPairedEnd = TRUE,
    nthreads = threads)

colnames(f_counts$counts) <- samples$SampleID

stats <- f_counts[["stat"]]
colnames(stats) <- c("Status", as.character(samples$SampleID))

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

write.table(
    as.data.frame(
        stats),
    file=snakemake@output[[2]],
    quote=FALSE,
    row.names=FALSE,
    sep="\t")

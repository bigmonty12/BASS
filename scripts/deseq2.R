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
    threads <- snakemake@threads
}

dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])
res <- results(
    dds,
    contrast=contrast,
    parallel=parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(
    dds,
    contrast=contrast,
    res=res)
# sort by p-value
res <- res[order(res$padj),]
res_ma <- res
# add peak ranges for annotation
r <- rownames(res)
r_split <- strsplit(r, ".", 2)

res$PeakID <- r
res$Chr <- sapply(r_split, "[", 1)
res$Start <- sapply(r_split, "[", 2)
res$End <- sapply(r_split, "[", 3)
res$Strand <- "."

res <- as.data.frame(res) %>%
        dplyr::select(
            Chr,
            Start,
            End,
            PeakID,
            baseMean,
            Strand,
            dplyr::everything()) %>%
        dplyr::filter(
            padj < snakemake@params[["pval"]])

# store results
svg(snakemake@output[["ma_plot"]])
plotMA(
    res_ma,
    ylim=c(-2,2))
dev.off()

write.table(
        res,
        row.names=F,
        col.names=F,
        quote=F,
        sep="\t",
        file=snakemake@output[["table"]])



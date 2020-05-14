log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)

# pca

pca.data <- DESeq2::plotPCA(
    counts,
    intgroup=snakemake@params[["pca_labels"]],
    returnData=TRUE)
percentVar <- round(100 * attr(pca.data, "percentVar"))

pca.vals <- pca.data[,1:2]
colnames(pca.vals) <- paste(
    colnames(pca.vals),
    paste(
        percentVar,
        '% variance',
        sep=""),
    sep=": ")
pca.vals <- cbind(
    sample = rownames(pca.vals),
    pca.vals)
write.table(
    pca.vals,
    file=snakemake@output[[1]],
    row.names=FALSE,
    col.names=TRUE,
    sep="\t",
    quote=FALSE)

# heatmap
samplesDists <- dist(t(assay(counts)))
sampleDistMatrix <- as.matrix(samplesDists)
print(snakemake@output[[1]])
print(snakemake@output[[2]])
write.table(
    cbind(
        sample = rownames(sampleDistMatrix),
        sampleDistMatrix),
    file=snakemake@output[[2]],
    row.names=FALSE,
    col.names=TRUE,
    sep="\t",
    quote=FALSE)

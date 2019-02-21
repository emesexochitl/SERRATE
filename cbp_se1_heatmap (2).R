select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
pheatmap(log2.norm.counts, cluster_rows=F, show_rownames=T, cluster_cols=T, main = "DESeq2 log2 normalized heatmap")

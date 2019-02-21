#general data load
countdata <- read.table("ars2_human_counts_all_s2_mm.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\_accepted_hits.[sb]am$", "", colnames(countdata))
coldata <- read.table("coldata.txt", header=T, row.names=1)

# specific variables
control <- c("wt_1", "wt_2", "wt_3")
treatment <- c("ars2_1", "ars2_2", "ars2_3")
control_name <- c("WT")
treatment_name <- c("mutant")
#samples <- c(control, treatment)


#specify comparison
#countdata2 <- countdata[ , which(names(countdata) %in% samples)]
#coldata2 <- coldata[ which(rownames(coldata) %in% samples), ]
#coldata2 <- read.table("coldata2.txt", header=T, row.names=1)

#do the differential expression
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds))>1, ] 
dds$condition <- relevel(dds$condition, ref="WT") 
dds <- DESeq(dds) 
res <- results(dds) 
res 
summary(res) 
summary(res, alpha = 0.05)
resOrdered <- res[order(res$padj), ]
resSig01 <- res[which(res$padj < 0.1), ]
resSig005 <- res[which(res$padj < 0.05), ]
rld <- rlog(dds, blind = F) 
data <- plotPCA(rld, intgroup="condition", returnData =T) 
percentVar <- round(100*attr(data, "percentVar")) 
ggplot(data, aes(PC1, PC2, color=condition, shape=condition)) + geom_point(size=9) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + ggtitle(paste0("PCA plot of WT and ", treatment_name,  " samples")) + theme(plot.title= element_text(lineheight = 1.6, face = "bold", size = 24)) 
plotMA(res[order(res$padj),], main=paste0("DESeq2 WT and ", treatment_name, " samples, padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
plotMA(res[order(res$padj),], main=paste0("DESeq2 WT and ", treatment_name, " samples, padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
hist(res$pvalue, breaks = 20, col="lightskyblue1", border = F, main = paste0("Histogram of pvals of DESeq2\nWT and ", treatment_name, " samples"), xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))

png(filename=paste0("maplot_WT_", treatment_name, "_01.png"))
plotMA(res[order(res$padj),], main=paste0("DESeq2 WT and ", treatment_name, " samples, padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
dev.off()

png(filename=paste0("maplot_WT_", treatment_name, "_005.png"))
plotMA(res[order(res$padj),], main=paste0("DESeq2 WT and ", treatment_name, " samples, padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
dev.off()

png(filename=paste0("hist_WT_", treatment_name, ".png"))
hist(res$pvalue, breaks = 20, col="lightskyblue1", border = F, main = paste0("Histogram of pvals of DESeq2\nWT and ", treatment_name, " samples"), xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))
dev.off()

png(filename=paste0("pca_drozi_ars2.png"))
ggplot(data, aes(PC1, PC2, color=condition, shape=condition)) + geom_point(size=9) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + ggtitle(paste0("PCA plot of WT and ", "\n", treatment_name, " samples")) + theme(plot.title= element_text(lineheight = 1.6, face = "bold", size = 24))
dev.off()

#sum(res$padj < 0.05, na.rm = T)

write.table((res[which(res$padj < 0.05), ]), file = paste0("deseq2_WT_", treatment_name, "_005.txt"), quote=F)
write.table((res[which(res$padj < 0.1), ]), file = paste0("deseq2_WT_", treatment_name, "_01.txt"), quote=F)

intdata <- read.table(file = "Homo_005_log2_1updown_intronnumber_iso1.txt", header = F)
table(intdata$V1)
intdata$V1 <- factor(intdata$V1, levels=c("Homo_genomic_all", "Homo_ars2_all", "Homo_ars2_down", "Homo_ars2_up"))
ggplot(data = intdata, aes(x = V1, y = V2, fill = intdata$V1)) + geom_boxplot(outlier.colour = NA) + ggtitle(label = "Intron numbers of Homo sapiens\nWT and ars2 RNA seq libraries(padj<0.05, 1st isoforms)") + theme_classic() + theme(plot.title=element_text(size = 20)) + xlab(label = "") + ylab(label = "Number of introns") + scale_x_discrete(labels = c("Genomic background", "All sig.\n unfiltered genes", "Downregulated\n (log2FC<=-1)", "Upregulated\n (log2FC>=1)")) + scale_fill_manual(values = rev(getPalette(4))) + geom_signif(comparisons = list(c("Homo_genomic_all", "Homo_ars2_all"), c("Homo_genomic_all", "Homo_ars2_down"), c("Homo_genomic_all", "Homo_ars2_up")), map_signif_level =TRUE, step_increase = 0.1, tip_length = 0.001, y_position=c(27,29,32)) + coord_cartesian(ylim=c(0,35))


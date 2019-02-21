#general data load
countdata <- read.table("counts_all.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\_accepted_hits.[sb]am$", "", colnames(countdata))
coldata <- read.table("coldata.txt", header=T, row.names=1)

# specific variables
treatment <- c("DRR006600", "DRR006709")
control <- c("DRR006654", "DRR006655")
samples <- c(treatment, control)
treatment_time <- c("12h")
treatment_typ <- c("osmotic")
tissue <- c("shoot")

#specify comparison
countdata2 <- countdata[ , which(names(countdata) %in% samples)]
coldata2 <- coldata[ which(rownames(coldata) %in% samples), ]

#do the differential expression
d
dds <- dds[rowSums(counts(dds))>1, ] 
dds$condition <- relevel(dds$condition, ref="control") 
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
ggplot(data, aes(PC1, PC2, color=condition, shape=condition)) + geom_point(size=9) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + ggtitle(paste0("PCA plot of control and ", "\n", treatment_typ, " treated ", tissue, " samples (", treatment_time, ")")) + theme(plot.title= element_text(lineheight = 1.6, face = "bold", size = 24)) 
plotMA(res[order(res$padj),], main=paste0("DESeq2 control and ", treatment_typ, " treated ", "\n", tissue, " samples (", treatment_time, "), padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
plotMA(res[order(res$padj),], main=paste0("DESeq2 control and ", treatment_typ, " treated ", "\n", tissue, " samples (", treatment_time, "), padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
hist(res$pvalue, breaks = 20, col="lightskyblue1", border = F, main = paste0("Histogram of pvals of DESeq2\ncontrol and ", treatment_typ, " treated ", tissue,  " samples (", treatment_time, ")"), xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))

png(filename=paste0("maplot_", treatment_typ,"_", tissue, "_", treatment_time,"_01.png"))
plotMA(res[order(res$padj),], main=paste0("DESeq2 control and ", treatment_typ, " treated ", "\n", tissue, " samples (", treatment_time, "), padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
dev.off()

png(filename=paste0("maplot_", treatment_typ, "_", tissue, "_", treatment_time,"_005.png"))
plotMA(res[order(res$padj),], main=paste0("DESeq2 control and ", treatment_typ, " treated ", "\n", tissue, " samples (", treatment_time, "), padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
dev.off()

png(filename=paste0("hist_", treatment_typ, "_", tissue, "_", treatment_time,".png"))
hist(res$pvalue, breaks = 20, col="lightskyblue1", border = F, main = paste0("Histogram of pvals of DESeq2\ncontrol and ", treatment_typ, " treated ", tissue,  " samples (", treatment_time, ")"), xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))
dev.off()

png(filename=paste0("pca_", treatment_typ, "_", tissue, "_", treatment_time,".png"))
ggplot(data, aes(PC1, PC2, color=condition, shape=condition)) + geom_point(size=9) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + ggtitle(paste0("PCA plot of control and ", "\n", treatment_typ, " treated ", tissue, " samples (", treatment_time, ")")) + theme(plot.title= element_text(lineheight = 1.6, face = "bold", size = 24)) 
dev.off()

#sum(res$padj < 0.05, na.rm = T)

write.table((res[which(res$padj < 0.05), ]), file = paste0("deseq2_", treatment_typ, "_", tissue, "_", treatment_time,"_005.txt"), quote=F) 
write.table((res[which(res$padj < 0.1), ]), file = paste0("deseq2_", treatment_typ, "_", tissue, "_", treatment_time,"_01.txt"), quote=F)
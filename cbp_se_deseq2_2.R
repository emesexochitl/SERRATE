# Import and edit read count file.
countdata <- read.table("counts_all_mod.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)] 
colnames(countdata) <- gsub("\\_accepted_hits.[sb]am$", "", colnames(countdata))

# Set up DESeq2 and run it.
mutant <- c("se1")

coldata <- read.table("coldata.txt", header=T, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds))>1, ] 
dds$condition <- relevel(dds$condition, ref="WT") # for double-check.
dds <- DESeq(dds) 
res <- results(dds) 
res <- results(dds, contrast=c("condition", mutant, "WT")) # Note that the control is always the third value!
res 
summary(res) 
summary(res, alpha = 0.05)

resOrdered <- res[order(res$padj), ]
resSig01 <- res[which(res$padj < 0.1), ]
resSig005 <- res[which(res$padj < 0.05), ]
#sum(res$padj < 0.05, na.rm = T)
#head(res[which(res$padj < 0.1), ])
# Visualization.
#plotMA(res[order(res$padj),], main=paste0("DESeq2 control (WT) and ", mutant, " samples, padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
#plotMA(res[order(res$padj),], main=paste0("DESeq2 control (WT) and ", mutant, " samples, padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
#hist(res$pvalue, breaks = 20, main = paste0("Histogram of pvals of DESeq2\ncontrol (WT) and ", mutant, " samples"), col="lightskyblue1", border = F, xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))

#png(filename=paste0("maplot_WT_", mutant, "_01.png"))
#plotMA(res[order(res$padj),], main=paste0("DESeq2 control (WT) and ", mutant, " samples, padj < 0.1"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)))
#dev.off()

#png(filename=paste0("maplot_WT_", mutant, "_005.png"))
#plotMA(res[order(res$padj),], main=paste0("DESeq2 control (WT) and ", mutant, " samples, padj < 0.05"), ylim = c(1.2*min(res$log2FoldChange),1.2*max(res$log2FoldChange)), alpha = 0.05)
#dev.off()

#png(filename=paste0("hist_WT_", mutant, ".png"))
#hist(res$pvalue, breaks = 20, main = paste0("Histogram of pvals of DESeq2\ncontrol (WT) and ", mutant, " samples"), col="lightskyblue1", border = F, xlab="p-value", ylim = c(0, 1.2*(sum(res$pvalue < 0.05, na.rm = T))))
#dev.off()

# Write significant results to txt files.
#write.table((res[which(res$padj < 0.05), ]), file = paste0("deseq2_WT_", mutant, "_005.txt"), quote=F) 
#write.table((res[which(res$padj < 0.1), ]), file = paste0("deseq2_WT_", mutant, "_01.txt"), quote=F)
#write.table(res[order(res$padj), ], file = paste0("deseq2_WT_", mutant, "_all.txt"), quote =F)
#write.table((res[which(res$padj < 0.05 & res$log2FoldChange> 0.59 | res$log2FoldChange < -0.59), ]), file = paste0("deseq2_WT_", mutant, "_fc1.5_005.txt"), quote=F)
#write.table((res[which(res$padj < 0.01 & res$log2FoldChange> 0.59 | res$log2FoldChange < -0.59), ]), file = paste0("deseq2_WT_", mutant, "_fc1_01.txt"), quote=F) 

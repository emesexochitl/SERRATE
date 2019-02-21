countdata <- read.table("counts_all_mod.txt", header=TRUE, row.names=1) #must change! 1
countdata <- countdata[ ,6:ncol(countdata)] 
#countdata <- countdata[ , -c(5)]
colnames(countdata) <- gsub("\\_accepted_hits.[sb]am$", "", colnames(countdata))

#condition <- factor(c(rep("WT", 3), rep("bdr5-1", 2), rep("bdr5-3", 3))) # must change! 1

#(coldata <- data.frame(row.names=colnames(countdata), condition)) 
coldata <- read.table("coldata.txt", header=T, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds))>1, ] 
dds$condition <- relevel(dds$condition, ref="WT") # must change! 2
dds <- DESeq(dds) 
res <- results(dds) 
res 
summary(res) 
summary(res, alpha = 0.05)

resmod <- results(dds, contrast=c("condition", "WT", "se1"))
resmod 
summary(resmod) 
summary(resmod, alpha = 0.05)

head(res[which(res$padj < 0.1), ])
head(resmod[which(res$padj < 0.1), ])
resOrdered <- res[order(res$padj), ]
rld <- rlog(dds, blind = F) 
data <- plotPCA(rld, intgroup="condition", returnData =T) 
percentVar <- round(100*attr(data, "percentVar")) 
ggplot(data, aes(PC1, PC2, color=group, shape=condition)) + geom_point(size=9) + geom_text(aes(label=condition), size =3, vjust=-1.75,colour = "black") + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + ggtitle("PCA plot of RNA-Seq samples") +theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 


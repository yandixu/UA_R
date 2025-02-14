if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

readcounts <- read.table("counts_gene.txt", header=TRUE)
row.names(readcounts) <- readcounts$gene_id
readcounts <- readcounts[,-c(13)]
readcounts <- readcounts[rowSums(readcounts>=12) > 6, ]

#names(readcounts) <- gsub("SRR12066","",names(readcounts))

dim(readcounts) 

condition = c("test","ctrl","test","ctrl","test","ctrl","test","ctrl","test","ctrl","test","ctrl")

sample_info <- data.frame (row.names=names(readcounts), condition = c("test","ctrl","test","ctrl","test","ctrl","test","ctrl","test","ctrl","test","ctrl"))

DESeq.ds <- DESeqDataSetFromMatrix (countData = readcounts, colData = sample_info, design = ~condition)
DESeq.ds <- estimateSizeFactors (DESeq.ds)
colData(DESeq.ds)
DESeq.rlog <- rlog (DESeq.ds, blind = TRUE)
rlog.counts <- assay(DESeq.rlog)
head(rlog.counts)

library (ggplot2)
plotPCA (DESeq.rlog)

write.table (rlog.counts, file = "rlog.table_1.txt",sep = "\t", quote = FALSE, row.names = FALSE)

colData(DESeq.ds)$condition <- relevel (colData(DESeq.ds)$condition, "ctrl")
DESeq.ds <- DESeq (DESeq.ds)
DGE.results <- results (DESeq.ds, independentFiltering = TRUE , alpha = 0.05)
summary (DGE.results) 

head (DGE.results, n=10)
table (DGE.results$padj < 0.05)
row.names (subset(DGE.results, padj < 0.05))

hist (DGE.results$padj)
hist (DGE.results$padj, 
      col='grey', breaks=20, 
      xlab = "adjusted p-value", ylab = "number of genes", main = "frequencies of p-values")

plotMA (DGE.results, alpha = 0.05 , main = "mutants", ylim = c(-10,10))

#install NMF
BiocManager::install("NMF")
library (NMF)

DGEgenes <- row.names(subset(DGE.results, padj < 0.05))

mat_DGEgenes <- rlog.counts[DGEgenes, ]

pdf("heatmap_3.pdf")
#aheatmap (mat_DGEgenes,color = colorRampPalette(c("yellow", "darkblue", "grey"))(50),Rowv = TRUE , Colv = TRUE ,distfun = "euclidean", hclustfun = "median",scale = "row") 
aheatmap (mat_DGEgenes,Rowv = TRUE , Colv = TRUE ,distfun = "euclidean", hclustfun = "median",scale = "row") 
dev.off()

DGE.dexpress <- DGE.results

upregulated.genes <- row.names(subset(DGE.dexpress, log2FoldChange >= 1))
upregulated.genes <- gsub("[.]([0-9]+)","",upregulated.genes)

write(upregulated.genes, file = "upregulatedGenes.txt", sep = "\t")

table(DGE.dexpress$log2FoldChange >=1)
table(DGE.dexpress$log2FoldChange <= -1)

downregulated.genes <- row.names(subset(DGE.dexpress, log2FoldChange <= -1))
downregulated.genes <- gsub("[.]([0-9]+)","",downregulated.genes)

write(upregulated.genes, file = "downregulatedGenes.txt", sep = "\t")


#########################Pipeline################################
#par( mfrow =c(2, 1))
#boxplot (readcounts, notch = TRUE, main = "untransformed read counts", ylab = "read counts")
#boxplot (rlog.counts, notch = TRUE, main = "rlog-transformed read counts", ylab = "rlog (read counts)")


#head(rlog.counts, n=10)
#cor(rlog.counts[,1], rlog.counts[,10])
#cor(rlog.counts[,1], rlog.counts[,10], method = "spearman")

#cor(rlog.counts, method = "spearman")

# clustering
#par( mfrow =c(2, 1))
#distance.readcounts <- as.dist (1 - cor(readcounts, method = "pearson"))
#plot(hclust(distance.readcounts)) 
#distance.rlog <- as.dist (1 - cor(rlog.counts, method = "pearson")) 
#plot(hclust(distance.rlog)) 



#######################################Annotation##############################
BiocManager :: install("org.Sc.sgd.db")
library (org.Sc.sgd.db)

keytypes (org.Sc.sgd.db)
columns (org.Sc.sgd.db)

#"CHR" doesn't work here.
anno <- select (org.Sc.sgd.db, keys = DGEgenes,keytype = "ALIAS", columns = c("GENENAME","PMID"))

subset (anno, GENENAME == "ctrl")

#plotCounts(dds = DESeq.ds, gene = "AC002398.1",transform = TRUE) 

out.df <- merge (as.data.frame(DGE.results), anno, by.x="row.names", by.y="ORF")

write.table (out.df, file = "result.txt",sep = "\t", quote = FALSE, row.names = FALSE)
row.names(out.df)


#install use org.Hs.eg.db for human gene
BiocManager :: install ("org.Hs.eg.db")
library (org.Hs.eg.db)
keytypes (org.Hs.eg.db)
columns (org.Hs.eg.db)
#anno <- select(org.Hs.eg.db,keys = c("E2F1","RB1"), keytype = "SYMBOL", columns = c("GENENAME", "MAP", "PMID"))
anno <- select (org.Hs.eg.db, keys = DGEgenes,keytype = "ENTREZID", columns = c("GENENAME", "MAP","PMID"))

DGE.results$

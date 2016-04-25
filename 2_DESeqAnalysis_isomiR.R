# Paul Wackers
# Perform a DEG analysis on isomiRs
# 30-10-2015
# R-3.2.1
# cn-03
# /zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1209-Paul_Lucassen/MAD1209-P001-brain_tissue/MAD1209-P001-E001_2014_FFPE_miRNASeq_svleeuw1

rm(list=ls())
options(stringsAsFactors = FALSE)
X11.options(type="Xlib")

design <- read.delim("./ExperimentInfo/SampleInformation/seqdesign.txt",stringsAsFactors=FALSE)
df <- read.delim("./Results/handOverFiles/countTables/isomiR_counttable.csv",stringsAsFactors=FALSE)

geneIDs <- df[,c(1:4)]
CountTable <- data.matrix(df[,-c(1:4)])

design$SampleID
colnames(CountTable)
### NIET OK!!!
colnames(CountTable)[order(colnames(CountTable))]
# Ok!
CountTable <- CountTable[, order(colnames(CountTable))]
colnames(CountTable)
# Ok!
write.table(CountTable,"./Results/isomiR_analysis/CountTable.raw.txt",sep="\t",row.names=FALSE,quote=FALSE)
tmp <- data.frame(geneIDs,CountTable)
write.table(tmp,"./Results/CountTable.isomiRs.raw.txt",sep="\t",row.names=FALSE,quote=FALSE)

sizeFactors.mad <- function (counts, locfunc = median){
    loggeomeans <- rowMeans(log(counts))
    apply(counts, 2, function(cnts) exp(locfunc((log(cnts) - 
        loggeomeans)[is.finite(loggeomeans)])))
}
sf <- sizeFactors.mad(CountTable)

CountTable.norm <- CountTable
for(i in 1:ncol(CountTable.norm)){CountTable.norm[,i] <- CountTable.norm[,i]/sf[i]}

boxplot(data.frame(log2(as.matrix(CountTable)+1)),pch=".")
X11()
boxplot(data.frame(log2(as.matrix(CountTable.norm)+1)),pch=".", main="gewone norm")

colour <- c("red", "blue")
for(i in 1:ncol(CountTable.norm)){
  if(i == 1){
    plot(density(log2(as.matrix(CountTable.norm[ ,i])+1)), ylim=c(0, 0.45), col=colour[design$Group[i]])
  }
  if(i > 1){
    lines(density(log2(as.matrix(CountTable.norm[ ,i])+1)), col=colour[design$Group[i]])
  }
}

library(annotate)
library(DESeq)

geneIDs.tmp <- c()
for(i in 1:nrow(geneIDs)){
  geneIDs.tmp <- c(geneIDs.tmp, paste(geneIDs[i,1], geneIDs[i,2], geneIDs[i,3], geneIDs[i,4],sep="_"))
}
rownames(CountTable) <- geneIDs.tmp
cds <- newCountDataSet(CountTable,as.factor(design$Group))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)
 
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(design$Group)[2], unique(design$Group)[1])
  
completePlotPath= paste0("./Results/isomiR_analysis/images/MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()
  
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

## plot Padj Histogram
completePlotPath <- paste0("./Results/isomiR_analysis/images/histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()
  
completeTablePath = paste0("./Results/isomiR_analysis/DESeqTable.txt")
res.write <- cbind(geneIDs, res[, -1])
all(res$id == geneIDs.tmp) # [1] TRUE
write.table(res.write, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

which(res$log2FoldChange == max(res$log2FoldChange))
CountTable.norm[3557,]
mean(CountTable.norm[3557,1:10])
mean(CountTable.norm[3557,11:20])
res[3557,]




library(ggplot2)
library(ballgown)
library(DESeq2)
library(limma)
setwd("/Users/salvatoreesposito/Desktop/RNAassessment")

# GUI 2439915


# Task 1: Genes

#Making the data files of the experiment and constructing the table
sampleNames <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
pData=data.frame(id=sampleNames, group = rep(c("Group A","Group B","Group C"), each=4))
write.csv(pData, "/Users/salvatoreesposito/Desktop/RNAassessment/Design-Exp.csv")
countTable <- read.csv("gene_count_matrix.csv",row.names=1)
new_order = c("s1.c2", "s2.c2", "s3.c2", "s4.c2","s5.c2","s6.c2","s7.c2", "s8.c2","s9.c2","s10.c2","s11.c2","s12.c2")
countTable = countTable[,new_order]
colTable <- read.csv("Design-Exp.csv",row.names=1)
deseqdata <- DESeqDataSetFromMatrix(countData=countTable,colData=colTable,design=~group)                    
deseqdata

# Remove all the 0 from the gene counts
removeZero <- (rowSums(counts(deseqdata))>0)
deseqdata <- deseqdata[removeZero,]
deseqdata <- DESeq(deseqdata)
head(deseqdata)
str(deseqdata)

#Making a dispersion plot
par(mfrow =c(1,1))
plotDispEsts(deseqdata)

# This function regularlises log transformation
rld = rlog(deseqdata) 
#log2 raw counts
lgc.raw <- log2(counts(deseqdata, normalized = F))
# log2 normalised counts
lgc.norm <- log2(counts(deseqdata, normalized = T))

#This function makes the boxplots from the results  
par(mfrow = c(1,3))
boxplot(lgc.raw,main= "Raw Count", )
boxplot(lgc.norm,main="Normalised Counts")
boxplot(assay(rld),main="Regularised Counts")

# This function makes the density plots
par(mfrow = c(1,3))
plotDensities(lgc.raw,group=deseqdata$group,legend="topright",main= "Raw Count")
plotDensities(lgc.norm,group=deseqdata$group,legend="topright",main="Normalised Count")
plotDensities(assay(rld),group=deseqdata$group,legend="topright",main="Regularised Count")

# This function makes a PCA plot
plotPCA(rld, intgroup="group")
colData(deseqdata)

resultsAB = results(deseqdata,contrast=c("group","Group A","Group B"))
resultsAC = results(deseqdata,contrast=c("group","Group A","Group C"))
resultsBC = results(deseqdata,contrast=c("group","Group B","Group C"))

# Summary of the statistics
summary(resultsAB,alpha=0.001)
summary(resultsAC,alpha=0.001)
summary(resultsBC,alpha=0.001)

# Order the results of the three groups
res_AB_sort = resultsAB[order(resultsAB$padj),]
res_AC_sort = resultsAC[order(resultsAC$padj),]
res_BC_sort = resultsBC[order(resultsBC$padj),]

#Subset the data from the three groups
res_AB_sig <- subset(res_AB_sort, alpha=0.001)
res_AC_sig <- subset(res_AC_sort, alpha=0.001)
res_BC_sig <- subset(res_BC_sort, alpha=0.001)

# write csv files 
write.csv(res_AB_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/Sig_Res_AB.csv")
write.csv(res_AC_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/Sig_Res_AC.csv")
write.csv(res_BC_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/Sig_Res_BC.csv")

# DESeq Results for genes
resAB.lfc2 <- results(deseqdata,contrast=c("group","Group A","Group B"),lfcThreshold=2)
resAC.lfc2 <- results(deseqdata,contrast=c("group","Group A","Group C"),lfcThreshold=2)
resBC.lfc2 <- results(deseqdata,contrast=c("group","Group B","Group C"),lfcThreshold=2)
resAB.lfc2_sort = resAB.lfc2[order(resAB.lfc2$padj),]
resAC.lfc2_sort = resAC.lfc2[order(resAC.lfc2$padj),]
resBC.lfc2_sort = resBC.lfc2[order(resBC.lfc2$padj),]
resAB.lfc2_sig <- subset(resAB.lfc2_sort, alpha=0.001)
resAC.lfc2_sig <- subset(resAC.lfc2_sort, alpha=0.001)
resBC.lfc2_sig <- subset(resBC.lfc2_sort, alpha=0.001)

write.csv(resAB.lfc2_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/resAB.lfc2_sig.csv")
write.csv(resAC.lfc2_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/resAC.lfc2_sig.csv")
write.csv(resBC.lfc2_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/resBC.lfc2_sig.csv")

par(mfrow =c(2,3))


# Differential expression mean of normalized counts MA plot of genes
DESeq2::plotMA(resultsAB, main="A vs B, lfc = 0")
DESeq2::plotMA(resultsAC, main ="A vs C, lfc = 0")
DESeq2::plotMA(resultsBC, main="B vs C, lfc = 0")
DESeq2::plotMA(resAB.lfc2, main="A vs B, lfc < 2")
DESeq2::plotMA(resAC.lfc2, main="A vs C, lfc < 2")
DESeq2::plotMA(resBC.lfc2, main="B vs C, lfc < 2")

# Task 2 Trancripts

#Making the data files of the experiment and constructing the table
sampleNames <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
pData=data.frame(id=sampleNames, group = rep(c("Group A","Group B","Group C"), each=4))
countTabletranscript <- read.csv("transcript_count_matrix.csv",row.names=1)
countTabletranscript = countTabletranscript[,new_order]
deseqdata <- DESeqDataSetFromMatrix(countData=countTabletranscript,colData=colTable,design=~group)                    
deseqdata

# Remove all the 0 from the gene counts
removeZerofromtranscripts <- (rowSums(counts(deseqdata))>0)
deseqdata <- deseqdata[removeZero,]
deseqdata <- DESeq(deseqdata)
head(deseqdata)
str(deseqdata)

#Making a dispersion plot
par(mfrow =c(1,1))
plotDispEsts(deseqdata)

# This function regularlises log transformation
rldt = rlog(deseqdata) # used blind to sample groups?
#log2 raw counts
lgc.rawt <- log2(counts(deseqdata, normalized = F))
# log2 normalised counts
lgc.normt <- log2(counts(deseqdata, normalized = T))

#This function makes the boxplots from the results  
par(mfrow = c(1,3))
boxplot(lgc.rawt,main= "Raw Count", )
boxplot(lgc.normt,main="Normalised Counts")
boxplot(assay(rldt),main="Regularised Counts")

# This function makes the density plots
par(mfrow = c(1,3))
plotDensities(lgc.rawt,group=deseqdata$group,legend="topright",main= "Raw Count")
plotDensities(lgc.normt,group=deseqdata$group,legend="topright",main="Normalised Count")
plotDensities(assay(rldt),group=deseqdata$group,legend="topright",main="Regularised Count")

# This function makes a PCA plot for the transcripts

plotPCA(rldt, intgroup="group")
colData(deseqdata)

resABt = results(deseqdata,contrast=c("group","Group A","Group B"))
resACt = results(deseqdata,contrast=c("group","Group A","Group C"))
resBCt = results(deseqdata,contrast=c("group","Group B","Group C"))

summary(resABt,alpha=0.001)
summary(resACt,alpha=0.001)
summary(resBCt,alpha=0.001)

resAB_sort_tra = resABt[order(resABt$padj),]
resAC_sort_tra = resACt[order(resACt$padj),]
resBC_sort_tra = resBCt[order(resBCt$padj),]

#Subset from the data
resABt_sig <- subset(resAB_sort_tra, alpha=0.001)
resACt_sig <- subset(resAC_sort_tra, alpha=0.001)
resBCt_sig <- subset(resBC_sort_tra, alpha=0.001)

# write the csv file
write.csv(resABt_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/Sig_Res_ABt_transcripts.csv")
write.csv(resACt_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/Sig_Res_ACt_transcripts.csv")
write.csv(resBCt_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/Sig_Res_BCt_transcripts.csv")

# DESeq Results for transcripts
resABt.lfc2 <- results(deseqdata,contrast=c("group","Group A","Group B"),lfcThreshold=2)
resACt.lfc2 <- results(deseqdata,contrast=c("group","Group A","Group C"),lfcThreshold=2)
resBCt.lfc2 <- results(deseqdata,contrast=c("group","Group B","Group C"),lfcThreshold=2)
resABt.lfc2_sort = resABt.lfc2[order(resABt.lfc2$padj),]
resACt.lfc2_sort = resACt.lfc2[order(resACt.lfc2$padj),]
resBCt.lfc2_sort = resBCt.lfc2[order(resBCt.lfc2$padj),]
resABt.lfc2_sig <- subset(resABt.lfc2_sort, alpha=0.001)
resACt.lfc2_sig <- subset(resACt.lfc2_sort, alpha=0.001)
resBCt.lfc2_sig <- subset(resBCt.lfc2_sort, alpha=0.001)

# write csv files
write.csv(resABt.lfc2_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/resABt.lfc2_sig.csv")
write.csv(resACt.lfc2_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/resACt.lfc2_sig.csv")
write.csv(resBCt.lfc2_sig, "/Users/salvatoreesposito/Desktop/RNAassessment/resBCt.lfc2_sig.csv")

# Differential expression mean of normalized counts MA plot of genes
par(mfrow =c(2,3))
DESeq2::plotMA(resABt, main="A vs B, lfc = 0")
DESeq2::plotMA(resACt, main ="A vs C, lfc = 0")
DESeq2::plotMA(resBCt, main="B vs C, lfc = 0")
DESeq2::plotMA(resABt.lfc2, main="A vs B, lfc < 2")
DESeq2::plotMA(resACt.lfc2, main="A vs C, lfc < 2")
DESeq2::plotMA(resBCt.lfc2, main="B vs C, lfc < 2")


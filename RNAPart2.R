
library(ggplot2)
library(amap)
library(clusterProfiler)
library(reshape2)
library(igraph)
BiocManager::install("clusterProfiler")
BiocManager::install("igraph")

# GUI 2439915


# Task 1

#Importing the files from the directory and assigning them variables 
expressionMatrix=read.table(file="EM.csv",header = TRUE, sep = "\t", row.names=1)
CC_HC=read.table(file="DE_CC_vs_HC.csv",header=TRUE,sep = "\t",row.names=1)
FLU_CC=read.table(file="DE_FLU_vs_CC.csv",header=TRUE,sep = "\t", row.names=1)
FLU_HC=read.table(file="DE_FLU_vs_HC.csv",header=TRUE,sep = "\t", row.names=1)
sample_sheet=read.table(file="sample_sheet.csv", header=TRUE, sep = "\t", row.names=1)
annotation = read.table("mart_export2.txt", sep = '\t', header = TRUE, row.names=1)

# This function merges the expression matrix and saves it as an csv file
annotationExpressionmatrix = merge(expressionMatrix, annotation, by = 0)
row.names(annotationExpressionmatrix)<-annotationExpressionmatrix$Row.names
annotationExpressionmatrix<-annotationExpressionmatrix[,-1]
write.table(annotationExpressionmatrix,
            "annotationExpressionmatrix.csv")

# this function introduces the mean to the annotation expression matrix
annotationExpressionmatrix$mean <- apply(annotationExpressionmatrix[, c(1:9)], 1, mean)

# this function adds columns to the annotation expression matrix by merging them
addcolumns<-function(data){ data <- merge(data, annotation, by = 0)[, c(1:4,9)]
  row.names(data) <- data$Row.names
  data <- data[, -1]
  # this function adds column for the -log10p
  data$log10p <- -log(data$p, 10)
  # this function adds significant columns
  data$significant <- as.factor(data$p.adj < 0.05 & ((abs(data$log2fold)) > 1.0))
  return(data)
}
  # datas from the three files are merged together
  CC_HC<-addcolumns(CC_HC)
  FLU_CC<-addcolumns(FLU_CC)
  FLU_HC<-addcolumns(FLU_HC)
  #this function renames the title of each column as shown below
  renamecolumns<-function(data,name){
  colnames(data) <- c(paste(name,"_log2fold", sep = ""), paste(name,"_p", sep = ""), paste(name,"_p_adjust", sep = ""),"Name",
  paste(name,"_log10p", sep = ""), paste(name,"_significant", sep = ""))
  data <- data[, c(4, 1:3,5:6)]
  return(data)
  }
  
  #this function renames the title of each column as shown below
  CC_HC<-renamecolumns(CC_HC,"CC_HC")
  FLU_CC<-renamecolumns(FLU_CC,"FLU_CC")
  FLU_HC<-renamecolumns(FLU_HC,"FLU_HC")
  
  # the tables from the three groups of data are merged into one table file called mastertable
  maintable <- merge(CC_HC, FLU_CC[2:6], by = 0)
  row.names(maintable) <- maintable$Row.names
  maintable <- maintable[, -1]
  # the tables from the three groups of data are merged into one table file called mastertable
  maintable <- merge(maintable, FLU_HC[2:6], by = 0)
  row.names(maintable) <- maintable$Row.names
  maintable <- maintable[, -1]
  # the tables from the three groups of data are merged into one table file called mastertable
  maintable <- merge(maintable, annotationExpressionmatrix[, c(1:46, 48:49)], by =0)
  row.names(maintable) <- maintable$Row.names
  maintable <- maintable[, -1]
  
  # This function makes a Volcano plot using the data from the mastertable
  CC_HC_volcanoPlot <- ggplot(maintable, aes(x=CC_HC_log2fold, 
  y=CC_HC_log10p, group=CC_HC_significant, color=CC_HC_significant )) + geom_point()
  # the volcano plot is saved as a png file
  png("Volcano_CC_HC.png") 
  CC_HC_volcanoPlot 
  dev.off()
  # This function makes a Volcano plot using the data from the mastertable
  FLU_CC_volcanoPlot <- ggplot(maintable, aes(x=FLU_CC_log2fold, 
   y=FLU_CC_log10p, group=FLU_CC_significant, color=FLU_CC_significant )) + geom_point()
  # the volcano plot is saved as a png file
  png("Volcano_FLU_CC.png")
  FLU_CC_volcanoPlot
  dev.off()
  # This function makes a Volcano plot using the data from the mastertable
  FLU_HC_volcanoPlot <- ggplot(maintable, aes(x=FLU_HC_log2fold, 
  y=FLU_HC_log10p, group=FLU_HC_significant, color=FLU_HC_significant )) + geom_point()
  # the volcano plot is saved as a png file
  png("Volcano_FLU_HC.png")
  FLU_HC_volcanoPlot
  dev.off()
  
  # This function makes an MAplot from the CC_HC data
  MAplot_CC_HC <- ggplot(maintable, aes(y = CC_HC_log2fold,
  x = log10(mean),colour = CC_HC_significant, group = CC_HC_significant)
  ) + geom_point() + theme_classic() + scale_color_manual(values = colours) +
  labs(title = "MAplot_CC_HC")
  # the MAplot is saved as a png file
  png("MAplot_CC_HC.png")
  MAplot_CC_HC
  dev.off()
  
  # This function makes an MAplot from the FLU_CC data
  MAplot_FLU_CC <- ggplot(maintable, aes(y = FLU_CC_log2fold,
  x = log10(mean), colour = FLU_CC_significant, group = FLU_CC_significant )
  ) + geom_point() + theme_classic() + scale_color_manual(values = colours) +
  labs(title = "MAplot_FLU_CC") 
  # the MAplot is saved as a png file
  png("MAplot_FLU_CC.png")
  MAplot_FLU_CC
  dev.off()
    
  # This function makes an MAplot from the FLU_hC data
  MAplot_FLU_HC <- ggplot(maintable,aes( y = FLU_HC_log2fold, 
  x = log10(mean),colour = FLU_HC_significant,group = FLU_HC_significant)
  ) + geom_point() + theme_classic() + scale_color_manual(values = colours) +
  labs(title = "MAplot_FLU_HC")
  # the MAplot is saved as a png file
  png("MAplot_FLU_HC.png")
  MAplot_FLU_HC
  dev.off()
  
  # This function makes PCA plots
  # colours are chosen for the plots
  colours3 = c("lightblue4", "lightpink4", "palegreen4")
  # this function makes numeric matrix for pca
  expressionnumericMatrix <- as.matrix(sapply(expressionMatrix, as.numeric))
  class(expressionnumericMatrix)
  #the pca analysis function is used to reduce dimensionality of the matrix
  pca <- prcomp(t(expressionnumericMatrix))
  #get coordinated of pc comparisons
  pca_coordinates <- data.frame(pca$x)
  
  # First PCA plot without carrying normalisation
  pca_plot <-ggplot(pca_coordinates,
  aes(x = PC1,y = PC2,colour = sample_sheet$GROUP,
  group = sample_sheet$GROUP,fill = sample_sheet$GROUP)
  ) +geom_point(size = 5) + labs(title = "PCA plot")
  pca_plot
  # this function saves the pca plot as a png file
  png("PCA1-2.png")
  pca_plot
  dev.off()

  # The second PCA plot uses z-scores to carry normalization of the data
  z_scores <- data.frame(t(scale(t(expressionMatrix))))
  z_expressionnumericMatrix <- as.matrix(sapply(z_scores, as.numeric))
  z_expressionnumericMatrix <- na.omit(z_expressionnumericMatrix)
  #PCA analysis of the matrix
  z_pca <- prcomp(t(z_expressionnumericMatrix))
  # this function gets the coordinates of the pc comparisons
  z_pca_coordinates <- data.frame(z_pca$x)
  
  z_pca_plot <- ggplot(z_pca_coordinates,
  aes(x = PC1,y = PC2,colour = sample_sheet$GROUP,
  group = sample_sheet$GROUP,fill = sample_sheet$GROUP)
  ) +geom_point(size = 5) + labs(title = "Z-scores PCA plot")
  z_pca_plot
  # this function saves the z-pca plot as a png file
  png("Z-score_PCA1-2.png")
  z_pca_plot
  dev.off()
  
  # Function to make the heatmap to see the similarity of gene expression
  #correlation coordinates between all samples
  cor_heatmap <- cor(z_expressionnumericMatrix, method = "spearman")
  #melt so all values in one column
  melted <- melt(cor_heatmap)
  colours_heatmap <- c("red", "white", "darkblue")
  correlation_heatmap <- ggplot(melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + scale_fill_gradientn(limits = c(-1, 1),
  colours = colorRampPalette(colours_heatmap)(100))
  correlation_heatmap
  # save correlation heatmap
  png("mastertable_correlation_heatmap.png")
  correlation_heatmap
  dev.off()
  
  # Full expression HeatMap 
  HM_gene_expression <-subset(maintable, ((
  abs(maintable$FLU_CC_log2fold) > 1 & maintable$FLU_CC_significant == TRUE)|
  (abs(maintable$FLU_HC_log2fold) > 1 & maintable$FLU_HC_significant == TRUE)|
  (abs(maintable$CC_HC_log2fold) > 1 & maintable$CC_HC_significant == TRUE)
  ))
  
  HM_expression <- expressionMatrix[row.names(HM_gene_expression),]
  HM_expression$Gene <- row.names(HM_expression)
  HM_melted <- melt(HM_expression, na.rm = T)
  
  # z-scores normalization
  HMz_expression <- data.frame(t(scale(t(HM_expression[,1:42]))))
  HMz_nexpression <- HMz_expression
  HMz_nexpression$Gene <- row.names(HMz_expression)
  HMz_melted <- melt(HMz_nexpression, na.rm = T)

  HMz_heatmap <- ggplot(HMz_melted, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() + theme_classic() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colours = colorRampPalette(colours_heatmap)(100)) +
  labs(x = "Sample", y = "Gene", title = "Gene expression heat map")
  HMz_heatmap
  # save heatmap as a png file
  png("z-score_gene_expression_heatmap.png")
  HMz_heatmap
  dev.off()
  
  #measure distances
  y.dist <- Dist(HMz_expression, method = "spearman")
  #use mean agglomeration
  y.cluster <- hclust(y.dist, method = "average")
  #extract dendrogram
  y.dd <- as.dendrogram(y.cluster)
  #reorder dendrogram
  y.dd.reordered <- reorder(y.dd, 0, FUN = "average")
  #re-order data (still need to order ggplot by levels when plotting)
  y.order <- order.dendrogram(y.dd.reordered)
  #change order of genes
  HMz_expression_ordered <- HMz_expression[y.order, ]
  #make new table to add Gene name column
  HMz_nexpression_ordered <- HMz_expression_ordered
  #add gene name column
  HMz_nexpression_ordered$Gene <- row.names(HMz_nexpression_ordered)
  #melt to use with ggplot
  HMz_ordered_melted <- melt(HMz_nexpression_ordered, na.rm = T)
  
  HMz_ordered_heatmap <- ggplot(HMz_ordered_melted, aes(
  x = variable, y = factor(Gene, level = HMz_nexpression_ordered$Gene),
  fill = value)) + geom_tile() + theme_classic() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colours = colorRampPalette(col_heatmap)(100)) +
  labs(x = "Sample", y = "Gene", title = " cluster y gene expression heat map")
  HMz_ordered_heatmap
  # save ordered expresssion of z-score heatmap
  png("z-score_ordered_exp_heatmap.png")
  HMz_ordered_heatmap
  dev.off()
  #order x axis
  # must first transform the data
  tHMz_expression <- t(HMz_expression)
  #measure distances
  x.dist <- Dist(tHMz_expression, method = "spearman")
  #use mean agglomeration
  x.cluster <- hclust(x.dist, method = "average")
  #extract dendrogram
  x.dd <- as.dendrogram(x.cluster)
  #reorder dendrogram
  x.dd.reordered <- reorder(x.dd, 0, FUN = "average")
  #re-order data (still need to order ggplot by levels when plotting)
  x.order <- order.dendrogram(x.dd.reordered)
  #change order of genes
  HMz_expression_ordered_xy <- HMz_expression_ordered[, x.order]
  #make new table to add Gene name column
  HMz_nexpression_ordered_xy <- HMz_expression_ordered_xy
  #add gene name column
  HMz_nexpression_ordered_xy$Gene <-
    row.names(HMz_nexpression_ordered_xy)
  #melt to use with ggplot
  HMz_ordered_xy_melted <- melt(HMz_nexpression_ordered_xy, na.rm = T)
  
  HMz_ordered_xy_heatmap <- ggplot(HMz_ordered_melted, aes(
  x = factor(variable, level = colnames(HMz_nexpression_ordered_xy)),
  y = factor(Gene, level = HMz_nexpression_ordered$Gene), fill = value
  )) + geom_tile() + theme_classic() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradientn(colours = colorRampPalette(col_heatmap)(100)) + 
  labs(x = "Sample", y = "Gene", title = "Gene expression heat map")
  HMz_ordered_xy_heatmap
  
  png("z-score_ordered_xy_geneexpression_heatmap.png")
  HMz_ordered_xy_heatmap
  dev.off()
  
  
  # Genontology analysis
  
  library(DOSE)
  library(org.Hs.eg.db)
  BiocManager::install("org.Hs.eg.db")
  
  # analysis for CC_HC data and generation of barplot
  G_O=subset(maintable, CC_HC_significant == "TRUE") 
  gen_ont=row.names(G_O)
  gene.df <- bitr(gen_ont, fromType = "ENSEMBL",
  toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  Gene_entre=gene.df$ENTREZID 
  edo <- enrichDGN(Gene_entre)
  barplot(edo, showCategory=20)
  # save barplot
  png("CC_HC_barplot.png")
  barplot(edo, showCategory=20)
  dev.off()
  
  
  # analysis for FLU_HC data and generation of barplot
  G_O_flu_hc=subset(maintable, FLU_HC_significant == "TRUE") 
  gen_ont_flu_hc=row.names(G_O_flu_hc)
  gene.df_flu_hc <- bitr(gen_ont_flu_hc, fromType = "ENSEMBL",
  toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  Gene_entre_flu_hc=gene.df_flu_hc$ENTREZID 
  edo_flu_hc <- enrichDGN(Gene_entre_flu_hc)
  barplot(edo_flu_hc, showCategory=20)
  # save barplot
  png("FLU_HC_barplot.png")
  barplot(edo_flu_hc, showCategory=20)
  dev.off()
  
  # analysis for FLU_CC data and generation of barplot
  G_O_flu_cc=subset(maintable, FLU_CC_significant == "TRUE") 
  gen_ont_flu_cc=row.names(G_O_flu_cc)
  gene.df_flu_cc <- bitr(gen_ont_flu_cc, fromType = "ENSEMBL",
  toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  Gene_entre_flu_cc=gene.df_flu_cc$ENTREZID 
  edo_flu_cc <- enrichDGN(Gene_entre_flu_cc)
  barplot(edo_flu_cc, showCategory=20)
  # save barplot
  png("FLU_CC_barplot.png")
  barplot(edo_flu_cc, showCategory=20)
  dev.off()
  
  
  
  
  

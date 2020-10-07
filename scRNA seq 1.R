


library("ggplot2")
library("reshape2")
library("gridExtra")
library("amap")
library("Seurat")
library("scran")
library("dplyr")
library("magrittr")
library("clusterProfiler")
library("igraph")
library("DOSE")
library("org.Hs.eg.db")
library("monocle")
library("readr")
library("readxl")
library("tidyverse")
library("M3Drop")


#TASK2

# expression table
expression = read.table("Expression.txt",header = TRUE)
# subset the expression for gene trem2
trem2 = subset(expression,Gene == "TREM2")
# make a boxplot
BoxPlot <- ggplot(trem2, aes(x=Condition, y=-log10(Expression))) + geom_boxplot()
BoxPlot
# make a violin plot
ViolinPlot = ggplot(trem2, aes(x=Condition, y=log(Expression))) + geom_violin(trim=FALSE, fill="lightblue") + geom_jitter(shape=16, position=position_jitter(0.05))+ stat_summary(color="Red", geom="pointrange")
ViolinPlot
# condition for tremRA and remission
trem2_A<-trem2[which(trem2$Condition=="ActiveRA"),]
trem2_R<-trem2[which(trem2$Condition=="Remission"),]
# make heatmap
heatmap_A <-ggplot(trem2_A,aes(x = row.names(trem2_A),y = Gene,fill = Expression)) +
geom_tile() + theme_classic()
heatmap_R <-ggplot(trem2_R,aes(x = row.names(trem2_R),y = Gene,fill = Expression)) +
geom_tile() + theme_classic()
CombinePlots(list(heatmap_A,heatmap_R))
t.test(trem2$Expression~trem2$Condition)

#TASK 3
# import ran and rem data
acRAN.data <- Read10X(data.dir = "acRAN_Seurat/GRCh38/")
REM.data <- Read10X(data.dir = "REM_Seurat/GRCh38/")
# create them as a seurat object
acRAN.data <- CreateSeuratObject(counts = acRAN.data, project = "acRAN.data", min.cells = 3, min.features = 200)
REM.data <- CreateSeuratObject(counts = REM.data, project = "REM.data", min.cells = 3, min.features = 200)
acRAN.data$sample <- "acRAN.data"
acRAN.data$group <- "activeRA"

REM.data$sample <- "REM.data"
REM.data$group <- "remission"

# Initiate the QC statistical analysis
acRAN.data[["percent.mt"]] <- PercentageFeatureSet(acRAN.data, pattern = "^MT-")
REM.data[["percent.mt"]] <- PercentageFeatureSet(REM.data, pattern = "^MT-")
# make a scatter plot and combine them for RAN
plot1 <- FeatureScatter(object = acRAN.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = acRAN.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# make a scatter plot and combine them for REM
plot3 <- FeatureScatter(object = REM.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(object = REM.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))
# subset the data for both including nfeature RNA and percent mt
acRAN.data <- subset(acRAN.data, subset = nFeature_RNA > 700 & nFeature_RNA < 3500 & percent.mt < 18)
REM.data <- subset(REM.data, subset = nFeature_RNA > 700 & nFeature_RNA < 3500 & percent.mt < 18)
# clear the data with the SCTransform
acRAN.data <- SCTransform(acRAN.data, vars.to.regress = "percent.mt", verbose = FALSE)
REM.data <- SCTransform(REM.data, vars.to.regress = "percent.mt", verbose = FALSE)
# PCA analysis
acRAN.data <- RunPCA(acRAN.data, features = VariableFeatures(object = acRAN.data))                 
REM.data <- RunPCA(REM.data, features = VariableFeatures(object = REM.data))                 
# Add elbow plot
ElbowPlot(object = acRAN.data, ndims = 50)                
ElbowPlot(object = REM.data, ndims = 50)   
# cluster the genes using KNN
dim = 17
TC.anchors <- FindIntegrationAnchors(object.list = list(acRAN.data,REM.data), dims = 1:dim)
Combined <- IntegrateData(anchorset = TC.anchors, dims = 1:dim)
DefaultAssay(Combined) <- "integrated"
Combined <- ScaleData(Combined, verbose = FALSE)
Combined <- RunPCA(Combined, npcs = dim, verbose = FALSE)
Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
Combined <- FindClusters(Combined, resolution = 0.1)
p1 <- DimPlot(Combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Combined, reduction = "umap", label = TRUE)
p1+p2
# make a feature plot to visualise the dataset
Combined[["UMI"]] <-  Combined$nCount_RNA 
Combined[["genes"]] <-  Combined$nFeature_RNA
FeaturePlot(Combined, features = "UMI")
# save the objects due to their long running time
saveRDS(acRAN.data, file = "acRA.rds")
saveRDS(REM.data, file = "REM.rds")
saveRDS(Combined, file = "Combined.rds")
# upload the objects
acRA = readRDS(file = "acRA.rds")
REM = readRDS(file = "REM.rds")
Combined = readRDS(file = "Combined.rds")
# find the markers and group them by clusters
Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# top20 and top 10 markers
top20 <- Combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top10 <- Combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# name the clusters by researching the articles
New.cluster.ids <- c("Macrophages_01", "Fibroblasts_CXCL12", "Fibroblasts_CD55", "TCells_IL32", "Platelets","T cell","Fibroblasts_PDPN","Macrophage_CD163","Mast cellTPSAB1","B Cell","Plasma Cell","Dendritic cell")
names(New.cluster.ids) <- levels(Combined)
Combined <- RenameIdents(Combined, New.cluster.ids)
# heatmap
DoHeatmap(Combined, features = top10$gene) + NoLegend()

View(top10)
View(top20)
# Subset Macrophages
All.macrophages <- subset(Combined, idents = c("Macrophages_01","Macrophage_CD163"), invert = FALSE)
DimPlot(All.macrophages, reduction = "umap")
# PCA analysis 
All.macrophages <- RunPCA(All.macrophages, features = VariableFeatures(object = All.macrophages))
plot_mac<-ElbowPlot(object = All.macrophages, ndims = 40)
#Differential Expression
All.macrophages$celltype.group <- paste(Idents(All.macrophages), All.macrophages$group, sep = "_")
All.macrophages$celltype <- Idents(All.macrophages)
DefaultAssay(object = All.macrophages) <- "integrated"
# Run the standard workflow for visualization and clustering
All.macrophages <- ScaleData(object = All.macrophages, verbose = FALSE)
All.macrophages <- RunPCA(object = All.macrophages, npcs = dim, verbose = FALSE)
# t-SNE and Clustering
All.macrophages <- RunUMAP(object = All.macrophages, reduction = "pca", dims = 1:dim)
All.macrophages <- FindNeighbors(object = All.macrophages, reduction = "pca", dims = 1:dim)
All.macrophages <- FindClusters(All.macrophages, resolution = 0.15)
# make a dimplot of the macrophages subset
DimPlot(All.macrophages, reduction = "umap", split.by = "sample")
DimPlot(object = All.macrophages, reduction = "umap")
Combined.markers02 <- FindAllMarkers(All.macrophages, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.markers02 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
All.macrophages$celltype.group <- paste(Idents(All.macrophages), All.macrophages$group, sep = "_")
All.macrophages$celltype <- Idents(All.macrophages)
# upregulated genes in remission found=0
n="0"
check_sanity = FALSE
# differential expression and make table of p-values for genes found
Idents(All.macrophages) <- "celltype.group"
REM.upregulated <- FindMarkers(All.macrophages, ident.1 = paste(n,"_activeRA", sep=""), ident.2 = paste(n,"_remission", sep=""), verbose = FALSE, min.pct = 0.25)
write.table(REM.upregulated,paste("Table of P-values.",n,".txt", sep="")) 


# Gene Ontology analysis
tablepvalue=read.table("Table of P-values.0.txt")
nameofthegenes=subset(tablepvalue, tablepvalue$avg_logFC>0 & tablepvalue$p_val_adj<0.01)
genes <- bitr(row.names(nameofthegenes), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
edo <- enrichDGN(genes[,2])
plot_GeneOntology<-ggplot(edo[c(1:20)],aes(x=Count,y=reorder(Description,-p.adjust),fill=p.adjust))+geom_col()+
  scale_fill_gradientn(trans='reverse',colours = colorRampPalette(c("blue","red"))(100))+
  labs(x = "Gene count", y = "", title = "GO enrichment bases on 3 upregulated genes")+theme_classic()
plot_GeneOntology

# Pseudo time-analysis Pre-processing

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(All.macrophages@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = All.macrophages@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct macrop cds
macrop_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,
expressionFamily = negbinomial.size())

filtered_only <- M3DropCleanData(data, is.counts=FALSE,min_detected_genes=1000)
cells <- colnames(as.matrix(filtered_only$data))
pd <- pd[cells,]
macrop_cds <- newCellDataSet(as.matrix(filtered_only$data), phenoData = pd)

#Cluster Data
# ESTIMATE THE SIZE FACTORS AND ESTIMATE THE DISPERSIONS
macrop_cds <- estimateSizeFactors(macrop_cds)
macrop_cds <- estimateDispersions(macrop_cds)
pData(macrop_cds)
# DETECT GENES THAT ARE EXPRESSED ABOVE 0.1 10+ CELLS
macrop_cds <- detectGenes(macrop_cds, min_expr = 0.1)
print(head(fData(macrop_cds)))
genes_expressed <- row.names(subset(fData(macrop_cds),num_cells_expressed >= 10))
# SELECTION GENES OF ADEQUATE EXPRESSION AND MARKED THEM AS ORDERING GENES.
dispersion_table <- dispersionTable(macrop_cds)
head(dispersion_table)
unsup_ordering_genes <- subset(dispersion_table, mean_expression >= 0.1)
macrop_cds <- setOrderingFilter(macrop_cds, unsup_ordering_genes$gene_id)
print(head(fData(macrop_cds)))
plot_ordering_genes(macrop_cds)
# Genes Identified
adequate_expression_genes <- row.names(subset(fData(macrop_cds), use_for_ordering == TRUE))
4962
# GENERATE A SCREE PLOT AND REDUCE DIMENSIONALITY
plot_pc_variance_explained(macrop_cds, return_all = F) 
macrop_cds <- reduceDimension(macrop_cds, max_components = 2, num_dim = 7, reduction_method = 'tSNE', verbose = F)
# CLUSTER THE CELLS AND THEN VISUALISE THIS CLUSTERING.
macrop_cds <- clusterCells(macrop_cds)
plot_cell_clusters(macrop_cds)

# IDENTIFY GENES THAT HAVE ADEQUATE EXPRESSION (>=0.5) AND GREATER DISPERSION THAN EXPECTED BY THE MODEL
ordering_genes <- subset(dispersion_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
macrop_cds <- setOrderingFilter(macrop_cds, ordering_genes)
print(head(fData(macrop_cds)))
# Genes Identified
adequate_expression_genes <- row.names(subset(fData(macrop_cds), use_for_ordering == TRUE))
2009
# PSEUDOTIME ANALYSIS
macrop_cds <- reduceDimension(macrop_cds, max_components = 2, num_dim = 7,method = 'tSNE')
macrop_cds <- orderCells(macrop_cds)
plot_cell_trajectory(macrop_cds)
# MODIFY THE plot_cell_trajectory() CALL TO COLOUR BY lopez
plot_cell_trajectory(macrop_cds, color_by = "group")
# MODIFY THE plot_cell_trajectory() CALL TO COLOUR BY THE CLUSTER
plot_cell_trajectory(macrop_cds, color_by = "Cluster")
# IDENTIFICATION OF A SMALLER GENE SET WITH HIGH EXPRESSION AND HIGH DISPERSION
high_expressed_genes <- subset(dispersion_table, mean_expression >= 10 & dispersion_empirical >= 5 * dispersion_fit)$gene_id
high_expressed_genes_ids <- as.array(high_expressed_genes)
differential_test_res <- differentialGeneTest(macrop_cds[high_expressed_genes_ids,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
gene_names <- row.names(subset(differential_test_res, qval < 0.1))
plot_pseudotime_heatmap(macrop_cds[gene_names,], num_clusters = 6, cores = 1, show_rownames = T)    
# PLOT THE GENE EXPRESSION AGAINST PSEUDOTIME FOR THE TOP 3 VARIABLE GENES
top3 =head(dispersion_table[order(-dispersion_table[,2]),],n=3)$gene_id
to_be_tested <-top3
cds_subset <- macrop_cds[to_be_tested,]
plot_genes_in_pseudotime(cds_subset, color_by = "Cluster")
# MODEL EACH GENE'S EXPRESSION ACROSS PSEUDOTIME
differential_test_res<- differentialGeneTest(macrop_cds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
differential_test_res <- differential_test_res[order(differential_test_res$qval), ]
genes <- row.names(differential_test_res)
top_200_genes <-  row.names(head(differential_test_res, 200))

plot_pseudotime_heatmap(macrop_cds[top_200_genes,], num_clusters = 6, cores = 1,
show_rownames = T)


                  
#Setting up environment.
library(tidyverse)
library(Seurat)
path <- "~/Desktop/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"
setwd("~/Desktop/Lung Cells")
list.files(path)
set.seed(42)

#Data import
nsclc<-Read10X_h5(filename = "~/Desktop/Lung Cells/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc)
cts<-nsclc$`Gene Expression`
#cts is a large sparse matrix, representing zeroes as dots. 

#Creating a Seurat Object (raw data).
nsclc_seu<-CreateSeuratObject(counts=cts, project = "NSCLC",min.cells = 3, min.features = 200 )
#min.cells = only genes that are expressed in at least 3 cells will be kept in the object. min. features= only cells expresseing this amount of genes will be kept in the object. 

#Quality control.
#Analyzing percent mitochondrial genes (distinguishes low-quality/ dying cells)
nsclc_seu[['percent_mt']]<-PercentageFeatureSet(nsclc_seu, pattern = '^MT-')
View(nsclc_seu@meta.data)
VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(nsclc_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
#nFeature_RNA = number of unique genes detected in a cell. ncount_RNA= number of molecules detected in the cell. 

#Filtering
nsclc_seu<- subset(nsclc_seu, nFeature_RNA > 200 & nFeature_RNA < 1000 & percent_mt < 5)
#Changed filtering restrictions from cells with <2500 nFeature_RNA to <1000 due to computer memory limits. 

#Normalization
nsclc_seu <- NormalizeData(nsclc_seu, normalization.method = 'LogNormalize', scale.factor = 10000)
nsclc_seu <- NormalizeData(nsclc_seu)
str(nsclc_seu)

#Identifying highly variable genes. 
nsclc_seu <- FindVariableFeatures(nsclc_seu, selection.method =  'vst', nfeatures = 2000)
#This will allow us to distinguish different cell types. 

#Identify the top 10 highly variable genes. 
top10 <- head(VariableFeatures(nsclc_seu), 10)
top10_plot <- VariableFeaturePlot(nsclc_seu)
LabelPoints(plot = top10_plot, points = top10, repel = TRUE)

#Scaling.
#We are going to scale data to get rid of unwanted variation. 
all_gene<- rownames(nsclc_seu)
nsclc_seu <- ScaleData(nsclc_seu, features = all_gene)
View(nsclc_seu@assays$RNA)

#PCA 
nsclc_seu<-RunPCA(nsclc_seu, features = VariableFeatures(nsclc_seu))
print(nsclc_seu[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc_seu, dims = 1, cells = 500, balanced = TRUE)
DimPlot(nsclc_seu, reduction = "pca") + NoLegend()

#Analyze dimensionality of the data.
ElbowPlot(nsclc_seu)

#Clustering.
nsclc_seu <- FindNeighbors(nsclc_seu, dims = 1:10)
#Changed dimensions from 1,15 to 1,10 based on reduced variation from smaller subset of cells. 
nsclc_seu <- FindClusters(nsclc_seu, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(nsclc_seu@meta.data)                          
   
DimPlot(nsclc_seu, group.by = 'RNA_snn_res.0.3', label = TRUE)                       
Idents(nsclc_seu) <- 'RNA_snn_res0.1' # set identity of clusters
#Low resolution due to reduced number of cells from data filtering. 

#UMAP
nsclc_seu <- RunUMAP(nsclc_seu, dims = 1:10)
DimPlot(nsclc_seu, reduction = 'umap')

saveRDS(nsclc_seu, file = 'nsclc_seu.RDS')

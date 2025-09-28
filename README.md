# Single Cell RNA Sequencing Preprocessing with Seurat (Biostats Squid Tutorial)

## Note:
Coding and workflow outlined below are from Biostats Squid's "Standard scRNAseq pre-processing workflow with Seurat" tutorial (https://biostatsquid.com/scrnaseq-preprocessing-workflow-seurat/#step1).
## Goal: 
Transform raw single cell RNA sequencing data into preprocessed data ready for downstream analysis. 

## Methods: 
Raw data was acquired from the 10x Genomics's " 40k Mixture of NSCLC DTCs from 7 donors, 3' HT v3.1 - Gene Expression - Feature / cell matrix HDF5 (raw)"  dataset. Upon data import, the Rstudio software was used to perform preproccessing functions from the Seurat and tidyverse packages. Violin plots and linear regression models were constructed for visualization of percent mitochondrial DNA and abnormal gene expression patterns. Cells exhibiting <5% mitochondrial DNA and 200-1000 detectable genes were retained. Due to computational constraints, less cells were retained than suggested from the tutorial. Data was then normalized using a logarithmic scale of a factor of 10,000. Highly variable genes were identified and used to guide Principal Component Analysis (PCA). A scree plot was then utilized to determine dimensionality of the data. Data was clustered with a resolution parameter of 0.3, and projected onto a UMAP.

## Results: 
Groups generated from clustering overlapped heavily, with unclear spatial boundaries for some clusters. The UMAP failed to generate subgroups, labelling all data points under the same cluster. 

Quality control: Identifying low-quality/dying cells.
![](https://github.com/sarutoor2002/scRNA-seq-Lung-Cell-Preprocessing/blob/main/QC%20Violin%20Plot.png)

![](https://github.com/sarutoor2002/scRNA-seq-Lung-Cell-Preprocessing/blob/main/Regression.png)

Identification of top ten most highly variable genes. 
![](https://github.com/sarutoor2002/scRNA-seq-Lung-Cell-Preprocessing/blob/main/Highly%20Variable%20Genes.png)

Principal Component Analysis.

Determining dimensionality with scree plot.
![](https://github.com/sarutoor2002/scRNA-seq-Lung-Cell-Preprocessing/blob/main/Scree%20Plot.png)

Clustering.
![](https://github.com/sarutoor2002/scRNA-seq-Lung-Cell-Preprocessing/blob/main/Clustering.png)

UMAP.
![]([https://github.com/sarutoor2002/scRNA-seq-Lung-Cell-Preprocessing/blob/main/UMAP.png)

## Conclusion: 
Lower cell counts were retained in the filtered dataset due to computational constraints. As a result, reduced cell counts likely contributed to lower levels of detectable variation, thereby imparing clustering ability and UMAP synthesis. To enhance data preprocessing, computers with a higher RAM should be considered.    

## Citation:
Biostatsquid. (2025, March 28). Standard scrnaseq pre-processing workflow with Seurat. biostatsquid.com - Easy bioinformatics and biostatistics. https://biostatsquid.com/scrnaseq-preprocessing-workflow-seurat/#step1 

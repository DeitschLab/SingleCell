#This code generates the seurat object and all the figures for the HIVE experiment 
#Figure 5 and Supplementary Figure 2

setwd("") #Set to the folder where all scripts and files are found

library(Seurat) #version 4.4.0
library(dplyr)
library(patchwork)
library(Matrix) #version 1.6.3
library(data.table)
library(scater)
library(ape)
library(viridis)
library(scCustomize)
lapply(list("umap","parallel","ggplot2","magrittr","reshape2","pheatmap","scales","splines","zoo","RColorBrewer"),require,character.only = TRUE)

#Upload gene lists used in the code
multicopy <- read.table("multicopy.txt")
multicopy <- c(multicopy[,1])
multicopy <- gsub("_", "-", multicopy)

#Read count matrices to make Seurat Object (done for every sample)
matHSA <- fread("HighSingleAHIVE.tsv")
setDF(matHSA)
rownames(matHSA) <- matHSA[[1]]
matHSA[[1]] <- NULL
#make a Seurat object
HSA_Seurat <- CreateSeuratObject(counts = matHSA)

matHSB <- fread("HighSingleBHIVE.tsv")
setDF(matHSB)
rownames(matHSB) <- matHSB[[1]]
matHSB[[1]] <- NULL
#make a Seurat object
HSB_Seurat <- CreateSeuratObject(counts = matHSB)

matLMA <- fread("LowManyAHIVE.tsv")
setDF(matLMA)
rownames(matLMA) <- matLMA[[1]]
matLMA[[1]] <- NULL
#make a Seurat object
LMA_Seurat <- CreateSeuratObject(counts = matLMA)

matLMB <- fread("LowManyBHIVE.tsv")
setDF(matLMB)
rownames(matLMB) <- matLMB[[1]]
matLMB[[1]] <- NULL
#make a Seurat object
LMB_Seurat <- CreateSeuratObject(counts = matLMB)

#Merge Seurat objects
MergeAB <- merge(x=HSA_Seurat, y=HSB_Seurat)
MergeABA <- merge(x=MergeAB, y=LMA_Seurat)
MergeAll <- merge(x=MergeABA, y=LMB_Seurat)

#PCA, UMAP and Plot
MergeAll <- SCTransform(MergeAll, verbose = T)
MergeAll <- RunPCA(MergeAll, verbose = T)
print(MergeAll[["pca"]], dims = 1:12, nfeatures = 5)

MergeAll <- RunUMAP(MergeAll, dims = 1:12, verbose = T)
MergeAll <- FindNeighbors(MergeAll, verbose = T)
MergeAll <- FindClusters(MergeAll, verbose = T, resolution = 0.5)

DimPlot(MergeAll, label = TRUE, label.box = TRUE, label.size=4)
DimPlot(MergeAll, reduction = "umap", group.by = "orig.ident")

#Feature plots in Figure 5F and 5G
pal <- viridis(n = 10, option = "C")
FeaturePlot_scCustom(MergeAll, features = "PF3D7-0711700", cols = pal)

#Identify Cluster Markers
AllMarkers <- FindAllMarkers(MergeAll, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25)
PositiveMarkers <- FindAllMarkers(MergeAll, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(PositiveMarkers, "PositiveClusterMarkers.csv")

#ViolinPlot in Supplementary Figure 2A
top_genes <- Extract_Top_Markers(AllMarkers, num_genes = 1,
                                 group_by = "cluster", rank_by = "avg_log2FC", named_vector = FALSE)
varibow_pal <- DiscretePalette_scCustomize(num_colors = 11, palette = "varibow")
Stacked_VlnPlot(MergeAll, features = top_genes, colors_use=varibow_pal)

#Clusters composition by sample in Supplememtary Figure 2B and 2C
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")
library(dittoSeq)

dittoBarPlot(MergeAll, var="orig.ident", group.by ="seurat_clusters", color.panel = c("#F8766D","#7CAE00","#00BFC4","#C77CFF"), scale = "percent", x.reorder=c(1,2,4,5,6,7,8,9,10,11,3))
dittoBarPlot(MergeAll, var="orig.ident", group.by ="seurat_clusters", color.panel = c("#F8766D","#7CAE00","#00BFC4","#C77CFF"), scale = "count", x.reorder=c(1,2,4,5,6,7,8,9,10,11,3))

#to remove clonally variant genes in the Seurat object for Figure 5E
MergeAll_MC <- MergeAll[!rownames(MergeAll) %in% multicopy,]

MergeAll_MC <- SCTransform(MergeAll_MC, verbose = T)
MergeAll_MC <- RunPCA(MergeAll_MC, verbose = T)
print(MergeAll_MC[["pca"]], dims = 1:12, nfeatures = 5)

MergeAll_MC <- RunUMAP(MergeAll_MC, dims = 1:12, verbose = T)
MergeAll_MC <- FindNeighbors(MergeAll_MC, verbose = T)
MergeAll_MC <- FindClusters(MergeAll_MC, verbose = T, resolution = 0.5)

DimPlot(MergeAll_MC, reduction = "umap", group.by = "orig.ident")
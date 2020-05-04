library(Seurat)##Seurat V2.3.4 was used
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(DoubletFinder)
rm(list = ls())

ags <- commandArgs(trailingOnly = T)
setwd(paste0("/data/home/liuyang/RAW_data/1-SingleCell/3-NPC/data/",ags[1]))

#Read in all input expression matrices
# Create and setup Seurat objects for each dataset 
way <- paste0(strsplit(ags[1],"_")[[1]][1],"_",strsplit(ags[1],"_")[[1]][2],"_",strsplit(ags[1],"_")[[1]][3])
TenXdat1 <- Read10X(data.dir = paste0("/data/home/hanbw/project/Single_cell_1808_NPC/",way,"/",ags[1],"/outs/filtered_gene_bc_matrices/refdata-cellranger-GRCh38_and_EBV_5p"))
TenXdat <- CreateSeuratObject(raw.data = TenXdat1, min.cells = 3, min.genes = 500, project = as.character(ags[1]))

mito.genes <- grep(pattern = "^MT-", x = rownames(x = TenXdat@data), value = TRUE)
percent.mito <- Matrix::colSums(TenXdat@raw.data[mito.genes, ])/Matrix::colSums(TenXdat@raw.data)

TenXdat <- AddMetaData(object = TenXdat, metadata = percent.mito, col.name = "percent.mito")
TenXdat <- FilterCells(object = TenXdat, subset.names = c("nGene", "percent.mito","nUMI"), 
    low.thresholds = c(200, -Inf,1000), high.thresholds = c(25000, 0.25,500000))

TenXdat <- NormalizeData(object = TenXdat, normalization.method = "LogNormalize", 
    scale.factor = 10000)	
TenXdat <- ScaleData(TenXdat, display.progress=F)
TenXdat <- FindVariableGenes(TenXdat, do.plot=F)
TenXdat <- RunPCA(object = TenXdat, features = VariableFeatures(object = TenXdat), verbose = FALSE)
# TenXdat <- ProjectDim(object = TenXdat)
TenXdat <- FindClusters(object = TenXdat, resolution = 0.5)
TenXdat <- RunTSNE(object = TenXdat,reduction.use = "pca", dims = 1:20)


#######add meta.data###########
TenXdat@meta.data$Patient <- paste0(strsplit((TenXdat@project.name),"_")[[1]][1],"_",strsplit((TenXdat@project.name),"_")[[1]][2],"_",strsplit((TenXdat@project.name),"_")[[1]][3])
TenXdat@meta.data$Barcode <- paste0(TenXdat@meta.data$Patient,"_",colnames(TenXdat@data))
TenXdat@meta.data$Source <- strsplit((TenXdat@project.name),"_")[[1]][4]

##########Find doublet######
sweep.res.list_TenXdat <- paramSweep(TenXdat, PCs = 1:20)
sweep.stats_TenXdat<- summarizeSweep(sweep.res.list_TenXdat, GT = FALSE)
bcmvn_TenXdat <- find.pK(sweep.stats_TenXdat)


annotations <- TenXdat@meta.data$res.0.5
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*length((TenXdat@meta.data$Barcode)))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

TenXdat <- doubletFinder(TenXdat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
colnames(TenXdat@meta.data[10]) <- "DF"
save(TenXdat, file = paste0(ags[1],"_seurat_2_vision",".RData"))

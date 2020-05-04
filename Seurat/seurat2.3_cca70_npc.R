
#######################################Seurat 2.3.4 subsetdata#######################
library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(DoubletFinder)
rm(list = ls())

#############merge samples###############
list_name = c("NPC_SC_1802_PBMC_cDNA",
              "NPC_SC_1802_Tumor_cDNA",
              "NPC_SC_1805_PBMC_cDNA",
              "NPC_SC_1805_Tumor_cDNA",
              "NPC_SC_1806_PBMC_cDNA",
              "NPC_SC_1806_Tumor_cDNA",
              "NPC_SC_1807_PBMC_cDNA",
              "NPC_SC_1807_Tumor_cDNA",
              "NPC_SC_1808_PBMC_cDNA",
              "NPC_SC_1808_Tumor_cDNA",
              "NPC_SC_1810_PBMC_cDNA",
              "NPC_SC_1810_Tumor_cDNA",
              "NPC_SC_1811_PBMC_cDNA",
              "NPC_SC_1811_Tumor_cDNA",
              "NPC_SC_1813_PBMC_cDNA",
              "NPC_SC_1813_Tumor_cDNA",
              "NPC_SC_1815_PBMC_cDNA",
              "NPC_SC_1815_Tumor_cDNA",
              "NPC_SC_1816_PBMC_cDNA",
              "NPC_SC_1816_Tumor_cDNA"
)
slist = list()
for (i in list_name){
  i <- load(paste0("/data/home/liuyang/liuyang/RAW_data/1-SingleCell/3-NPC/data/",i,"/",i,"_seurat_2_vision_minimal_1%_cells.RData"))
  i <- TenXdat
  slist <- append(slist,i)
} 

################rm doublet###############
for (i in 1:20){
  slist[[i]] <- SubsetData(object = slist[[i]], cells.use = rownames(slist[[i]]@meta.data[slist[[i]]@meta.data$DF == "Singlet",]))
}
###############rename cell names##########
sa <- slist[[1]]
sb <- slist[[2]]
sc <- slist[[3]]
sd <- slist[[4]]
se <- slist[[5]]
sf <- slist[[6]]
sg <- slist[[7]]
sh <- slist[[8]]
si <- slist[[9]]
sj <- slist[[10]]
sk <- slist[[11]]
sl <- slist[[12]]
sm <- slist[[13]]
sn <- slist[[14]]
so <- slist[[15]]
sp <- slist[[16]]
sq <- slist[[17]]
sr <- slist[[18]]
ss <- slist[[19]]
st <- slist[[20]]

sa <- RenameCells(sa, add.cell.id = sa@meta.data$orig.ident )
sb <- RenameCells(sb, add.cell.id = sb@meta.data$orig.ident )
sc <- RenameCells(sc, add.cell.id = sc@meta.data$orig.ident)
sd <- RenameCells(sd, add.cell.id = sd@meta.data$orig.ident)
se <- RenameCells(se, add.cell.id = se@meta.data$orig.ident)
sf <- RenameCells(sf, add.cell.id = sf@meta.data$orig.ident)
sg <- RenameCells(sg, add.cell.id = sg@meta.data$orig.ident)
sh <- RenameCells(sh, add.cell.id = sh@meta.data$orig.ident)
si <- RenameCells(si, add.cell.id = si@meta.data$orig.ident)
sj <- RenameCells(sj, add.cell.id = sj@meta.data$orig.ident)
sk <- RenameCells(sk, add.cell.id = sk@meta.data$orig.ident)
sl <- RenameCells(sl, add.cell.id = sl@meta.data$orig.ident)
sm <- RenameCells(sm, add.cell.id = sm@meta.data$orig.ident)
sn <- RenameCells(sn, add.cell.id = sn@meta.data$orig.ident)
so <- RenameCells(so, add.cell.id = so@meta.data$orig.ident)
sp <- RenameCells(sp, add.cell.id = sp@meta.data$orig.ident)
sq <- RenameCells(sq, add.cell.id = sq@meta.data$orig.ident)
sr <- RenameCells(sr, add.cell.id = sr@meta.data$orig.ident)
ss <- RenameCells(ss, add.cell.id = ss@meta.data$orig.ident)
st <- RenameCells(st, add.cell.id = st@meta.data$orig.ident)

#####calcullate variable genes of all samples#####

com_3000 <- Reduce(intersect, list(v1=head(rownames(sa@hvg.info), 3000),
v2=head(rownames(sb@hvg.info), 3000), v3=head(rownames(sc@hvg.info), 3000),
v4=head(rownames(sd@hvg.info), 3000), v5=head(rownames(se@hvg.info), 3000),
v6=head(rownames(sf@hvg.info), 3000), v7=head(rownames(sg@hvg.info), 3000),
v8=head(rownames(sh@hvg.info), 3000), v9=head(rownames(si@hvg.info), 3000),
v10=head(rownames(sj@hvg.info), 3000),v11=head(rownames(sk@hvg.info), 3000),
v12=head(rownames(sl@hvg.info), 3000),v13=head(rownames(sm@hvg.info), 3000),
v14=head(rownames(sn@hvg.info), 3000),v15=head(rownames(so@hvg.info), 3000),
v16=head(rownames(sp@hvg.info), 3000),v17=head(rownames(sq@hvg.info), 3000),
v18=head(rownames(sr@hvg.info), 3000),v19=head(rownames(ss@hvg.info), 3000),
v20=head(rownames(st@hvg.info), 3000)));

##################Run CCA###########
alist <- list(sa,sb,sc,sd,se,sf,sg,sh,si,sj,sk,sl,sm,sn,so,sp,sq,sr,ss,st)
npc <- RunMultiCCA(object.list=alist, genes.use=com_3000, num.cc=70)

#################Plot###############
p1 <- DimPlot(npc, reduction.use="cca", group.by="Patient", pt.size=0.5, do.return=T)
p2 <- VlnPlot(npc, features.plot="CC1", group.by="Patient", do.return=TRUE)
p3 <- MetageneBicorPlot(npc, grouping.var="Patient", dims.eval=1:70); dev.off()

#########Standardized processes#########
npc <- AlignSubspace(npc,reduction.type="cca",grouping.var="Patient",dims.align=1:20)
npc <- RunUMAP(npc,reduction.use="cca.aligned",n_neighbors=5L,dims.use=1:20,min_dist=0.1)
npc <- FindClusters(npc,reduction.type="cca.aligned",resolution=1,dims.use=1:20)

save(npc,file = "npc.RData")

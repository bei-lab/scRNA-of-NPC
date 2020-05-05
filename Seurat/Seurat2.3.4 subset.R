library(Seurat)
library(Matrix)
library(stringr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(ggsci)

###########---------load data--------############
load("subset.RData")

#####----------find variable gene--------------###
subset <- FindVariableGenes(subset,
                        x.low.cutoff=0.15,
                        x.high.cutoff=3,
                        y.cutoff=0.4,
                        y.high.cutoff = 10,
                        mean.function=ExpMean,
                        dispersion.function=LogVMR)
length(subset@var.genes)

#rm genes
mit=grep(pattern="^MT-",x=subset@var.genes,value=T); com=setdiff(subset@var.genes,mit)
rps <- grep(pattern="^RPS", x=subset@var.genes, value=T); com <- setdiff(com, rps)
rpl <- grep(pattern="^RPL", x=subset@var.genes, value=T); com <- setdiff(com, rpl)

subset@var.genes <- com; length(subset@var.genes)

#select dims
png("subset_MetageneBicorPlot.png",width = 720*12, height = 720*10, res = 720);
MetageneBicorPlot(subset, grouping.var="Patient", dims.eval=1:30);
dev.off()

#scale
subset <- AlignSubspace(subset,reduction.type="cca",
                    grouping.var="Patient",
                    genes.use=subset@var.genes,
                    dims.align=1:20,
                    do.par = T,num.cores = 16); 

subset <- RunUMAP(subset,reduction.use="cca.aligned",
              dims.use=1:20,
              n_neighbors=30L,
              min_dist=0.2);

subset <- FindClusters(subset,dims.use=1:20,
                   reduction.type="cca.aligned",
                   resolution=2,
                   force.recalc=T)

###########--------dimplot-----------########
cols <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_lancet()(9))[-8]
p1 = DimPlot(subset,reduction.use="umap",no.legend=F,do.return=T,pt.size=1,cols.use = cols,
             do.label=T,label.size=6); 
p2=DimPlot(subset,do.label=T,reduction.use="umap",
           no.legend=F,pt.size=0.1,do.return=T,
           group.by="Patient")
p3=DimPlot(subset,do.label=T,reduction.use="umap",
           no.legend=F,pt.size=0.1,do.return=T,
           group.by="Source")
png("subset_dimplot.png",width = 720*12, height = 720*10, res = 720);
plot_grid(p1); dev.off()
plot_grid(p1,p2,p3,ncol = 3); dev.off()

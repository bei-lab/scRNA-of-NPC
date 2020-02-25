library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(monocle)
rm(list = ls())

load("monocle.RData")

#######seurat to monocle######
monocle_data <- importCDS(monocle_data, import_all = TRUE)

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
monocle_data <- estimateSizeFactors(monocle_data)
monocle_data <- estimateDispersions(monocle_data)

###########find variable genes#########
diff_test_res <- differentialGeneTest(monocle_data,
                                      fullModelFormulaStr = "~cell_cluster",cores = 16)

save(diff_test_res,file = "monocle_data_pseudotime_all_gene.RData")

###############select variable genes###
ordering_genes <- row.names (subset(diff_test_res, qval < 0.00000000000000000000000000000000000000000000001))
length(ordering_genes)
RPL_genes <- grep(pattern = c("^RPL"), x = ordering_genes, value = TRUE)
RPS_genes <- grep(pattern = c("^RPS"), x = ordering_genes, value = TRUE)
IGH_genes <- grep(pattern = c("^IGH"), x = ordering_genes, value = TRUE)
IGL_genes <- grep(pattern = c("^IGL"), x = ordering_genes, value = TRUE)
IGK_genes <- grep(pattern = c("^IGK"), x = ordering_genes, value = TRUE)
MT_genes <- grep(pattern = c("^MT-"), x = ordering_genes, value = TRUE)
rm_genes <- c(RPL_genes,RPS_genes,MT_genes,IGK_genes,IGH_genes,IGL_genes)
ordering_genes <- ordering_genes[!(ordering_genes %in% rm_genes)]
length(ordering_genes)

monocle_data_select <- setOrderingFilter(monocle_data, ordering_genes)
monocle_data_select <- reduceDimension(monocle_data_select, max_components = 2,num_dim=10, norm_method="log", 
                             reduction_method="DDRTree",residualModelFormulaStr="~Patient")
monocle_data_select <- orderCells(monocle_data_select)

#######set color############
cell_type_color <- unique((pData(monocle_data_select))$cell_cluster)
cols = c("#BC3C29FF" ,"#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF",
         "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2" ,"#8491B4B2", "#91D1C2B2", "#DC0000B2" ,"#7E6148B2",
         "#00468BFF" ,"#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF" ,"#FDAF91FF", "#AD002AFF" ,"#ADB6B6FF" ,"#1B1919FF")
names(cols) <- cell_type_color

#######################plot###################
plot_cell_trajectory(monocle_data_select,color_by = "cell_cluster",cell_size = 2) +
  scale_color_manual(values = cols)+ 
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))

plot_genes_branched_heatmap(monocle_data_select[TFs,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)

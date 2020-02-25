library(pheatmap)

pheatmap(plot_matrix,cluster_row = F, cluster_col = F,scale = "row",cellwidth = 15,cellheight = 15,
         color = (inferno(20)),border = F,gaps_row = c(4,8,13,21),
         gaps_col =c(),treeheight_col = 0,border_color = "white",angle_col = 45)

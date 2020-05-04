# scRNA-of-NPC
Environment & Dependencies

Red Hat Enterprise Linux Server release 6.6	

R version 3.6.1	

Python version 3.7.4	

Java version 1.8.0_60	


#install R package Seurat v2.3.4 	

source("https://z.umn.edu/archived-seurat")


#install R package DoubletFinder v2.0

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')


#install R package Monocle v2.8 	

source("http://bioconductor.org/biocLite.R") 

biocLite("monocle")	


#install R package STARTRAC	

install.packages("devtools")	

devtools::install_github("Japrin/STARTRAC")	


#install software GSEA v3.0	

download from https://www.gsea-msigdb.org/gsea/downloads.jsp	


#Online software Metascape	

https://metascape.org/gp/index.html	


#insatll software ARACNe-AP 

download from https://github.com/califano-lab/ARACNe-AP


#insatll software CellPhoneDB 

dowload from https://github.com/Teichlab/cellphonedb

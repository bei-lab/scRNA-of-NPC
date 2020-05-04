# scRNA-of-NPC
## Environment 

Red Hat Enterprise Linux Server release 6.6	

R version 3.6.1	

Python version 3.7.4	

Java version 1.8.0_60	

## Install software
### Install R package Seurat v2.3.4 	

\tsource("https://z.umn.edu/archived-seurat")



### Install R package DoubletFinder v2.0
\tremotes::install_github('chris-mcginnis-ucsf/DoubletFinder')


### Install R package Monocle v2.8 	
\tsource("http://bioconductor.org/biocLite.R") 
\tbiocLite("monocle")	


### Install R package STARTRAC v0.1
\tinstall.packages("devtools")	
\tdevtools::install_github("Japrin/STARTRAC")	


### Install software GSEA v3.0	

download from https://www.gsea-msigdb.org/gsea/downloads.jsp	


### Online software Metascape	

https://metascape.org/gp/index.html	


### Insatll software ARACNe-AP 

download from https://github.com/califano-lab/ARACNe-AP


### Insatll software CellPhoneDB v2.06

dowload from https://github.com/Teichlab/cellphonedb

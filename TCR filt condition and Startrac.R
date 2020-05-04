#######################
library("Startrac")
library("ggplot2")
library("RColorBrewer")
library("dplyr")

#############################
TCR <- read.csv("filtered_contig_annotations.csv")
load("tc.RData") #tc Seurat data

#############------------------transform the barcode of cells---------------###########
T_cell_meta.data <- tc@meta.data
a <- T_cell_meta.data$Barcode
barc <- strsplit(a,split = "_")
aa <- list()
for (i in 1:length(barc)){
  aa <- append(aa,barc[[i]][4])
}
T_cell_meta.data$Barcode <- unlist(aa)
T_cell_meta.data$BC <- paste0(T_cell_meta.data$Patient,"-",T_cell_meta.data$Source,"-",T_cell_meta.data$Barcode)
number_to_sample <- data.frame(sort(unique(paste0(tc@meta.data$Patient,"-",tc@meta.data$Source))),c(1:20))
colnames(number_to_sample) <- c("sample","number")
number_cor_tissue_T <- str_split(as.character(TCR$barcode), pattern = "-", simplify = T)
number_cor_tissue_T[, 2] <- mapvalues(number_cor_tissue_T[, 2], from = as.character(number_to_sample$number), to = as.character(number_to_sample$sample))
TCR$Sample <- number_cor_tissue_T[, 2]
TCR$Barcode <- number_cor_tissue_T[, 1]
TCR$cell_barcode <- (with(TCR, paste0(Sample, "-",Barcode)))

##############
TCR <- TCR[!(table(TCR$barcode) == 1),]
###############remove contaminated cells with BCR###########
contamination_cells_T <- unique(TCR[(grepl(pattern = "IG", x = with(TCR, paste0(chain, v_gene, d_gene, j_gene, c_gene)))), "cell_barcode"])
TCR <- TCR[!(TCR$cell_barcode %in% contamination_cells_T), ]

muliti_cells_T <- unique(TCR$cell_barcode[TCR$chain == "Multi"])
TCR <- TCR[!(TCR$cell_barcode %in% muliti_cells_T), ]
TCR_none_rm <- TCR[!TCR$raw_clonotype_id == "None", ]

###remove non T cells TCR##########
TCR_uniq <- TCR_none_rm
TCR_uniq <- filter(TCR_uniq, cell_barcode %in% (T_cell_meta.data$BC))

######UMIS >1##########
TCR_uniq <- TCR_uniq[TCR_uniq$umis > 1,]

#####productive == True########
TCR_uniq <- TCR_uniq[TCR_uniq$productive == "True",]
TCR_uniq$cluster <- mapvalues(x = TCR_uniq$cell_barcode,from = tc@meta.data$BC,to = tc@meta.data$id)

#######select paired TCR########
TCR_split <- by(data = TCR_uniq,INDICES = aa$chain,FUN = `[`)
tcr_bc <- intersect(TCR_split$TRA$cell_barcode,aaa$TRB$cell_barcode)
TRA <- TCR_split$TRA
TRB <- TCR_split$TRB

TRA_rm <- TRA$cell_barcode[TRA$cell_barcode %in% tcr_bc]
TRA <- TRA[TRA$cell_barcode %in% tcr_bc,]
TRA_NAME <- names(table(TRA_rm))[unname(table(TRA_rm)) >1]

TRB_rm <- TRB$cell_barcode[TRB$cell_barcode %in% tcr_bc]
TRB <- TRB[TRB$cell_barcode %in% tcr_bc,]
TRB_NAME <- names(table(TRB_rm))[unname(table(TRB_rm)) >1]
bc_rm <- unique(c(TRA_NAME,TRB_NAME))

TRA_data <- TRA[!TRA$cell_barcode %in% bc_rm,]
TRB_data <- TRB[!TRB$cell_barcode %in% bc_rm,]

############Startrac data###############
#Cell_Name	clone.id	clone.status	patient	sampleType	stype	majorCluster	loc
TCR_uniq<- TRA_data
Cell_Name <- as.character(TCR_uniq$cell_barcode)
clone.id <- as.character(TCR_uniq$raw_clonotype_id)

#############clone cells#############
duplicated_clones <- names(table(TCR_uniq$raw_clonotype_id)[unname(table(TCR_uniq$raw_clonotype_id))>1])
clone.status <- as.character(TCR_uniq$raw_clonotype_id)
clone.status[!clone.status %in% duplicated_clones] <- "NoClonal"
clone.status[clone.status %in% duplicated_clones]  <-"Clonal"

#patient
patient <- as.character(mapvalues(as.character(TCR_uniq$cell_barcode),as.character(T_cell_meta.data$BC),
                                  as.character(T_cell_meta.data$Patient)))

##cluster
majorCluster <- as.character(mapvalues(as.character(TCR_uniq$cell_barcode),as.character(T_cell_meta.data$BC),
                                       as.character(T_cell_meta.data$id)))

##loc
loc <- as.character(mapvalues(as.character(TCR_uniq$cell_barcode),as.character(T_cell_meta.data$BC),
                              as.character(T_cell_meta.data$Source)))

TCR_startrac <- data.frame(Cell_Name,clone.id,clone.status,patient,majorCluster,loc,stringsAsFactors = FALSE)
save(TCR_startrac,file = "TCR_startrac.RData")

######Startrac.run##########
tic("Startrac.run")
out <- Startrac.run(TCR_startrac, proj="CRC", cores=NULL,verbose=F)
toc()

###Startrac.plot###########

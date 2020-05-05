library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)

#####-------CellphoneDB results merge and filter------######

pvalues <-read.table(paste0("./out/pvalues.txt"),sep = "\t",row.names = 1, header = T, stringsAsFactors = F)
signification_means <- read.table(paste0("./out/means.txt"),sep = "\t",row.names = 1, header = T, stringsAsFactors = F)

#Merge
pvalues_npc <- pvalues[,c(1,11:length(pvalues))]
pvalues_npc_1 <- melt(pvalues_npc,id.var=colnames(pvalues_npc)[1])
colnames(pvalues_npc_1) <- c("interacting_pair","cell_type","pvalues")

signification_means_npc <- signification_means[,c(1,11:length(signification_means))]
signification_means_npc_1 <- melt(signification_means_npc,id.var=colnames(signification_means_npc)[1])
colnames(signification_means_npc_1) <- c("interacting_pair","cell_type","signification_means")

cp <- (pvalues_npc_1)[,1]
cell_type <- pvalues_npc_1[,2]
pvalues <- pvalues_npc_1[,3]
signification_means <- signification_means_npc_1[,3]
cp_DB <- data.frame(cp,cell_type,pvalues,signification_means)
write.table(cp_DB,file =paste0("npc_project_merge.txt") , sep = "\t", row.names = F, quote = F)

#Filter P < 0.05
cp_DB <- na.omit(cp_DB)
cp_DB <- cp_DB[cp_DB$pvalues < 0.05,]
write.table(cp_DB,file =paste0("npc_project_filt.txt") , sep = "\t", row.names = F, quote = F)


########-----------Plot-------##########

#Read data
npc_project_merge <- read.table(file = "npc_project_merge.txt",header = T,sep = "\t",stringsAsFactors = F)

#Preprocess
dcat_dt <- dcast(data = npc_project_merge,cp ~ cell_type, value.var = "signification_means")
melt_dt <- melt(dcat_dt,value.name = "cp")
colnames(melt_dt)[3] <- "signification_means"
melt_dt$variable <- as.character(melt_dt$variable)
melt_dt$cp_cell_type <- paste0(melt_dt$cp,melt_dt$variable)
a$cp_cell_type <- paste0(a$cp,a$cell_type )
melt_dt$p_value <- 1
melt_dt$p_value <- mapvalues(melt_dt$cp_cell_type , from = a$cp_cell_type, to = a$pvalues)
melt_dt$p_value[melt_dt$p_value %in% melt_dt$cp_cell_type] <- 1
melt_dt$p_value <- as.numeric(melt_dt$p_value)
melt_dt$p_value[melt_dt$p_value == 0] <- 0.01
melt_dt$signification_means[melt_dt$signification_means > 2.5] <- 2.5

#Select data
Select_data <- read.table("Select_data.txt",header=T,sep='\t',stringsAsFactors = F)
cell_type <- unique(Select_data$cell_type)
cp <- unique(Select_data$ligand_receptor)

Plot_data <- melt_dt[melt_dt$cp %in% cp,]
Plot_data <- Plot_data[Plot_data$variable %in% cell_type, ]
Plot_data$cp <- factor(Plot_data$cp)

#Plot
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(20) %>% rev()
ggplot(Plot_data, aes(x = variable, y = cp)) + geom_point(aes(size = -log10(p_value),color = signification_means)) +
  scale_colour_gradientn(colours = cols) +
  coord_flip()+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle = 90,hjust = 1))

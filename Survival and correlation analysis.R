library(pROC)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(ggsignif)
library(ggpubr)
rm(list = ls())

###############---------Survival analysis------------############

#Read data
NPC_RNAseq <- read.table("NPC_TPM.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
clinical_features <- read.table("SRR_survival_statue.txt", row.names = 1, header = T, sep = "\t", stringsAsFactors = F)

#RNA data overlap with clinical information
names(NPC_RNAseq) <- do.call(rbind, strsplit(names(NPC_RNAseq), split = "_"))[,1]
samples_of_RNA_seq <- names(NPC_RNAseq)
overlap_samples <- intersect(row.names(clinical_features), samples_of_RNA_seq)
Rna_seq <- NPC_RNAseq[, overlap_samples]
clinical_features_overlaps <- clinical_features[overlap_samples, ]
clinical_features_overlaps$XIST <- ""
clinical_features_overlaps$XIST <- as.numeric(Rna_seq["XIST", ])> 1
clinical_features_overlaps[clinical_features_overlaps$XIST, "XIST"] <- "Female"
clinical_features_overlaps[clinical_features_overlaps$XIST==FALSE, "XIST"] <- "Male"
strs <- c("Male", "Female")
names(strs) <- c("TRUE", "FALSE")

#Select signature
genes <- "CCR8_FOXP3"
expression_of_NLRC5 <- with(NPC_RNAseq, Rna_seq[row.names(Rna_seq)== genes, ])

#Calculate ROC
roc.module <- roc(clinical_features_overlaps$DFS, as.numeric(expression_of_NLRC5))
coords.modult <- coords(roc.module, "b", input=c("threshold", "specificity",
                                                 "sensitivity"), ret=c("threshold", "specificity", "sensitivity"),
                        as.list=FALSE, drop=TRUE, best.method=c("closest.topleft"))

#Set cut-off
binary_expression <- as.numeric((expression_of_NLRC5 > coords.modult["threshold"])+0)
clinical_features_overlaps$genes_exp <- binary_expression
km.as.gene <- with(clinical_features_overlaps, survfit(Surv(DFS_time, DFS) ~ genes_exp, data = clinical_features_overlaps, conf.type = "log-log"))

#Plot
ggsurvplot(km.as.gene, conf.int=F, pval=TRUE, risk.table=TRUE,
           legend.labs=c("low", "high"), legend.title=genes,
           palette=c("blue", "red"),
           main="Kaplan-Meier Curve for NPC Survival",
           risk.table.height=.20)


##########--------Correlation analysis-------############

#Read signature genes
gene1 <-read.table("siganture1.txt",sep='\t',stringsAsFactors = F,header = F)
gene2 <-read.table("siganture2.txt",sep='\t',stringsAsFactors = F,header = F)
Rna_seq <- Rna_seq[c(1:27210),]

#Scale
new_matrix <- scale(log2(Rna_seq+1))

#Calculate signature values
gene1_signature <- new_matrix[gene1,] %>% na.omit() %>% colMeans()
gene2_signature <- new_matrix[gene2,] %>% na.omit() %>% colMeans()

#Plot
dat_signature <- data.frame(gene1_signature = gene1_signature, gene2_signature = gene2_signature)
ggplot(data = dat_signature, mapping = aes(x = gene1_signature,
                                       y = gene2_signature)) + 
  geom_point( size = 2) +  geom_smooth(method = lm)+
  stat_cor(method = "pearson" )+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))


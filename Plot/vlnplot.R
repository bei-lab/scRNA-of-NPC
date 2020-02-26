#vlnplot###
library(ggplot2)
library(ggsci)

load("vlnplot_data.RData")

ggplot(vlnplot_data, aes(x=id, y=proliferating,fill=id)) + 
  geom_violin(trim=F,color="white",size = 1 ,scale = "width") 
  geom_boxplot(width=0.1,position=position_dodge(0.9),fill="lightgrey")+ 
  scale_fill_manual(values = cols)+ 
  theme_bw()+ 
  theme(axis.text.x=element_text(angle=30,hjust = 1,colour="black",size=16),
        axis.text.y=element_text(size=12,face="plain"), 
        axis.title.y=element_text(size = 20,face="plain"), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), 
        legend.text=element_text(face="italic",  colour="black",  size=16),
        legend.title=element_text(face="italic", colour="black",  size=16),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
        ylab("Proliferating")+xlab("")
library(ggplot2)

load("barplot_data.RData")

ggplot(data=barplot_data, aes(x=majorCluster, fill =clone.status)) +
  geom_bar(stat="count") +
  theme(axis.text.x = element_text(angle =30, hjust = 1, vjust = 0.95,size= 7))+
  labs(x="cluster",y = "number of clonetypes")+
  scale_fill_manual(values=c( "#8491B4FF" ,"#91D1C2FF","#F39B7FFF"))+
  coord_flip()+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

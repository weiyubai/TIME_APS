diff1 <- read.csv("./output/diff1.csv",header = T,row.names = 1)
diff2 <- read.csv("./output/diff2.csv",header = T,row.names = 1)
TIME_DEGs <- rownames(diff1)
Chose_pathway_DEGs <- rownames(diff2)
overlapping_Degs <- intersect(TIME_DEGs,Chose_pathway_DEGs)
library(ggvenn)
x <- list(TIME_DEGs=TIME_DEGs,Chose_pathway_DEGs=Chose_pathway_DEGs)

ggvenn(x,show_percentage = F,text_size = 4,set_name_size = 3,stroke_color = "white",fill_color = c("#DF8F4499","#00A1D599"),
       set_name_color = c("#DF8F4499","#00A1D599"))
ggsave(file= "./output/vn.pdf",height=4, width=4 )
write.csv(overlapping_Degs,"output/Hub_genes.csv")
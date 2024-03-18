
library(clusterProfiler)

library(GOplot)

library(stringr)
gene<-bitr(rownames(diff1),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db') 
KEGG<-enrichKEGG(
  gene$ENTREZID,
  organism = "hsa",minGSSize =5,
  maxGSSize =500,
  keyType = "kegg",
  use_internal_data = FALSE
)
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')+scale_color_gradient(low="green",high = "red")
write.csv(KEGG,"./output/kegg.csv")
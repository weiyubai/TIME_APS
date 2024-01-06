
library(pheatmap)
library(Seurat)
library(ggplot2)


group <- ifelse(proj@meta.data$Immune >= median(proj@meta.data$Immune),'Immune_high','Immune_low')
table(group)
proj@meta.data$group <- group
diff_dat <- FindMarkers(proj,ident.1="Immune_high",ident.2="Immune_low",
                        group.by='group')
diff1 <- diff_dat[diff_dat$p_val_adj<0.05 & abs(diff_dat$avg_log2FC)>=1,]
write.csv(diff1,"./output/diff1.csv")


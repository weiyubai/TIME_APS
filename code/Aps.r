group1 <- ifelse(SC_data@meta.data$hsa04979>=median(SC_data@meta.data$hsa04979),'hsa04979_high','hsa04979_low')
table(group1)
Sc_data@meta.data$group1 <- group1
diff_dat1 <- FindMarkers(Sc_data,ident.1="hsa04979_high",ident.2="hsa04979_low",
                        group.by='group1')
diff2 <- diff_dat1[diff_dat1$p_val_adj<0.05 & abs(diff_dat1$avg_log2FC)>=1,]
write.csv(diff2,"./output/diff2.csv")
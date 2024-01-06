
folder_path <- "./ESTIMAT"
if (!file.exists(folder_path)) {
  dir.create("./ESTIMAT")
}
myfilterCommonGenes <- function(input.f, output.f, id = c("GeneSymbol", "EntrezID"))
{
  
  id <- match.arg(id)
  input.df <- input.f
  merged.df <- merge(common_genes, input.df, by.x = id, by.y = "row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).",
                nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
  outputGCT(merged.df, output.f)
}

environment(myfilterCommonGenes) <-  environment(estimate::filterCommonGenes)

library("Seurat")
library("estimate")
DefaultAssay(proj) <- "RNA"
myfilterCommonGenes(input.f= as.matrix(proj@assays$RNA@data), output.f=paste0("matrixgenes.gct"), id="GeneSymbol")

estimateScore(paste0("matrixgenes.gct"),paste0("estimate_score.gct"), platform= "affymetrix")

estimate_score <- read.table(paste0("estimate_score.gct"),header = F,skip=2,stringsAsFactors = FALSE)

estimate_score <- t(estimate_score)
rownames(estimate_score) <- estimate_score[,1]
colnames(estimate_score) <- estimate_score[2,]
estimate_score <- estimate_score[-c(1:2),]
estimate_score[1:4,1:4]
estimate_score <- as.data.frame(estimate_score)

proj@meta.data$Stromal <- as.numeric(estimate_score$StromalScore)
proj@meta.data$Immune <- as.numeric(estimate_score$ImmuneScore)
proj@meta.data$TumorPurity <- as.numeric(estimate_score$TumorPurity)
proj@meta.data$ESTIMATE <- as.numeric(estimate_score$ESTIMATEScore)

p1<- RidgePlot(proj,features = "Stromal",group.by="Type")
p2<- RidgePlot(proj,features = "Immune",group.by="Type")
p3<- RidgePlot(proj,features = "TumorPurity",group.by="Type")
p4<- RidgePlot(proj,features = "ESTIMATE",group.by="Type")
# p4<- Seurat::CombinePlots(c(p1  ,p2,p3,p4))
library(gridExtra)
pdf(file= "./output/RidgePlot.pdf",height=4, width=10 )
grid.arrange(p1,p2,p3,p4,ncol = 2 )
dev.off()
library(ggplot2)
mydata1<- FetchData(proj,vars = c("tSNE_1","tSNE_2","Stromal","Immune","TumorPurity","ESTIMATE"))
a1 <- ggplot(mydata1,aes(x = tSNE_1,y =tSNE_2,colour = TumorPurity))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
colours = c('#0425D4',"#0085FF","#00B9D1","#66D404","#EF6909","red","#74120C"))+theme_bw()
a2 <- ggplot(mydata1,aes(x = tSNE_1,y =tSNE_2,colour = Stromal))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#0425D4',"#0085FF","#00B9D1","#66D404","#EF6909","red","#74120C"))+theme_bw()
a3 <- ggplot(mydata1,aes(x = tSNE_1,y =tSNE_2,colour = Immune))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#0425D4',"#0085FF","#00B9D1","#66D404","#EF6909","red","#74120C"))+theme_bw()
a4 <- ggplot(mydata1,aes(x = tSNE_1,y =tSNE_2,colour = ESTIMATE))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#0425D4',"#0085FF","#00B9D1","#66D404","#EF6909","red","#74120C"))+theme_bw()
pdf(file= "./output/HeatMap.pdf",height=4, width=7 )
grid.arrange(a1,a2,a3,a4,ncol = 2 )
dev.off()
file.remove(c("./estimate_score.gct","matrixgenes.gct"))
write.csv(estimate_score,"./output/estimate_score.csv")
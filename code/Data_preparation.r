rm(list=ls())
set.seed(123456)
library(Seurat)
library(ggplot2)
library(data.table)
samples=list.files('./set')

sceList = lapply(samples, function(pro){ 
  folder=file.path('./set', pro) #file.path(a,b) means the path: a/b
  sce=CreateSeuratObject(counts = Read10X(folder),
                         project = pro)
  return(sce)
})

sceList <- lapply(sceList, FUN = function(x) {
  x <- NormalizeData(x)
  
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})


sceall = merge(sceList[[1]], y = sceList[[2]], 
              add.cell.ids = samples, merge.data = TRUE)
#
#library(harmony)
#sceall <- NormalizeData(sceall, normalization.method = "RC", scale.factor = 1e4) 
#sceall <- FindVariableFeatures(sceall)
#sceall <- ScaleData(sceall)
#sceall <- RunPCA(sceall, features = VariableFeatures(object = sceall))
#sceall <- RunHarmony(sceall, group.by.vars = "orig.ident") 

sample1 <- read.table(gzfile("./set/set1/samples.txt.gz"),sep = "\t",header = T)
rownames(sample1) <-sample1$Cell.Barcode
sample2 <- read.table(gzfile("./set/set2/samples.txt.gz"),sep = "\t",header = T)
rownames(sample2) <-sample2$Cell.Barcode
sample <- rbind(sample1,sample2)
sceall@meta.data$sample <- sample$Sample
sceall@meta.data$Type <- sample$Type
Idents(sceall) <- sceall@meta.data$sample
LIHC <- subset(sceall,idents =c("S16_P10_LCP18","S02_P01_LCP21","S10_P05_LCP23","S07_P02_LCP28",
                                "S12_P07_LCP30","S21_P13_LCP37","S15_P09_LCP38","S351_P10_LCP34",
                                "S364_P21_LCP65"),invert = FALSE)
LIHC[["percent.mt"]] <- PercentageFeatureSet(LIHC, pattern = "^MT-")
LIHC <- subset(LIHC, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
LIHC <- NormalizeData(LIHC, normalization.method = "RC", scale.factor = 10000)
LIHC <- FindVariableFeatures(LIHC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(LIHC)
LIHC <- ScaleData(LIHC, features = all.genes)
LIHC <- ScaleData(LIHC, vars.to.regress = "percent.mt")
LIHC <- RunPCA(LIHC, features = VariableFeatures(object = LIHC))
LIHC<- RunTSNE(LIHC, dims = 1:10)


DimPlot(LIHC, group.by = "Type",reduction = "tsne",pt.size=0.5,label = FALSE)
ggsave("output/ScType.pdf",width = 4,height = 2.5)






---
title: "TIME-ASP"
author: "Weiyu Bai"
date: "2024-5-10"
---
# Set the RStudio working directory to the project folder.
```{r}
setwd("D:/TIME_ASP ")
```
# Prepare Seurat object format single-cell sequencing data of interest. 
# Option 1: User’s own or public database single-cell sequencing data.
```{r}
library(Seurat)
library(ggplot2)
#a.	Put the three files ‘barcodes.tsv’, ‘genes.tsv’, and ‘matrix.mtx’ into the 'folder' folder
#b.	Build Seurat object.

folder <- Read10X(data.dir = './folder')
Sc_data <- CreateSeuratObject(counts = folder)
#c.	Calculate the proportion of mitochondrial genes.
# The variable in this step is defined as "Sc_data" in order to connect with the analysis steps of subsequent sample data.
Sc_data[["percent.mt"]] <- PercentageFeatureSet(Sc_data, pattern = "^MT-")
#d.	Quality control.
#i.	Filter cells with UMI (Unique Molecular Identifier) numbers greater than 2500. 
#ii.	Filter cells with UMI numbers less than 200.
#iii.	Filter cells with a percentage of mitochondrial genes greater than 5% of the total number of genes
Sc_data <- subset(Sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#e.	Data standardization.
Sc_data <- NormalizeData(Sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
#The default data standardization method used is LogNormalize, where the total expression level of each cell is standardized to 10000, and then log is taken as the logarithm; The results are stored in Sc_data["RNA"] @ data.
#f.	Perform cell type annotation and linear dimensionality reduction.
Sc_data <- FindVariableFeatures(Sc_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Sc_data)
Sc_data <- ScaleData(Sc_data, features = all.genes)
Sc_data <- ScaleData(Sc_data, vars.to.regress = "percent.mt")
#Using SingleR for cell annotation
library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
lihcSingleR <- GetAssayData(Sc_data, slot = "data")
lihc.hesc <- SingleR(test=lihcSingleR, ref = hpca.se, labels =hpca.se$label.main)
table(lihc.hesc$labels)
Sc_data @meta.data$labels <- lihc.hesc$labels
Sc_data <- RunPCA(Sc_data, features = VariableFeatures(object = Sc_data))
Sc_data<- RunTSNE(Sc_data, dims = 1:10)
DimPlot(Sc_data, group.by ="labels",reduction = "tsne",label = TRUE)
ggsave("ScType.pdf",width = 4,height = 2.5)
save(Sc_data,file="./Sc_data.Rdata")
#Save single-cell sequencing expression matrix.
```
# Option 2: Using sample data.
```{r}
rm(list=ls())
#Clear variables
set.seed(123456)
# Set random seeds
library(Seurat)
library(ggplot2)
library(data.table)
samples=list.files('./set')
sceList = lapply(samples, function(pro){ 
folder=file.path('./set', pro) 
sce=CreateSeuratObject(counts = Read10X(folder),
project = pro)
return(sce)
})
#file.path(a,b) means the path: a/b
sceList <- lapply(sceList, FUN = function(x) {
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
#Build a list based on two datasets
sceall = merge(sceList[[1]], y = sceList[[2]], 
add.cell.ids = samples, merge.data = TRUE) 
# Merge data 
sample1 <- read.table(gzfile("./set/set1/samples.txt.gz"),sep = "\t",header = T)
#Import sample information from dataset1
rownames(sample1) <-sample1$Cell.Barcode
sample2 <- read.table(gzfile("./set/set2/samples.txt.gz"),sep = "\t",header = T)
#Import sample information from dataset2.
rownames(sample2) <-sample2$Cell.Barcode
sample <- rbind(sample1,sample2)
sceall@meta.data$sample <- sample$Sample
sceall@meta.data$Type <- sample$Type
#Import cell annotation information
Idents(sceall) <- sceall@meta.data$sample
Sc_data <-subset(sceall,idents = c("S16_P10_LCP18","S02_P01_LCP21","S10_P05_LCP23","S07_P02_LCP28","S12_P07_LCP30","S21_P13_LCP37","S15_P09_LCP38","S351_P10_LCP34","S364_P21_LCP65"),invert = FALSE)
#Filter data for specified samples (hepatocellular carcinoma)

Sc_data[["percent.mt"]] <- PercentageFeatureSet(Sc_data, pattern = "^MT-")
Sc_data <- subset(Sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
Sc_data <- NormalizeData(Sc_data, normalization.method = "RC", scale.factor = 10000)
Sc_data <- FindVariableFeatures(Sc_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Sc_data)
Sc_data <- ScaleData(Sc_data, features = all.genes)
Sc_data <- ScaleData(Sc_data, vars.to.regress = "percent.mt")
Sc_data <- RunPCA(Sc_data, features = VariableFeatures(object = Sc_data))
Sc_data<- RunTSNE(Sc_data, dims = 1:10)
DimPlot(Sc_data, group.by = "Type",reduction = "tsne",pt.size=0.5,label = FALSE)
ggsave("output/ScType.pdf",width = 4,height = 2.5)
#Note: A file named ‘ScType.pdf’ will be generated in the ‘output’ folder. If you need to eliminate batch effects from your data, you can utilize the Integration method from the Seurat package or the harmony30 package. This step uses data from liver cancer patients in the GSE125449 dataset in the GEO database, with patient numbers listed as follows: "S16-P10-LCP18", "S02-P01-LCP21", "S10-P05-LCP23", "S07-P02-LCP28", "S12-007-LCP30", "S21P13-LCP37", "S15-P09-LCP38", "S351-P10-LCP34", "S364-P21-LCP65"
```
# Calculate immune score and tumor purity
```{r}
#Here, we describe the steps for calculating immune score, stromal score, tumor purity score, and estimate score based on the ESTIMATE R package.
#Compute the immune score, tumor purity score, stromal score, and ESTIMATE score. 
proj <- Sc_data
source("./code/AddModuleScore.r")
```
# Grouping all cells based on immune scores
```{r}
#Here, we describe the steps of grouping all cancer cells based on immune scores. Grouping is performed based on the calculated immune scores from the previous step (or based on tumor purity), and the FindMarkers() function is used for differential analysis. 
source("./code/TIME.r")
#Note: If you want to use tumor purity as the grouping condition, please open the ‘Estimate.r’ file and modify the grouping code as follows:
#group <- ifelse(proj@meta.data$TumorPurity>=median(proj@meta.data$ TumorPurity),' TumorPurity_high',' TumorPurity_low')

```
# Signal pathway enrichment
```{r}
#Here, we describe the steps of pathway enrichment using differentially expressed genes obtained through immune scoring grouping.
source(Kegg.r)
#If the gene species is not human, please modify the 'organism' parameter in the Kegg.r code according to the abbreviation of the queried species. Please refer to 'https://www.genome.jp/kegg/catalog/org_list.html' for the abbreviations of the specific species. For detailed explanations of other additional parameters, please refer to the usage guide of the 'clusterProfiler' R package.
```
# Evaluation of signaling pathway enrichment level
```{r}
#Here, we describe the steps of quantitatively analyzing signal pathways using the 'AddModuleScore' function.
#Computational analysis of pathway expression across different cell types.
#a.	Specify the identifier for the signal pathway that needs to be computed. Using ‘Cholesterol metabolism’ as an example, input the identifier for Cholesterol metabolism.
ID = "hsa04979"
#Note: KEGG ID can be queried from the following two websites: https://www.kegg.jp/ Or http://lmmd.ecust.edu.cn/netinfer/get_results.php?id=6248
#b.	Running the calculation program will generate a file in Seurat object format named Inscore, and the calculated scores will be named after the first word of the signal pathway, contained within the Inscore file.
source("./code/Chose_pathway.r")

#Display multiple signal pathways in the form of heatmaps, and repeat steps 6a-6d. Here, we will use signal pathways such as hsa05200, hsa04510, hsa05208, hsa04932, hsa04151, hsa05205, hsa04610, hsa00190, hsa04512, hsa04979, hsa04066, and hsa04668 as examples to draw a heat map.
source("./code/All_pathways.r")
#Note: In this step, the analysis process from a to k remains exactly the same, except for the different IDs of the signaling pathways. If you need to analyze a specific signaling pathway, simply modify the ID of the pathway. For example, if you need to analyze the ‘Adherens junction signaling pathway (ID hsa04520):
#gsInfo = keggGet(‘hsa04520’)[[1]]; 
#geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
#geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
#names(geneSet) = gsInfo$NAME
#hsa04520  = list(geneSet$`Adherens junction signaling - Homo sapiens (human)`) 
#Inscore <- AddModuleScore(Sc_data, features = hsa04520, ctrl = 10, name = "hsa0452")
#Sc_data@meta.data$hsa04520  <- as.numeric(Inscore$hsa045201)# It should be noted that the algorithm for this line of code will default to adding a sequence number 1 after the ID number. 
```
#Differential analysis of abnormal signaling pathways.
```{r}
#Here, we describe the steps of screening differentially expressed genes based on abnormal signaling pathway grouping.
#Calculate differential genes in single-cell sequencing data by grouping a certain signaling pathway. For example, using the ‘hsa04979’ signaling pathway, all cells are divided into two groups: high and low.
source("./code/Aps.r")
```
# Hub gene screening
```{r}
#Here, we describe the steps for screening hub genes based on differential analysis of immune and abnormal signaling pathways.
#Identify and overlap the differential genes associated with the tumor immune microenvironment and the abnormal signaling pathways for screening (Table S4), and visually represent this information using a Venn diagram.
source("./code/hub_genes.r")
```






#kegg

#Cholesterol metabolism - Homo sapiens
library("KEGGREST")
gsInfo = keggGet('hsa04979')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04979  <- list(geneSet$`Cholesterol metabolism - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04979,
                          ctrl = 10,
                          name = "hsa04979")
SC_data@meta.data$hsa04979  <- as.numeric(Inscore$hsa049791)



#Proteoglycans in cancer
gsInfo = keggGet('hsa05205')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa05205  <- list(geneSet$`Proteoglycans in cancer - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa05205,
                          ctrl = 10,
                          name = "hsa05205")
SC_data@meta.data$hsa05205  <- as.numeric(Inscore$hsa052051)




#Oxidative phosphorylation
gsInfo = keggGet('hsa00190')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa00190  <- list(geneSet$`Oxidative phosphorylation - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa00190,
                          ctrl = 10,
                          name = "hsa00190")
SC_data@meta.data$hsa00190  <- as.numeric(Inscore$hsa001901)



#TNF

gsInfo = keggGet('hsa04668')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04668  <- list(geneSet$`TNF signaling pathway - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04668,
                          ctrl = 10,
                          name = "hsa04668")
SC_data@meta.data$hsa04668  <- as.numeric(Inscore$hsa046681)





#oxygen
gsInfo = keggGet('hsa05208')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa05208  <- list(geneSet$`Chemical carcinogenesis - reactive oxygen species - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa05208,
                          ctrl = 10,
                          name = "hsa05208")
SC_data@meta.data$hsa05208  <- as.numeric(Inscore$hsa052081)



#Pathways in cancer
gsInfo = keggGet('hsa05200')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa05200  <- list(geneSet$`Pathways in cancer - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa05200,
                          ctrl = 10,
                          name = "hsa05200")
SC_data@meta.data$hsa05200  <- as.numeric(Inscore$hsa052001)



#Focal adhesion
gsInfo = keggGet('hsa04510')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04510  <- list(geneSet$`Focal adhesion - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04510,
                          ctrl = 10,
                          name = "hsa04510")
SC_data@meta.data$hsa04510  <- as.numeric(Inscore$hsa045101)




#Non-alcoholic fatty liver diseas
gsInfo = keggGet('hsa04932')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04932  <- list(geneSet$`Non-alcoholic fatty liver disease - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04932,
                          ctrl = 10,
                          name = "hsa04932")
SC_data@meta.data$hsa04932  <- as.numeric(Inscore$hsa049321)

#hsa04151:PI3K-Akt signaling pathway
gsInfo = keggGet('hsa04151')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04151  <- list(geneSet$`PI3K-Akt signaling pathway - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04151,
                          ctrl = 10,
                          name = "hsa04151")
SC_data@meta.data$hsa04151  <- as.numeric(Inscore$hsa041511)


#hsa04610:Complement and coagulation cascades
gsInfo = keggGet('hsa04610')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04610  <- list(geneSet$`Complement and coagulation cascades - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04610,
                          ctrl = 10,
                          name = "hsa04610")
SC_data@meta.data$hsa04610  <- as.numeric(Inscore$hsa046101)

#ECM
#hsa04512:ECM-receptor interaction
gsInfo = keggGet('hsa04512')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04512  <- list(geneSet$`ECM-receptor interaction - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04512,
                          ctrl = 10,
                          name = "hsa04512")
SC_data@meta.data$hsa04512 <- as.numeric(Inscore$hsa045121)

#HIF
#hsa04066:HIF-1 signaling pathway
gsInfo = keggGet('hsa04066')[[1]]
names(gsInfo)
geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet
hsa04066  <- list(geneSet$`HIF-1 signaling pathway - Homo sapiens (human)`)

Inscore <- AddModuleScore(SC_data,
                          features = hsa04066,
                          ctrl = 10,
                          name = "hsa04066")
SC_data@meta.data$hsa04066  <- as.numeric(Inscore$hsa040661)


mydata1<- FetchData(SC_data,vars = c("Type","sample","hsa05200", "hsa04510", "hsa05208", "hsa04932", "hsa04151", "hsa05205", "hsa04610", "hsa00190", "hsa04512", "hsa04979", "hsa04066", "hsa04668"))

ann_col <- data.frame(mydata1[,1:2])
ann_col$id <-rownames(ann_col)
ann_col <- ann_col[order(ann_col[,1]),]
ac <- data.frame(ann_col$Type)
rownames(ac) <- rownames(ann_col)
colnames(ac) <- c("cell type")
ac$sample <- ann_col$sample
mydata <- mydata1[row.names(ac),3:14]


mydata <- t(mydata)
library(pheatmap)
pdf(file= "./output/KEGG_HeatMap.pdf",height=5.7, width=5 )
pheatmap(mydata, show_rownames = T, show_colnames = F,cluster_cols = FALSE, cluster_rows = F,fontsize=6,
         annotation_col = ac)
dev.off()
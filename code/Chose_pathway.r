library(KEGGREST)
library(Seurat)
gsInfo = keggGet(ID)[[1]]
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
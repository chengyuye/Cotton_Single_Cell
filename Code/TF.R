library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(data.table)
library(readxl)
library(monocle)
library(ggsci)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
## TFN using SCODE
#first extract gene expression matri
## read in the 142 TF
TF <- read.table("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/filtered/TF_sorted.txt")
colnames(TF) <- c('TF_ID','Gene_ID','XP_ID','LOC_ID',	'Family','Sort_index')

TF_ID <- TF$LOC_ID

### read in SCODE result
scode <- read.table("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/meanA.txt")
### rename the rownames and colnames
colnames(scode) <- TF_ID
rownames(scode) <- TF_ID

### view scode results
scode[1:5, 1:5]

###save the scode result
write.table(scode, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/scode.txt")

scode <- as.matrix(scode)
library(igraph)
# Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph_from_adjacency_matrix(scode)

m <- data.frame()
TF_family <- TF$Family
source <- colnames(scode)
target <- rownames(scode)
for (i in 1:115) {
  for (j in 1:115) {
    n <- data.frame('TF.id' = source[i],
                  'Target.id' = target[j],
                  'Scode_values' = scode[i,j])
    n$Source <- TF_family[i]
    n$Target <- TF_family[j]
    
    m <- rbind(m, n)
  }
  
}

m$Regulation <- ifelse(m$Scode_values<0,"R","A")
head(m)

### CUT-OFF 0.1
scode_new <- m[,c(4,5,3,6)]
### cutoff-0.1
scode_new <- scode_new %>% filter(abs(scode_new$Scode_values) >= 0.1)
###save the scode result
write.csv(scode_new, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/scode_0.1.csv")

### CUT-OFF 0.5
scode_new <- m[,c(4,5,3,6)]
### cutoff-0.5
scode_new <- scode_new %>% filter(abs(scode_new$Scode_values) >= 0.5)
###save the scode result
write.csv(scode_new, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/scode_0.5.csv")



### CUT-OFF 1.5
scode_new <- m[,c(4,5,3,6)]
### cutoff-1.5
scode_new <- scode_new %>% filter(abs(scode_new$Scode_values) >= 1.5)
###save the scode result
write.csv(scode_new, row.names  = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/scode_1.5.csv")


### CUT-OFF 1.0
scode_new <- m[,c(4,5,3,6)]
### cutoff-1.0
scode_new <- scode_new %>% filter(abs(scode_new$Scode_values) >= 1.0)
###save the scode result
write.csv(scode_new, row.names  = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/scode_1.0.csv")



#### save the results
###save the scode result
write.csv(m, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/SCODE/average_output/scode.csv")


###GENIE3----part 1
#candidate regulators
library(GENIE3)
## read in the 137 TF
TF <- read.table("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/crc_tf.txt")
colnames(TF) <- c('TF_ID','Gene_ID','Family','XP_ID','LOC_ID')

TF_ID <- TF$LOC_ID
TF_family <- TF$Family


regulators <- TF_ID

### read in genes that are asscocated with gossipol pathway
### Calculate each gene expression level
## read in LOC files first
LOC_list <- read.csv("C:/BGI/mianhua/CK_harmony/dim50_annotation/pathway_CKN.csv")
LOC_list <- LOC_list$NCBI基因号
#LOC_list <- as.list(LOC_list)
### read in root cap cell
LRC <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds')
candidate_genes <- c(LOC_list,regulators)

LRC <- subset(LRC, features = candidate_genes)

LOC_list <- intersect(LOC_list, rownames(LRC))
### get expression matrix
expr_matrix <- GetAssayData(LRC, slot = 'data')
expr_matrix <- as.matrix(expr_matrix)

# For reproducibility of results
set.seed(123)
weightMat <- GENIE3(expr_matrix, regulators=regulators, targets = LOC_list,
                    nCores=8, verbose=TRUE)

weightMat[1:5,1:5]

dim(weightMat)


###Get all the regulatory links
linkList <- getLinkList(weightMat)
dim(linkList)

## replace the gene name of source to tf-family name
library(stringr)
for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}

### CUT-OFF 0
###save the  result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3.csv")


### CUT-OFF 0.01
###save the  result
linkList <- getLinkList(weightMat, threshold=0.01)
for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}

write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.01.csv")

### CUT-OFF 0.02
###save the  result
linkList <- getLinkList(weightMat, threshold=0.02)
for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.02.csv")

### CUT-OFF 0.03
###save the  result
linkList <- getLinkList(weightMat, threshold=0.03)
for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.03.csv")

### CUT-OFF 0.03
###save the  result
linkList <- getLinkList(weightMat, threshold=0.05)
for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.05.csv")


### replace the gene name to pathway name
### CUT-OFF 0.01
###save the  result
linkList <- getLinkList(weightMat, threshold=0.01)


LOC <- read.csv("C:/BGI/mianhua/CK_harmony/dim50_annotation/pathway_CKN.csv")
LOC_family <- LOC$对应基因名称.棉酚通路.
loc_id <- LOC$NCBI基因号
for (i in 1:99) {
  LOC_family[i] <- strsplit(LOC_family[i], " ")[[1]][2]
  
}
### replace the target gene id
for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}
### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.01_pathway.csv")




### replace the gene name to pathway name
### CUT-OFF 0.02
###save the  result
linkList <- getLinkList(weightMat, threshold=0.02)

for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
### replace the target gene id

for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}
### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.02_pathway.csv")


### replace the gene name to pathway name
### CUT-OFF 0.05
###save the  result
linkList <- getLinkList(weightMat, threshold=0.05)

for (i in 1:137) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
### replace the target gene id

for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}

### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/CRC_DE/genie3_0.05_pathway.csv")




###GENIE3----part 2
## read in the 137 TF
TF <- read.table("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/9TF.txt")
colnames(TF) <- c('TF_ID','Gene_ID','Family','XP_ID','LOC_ID')

TF_ID <- TF$LOC_ID
TF_family <- TF$Family


regulators <- TF_ID

### read in genes that are asscocated with gossipol pathway
### Calculate each gene expression level
## read in LOC files first
LOC_list <- read.csv("C:/BGI/mianhua/CK_harmony/dim50_annotation/pathway_CKN.csv")
LOC_list <- LOC_list$NCBI基因号
#LOC_list <- as.list(LOC_list)
### read in root cap cell
LRC <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds')
candidate_genes <- c(LOC_list,regulators)

LRC <- subset(LRC, features = candidate_genes)

LOC_list <- intersect(LOC_list, rownames(LRC))
regulators <- intersect(regulators, rownames(LRC))
### get expression matrix
expr_matrix <- GetAssayData(LRC, slot = 'data')
expr_matrix <- as.matrix(expr_matrix)

# For reproducibility of results
set.seed(123)
weightMat <- GENIE3(expr_matrix, regulators=regulators, targets = LOC_list,
                    nCores=8, verbose=TRUE)

weightMat[1:5,1:5]

dim(weightMat)


###Get all the regulatory links
linkList <- getLinkList(weightMat)
dim(linkList)

## replace the gene name of source to tf-family name
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}

### CUT-OFF 0
###save the  result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3.csv")


### CUT-OFF 0.01
###save the  result
linkList <- getLinkList(weightMat, threshold=0.01)
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}

write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.01.csv")

### CUT-OFF 0.02
###save the  result
linkList <- getLinkList(weightMat, threshold=0.02)
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}

write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.02.csv")

### CUT-OFF 0.03
###save the  result
linkList <- getLinkList(weightMat, threshold=0.03)
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.03.csv")

### CUT-OFF 0.05
###save the  result
linkList <- getLinkList(weightMat, threshold=0.05)
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.05.csv")

### CUT-OFF 0.1
###save the  result
linkList <- getLinkList(weightMat, threshold=0.1)
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.1.csv")


### CUT-OFF 0.2
###save the  result
linkList <- getLinkList(weightMat, threshold=0.2)
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.2.csv")



### replace the gene name to pathway name
### CUT-OFF 0.01
###save the  result
linkList <- getLinkList(weightMat, threshold=0.01)


LOC <- read.csv("C:/BGI/mianhua/CK_harmony/dim50_annotation/pathway_CKN.csv")
LOC_family <- LOC$对应基因名称.棉酚通路.
loc_id <- LOC$NCBI基因号
for (i in 1:99) {
  LOC_family[i] <- strsplit(LOC_family[i], " ")[[1]][2]
  
}
## replace tf name
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
### replace the target gene id
for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}
### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.01_pathway.csv")




### replace the gene name to pathway name
### CUT-OFF 0.05
###save the  result
linkList <- getLinkList(weightMat, threshold=0.05)

## replace tf name
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
### replace the target gene id
for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}
### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.05_pathway.csv")



### replace the gene name to pathway name
### CUT-OFF 0.1
###save the  result
linkList <- getLinkList(weightMat, threshold=0.1)

## replace tf name
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
### replace the target gene id
for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}
### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.1_pathway.csv")



### replace the gene name to pathway name
### CUT-OFF 0.2
###save the  result
linkList <- getLinkList(weightMat, threshold=0.2)

## replace tf name
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}
### replace the target gene id
for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}
### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.2_pathway.csv")


### replace the gene name to pathway name
### CUT-OFF 0.5
###save the  result
linkList <- getLinkList(weightMat, threshold=0.5)

## replace tf name
for (i in 1:9) {
  linkList$regulatoryGene <- str_replace_all(linkList$regulatoryGene, TF_ID[i], TF_family[i])
}

### replace the target gene id
for (i in 1:99) {
  linkList$targetGene <- str_replace_all(linkList$targetGene, loc_id[i], LOC_family[i])
}

### save the result
write.csv(linkList, row.names = F, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/9TF/genie3_0.5_pathway.csv")

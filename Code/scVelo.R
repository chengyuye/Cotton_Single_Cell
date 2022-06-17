## import library
library(Seurat)
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(Seurat)
library(kableExtra)
library(devtools)
library(stringr)


LRC <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds')
# save metadata table:
###rename barcode
LRC$barcode <- paste(LRC$Condition, colnames(LRC), sep = "-")
LRC$barcode <- str_sub(LRC$barcode,1,nchar(LRC$barcode)-2)


LRC$UMAP_1 <- LRC@reductions$umap@cell.embeddings[,1]
LRC$UMAP_2 <- LRC@reductions$umap@cell.embeddings[,2]
write.csv(LRC@meta.data, file='C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/metadata.csv', 
          quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(LRC, assay='RNA', slot='counts')
writeMM(counts_matrix, file='C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(LRC@reductions$pca@cell.embeddings, file='C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/pca.csv', 
          quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/gene_names.csv',
  quote=F,row.names=F,col.names=F
)



###############################################
######### method 2 ############
table(LRC$celltype)

# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
write.csv(Embeddings(LRC, reduction = "umap"), file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/cell_embeddings.csv")
# 获取每个细胞的barcode
barcode <- paste(LRC$Condition, Cells(LRC), sep = "-")
barcode <- str_sub(barcode,1,nchar(LRC$barcode)-2)
write.csv(barcode, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/cellID_obs.csv", row.names = FALSE)

# 提取每个细胞的cluster信息
write.csv(LRC@meta.data[, 'seurat_clusters', drop = FALSE], file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/cell_clusters.csv")
# 提取每个细胞的celltype信息
write.csv(LRC@meta.data[, 'celltype', drop = FALSE], file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/cell_celltype.csv")
# 获取celltype的颜色信息
hue_pal()(length(levels(LRC$celltype)))
# 获取cluster的颜色信息
hue_pal()(length(levels(LRC$cluster)))




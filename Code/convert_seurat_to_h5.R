## import library
library(Seurat)
library(SeuratDisk)

LRC1 <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds')
# rename the cell 
colnames(mtx$Ara) <- paste("Ara", colnames(mtx$Ara), sep = "-")
LRC1@meta.data$celltype = paste(LRC@meta.data$celltype) # user adaptation required on own dataset
SaveH5Seurat(LRC1, filename = "LRC.h5Seurat")
Convert("LRC.h5Seurat", dest = "h5ad")


### Load CK.harmony
CK.harmony <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds')
CK.harmony@meta.data$celltype = paste(CK.harmony@meta.data$celltype) # user adaptation required on own dataset
SaveH5Seurat(CK.harmony, filename = "CK_annotated.h5Seurat")
Convert("CK_annotated.h5Seurat", dest = "h5ad")

### Load Ara_combined
Ara <- readRDS('C:/BGI/root_integration/data/Ara/combined_Ara_annotated_doubletremoved.rds')
Ara@meta.data$celltype = paste(Ara@meta.data$celltype) # user adaptation required on own dataset
SaveH5Seurat(Ara, filename = "Ara_annotated.h5Seurat")
Convert("Ara_annotated.h5Seurat", dest = "h5ad")

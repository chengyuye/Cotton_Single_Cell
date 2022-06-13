library(CytoTRACE)
library(future)#加载包

# 或者使用“ncores”（默认值 = 1）进行多线程
#或使用“subsamplingsize”（默认值 = 1,000 个单元格）指示子采样大小
#比如：使用 8 个内核和 1,000 个子样本在快速模式下运行以下数据集
## LRC Expression matrix
expr <- as.matrix(LRC@assays$RNA@counts)


table(is.na(LRC$celltype))
####提取表型文件
table(LRC$celltype)
phe <- LRC$celltype
phe = as.character(phe)
names(phe) <- rownames(LRC@meta.data)

phe <- cds$celltype
phe = as.character(phe)
names(phe) <- rownames(LRC@meta.data)


#### extracting umap infor for figures
umap = LRC@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame()

#### extracting DDRTree infor for figures
DDRTree = reducedDimS(cds) %>%  #坐标信息
  as.data.frame()

DDRTree <- t(DDRTree)
colnames(DDRTree) <- c('Component 1', 'Component 2')
DDRTree <- as.data.frame(DDRTree)
## Run cytotrace
results <- CytoTRACE(expr)

# visulize
plotCytoTRACE(results, phenotype = phe,emb = umap)
plotCytoTRACE(results, phenotype = phe,emb = DDRTree)
plotCytoGenes(results, numOfGenes = 10)


## CK.harmony Expression matrix
expr <- as.matrix(CK.harmony@assays$RNA@counts)


table(is.na(CK.harmony$celltype))
####提取表型文件
table(CK.harmony$celltype)
phe <- CK.harmony$celltype
phe = as.character(phe)
names(phe) <- rownames(CK.harmony@meta.data)


#### extracting umap infor for figures---- celltypes
umap = CK.harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame()


## Run cytotrace
results <- CytoTRACE(expr)

# visulize
plotCytoTRACE(results, phenotype = phe, emb = umap)
plotCytoGenes(results, numOfGenes = 10)

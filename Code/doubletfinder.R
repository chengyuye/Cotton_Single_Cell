## library
library(DoubletFinder)

### path variables ###
CK.path <- "C:/BGI/mianhua/CK_harmony_updated_test/"

## read in data 
### load data ###
CKN.filtered <- readRDS("C:/BGI/mianhua/CKN_filtered.rds") # 10402 cells
CKG.filtered <- readRDS('C:/BGI/mianhua/CKG_filtered.rds') # 21842 cells


#### Normalization ####
CKN.filtered <- NormalizeData(CKN.filtered, normalization.method = "LogNormalize", scale.factor = 10000)
CKG.filtered <- NormalizeData(CKG.filtered, normalization.method = "LogNormalize", scale.factor = 10000)


### feature selection: HVGs ### -----method 2
CKN.filtered <- FindVariableFeatures(CKN.filtered, selection.method = "vst", nfeatures = 2000)
CKG.filtered <- FindVariableFeatures(CKG.filtered, selection.method = "vst", nfeatures = 2000)


### scale data ###
CKN.filtered <- ScaleData(CKN.filtered) # uncorrected
CKG.filtered <- ScaleData(CKG.filtered)

### dimensionality reduction: PCA ###
CKN.filtered <- RunPCA(CKN.filtered)
CKG.filtered <- RunPCA(CKG.filtered)


## elbow plot - CKN ##
pca.elbow.plot <- ElbowPlot(CKN.filtered, ndims = 50, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckn.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## elbow plot - CKG ##
pca.elbow.plot <- ElbowPlot(CKG.filtered, ndims = 50, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckg.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## run umap
CKN.filtered <- RunUMAP(CKN.filtered, reduction = "pca", dims = 1:50)
CKG.filtered <- RunUMAP(CKG.filtered, reduction = "pca", dims = 1:50)

## FindNeighbors
CKN.filtered <- FindNeighbors(CKN.filtered, reduction = "pca", dims = 1:50)
CKG.filtered <- FindNeighbors(CKG.filtered, reduction = "pca", dims = 1:50)

## FindClusters
CKN.filtered <- FindClusters(CKN.filtered, resolution = 0.5)
CKG.filtered <- FindClusters(CKG.filtered, resolution = 0.5)

## 寻找最优pK值 CKN
sweep.res.list <- paramSweep_v3(CKN.filtered, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = ncol(CKN.filtered)*8*1e-6                    # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(CKN.filtered$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(CKN.filtered)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
CKN.filtered <- doubletFinder_v3(CKN.filtered, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F)

## 结果展示，分类结果在pbmc@meta.data中
p <- DimPlot(CKN.filtered, reduction = "umap", group.by = "DF.classifications_0.25_0.26_769")
ggsave('CKN_doublet_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/",
       width = 7, height = 5)

## remove doublets ----9633 cells remained
CKN.filtered <- subset(CKN.filtered, subset = DF.classifications_0.25_0.26_769 == 'Singlet')

### save data ###
saveRDS(CKN.filtered, paste0(CK.path, "CKN_doublet_filtered.rds"))

## 寻找最优pK值 CKG
sweep.res.list <- paramSweep_v3(CKG.filtered, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = ncol(CKG.filtered)*8*1e-6                    # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(CKG.filtered$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(CKG.filtered)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
CKG.filtered <- doubletFinder_v3(CKG.filtered, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                                 nExp = nExp_poi.adj, reuse.pANN = F)

## 结果展示，分类结果在pbmc@meta.data中
p <- DimPlot(CKG.filtered, reduction = "umap", group.by = "DF.classifications_0.25_0.12_998")
ggsave('CKG_doublet_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/",
       width = 7, height = 5)

## remove doublets ----10672 cells remained
CKG.filtered <- subset(CKG.filtered, subset = DF.classifications_0.25_0.12_998 == 'Singlet')

### save data ###
saveRDS(CKG.filtered, paste0(CK.path, "CKG_doublet_filtered.rds"))

## re-run finderallmarkers 
## compute cluster marker genes ###
CK.markers <- FindAllMarkers(CK.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
CK.top10.markers <- CK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(CK.markers, file = "C:/BGI/mianhua/CK_harmony/dim50_annotation/CK_Marker_genes_afterdoublet.csv")
write.csv(CK.top10.markers, file = "C:/BGI/mianhua/CK_harmony/dim50_annotation/top10_marker_genes_afterdoublet.csv")

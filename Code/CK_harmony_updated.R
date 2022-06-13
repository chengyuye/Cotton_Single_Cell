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
library(ggsci)

### path variables ###
CK.path <- "C:/BGI/mianhua/CK_harmony_updated_test/"

### load data ###
CKN.filtered <- readRDS("C:/BGI/mianhua/CK_harmony_updated_test/CKN_doublet_filtered.rds") # 21842 cells
CKG.filtered <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/CKG_doublet_filtered.rds') # 21842 cells


### Merge CK together ###
CK <- merge(CKN.filtered, CKG.filtered)

#### Normalization ####
CK <- NormalizeData(CK, normalization.method = "LogNormalize", scale.factor = 10000)


### feature selection: HVGs ### -----method 2
CK <- FindVariableFeatures(CK, selection.method = "vst", nfeatures = 2000)


### scale data ###
CK <- ScaleData(CK) # uncorrected

### dimensionality reduction: PCA ###
CK <- RunPCA(CK, features = VariableFeatures(object = CK),npcs = 100)


## elbow plot - CK ##
pca.elbow.plot <- ElbowPlot(CK, ndims = 100, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.condition.plot <- DimPlot(object = CK, reduction = "pca", pt.size = .1, group.by = "Condition")
png(paste0(CK.path,"PC1_2.condition.png"), width=1000,height=1000,units="px")
print(PC1_2.condition.plot)
dev.off()

# ## jackstraw ---- take too much time----withdraw
# CK <- JackStraw(object = CK, num.replicate = 100, dims=100); # takes around 4 minutes
# CK <- ScoreJackStraw(object = CK, dims = 1:100)
# jackstraw_plot <- JackStrawPlot(object = CK, dims = 1:100 )
# png(paste0(CK.path,"jackstraw_100_PC.png"), width=1000,height=1000,units="px")
# print(jackstraw_plot)
# dev.off()


# save genes making up the PCs  
sink(paste0(CK.path, "PC_genes.txt"))
print(CK[["pca"]], dims = 1:100, nfeatures = 20)
sink()


### integration with harmony ###
# harmonize samples
CK.harmony <- CK %>% RunHarmony("Condition", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) 

# # ## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(CK.harmony, ndims = 100, reduction = "harmony_theta2")
png(paste0(CK.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()
# #
## explore harmony coordinates ##
## visualize PCs ##

harmony.PC1_2.condition.plot <- DimPlot(object = CK.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "Condition")
png(paste0(CK.path,"harmony_theta2.PC1_2.condition.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.condition.plot)
dev.off()
# #
# # harmony.PC1_2.studies.plot <- DimPlot(object = kidney.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
# # png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
# # print(harmony.PC1_2.studies.plot)
# # dev.off()
# #
# # ### save genes making up the PCs ###
sink(paste0(CK.path, "harmony_PC_genes.txt"))
print(CK.harmony[["harmony_theta2"]], dims = 1:100, nfeatures = 20)
sink()
# 
# ### save data ###
# #saveRDS(kidney.harmony, paste0(harmony.samples.path, "all_kidney_harmony.rds"))

### explore different numbers of harmony PCs ###

# visualize marker gene expression
### explore different numbers of harmony PCs ###
# Lateral root cap genes
LTC.genes <- c("LOC107891486","LOC107953062")

# Protophloem genes 
protophloem.genes <- c("LOC107887822")

# Pericycle genes
pericycle <- c('LOC107891006')

# Stele
stele.gene <- c('LOC107944462')

# Root cortex genes
RC.genes <- c("LOC107923064")

# Non???hair root epidermal cell genes
Nh_RE.genes <- c('LOC107951031')

# Columella
columella.genes <- c('LOC107934072')

# Root endodermis genes 
RE.genes <- c("LOC107930672")

# Quiescent center
QC.genes <- c('LOC107907772')

# Trichoblast
trichoblast.genes <- c('LOC107959956')

# Columella root cap cell
CRC.genes <- c('LOC107954183')

# Xylem
xylem.genes <- c('LOC107891487')

# Phloem
phloem.genes <- c("LOC107895190")

markers <- c(LTC.genes, protophloem.genes,pericycle, 
             stele.gene,RC.genes,Nh_RE.genes,columella.genes,
             RE.genes, QC.genes, trichoblast.genes,
             CRC.genes,xylem.genes,phloem.genes)

library(MySeuratWrappers)  


#颜色：和UMAP图一致
col <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
p <- VlnPlot(CK.harmony, features = markers,  
             stack=T, pt.size=0,  
             #cols = col,#颜色  
             direction = "horizontal", #水平作图
             x.lab = '', y.lab = '')+#横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度
p <- p + scale_fill_d3('category20')


ggsave('CK_stacked_markers.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/umap/",
       width = 8, height = 6)
## colors 
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
StackedVlnPlot(CK, markers, pt.size=0, cols=my36colors)

### For CK ###
dims <- c(50,60,70)
for(d in dims){
  
  # create folder
  dir.create(paste0(CK.path, "dim", d, "_annotation"))
  dim.path <- paste0(CK.path, "dim", d, "_annotation/")
  
  # run UMAP
  CK.harmony <- RunUMAP(CK.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(CK.harmony, reduction = "umap", pt.size = 0.1)
  condition.batch.plot <- DimPlot(CK.harmony, reduction = "umap", group.by = "Condition", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=900, units="px")
  print(batch1.plot)
  dev.off()
  
  
  # # no.6 - colored by condition
  # batch6.plot <- DimPlot(CK.harmony, reduction = "umap", group.by = "Condition", pt.size = 0.01)
  # png(paste0(dim.path, "UMAP_dim", d, ".condition.batch1.png"), width=1600, height=1000, units="px")
  # print(batch6.plot)
  # dev.off()
  
  # no.7 - split by condition
  batch7.plot <- DimPlot(CK.harmony, reduction = "umap", split.by = "Condition", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch2.png"), width=1800, height=900, units="px")
  print(batch7.plot)
  dev.off()
  
  # # no.8 - colored by lineage
  # batch8.plot <- DimPlot(CK.harmony, reduction = "umap", group.by = "lineage", pt.size = 0.01)
  # png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch1.png"), width=1600, height=1000, units="px")
  # print(batch8.plot)
  # dev.off()
  
  # # no.9 - colored by annotation
  # batch9.plot <- DimPlot(CK.harmony, reduction = "umap", group.by = "ident", label=TRUE, pt.size = 0.01)
  # png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch2.png"), width=1600, height=1000, units="px")
  # print(batch9.plot)
  # dev.off()
  
  # # no.10 - lineage overview
  # lineage.overview <- condition.batch.plot + batch8.plot
  # png(paste0(dim.path, "UMAP_dim", d, ".lineage.overview.png"), width=2000, height=1000, units="px")
  # print(lineage.overview)
  # dev.off()
  
  ## compare CKN vs. CKG ##
  compare <- DimPlot(CK.harmony, reduction = "umap", split.by = "Condition", group.by = "Condition", 
                     pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
  png(paste0(dim.path, "UMAP_dim", d, "CKN_CKG_1.png"), width=1800, height=900, units="px")
  print(compare)
  dev.off()
  
  # ## compare original annotation ##
  # annotation <- DimPlot(CK.harmony, reduction = "umap", split.by = "Condition", group.by = "celltype",
  #                       label = TRUE, repel = TRUE, pt.size = 0.01) + NoLegend()
  # png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_2.png"), width=1800, height=900, units="px")
  # print(annotation)
  # dev.off()
  
  # ## compare original annotation in healthy vs. fibrotic data ##
  # plot <- compare / annotation
  # png(paste0(dim.path,"UMAP_dim", d, "healthy_fibrotic_3.png"), width=2000, height=2000, units="px")
  # print(plot)
  # dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(CK.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # # overview feature plot
  # lineage.overview <- FeaturePlot(CK.harmony, features = lineage.genes, pt.size = 0.5, ncol = 4) & 
  #   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  # png(paste0(dim.path,"UMAP_dim", d, "overview.markers.png"), width=1600,height=500,units="px")
  # print(lineage.overview )
  # dev.off()
  
}

### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
res <- seq(0.2,0.6,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(CK.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  CK.harmony <- RunUMAP(CK.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # create kNN graph
  CK.harmony <- FindNeighbors(CK.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  for (r in res) {
    
    CK.harmony <- FindClusters(CK.harmony, reduction = "harmony_theta2", resolution = r)
    umap.plot <- DimPlot(CK.harmony, reduction = "umap", repel = T, label = TRUE, pt.size = 0.1)
    
    # create eval plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=800, height=600, units="px")
    print(umap.plot)
    dev.off()
    
    ## clustering by condition
    condition.clustering <- DimPlot(CK.harmony, split.by = "Condition", pt.size = 0.1)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r,".condition.png"), width=1000,height=600,units="px")
    print(condition.clustering)
    dev.off()
    
  }
}

### preliminary clustering ###
CK.harmony <- FindNeighbors(CK.harmony, reduction = "harmony_theta2", dims = 1:50)
CK.harmony <- FindClusters(CK.harmony, reduction = "harmony_theta2", resolution = 0.4)
# run UMAP
CK.harmony <- RunUMAP(CK.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

# run tsne
CK <- RunTSNE(CK, reduction = "harmony_theta2", dims = 1:50, seed.use = 1)
### save data ###
saveRDS(CK.harmony, paste0(CK.path, "CK_harmony.rds"))


## compute cluster marker genes ###
CK.markers <- FindAllMarkers(CK.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
CK.top10.markers <- CK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(CK.markers, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/CK_Marker_genes_res0.4.csv")
write.csv(CK.top10.markers, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/top10_marker_genes_res0.4.csv")

### cell type annotation ###
cluster.annotation <- c("Unknown 1", 'Lateral root cap 1', 'Protophloem',
                        'Pericycle', 'Stele','Root cortex', 
                        'Non-hair root epidermal cell','Columella',  
                        'Root endodermis', 'Quiescent center', 
                        'Trichoblast', 'Columella root cap cell', 'Xylem',
                        'Lateral root cap 2', 'Phloem', 'Unknown 2')


names(cluster.annotation) <- levels(CK.harmony)
CK.harmony <- RenameIdents(CK.harmony, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(CK.harmony),
                        celltype = Idents(CK.harmony))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
CK.harmony <- AddMetaData(CK.harmony, cell.data, col.name = "celltype")

p <- DimPlot(CK.harmony, group.by = 'celltype', label = T, label.size = 3,pt.size = 0.01, repel = T)
p <- p + scale_color_d3('category20')

ggsave('CK_annotated_umap_res0.4.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 8, height = 5)


### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(CK.harmony),
                        celltype = Idents(CK.harmony))

lineage.annotation <- c("Unknown","Root cap cell (RCC)","Protophloem",
                        "Pericycle","Stele","Cortex","Epidermal cell/non root hair (NRH)",
                        "Columella","Endodermis","QC",
                        "Trichoblast","Root cap cell (RCC)","Xylem (X)",
                        "Root cap cell (RCC)","Phloem (P)", "Unknown")

# lineage.annotation <- c("Unknown","Root cap cell (RCC)","Protophloem",
#                         "Pericycle","Stele","Cortex","Epidermal cell/non root hair (NRH)",
#                         "Columella","Endodermis","QC",
#                         "Trichoblast","Root cap cell (RCC)","Stele",
#                         "Root cap cell (RCC)","Stele", "Unknown")

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
CK.harmony <- AddMetaData(CK.harmony, meta.data, col.name = "lineage")

p <- DimPlot(CK.harmony, group.by = 'lineage', label = T, label.size = 3, pt.size = 0.01, repel = T)
p <- p + scale_color_d3('category20')

ggsave('CK_annotated_umap_lineage_doubletremoved_updated.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 8, height = 5)

### save data ### 19279 cells remained
saveRDS(CK.harmony, "C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds")

## remove doublets
CK.harmony <- subset(CK.harmony, subset = celltype != 'Doublet 1')
CK.harmony <- subset(CK.harmony, subset = celltype != 'Doublet 2')

## re-umap
# run UMAP
CK.harmony <- RunUMAP(CK.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

## plot umap 
p <- DimPlot(CK.harmony, group.by = 'celltype', label = T, label.size = 3,pt.size = 0.01, repel = T)
p <- p + scale_color_d3('category20')

ggsave('CK_annotated_res0.4_doubletremoved.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 8, height = 5)
# run tsne
CK.harmony <- RunTSNE(CK.harmony, reduction = "harmony_theta2", dims = 1:75, seed.use = 1)


##colors
# SAVE FIGURES 
library(ggsci)
p <- DimPlot(CK.harmony, reduction = "umap", label = T,   
             cols= allcolour, #设置颜色  
             pt.size = 0.01,#设置点的大小  
             repel = T)+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


p <- p+ scale_color_d3("category20")

ggsave('CK_annotated_umap.pdf', path = "C:/BGI/mianhua/CK_harmony_updated_test/dim75_annotation/",
       width = 8, height = 5)

##
p <- DimPlot(CK.harmony, reduction = "umap", label = T,
             group.by = 'lineage',
             pt.size = 0.01,#设置点的大小  
             repel = T)+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


p <- p+ scale_color_d3("category20")

ggsave('CK_annotated_umap_lineage.pdf', path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 8, height = 5)

p <- DimPlot(CK.harmony, reduction = "umap", label = T,   
             cols= col, #设置颜色  
             pt.size = 0.5,#设置点的大小  
             repel = T,
             split.by = "Condition")+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
ggsave('CK_annotated_umap_split_by_condition.pdf', path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/",
       width = 12, height = 5)

### extracting umap infor for figures---- celltypes
umap = CK.harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = CK.harmony@meta.data$celltype) # 注释后的label信息 ，改为cell_type

### extracting umap infor for figures -----cell lineage
umap = CK.harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_lineage = CK.harmony@meta.data$lineage) # 注释后的label信息 ，改为cell_type

head(umap)
## colors
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
## generate figures
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 0.01 , alpha =1 ) + 
  scale_color_manual(values = allcolour)

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_lineage)) +
  geom_point(size = 0.01 , alpha =1 ) + 
  scale_color_manual(values = allcolour)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))

p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=12), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=2))) #设置legend中 点的大小 


p4 <- p3 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p4 <- p4+scale_color_d3("category20")

p4 <- p4+		scale_color_d3('category20c')

##celltype
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

## cell lineage
cell_lineage_med <- umap %>%
  group_by(cell_lineage) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

library(ggrepel)
##lineage
p4 <- p4 + geom_text_repel(aes(label= cell_lineage), fontface="bold", 
                           data = cell_lineage_med,
                     point.padding=unit(1.6, "lines"))

##cell type
p4 <- p4 + geom_text_repel(aes(label= cell_type), fontface="bold", 
                           data = cell_type_med,
                           point.padding=unit(1.6, "lines"))

##
ggsave('CK_annotated_umap_doublet_removed_refined.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)
##
ggsave('CK_annotated_lineage_doublet_removed_refined.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)

##
ggsave('CK_annotated_lineage_doublet_removed_updated.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)
##
ggsave('CK_annotated_lineage_doublet_removed_labeled.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)

##
ggsave('CK_annotated_umap_doublet_removed_labeled.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)
### save data ###
saveRDS(CK.harmony, "C:/BGI/mianhua/CK_harmony/combined_CK_annotated.rds")


### pathway expression function define
path_expr <- function(object, gene_list, path_name){
  # calculate the module score first
  object <- AddModuleScore(object, features = list(gene_list), name = path_name)
  ## feature plot to see the expression of the whole pathway
  #1. not split by condition
  p <- FeaturePlot(object,
                   features = paste0(path_name, '1'), label = TRUE, repel = TRUE,
                   cols=c("grey","yellow","red","brown"),pt.size = 0.01,label.size = 3
  )

  # save the plot
  ggsave(paste0(path_name,"_expression.pdf"), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/pathway_expression/",
         width = 6, height = 5) 
  #2. split by condition
  p <- FeaturePlot(object,
                   features = paste0(path_name, '1'), 
                   label = TRUE, repel = TRUE, 
                   #min.cutoff = "q10", max.cutoff = "q90",
                   split.by = 'Condition',
                   cols=c("grey","yellow","red","brown"),pt.size = 0.01,label.size = 3)
    #&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  # save the plot
  ggsave(paste0(path_name,"_expression_by_condition.pdf"), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/pathway_expression/",
         width = 12, height = 5)
}


#2-ODD-1 pathway
ODD <- c('LOC107938034','LOC107938040','LOC107938601',
         'LOC107938602','LOC107938615','LOC107940602',
         'LOC107941626','LOC107941627','LOC107941629',
         'LOC107942041','LOC107942044','LOC107944737')

path_expr(CK.harmony,ODD,'ODD')

### MVK Pathway
MVK <- c('LOC107930932')
p <- FeaturePlot(CK.harmony, features = MVK, label = T, repel = T,
                 cols=c("grey","yellow","red","brown"),pt.size = 0.01,label.size = 3)
ggsave("MVK_expression.pdf", p, 
       path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/pathway_expression/",
       width = 6, height = 5) 

p <- FeaturePlot(CK.harmony, features = MVK, label = T, 
                 repel = T,split.by = 'Condition',
                 cols=c("grey","yellow","red","brown"),pt.size = 0.01,label.size = 3)

ggsave("MVK_expression_by_condition.pdf", p, 
       path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/pathway_expression/",
       width = 12, height = 5)

## HMGR pathway   ###　LOC107921012, LOC121203148　not exist
HMGR <- c('LOC107894010','LOC107898896','LOC107900974',
          'LOC107908970','LOC107908971','LOC107908972',
          'LOC107908973','LOC107920223','LOC107921012',
          'LOC107928879','LOC107929696','LOC107939227',
          'LOC107939235','LOC107942380','LOC107949229',
          'LOC121203148') ### 

path_expr(CK.harmony,HMGR,'HMGR')


### HMGS pathway
HMGS <- c('LOC107897938','LOC107909098','LOC107926921',
          'LOC107928371','LOC107929285','LOC107929989',
          'LOC107936427','LOC107936610','LOC107946991',
          'LOC107963679')
path_expr(CK.harmony,HMGS,'HMGS')


## ACAT pathway
ACAT <- c('LOC107912464','LOC107924144',
          'LOC107960146','LOC107963730')
path_expr(CK.harmony,ACAT,'ACAT')


### CYP71BE79 
CYP71BE79 <- c('LOC107897460','LOC107898600')
path_expr(CK.harmony, CYP71BE79, 'CYP71BE79')

### CYP82D113
CYP82D113 <- c('LOC107887028','LOC107903956','LOC107905382',
               'LOC107905383','LOC107935587','LOC107942683',
               'LOC107944158','LOC107944380','LOC107959934')
path_expr(CK.harmony, CYP82D113, 'CYP82D113')

###  DH1 pathway
DH1 <- c('LOC107925529',   ###LOC107928169 not exist
         'LOC107925530',
         'LOC107925589',
         'LOC107925690',
         'LOC107928123',
         'LOC107928135',
         'LOC107928149',
         'LOC107928169'
)
path_expr(CK.harmony, DH1, 'DH1')

### CYP706B1
CYP706B1 <- c('LOC107920158',
              'LOC107963672'
)
path_expr(CK.harmony, CYP706B1, 'CYP706B1')

### CDNC
CDNC <- c('LOC107895206', ### LOC107903509, LOC107918736, LOC107919150 not exist
          'LOC107903499',
          'LOC107903504',
          'LOC107903509',
          'LOC107918736',
          'LOC107918752',
          'LOC107919150',
          'LOC107920688',
          'LOC107920806',
          'LOC107925552',
          'LOC107928183'
)
path_expr(CK.harmony, CDNC, 'CDNC')

### FPS
FPS <- c('LOC107905701',
         'LOC107905737',
         'LOC107922380',
         'LOC107944620'
)
path_expr(CK.harmony, FPS, 'FPS')

### PMD
PMD <- c('LOC107899945',
         'LOC107901154',
         'LOC107901410',
         'LOC107903883',
         'LOC107914631',
         'LOC107915093',
         'LOC107919001',
         'LOC107933302',
         'LOC107958617'
)
path_expr(CK.harmony, PMD, 'PMD')

### MVP
MVP <- c('LOC107946659',
         'LOC107946782'
)
path_expr(CK.harmony, MVP, 'MVP')


### Calculate each gene expression level
## read in LOC files first
LOC_list <- read.csv("C:/BGI/mianhua/CK_harmony/dim50_annotation/pathway_CKN.csv")
LOC_list <- LOC_list$NCBI基因号
celltypes <- as.vector(unique(CK.harmony@meta.data$celltype))


### obtain the gene expression level of each gene in each cluster
gene_exp <- function(seu_object,cluster,gene){
  tmp <- subset(seu_object, subset = celltype == cluster)
  expr <- GetAssayData(object = tmp, slot = "data")
  expr_indi <- as.numeric(Matrix::colSums(FetchData(object = tmp, 
                                                    vars = c(gene))))
  return(expr_indi)
}

### for CKG 
CKG <- subset(CK.harmony, subset = Condition == 'CKG')


###  obtain a list of LOC genes that have zero expression 
zeros <- which(Matrix::rowSums(CKG) == 0)
zeros_list <- names(zeros)
loc_zero_ckg <- intersect(x = LOC_list, y = zeros_list)
### obtain a list of LOC genes that are not present in the Seurat object
whole_list <- rownames(CKG)
absent_list_ckg <- setdiff(x = LOC_list, y = whole_list)
present_list_ckg <- intersect(LOC_list, whole_list)


### run for all genes for all clusters
loc_expr <- data.frame(matrix(NA, ncol=0, nrow = 82))
#colnames(loc_expr) <- celltypes
rownames(loc_expr) <- present_list_ckg
#loc_expr <- as_tibble(loc_expr)
for (i in celltypes) {
  res <- lapply(present_list_ckg, function(x) gene_exp(CKG, i, x))
  names(res) <- present_list_ckg
  res <- as.data.frame(res)
  res <- t(res)
  colnames(res) <- i
  loc_expr <- cbind(loc_expr,res)
  
}

### save the results
write.csv(loc_expr,file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/CKG_path_expression.csv")

## for CKN
CKN <- subset(CK.harmony, subset = Condition == 'CKN')

###  obtain a list of LOC genes that have zero expression 
zeros <- which(Matrix::rowSums(CKN) == 0)
zeros_list <- names(zeros)
loc_zero_ckn <- intersect(x = LOC_list, y = zeros_list)
### obtain a list of LOC genes that are not present in the Seurat object
whole_list <- rownames(CKN)
absent_list_ckn <- setdiff(x = LOC_list, y = whole_list)
present_list_ckn <- intersect(LOC_list, whole_list)


### run for all genes for all clusters
loc_expr <- data.frame(matrix(NA, ncol=0, nrow = 82))
#colnames(loc_expr) <- celltypes
rownames(loc_expr) <- present_list_ckn
#loc_expr <- as_tibble(loc_expr)
for (i in celltypes) {
  res <- lapply(present_list_ckn, function(x) gene_exp(CKN, i, x))
  names(res) <- present_list_ckn
  res <- as.data.frame(res)
  res <- t(res)
  colnames(res) <- i
  loc_expr <- cbind(loc_expr,res)
  
}

### save the results
write.csv(loc_expr,file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/CKN_path_expression.csv")





CK.harmony <- MetaFeature(object = CK.harmony, features = ACAT, meta.name = "ACAT") 


CK.harmony <- MetaFeature(object = CK.harmony, features = ACAT, meta.name = "ACAT") 
data = object@meta.data$oxidation_reduction
threshold = sort(data, decreasing = T)[0.2*length(data)]
object@meta.data$oxidation_reduction = replace(data,data<threshold,0)
FeaturePlot(object = object, features = "oxidation_reduction",max.cutoff = "q80")

markers_percent <- Matrix::colSums(GetAssayData(object = CK.harmony, slot = "counts")[ACAT, ])/Matrix::colSums(GetAssayData(object = CK.harmony, slot = "counts"))
CK.harmony <- AddMetaData(object = CK.harmony, metadata = markers_percent , col.name = "markers_percent")
### save data ###
saveRDS(CK.harmony, "C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated.rds")



### TESTING- deleting the Unknown 1
CK.subset <- subset(CK.harmony, subset = celltype != 'Unknown 1')
CK.subset <- subset(CK.subset, subset = celltype != 'Unknown 2')
ckn <- subset(CK.subset, subset = Condition != 'CKN')#8110
ckg <- subset(CK.subset, subset = Condition != 'CKG')#8288
## re-runing umap
### preliminary clustering ###
CK.subset <- FindNeighbors(CK.subset, reduction = "harmony_theta2", dims = 1:50)
CK.subset <- FindClusters(CK.subset, reduction = "harmony_theta2", resolution = 0.5)

# run UMAP
CK.subset <- RunUMAP(CK.subset, reduction = "harmony_theta2", dims=1:50, seed.use=1)


### save data ###
saveRDS(CK.subset, "C:/BGI/mianhua/CK_harmony/sub_CK.rds")


# SAVE FIGURES 
p <- DimPlot(CK.subset, reduction = "umap", label = T,   
             cols= allcolour, #设置颜色  
             pt.size = 0.001,#设置点的大小  
             repel = T)+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p + scale_color_d3("category20c")
ggsave("ck_subset.pdf",p, width = 8,height = 5)

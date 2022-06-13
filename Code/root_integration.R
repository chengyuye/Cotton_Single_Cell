library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(Seurat)
library(kableExtra)
library(devtools)
#install_github('welch-lab/liger')
library(rliger)
library(stringr)

# ara data
Ara <-  readRDS(file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/combined_Ara_annotated_doubletremoved.rds")
head(rownames(Ara))
# cotton data
Cotton <- readRDS(file = "C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds")
head(rownames(Cotton))
datasets <- c("Ara", "Cotton")

### rename the lineage annotation of cotton
### cell type annotation ###
cluster.annotation <- c("Unknown 1", 'Lateral root cap 1', 'Protophloem',
                        'Pericycle', 'Stele','Root cortex', 
                        'Non-hair root epidermal cell','Columella',  
                        'Root endodermis', 'Quiescent center', 
                        'Trichoblast', 'Columella root cap cell', 'Xylem',
                        'Lateral root cap 2', 'Phloem', 'Unknown 2')


names(cluster.annotation) <- levels(Cotton)
Cotton <- RenameIdents(Cotton, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(Cotton),
                        celltype = Idents(Cotton))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
Cotton <- AddMetaData(Cotton, cell.data, col.name = "celltype")


### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(Cotton),
                        celltype = Idents(Cotton))

lineage.annotation <- c("Unknown","Root cap cell (RCC)","VT",
                        "Stele","Stele","Cortex","Non root hair (NRH)",
                        "Columella","Endodermis","Root Meristem",
                        "Trichoblast","Root cap cell (RCC)","Stele(X)",
                        "Root cap cell (RCC)","Stele(P)", "Unknown")


lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
Cotton <- AddMetaData(Cotton, meta.data, col.name = "lineage")

### extracting umap infor for figures -----cell lineage
umap = Cotton@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_lineage = Cotton@meta.data$lineage) # 注释后的label信息 ，改为cell_type

head(umap)
## generate figures

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_lineage)) +
  geom_point(size = 0.1, alpha =1 ) 

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

p4 <- p4+	scale_color_d3('category20c')
##
ggsave('CK_annotated_lineage.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/物种整合part2/",
       width = 12, height = 8)

### rename the annotation of ara
### cell type annotation ###
cluster.annotation <- c("Root cap cell 1", 'Stele 1', 'Non-hair root epidermal cell 1',
                        'Root hair cell 1', 'Root cortex','Stele 2', 
                        'Unknown','Quiescent Center', 'Phloem 1', 
                        'Root hair cell 2', 'Non-hair root epidermal cell 2', 'Root cap cell 2',
                        'Meristem 1', 'Xylem 1', 'Root endodermis',
                        'Meristem 2', 'Phloem 2','Root cap cell 3', 
                        'Xylem 2')


names(cluster.annotation) <- levels(Ara)
Ara <- RenameIdents(Ara, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(Ara),
                        celltype = Idents(Ara))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode) 
cell.data$barcode <- NULL
Ara <- AddMetaData(Ara, cell.data, col.name = "celltype")



### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(Ara),
                        celltype = Idents(Ara))

lineage.annotation <- c("Root cap cell (RCC)", "Stele", "Non root hair (NRH)",
                        "Root hair (RH)", "Cortex", "Stele",
                        "Unknown", "Root Meristem","Stele(P)",
                        "Root hair (RH)","Non root hair (NRH)","Root cap cell (RCC)",
                        "Root Meristem", "Stele(X)", 'Endodermis',
                        'Root Meristem', "Stele(P)", 'Root cap cell (RCC)',
                        "Stele(X)")


lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
Ara <- AddMetaData(Ara, meta.data, col.name = "lineage")

### extracting umap infor for figures -----cell lineage
umap = Ara@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_lineage = Ara@meta.data$lineage) # 注释后的label信息 ，改为cell_type

head(umap)
## generate figures

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_lineage)) +
  geom_point(size = 0.05 , alpha =1 ) 

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

p4 <- p4+	scale_color_d3('category20c')
##
ggsave('Ara_annotated_lineage.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/物种整合part2/",
       width = 12, height = 8)

# overview of samples
data.frame(
  "Species" = c("Ara", "Cotton"),
  "Clusters" = map_int(datasets, ~ length(unique(get(.x)@active.ident)))
) %>%
  kable(caption = "Overview") %>%
  kable_styling(latex_options = "hold_position")



# plot_settings <- list(
#   theme_void(),
#   theme(legend.position = "none",
#         plot.title = element_text(hjust = 0.5, size = 10))
# )
# cowplot::plot_grid(plotlist = list(
#   DimPlot(mouse, label = TRUE) +
#     scale_color_manual(values = mouse@misc$cluster_colors) +
#     ggtitle("Mouse") +
#     plot_settings,
#   DimPlot(nhp, label = TRUE) +
#     ggtitle("Macaque") +
#     scale_color_manual(values = nhp@misc$cluster_colors) +
#     plot_settings
# ))

# keep metadata for future use
metadata <- map(
  datasets,
  ~ data.frame(
    "cell" = rownames(get(.x)@meta.data),
    "Cluster" = get(.x)@active.ident)
)
names(metadata) <- datasets



# Preprocessing

#Before I integrate the data, 
#I'm going to keep only the genes that are shared in both datasets.

# isolate just the raw count matrix files by species
mtx <- map(datasets, ~ get(.x)@assays$RNA@counts)
names(mtx) <- datasets


#print("Reading in ortholog info from Ensembl ...")
# read in ortholog info and only keep orthologs with 1:1 mapping
orthologs <- read.csv('C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/single-copy_orthologs.csv', header = F)

## filter na
orthologs <- orthologs %>% filter(V3 != '#N/A')
# get unique genes by species
unique_ara <- orthologs$V2
unique_ara <- str_sub(unique_ara,1,nchar(unique_ara)-2)
unique_cotton <- orthologs$V3
gene_num <- length(unique_ara) %>% as.numeric()

## replace the cotton name to ara name
for(i in 1:gene_num){
  rownames(mtx$Cotton) <- gsub(unique_cotton[i], unique_ara[i],rownames(mtx$Cotton))
  }

# rename the cell to avoid errors
colnames(mtx$Ara) <- paste("Ara", colnames(mtx$Ara), sep = "_")
colnames(mtx$Cotton) <- paste("Cot", colnames(mtx$Cotton), sep = "_")

###using liger for integration
species.liger <- createLiger(list(ara = mtx$Ara, Cotton = mtx$Cotton))

###Normalize the datasets. 
#The normalization is applied to the datasets in their entirety.
species.liger <- rliger::normalize(species.liger)
  
#Select shared, homologous genes between the two species, 
#as well as unshared, non-homologous genes
species.liger <- selectGenes(species.liger, var.thres= 0.1,unshared = TRUE, 
                             unshared.datasets = list(1,2), unshared.thresh= c(0.3,0.3),do.plot = T)

#scale, but do not center, the data.
species.liger <- scaleNotCenter(species.liger)


####Joint Matrix Factorization
species.liger <- optimizeALS(species.liger,  lambda = 5, 
                             use.unshared = TRUE, thresh=1e-10, k =30)

#Quantile Normalization and Joint Clustering
species.liger <- quantile_norm(species.liger, ref_dataset = "ara")
species.liger <- louvainCluster(species.liger, resolution = 0.1)

##Visualizations and Downstream processing
#UMAP
species.liger <- runUMAP(species.liger)

umap_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE) 
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/Ara_cotton_integrated_atlas_UMAP.pdf',
    width = 7, height = 5.5)
umap_plots[[1]] 
dev.off()

umap_plots[[2]] 

###using original cell annotation
Ara_anno <- Ara$lineage
Cotton_anno <- Cotton$lineage


names(Ara_anno) <- paste("Ara", names(Ara_anno), sep = "_")
names(Cotton_anno) <- paste("Cot", names(Cotton_anno), sep = "_")
anno <- append(Ara_anno, Cotton_anno)

umap_plots <- plotByDatasetAndCluster(species.liger, pt.size = 0.6, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE, 
                                      clusters = Ara_anno)

umap_plots <- plotByDatasetAndCluster(species.liger, pt.size = 0.6, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE, 
                                      clusters = anno)
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/atlas_annotated_by_Ara_UMAP.pdf',
    width = 10, height = 5.5)
umap_plots[[2]] 
dev.off()

umap_plots <- plotByDatasetAndCluster(species.liger, pt.size = 0.6, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE, 
                                      clusters = Cotton_anno)

pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/atlas_annotated_by_cotton_UMAP.pdf',
    width = 10, height = 6)
umap_plots[[2]] 
dev.off()
##TSNE
species.liger <- runTSNE(species.liger)

tsne_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("TSNE1","TSNE2"), 
                                      return.plots = TRUE) 
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/Ara_lot_rice_integrated_atlas_TSNE.pdf',
    width = 7, height = 5.5)
tsne_plots[[1]] 
dev.off()

##### remove Schiefelbein data for testing 
Ara <- subset(Ara, subset = Study != 'Schiefelbein')

## remove Cuperus and Rybel data for testing
# Ara <- subset(Ara, subset = Study != 'Cuperus')
# Ara <- subset(Ara, subset = Study != 'Rybel')


##### remove Rybel data for testing 
# Ara <- subset(Ara, subset = Study != 'Rybel') ### not good; should remain

##re- run umap for ara
Ara <- RunUMAP(Ara, reduction = "harmony_theta2", dims=1:50, seed.use=1)

##colors
# SAVE FIGURES  ----group by study and samples
library(ggsci)
p <- DimPlot(Ara, reduction = "umap", label = T,   
             #cols= allcolour, #??????色  
             pt.size = 0.01,#???玫??拇?小  
             group.by = 'Study',
             repel = T)+
  theme(axis.text.y = element_blank(),   #去??y???潭?注??
        axis.ticks.y = element_blank(),    #去??y???潭?
        axis.text.x = element_blank(),   #去??x???潭?注??
        axis.ticks.x = element_blank())+  #去??x???潭?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


#p <- p+ scale_color_npg()

ggsave('Ara_group_by_study_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 8, height = 6)

## group by celltype but split by study
p <- DimPlot(Ara, reduction = "umap", label = F,   
             #cols= allcolour, #??????色  
             pt.size = 0.01,#???玫??拇?小  
             group.by = 'celltype',
             split.by = 'Study',
             ncol = 3)+
  theme(axis.text.y = element_blank(),   #去??y???潭?注??
        axis.ticks.y = element_blank(),    #去??y???潭?
        axis.text.x = element_blank(),   #去??x???潭?注??
        axis.ticks.x = element_blank())+  #去??x???潭?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


p <- p+ scale_color_d3('category20')

ggsave('Ara_split_by_study_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 15, height = 8)

## group by samples
p <- DimPlot(Ara, reduction = "umap", label = F,   
             #cols= allcolour, #??????色  
             pt.size = 0.01,#???玫??拇?小  
             group.by = 'Sample',
             repel = T)+
  theme(axis.text.y = element_blank(),   #去??y???潭?注??
        axis.ticks.y = element_blank(),    #去??y???潭?
        axis.text.x = element_blank(),   #去??x???潭?注??
        axis.ticks.x = element_blank())+  #去??x???潭?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


#p <- p+ scale_color_lancet()

ggsave('Ara_group_by_sample_umap.pdf', p,  path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 8, height = 6)

## group by celltype but split by sample
p <- DimPlot(Ara, reduction = "umap", label = F,   
             #cols= allcolour, #??????色  
             pt.size = 0.01,#???玫??拇?小  
             group.by = 'celltype',
             split.by = 'Sample',
             ncol = 3)+
  theme(axis.text.y = element_blank(),   #去??y???潭?注??
        axis.ticks.y = element_blank(),    #去??y???潭?
        axis.text.x = element_blank(),   #去??x???潭?注??
        axis.ticks.x = element_blank())+  #去??x???潭?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


p <- p+ scale_color_d3('category20')

ggsave('Ara_split_by_sample_umap.pdf', path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 15, height = 12)


### extracting umap infor for figures---- celltypes
umap = Ara@reductions$umap@cell.embeddings %>%  #??????息
  as.data.frame() %>% 
  cbind(cell_type = Ara@meta.data$celltype) # 注?秃???label??息 ????为cell_type

### extracting umap infor for figures -----cell lineage
umap = Ara@reductions$umap@cell.embeddings %>%  #??????息
  as.data.frame() %>% 
  cbind(cell_lineage = Ara@meta.data$lineage) # 注?秃???label??息 ????为cell_type

head(umap)
## colors
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
## generate figures
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 0.01 , alpha =1 ) 

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_lineage)) +
  geom_point(size = 0.01 , alpha =1 ) 

p2 <- p  +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_blank(), #?呖?
        axis.title = element_blank(),  #??????
        axis.text = element_blank(), # ?谋?
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #????色
        plot.background=element_rect(fill="white"))

p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去??legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=12), #????legend??签?拇?小
    legend.key.size=unit(1,'cm') ) +  # ????legend??签之???拇?小
  guides(color = guide_legend(override.aes = list(size=2))) #????legend?? ???拇?小 


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
ggsave('Ara_annotated_umap_refined.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 12, height = 8)
##
ggsave('Ara_annotated_lineage_refined.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 12, height = 8)




# # remove unknown cluster
# Cotton <- subset(Cotton, subset = celltype != 'Unknown 1')
# Cotton <- subset(Cotton, subset = celltype != 'Unknown 2')
# #Cotton <- subset(Cotton, subset = celltype != 'Protophloem')
# Ara <- subset(Ara, subset = celltype != 'Unknown')


# isolate just the raw count matrix files by species
mtx <- map(datasets, ~ get(.x)@assays$RNA@counts)
names(mtx) <- datasets

## replace the cotton name to ara name
for(i in 1:gene_num){
  rownames(mtx$Cotton) <- gsub(unique_cotton[i], unique_ara[i],rownames(mtx$Cotton))
}

# rename the cell to avoid errors
colnames(mtx$Ara) <- paste("Ara", colnames(mtx$Ara), sep = "_")
colnames(mtx$Cotton) <- paste("Cot", colnames(mtx$Cotton), sep = "_")

###using liger for integration
species.liger <- createLiger(list(Ara = mtx$Ara, Cotton = mtx$Cotton))

###Normalize the datasets. 
#The normalization is applied to the datasets in their entirety.
species.liger <- rliger::normalize(species.liger)

#Select shared, homologous genes between the two species, 
#as well as unshared, non-homologous genes
species.liger <- selectGenes(species.liger, var.thres = 0.2,unshared = TRUE, 
                             unshared.datasets = list(2), unshared.thresh= 0.3,do.plot = T)

#scale, but do not center, the data.
species.liger <- scaleNotCenter(species.liger)


####Joint Matrix Factorization
species.liger <- optimizeALS(species.liger,  lambda = 5, 
                             use.unshared = TRUE, thresh=1e-10, k =30)

#Quantile Normalization and Joint Clustering
species.liger <- quantile_norm(species.liger, ref_dataset = "Ara")
species.liger <- louvainCluster(species.liger, resolution = 0.2)

##Visualizations and Downstream processing
#UMAP
species.liger <- runUMAP(species.liger)

umap_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE) 
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/Ara_cotton_integrated_atlas_withoutSchiefelbein_UMAP.pdf',
    width = 7, height = 5.5)
umap_plots[[1]] 
dev.off()

### cluster
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/Ara_cotton_integrated_atlas_withoutSchiefelbein_UMAP_cluster.pdf',
    width = 7, height = 5.5)
umap_plots[[2]] 
dev.off()

###using original cell annotation
Ara_anno <- Ara$lineage
Cotton_anno <- Cotton$lineage

names(Ara_anno) <- paste("Ara", names(Ara_anno), sep = "_")
names(Cotton_anno) <- paste("Cot", names(Cotton_anno), sep = "_")

## merge the annotation 
annotation <- append(Ara_anno, Cotton_anno)

umap_plots <- plotByDatasetAndCluster(species.liger, pt.size = 0.6, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE, 
                                      clusters = Ara_anno)

pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/atlas_annotated_by_Ara_UMAP.pdf',
    width = 10, height = 5.5)
umap_plots[[2]] 
dev.off()

umap_plots <- plotByDatasetAndCluster(species.liger, pt.size = 0.6, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE, 
                                      clusters = Cotton_anno)

pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/atlas_annotated_by_cotton_UMAP.pdf',
    width = 10, height = 6)
umap_plots[[2]] 
dev.off()

##all annotation
umap_plots <- plotByDatasetAndCluster(species.liger, pt.size = 0.6, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE, 
                                      clusters = annotation)

pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/atlas_annotated_UMAP.pdf',
    width = 10, height = 6)
umap_plots[[2]] 
dev.off()

##TSNE
species.liger <- runTSNE(species.liger)
species.liger <- louvainCluster(species.liger, resolution = 0.2)

tsne_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("TSNE1","TSNE2"), 
                                      return.plots = TRUE) 
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/Ara_cotton_integrated_atlas_withoutSchiefelbein_TSNE.pdf',
    width = 7, height = 5.5)
tsne_plots[[1]] 
dev.off()

### cluster
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/Ara_cotton_integrated_atlas_withoutSchiefelbein_TSNE_cluster.pdf',
    width = 7, height = 5.5)
tsne_plots[[2]] 
dev.off()
##annotation
tsne_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("TSNE1","TSNE2"), 
                                      return.plots = TRUE,
                                      clusters = annotation) 
pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/Ara_cotton_integrated_atlas_annotated_TSNE.pdf',
    width = 8, height = 5.5)
tsne_plots[[2]] 
dev.off()
### riverplot

pdf(file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/cluster_sankeyplot.pdf',
    width = 18, height = 8) 
makeRiverplot(species.liger, as.factor(Ara_anno), as.factor(Cotton_anno),river.usr = c(0.5,0.75,-0.6,1.6))
dev.off()

### findmarkers
new_markers <- getFactorMarkers(species.liger, dataset1 = 'Ara', dataset2 = 'Cotton', num.genes = 10)

## word clouds
word_clouds <- plotWordClouds(species.liger,num.genes = 10,do.spec.plot = F, return.plots = T)

# Create Seurat object from liger object, keeping liger highly variable genes
seurat_obj = ligerToSeurat(species.liger, use.liger.genes = T)

### extracting umap infor for figures -----cell lineage
umap = seurat_obj@reductions$tsne@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = annotation) # 注释后的label信息 ，改为cell_type

head(umap)
## generate figures
p <- ggplot(umap,aes(x= tSNE_1 , y = tSNE_2 ,color = cell_type)) +
  geom_point(size = 0.01 , alpha =1 ) 


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
  geom_segment(aes(x = min(umap$tSNE_1) , y = min(umap$tSNE_2) ,
                   xend = min(umap$tSNE_1) +3, yend = min(umap$tSNE_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$tSNE_1)  , y = min(umap$tSNE_2)  ,
                   xend = min(umap$tSNE_1) , yend = min(umap$tSNE_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$tSNE_1) +1.5, y = min(umap$tSNE_2) -1, label = "tSNE_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$tSNE_1) -1, y = min(umap$tSNE_2) + 1.5, label = "tSNE_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
#p4 <- p4+scale_color_d3("category20")

#p4 <- p4+		scale_color_d3('category20c')

##
ggsave('refined_Ara_cotton_integrated_atlas_annotated_TSNE.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 12, height = 8)

##another plot
seurat_obj$celltype <- annotation
Idents(seurat_obj) <- seurat_obj$celltype


##by cluster
p <- DimPlot(seurat_obj, reduction = "tsne", 
             label = T,  
             pt.size = 0.04,#设置点的大小  
             repel = T)+
  scale_color_manual(values = allcolour)+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

p <- p +scale_color_d3('category20')

ggsave('Ara_cotton_integrated_atlas_TSNE_by_cluster.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 8, height = 6)
## by celltype
p <- DimPlot(seurat_obj, reduction = "tsne", 
             label = F, group.by = 'celltype',  
             pt.size = 0.01,#设置点的大小  
             repel = T)+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))



ggsave('refined_Ara_cotton_integrated_atlas_annotated_TSNE.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/",
       width = 8, height = 6)

## compute cluster marker genes ###
integratied.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
top10.markers <- integratied.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10.markers, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/top10_marker_genes.csv")

### cell type annotation ###
cluster.annotation <- c("stele 1", 'columella', 'unknown 1',
                        'root cortex 1', 'root hair cell','lateral root cap', 
                        'Xylem 1','pericycle', 'Protophloem', 
                        'trichoblast', 'stele 2', 'root cortex 2',
                        'Phloem', 'Xylem 2', 'quiescent center',
                        'unknown 2', 'root cortex 3','Root endodermis')


names(cluster.annotation) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(seurat_obj),
                        celltype = Idents(seurat_obj))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
seurat_obj <- AddMetaData(seurat_obj, cell.data, col.name = "celltype")

p <- DimPlot(seurat_obj, group.by = 'celltype', label = T, label.size = 3,pt.size = 0.03, repel = T)
p <- p + scale_color_d3('category20')

ggsave('Ara_annotated_umap_res0.5.pdf', p, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 8, height = 5)
### save data ### 34747 cells remained
saveRDS(seurat_obj, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/integrated_ligertoseurat.rds")

## extract stele and meristem
stele <- WhichCells(seurat_obj, idents = 'Stele')
meristem <- WhichCells(seurat_obj, idents = 'Root Meristem')

DimPlot(seurat_obj, label=T, group.by="celltype", pt.size = 0.01,
        cells.highlight= list(meristem), cols.highlight = c("darkblue"), cols= "grey")



## subset orginal liger object
species.liger_subset <- subsetLiger(species.liger, clusters.use = c(0,1,2,3,5,6,7,8,
                                                                    9,10))
## compute cluster marker genes ###
integratied.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
top10.markers <- integratied.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10.markers, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/top10_marker_genes.csv")



# Create Seurat object from liger object again
seurat_obj = ligerToSeurat(species.liger_subset, use.liger.genes = T)

### cell type annotation ###
cluster.annotation <- c("Root Meristem", 'Cortex', 'Root cap cell', 
                        'Xylem','Phloem','Trichoblast',
                        'Endodermis', 'Non root hair cell', 
                        'Unknown 1', 'Unknown 2')


names(cluster.annotation) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(seurat_obj),
                        celltype = Idents(seurat_obj))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
seurat_obj <- AddMetaData(seurat_obj, cell.data, col.name = "celltype")

## remove the 'removed'
seurat_obj <- subset(seurat_obj, subset = celltype != 'removed')

p <- DimPlot(seurat_obj, group.by = 'celltype', label = T, label.size = 4,pt.size = 0.05, repel = T)
p <- p + scale_color_d3('category20')+ xlab('UMAP_1')+ylab('UMAP_2')
ggsave('annotated_integrated_atlas_res0.1.pdf', p, path = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/',
       width = 8, height = 6)

### splity by species
p <- DimPlot(seurat_obj, split.by = 'orig.ident', label = T, 
             label.size = 4, pt.size = 0.05, repel = T)

### umap
p <- DimPlot(seurat_obj, split.by = 'orig.ident', 
             group.by = 'orig.ident',pt.size = 0.05) +
  scale_color_npg()
ggsave('integrated_atlas_splitby_species.pdf', p, path = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/results_version2/',
       width = 12, height = 6)

###TSNE
p <- DimPlot(seurat_obj, split.by = 'orig.ident', 
             group.by = 'orig.ident',pt.size = 0.05) +
  scale_color_npg()
ggsave('integrated_atlas_splitby_species_TSNE.pdf', p, path = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version1/',
       width = 12, height = 6)


### check ara marker genes in intergrated object
Idents(seurat_obj) <- seurat_obj$celltype

# root cap cell (RCC) genes
RCC.genes <- c("AT1G33280","AT1G79580")

# root meristematic cell (RMC) genes 
RMC.genes <- c("AT1G73590",'AT3G20840','AT2G04025','AT3G25980','AT5G13840')

# QC
QC.genes <- c('AT1G02870')

# epidermal cell/root hair (RH) genes
RH.genes <- c('AT5G49270','AT4G22080','AT2G03720')

# epidermal cell/non root hair (NRH) genes
NRH.gene <- c('AT3G02240','AT2G37260','AT1G65310','AT1G79840')

# cortex (Co) genes
RC.genes <- c("AT1G62510",'AT5G02000','AT5G64620')

# endodermis (CS+) genes
RE.genes <- c('AT5G57620','AT4G20140','AT5G06200','AT3G11550','RALF1')

# phloem (P)
phloem.genes <- c('AT4G19840','AT1G05760','AT2G15310','AT1G79430')

# stele/vacular cell genes 
VC.genes <- c('AT3G25710','AT5G48070','AT2G31083',
              'AT1G22710','AT3G23430','AT1G32450')

# xylem (X)
xylem.genes <- c('AT1G68810','AT1G20850',
                 'AT5G44030','AT3G16920','AT2G37090')

#trichoblast
trichoblast.gemes <- c('AT5G49270','AT1G63450','AT2G39690','AT1G07795')

##procabium
procabium.genes <- c('AT5G57130', 'AT4G32880')

###pericyle
pericyle.genes <- c('AT1G32450','AT5G01740')

markers <- c(RCC.genes, RMC.genes,  RH.genes, 
             NRH.gene, RC.genes, RE.genes, phloem.genes,
             VC.genes, xylem.genes, QC.genes)

## check which markers are in orthologs
shared_markers <- intersect(markers, unique_ara)

###dot plot
p <- DotPlot(seurat_obj, features = shared_markers) +
  coord_flip() +
  theme_bw()+ #  
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+ 
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3)) 

ggsave('Ara_dotplot_markers.pdf', p, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 6, height = 8)

### get shared orthologs and unshared genes together
unshared_genes <- species.liger@var.unshared.features$Cotton
all_genes <- append(unique_ara, unshared_genes)

# keep orthologs
# orthologs <- orthologs %>% filter(
#   `Gene name` %in% unique_mouse & `Gene name_1` %in% unique_nhp
# )

### only keep these genes
seurat_obj <- subset(seurat_obj, features = unique_ara)

###dot plot
p <- DotPlot(seurat_obj, features = shared_markers) +
  coord_flip() +
  theme_bw()+ #  
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+ 
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3)) 


seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@meta.data$celltype, 
                                                seurat_obj$orig.ident,sep = "_") # user adaptation required on own dataset
### caculate average expression level
ara_cot_Exp = AverageExpression(seurat_obj, group.by = "celltype_aggregate")
ara_cot_Exp = ara_cot_Exp[[1]]
ara_cot_Exp = ara_cot_Exp[rownames(ara_cot_Exp)%in%unique_ara,]
head(ara_cot_Exp)

### split into cotton and ara
ara_exp <- ara_cot_Exp[,c(2,4,6,8,10,11,13,15,17,20)]
ara_cot <- ara_cot_Exp[,c(1,3,5,7,9,12,14,16,18,19,21,22)]

#Prune Expression Tables & Remove rows with no expression
#Sp1 = ExpressionTableSpecies1[rownames(ExpressionTableSpecies1) %in% DEgenesSpecies1,]
ara_exp = ara_exp[rowSums (ara_exp)!=0,]
#Sp2 = ExpressionTableSpecies2[rownames(ExpressionTableSpecies2) %in% DEgenesSpecies2,]
ara_cot = ara_cot[rowSums (ara_cot)!=0,]

#Scale Expression Tables by gene average
avg = rowMeans(ara_exp)
ara_exp = sweep(ara_exp,1,avg,"/")
rm(avg)
avg = rowMeans(ara_cot)
ara_cot = sweep(ara_cot,1,avg,"/")
rm(avg)

#Step6: Merge Expression Tables
geTable = merge(ara_exp,ara_cot, by='row.names', all=F)
rownames(geTable) = geTable$Row.names
geTable = geTable[,2:ncol(geTable)]

# # calculate correlation
# cor_ara_cot_Exp <- cor(ara_cot_Exp)

#Step7:  Correlation
#7a:  Correlation
corr.method <- 'pearson'
Corr.Coeff.Table = cor(geTable,method=corr.method)

#7b:  Shuffle data
shuffled.cor.list = list()
pb <- txtProgressBar(1, 100, style=3)

nPermutations <- 1000
for (i in 1:nPermutations){
  shuffled = apply(geTable[,1:ncol(ara_exp)],1,sample)
  shuffled2 = apply(geTable[,(ncol(ara_exp)+1):ncol(geTable)],1,sample)
  shuffled = cbind(t(shuffled),t(shuffled2))
  shuffled.cor = cor(shuffled,method=corr.method)
  shuffled.cor.list[[i]] = shuffled.cor
  rm(list=c('shuffled','shuffled2','shuffled.cor'))
  if ((i %% 100) ==0){
    setTxtProgressBar(pb, (i*100)/nPermutations)
  }
}

p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
rownames(p.value.table) = colnames(geTable)
colnames(p.value.table) = colnames(geTable)

shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
rownames(shuffled.mean.table) = colnames(geTable)
colnames(shuffled.mean.table) = colnames(geTable)

a = combn(1:ncol(geTable),2)
for (i in 1:ncol(a)){
  cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
  shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
  shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
  p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
  p.value.table[a[1,i],a[2,i]] = p.value
  p.value.table[a[2,i],a[1,i]] = p.value
  rm(list=c('cor.scores','p.value'))
  setTxtProgressBar(pb, (i*100)/ncol(a))
}

neg.log10.p = -log10(p.value.table)

#step8 "Overlap in Markers"
#for all pairs of cell-types generate list of genes that are at least 1.5x avg in both cells

#from above a = combn(1:ncol(geTable),2)
marker.overlap.list = list()
for (i in 1:ncol(a)){
  datasubset = cbind(geTable[,a[1,i]],geTable[,a[2,i]])
  markers = rownames(geTable[datasubset[,1]>1.5 & datasubset[,2]>1.5,])
  marker.overlap.list[[i]] = markers
  names(marker.overlap.list)[i] = paste(colnames(geTable)[a[1,i]], colnames(geTable)[a[2,i]],sep='_')
  rm(list=c('datasubset','markers'))
}

list.to.return = list(Corr.Coeff.Table,shuffled.mean.table,p.value.table,neg.log10.p,DEgenes,rownames(geTable),length(DEgenes),length(rownames(geTable)),nDESp1,nDESp2,geTable,marker.overlap.list)
names(list.to.return) = c('corr.coeff','shuffled_correlation_score_means','p.value','negative_log10_p.value','DEgenes_intersect','DEgenes_in_analysis','nDEgenes_intersect','nDEgenes_in_analysis','nDEgenes_Sp1','nDEgenes_Sp2','scaled_table','overlapping_markers')
return(list.to.return)

### make plots
comp_table.species1 <- t(Corr.Coeff.Table[1:ncol(ara_exp),
                                            (ncol(ara_exp)+1):nrow(Corr.Coeff.Table)])

p_table.species1 <- t(p.value.table[1:ncol(ara_exp),(ncol(ara_exp)+1):nrow(Corr.Coeff.Table)])
p_table.species1 <- (-p_table.species1)


col1 <- colorRampPalette(c("darkblue", "white","darkred"))

pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/AT_GH_Correlation.pdf",width = 6, height = 6)
corrplot(comp_table.species1, order="original",tl.pos="lt", 
         method="color", tl.col="black",
         cl.lim=c(min(comp_table.species1),max(comp_table.species1)), 
         is.corr=F,tl.cex=0.7, sig.level=(-0.05),insig="pch", pch=19, 
         pch.cex=0.25,pch.col="black",p.mat=p_table.species1, col=col1(200), 
         #main= paste(Method3,",",comp.species1[[8]], "genes", sep=" "),
         mar=c(3,1,5,1),cl.align.text="l")  
dev.off()


### find conserved markers between ara and cotton for root meristem
seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@active.ident, 
                                                seurat_obj$orig.ident,
                                                sep = "_") # user adaptation required on own dataset

### meristem DEG between ara and cotton
meristem_deg <- FindMarkers(seurat_obj, ident.1 = '5_Ara', ident.2 = '5_Cotton',
                            group.by = 'celltype_aggregate')

meristem_deg$cluster <- 'Meristem'
meristem_deg$geneID <- rownames(meristem_deg)
###


head(seurat_obj$celltype_aggregate)
seurat_obj1 <- subset(seurat_obj, celltype_aggregate == c('5_Ara','5_Cotton'))
Meristem_Exp = AverageExpression(seurat_obj1, group.by = "celltype_aggregate")
Meristem_Exp = Meristem_Exp[[1]] %>% as.data.frame()
Meristem_Exp_conserved = Meristem_Exp[rownames(meristem_marker),]
Meristem_Exp_conserved = as.matrix(Meristem_Exp_conserved)
ac=data.frame(cluster=new_celltype)
rownames(ac)=colnames(cs_data)
##complexheatmap
library(ComplexHeatmap)
library(paletteer)  
color = paletteer_d("ggsci::nrc_npg")[c(1,3)]
names(color) <- unique(new_celltype)
top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=color),
                                                 labels = unique(new_celltype),
                                                 labels_gp = gpar(cex=0.5,color='white',fontsize = 18)))
col_fun = colorRamp2(c(-2, 1, 4), c("#377EB8", "white", "#E41A1C"))
col_fun2 = colorRamp2(c(-2, 1, 4), c("#92b7d1", "white", "#d71e22"))
Heatmap(Meristem_Exp_conserved,
        #km = 4,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = T,
        column_split = c('5_Ara','5_Cotton'),
        top_annotation = top_anno, 
        column_title = NULL,
        heatmap_legend_param = list(
          title='Expression',
          title_position='leftcenter-rot'),
        #col = col_fun, 
          col = colorRampPalette(c("white","#66CCFF","#333366"))(100),
          border = "black",
        row_names_gp = gpar(fontsize = 8))

#### meristem deg
Meristem_Exp_deg = Meristem_Exp[rownames(meristem_deg),]
Meristem_Exp_deg = as.matrix(Meristem_Exp_conserved)


##### markers for all clusters
## compute cluster marker genes ###
head(seurat_obj$celltype_aggregate)

seurat_obj1 <- subset(seurat_obj, celltype_aggregate == c('5_Ara','5_Cotton'))


### find all markers
seurat_obj2 <- seurat_obj
### find conserved markers between ara and cotton for root meristem
seurat_obj2$celltype_aggregate = paste(seurat_obj$celltype, 
                                                seurat_obj$orig.ident,
                                                sep = "_") # user adaptation required on own dataset

Idents(seurat_obj2) <- seurat_obj2$celltype_aggregate

### remove the cluster that only cotton has
seurat_obj2 <- subset(seurat_obj2, subset = celltype != 'VT')
seurat_obj2 <- subset(seurat_obj2, subset = celltype != 'Trichoblast')
seurat_obj2 <- subset(seurat_obj2, subset = celltype != 'Columella')
seurat_obj2 <- subset(seurat_obj2, subset = celltype != 'Root hair (RH)')
##################
## determine the downsample number
nsample <- min(table(seurat_obj2$celltype_aggregate))
### sampling
seurat_obj3 <- subset(seurat_obj2, downsample = nsample)
## normalize and scaling
seurat_obj3 <- NormalizeData(seurat_obj3)
seurat_obj3 <- ScaleData(seurat_obj3,features = rownames(seurat_obj3))


Idents(seurat_obj3) <- seurat_obj3$celltype
ara_cot.markers <- FindAllMarkers(seurat_obj3,only.pos = TRUE, 
                                      min.pct = 0.25, logfc.threshold = 0.25) 
##################
exp <- GetAssayData(seurat_obj3, slot = 'scale.data')
new_celltype <- sort(seurat_obj3$celltype_aggregate)
#new_celltype <- new_celltype[c('5_Ara','5_Cotton')%in%new_celltype]
head(new_celltype)
cs_data <- as.matrix(exp[ara_cot.markers$gene,names(new_celltype)])
#cs_data <- as.matrix(exp[rownames(meristem_marker),names(new_celltype)])
#head(cs_data)
mean <- matrix(nrow = 36, ncol = 2)
rownames(mean) <- rownames(meristem_marker)
colnames(mean) <- c('5_Ara', '5_Cotton')
mean_value <- as.data.frame(rowMeans(cs_data))


ac=data.frame(cluster=new_celltype)
rownames(ac)=colnames(cs_data)
library(pheatmap)

pheatmap(cs_data,
         show_colnames =F,
         show_rownames = F,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col=ac,
         border_color = NA)
##complexheatmap


pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/Ara_Cot_deg_heatmap_common_clusters.pdf",
    width = 7, height = 5)
plot_heatmap(dataset = seurat_obj3, 
             markers = ara_cot.markers$gene,
             row_font_size = 0,
             sort_var = c("celltype","orig.ident"),
             anno_var = c("celltype","orig.ident"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("#0570b0", "#d7301f"), # color vector
                                "Reds"))
             #hm_colors = c('#4575B4','#FFFFBF','#D73027'))
dev.off()

##complexheatmap
library(ComplexHeatmap)
library(paletteer)  
color = paletteer_d("ggsci::category20")[c(-19,-20)]
names(color) <- unique(new_celltype)
top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=color),
                                                 labels = unique(new_celltype),
                                                 labels_gp = gpar(cex=0.5,color='white',fontsize = 18)))
col_fun = colorRamp2(c(-2, 1, 4), c("#377EB8", "white", "#E41A1C"))
col_fun2 = colorRamp2(c(-2, 1, 4), c("#92b7d1", "white", "#d71e22"))
Heatmap(cs_data,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = F,
        column_split = c('5_Ara','5_Cotton'),
        top_annotation = top_anno, 
        column_title = NULL,
        heatmap_legend_param = list(
          title='Expression',
          title_position='leftcenter-rot'),
        #col = col_fun, 
        col = colorRampPalette(c("white","#66CCFF","#333366"))(100),
        border = "black",
        row_names_gp = gpar(fontsize = 8))




###
celltypes <- unique(seurat_obj3$celltype)
diffgenes <- data.frame()
Idents(seurat_obj3) <- seurat_obj3$celltype_aggregate
for (i in celltypes){ #or however many clusters you have
    ident1 <- paste0(i,"_Cotton")
    ident2 <- paste0(i,"_Ara")
    ara_cot.diffgenes <- FindMarkers(seurat_obj3, ident.1 = ident1, ident.2=ident2, 
                                     min.pct = 0.1,logfc.threshold=2.3)
    ara_cot.diffgenes$Cluster <- as.character(i)
    diffgenes <- rbind(diffgenes, ara_cot.diffgenes)
}

diffgenes$gene <- rownames(diffgenes)


#### calculate the gene expression difference
all_Exp = AverageExpression(seurat_obj3, group.by = "celltype_aggregate")
all_Exp = all_Exp[[1]] %>% as.data.frame()
ara_exp <- all_Exp[,c(1,3,5,7,9,11,13,15,17)]
cot_exp <- all_Exp[,c(2,4,6,8,10,12,14,16,18)]

diff_exp <- abs(cot_exp - ara_exp) 


diff_exp <- diff_exp[diffgenes$gene,] %>% as.matrix()
diff_exp <- na.omit(diff_exp)

colnames(diff_exp) <- c("Cortex","Endodermis","Non root hair (NRH)",
                        "Root cap cell (RCC)","Root Meristem","Stele(P)",          
                        "Stele(X)","Stele","Unknown")

write.csv(diff_exp, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/DEG_diff.csv')

pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/Ara_Cot_deg_expression_change_heatmap.pdf",
    width = 6, height = 5)
Heatmap(diff_exp,
        show_row_names = F,
        col =col_fun,
        heatmap_legend_param = list(
          title='Gene Expression Difference',
          title_position='leftcenter-rot'))
dev.off()

### kegg of these divergent genes
library(clusterProfiler)
library(org.Hs.eg.db)

divergent_genes <- rownames(diff_exp)

## KEGG ----- up-regulated genes
divergent_genes <- unique(divergent_genes)
kegg <- enrichKEGG(gene = divergent_genes,
                    organism = "ath",
                    pvalueCutoff = 0.5,
                    qvalueCutoff = 0.5)
write.csv(kegg, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/deg_KEGG.csv')
p <- clusterProfiler::dotplot(kegg,showCategory = 10)
ggsave('deg_KEGG_dotplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/",
       width = 8, height = 5)

p <- clusterProfiler::barplot(kegg,,showCategory = 10)
ggsave('deg_KEGG_barplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/",
       width = 8, height = 5)

#### GO ENRICHMENT
GO = enrichGO(divergent_genes, OrgDb = "org.At.tair.db", 
              keyType="TAIR", 
              pvalueCutoff = 0.5,
              qvalueCutoff = 0.5)

write.csv(GO, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/deg_GO.csv')
p <- clusterProfiler::dotplot(GO,showCategory = 10)
ggsave('deg_GO_dotplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/",
       width = 12, height = 5)

p <- barplot(GO,,showCategory = 10)
ggsave('deg_GO_barplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparision_with_AT/result_version3/",
       width = 12, height = 5)

#####
heatmap(exp_dif)

Meristem_Exp_conserved = Meristem_Exp[rownames(meristem_marker),]
Meristem_Exp_conserved = as.matrix(Meristem_Exp_conserved)







#######
DoHeatmap(seurat_obj, features = rownames(meristem_deg),assay = 'RNA',
          group.by = 'orig.ident')

### stele DEG between ara and cotton
stele_deg <- FindMarkers(seurat_obj, ident.1 = '1_Ara', ident.2 = '1_Cotton',
                            group.by = 'celltype_aggregate')


DoHeatmap(seurat_obj, features = rownames(stele_deg),assay = 'RNA',
          group.by = 'orig.ident') + 
  scale_fill_gradientn(colours = c('navy','white','firebrick3'))


### meristem conserved markers
meristem_marker <- FindConservedMarkers(seurat_obj,grouping.var = 'orig.ident',
                                        ident.1 = 5, logfc.threshold = 0.1)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))


library(circlize)
DoHeatmap(seurat_obj, features = rownames(meristem_marker),assay = 'RNA',
          group.by = 'orig.ident') + 
  scale_fill_gradientn(colours = colorRamp2(c(-2,1,4),c('#377EB8','white','#E41A1C')))

### stele conserved markers
stele_marker <- FindConservedMarkers(seurat_obj,grouping.var = 'orig.ident',
                                     ident.1 = 1, logfc.threshold = 0.1)

DoHeatmap(seurat_obj, features = rownames(stele_marker),assay = 'RNA',
          group.by = 'orig.ident')


##更改主题
#添加相关系数
pheatmap::pheatmap(cor_ara_cot_Exp, cluster_rows = F, cluster_cols = F)
##
ggcorrplot(cor_ara_cot_Exp,hc.order=TRUE,outline.color="white",
           type="lower",colors = c("#6D9EC1", "white", "#E46726"),
           lab = FALSE)

# only keep genes that are in the dataset
unique_ara <- unique_ara[unique_ara %in% rownames(Ara)]
unique_cotton <- unique_cotton[unique_cotton %in% rownames(Cotton)]

# keep orthologs
orthologs <- orthologs %>% filter(
  `Gene name` %in% unique_mouse & `Gene name_1` %in% unique_nhp
)


p <- plot_heatmap(dataset = seurat_obj1, 
                  markers = rownames(meristem_marker),
                  sort_var = c("celltype_aggregate","orig.ident"),
                  anno_var = c("celltype_aggregate","orig.ident"),
                  anno_colors = list("Set2",                                             # RColorBrewer palette
                                     c("#0570b0", "#d7301f"), # color vector
                                     "Reds"))










































orthofinder = read.csv('./single-copy_orthologs.csv', header = F)
orthofinder = orthofinder[orthofinder[,5]!="",]
Ara_filter <- subset(x = Ara_root, features = orthofinder$X1)
length(rownames(Ara_filter))

# Lj??????
Lj_root = readRDS(file = "D:/singlecell-liang/write/result_new_2.rds")
orthofinder$X3 = gsub('_',"-",orthofinder$X3)
orthofinder$X4 = gsub('_',"-",orthofinder$X4)
Lj_filter <- subset(x = Lj_root, features = orthofinder$X3)
length(rownames(Lj_filter))

result=orthofinder[0,]
for (i in 1:length(Lj_filter@assays$RNA@counts@Dimnames[[1]])){
  results = subset(orthofinder,X3==Lj_filter@assays$RNA@counts@Dimnames[[1]][i])
  result=rbind(result,results)
}

Lj_filter@assays$RNA@counts@Dimnames[[1]] <- result$X1
Lj_filter@assays$RNA@data@Dimnames[[1]] <- result$X1
Lj_filter@assays$RNA@meta.features <- data.frame (row.names = rownames (Lj_filter@assays$RNA))

# 水????????
rice_root = readRDS(file = "D:/singlecell-gu/Nip_data/result/result.rds")
rice_filter <- subset(x = rice_root, features = orthofinder$X4)
length(rownames(rice_filter))

result=orthofinder[0,]
for (i in 1:length(rice_filter@assays$RNA@counts@Dimnames[[1]])){
  results = subset(orthofinder,X4==rice_filter@assays$RNA@counts@Dimnames[[1]][i])
  result=rbind(result,results)
}

rice_filter@assays$RNA@counts@Dimnames[[1]] <- result$X1
rice_filter@assays$RNA@data@Dimnames[[1]] <- result$X1
rice_filter@assays$RNA@meta.features <- data.frame (row.names = rownames (rice_filter@assays$RNA))



# RenameGenesSeurat <- function(obj, newnames ) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
#   print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
#   RNA <- obj@assays$RNA
#   
#   if (nrow(RNA) == length(newnames)) {
#     if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
#     if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
#   } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
#   obj@assays$RNA <- RNA
#   return(obj)
# }
# RenameGenesSeurat(obj = rice_filter, newnames = result$Arabidopsis_thaliana.TAIR10.pep.all)


merge1 = SplitObject(Ara_filter, split.by = "orig.ident")
merge2 = SplitObject(Lj_filter, split.by = "orig.ident")
merge3 = SplitObject(rice_filter, split.by = "orig.ident")
merge.list <- c(merge1$tair_root, merge2$`Cro-root`,merge3$Nip_root)


for (i in 1:length(merge.list)) {
  merge.list[[i]] <- NormalizeData(merge.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  merge.list[[i]] <- FindVariableFeatures(merge.list[[i]], selection.method = "vst",
                                             nfeatures = 1000, verbose = FALSE)
}

features = SelectIntegrationFeatures(object.list = merge.list,nfeatures=500)
#?FindIntegrationAnchors
pancreas.anchors <- FindIntegrationAnchors(object.list = merge.list , anchor.features = features, dims = 1:30 )  ## dims = 1:20,k.anchor = 5,k.filter = 30
#?IntegrateData
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:20)  ## dims = 1:20

DefaultAssay(pancreas.integrated) <- "integrated"        
pancreas.integrated <- ScaleData(pancreas.integrated,features=VariableFeatures(pancreas.integrated)) 
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))
ElbowPlot(pancreas.integrated)
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:50)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:50)

# pancreas.integrated <- RunTSNE(pancreas.integrated, dims = 1:16)
DimPlot(pancreas.integrated, reduction = "umap",label = F,group.by="orig.ident")
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident") 
DimPlot(pancreas.integrated, reduction = "umap") 

Idents(pancreas.integrated) = "celltype"
# pancreas.integrated@meta.data[pancreas.integrated$celltype == "Metaxylem",]$celltype = "xylem"
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident") 
DimPlot(pancreas.integrated, reduction = "umap",label = T, split.by  = "orig.ident") 

DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = c("Epidermis","Epidermis (near root hair)")))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = "xylem"))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = c("Root_cap","Root_hair","Cap&Root hair")))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = "Stele"))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = c("Root_cap")))



saveRDS(pancreas.integrated, file = "D:/singlecell-liang/write/cro_ara_merge_result/ara_Lj_merge_result.rds")
#pancreas.integrated = readRDS(file = "D:/singlecell-liang/write/cro_ara_merge_result/ara_Lj_merge_result.rds")

pbmc.markers <- FindAllMarkers(pancreas.integrated, 
                               only.pos = TRUE, 
                               min.pct = 0.5, 
                               min.diff.pct = 0.3,
                               logfc.threshold = 0.25) 


DefaultAssay(pancreas.integrated) <- "RNA" 
FeaturePlot(pancreas.integrated, features = c("AT5G05500", "AT1G30870", "AT1G12560"),max.cutoff = 'q80',keep.scale = "all")

################ ??群??????????图
pancreas.integrated = readRDS(file = "D:/singlecell-liang/write/cro_ara_merge_result/ara_Lj_merge_result.rds")
orthofinder = read.csv('D:/singlecell-liang/write/Orthogroup_result.csv')
orthofinder = orthofinder[orthofinder[,5]!="",]
list = list(orthofinder[1:50,"X1"])
object <- AddModuleScore(object = pancreas.integrated, features = list, name = "gene_type")
FeaturePlot(object = object, features = "gene_type1")
### top 20% of the cells with high expression levels           
object <- MetaFeature(object = pancreas.integrated, features = orthofinder[1:50,"X1"], meta.name = "Aggregate_Feature") 
data = object@meta.data$Aggregate_Feature
threshold = sort(data, decreasing = T)[0.2*length(data)]
object@meta.data$Aggregate_Feature = replace(data,data<threshold,0)
FeaturePlot(object = object, features = "Aggregate_Feature",max.cutoff = "q80")
################

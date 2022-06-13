###
library(Seurat)
library(patchwork)
library(ggplot2)
#devtools::install_github("hannet91/ggcor")
library(ggcor)
library(corrplot)
library(ggcorrplot)
library(tidyverse)
library(grid)
library(VennDiagram)
library(ggVennDiagram)
library("gridExtra")   

######## Load CK.harmony
CK.harmony <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds')
CK.harmony@meta.data$Sample = CK.harmony@meta.data$Condition # user adaptation required on own dataset
CK.harmony@meta.data$celltype_aggregate = paste(CK.harmony@meta.data$celltype, 
                                                CK.harmony@meta.data$Condition,sep = "_") # user adaptation required on own dataset
### save data ###
saveRDS(CK.harmony, "C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds")

### umap
### extracting umap infor for figures---- SAMPLE
umap = CK.harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(Sample = CK.harmony@meta.data$Sample) # 

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = Sample)) +
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
p4 <- p4 + scale_color_nejm()

### extracting umap infor for figures---- celltype
umap = CK.harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(celltype = CK.harmony@meta.data$celltype) # 

p5 <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = celltype)) +
  geom_point(size = 0.01 , alpha =1 )


p5 <- p5  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))

p5 <- p5 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=12), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=2))) #设置legend中 点的大小 


p5 <- p5 + 
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
p5 <- p5 + scale_color_d3('category20b')

p <- p4 | p5

##
ggsave('CK_annotated_umap_sample.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/",
       width = 16, height = 7)



### correlation of marker genes between CKG/CKN
CKG = subset(CK.harmony, Sample == "CKG")
CKN = subset(CK.harmony, Sample == "CKN")

CKG_exp = AverageExpression(CKG, group.by = "celltype_aggregate")
CKG_exp = CKG_exp$RNA

CKN_exp = AverageExpression(CKN, group.by = "celltype_aggregate")
CKN_exp = CKN_exp$RNA

CK_Exp = cbind(CKG_exp,CKN_exp)
CK_Exp = AverageExpression(CK.harmony, group.by = "celltype_aggregate")
CK_Exp = CK_Exp[[1]]
head(CK_Exp)


### corrplot
CK_Exp_mat <- as.matrix(CK_Exp)

# calculate correlation
cor_CK_Exp <- cor(CK_Exp, method = 'spearman')

##计算p值
cor_p <- cor_pmat(cor_CK_Exp)


##更改主题
#添加相关系数
ggcorrplot(cor_CK_Exp,hc.order=TRUE,outline.color="white",
           type="lower",colors = c("#6D9EC1", "white", "#E46726"),
           lab = FALSE)
ggcorrplot(cor_CK_Exp,
           hc.order = TRUE, type = "lower",
           outline.color = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726")
)

p <- quickcor(CK_Exp,type = "lower") + geom_square()+
  scale_fill_gradient2(midpoint = 0.8, low = 'blue', mid = 'white', high = 'red', space = 'Lab')+
  guides(fill=guide_colorbar(title = "Pearson's r", order = 3))
##
ggsave('CKN_CKG_Correlation_heatmap.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/",
       width = 10, height = 7)


pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/CKG_CKN_Correlation_spearman.pdf",width = 8, height = 6)
p <- pheatmap::pheatmap(cor(CK_Exp, method = 'spearman'), cluster_rows = F, cluster_cols = F)
##
dev.off()



### venn diagram for gene overlap between CKG/CKN
##remove genes with 0 expression -----CKN
zeros <- which(Matrix::rowSums(CKN) == 0)
CKN_genes.use <- rownames(CKN[-zeros,])

##remove genes with 0 expression -----CKG
zeros <- which(Matrix::rowSums(CKG) == 0)
CKG_genes.use <- rownames(CKG[-zeros,])


## CKG_exp
CKG_exp <- as.data.frame(CKG_exp)

## CKN_exp
CKN_exp <- as.data.frame(CKN_exp)

### CRC as a test
CRC_CKN <- CKN_exp %>% dplyr::select(1) %>% dplyr::filter(`Columella root cap cell_CKN`>0)
CRC_CKG <- CKG_exp %>% dplyr::select(1) %>% dplyr::filter(`Columella root cap cell_CKG`>0)
LRC1_CKN <- CKN_exp %>% dplyr::select(3) %>% dplyr::filter(`Lateral root cap 1_CKN`>0)
LRC1_CKG <- CKG_exp %>% dplyr::select(3) %>% dplyr::filter(`Lateral root cap 1_CKG`>0)

## OBTAIN THE NUMBER
CRC_CKN.genes <- rownames(CRC_CKN) 
CRC_CKG.genes <- rownames(CRC_CKG)
intersce <- intersect(CRC_CKN.genes, CRC_CKG.genes)

## convert to number
CRC_CKN.genes <- length(CRC_CKN.genes) %>% as.numeric()
CRC_CKG.genes <- length(CRC_CKG.genes) %>% as.numeric()
intersce <- length(intersce) %>% as.numeric()

## test which number is larger
if(CRC_CKN.genes > CRC_CKG.genes ){
  a <- CRC_CKN.genes
  b <- CRC_CKG.genes
} else{
  a <- CRC_CKG.genes
  b <- CRC_CKN.genes
}

## caculate p value using hypergeometric test
p <- phyper(intersce-1, a, all_genes_num - a, b, lower.tail = F)


## VENN DIAGRAM
pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/CRC_venn.pdf",width = 4, height = 4)
p <- draw.pairwise.venn(area1 = CRC_CKN.genes,
                   area2 = CRC_CKG.genes,
                   cross.area = intersce,
                   category.names = c("Set 1" , "Set 2"),
                   col  = c("#3C5488FF", "#DC0000FF"), 
                   alpha = c(0.5, 0.5),
                   ind = FALSE,
                   scaled = FALSE,
                   euler.d = FALSE)
gridExtra::grid.arrange(gTree(children=p), bottom = "Columella root cap")
#grid.draw(p)
dev.off()


### define a function
marker_venn <- function(cell_CKN){
  
  celltype <- gsub('_CKN', '', cell_CKN)
  cell_CKG <- gsub('_CKN', '_CKG', cell_CKN) 
  
  ###
  gene_CKN <- CKN_exp %>% dplyr::select(cell_CKN)
  colnames(gene_CKN) <- 'celltype_ckn'
  # caculate the number of all the genes
  all_genes_num <- length(gene_CKN) %>% as.numeric()
  # remove CKN =0 
  gene_CKN <- gene_CKN %>% dplyr::filter(celltype_ckn > 0)
  
  gene_CKG <- CKG_exp %>% dplyr::select(cell_CKG)
  colnames(gene_CKG) <- 'celltype_ckg'
  gene_CKG <- gene_CKG %>% dplyr::filter(celltype_ckg > 0)
  
  ## OBTAIN THE NUMBER
  gene_CKN.use <- rownames(gene_CKN) 
  gene_CKG.use <- rownames(gene_CKG) 
  intersce <- intersect(gene_CKN.use, gene_CKG.use)
  
  ## convert to number
  gene_CKN.use <- length(gene_CKN.use) %>% as.numeric()
  gene_CKG.use <- length(gene_CKG.use) %>% as.numeric()
  intersce <- length(intersce) %>% as.numeric()
  
  ## test which number is larger
  if(gene_CKN.use > gene_CKG.use ){
    a <- gene_CKN.use
    b <- gene_CKG.use
  } else{
    a <- gene_CKG.use
    b <- gene_CKN.use
  }
  
  ## caculate p value using hypergeometric test
  p <- phyper(intersce-1, a, all_genes_num-a, b, lower.tail = F) %>% as.character()
  
  ## VENN DIAGRAM
  pdf(paste0("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/", celltype,"_venn.pdf"),
      width = 4, height = 4)
  p <- draw.pairwise.venn(area1 = gene_CKN.use,
                          area2 = gene_CKG.use,
                          cross.area = intersce,
                          category.names = c("Set 1" , "Set 2"),
                          col  = c("#3C5488FF", "#DC0000FF"), 
                          alpha = c(0.5, 0.5),
                          ind = FALSE,
                          scaled = FALSE,
                          inverted = T)
  gridExtra::grid.arrange(gTree(children = p), # Add title & subtitle
                          top=celltype,
                          bottom = p)
  #grid.draw(p)
  
  dev.off()
}

### celltype_cnk
celltypes_CKN <- unique(CKN$celltype_aggregate)
celltypes_CKN <- celltypes_CKN[which(!is.na(celltypes_CKN))]
celltypes_CKN <- setdiff(celltypes_CKN, c("red_blood_cells", "low_quality_cells", "Gamma-Delta_T_cells"))

#run for all celltypes
lapply(celltypes_CKN, function(x) marker_venn(x))


### KEGG for all cluster together using shared gene
### define a function
shared_gene <- function(cell_CKN){
  
  celltype <- gsub('_CKN', '', cell_CKN)
  cell_CKG <- gsub('_CKN', '_CKG', cell_CKN) 
  
  ###
  gene_CKN <- CKN_exp %>% dplyr::select(cell_CKN)
  colnames(gene_CKN) <- 'celltype_ckn'

  # remove CKN =0 
  gene_CKN <- gene_CKN %>% dplyr::filter(celltype_ckn > 0)
  
  gene_CKG <- CKG_exp %>% dplyr::select(cell_CKG)
  colnames(gene_CKG) <- 'celltype_ckg'
  gene_CKG <- gene_CKG %>% dplyr::filter(celltype_ckg > 0)
  
  ## OBTAIN the shared gene
  gene_CKN.use <- rownames(gene_CKN) 
  gene_CKG.use <- rownames(gene_CKG) 
  intersce <- intersect(gene_CKN.use, gene_CKG.use)
  
  ##returen
  intersce
}

#run for all celltypes
shared_gene_list <- lapply(celltypes_CKN, function(x) shared_gene(x))

## change the name
celltypes <- gsub('_CKN', '', celltypes_CKN)
names(shared_gene_list) <- celltypes

### modify the gene ID
delete_LOC <- function(list){
  list <- gsub('[LOC]', '', list)
}
  
## run for all the shared_gene_list
shared_gene_list <- lapply(shared_gene_list, function(x) delete_LOC(x))

### run grouped KEGG
### GO/KEGG analysis
library(clusterProfiler)
devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
library(clusterProfiler.dplyr)
kk <- compareCluster(shared_gene_list, 
                     fun="enrichKEGG",
                     organism="ghi", 
                     pvalueCutoff=0.05)

table(kk@compareClusterResult$Cluster) #每个基因集富集个数
head(as.data.frame(kk)) #查看完整结果
#Save the result
write.csv(kk, 
          file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/CKG_CKN_shared_gene_KEGG.csv')

## dotplot
sig.kegg<-filter(kk, p.adjust<0.01, qvalue < 0.05)

p <- clusterProfiler::dotplot(sig.kegg)
p <- p + scale_x_discrete(guide = guide_axis(angle = 45))
ggsave('CRC_KEGG_shared_gene_dotplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/comparison_CKN_CKG/",
       width = 12, height = 8)

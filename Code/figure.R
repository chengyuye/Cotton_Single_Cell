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

##### Figure 4 #####
## colors 
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
## add cultivar information to metadata 
CK.harmony$Cultivars <- CK.harmony$Condition
library(paletteer)
col <- c("darkolivegreen2","lightpink","lightblue2","#00F5FF","#FFA500","plum2","#FF6A6A","#7FFFD4", "#AB82FF")
col <- paletteer_d("ggsci::nrc_npg") #有几群细胞需要标记就选几种颜色  
p <- DimPlot(CK.harmony, reduction = "umap", group.by = 'Cultivars',
             label = T, 
             #cols= col, #设置颜色  
             pt.size = 0.1,#设置点的大小  
             repel = T)+
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())  #去掉x轴刻度
  #theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
ggsave('CK_cultivar_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/",
       width = 8, height = 6)





### Figure 6 Differentiation trajectories of lateral root cap
LRC <- subset(CK.harmony, subset = celltype == c('Lateral root cap 1',
                                                 'Lateral root cap 2',
                                                 'Lateral root cap 3'))###1126 cells

# recluster the lateral root cap
LRC <- NormalizeData(LRC, normalization.method = "LogNormalize", scale.factor = 1e4) 
LRC <- FindVariableFeatures(LRC, selection.method = 'vst', nfeatures = 2000)
LRC <- ScaleData(LRC)
LRC <- RunPCA(LRC, features = VariableFeatures(object = LRC)) 

## elbow plot - CK ##
pca.elbow.plot <- ElbowPlot(LRC, ndims = 50, reduction = "pca")
pca.elbow.plot

LRC <- FindNeighbors(LRC, dims = 1:50)
LRC <- FindClusters(LRC, resolution = 0.2 )
# run UMAP
LRC <- RunUMAP(LRC, dims=1:50, seed.use=1)

# visulization
#colors 
col <- c(brewer.pal(3, "Set1"))

p <- DimPlot(LRC, reduction = 'umap', 
             pt.size = 1.5, group.by = 'celltype',
             cols = col)
ggsave('LRC_recluster_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/figure6/",
       width = 8, height = 6)

## save data 
saveRDS(LRC, "C:/BGI/mianhua/CK_harmony/CK_LRC_annotated.rds")
## monocle 2
# read in data 
LRC <- readRDS("C:/BGI/mianhua/CK_harmony/CK_LRC_annotated.rds")
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(LRC@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- LRC@meta.data 
p_data$celltype <- LRC$celltype  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(LRC),row.names = row.names(LRC))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- as.CellDataSet(LRC)

#Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## ordering by marker gene per cluster
##使用monocle选择的高变基因??????
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1) 
cds <- setOrderingFilter(cds, disp.genes$gene_id)
plot_ordering_genes(cds)


## dimension reduciton
cds <- monocle::reduceDimension(cds, max_components = 2, method = 'DDRTree')

## ordering cells
cds <- orderCells(cds)

color1 <- c(brewer.pal(6, "Set1"))

# 6.1 celltype
p <- plot_cell_trajectory(cds,color_by="celltype", size=1,
                     show_branch_points = FALSE, show_backbone=TRUE)+
  scale_colour_manual(values = col)

ggsave('LRC_trajectory_by_celltype.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/figure6/",
       width = 6, height = 5)
colour=c('#91D1C2FF', '#E39A35', 'black')

p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "celltype")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 
ggsave('LRC_trajectory_by_celltype_treeplot.pdf', p2, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/figure6/",
       width = 6, height = 5)
#6.2 Pseudotime
p <- plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE,) + 
  scale_color_viridis_c()

ggsave('LRC_trajectory_by_pseudotime.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/figure6/",
       width = 6, height = 5)

#6.3 states
p <- plot_cell_trajectory(cds,color_by="State", size= 0.01,
                          show_branch_points = FALSE, show_backbone=TRUE)+
  scale_colour_manual(values = col)

ggsave('LRC_trajectory_by_state.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/figure6/",
       width = 6, height = 5)

#6.4 density
library(ggpubr)
df <- cds@phenoData@data
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(df)
p <- ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()

ggsavee('LRC_trajectory_cell_density.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
       width = 8, height = 5)


##6.5 color by genes
##2-ODD-1
colnames(pData(cds))
pData(cds)$LOC107942044 = log2( exprs(cds)['LOC107942044',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107942044",show_branch_points = FALSE)  +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

pData(cds)$LOC107942041 = log2(exprs(cds)['LOC107942041',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107942041",show_branch_points = FALSE)    +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
library(patchwork)
p3 <- p1+p2
ggsave('LRC_trajectory_2-ODD-1.pdf', p3, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
        width = 12, height = 6)

##DH1
pData(cds)$LOC107928135 = log2( exprs(cds)['LOC107928135',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107928135",show_branch_points = FALSE)  +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

pData(cds)$LOC107925530 = log2(exprs(cds)['LOC107925530',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107925530",show_branch_points = FALSE)    +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- p1+p2

ggsave('LRC_trajectory_DH1.pdf', p3, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
       width = 12, height = 6)

##CYP82D113
pData(cds)$LOC107903956 = log2( exprs(cds)['LOC107903956',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107903956",show_branch_points = FALSE)  +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

pData(cds)$LOC107944158 = log2(exprs(cds)['LOC107944158',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107944158",show_branch_points = FALSE)    +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- p1+p2

ggsave('LRC_trajectory_CYP82D113.pdf', p3, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
       width = 12, height = 6)

##CYP71BE79
pData(cds)$LOC107898600 = log2( exprs(cds)['LOC107898600',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107898600",show_branch_points = FALSE)  +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

pData(cds)$LOC107897460 = log2(exprs(cds)['LOC107897460',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107897460",show_branch_points = FALSE)    +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- p1+p2

ggsave('LRC_trajectory_CYP71BE79.pdf', p3, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
       width = 12, height = 6)

##CDNC
pData(cds)$LOC107903504 = log2( exprs(cds)['LOC107903504',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107903504",show_branch_points = FALSE)  +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

pData(cds)$LOC107903499 = log2(exprs(cds)['LOC107903499',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107903499",show_branch_points = FALSE)    +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- p1+p2

ggsave('LRC_trajectory_CDNC.pdf', p3, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
       width = 12, height = 6)

##CYP706B1
pData(cds)$LOC107920158 = log2( exprs(cds)['LOC107920158',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107920158",show_branch_points = FALSE) + scale_colour_gradient(low = "grey", high = "red")

pData(cds)$LOC107963672 = log2(exprs(cds)['LOC107963672',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107963672",show_branch_points = FALSE)    +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- p1+p2

ggsave('LRC_trajectory_CYP706B1.pdf', p3, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/trajectory/",
       width = 12, height = 6)


##gene expression 
#指定基因
s.genes <- c("LOC107942044",	"LOC107942041",
             "LOC107928135",	"LOC107925530",
             "LOC107903956",	"LOC107944158",
             'LOC107898600',	'LOC107897460',
             'LOC107903504',	'LOC107903499',
             'LOC107920158',	'LOC107963672')
cds_subset <- cds[s.genes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 16)

p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 16)

### branch point heatmap
BEAM_res <- BEAM(cds, branch_point = 3, cores = 2) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res, "BEAM_res.csv", row.names = F)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 3, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#


#选前100个基因可视化
BEAM_genes <- dplyr::top_n(BEAM_res, n = 100, dplyr::desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 3, 
                                 num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave("BEAM_heatmap.pdf", p$ph_res, width = 6.5, height = 10)

#显著差异基因(top100)按热图结果排序并保存
##如果要所有的差异基因，就把前面所632个基因的热图存为p
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig, "BEAM_sig.csv", row.names = F)




#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds, cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
sig_gene_names <- row.names(subset(Time_diff, qval < 0.05))
p=plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)


Time_diff_100 <- dplyr::top_n(Time_diff, n = 100, dplyr::desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_diff_100,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapTop100.pdf", p, width = 5, height = 10)

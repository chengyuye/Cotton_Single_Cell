library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
#library(dplyr)

### path variables ###
CK.path <- "C:/BGI/mianhua/CK_harmony_updated_test/CK.annotated_doubletremoved.h5ad"
LRC.path <- 'C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds'

### read files ###
CK.harmony <- readRDS(CK.path)


### Monocle3 ###
##创建CDS对象并预处理数据
data <- GetAssayData(CK.harmony, assay = 'RNA', slot = 'counts')
cell_metadata <- CK.harmony@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


### one setp to cds (monocle3)
cds <- as.cell_data_set(CK.harmony)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(CK.harmony, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p = p1|p2
ggsave("Reduction_Compare.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 10, height = 5)

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) +
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)


## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE,
               group_label_size = 4)+ scale_color_d3('category20')

ggsave("Trajectory.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 7, height = 6)

p <- plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
                label_leaves=TRUE, label_branch_points=TRUE,
                group_label_size = 4, graph_label_size = 2)+ scale_color_d3('category20')

ggsave("Trajectory_with_label.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 7, height = 6)



##细胞按拟时排序
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)

ggsave("Trajectory_Pseudotime.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 8, height = 6)

p = plot_cells(cds, color_cells_by = "pseudotime",
               show_trajectory_graph = FALSE,
               label_cell_groups = FALSE, 
               label_leaves = FALSE,  
               label_branch_points = FALSE)

ggsave("Trajectory_Pseudotime_no_line.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 8, height = 6)

saveRDS(cds, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/cds.rds")

p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(0,1,0.25))
embed <- data.frame(Embeddings(scRNAsub, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -6.75 & UMAP_1 < -6.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
# 解决order_cells(cds)报错"object 'V1' not found"
# rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)


p = plot_cells(cds, color_cells_by = "pseudotime", 
               label_cell_groups = FALSE, 
               label_leaves = FALSE,  
               label_branch_points = FALSE)
ggsave("Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)
saveRDS(cds, file = "cds.rds")


#5. 差异表达分析
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)



###Analyzing branches in single-cell trajectories
cds_subset <- choose_cells(cds)

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)



agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


### using subsetted LRC
### read files ###
LRC <- readRDS(LRC.path)


### Monocle3 ###
##创建CDS对象并预处理数据
data <- GetAssayData(LRC, assay = 'RNA', slot = 'counts')
cell_metadata <- LRC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


### one setp to cds (monocle3)
cds <- as.cell_data_set(LRC)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", 
                 cell_size = 0.8,
                 color_cells_by="seurat_clusters") + ggtitle('cds.umap')

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(LRC, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP",
                 cell_size = 0.5,
                 color_cells_by="celltype") + ggtitle('int.umap')
p = p1|p2
ggsave("LRC_Reduction_Compare.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 10, height = 5)

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) +
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)


## 识别轨迹
## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE,
               group_label_size = 4)+ scale_color_d3('category20')

ggsave("LRC_Trajectory.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 8, height = 6)

p <- plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
                label_leaves=TRUE, label_branch_points=TRUE,
                group_label_size = 4, graph_label_size = 2)+ scale_color_d3('category20')

ggsave("LRC_Trajectory_with_label.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 8, height = 6)

##细胞按拟时排序
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)

ggsave("Trajectory_Pseudotime.pdf", path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/",
       plot = p, width = 8, height = 6)



saveRDS(cds, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/monocle3/cds.rds")

p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(0,1,0.25))
embed <- data.frame(Embeddings(scRNAsub, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -6.75 & UMAP_1 < -6.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
# 解决order_cells(cds)报错"object 'V1' not found"
# rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)


p = plot_cells(cds, color_cells_by = "pseudotime", 
               label_cell_groups = FALSE, 
               label_leaves = FALSE,  
               label_branch_points = FALSE)
ggsave("Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)
saveRDS(cds, file = "cds.rds")


#5. 差异表达分析
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)












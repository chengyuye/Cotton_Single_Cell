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

### read in data 
LRC <- readRDS('C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds')
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
col <- paletteer_d("ggsci::nrc_npg") #?м?Ⱥϸ????Ҫ???Ǿ?ѡ??????ɫ  
p <- DimPlot(CK.harmony, reduction = "umap", group.by = 'Cultivars',
             label = T, 
             #cols= col, #??????ɫ  
             pt.size = 0.1,#???õ??Ĵ?С  
             repel = T)+
  theme(axis.text.y = element_blank(),   #ȥ??y???̶?ע??
        axis.ticks.y = element_blank(),    #ȥ??y???̶?
        axis.text.x = element_blank(),   #ȥ??x???̶?ע??
        axis.ticks.x = element_blank())  #ȥ??x???̶?
  #theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
ggsave('CK_cultivar_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony/dim50_annotation/",
       width = 8, height = 6)





### Differentiation trajectories of lateral root cap
## check the celltype
unique(CK.harmony$celltype)
LRC <- subset(CK.harmony, subset = celltype == c('Lateral root cap 1',
                                                 'Lateral root cap 2',
                                                 'Columella root cap cell'))###1126 cells
DimPlot(LRC, group.by = "celltype")

# recluster the lateral root cap
LRC <- NormalizeData(LRC, normalization.method = "LogNormalize", scale.factor = 1e4) 
LRC <- FindVariableFeatures(LRC, selection.method = 'vst', nfeatures = 2000)
LRC <- ScaleData(LRC)
LRC <- RunPCA(LRC, features = VariableFeatures(object = LRC)) 

## elbow plot - CK ##
pca.elbow.plot <- ElbowPlot(LRC, ndims = 50, reduction = "pca")
pca.elbow.plot

LRC <- FindNeighbors(LRC, dims = 1:50)
LRC <- FindClusters(LRC, resolution = 0.5 )
# run UMAP
LRC <- RunUMAP(LRC, dims=1:50, seed.use=1)

##dimplot
DimPlot(LRC)

# visulization
#colors 
col <- c(brewer.pal(3, "Set1"))

p <- DimPlot(LRC, reduction = 'umap', 
             pt.size = 1.5, group.by = 'celltype') + scale_color_npg()

ggsave('LRC_recluster_umap.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6.5, height = 5)
### extracting umap infor for figures---- celltypes
umap = LRC@reductions$umap@cell.embeddings %>%  #??????Ϣ
  as.data.frame() %>% 
  cbind(cell_type = LRC@meta.data$celltype) # ע?ͺ???label??Ϣ ????Ϊcell_type

head(umap)

## generate figures
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 1.5 , alpha =1 ) + 
  scale_color_manual(values = allcolour)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_blank(), #?߿?
        axis.title = element_blank(),  #??????
        axis.text = element_blank(), # ?ı?
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background=element_rect(fill="white"))

p3 <- p2 +         
  theme(
    legend.title = element_blank(), #ȥ??legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=12), #????legend??ǩ?Ĵ?С
    legend.key.size=unit(1,'cm') ) +  # ????legend??ǩ֮???Ĵ?С
  guides(color = guide_legend(override.aes = list(size=2))) #????legend?? ???Ĵ?С 


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
p4 <- p4+scale_color_npg()


ggsave('LRC_recluster_umap.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 7, height = 5)

####
pal= pal_npg()(10)
p <- ggplot(umap, aes(x= UMAP_1 , y = UMAP_2 ,fill = cell_type))+
  geom_point(size = 2.5,colour="grey35",shape=21) +
  scale_fill_manual(values = c("#DC0000FF","#4DBBD5FF","#00A087FF"))+ 
  theme_classic()

ggsave('LRC_recluster_umap_modified.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6.5, height = 5)

p <- DimPlot(LRC, reduction = 'umap', 
             pt.size = 1.0) + scale_color_nejm()

ggsave('LRC_recluster_umap_cluster.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6, height = 6)
## save data 
saveRDS(LRC, "C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds")

## monocle 2
# read in data 
LRC <- readRDS("C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds")
##??ȡ??????Ϣ--ϸ????Ϣ(????????ϸ???ľ???????ϸ?????ͼ?????Ϣ??ʵ??????????Ϣ)
expr_matrix <- as(as.matrix(LRC@assays$RNA@counts), 'sparseMatrix')
##??ȡ??????Ϣ??p_data(phenotype_data)???? 
p_data <- LRC@meta.data 
##??ȡ??????Ϣ ?????????͡?gc??��??
f_data <- data.frame(gene_short_name = row.names(LRC),row.names = row.names(LRC))
##expr_matrix????????f_data????????ͬ(gene number), expr_matrix????????p_data????????ͬ(cell number)

#????CDS????
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#??p_data??f_data??data.frameת??AnnotatedDataFrame??????
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())
# # ONE STEP FOR CONVERSION
# cds <- as.CellDataSet(LRC)

#Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

##filter
cds <- detectGenes(cds, min_expr = 0.1) #??һ????????fData(cds)??????һ??num_cells_expressed
print(head(fData(cds)))#z
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #??? ordering by marker gene per cluster
##ʹ??monocleѡ???ĸ߱?????
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table, mean_expression >= 0.1 )$gene_id 
# cds <- setOrderingFilter(cds, disp.genes)
# plot_ordering_genes(cds)

###dpfeature--- recommended?
diff_test_res_test <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~celltype")
#diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~celltype")

## select genes
#ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
ordering_genes <- row.names(diff_test_res_test)[order(diff_test_res_test$qval)][1:1000] #1321 genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

## dimension reduciton
cds <- monocle::reduceDimension(cds, max_components = 2, method = 'DDRTree')

## ordering cells
cds <- orderCells(cds)

#color1 <- c(brewer.pal(6, "Set1"))

# 6.1 celltype
p1 <- plot_cell_trajectory(cds,color_by="celltype",
                     show_branch_points = FALSE, show_backbone=TRUE)+
  scale_color_lancet()

ggsave('LRC_trajectory_by_celltype.pdf', p1, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6, height = 5)

### splity by CKG/CKN
p1 <- plot_cell_trajectory(cds,color_by="celltype",
                           show_branch_points = FALSE, show_backbone=TRUE) + 
  facet_wrap(~Condition, ncol  = 2)+
  scale_color_lancet()

ggsave('LRC_trajectory_celltype_CKG_CKN.pdf', p1, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 9, height = 4.5)


#colour=c('#91D1C2FF', '#E39A35', 'black')
p_2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "celltype")+
  scale_color_lancet() +
  theme(legend.title = element_blank()) 
ggsave('LRC_trajectory_by_celltype_treeplot.pdf', p_2, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6, height = 5)

#6.2 Pseudotime
###??ɫ1
temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)

p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE) + 
  scale_colour_gradientn(colours = rev(rbPal(50)), 
                         guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                          ticks.linewidth = 1, 
                                                          frame.colour = "black"))
### ??ɫ2
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE) + 
  scale_color_viridis_c()

ggsave('LRC_trajectory_by_pseudotime.pdf', p2, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 5.5, height = 4.5)

## split by CKG/CKN
###??ɫģʽ1

p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE) + 
  scale_color_viridis_c()+
  facet_wrap(~Condition, ncol  = 2)

###??ɫģʽ2
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE) + 
  scale_colour_gradientn(colours = rev(rbPal(50)), 
                       guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                        ticks.linewidth = 1, 
                                                        frame.colour = "black"))+
  facet_wrap(~Condition, ncol  = 2)

ggsave('LRC_trajectory_pseudotime_split_CKG_CKN.pdf', p2, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 9, height = 4.5)

#6.3 states
p3 <- plot_cell_trajectory(cds,color_by="State", cell_size=1,
                          show_branch_points = FALSE, show_backbone=TRUE)+
  scale_color_lancet()

ggsave('LRC_trajectory_by_state.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6, height = 5)

## split by CKG/CKN
p3 <- plot_cell_trajectory(cds,color_by="State",
                           show_branch_points = FALSE, 
                           show_backbone=TRUE)+
  scale_color_npg()+
  facet_wrap(~Condition, ncol = 2)

ggsave('LRC_trajectory_state_CKG_CKN.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 9, height = 4.5)

##?ϲ???ͼ
plotc <- p1|p2|p3
ggsave('LRC_trajectory_combined_plot.pdf', plotc, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 15, height = 5)

#6.4 density

df <- cds@phenoData@data
## pData(cds)ȡ??????cds??????cds@phenoData@data??????
View(df)
p <- ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
    geom_density(bw=0.5,size=1,alpha = 0.5) + theme_classic2() 

ggsave('LRC_trajectory_cell_density.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 8, height = 5)

df <- pData(cds)
ClusterName_color_panel <- c(
  "Lateral root cap 1" = "#6495ED", "Lateral root cap 2" = "#339933", 
  "Columellar root cap cell" = "#FF4500"
)
p <-ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ 
  scale_fill_manual(name = "", values = ClusterName_color_panel)+
  scale_color_manual(name = "", values = ClusterName_color_panel)

ggsave('LRC_trajectory_cell_density.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 8, height = 5)

##6.5 color by genes

## define a function
traj_gene <- function(cds, gene_list, path_name){
  # gradient_color <- c('#ffffcc','#ffeda0','#fed976','#feb24c',
  #                     '#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
  gene1 <- gene_list[1] %>% as.character()
  gene2 <- gene_list[2] %>% as.character()
  
  A <- paste0('pData(cds)','$', gene_list[1])
  A <-  log2(exprs(cds)[gene1,]+1)
  
  p1 <- plot_cell_trajectory(cds, color_by = gene1, show_branch_points = FALSE) +
    scale_color_gradient(low='grey', high = 'blue')+
    facet_wrap(~Condition, ncol  = 2)
  
  B <- paste0('pData(cds)','$', gene_list[2])
  B <-  log2(exprs(cds)[gene2,]+1)
  p2 <- plot_cell_trajectory(cds, color_by = gene2, show_branch_points = FALSE) +
    scale_color_gradient(low='grey', high = 'blue')+
    facet_wrap(~Condition, ncol  = 2)
  
  p3 <- p1+p2
  
  ggsave(paste0('LRC_trajectory_', path_name,'_CKG_CKN.pdf'), p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
         width = 12, height = 8)
}
##2-ODD-1
ODD <- c('LOC107942044','LOC107942041')
traj_gene(cds, ODD, '2-ODD-1')

### not using function
colnames(pData(cds))
pData(cds)$LOC107942044 = log2( exprs(cds)['LOC107942044',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107942044",show_branch_points = FALSE) +
  scale_colour_gradient(low='grey', high = 'blue')

pData(cds)$LOC107942041 = log2(exprs(cds)['LOC107942041',]+1)
p2 <- plot_cell_trajectory(cds, color_by = "LOC107942041",show_branch_points = FALSE) +
  scale_colour_gradient(low='grey', high = 'blue')


p3 <- p1+p2
ggsave('LRC_trajectory_2-ODD-1.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
        width = 12, height = 6)

### SPLIT BY CKN/CKG
pData(cds)$LOC107942044 = log2( exprs(cds)['LOC107942044',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107942044",show_branch_points = FALSE) +
  scale_colour_gradient(low='grey', high = "blue")+
  facet_wrap(~Condition, ncol  = 2)

pData(cds)$LOC107942041 = log2(exprs(cds)['LOC107942041',]+1)
p2 <- plot_cell_trajectory(cds, color_by = "LOC107942041",show_branch_points = FALSE) +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

p3 <- p1+p2
ggsave('LRC_trajectory_2-ODD-1_CKG_CKN.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 14, height = 5.5)

##DH1
DH1 <- c('LOC107925530','LOC107928135')
traj_gene(cds, DH1, 'DH1')

## not using function
pData(cds)$LOC107928135 = log2( exprs(cds)['LOC107928135',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107928135",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')

pData(cds)$LOC107925530 = log2(exprs(cds)['LOC107925530',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107925530",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')
p3 <- p1+p2

ggsave('LRC_trajectory_DH1.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 12, height = 6)

## split by CKG/CKN
pData(cds)$LOC107928135 = log2( exprs(cds)['LOC107928135',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107928135",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

pData(cds)$LOC107925530 = log2(exprs(cds)['LOC107925530',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107925530",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

p3 <- p1+p2

ggsave('LRC_trajectory_DH1_CKG_CKN.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 14, height = 5.5)

##CYP82D113
CYP82D113 <- c('LOC107903956','LOC107944158')
traj_gene(cds, CYP82D113, 'CYP82D113')

## not using function
pData(cds)$LOC107903956 = log2( exprs(cds)['LOC107903956',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107903956",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')

pData(cds)$LOC107944158 = log2(exprs(cds)['LOC107944158',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107944158",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')
p3 <- p1+p2

ggsave('LRC_trajectory_CYP82D113.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 12, height = 6)


### SPLIT BY CKG/CKN
pData(cds)$LOC107903956 = log2( exprs(cds)['LOC107903956',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107903956",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

pData(cds)$LOC107944158 = log2(exprs(cds)['LOC107944158',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107944158",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

p3 <- p1+p2

ggsave('LRC_trajectory_CYP82D113_CKN_CKN.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 14, height = 5.5)

##CYP71BE79
CYP71BE79 <- c('LOC107898600','LOC107897460')
traj_gene(cds, CYP71BE79, 'CYP71BE79')


## not using function
pData(cds)$LOC107898600 = log2( exprs(cds)['LOC107898600',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107898600",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')

pData(cds)$LOC107897460 = log2(exprs(cds)['LOC107897460',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107897460",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')
p3 <- p1+p2

ggsave('LRC_trajectory_CYP71BE79.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 12, height = 6)

### SPLIT BY CKG/CKN
pData(cds)$LOC107898600 = log2( exprs(cds)['LOC107898600',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107898600",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

pData(cds)$LOC107897460 = log2(exprs(cds)['LOC107897460',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107897460",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

p3 <- p1+p2

ggsave('LRC_trajectory_CYP71BE79_CKG_CKN.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 14, height = 5.5)

##CDNC
pData(cds)$LOC107903504 = log2( exprs(cds)['LOC107903504',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107903504",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')

pData(cds)$LOC107903499 = log2(exprs(cds)['LOC107903499',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107903499",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')
p3 <- p1+p2

ggsave('LRC_trajectory_CDNC.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 12, height = 6)

### SPLIT BY CKG/CKN
pData(cds)$LOC107903504 = log2( exprs(cds)['LOC107903504',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107903504",show_branch_points = FALSE)  +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

pData(cds)$LOC107903499 = log2(exprs(cds)['LOC107903499',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107903499",show_branch_points = FALSE)    +
  scale_colour_gradient(low='grey', high = 'blue')+
  facet_wrap(~Condition, ncol  = 2)

p3 <- p1+p2

ggsave('LRC_trajectory_CDNC_CKN_CKG.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 14, height = 5.5)

##CYP706B1
pData(cds)$LOC107920158 = log2( exprs(cds)['LOC107920158',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107920158",show_branch_points = FALSE) + 
  scale_colour_gradient(low = "grey", high = "blue")

pData(cds)$LOC107963672 = log2(exprs(cds)['LOC107963672',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107963672",show_branch_points = FALSE)    +
  scale_colour_gradient(low = "grey", high = "blue")
p3 <- p1+p2

ggsave('LRC_trajectory_CYP706B1.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 12, height = 6)


##SPLIT BY CKG/CKN
pData(cds)$LOC107920158 = log2( exprs(cds)['LOC107920158',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "LOC107920158",show_branch_points = FALSE) + 
  scale_colour_gradient(low = "grey", high = "blue")+
  facet_wrap(~Condition, ncol  = 2)

pData(cds)$LOC107963672 = log2(exprs(cds)['LOC107963672',]+1)
p2=plot_cell_trajectory(cds, color_by = "LOC107963672",show_branch_points = FALSE)    +
  scale_colour_gradient(low = "grey", high = "blue")+
  facet_wrap(~Condition, ncol  = 2)

p3 <- p1+p2

ggsave('LRC_trajectory_CYP706B1_CKN_CKG.pdf', p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 14, height = 5.5)


##gene expression 
#ָ??????
s.genes <- c("LOC107942044",	"LOC107942041",
             "LOC107928135",	"LOC107925530",
             "LOC107903956",	"LOC107944158",
             'LOC107898600',	'LOC107897460',
             'LOC107903504',	'LOC107903499',
             'LOC107920158',	'LOC107963672')
cds_subset <- cds[s.genes,]
##???ӻ?????state/celltype/pseudotime????
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3

ggsave("Genes_pseudotimeplot.pdf", plot = p3, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/Trajectory_CKG_CKN/",
       width = 6, height = 16)
ggsave("Genes_pseudotimeplot.pdf", plot = plotc, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 16, height = 16)

p1 <- plot_genes_jitter(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "celltype", color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 18, height = 16)



### branch point heatmap
BEAM_res <- BEAM(cds, branch_point = 1, cores = 1) 
#?????õ???ordergene??Ҳ???ǵ?????dpFeature?ҳ?��?Ļ?????
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #??2829??????????????????????
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/BEAM_res.csv", row.names = F)
### heatmap---all the genes
pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/BEAM_pseudotime_heatmap.pdf",width = 8, height = 8)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                            branch_point = 1, #???Ƶ????ĸ???֧
                            num_clusters = 4, #?ֳɼ???cluster????????Ҫ????
                            show_rownames = F)#

tmp1 <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 3, #??Щ???򱻷ֳɼ???group
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #Ĭ??ֵ
                                 #hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2?ֱ???ʲô??ɫ
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #?Ƿ񷵻?һЩ??Ҫ??Ϣ
)

pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/Trajectory_CKG_CKN/Branched_pseudotime_heatmap_3cluster.pdf",width = 5, height = 6)
tmp1$ph_res
dev.off()
### select cluster for kegg analysis
gene_group <- tmp1$annotation_row
gene_group$gene <- rownames(gene_group)


allcluster_kegg <- data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group <- filter(gene_group,gene_group$Cluster==i)
  
### modify the gene id
  gene_name <- gsub('LOC', '', small_gene_group$gene)
  
  kegg <- enrichKEGG(gene = gene_name,
                      organism = "ghi",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  # p <- clusterProfiler::dotplot(kegg,showCategory = 10)
  # ggsave(paste0('Branched_cluster', i, '_KEGG_dotplot.pdf'), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/Trajectory_CKG_CKN/",
  #        width = 10, height = 5)
  kegg_res <- kegg@result
  if (dim(kegg_res)[1] != 0) {
    kegg_res$cluster <- i
    allcluster_kegg <- rbind(allcluster_kegg,kegg_res)
  }
}

head(allcluster_kegg[,c("ID","Description","qvalue","cluster")])
### save the result 
write.csv(allcluster_kegg, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/Trajectory_CKG_CKN/Branched_3cluster_KEGG.csv')

p <- clusterProfiler::dotplot(allcluster_kegg$cluster,showCategory = 10)

#ѡǰ100?????????ӻ?
BEAM_genes <- dplyr::top_n(BEAM_res, n = 100, dplyr::desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  
                                 branch_point = 1, 
                                 num_clusters = 3, 
                                 show_rownames = T, 
                                 return_heatmap = T)

ggsave("BEAM_heatmap.pdf", p$ph_res, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 6.5, height = 10)

#????????????(top100)????ͼ???????򲢱???
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/BEAM_sig_top100.csv", row.names = F)




#?????ǰ
entialGeneTest(cds[ordering_genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,
 "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/Time_diff_all.csv", row.names = F)
### p < 0.001 heatmap
sig_gene_names <- row.names(subset(Time_diff, qval < 0.001))
p <- plot_pseudotime_heatmap(cds[sig_gene_names,], 
                             num_clusters=3, 
                             show_rownames=F, 
                             return_heatmap=T)
ggsave("Time_heatmap_p_less_than_0.001.pdf", p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/", width = 5, height = 10)
### selected pathway genes for heatmap
p <- plot_pseudotime_heatmap(cds[s.genes,], 
                             num_clusters=1, 
                             show_rownames=T, 
                             return_heatmap=T)
ggsave("Time_heatmap_p_selected_pathway_genes.pdf", p, 
       path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/", 
       width = 5, 
       height = 5)


Time_diff_100 <- dplyr::top_n(Time_diff, n = 100, dplyr::desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <-  plot_pseudotime_heatmap(cds[Time_diff_100,], 
                              num_clusters=3, 
                              show_rownames=T, 
                              return_heatmap=T)
ggsave("Time_heatmapTop100.pdf", p, 
       path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/",
       width = 5, 
       height = 8)
### transcriptional factor
#?????л???��??????????ʱ??????????
Time_diff <- differentialGeneTest(cds, cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

##top2000
top2000_time_diff_gene <- Time_diff$gene_short_name[1:2000]
write.csv(top2000_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/top2000_time_diff_gene.csv", row.names = F)

#top3000
top3000_time_diff_gene <- Time_diff$gene_short_name[1:3000]
write.csv(top3000_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/trajectory/top3000_time_diff_gene.csv", row.names = F)


### ʹ?ò?Ϊ0?Ļ???��??????ʱ??????????
#?????л???��??????????ʱ??????????
expr <- GetAssayData(object = LRC, slot = "counts")
zeros <- which(Matrix::rowSums(expr) == 0)
expr <- data.matrix(expr[-zeros,])
### obtain the gene names to use
genes.use <- rownames(expr)

## use these non-zero genes to infer the relationship with time
Time_diff <- differentialGeneTest(cds[genes.use,], cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

##top2000
top2000_time_diff_gene <- Time_diff$gene_short_name[1:2000]
write.csv(top2000_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/top2000_time_diff_non_zero_gene.csv", row.names = F)

##top2500
top2500_time_diff_gene <- Time_diff$gene_short_name[1:2500]
write.csv(top2500_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/top2500_time_diff_non_zero_gene.csv", row.names = F)


#top3000
top3000_time_diff_gene <- Time_diff$gene_short_name[1:3000]
write.csv(top3000_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/top3000_time_diff_non_zero_gene.csv", row.names = F)

### ʹ??֮ǰ???˵Ļ?????expr>0.1, mincell=10) ��??????ʱ??????????
### obtain the gene names to use
genes.use <- rownames(expr)

## use these previously filtered genes to infer the relationship with time
Time_diff <- differentialGeneTest(cds[expressed_genes,], cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

##top2000
top2000_time_diff_gene <- Time_diff$gene_short_name[1:2000]
write.csv(top2000_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/top2000_time_diff_filtered_gene.csv", row.names = F)

##top2500
top2500_time_diff_gene <- Time_diff$gene_short_name[1:2500]
write.csv(top2500_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/top2500_time_diff_filtered_gene.csv", row.names = F)


#top3000
top3000_time_diff_gene <- Time_diff$gene_short_name[1:3000]
write.csv(top3000_time_diff_gene, "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/top3000_time_diff_filtered_gene.csv", row.names = F)

## identified 142 TF using top3000
## TFN using SCODE
#first extract gene expression matri
## read in the 142 TF
TF <- read.table("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/filtered/TF_sorted.txt")
colnames(TF) <- c('TF_ID','Gene_ID','XP_ID','LOC_ID',	'Family','Sort_index')

TF_ID <- TF$LOC_ID

F
p <-  plot_pseudotime_heatmap(cds[TF_ID,], 
                              num_clusters=5, 
                              show_rownames=T, 
                              return_heatmap=T)
ggsave("TF_Pseudotime.pdf", p, 
       path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/",
       width = 6, 
       height = 8)


## generate heatmap based on branch
p <- plot_genes_branched_heatmap(cds[TF_ID,],
                                    branch_point = 1,
                                    num_clusters = 5, #??Щ???򱻷ֳɼ???group
                                    cores = 1,
                                    branch_labels = c("Cell fate 1", "Cell fate 2"),
                                    #hmcols = NULL, #Ĭ??ֵ
                                    #hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                    branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2?ֱ???ʲô??ɫ
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T #?Ƿ񷵻?һЩ??Ҫ??Ϣ
)

pdf("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/Branched_TF_pseudotime_heatmap.pdf",width = 6, height = 8)
p$ph_res
dev.off()



### generate plots to see where 
TFime_point <- data.frame(
  cell_index = cds$celltype,
  pseudotime = cds$Pseudotime)

## change the celltype
pseudotime_point$cell_index <- gsub("Lateral root cap 1", as.numeric(1), 
                                    pseudotime_point$cell_index)

pseudotime_point$cell_index <- gsub("Lateral root cap 2", as.numeric(2), 
                                    pseudotime_point$cell_index)

pseudotime_point$cell_index <- gsub("Columella root cap cell", as.numeric(3), 
                                    pseudotime_point$cell_index)
## normalize the pseudotime value
p <- pseudotime_point$pseudotime
pseudotime_point$pseudotime <- (p-min(p))/(max(p)-min(p))

##export pseudotime file
pseudotime_point <- as.matrix(pseudotime_point)
write.table(pseudotime_point,
            file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/pseudotime_point.txt", sep = "\t",
            row.names = FALSE)


###GENIE3
# Gencol.names = F, es that are used as # genSmoothCurves to Fit smooth spline curves and return the response matrix
cds_subset <- cds[TF_ID,]
newdata <- data.frame(Pseudotime = cds_subset$Pseudotime)
newdata <- (newdata-min(newdata))/(max(newdata)-min(newdata))

TF_Matrix <- genSmoothCurves(cds_subset, new_data = newdata, cores = 4)
TF_Matrix <- as.matrix(TF_Matrix)
#export matrx
write.table(TF_Matrix, col.names = FALSE, row.names = FALSE, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/TF_Analysis/TF_matrix.txt',sep = "\t")
#candidate regulators
library(GENIE3)
regulators <- TF_ID

# For reproducibility of results
set.seed(123)
weightMat <- GENIE3(expr_matrix, regulators=regulators, nCores=8, verbose=TRUE)


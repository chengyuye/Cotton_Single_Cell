### DA Analysis
library(plyr)
library(ggplot2)
library(scater)
library(ggpubr)
library(ggsci)
#rm(list=ls())

###read in data
CK.harmony <- readRDS("C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds")
## umap difference
### extracting umap infor for figures---- celltypes
umap = CK.harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = CK.harmony@meta.data$celltype,
        group = CK.harmony$Condition) # 注释后的label信息 ，改为cell_type

## generate figures
p <- ggplot(umap,aes(x= UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.1 , alpha =1) + 
  scale_color_manual(values = allcolour)+
  facet_wrap(~group,nrow=1)

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
p4 <- p4+ scale_color_d3("category20")

p4 <- p4+	scale_color_d3('category20c')


library(ggrepel)

p4 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                      point.padding=unit(0.5, "lines"))
##
ggsave('CK_annotated_umap_split_by_group.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DA/",
       width = 15, height = 8)



#theme函数设置
theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          # axis.ticks = element_line(color='black'),#坐标轴刻度线
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),#去除图例标题
          # legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
          legend.position=c(0.28, 0.9),#图例在绘图区域的位置
          # legend.position='top',#图例放在顶部
          legend.direction = "horizontal",#设置图例水平放置
          # legend.spacing.x = unit(2, 'cm'),
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
          legend.background = element_rect( linetype="solid",colour ="black")
          # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
          # legend.box.margin =margin(-10,0,0,0)
    )
  
}



metadata <- CK.harmony@meta.data
## colors
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
table(CK.harmony$Condition)#查看各组细胞数
# CKG   CKN 
# 10477  8802
prop.table(table(Idents(CK.harmony)))
table(Idents(CK.harmony), CK.harmony$Condition)#各组不同细胞群细胞数
#                               CKG  CKN
# Unknown 1                    1301 1436
# Lateral root cap 1            387 2013
# Protophloem                  1309  681
# Pericycle                    1433  552
# Stele                        1323  575
# Root cortex                   433 1137
# Non-hair root epidermal cell 1003  475
# Columella                    1048  289
# Root endodermis               551  355
# Quiescent center              437  326
# Trichoblast                   322  332
# Columella root cap cell       417  147
# Xylem                         308  242
# Lateral root cap 2             64  103
# Phloem                         91   75
# Unknown 2                      50   64
Cellratio <- prop.table(table(Idents(CK.harmony), CK.harmony$Condition), margin = 2)#计算各组样本不同细胞群比例
Cellratio
#                                      CKG         CKN
# Unknown 1                    0.124176768 0.163144740
# Lateral root cap 1           0.036938055 0.228698023
# Protophloem                  0.124940346 0.077368780
# Pericycle                    0.136775795 0.062713020
# Stele                        0.126276606 0.065326062
# Root cortex                  0.041328625 0.129175187
# Non-hair root epidermal cell 0.095733512 0.053965008
# Columella                    0.100028634 0.032833447
# Root endodermis              0.052591391 0.040331743
# Quiescent center             0.041710413 0.037037037
# Trichoblast                  0.030733989 0.037718700
# Columella root cap cell      0.039801470 0.016700750
# Xylem                        0.029397728 0.027493751
# Lateral root cap 2           0.006108619 0.011701886
# Phloem                       0.008685692 0.008520791
# Unknown 2                    0.004772358 0.007271075
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1] <- 'celltype'
colnames(Cellratio)[2] <- 'sample'
## save the cell ratio result
write.csv(Cellratio, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DA/cell_ratio.csv")
library(ggplot2)
p <- ggplot(data = Cellratio,aes(x =celltype, y= Freq, fill = sample)) + 
  geom_bar(stat = "identity",position = "dodge", width = 0.8)+ 
  theme_classic() +
  labs(x='Cell Type',y = 'Ratio')+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme_bar()
  #coord_flip()+
  #theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p <- p +  scale_fill_npg()
ggsave('cellratio_barplot_celltype.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DA/",
       width = 12, height = 8)


### box plot
table(CK.harmony$Condition)#查看各组细胞数
prop.table(table(Idents(CK.harmony)))
cell_count <- table(Idents(CK.harmony), CK.harmony$Condition)#各组不同细胞群细胞数
write.csv(cell_count, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DA/cell_count.csv")
Cellratio <- prop.table(table(Idents(CK.harmony), CK.harmony$Condition), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###添加分组信息
sample <- c("CKN1","CKG1")
group <- c("CKN","CKG")
samples <- data.frame( group)#创建数据框

rownames(samples)=samples$group
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列

###作图展示
pplist = list()
sce_groups = c(CK.harmony$celltype %>% unique() %>% as.vector())
library(ggplot2)
library(dplyr)
library(ggpubr)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('group',group_)))#选择一组数据
  colnames(cellper_) = c('group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("CKN", "CKG") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

library(cowplot)
p <- plot_grid(pplist[['Stele']],
          pplist[['Root cortex']],
          pplist[['Protophloem']],
          pplist[['Columella']],
          pplist[['Trichoblast']],
          pplist[['Unknown 1']],
          pplist[['Lateral root cap 1']],
          pplist[['Xylem']],
          pplist[['Root endodermis']],
          pplist[['Quiescent center']],
          pplist[['Non-hair root epidermal cell']],
          pplist[['Pericycle']],
          pplist[['Lateral root cap 2']],
          pplist[['Columella root cap cell']],
          pplist[['Phloem']],
          pplist[['Unknown 2']])

ggsave('cellratio_all.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DA/",
       width = 12, height = 12)


if(T){
  library(ggrepel)
  require(qdapTools)
  require(REdaS)
  # 
  metadata <- CK.harmony@meta.data
  # Keep only cells from tissues that are not brain or pleura 
  metadata <- metadata[-which(metadata$Condition=="CKG" | metadata$Condition=="CKN"),]
  # Convert to factor with logical order 
  metadata$sample <- factor(metadata$Condition, levels = c("CKG", "CKN"))
  # Create table and keep selected cell types 
  meta.temp <- metadata[,c("celltype", "sample")]
  # Loop over treatment response categories 
  # Create list to store frequency tables 
  prop.table.error <- list()
  for(i in 1:length(unique(meta.temp$analysis))){
    vec.temp <- meta.temp[meta.temp$analysis==unique(meta.temp$analysis)[i],"immuSub"]
    # Convert to counts and calculate 95% CI 
    # Store in list 
    table.temp <- freqCI(vec.temp, level = c(.95))
    prop.table.error[[i]] <- print(table.temp, percent = TRUE, digits = 3)
    # 
  }
  # Name list 
  names(prop.table.error) <- unique(meta.temp$analysis)
  # Convert to data frame 
  tab.1 <- as.data.frame.array(do.call(rbind, prop.table.error))
  # Add analysis column 
  b <- c()
  a <- c()
  for(i in names(prop.table.error)){
    a <- rep(i,nrow(prop.table.error[[1]]))
    b <- c(b,a)
  }
  tab.1$analysis <- b
  # Add common cell names 
  aa <- gsub(x = row.names(tab.1), ".1", "")
  aa <- gsub(x = aa, ".2", "")
  tab.1$cell <- aa
  # 
  # Resort factor analysis 
  tab.1$analysis <- factor(tab.1$analysis, levels = c("naive", "grouped_pr", "grouped_pd"))
  # Rename percentile columns 
  colnames(tab.1)[1] <- "lower"
  colnames(tab.1)[3] <- "upper"
  # 
  p<- ggplot(tab.1, aes(x=analysis, y=Estimate, group=cell)) +
    geom_line(aes(color=cell))+
    geom_point(aes(color=cell)) + facet_grid(cols =  vars(cell)) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=0.5), legend.position="bottom") + 
    xlab("") + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
  p1<- ggplot(tab.1, aes(x=analysis, y=Estimate, group=cell)) +
    geom_bar(stat = "identity", aes(fill=cell)) + facet_grid(cols =  vars(cell)) + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=0.5), legend.position= "none") + 
    xlab("") + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05)) 
  p1
}
library(gplots)
tab.1=table(CK.harmony$Condition,CK.harmony$celltype) 
head(tab.1)
balloonplot(tab.1)





### 带比例
PropPlot <- function(object, groupBy){
  # (1)获取绘图数据
  plot_data = object@meta.data %>% 
    dplyr::select(orig.ident, {{groupBy}}) %>% 
    dplyr::rename(group = as.name(groupBy))
  
  # (2)绘图
  figure = ggbarstats(data = plot_data, 
                      x = group, y = orig.ident,
                      package = 'ggsci',
                      palette = 'category20c_d3',
                      results.subtitle = FALSE,
                      bf.message = FALSE,
                      proportion.test = FALSE,
                      label.args = list(size = 2, 
                                        fill = 'white', 
                                        alpha = 0.85,
                                        family = 'Arial',
                                        fontface = 'bold'),
                      perc.k = 2,
                      title = '',
                      xlab = '',
                      legend.title = 'Seurat Cluster',
                      ggtheme = ggpubr::theme_pubclean()) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = 'black', lineend = 'round'),
          legend.position = 'right',
          axis.text.x = element_text(size = 15, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 15, color = 'black', family = 'Arial'),
          legend.text = element_text(family = 'Arial', size = 10, color = 'black'),
          legend.title = element_text(family = 'Arial', size = 13, color = 'black')) 
  
  # (3)去除柱子下面的样本量标识：
  gginnards::delete_layers(x = figure, match_type = 'GeomText')
}


PropPlot(object = CK.harmony, groupBy = 'celltype')











patient_summary <- metadata[,c("cancer","patient")]
patient_summary <- patient_summary[!duplicated(patient_summary),]
table(patient_summary$cancer)

metadata$sample <- paste0(metadata$Condition)
sample_summary <- metadata[,c("sample")]
sample_summary <- sample_summary[!duplicated(sample_summary),]
table(sample_summary$cancer)

metadata_filt <- metadata[metadata$celltype,]
#metadata_filt <- metadata_filt[metadata_filt$tissue!='L',]
metadata_filt$sample <- factor(metadata_filt$sample, levels=c("CKN","CKG"))
#metadata_filt$cancer <- factor(metadata_filt$cancer, levels=c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L"))

#metadata_filt$celltype <- factor(metadata_filt$Cluster, levels = c("Mast","pDC",  "cDC", "Mo/Mq"))

## cluster distribution
ggplot(metadata_filt, aes(factor(sample)))+ 
  geom_bar(aes(fill = celltype), position = "fill")+ 
  xlab("")+ylab("Proportion")+facet_wrap(~sample,ncol=5,scales='free_x')+
  theme(legend.title=element_blank(),strip.background.x = element_blank())+ 
  #scale_fill_manual(values=c("#6495ED","#FFA500","#FF4500","#999999"))+
  theme_classic2()


## MajorCluster comparison
## only include N and T
metadata_filt_N_T <- metadata_filt[metadata_filt$tissue %in% c("N","T"),]

## Mast
metadata_filt_N_T <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("BC","CRC","ESCA","Gastric","HCC","Lung","OV","PACA","RC","THCA","UCEC"),]
metadata_filt_N_T$tissue <- factor(metadata_filt_N_T$tissue)
cell_num <- ddply(metadata_filt_N_T, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_Mast <- metadata_filt_N_T[metadata_filt_N_T$Cluster %in% c("Mast"),]
cell_num_cluster <- ddply(metadata_filt_Mast, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

## paired wilcox test
library(ggpubr)
cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cellper, x = "variable", y = "value",
              color='variable',palette = "jco", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = TRUE)+facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of mast cells")+ xlab("")
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "less"))
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "greater"))

## cDC
metadata_filt_N_T <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("BC","CRC","ESCA","Gastric","HCC","Lung","OV","PACA","RC","THCA","UCEC"),]
metadata_filt_N_T$tissue <- factor(metadata_filt_N_T$tissue)
cell_num <- ddply(metadata_filt_N_T, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_cDC <- metadata_filt_N_T[metadata_filt_N_T$Cluster %in% c("cDC"),]
cell_num_cluster <- ddply(metadata_filt_cDC, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=100,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

## paired wilcox test
library(ggpubr)
cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
              color='variable',palette = "jco", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = TRUE)+facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of cDCs")+ xlab("")
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "less"))
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "greater"))



## DC subset comparison
## DC subset cluster distribution
metadata_filt_DC <- metadata_filt[metadata_filt$Cluster %in% c("cDC"),]
metadata_filt_DC$DC_type <- unlist(lapply(strsplit(metadata_filt_DC$MajorCluster,"_"),function(x){x[2]}))

ggplot(metadata_filt_DC, aes(factor(tissue)))+ geom_bar(aes(fill = DC_type), position = "fill")+ xlab("")+ylab("Proportion")+facet_wrap(~cancer,ncol=5,scales='free_x')+theme(legend.title=element_blank(),strip.background.x = element_blank())+scale_fill_brewer(palette="Set1")

metadata_filt_DC$DC_type <- factor(metadata_filt_DC$DC_type, levels=c("pDC",'cDC1','cDC2','cDC3'))
current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt_DC$cancer <- plyr::mapvalues(x = metadata_filt_DC$cancer, from = current.cluster.ids, to = new.cluster.ids)
ggplot(metadata_filt_DC, aes(factor(tissue)))+ geom_bar(aes(fill = DC_type), position = "fill")+ xlab("")+ylab("Proportion")+facet_wrap(~cancer,ncol=5,scales='free_x')+theme(legend.title=element_blank(),strip.background.x = element_blank())+ scale_fill_manual(values=c("#26933B","#F3951B","#E71638"))+theme_classic2()

## LAMP3 DC
metadata_filt_DC_N_T <- metadata_filt_DC[metadata_filt_DC$tissue %in% c("N","T"),]
current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt_DC_N_T$cancer <- plyr::mapvalues(x = metadata_filt_DC_N_T$cancer, from = current.cluster.ids, to = new.cluster.ids)

metadata_filt_DC_N_T <- metadata_filt_DC_N_T[metadata_filt_DC_N_T$cancer %in% c("BRCA","CRC","ESCA","STAD","HCC","LUNG","OV-FTC","PAAD","KIDNEY","THCA","UCEC"),]
cell_num <- ddply(metadata_filt_DC_N_T, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_DC_N_T_LAMP3 <- metadata_filt_DC_N_T[metadata_filt_DC_N_T$DC_type %in% c("cDC3"),]
cell_num_cluster <- ddply(metadata_filt_DC_N_T_LAMP3, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
              color='variable',palette = "jco", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = TRUE) +facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of LAMP3+ cDCs in cDCs")+ xlab("")

p + stat_compare_means(paired = TRUE,method.args = list(alternative = "less"))


### plot macrophage subset distribution for lung cancer and UCEC
## macrophage from Lung cancer and UCEC
metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("Lung","UCEC","THCA"),]
metadata_filt_N_T_cancer$tissue <- factor(metadata_filt_N_T_cancer$tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_SPP1 <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$MajorCluster %in% c("M10_Macro_SPP1","M09_Macro_SPP1"),]
cell_num_cluster <- ddply(metadata_filt_SPP1, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

## paired wilcox test
library(ggpubr)
cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
              color='variable',palette = "jco", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = TRUE) +facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of Macro_SPP1")+ xlab("")
p + stat_compare_means(paired = TRUE)



########### only consider composition of myeloid cell in tumor tissues from different tumor types
setwd("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate")
metadata <- read.csv("metadata.csv")
metadata_filt <- metadata[metadata$Cluster!='Myeloid',]
metadata_filt <- metadata_filt[metadata_filt$tissue=='T',]

current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt$cancer <- plyr::mapvalues(x = metadata_filt$cancer, from = current.cluster.ids, to = new.cluster.ids)


cell_num <- ddply(metadata_filt, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(metadata_filt, .(patient,cancer,Cluster), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer+Total~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer","Total"))

cell_num_cluster_summary_formated_melt <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$Total >=100,]
cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)

# Color panel -----------------
c54 <- c('dodgerblue2','green4','#E31A1C','#6A3D9A','#FF7F00',
         '#FB9A99','#CAB2D6','khaki2','deeppink1','blue1',      
         'steelblue4','green1','yellow4','yellow3','forestgreen',
         'red2','orange','cornflowerblue', 'magenta','darkolivegreen4',
         'indianred1','tan4','darkblue','mediumorchid1','firebrick4',
         'yellowgreen','lightsalmon','tan3','tan1','darkgray',
         'wheat4','#DDAD4B','chartreuse','seagreen1','moccasin',
         'mediumvioletred','seagreen','cadetblue1','darkolivegreen1','tan2',
         'tomato3','#7CE3D8','gainsboro','gold1','skyblue2',
         'palegreen2','#FDBF6F','gray70','darkorange4','orchid1',
         'darkturquoise','maroon','brown','black')


c54 <- c("BRCA" = 'dodgerblue2',"ESCA"='green4',"STAD"='#E31A1C',"LUNG"='#6A3D9A',"OV-FTC"='#FF7F00',
         "CRC"='#FB9A99',"PAAD"='#CAB2D6',"KIDNEY"='khaki2',"THCA"='deeppink1',"UCEC"='blue1',      
         "HCC"='steelblue4',"NPC"='green1',"MEL"='yellow4',"MYE"='yellow3',"LYM"='forestgreen')

type <- "Mast"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p1 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 70, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")

type <- "pDC"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p2 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 50, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")
type <- "cDC"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p3 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 50, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")


type <- "Mo/Mq"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p4 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")

library(scater)
pdf("/data2/csj/Pan_Myeloid/A20191105/P1_fraction.pdf",width=4.78,height=6.23)
multiplot(p1,p2,cols=1)
dev.off()
pdf("/data2/csj/Pan_Myeloid/A20191105/P2_fraction.pdf",width=4.78,height=6.23)
multiplot(p3,p4,cols=1)
dev.off()

multiplot(p4,p2,p3,p1,cols=1)


#### only consider DC cells
## DC subset comparison
## DC subset cluster distribution
metadata_filt_DC <- metadata_filt[metadata_filt$Cluster %in% c("cDC"),]
metadata_filt_DC$DC_type <- unlist(lapply(strsplit(metadata_filt_DC$MajorCluster,"_"),function(x){x[2]}))

cell_num <- ddply(metadata_filt_DC, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(metadata_filt_DC, .(patient,cancer,DC_type), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=50,]

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+Total+cancer~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","Total","cancer"))

cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)

type <- "cDC1"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p1 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")
type <- "cDC2"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p2 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")

type <- "cDC3"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p3 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
                  legend.position = "null",
                  plot.title = element_text(
                    size = 16,
                    face = "bold",
                    hjust = 0.5
                  ),
                  text = element_text(size = 10),
                  plot.margin = unit(c(1, 1, 1, 1), "char"),
                  axis.text.x = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                  ),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 15)
                )+xlab("")
library(scater)
multiplot(p1,p2,p3,cols=1)
pdf("/data2/csj/Pan_Myeloid/A20191105/cDC_fraction.pdf",width=4.78,height=9)
multiplot(p1,p2,p3,cols=1)
dev.off()

########### only consider composition of myeloid cell in tumor tissues from different tumor types and calcute anogigenesis and phagocytosis ratio

setwd("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate")
metadata <- read.csv("metadata.csv")
metadata_filt <- metadata[metadata$Cluster!='Myeloid',]
metadata_filt <- metadata_filt[metadata_filt$tissue=='T',]

current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt$cancer <- plyr::mapvalues(x = metadata_filt$cancer, from = current.cluster.ids, to = new.cluster.ids)

metadata_filt$class <- "other"
metadata_filt[metadata_filt$cancer=='BRCA' & metadata_filt$MajorCluster =='M09_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='BRCA' & metadata_filt$MajorCluster =='M08_Macro_CLEC5A',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='BRCA' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='OV-FTC' & metadata_filt$MajorCluster =='M08_Macro_MARCO',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='OV-FTC' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='PAAD' & metadata_filt$MajorCluster =='M08_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='PAAD' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='KIDNEY' & metadata_filt$MajorCluster =='M08_Macro_FN1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='KIDNEY' & metadata_filt$MajorCluster =='M11_Macro_SLC40A1',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='STAD' & metadata_filt$MajorCluster =='M08_Macro_IL1B',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='STAD' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='THCA' & metadata_filt$MajorCluster =='M09_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='THCA' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='UCEC' & metadata_filt$MajorCluster =='M10_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='UCEC' & metadata_filt$MajorCluster =='M11_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='CRC' & metadata_filt$MajorCluster =='M11_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='CRC' & metadata_filt$MajorCluster =='M12_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='ESCA' & metadata_filt$MajorCluster =='M09_Macro_IDO1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='ESCA' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='HCC' & metadata_filt$MajorCluster =='M08_Macro_FCN1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='HCC' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='LUNG' & metadata_filt$MajorCluster =='M09_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='LUNG' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='LYM' & metadata_filt$MajorCluster =='M06_Macro_IL1B',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='LYM' & metadata_filt$MajorCluster =='M07_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='MEL' & metadata_filt$MajorCluster =='M07_Macro_VCAN',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='MEL' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='MYE' & metadata_filt$MajorCluster =='M08_Macro_IL1B',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='MYE' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'


metadata_filt[metadata_filt$cancer=='NPC' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'anogigenesis'

cell_num_cluster <- ddply(metadata_filt, .(patient,cancer,class), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_formated <- dcast(cell_num_cluster,patient+cancer~cluster,value.var=c('number'))
cell_num_cluster_formated[is.na(cell_num_cluster_formated)] <- 0
cell_num_cluster_formated$Total <- cell_num_cluster_formated$anogigenesis + cell_num_cluster_formated$phagocytosis + cell_num_cluster_formated$other
cell_num_cluster_formated <- cell_num_cluster_formated[cell_num_cluster_formated$Total >=100,]
cell_num_cluster_formated$Ratio <- round(log2(cell_num_cluster_formated$anogigenesis/cell_num_cluster_formated$phagocytosis),2)


median_table <- ddply(cell_num_cluster_formated,.(cancer), function(x){median(x$Ratio)})
median_table_sorted <- median_table[order(median_table$V1),]
cell_num_cluster_formated$cancer <- factor(cell_num_cluster_formated$cancer, levels=median_table_sorted$cancer)


ggboxplot(cell_num_cluster_formated, x = "cancer", y = "Ratio",
          color = "cancer", palette =c54,
          add = "jitter") + stat_compare_means(label.y = 6, label.x = 2) + ggtitle(paste0("Anogigenesis/phagocytosis macrophages"))+ theme(
            legend.position = "null",
            plot.title = element_text(
              size = 16,
              face = "bold",
              hjust = 0.5
            ),
            text = element_text(size = 10),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            axis.text.x = element_text(
              size = 12,
              angle = 45,
              hjust = 1
            ),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 15)
          )

## calcutate the anogigenesis fraction in all myeloid cell
cell_num <- ddply(metadata_filt, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(metadata_filt, .(patient,cancer,class), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+Total+cancer~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","Total","cancer"))

cell_num_cluster_summary_formated_melt <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$Total >=100,]
cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)


type <- "anogigenesis"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

pdf("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/angiogenesis-associated.pdf",width=6.19,height=3.51)
ggboxplot(df, x = "cancer", y = "Proportion",
          color = "cancer", palette =c54,
          add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of angiogenesis-associated macrophages"))+ theme(
            legend.position = "null",
            plot.title = element_text(
              size = 16,
              face = "bold",
              hjust = 0.5
            ),
            text = element_text(size = 10),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            axis.text.x = element_text(
              size = 12,
              angle = 45,
              hjust = 1
            ),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 15)
          )
dev.off()	  
## whole proportion

colorder = c("Ctrl-1","Ctrl-2","Ctrl-3",
             "C2-1","C2-2","C2-3",
             "C5-1","C5-2", "C5-3","C5-4")

x <- table(CK.harmony$celltype,CK.harmony$Condition)
x <- x[, colorder]
x3= t(t(x)/rowSums(t(x)))

x4 = as.data.frame(as.table(t(x3)))
colnames(x4) = c("sample","celltype","Freq")
x4$group = x4$sample %>% str_replace("-.*","")
x4$group = factor(x4$group, levels = c("Ctrl","C2","C5"))


write.csv(x4, "Fig1e.proportion_each_cluster.csv")

dose<-read.csv( "Fig1e.proportion_each_cluster.csv")
which(dose$celltype=="RBC")
which(dose$celltype=="Doublet")
dose<-dose[-c(121:140),]
top<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottom<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
dose_Ctrl<-dose[which(dose$group=="Ctrl"),]
dose_C2<-dose[which(dose$group=="C2"),]
dose_C5<-dose[which(dose$group=="C5"),]

ggplot(data=dose_Ctrl,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_Ctrl,aes(celltype,Freq),size=3,pch=19) 


ggplot(data=dose_C2,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_C2,aes(celltype,Freq),size=3,pch=19) 


ggplot(data=dose_C5,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_C5,aes(celltype,Freq),size=3,pch=19) 
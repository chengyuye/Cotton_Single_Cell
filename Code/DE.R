### import libraries
#???ذ???????
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(MAST)

#read in data
CK.harmony <- readRDS("C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds")

### seperate group
CK.harmony@meta.data$celltype_aggregate = paste(CK.harmony@meta.data$celltype, 
                                                CK.harmony@meta.data$Condition,sep = "_") # user adaptation required on own dataset
DimPlot(CK.harmony, group.by = "celltype_aggregate")


#DEGs ---LRC1
diff_LRC1 <- FindMarkers(CK.harmony, 
                         group.by = "celltype_aggregate",
                         ident.1 ="Lateral root cap 1_CKG",
                         ident.2="Lateral root cap 1_CKN")
#DEGs ---LRC2
diff_LRC2 <- FindMarkers(CK.harmony, 
                         group.by = "celltype_aggregate",
                         ident.1 ="Lateral root cap 2_CKG",
                         ident.2="Lateral root cap 2_CKN")

#DEGs ---CRC
diff_CRC <- FindMarkers(CK.harmony, 
                         group.by = "celltype_aggregate",
                         ident.1 ="Columella root cap cell_CKG",
                         ident.2="Columella root cap cell_CKN")

#diff_LRC1 <- as.data.frame(diff_LRC1)
diff_LRC1$cluster <- 'LRC 1'
diff_LRC1$geneID <- rownames(diff_LRC1)
rownames(diff_LRC1) <- c(1:dim(diff_LRC1)[1])
diff_LRC1 %>% head()

diff_LRC2$cluster <- 'LRC 2'
diff_LRC2$geneID <- rownames(diff_LRC2)
rownames(diff_LRC2) <- c(1:dim(diff_LRC2)[1])
diff_LRC2 %>% head()

diff_CRC$cluster <- 'CRC'
diff_CRC$geneID <- rownames(diff_CRC)
rownames(diff_CRC) <- c(1:dim(diff_CRC)[1])
diff_CRC %>% head()

write.csv(diff_LRC1, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/Diff_LRC1.csv")
write.csv(diff_LRC2, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/Diff_LRC2.csv")
write.csv(diff_CRC, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/Diff_CRC1.csv")

## merge them together
diff_LRC <- rbind(diff_LRC1, diff_LRC2)
diff_LRC <- rbind(diff_LRC,diff_CRC)
head(diff_LRC)

library(Seurat)
library(patchwork)
library(clusterProfiler)
library(tidyverse)

### read in possypol related gene list
gossipol <- read.table("C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/gossypol_gene.txt",col.names = 'Gene_ID')

##finder overlapped gene for LRC 1 2 and CRC
overlap_lrc1 <- intersect(gossipol$Gene_ID, diff_LRC1$geneID)
overlap_lrc2 <- intersect(gossipol$Gene_ID, diff_LRC2$geneID)
overlap_crc <- intersect(gossipol$Gene_ID, diff_CRC$geneID)

#diff_LRC1$names <- rownames(diff_LRC1)
#sig_dge.all <- subset(object.markers, p_val_adj<0.05&abs(avg_log2FC)>0.15) #???в???????
#View(sig_dge.all)
library(dplyr)
diff_LRC1 <- diff_LRC1 %>%
  mutate(Difference = pct.1 - pct.2)
library(ggplot2)
library(ggrepel)

# ggplot(diff_LRC1, aes(x=Difference, y=avg_log2FC)) + 
#   geom_point(size=0.5,aes(color=group)) + 
#   scale_color_manual(values=c('blue','grey','red'))+
#   geom_label_repel(data=subset(object.markers, group !='no'), aes(label=names), segment.size = 0.25, size=2.5)+
#   geom_vline(xintercept = 0.0,linetype=2)+
#   geom_hline(yintercept = 0,linetype=2)+
#   theme_classic()
# ggsave("TopMarkerVol2.pdf", height=8, width=8)
# ggsave("TopMarkerVol1.pdf", height=8, width=8)



ggplot(diff_LRC1,aes(x=avg_log2FC,y= -log10(p_val),fill = threshold))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_size(limits  = c(2, 16))+ #控制最大气泡和最小气泡，调节气泡相对大小
  scale_fill_manual(values=c("#fe0000","#13fc00","#bdbdbd"))+#确定点的颜色
  geom_text_repel(
    data = subset(diff_LRC1, diff_LRC1$geneID == overlap_lrc1),
    aes(label = geneID),
    size = 3.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  ylab('-log10 (Pvalue)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.25,0.25),lty=2,col="black",lwd=0.5) +#添加横线|logFoldChange|>0.25
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +#添加竖线padj<0.05
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(fill=guide_legend(override.aes = list(size=5)))+ #图例圈圈大小
  theme(axis.title.x = element_text(size = 10, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 10,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 7, 
                                   face = "bold")
  )

ggplot(data = diff_LRC1, aes(x = avg_log2FC, y = -log10(p_val), colour=threshold,label = geneID)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-11.5, 11.5)) +
  geom_vline(xintercept=c(-3,3),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.001),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential metabolites") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
  geom_text_repel(
    data = subset(diff_LRC1, diff_LRC1$geneID == overlap_lrc1),
    aes(label = geneID),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )


### enhancedvolcano lrc1
library(EnhancedVolcano) 
#确定是上调还是下调，用于给图中点上色
diff_LRC1$threshold = factor(ifelse(diff_LRC1$p_val < 0.05 & abs(diff_LRC1$avg_log2FC) >= 0.25, 
                                    ifelse(diff_LRC1$avg_log2FC >= 0.25 ,'Up','Down'),'NoSignifi'),
                             levels=c('Up','Down','NoSignifi'))

diff_LRC1_up <- subset(diff_LRC1, diff_LRC1$avg_log2FC >= 0.25 & diff_LRC1$p_val <= 0.05)
overlap_lrc1 <- intersect(gossipol$Gene_ID, diff_LRC1_up$geneID)

lrc1 <- EnhancedVolcano(diff_LRC1,
                lab = diff_LRC1$geneID,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Lateral root cap 1',
                FCcutoff = 0.25,
                subtitle = NULL,
                pointSize = 1.5,
                labSize  = 3.0,
                selectLab = overlap_lrc1,
                #boxedLabels  = TRUE,
                colAlpha = 4/5,
                legendLabels =c('NS','Log (base 2) fold-change','P value','P value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                maxoverlapsConnectors = Inf)
ggsave('LRC1_DEG_volcano_highlight.pdf', lrc1, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 12, height = 7)

ggsave('LRC1_DEG_volcano.pdf', lrc1, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 12, height = 6)

### LRC2
lrc2 <- EnhancedVolcano(diff_LRC2,
                        lab = diff_LRC2$geneID,
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = 'Lateral root cap 2',
                        subtitle = NULL,
                        pointSize = 1.5,
                        labSize  = 3.0,
                        # ?ӿ?
                        #boxedLabels  = TRUE,
                        colAlpha = 4/5,
                        legendLabels =c('NS','Log (base 2) fold-change','P value','P value & Log (base 2) fold-change'),
                        legendPosition = 'right',
                        legendLabSize = 10,
                        legendIconSize = 3.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black')

ggsave('LRC2_DEG_volcano.pdf', lrc2, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 10, height = 6)


### CRC
crc <- EnhancedVolcano(diff_CRC,
                        lab = diff_CRC$geneID,
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = 'Columella root cap',
                        subtitle = NULL,
                        pointSize = 1.5,
                        labSize  = 3.0,
                        # ?ӿ?
                        #boxedLabels  = TRUE,
                        colAlpha = 4/5,
                        legendLabels =c('NS','Log (base 2) fold-change','P value','P value & Log (base 2) fold-change'),
                        legendPosition = 'right',
                        legendLabSize = 10,
                        legendIconSize = 3.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        colConnectors = 'black')

ggsave('CRC_DEG_volcano.pdf', crc, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 10, height = 6)

#### 
#?????????Ա?ǩ??
diff_LRC$label <- ifelse(diff_LRC$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
head(diff_LRC)

#??ȡÿ??cluster?б???????????????10????????
top10sig1 <- filter(diff_LRC,cluster=="LRC 1") %>% 
  distinct(geneID,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))
head(top10sig1)


top10sig2 <- filter(diff_LRC,cluster=="LRC 2") %>% 
  distinct(geneID,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))
head(top10sig2)

top10sig3 <- filter(diff_LRC,cluster=="CRC") %>% 
  distinct(geneID,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))
head(top10sig3)

#????ȡ????cluster??Top10?????????ϲ???
top10sig <- rbind(top10sig1,top10sig2,top10sig3)

#????һ?У???Top10?Ĳ???????????Ϊ2???????ı???Ϊ1??
diff_LRC$size <- case_when(!(diff_LRC$geneID %in% top10sig$geneID)~ 1,
                           diff_LRC$geneID %in% top10sig$geneID ~ 2)

#??ȡ??Top10?Ļ?????????
dt <- filter(diff_LRC,size==1)
head(dt)

#????ÿ??Cluster Top10??????????ɢ????ɽͼ??
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)
p

#????ÿ??Cluster Top10????ɢ??(??ɢ???ʵ??Ŵ?ǿ??????
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

#????ͼp??log2FC????ȷ???????????ȣ?
dfbar <-data.frame(x=c("CRC","LRC 1", "LRC 2"),
                  y=c(2.5,2.5,3))
dfbar1 <-data.frame(x=c("CRC","LRC 1", "LRC 2"),
                   y=c(-2.8,-2.2,-2.4))
#???Ʊ???????
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

#??ɢ????ɽͼ???ӵ????????ϣ?
p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p2


#????X????clusterɫ????ǩ??
dfcol<-data.frame(x=c("CRC","LRC 1", "LRC 2"),
                  y=0,
                  label=c("CRC","LRC 1", "LRC 2"))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3



#??ÿ??Cluster????????ǰTop10???????ϱ?ǩ??
options(ggrepel.max.overlaps = Inf)
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=avg_log2FC,label=geneID),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4

#ɢ????ɫ??????
p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5

#?޸?X/Y????????????cluster???֣?
p6 <- p5+
  labs(x="Cluster",y="Average log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =6,
            color ="white")

p6

#?Զ????????��???
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7

ggsave('LRC_DEG_volcano_all.pdf', p7, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 8, height = 10)

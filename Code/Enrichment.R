### GO/KEGG analysis
library(clusterProfiler)

## read in gene lists
# #CRC
# CRC_gene <- read.table('C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/DIFF_CRC1_idlist.txt')
# diff_CRC <- CRC_gene$V2
## 获取上下调基因 ---LRC 1
gene_up <-  diff_CRC %>% filter(avg_log2FC > 0) %>% 
  select(geneID)
gene_up <- gene_up$geneID

gene_down <-  diff_CRC %>% filter(avg_log2FC < 0) %>% 
  select(geneID)
gene_down <- gene_down$geneID

### modify the gene name ###
### delete the LOC ###
gene_up <- gsub('[LOC]', '', gene_up)
gene_down <- gsub('[LOC]', '', gene_down)

## KEGG ----- up-regulated genes
gene_up <- unique(gene_up)
kk.up <- enrichKEGG(gene = gene_up,
                    organism = "ghi",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
write.csv(kk.up, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/CRC_KEGG_up_regulated.csv')
p <- clusterProfiler::dotplot(kk.up,showCategory = 10)
ggsave('CRC_KEGG_Up_dotplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 10, height = 5)
p <- barplot(kk.up,,showCategory = 10)
ggsave('CRC_KEGG_Up_barplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 10, height = 5)

## KEGG ----- down-regulated genes
gene_down <- unique(gene_down)
kk.down <- enrichKEGG(gene = gene_down,
                    organism = "ghi",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
write.csv(kk.down, file = 'C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/CRC_KEGG_down_regulated.csv')
p <- clusterProfiler::dotplot(kk.down,showCategory = 10)
ggsave('CRC_KEGG_down_dotplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 10, height = 5)
p <- barplot(kk.down,,showCategory = 10)
ggsave('CRC_KEGG_down_barplot.pdf', p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
       width = 10, height = 5)



### saving as a function 
LRC_KEGG <- function(gene_list, name){
  ## 获取上下调基因 ---
  gene_up <-  gene_list %>% filter(avg_log2FC > 0) %>% 
    select(geneID)
  gene_up <- gene_up$geneID
  
  gene_down <-  gene_list %>% filter(avg_log2FC < 0) %>% 
    select(geneID)
  gene_down <- gene_down$geneID
  
  ### modify the gene name ###
  ### delete the LOC ###
  gene_up <- gsub('[LOC]', '', gene_up)
  gene_down <- gsub('[LOC]', '', gene_down)
  
  ## KEGG ----- up-regulated genes
  gene_up <- unique(gene_up)
  kk.up <- enrichKEGG(gene = gene_up,
                      organism = "ghi",
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  write.csv(kk.up, file = paste0('C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/',name, '_KEGG_up_regulated.csv'))
  p <- dotplot(kk.up,showCategory = 10)
  ggsave(paste0(name,'_KEGG_Up_dotplot.pdf'), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
         width = 10, height = 5)
  p <- barplot(kk.up,,showCategory = 10)
  ggsave(paste0(name,'_KEGG_Up_barplot.pdf'), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
         width = 10, height = 5)
  
  ## KEGG ----- down-regulated genes
  gene_down <- unique(gene_down)
  kk.down <- enrichKEGG(gene = gene_down,
                        organism = "ghi",
                        pvalueCutoff = 0.9,
                        qvalueCutoff = 0.9)
  write.csv(kk.down, file = paste0('C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/',name, '_KEGG_down_regulated.csv'))
  p <- clusterProfiler::dotplot(kk.down,showCategory = 10)
  ggsave(paste0(name,'_KEGG_down_dotplot.pdf'), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
         width = 10, height = 5)
  p <- barplot(kk.down,,showCategory = 10)
  ggsave(paste0(name,'_KEGG_down_barplot.pdf'), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DE/",
         width = 10, height = 5)
}
## CRC
LRC_KEGG(diff_CRC,'CRC')


### LRC 1
LRC_KEGG(diff_LRC1,'LRC1')


### LRC 2
LRC_KEGG(diff_LRC2,'LRC2')

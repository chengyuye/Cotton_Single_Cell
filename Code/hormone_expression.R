# read in data 
LRC <- readRDS("C:/BGI/mianhua/CK_harmony_updated_test/CK_LRC_annotated.rds")
CK.harmony <- 
Idents(LRC) <- LRC$celltype

### pathway expression function define
path_expr <- function(object, gene_list, path_name){
  # calculate the module score first
  object <- AddModuleScore(object, features = list(gene_list), name = path_name)
  ## feature plot to see the expression of the whole pathway
  #1. not split by condition
  p <- FeaturePlot(object,
                   features = paste0(path_name, '1'), label = TRUE, repel = TRUE,
                   pt.size = 0.5,label.size = 3
  )
  
  # save the plot
  ggsave(paste0(path_name,"_expression.pdf"), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/hormone_pathway/Root cap/",
         width = 6, height = 5) 
  #2. split by condition
  p <- FeaturePlot(object,
                   features = paste0(path_name, '1'), 
                   label = TRUE, repel = TRUE, 
                   #min.cutoff = "q10", max.cutoff = "q90",
                   split.by = 'Condition',
                   pt.size = 0.5,label.size = 3)
  #&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  # save the plot
  ggsave(paste0(path_name,"_expression_by_condition.pdf"), p, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/hormone_pathway/Root cap/",
         width = 12, height = 5)
}

### pathway expression function define using scillus
path_expr_scillus <- function(object, gene_list, path_name){
  path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/hormone_pathway/"
  pdf(paste0(path, path_name,"_expression_heatmap.pdf"),
      width = 10, height = 6)
  p <- plot_heatmap(dataset = object, 
                    markers = gene_list,
                    sort_var = c("celltype","Sample"),
                    anno_var = c("celltype","Sample"),
                    anno_colors = list("Set2",                                             # RColorBrewer palette
                                       c("#0570b0", "#d7301f"), # color vector
                                       "Reds"),
                    hm_colors = c("black","white","purple"))
  dev.off()
  
}

##scale all genes
CK.harmony <- ScaleData(CK.harmony, features = rownames(CK.harmony))
DoHeatmap(CK.harmony, features = BBB,assay = 'RNA',group.by = 'celltype')


#BAK1-BRI1-Brassinosteroids pathway  -----LOC107910843, LOC121210267 absent
BBB <- c('LOC107890835','LOC107899559',
         'LOC107922051','LOC107950600')

path_expr(CK.harmony,BBB,'BBB')
path_expr(LRC,BBB,'BBB')
##sciluus
path_expr_scillus(CK.harmony, BBB, 'BBB')


###
exp <- GetAssayData(CK.harmony, slot = "counts")
exp <- log10(exp + 1)
head(CK.harmony$celltype)
new_celltype <- sort(CK.harmony$celltype)
head(new_celltype)
cs_data <- as.matrix(exp[AUX, names(new_celltype)])


###PYR-PYL-Abscisic acid  ------LOC107903191, LOC107903272, LOC107929906, LOC107953883, LOC121204798 absent
PPA <- c('LOC107890765','LOC107890941','LOC107900548',
         'LOC107900942','LOC107902420',
         'LOC107904257','LOC107904467',
         'LOC107909797','LOC107911088','LOC107912184',
         'LOC107924284','LOC107924666','LOC107924725',
         'LOC107924990','LOC107929605',
         'LOC107930169','LOC107931756','LOC107932889',
         'LOC107933250','LOC107937038','LOC107937632',
         'LOC107938971','LOC107939394','LOC107946749',
         'LOC107946830','LOC107947093','LOC107947117',
         'LOC107947258','LOC107958082',
         'LOC107961042','LOC107961366','LOC107962125'
         )

path_expr(CK.harmony,PPA,'PPA')
path_expr(LRC,PPA,'PPA')
DoHeatmap(CK.harmony, features = PPA,assay = 'RNA',group.by = 'celltype')
##sciluus
path_expr_scillus(CK.harmony, PPA, 'PPA')


###AUX1-Auxins ---LOC121207860 absent
AUX <- c('LOC107889766','LOC107890973','LOC107891563',
         'LOC107899523','LOC107901395','LOC107908661',
         'LOC107916752','LOC107918785','LOC107920883',
         'LOC107927876','LOC107930145','LOC107934930',
         'LOC107941173','LOC107943649','LOC107955165',
         'LOC107958047','LOC107962491')

path_expr(CK.harmony,AUX,'AUX')
path_expr(LRC,AUX,'AUX')
DoHeatmap(CK.harmony, features = AUX,assay = 'RNA',group.by = 'celltype')
##sciluus
path_expr_scillus(CK.harmony, AUX, 'AUX')


###CRE1-Cytokinine ----LOC121210995 absent
CRE <- c('LOC107886429','LOC107909856','LOC107912602',
         'LOC107913518','LOC107918097','LOC107919969',
         'LOC107922254','LOC107924440','LOC107934511',
         'LOC107938547','LOC107946265','LOC107952812',
         'LOC107954738')

path_expr(CK.harmony,CRE,'CRE')
path_expr(LRC,CRE,'CRE')
##sciluus
path_expr_scillus(CK.harmony, CRE, 'CRE')


####ETR-Ethylene  ----LOC121214464 absent
ETR <- c('LOC107889254','LOC107903992','LOC107912356',
         'LOC107923931','LOC107925650','LOC107926811',
         'LOC107927738','LOC107928557','LOC107929124',
         'LOC107932670','LOC107936175','LOC107946055',
         'LOC107948437')

path_expr(CK.harmony,ETR,'ETR')
path_expr(LRC,ETR,'ETR')
##sciluus
path_expr_scillus(CK.harmony, ETR, 'ETR')

###GID1-Gibberellins
GID <- c('LOC107887986','LOC107888297','LOC107888298',
         'LOC107893263','LOC107893264','LOC107894563',
         'LOC107899135','LOC107912396','LOC107924318',
         'LOC107930165','LOC107936899','LOC107941282',
         'LOC107943945')

path_expr(CK.harmony,GID,'GID')
path_expr(LRC,GID,'GID')
##sciluus
path_expr_scillus(CK.harmony, GID, 'GID')


###JAR1-Jasmonates ----LOC107909788 absent
JAR <- c('LOC107900937','LOC107909787',
         'LOC107923052','LOC107923326','LOC107933526',
         'LOC107952940','LOC107954506')

path_expr(CK.harmony,JAR,'JAR')
path_expr(LRC,JAR,'JAR')
##sciluus
path_expr_scillus(CK.harmony, JAR, 'JAR')


### NPR1-Salicylic acid ---LOC107889382, LOC107891351, LOC121220834 absent
NPR <- c('LOC107888467',
         'LOC107891749','LOC107897478','LOC107916443',
         'LOC107921109','LOC107932804','LOC107939029',
         'LOC107940188')

path_expr(CK.harmony,NPR,'NPR')
path_expr(LRC,NPR,'NPR')
##sciluus
path_expr_scillus(CK.harmony, NPR, 'NPR')

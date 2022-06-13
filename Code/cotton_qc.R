library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
CKG.path <- "C:/BGI/mianhua/CKG"
CKN.path <- "C:/BGI/mianhua/CKN"
ZMSG.path <- "C:/BGI/mianhua/ZMSG"
ZMSN.path <- "C:/BGI/mianhua/ZMSN"
sample.ID <- c('CKG','CKN','ZMSG','ZMSN')
cotton.path <- "C:/BGI/mianhua/"

### load data ###
# load samples samples separately #
CKG.counts <- Read10X(data.dir = CKG.path)
CKN.counts <- Read10X(data.dir = CKN.path)
ZMSG.counts <- Read10X(data.dir = ZMSG.path)
ZMSN.counts <- Read10X(data.dir = ZMSN.path)

# initialize Seurat objects #
for(i in 1:length(sample.ID)){
  data <- sample.ID[i]
  count_data <- paste0(sample.ID[i],".counts")
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count_data,", project = \"Cotton\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count_data, ")")))
  
  # add metadata: study/cohort
  #eval(parse(text=paste0(data,"$study <- \"liao_kidney\"")))
  #eval(parse(text=paste0(data,"$cohort <- \"Guangxi\"")))
  
  # add metadata: condition
  #eval(parse(text=paste0(data,"$condition <- \"healthy\"")))
  
  
  # add meta data: calculate mitochondrion
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^ACQ\")")))
  
  # add meta data: calculate chloroplast
  eval(parse(text=paste0(data,"[[\"percent.pt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^GohiC\")")))
}

#### add condition to metadata 
CKG$Condition = "CKG"
CKN$Condition = "CKN"
ZMSG$Condition = "ZMSG"
ZMSN$Condition = "ZMSN"



#### Merge CKG & CKN into CK, ZMSG & ZMSN into ZMS
CK <- merge(CKG, CKN)
ZMS <- merge(ZMSG, ZMSN)

### perform QC ###
### functions ###
# scatter plots
QC.scatter <- function(data){
  scatter.plot <- ggplot(data[[]], aes( x = nCount_RNA, y = nFeature_RNA)) + 
    geom_point(aes(colour = percent.mt), size = 1) + 
    coord_cartesian(xlim = c(0.0 , 50000), ylim = c(0.0 , 10000)) +
    labs(title = "Overall QC", x  ="Count depth", y = "Unique Genes") + 
    theme(
      plot.title = element_text(color = "black", size = 20 , face = "bold"),
      axis.title.x = element_text(color = "black", size = 20, face = "bold"),
      axis.title.y = element_text(color = "black", size = 20, face = "bold"),
      legend.title = element_text(color = "black", size = 16, face = "bold", angle = 90)
    ) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                           guide = guide_colourbar("Mitochondrial fraction", title.position = "right", title.vjust = 0.5, title.hjust = 0.5, barwidth = 1.0, barheight = 60))
  
  return(scatter.plot)
}

# histograms
QC.histograms <- function(data){
  histograms <- list()
  
  # distribution of genes per cell
  hist1 <- qplot(x =data[["nFeature_RNA"]]$nFeature_RNA , fill=..count.., geom="histogram", binwidth = 100,
                 xlab = "Unique genes per cell",
                 ylab = "Frequency",
                 main = "Gene Count Distribution")+scale_fill_gradient(low="lightblue", high="darkblue")
  histograms[[1]] <- hist1
  
  # distribution of count depth
  hist2 <- qplot(x =data[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 500,
                 xlab = "Count depth per cell",
                 ylab = "Frequency",
                 main = "Transcript Count Distribution")+scale_fill_gradient(low="orange", high="red")
  histograms[[2]] <- hist2
  
  # distribution of mitochondrial gene fraction
  hist3 <- qplot(x =data[["percent.mt"]]$percent.mt, fill=..count.., geom="histogram", binwidth = 0.1,
                 xlab = "Mitochondrial fraction per cell",
                 ylab = "Frequency",
                 main = "Mitochondrial Reads Fraction")+scale_fill_gradient(low="lightgreen", high="darkgreen")
  histograms[[3]] <- hist3
  return(histograms)
}

### QC ###
# CKG_ unique genes per cell # unique transcripts per cell #mitochondrial fraction
nFeature.plot1 <- VlnPlot(CKG, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.raw.CKG.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(CKG, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.raw.CKG.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(CKG, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.raw.CKG.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

### CKN
nFeature.plot1 <- VlnPlot(CKN, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.raw.CKN.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(CKN, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.raw.CKN.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(CKN, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.raw.CKN.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

#ZMSG
nFeature.plot1 <- VlnPlot(ZMSG, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.raw.ZMSG.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(ZMSG, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.raw.ZMSG.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(ZMSG, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.raw.ZMSG.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

#ZMSN
nFeature.plot1 <- VlnPlot(ZMSN, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.raw.ZMSN.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(ZMSN, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.raw.ZMSN.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(ZMSN, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.raw.ZMSN.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

#CK overall
nFeature.plot1 <- VlnPlot(CK, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.raw.CK.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(CK, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.raw.CK.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()


#ZMS overall
nFeature.plot1 <- VlnPlot(ZMS, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.raw.ZMS.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(ZMS, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.raw.ZMS.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

## QC histograms----CKG ##
QC.histograms.raw <- QC.histograms(CKG)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(cotton.path,"QC.raw.histogram.CKG.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot----CKG ##
QC.scatter.raw <- QC.scatter(CKG)
png(paste0(cotton.path,"QC.raw.scatter.CKG.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

## QC histograms----CKN ##
QC.histograms.raw <- QC.histograms(CKN)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(cotton.path,"QC.raw.histogram.CKN.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot----CKN ##
QC.scatter.raw <- QC.scatter(CKN)
png(paste0(cotton.path,"QC.raw.scatter.CKN.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

## QC histograms----ZMSG ##
QC.histograms.raw <- QC.histograms(ZMSG)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(cotton.path,"QC.raw.histogram.ZMSG.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot----ZMSG ##
QC.scatter.raw <- QC.scatter(ZMSG)
png(paste0(cotton.path,"QC.raw.scatter.ZMSG.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

## QC histograms----ZMSN ##
QC.histograms.raw <- QC.histograms(ZMSN)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(cotton.path,"QC.raw.histogram.ZMSN.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot----ZMSN ##
QC.scatter.raw <- QC.scatter(ZMSN)
png(paste0(cotton.path,"QC.raw.scatter.ZMSN.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()


### save data: 25279 cells ###
saveRDS(kidney, paste0(cotton.path, "liao_kidney.rds"))


### remove low quality cells ###
### CKN ###
### cells with < 5000 features
### initial try to determine the min features as the qc plots are vague
CKN.filtered <- subset(CKN, subset = nFeature_RNA < 10000 & nFeature_RNA > 400 & nCount_RNA < 15000)
CKG.filtered <- subset(CKG, subset = nFeature_RNA < 15000 & nFeature_RNA > 400 & nCount_RNA < 20000)
ZMSN.filtered <- subset(ZMSN, subset = nFeature_RNA < 15000 & nFeature_RNA > 400 & nCount_RNA < 20000)
ZMSG.filtered <- subset(ZMSG, subset = nFeature_RNA < 3000 & nFeature_RNA > 400)


CK.filtered <- subset(CK, subset = nFeature_RNA <5000 & nFeature_RNA > 400 & nCount_RNA < 15000)
ZMS.filtered <- subset(ZMS, subset = nFeature_RNA <4000 & nFeature_RNA > 400 & nCount_RNA < 15000)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(CK.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(cotton.path,"QC.filtered.histogram.CK.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(CK.filtered)
png(paste0(cotton.path,"QC.filtered.scatter.CK.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()


#### ZMS
QC.histograms.filtered <- QC.histograms(ZMS.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(cotton.path,"QC.filtered.histogram.ZMS.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(ZMS.filtered)
png(paste0(cotton.path,"QC.filtered.scatter.ZMS.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()


## post filtering QC  ##
# CKG_ unique genes per cell # unique transcripts per cell #mitochondrial fraction
nFeature.plot1 <- VlnPlot(CKG.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.filtered.CKG.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(CKG.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.filtered.CKG.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(CKG.filtered, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.filtered.CKG.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

### CKN
nFeature.plot1 <- VlnPlot(CKN.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.filtered.CKN.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(CKN.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.filtered.CKN.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(CKN.filtered, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.filtered.CKN.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

#ZMSG
nFeature.plot1 <- VlnPlot(ZMSG.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.filtered.ZMSG.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(ZMSG.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.filtered.ZMSG.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(ZMSG.filtered, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.filtered.ZMSG.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

#ZMSN
nFeature.plot1 <- VlnPlot(ZMSN.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.filtered.ZMSN.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(ZMSN.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.filtered.ZMSN.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

nFeature.plot3 <- FeatureScatter(ZMSN.filtered, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
png(paste0(cotton.path,"QC.filtered.ZMSN.3.png"), width=600,height=500,units="px")
print(nFeature.plot3)
dev.off()

#CK overall
nFeature.plot1 <- VlnPlot(CK.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.filtered.CK.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(CK.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.filtered.CK.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()


#ZMS overall
nFeature.plot1 <- VlnPlot(ZMS.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"),  pt.size = 0.5)
png(paste0(cotton.path,"QC.filtered.ZMS.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(ZMS.filtered, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0)
png(paste0(cotton.path,"QC.filtered.ZMS.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

### save CK data: 21842 healthy cells ###
saveRDS(CK.filtered, paste0(cotton.path, "CK_filtered.rds"))

### save ZMS data: 342267 healthy cells ###
saveRDS(ZMS.filtered, paste0(cotton.path, "ZMS_filtered.rds"))

### save CKN data: 10402 cells ###
saveRDS(CKN.filtered, paste0(cotton.path, "CKN_filtered.rds"))

### save CKG data: 11670 cells ###
saveRDS(CKG.filtered, paste0(cotton.path, "CKG_filtered.rds"))

### save ZMSN data: 12862 cells ###
saveRDS(ZMSN.filtered, paste0(cotton.path, "ZMSN_filtered.rds"))

### save ZMSG data: 329985 cells ###
saveRDS(ZMSG.filtered, paste0(cotton.path, "ZMSG_filtered.rds"))

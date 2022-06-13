# R code
# Lukas Simon
# This code will generate panels for Figure 3

# Load R libraries ####
library(biomaRt)
library(limma)
library(Seurat)
library(pheatmap)
library(preprocessCore)
library(corrplot)
library(ggplot2)
library(DESeq2)
library(mcr)

# Load whole lung bulk expression ####
fc <- read.csv('C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/bulk_scRNA/Bulk_Count/gene_count.csv',na.strings = "-")
counts <- fc[, c(11:13,23:26)]
counts <- counts[grep(pattern="LOC",counts[,7]),]
rownames(counts) <- counts$gene_name
counts <- counts[,-7]
counts <- counts[which(rowSums(counts) > 0),]
bulk <- counts
CK_bulk <- c(rep("CKG", 3), rep("CKN", 3))

# Perform DE analysis on whole lung bulk data ###
des <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(CK_bulk), design = ~ CK_bulk)
des <- DESeq(des)
res <- results(des, contrast = c('CK_bulk', 'CKG', 'CKN')) 

# Get gene symbols from biomart ####
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
t2g <- t2g[which(t2g$ensembl_gene_id %in% rownames(bulk)),]

# Save DE result table ####
res$symbol <- t2g$external_gene_name[match(rownames(res), t2g$ensembl_gene_id)]
final <- data.frame(res, counts(des, normalized = T))

# Perform bulk deconvolution on whole lung results ####
markers <- read.delim('../data/Table S1_AllMarkersCelltypes.txt')
tmp <- markers[which(markers$avg_logFC > 1 & markers$p_val_adj < 0.1),]
schiller_genes <- split(as.character(tmp$gene), tmp$cluster)

tmp <- do.call(rbind, lapply(schiller_genes, function(x){
  ok <- intersect(x, res$symbol)
  set <- res$log2FoldChange[match(ok, res$symbol)]
  rest <- res$log2FoldChange[-match(ok, res$symbol)]
  pval <- ks.test(set, rest)$p.value
  coef <- mean(set) - mean(rest)
  c(fc_diff = coef, pval = pval)
}))
tmp[which(tmp[,2] < 1e-50), 2] <- 1e-50
tmp <- data.frame(tmp)
tmp$celltype <- rownames(tmp)

# Generate Fig 4e ####
ggplot(tmp, aes(fc_diff, -log10(pval))) + geom_point() +
  geom_text_repel(data = tmp, aes(label = celltype), color = "black", size = 3) +
  geom_vline(xintercept = 0)

# Generate Fig 4f ####
genECDFplot2 <- function(celltype){
  ok <- intersect(schiller_genes[[celltype]], res$symbol)
  type <- rep('Rest', nrow(res))
  type[match(ok, res$symbol)] <- celltype 
  df <- data.frame(foldchange = res$log2FoldChange, type)
  ggplot(df, aes(foldchange, colour = type)) + scale_color_manual(values=c("red", "black")) + stat_ecdf() + scale_x_continuous(limits = c(-10, 10)) 
}
genECDFplot2('Ciliated_cells')

# Load and generate in silico bulk data ####
Cotton <- readRDS("C:/BGI/mianhua/CK_harmony_updated_test/combined_CK_annotated_doubletremoved.rds")

CKG <- read_count_output("C:/BGI/mianhua/CKG/")
CKG <- cbind(apply(CKG,1,sum))

CKN <- read_count_output('C:/BGI/mianhua/CKN/')
CKN <- cbind(apply(CKN,1,sum))

asplit <- split(1:nrow(Cotton@meta.data), Cotton@meta.data$Sample)
tmp <- do.call(cbind, lapply(asplit, function(x) Matrix::rowSums(Cotton@assays$RNA@counts[,x])))
tmp <- tmp[which(rowSums(tmp) > 0),]
insilico <- tmp
CK_insilico <- as.character(Cotton@meta.data$Sample[match(colnames(insilico), Cotton@meta.data$Sample)])
# CK_insilico <- gsub("10Days", "young", CK_insilico)
# CK_insilico <- gsub("6Months", "old", CK_insilico)


# Show correspondence between bulk and insilico ####
ok <- intersect(rownames(insilico), rownames(bulk))
ok <- intersect(names(CKN), rownames(bulk))
tmp <- bulk[rownames(bulk)[match(ok, rownames(bulk))], ]

# Generate Fig 3b - part 1 ####
tmp2 <- cbind(tmp, insilico[ok, ])
correl <- cor(tmp2, method = 'spearman')
corrplot(correl, method = 'ellipse', type = 'upper', cl.lim = c(0,1))

# Generate Fig 3b - part 2 ####
## compare ckn vs ckn bulk
insilico_CKN <- insilico[,2] %>% as.matrix()
colnames(insilico_CKN) <- 'CKN'
insilico_CKN_means <- insilico_CKN[ok,]
bulk_CKN <- tmp[,c(4:6)]
bulk_CKN_means <- rowMeans(bulk_CKN)

plot(log(insilico_CKN_means), log(bulk_CKN_means), col = rgb(0, 0, 0, 0.5), pch = 19, xlab = 'Average in silico bulk [log]', ylab = 'Average whole lung bulk [log]')
cor.test(log(insilico_CKN_means), log(bulk_CKN_means), method = 'spearman',exact = FALSE)
dem.reg <- mcreg(log(insilico_CKN_means), log(bulk_CKN_means), method.reg = "Deming")
abline(dem.reg@para[1:2], col = "red", lwd = 2)
abline(0, 1)

## compare ckg vs ckg bulk
insilico_CKG <- insilico[,1] %>% as.matrix()
colnames(insilico_CKG) <- 'CKG'
insilico_CKG_means <- insilico_CKG[ok,]
bulk_CKG <- tmp[,c(1:3)]
bulk_CKG_means <- rowMeans(bulk_CKG)

plot(log(insilico_CKG_means), log(bulk_CKG_means), col = rgb(0, 0, 0, 0.5), pch = 19, xlab = 'Average in silico bulk [log]', ylab = 'Average whole lung bulk [log]')
cor.test(log(insilico_CKG_means), log(bulk_CKG_means), method = 'spearman',exact = FALSE)
dem.reg <- mcreg(log(insilico_CKG_means), log(bulk_CKG_means), method.reg = "Deming")
abline(dem.reg@para[1:2], col = "red", lwd = 2)
abline(0, 1)
insilico_means <- rowMeans(insilico[ok,])
bulk_means <- rowMeans(tmp)
plot(log(insilico_means), log(bulk_means), col = rgb(0, 0, 0, 0.5), pch = 19, xlab = 'Average in silico bulk [log]', ylab = 'Average whole lung bulk [log]')
cor.test(log(insilico_means), log(bulk_means), method = 'spearman',exact = FALSE)
dem.reg <- mcreg(log(insilico_means), log(bulk_means), method.reg = "Deming")
abline(dem.reg@para[1:2], col = "red", lwd = 2)
abline(0, 1)


# Match genes for all three matrices ####
ok <- intersect(rownames(protein), rownames(insilico))
#ok <- intersect(ok, t2g$external_gene_name)
protein_ok <- protein[ok,]
write.csv(protein_ok, "C:/BGI/转录噪声/Proteome/Overlapped_protein_abundance.csv")
#bulk_ok <- bulk[t2g$ensembl_gene_id[match(ok, t2g$external_gene_name)], ]
insilico_ok <- insilico[ok, ]
write.csv(insilico_ok, "C:/BGI/转录噪声/Proteome/Overlapped_gene_sum.csv")
rownames(protein_ok) <- rownames(insilico_ok) <- ok

# Normalize and merge all three matrices ####
#bulk_ok <- voom(bulk_ok)$E
insilico_ok <- voom(insilico_ok)$E
#merged <- cbind(bulk_ok, insilico_ok, protein_ok)
merged <- cbind(insilico_ok, protein_ok)
merged <- normalize.quantiles(merged)
rownames(merged) <- rownames(insilico_ok)
merged <- merged[-46,]

# Define sample attributes ####
batch <- c(rep("insilico", ncol(insilico_ok)), rep("protein", ncol(protein_ok)))
age <- c(age_insilico, age_prot)

farben <- c("blue", "red")
names(farben) <- c('young', 'old')

shape <- c( 2, 3)
names(shape) <- c( "insilico", "protein")

# Calculate PCA ####
pca <- prcomp(t(merged))
# Generate Fig 3c ####
par(mfrow = c(1, 2))
par(oma=c(3,3,3,3)) 
par(mar=c(6,5,4,3) + 0.1) 
plot(pca$x[,1:2], col = farben[age], pch = shape[batch])
#legend("bottomright", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c("blue", "red"))
abline(h = 0,v = 0,lty = 3)
plot(pca$x[,2:3], col = farben[age], pch = shape[batch])
abline(h = 0,v = 0,lty = 3)
legend('topright',c(names(farben), names(shape)), xpd= FALSE, pch = c(16, 16, shape), bty = "n", col = c("blue", "red"))

### using ggplot2 for PCA ###



# Generate Fig 3f ####
col4a_genes <- c('Col4a1', 'Col4a2', 'Col4a3', 'Col4a4', 'Col4a5', 'Col4a6')
norm1 <- function(matr, treat){
  
  t(apply(matr, 1, function(x) (x - mean(x))/sd(x)))
}

batch <- c(rep("insilico", ncol(insilico_ok)), rep("bulk", ncol(bulk_ok)),rep("protein", ncol(protein_ok)))
age <- c(age_insilico, age_bulk, age_prot)

final <- cbind(norm1(insilico_ok[col4a_genes,]),
               norm1(bulk_ok[col4a_genes,]),
               norm1(protein_ok[col4a_genes,]))

tmp <- do.call(c, apply(final, 1, function(x) split(x, paste(batch, age))))
ord <- c(seq(from = 1, to = 36, by = 6),
         seq(from = 2, to = 36, by = 6),
         seq(from = 3, to = 36, by = 6),
         seq(from = 4, to = 36, by = 6),
         seq(from = 5, to = 36, by = 6),
         seq(from = 6, to = 36, by = 6))
tmp <- tmp[ord]
order4plot <- c(paste(col4a_genes, 'insilico young', sep = '.'),
                paste(col4a_genes, 'insilico old', sep = '.'),
                paste(col4a_genes, 'bulk young', sep = '.'),
                paste(col4a_genes, 'bulk old', sep = '.'),
                paste(col4a_genes, 'protein young', sep = '.'),
                paste(col4a_genes, 'protein old', sep = '.'))

boxplot(tmp[order4plot], las = 2, ylab = 'Expression z-score', outline = F)
abline(h = 0)


###
read_count_output <- function(dir) {
  #dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "matrix.mtx.gz"))
  genes <- read.table(paste0(dir, "/", 'features.tsv.gz'), stringsAsFactors = F,sep='\t',header = F)$V2
  #barcodes <- readLines(file(paste0(dir, "/", "barcodes.tsv.gz")))
  #colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

CKG <- read_count_output("C:/BGI/mianhua/CKG/")
CKG <- cbind(apply(CKG,1,sum))

CKN <- read_count_output('C:/BGI/mianhua/CKN/')
CKN <- cbind(apply(CKN,1,sum))

# Load whole lung bulk expression ####
fc <- read.csv('C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/bulk_scRNA/Bulk_Count/gene_count.csv',na.strings = "-")
counts <- fc[, c(11:13,23:26)]
counts <- counts[grep(pattern="LOC",counts[,7]),]
rownames(counts) <- counts$gene_name
counts <- counts[,-7]
bulk <- counts

sharedGenes <- intersect(rownames(CKG), rownames(bulk))


bulk <- bulk[sharedGenes,]
CKN <- CKN[sharedGenes,]
CKG <- CKG[sharedGenes,]


all <- cbind(bulk)
#all <- all/gL[rownames(all)]
## RPM
all <- t(t(all)/(apply(all,2,sum)/1e6))
## RPM +1 
all <- all + 1
all <- cbind(all,CKN,CKG)
all <- cbind(all,CKN+1,CKG+1)
colnames(all)[7:8] <- c('CKN', 'CKG')
all <- all[rowSums(all == 0) == 0,]


all <- log(all)
#all[is.na(all) | all == "-Inf"] <- NA

cor(all,method = "sp")
all <- as.data.frame(all)
all$CKG_Average <- (all$CKG_DR1 + all$CKG_DR2+all$CKG_DR3)/3
all$CKN_Average <- (all$CKN_DR1 + all$CKN_DR2+all$CKN_DR3)/3

## Compare with bulk seq
over_gene_ckg <- intersect(rownames(CKG), rownames(bulk))
over_gene_ckn <- intersect(rownames(CKN), rownames(bulk))
sub_ckg  <- CKG[over_gene_ckg, ]
sub_ckg_means <- cbind(apply(sub_ckg,1,sum))
sub_ckg_means <- sub_ckg_means[which(rowSums(sub_ckg_means) > 0),]

sub_ckn  <- CKN[over_gene_ckn, ]
sub_ckn_means <- cbind(apply(sub_ckn,1,sum))
sub_ckn_means <- sub_ckn_means[which(rowSums(sub_ckn_means) > 0),]

##bulk
bulk_CKG <- bulk[1:3]
bulk_CKG <- bulk_CKG[over_gene_ckg,]
bulk_CKG_means <- rowSums(bulk_CKG)


bulk_CKN <- bulk[4:6]
bulk_CKN <- bulk_CKN[over_gene_ckn,]
bulk_CKN_means <- rowSums(bulk_CKN)


pdf("Bulk_scRNA_CKG.pdf", width = 3.3*3, height = 2.5*2)
layout(matrix(c(1,2,5,6,3,4,7,8), nrow = 2, ncol = 4, byrow = TRUE))
par(mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
#plot(log(sub_ckn_means), log(bulk_CKN_means), col = rgb(0, 0, 0, 0.5), pch = 19, xlab = 'Average in silico bulk [log]', ylab = 'Average whole lung bulk [log]')
plot(all[,"CKG"],all[,"CKG_DR1"], pch = 16, 
     col=densCols(cbind(all[,"CKG"],all[,"CKG_DR1"])), 
     xlab="", ylab="",las=1)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("CKG_sc vs. CKG_bulk1", side = 3, line = 0.1, cex = 0.7)
corV <- round(cor(all[,"CKG"],all[,"CKG_DR1"], method = 'sp'),2)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKG_DR1"]~all[,"CKG"]), col = "red")
#abline(0, 1)

plot(all[,"CKG"],all[,"CKG_DR2"], pch = 16, 
     col=densCols(cbind(all[,"CKG"],all[,"CKG_DR2"])), 
     xlab="", ylab="",las=1)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("CKG_sc vs. CKG_bulk2", side = 3, line = 0.1, cex = 0.7)
corV <- round(cor(all[,"CKG"],all[,"CKG_DR2"], method = 'sp'),2)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKG_DR2"]~all[,"CKG"]), col = "red")

plot(all[,"CKG"],all[,"CKG_DR3"], pch = 16, 
     col=densCols(cbind(all[,"CKG"],all[,"CKG_DR3"])), 
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"CKG"],all[,"CKG_DR3"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('RPM+1'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("CKG_sc vs. CKG_bulk3", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKG_DR3"]~all[,"CKG"]), col = "red")

plot(all[,"CKG"],all[,"CKG_Average"], pch = 16, 
     col=densCols(cbind(all[,"CKG"],all[,"CKG_Average"])),
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"CKG"],all[,"CKG_Average"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('mean RPM+1'))"), side = 2, line = 1.5)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("CKG_sc vs. CKG_bulk_average", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKG_Average"]~all[,"CKG"]), col = "red")


#CKN
plot(all[,"CKN"],all[,"CKN_DR1"], pch = 16, 
     col=densCols(cbind(all[,"CKN"],all[,"CKN_DR1"])), 
     xlab="", ylab="",las=1)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("CKN_sc vs. CKN_bulk1", side = 3, line = 0.1, cex = 0.7)
corV <- round(cor(all[,"CKN"],all[,"CKN_DR1"], method = 'sp'),2)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKN_DR1"]~all[,"CKN"]), col = "red")
#abline(0, 1)

plot(all[,"CKN"],all[,"CKN_DR2"], pch = 16, 
     col=densCols(cbind(all[,"CKN"],all[,"CKN_DR2"])), 
     xlab="", ylab="",las=1)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("CKN_sc vs. CKN_bulk2", side = 3, line = 0.1, cex = 0.7)
corV <- round(cor(all[,"CKN"],all[,"CKN_DR2"], method = 'sp'),2)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKN_DR2"]~all[,"CKN"]), col = "red")

plot(all[,"CKN"],all[,"CKN_DR3"], pch = 16, 
     col=densCols(cbind(all[,"CKN"],all[,"CKN_DR3"])), 
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"CKN"],all[,"CKN_DR3"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('RPM+1'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("CKN_sc vs. CKN_bulk3", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKG_DR3"]~all[,"CKN"]), col = "red")

plot(all[,"CKN"],all[,"CKN_Average"], pch = 16, 
     col=densCols(cbind(all[,"CKN"],all[,"CKN_Average"])),
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"CKN"],all[,"CKN_Average"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('RPM+1'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('mean RPM+1'))"), side = 2, line = 1.5)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("CKN_sc vs. CKN_bulk_average", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CKN_Average"]~all[,"CKN"]), col = "red")
dev.off()




ggplot(data.frame(corrlation='corrlation',value=unlist(outcor[,1])),aes(x=corrlation,y=value))+geom_boxplot()+xlab('')+
  theme(axis.title = element_text(size=16),axis.text = element_text(size=16))


# 单细胞基因表达量与bulk基因表达量相关性分析
mat <- as.data.frame(t(as.matrix(GetAssayData(Cotton, assay = "RNA", slot = "count"))))
mat <- aggregate(mat, by=list(Cotton@meta.data[["Sample"]]), FUN="sum") 
mat = t(mat)
exp_sc = data.frame(mat[,-1])
sum(as.integer(exp_sc))
exp_sc['exp_sc']=log2(as.integer(exp_sc$mat..1...)/239.69+1)
exp_sc['gene_name'] = rownames(exp_sc)
write.csv(exp_sc,file = "D:/singlecell-liang/genome-baimaigen/rsem_result/root_exp_sc.csv")
# library(tidyverse)
# exp_sc = rename(exp_sc,log1p=mat..1...)
# exp = read.csv(file = "D:/singlecell-liang/genome-baimaigen/rsem_result/root_exp.csv")
# exp['log1p']= log1p(exp['RPM_final'])
# write.csv(exp,file = "D:/singlecell-liang/genome-baimaigen/rsem_result/root_exp.csv")
exp = read.csv(file = "D:/singlecell-liang/genome-baimaigen/rsem_result/root_exp.csv")
library(ggplot2)
library(ggpubr)
library(ggpmisc)
result=lm(formula =exp$unprotoplasted.Bulk.RNA.Expression.log2.rpm.1. ~ exp$ScRNA.Expression.log2.rpm.1., data = exp)
summary(result)
library(ggplot2)
theme_set(theme_bw())
a=summary(result)$r.squared
r = signif(sqrt(a),2)
p = signif(summary(result)$coefficients[2,4],2)


theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
b <- ggplot(exp, aes(x = ScRNA.Expression.log2.rpm.1., y = unprotoplasted.Bulk.RNA.Expression.log2.rpm.1.))
# Scatter plot with regression line
p1 = b + geom_point()+
  geom_smooth(method = "lm", color = "red", fill = "lightgray") 
p2=p1+labs(title=paste("r=",r,"p-value=",p))
p2

result=lm(formula =exp$protoplasted.Bulk.RNA.Expression.log2.rpm.1. ~ exp$ScRNA.Expression.log2.rpm.1., data = exp)
summary(result)
library(ggplot2)
theme_set(theme_bw())
a=summary(result)$r.squared
r = signif(sqrt(a),2)
p = signif(summary(result)$coefficients[2,4],2)


theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
b <- ggplot(exp, aes(x = ScRNA.Expression.log2.rpm.1., y = protoplasted.Bulk.RNA.Expression.log2.rpm.1.))
# Scatter plot with regression line
p1 = b + geom_point()+
  geom_smooth(method = "lm", color = "red", fill = "lightgray") 
p2=p1+labs(title=paste("r=",r,"p-value=",p))
p2

library(Matrix)
library(geuvPack)
library(preprocessCore)
library(limma)


popInfo <- read.csv("populationInfo/GD660.GeneQuantCount.txt.gz", sep = "\t")
gNames <- popInfo[,2]
popInfo <- popInfo[gNames %in% geneInfo[,1],]
gNames <- gNames[gNames %in% geneInfo[,1]]
gNames <- geneInfo[as.vector(gNames),2]
colnames(popInfo) <- unlist(lapply(strsplit(colnames(popInfo), "\\."), function(X){X[1]}))
popInfo <- popInfo[,colnames(geuFPKM)]
rownames(popInfo) <- make.unique(gNames)

CEU <- popInfo[,geuFPKM@colData$popcode == "CEU",]
CEU <- CEU[rownames(CEU) %in% geneLength[,1],]
CEU <- CEU/gL[rownames(CEU)]
CEU <- t(t(CEU)/(colSums(CEU)/1e6))
CEU <- cbind(rowMeans(CEU))

YRI <- popInfo[,geuFPKM@colData$popcode == "YRI",]
YRI <- YRI[rownames(YRI) %in% geneLength[,1],]
YRI <- YRI/gL[rownames(YRI)]
YRI <- t(t(YRI)/(colSums(YRI)/1e6))
YRI <- cbind(rowMeans(YRI))

bulkRNAseq <- read.csv("bulkFQ/RNAseq.tsv", sep = "\t")
bulkRNAseq <- bulkRNAseq[rownames(bulkRNAseq) %in% geneLength[,1],]

GM12878 <- readMM("../GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("../GM12878/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM12878 <- cbind(apply(GM12878,1,sum))

GM18502 <- readMM("../GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("../GM18502/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM18502 <- cbind(apply(GM18502,1,sum))

sharedGenes <- table(c(make.unique(rownames(GM12878)),
                       make.unique(rownames(GM18502)),
                       make.unique(rownames(CEU)),
                       make.unique(rownames(YRI)),
                       make.unique(rownames(bulkRNAseq))))
sharedGenes <- names(sharedGenes[sharedGenes == 5])

bulkRNAseq <- bulkRNAseq[sharedGenes,]
GM18502 <- GM18502[sharedGenes,]
GM12878 <- GM12878[sharedGenes,]
CEU <- CEU[sharedGenes,]
YRI <- YRI[sharedGenes,]

all <- cbind(bulkRNAseq)
all <- all/gL[rownames(all)]
all <- t(t(all)/(apply(all,2,sum)/1e6))
all <- cbind(all,GM12878,GM18502,CEU,YRI)
all <- all[rowSums(all == 0) == 0,]

colnames(all) <- c("bGM12878", "bGM18502", "scGM12878", "scGM18502", "CEU", "YRI")
all <- log(all)

cor(all,method = "sp")
all <- as.data.frame(all)

pdf("Figure5.pdf", width = 3.3*3, height = 2.5*2)
layout(matrix(c(1,2,5,5,3,4,5,5), nrow = 2, ncol = 4, byrow = TRUE))
par(mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
plot(all[,"scGM12878"],all[,"bGM12878"], pch = 16, 
     col=densCols(cbind(all[,"scGM12878"],all[,"bGM12878"])), 
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"scGM12878"],all[,"bGM12878"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('TPM'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("GM12878 vs. GM12878", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"bGM12878"]~all[,"scGM12878"]), col = "red")

plot(all[,"scGM18502"],all[,"bGM18502"], pch = 16, 
     col=densCols(cbind(all[,"scGM18502"],all[,"bGM18502"])), 
     xlab="", ylab="",las=1)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('TPM'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("GM18502 vs. GM18502", side = 3, line = 0.1, cex = 0.7)
corV <- round(cor(all[,"scGM18502"],all[,"bGM18502"], method = 'sp'),2)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"bGM18502"]~all[,"scGM18502"]), col = "red")

plot(all[,"scGM12878"],all[,"CEU"], pch = 16, 
     col=densCols(cbind(all[,"scGM12878"],all[,"CEU"])), 
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"scGM12878"],all[,"CEU"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('average TPM'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("GM12878 vs. CEU Population", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CEU"]~all[,"scGM12878"]), col = "red")

plot(all[,"scGM18502"],all[,"YRI"], pch = 16, 
     col=densCols(cbind(all[,"scGM18502"],all[,"YRI"])),
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"scGM18502"],all[,"YRI"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('average TPM'))"), side = 2, line = 1.5)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("GM18502 vs. YRI Population", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"YRI"]~all[,"scGM18502"]), col = "red")

geneLength <- read.csv("hg38_genes_length.tsv", sep = "\t", stringsAsFactors = FALSE)
gL <- geneLength[,4]
names(gL) <- geneLength[,1]
geneInfo <- cbind(geuFPKM@rowRanges$gene_id,geuFPKM@rowRanges$gene_name)
rownames(geneInfo) <- geneInfo[,1]
popInfo <- read.csv("populationInfo/GD660.GeneQuantCount.txt.gz", sep = "\t")
gNames <- popInfo[,2]
popInfo <- popInfo[gNames %in% geneInfo[,1],]
gNames <- gNames[gNames %in% geneInfo[,1]]
gNames <- geneInfo[as.vector(gNames),2]
colnames(popInfo) <- unlist(lapply(strsplit(colnames(popInfo), "\\."), function(X){X[1]}))
popInfo <- popInfo[,colnames(geuFPKM)]
rownames(popInfo) <- make.unique(gNames)

CEU <- popInfo[,geuFPKM@colData$popcode == "CEU",]
CEU <- CEU[rownames(CEU) %in% geneLength[,1],]

YRI <- popInfo[,geuFPKM@colData$popcode == "YRI",]
YRI <- YRI[rownames(YRI) %in% geneLength[,1],]

bulkRNAseq <- read.csv("bulkFQ/RNAseq.tsv", sep = "\t")
bulkRNAseq <- bulkRNAseq[rownames(bulkRNAseq) %in% geneLength[,1],]

GM12878 <- readMM("../GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("../GM12878/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM12878 <- cbind(apply(GM12878,1,sum))

GM18502 <- readMM("../GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("../GM18502/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM18502 <- cbind(apply(GM18502,1,sum))

sharedGenes <- table(c(make.unique(rownames(GM12878)),
                       make.unique(rownames(GM18502)),
                       make.unique(rownames(CEU)),
                       make.unique(rownames(YRI)),
                       make.unique(rownames(bulkRNAseq))))
sharedGenes <- names(sharedGenes[sharedGenes == 5])

bulkRNAseq <- bulkRNAseq[sharedGenes,]
GM18502 <- GM18502[sharedGenes,]
GM12878 <- GM12878[sharedGenes,]

CEU <- CEU[sharedGenes,]
YRI <- YRI[sharedGenes,]

colnames(CEU) <- paste0("CEU_",colnames(CEU))
colnames(YRI) <- paste0("YRI_", colnames(YRI))

all <- cbind(bulkRNAseq,GM12878,GM18502, CEU, YRI)
dType <- c(rep("b",2), rep("sc",2), rep("geu", ncol(CEU)), rep("geu", ncol(YRI)))
dType <- as.factor(dType)
dType <- relevel(dType,ref="sc")
all <- removeBatchEffect(all,batch = dType)

all <- all/gL[rownames(all)]
all <- round(t(t(all)/(apply(all,2,sum))) * 1e6)
all <- all[rowSums(all == 0) == 0,]

colnames(all)[1:4] <- c("bGM12878", "bGM18502", "scGM12878", "scGM18502")
allN <- normalize.quantiles(all)
colnames(allN) <- colnames(all)
rownames(allN) <- rownames(all)
all <- round(allN)

cList <- rep("gray", ncol(all))
cList[grepl("YRI", colnames(all))] <- rgb(1,0,0,0.5)
cList[grepl("scGM18502", colnames(all))] <- rgb(1,0,0,1)
cList[grepl("bGM18502", colnames(all))] <- rgb(1,0,0,1)
cList[grepl("CEU", colnames(all))] <- rgb(0,0,1,0.5)
cList[grepl("scGM12878", colnames(all))] <- rgb(0,0,1,1)
cList[grepl("bGM12878", colnames(all))] <- rgb(0,0,1,1)

par(mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
outData <- prcomp(scale(t(all)))$x[,1:2]
colnames(outData) <- c("PC1","PC2")
plot(outData, pch=16, col=cList, las=1, xlab = "", ylab = "", cex = 1.5)
mtext(parse(text = "PC1"), side = 1,  line = 2)
mtext(parse(text = "PC2"), side = 2,  line = 1.5)
mtext("PCA", side = 3, line = 1, cex = 1,font = 2)
mtext("ALL SAMPLES", side = 3, line = 0.1, cex = 0.7)
text(outData[1,1],outData[1,2]+10,"BULK-GM12878")
arrows(outData[1,1],outData[1,2]+2, outData[1,1],outData[1,2]+7, length = 0)
text(outData[2,1],outData[2,2]-10,"BULK-GM18502")
arrows(outData[2,1],outData[2,2]-2, outData[2,1],outData[2,2]-7, length = 0)
text(outData[3,1]+15,outData[3,2]+10,"SC-GM12878")
arrows(outData[3,1]+0.5,outData[3,2]+2, outData[3,1]+2,outData[3,2]+7, length = 0)
text(outData[4,1]+30,outData[4,2]-10,"SC-GM18502")
arrows(outData[4,1]+2,outData[4,2]-2, outData[4,1]+11,outData[4,2]-9, length = 0)
abline(h=0, lty=2, col=rgb(0,0,0,0.1))
legend("bottomright", legend = c("CEU", "YRI"), bty = "n", horiz = TRUE, pch = 16, col = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)))
dev.off()

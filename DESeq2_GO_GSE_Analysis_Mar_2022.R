#!/usr/local/bin/Rscript
# R 4.0.2
#---------------------------------------------------------------------------------
# Human BeWo cells cultured ± thapsigargin RNASeq analysis,  25.06.2021
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/Yung_Charnock-Jones
#
#
# Analysis Performed by Xiaohui Zhao
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#---------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+--------------------------------------------------------------------+")
message("+------        Basic settings and libraries           ---------------+")
message("+--------------------------------------------------------------------+")

suppressPackageStartupMessages({
  library("tidyr")
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("Matrix")
  library("matrixStats")
  library("useful")
  library("reshape")
  library("reshape2")
  library("DESeq2")
  library("biomaRt")
  library("ggforce")
  library("pheatmap")
  library('RColorBrewer')
  library("scales")
  library("ggbeeswarm")
  library("BiocParallel")
  library("ggalt")
  library("ComplexHeatmap")
  library("apeglm")
  library("openxlsx")
})

NUMCORES      <- 3
register(MulticoreParam(NUMCORES))

Project <- "CTR_gjb2_0002"
baseDir <- "/Users/xz289/Documents/CTR_gjb2_0002/Figures_Tables_Ensembl_ID_98_84/deseq2_qc"
setwd(baseDir)

TOPNUM       <- 2000
l2fc         <- 1
significance <- 0.05
elementTextSize <- 6

message("+-------------------------------------------------------------------------------+")
message("+                               Use ensEMBL Annotations                         +")
message("+-------------------------------------------------------------------------------+")

ensembl    =  useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host = 'ensembl.org')
listEnsembl()
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description','chromosome_name'),  mart = ensembl, useCache = FALSE) 
ensEMBL2id$description <- gsub("..Source.*", "", ensEMBL2id$description)
ensEMBL2id <- subset(ensEMBL2id, chromosome_name=="1"| chromosome_name=="2"| chromosome_name=="3"|
                       chromosome_name=="4"| chromosome_name=="5"| chromosome_name=="6"|
                       chromosome_name=="7"| chromosome_name=="8"| chromosome_name=="9"|
                       chromosome_name=="10"| chromosome_name=="11"| chromosome_name=="12"|
                       chromosome_name=="13"| chromosome_name=="14"| chromosome_name=="15"|
                       chromosome_name=="16"| chromosome_name=="17"| chromosome_name=="18"|
                       chromosome_name=="19"| chromosome_name=="20"| chromosome_name=="21"|
                       chromosome_name=="22"| chromosome_name=="MT"| chromosome_name=="X"|
                       chromosome_name=="Y")
head(ensEMBL2id)
nrow(ensEMBL2id)
## 67128--->60605
save(ensEMBL2id, file = "Ensembl_hsapiens_ID_Name_Des_Chr.RData")

message("+-------------------------------------------------------------------------------+")
message("+    DESeq2 Analaysis---->PCA plot                                              +")
message("+-------------------------------------------------------------------------------+")

load("deseq2.dds.RData")  ## 60504
cts                  <- assay(dds)
coldata              <- colData(dds)
coldata$Pairs        <- as.factor(rep(c(1:5), length=10))
coldata$condition    <- as.factor(ifelse(coldata$Group1=="control", "Control", "Tg"))
coldata$Individual   <- paste0("A0", c("03", "10", "18", "21", "25", "01", "09", "16", "20", "23"))

message("+----          DESeq2 model with paired information in the design           ----+")

dds.pair             <- DESeqDataSetFromMatrix(countData = cts,
                                               colData = coldata,
                                               design= ~ Pairs + condition)
dds.pair             <- estimateSizeFactors(dds.pair)
sizeFactors(dds.pair)

dds.pair             <- DESeq(dds.pair, parallel=TRUE)
vsd.pair             <- vst(dds.pair,     blind=F)
resultsNames(dds.pair)
colData(vsd.pair)


message("+----                     DESeq2 PCA plot                                   -------+")

customPCA     <- function(sampleTBL, RLD, TOPNUM, model) {
  
  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sample, pca$x, individual=sampleTBL$Individual, 
                          condition=sampleTBL$condition, pairs=sampleTBL$Pairs)
  
  
  plt.pca.new <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, shape=pairs, label=individual) ) +
    geom_text_repel(aes(label=individual), show.legend = FALSE, size=4) +
    geom_point(size=4, alpha=0.75) +
    scale_color_manual(values=c("blue", "red")) +
    scale_shape_manual(values = c(15,16,17,18,20)) +
    xlab(pc1lab) +
    ylab(pc2lab) + 
    coord_fixed() +
    theme_bw() +
    theme(text = element_text(size=12), 
          plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white",colour = "white"),
          legend.position='right', 
          aspect.ratio=1,
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  return(list(plt.pca.new))
  
}
sampleTBL            <- as.data.frame(colData(dds.pair))
RLD.pair             <- assay(vsd.pair)

TOPNUM               <- 2000
models               <- "vsd.pair"

## customised PCA function

pca.pair      <- customPCA(sampleTBL, RLD.pair, TOPNUM, models)

pdf("SFig1_PCA_nofilter.pdf")
print(pca.pair[[1]])
dev.off()

message("+-------------------------------------------------------------------------------+")
message("+    DESeq2 Analaysis with paired design                                        +")
message("+-------------------------------------------------------------------------------+")

res.pair        <- lfcShrink(dds=dds.pair, coef="condition_Tg_vs_Control", type="apeglm", parallel=TRUE)
res.pair.dat    <- as.data.frame(res.pair)
res.pair.dat$ensembl_gene_id <- rownames(res.pair.dat)

res.pair.datEM    <- merge(res.pair.dat, ensEMBL2id, by ="ensembl_gene_id") ## 56681, 
res.pair.datEM    <- subset(res.pair.datEM, !is.na(padj))  ## 21076
## remove padj missing, left 21076, 17743 with external_gene_name
## 3333 missing external_gene_names, and 23 has duplicated names
## duplicated gene names 
## [1] "Y_RNA"  (10)    "SPATA13"    "RGS5"       "LINC00484"  "DNAJC9-AS1"
## [6] "ALG1L9P"    "LINC01238"  "TMSB15B"    "DUXAP8"     "BMS1P4"    
## [11] "CYB561D2"   "POLR2J4"    "SNORD22"    "MATR3"  

res.pair.datEM    <- res.pair.datEM[order(res.pair.datEM$ensembl_gene_id),
                                    c("ensembl_gene_id", "external_gene_name", "chromosome_name",
                                      "description", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
res.mat           <- as.data.frame(assay(vsd.pair))
res.mat$ensembl_gene_id <- rownames(res.mat)
res.matMer        <- merge(res.pair.datEM, res.mat, by = "ensembl_gene_id", all.x=T)
res.matMer        <- res.matMer[order(res.matMer$ensembl_gene_id), 
                                c("ensembl_gene_id", "external_gene_name", "chromosome_name",
                                  "description", "baseMean", "log2FoldChange", "lfcSE",
                                  "pvalue", "padj", 
                                  paste0(rep(c("Tg", "Control"), each=5), "_REP", rep(c(1:5), length=10)))]
colnames(res.matMer)[10:19] <- c( "Tg-A003","Tg-A010",  "Tg-A018", "Tg-A021", "Tg-A025",
                                  "Control-A001", "Control-A009", "Control-A016","Control-A020","Control-A023")


save(dds.pair, vsd.pair, res.pair, res.pair.datEM, res.matMer, 
     file = paste0(Project, "-DESeq2_condition_pair_dds_vsd_res_nofilter.RData"))

message("+---- filter out the gene with duplicated gene names and no names -------------+")

res.matMer.mind  <- subset(res.matMer, external_gene_name!="") 
res.matMer.mind  <- res.matMer.mind[-which(duplicated(res.matMer.mind$external_gene_name)==T),] ## 17275


write.table(res.matMer, file = paste0(Project, "-DESeq2_res_all_normvsd_N", dim(res.matMer)[1], "_summaryTable_nofilter.txt"), row.names=F) 
write.table(res.matMer.mind, file = paste0(Project, "-DESeq2_res_all_normvsd_N", dim(res.matMer.mind)[1], "_rmNA_rmDup_summaryTable_nofilter.txt"), row.names=F) 

## select significant DEGs, padj <- 0.05 & abs(l2fc) >=1

res.pair.sig    <- subset(res.matMer, abs(log2FoldChange) >=l2fc & padj < significance) 
print(dim(res.pair.sig)) 
res.sig.dim     <- dim(res.pair.sig)
res.pair.sig    <- res.pair.sig[order(-res.pair.sig$log2FoldChang), ]

res.pair.sig.mind    <- subset(res.matMer.mind, abs(log2FoldChange) >=l2fc & padj < significance) 
print(dim(res.pair.sig.mind)) 
res.sig.dim.mind     <- dim(res.pair.sig.mind)
res.pair.sig.mind    <- res.pair.sig.mind[order(-res.pair.sig.mind$log2FoldChang), ]


write.table(res.pair.sig, file =  paste0(Project, "-DESeq2_condition_pair_res_l2fc1_padj05_N",res.sig.dim[1],"_nofilter.txt"), row.names=F)
write.table(res.pair.sig.mind, file =  paste0(Project, "-DESeq2_condition_pair_res_l2fc1_padj05_N",res.sig.dim[1],"_rmDup_rmName_nofilter.txt"), row.names=F)


message("+----------    DESeq2 Analaysis heatmap plot                       ------------+")

## Heatmap data prepare and plot functions
Heatdat.fn  <- function(resmat, selcols, topN){
  rvmat              <- resmat[,selcols]
  #rvmat              <- subset(rvmat, external_gene_name!="" & !is.na(external_gene_name))
  #rvmat              <- rvmat[-which(duplicated(rvmat[,1])==T),]
  rownames(rvmat)    <- rvmat[,1]
  rvmat              <- rvmat[,-1]
  rv                 <- rowVars(as.matrix(rvmat))
  select             <- order(rv, decreasing = TRUE)[seq_len(min(topN, length(rv)))]
  heatvsd            <- rvmat[select,]
  heatvsd
}

Heatmap_Top_fn  <- function(heatmat, coldeg, gwhclust, ord.col, ha){
  pht_heat = Heatmap(as.matrix(heatmat), col = coldeg, 
                     name = "Tg vs Control", show_row_names=T,
                     show_column_names = T, width = unit(6, "cm"),
                     heatmap_legend_param =list(title = "Expression \n(Normalised Counts)",
                                                title_position = "topleft",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     top_annotation = ha,
                     cluster_rows = T, show_row_dend = T,
                     cluster_columns = reorder(as.dendrogram(gwhclust[[1]]), ord.col), show_column_dend = T,
                     column_title=NULL,
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 6),
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 6, col = rep(c("blue","red"), each=5)),
                     column_labels = unlist(lapply(colnames(heatmat), function(x) strsplit(x, split="-")[[1]][2]))
  )
  pht_heat
}

message("+----------------  Heatmap with top 50, 100 and 200 most variable genes ------------------+")

model     <- "TgvsControl"
condition <- rep(c("Control", "Tg"), each=5)
library(seriation)
ha  = HeatmapAnnotation(Condition = condition,  col = list(Condition = c("Control" = "blue", "Tg"= "red")), annotation_name_side = "left")

heat.matrix.50          <- Heatdat.fn(res.matMer.mind, c(2,10:19), 50)
o2.50                   <- seriate(dist(t(heat.matrix.50)), method = "GW")

heat.matrix.100         <- Heatdat.fn(res.matMer.mind, c(2,10:19), 100)
o2.100                  <- seriate(dist(t(heat.matrix.100)), method = "GW")

heat.matrix.200         <- Heatdat.fn(res.matMer.mind, c(2,10:19), 200)
o2.200                  <- seriate(dist(t(heat.matrix.200)), method = "GW")

breaksList.deg          <- seq(8, 20, by = 1)
ScaleCols.deg           <- colorRampPalette(colors = c("purple4", "white", "darkgreen"))(length(breaksList.deg))


heatmap.T50             <- Heatmap_Top_fn(heat.matrix.50, ScaleCols.deg, o2.50, 1:10, ha)
heatmap.T100            <- Heatmap_Top_fn(heat.matrix.100, ScaleCols.deg, o2.100, 1:10, ha)
heatmap.T200            <- Heatmap_Top_fn(heat.matrix.200, ScaleCols.deg, o2.200, 1:10, ha)


pdf(paste0(Project, "-", model, "_top50_Fig_nofilter.Heatmap.pdf"), width=6, height= 8)
draw(heatmap.T50, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_top100_Fig_nofilter.Heatmap.pdf"), width=6, height= 10)
draw(heatmap.T100, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_top200_Fig_nofilter.Heatmap.pdf"), width=6, height= 16)
draw(heatmap.T200, heatmap_legend_side = "right", merge_legend=T)
dev.off()

## Hierarchical clustering  
# install.packages("ggdendro")
library("ggdendro")
sample_distances <- dist(t(assay(vsd.pair)))
hclust_Fig       <- ggdendrogram(hclust(sample_distances), rotate = FALSE, segments = TRUE)
pdf(paste0(Project, "-Hierarchical_Clustering_nofilter_plot.pdf"))
print(hclust_Fig)
dev.off()

message("+----Heatmap with the top 100/200 genes and highlight the genes relating to the pathway GO:0034976---+")

load(paste0(Project, "-DESeq2_condition_pair_dds_vsd_res_nofilter.RData"))
respair.all      <- read.table("CTR_gjb2_0002-DESeq2_res_all_normvsd_N20701_summaryTable.txt", header=T)
respair.sig      <- read.table("CTR_gjb2_0002-DESeq2_condition_pair_res_l2fc1_padj05_N3194.txt", header=T)
respair.sig.new  <- respair.sig[,-c(10:19)]
normCounts.pair  <- as.data.frame(counts(dds.pair, normalized=TRUE))
normCounts.pair$ensembl_gene_id <- rownames(normCounts.pair) 
respair.all.newM <- merge(respair.all, normCounts.pair, by = "ensembl_gene_id")
respair.all.newM <- respair.all.newM[,c(1:9,15:19, 10:14)]
respair.all.newM[,10:19]  <- log2(respair.all.newM[,10:19]+1)
colnames(respair.all.newM)[10:19] <- colnames(respair.all)[c(15:19, 10:14)]
colnames(respair.all.newM)[10:19] <- gsub("[.]", "-", colnames(respair.all.newM)[10:19])

respair.sig.newM <- merge(respair.sig.new, normCounts.pair, by = "ensembl_gene_id")
respair.sig.newM <- respair.sig.newM[,c(1:9,10:14,15:19)]
colnames(respair.sig.newM)[10:19] <- colnames(respair.sig)[c(15:19, 10:14)]
colnames(respair.sig.newM)[10:19] <- gsub("[.]", "-", colnames(respair.sig.newM)[10:19])
min(respair.sig.newM[,10:19])
respair.sig.newM[,10:19]  <- log2(respair.sig.newM[,10:19]+1)
respair.sig.newM <- respair.sig.newM[order(respair.sig.newM$padj),]


write.table(respair.all.newM, file = "DESeq2_condition_pair_res_all_l2normCounts_N20701_Oct_2021.txt", row.names=F)
write.table(respair.sig.newM, file = "DESeq2_condition_pair_res_l2fc1_padj05_l2normCounts_N3194_Oct_2021.txt", row.names=F)

sel.respair.sig.newM <- subset(respair.sig.newM, external_gene_name!="") 
pair1.l2fc      <- sel.respair.sig.newM[,"Tg-A003"]-sel.respair.sig.newM[,"Control-A001"]
pair2.l2fc      <- sel.respair.sig.newM[,"Tg-A010"]-sel.respair.sig.newM[,"Control-A009"]
pair3.l2fc      <- sel.respair.sig.newM[,"Tg-A018"]-sel.respair.sig.newM[,"Control-A016"]
pair4.l2fc      <- sel.respair.sig.newM[,"Tg-A021"]-sel.respair.sig.newM[,"Control-A020"]
pair5.l2fc      <- sel.respair.sig.newM[,"Tg-A025"]-sel.respair.sig.newM[,"Control-A023"]

selTop100       <- cbind(pair1.l2fc[1:100],pair2.l2fc[1:100],pair3.l2fc[1:100],pair4.l2fc[1:100],pair5.l2fc[1:100])
selTop200       <- cbind(pair1.l2fc[1:200],pair2.l2fc[1:200],pair3.l2fc[1:200],pair4.l2fc[1:200],pair5.l2fc[1:200])
rownames(selTop100) <- sel.respair.sig.newM[1:100, "external_gene_name"]
rownames(selTop200) <- sel.respair.sig.newM[1:200, "external_gene_name"]
colnames(selTop100) <- c("A003/A001", "A010/A009", "A018/A016", "A021/A020", "A025/A023")
colnames(selTop200) <- c("A003/A001", "A010/A009", "A018/A016", "A021/A020", "A025/A023")


l2fcTop100      <- as.matrix(sel.respair.sig.newM[1:100,"log2FoldChange"], ncol=1)
colnames(l2fcTop100) <- "L2FC"
rownames(l2fcTop100) <- rownames(selTop100)

l2fcTop200      <- as.matrix(sel.respair.sig.newM[1:200,"log2FoldChange"], ncol=1)
colnames(l2fcTop200) <- "L2FC"
rownames(l2fcTop200) <- rownames(selTop200)

library(seriation)
library(circlize)
o2Top100                   <- seriate(dist(t(selTop100)), method = "GW")
o2Top200                   <- seriate(dist(t(selTop200)), method = "GW")

### Heatmap function version_0 is used to plot by ranking the sigGenes with smallest p-values. 
### and plot pair l2fc, highlight selected GO genes
Heatmap_Top_fn_v0  <- function(heatmat, col.deg1, gwhclust, ord.col, selGO.GN){
  pht_heat = Heatmap(as.matrix(heatmat), col = col.deg1,  
                     name = "Heatmap", show_row_names=T,
                     show_column_names = T, width = unit(6, "cm"),
                     heatmap_legend_param =list(title = "Log2FoldChange",
                                                title_position = "leftcenter-rot",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     cluster_rows = T, 
                     cluster_columns = reorder(as.dendrogram(gwhclust[[1]]), ord.col), 
                     show_column_dend = T,
                     column_title=NULL,
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side="left",
                     show_row_dend=T,
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 8, fontface="bold")) +
    rowAnnotation(link = anno_mark(at = which(rownames(heatmat)%in%selGO.GN==T), 
                                   labels = rownames(heatmat)[which(rownames(heatmat)%in%selGO.GN==T)], 
                                   labels_gp = gpar(fontsize = 8, col="red"), padding = unit(1, "mm")))
  
  
  
  return(pht_heat)
}
modelv0 <- "Pair-l2fc"
cols.fn <- colorRamp2(c(-8, 0, 8), c("purple4", "white", "darkgreen"))
load("CTR_gjb2_0002-Ensembl_merge_resall_GOID.RData")
selGO.GID  <- ensembl_nameID[ensembl_nameID[,"go_id"]=="GO:0034976","ensembl_gene_id"]
selGO.GN   <- respair.sig.newM[respair.sig.newM$ensembl_gene_id%in%unique(selGO.GID), "external_gene_name"] ## 24
selTop100.hplt      <- Heatmap_Top_fn_v0(selTop100, cols.fn, o2Top100, 1:5, selGO.GN)
selTop200.hplt      <- Heatmap_Top_fn_v0(selTop200, cols.fn, o2Top200, 1:5, selGO.GN)

pdf(paste0(Project, "-", modelv0, "_top100_Fig_Heatmap.pdf"), width=6, height= 8)
draw(selTop100.hplt, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", modelv0, "_top200_Fig_Heatmap.pdf"), width=6, height= 16)
draw(selTop200.hplt, heatmap_legend_side = "right", merge_legend=T)
dev.off()

### Heatmap v1, plot log2NormCounts across all samples, and highlight selected GO genes. 
### plus a column with log2FC between Tg and Ctrl

Heatmap_Top_fn_v1  <- function(heatmat, col.deg1, gwhclust, col.deg2, ord.col, basemeans.plt, selGO.GN){
  pht_heat = Heatmap(as.matrix(heatmat), col = col.deg1,  
                     name = "Heatmap", show_row_names=T,
                     show_column_names = T, width = unit(6, "cm"),
                     heatmap_legend_param =list(title = "Log2(Norm Counts+1)",
                                                title_position = "leftcenter-rot",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     top_annotation = ha,
                     cluster_rows = T, 
                     cluster_columns = reorder(as.dendrogram(gwhclust[[1]]), ord.col), 
                     show_column_dend = T,
                     column_title=NULL,
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side="left",
                     show_row_dend=T,
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 8, col = rep(c("blue","red"), each=5)),
                     column_labels = unlist(lapply(colnames(heatmat), function(x) strsplit(x, split="-")[[1]][2])))+
    Heatmap(basemeans.plt, name = "l2FC", col=col.deg2, width = unit(5, "mm"),column_names_gp = gpar(fontsize = 8, fontface="bold"),
            heatmap_legend_param = list(title = "Log2FoldChange")) +
    rowAnnotation(link = anno_mark(at = which(rownames(heatmat)%in%selGO.GN==T), 
                                   labels = rownames(heatmat)[which(rownames(heatmat)%in%selGO.GN==T)], 
                                   labels_gp = gpar(fontsize = 8, col="red"), padding = unit(1, "mm")))
  
  return(pht_heat)
}
modelv1      <- "AllSample-l2Norm"
breaksListv1 <- seq(4, 20, by = 1)
cols.fnv1    <- colorRampPalette(colors = c("lightcyan", "darkblue"))(length(breaksListv1))
###
condition <- rep(c("Control", "Tg"), each=5)
ha  = HeatmapAnnotation(Condition = condition,  col = list(Condition = c("Control" = "blue", "Tg"= "red")), annotation_name_side = "left")
selTop100.v1 <- sel.respair.sig.newM[1:100,c(15:19,10:14)]
rownames(selTop100.v1) <- sel.respair.sig.newM[1:100, "external_gene_name"]
selTop200.v1 <- sel.respair.sig.newM[1:200,c(15:19,10:14)]
rownames(selTop200.v1) <- sel.respair.sig.newM[1:200, "external_gene_name"]
###
o2Top100.v1     <- seriate(dist(t(selTop100.v1)), method = "GW")
o2Top200.v1     <- seriate(dist(t(selTop200.v1)), method = "GW")

selTop100.hplt1      <- Heatmap_Top_fn_v1(selTop100.v1, cols.fnv1, o2Top100.v1, cols.fn, 1:10, l2fcTop100, selGO.GN)
selTop200.hplt1      <- Heatmap_Top_fn_v1(selTop200.v1, cols.fnv1, o2Top200.v1, cols.fn, 1:10, l2fcTop200, selGO.GN)

## SFig2B
pdf("SFig2B-AllSample-l2Norm_top100_Fig_Heatmap.pdf", width=6, height= 10)
draw(selTop100.hplt1, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", modelv1, "_top200_Fig_Heatmap.pdf"), width=6, height= 18)
draw(selTop200.hplt1, heatmap_legend_side = "right", merge_legend=T)
dev.off()

#### Figure 2C: GO:00006486 with FoldChange >= 1.5
ensembl    =  useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host = 'ensembl.org')
ensembl_nameID <- getBM(attributes=c('ensembl_gene_id', "go_id", "name_1006"),  mart = ensembl, useCache = FALSE) 

selGO.GID2C  <- ensembl_nameID[ensembl_nameID[,"go_id"]=="GO:0006486","ensembl_gene_id"]
selGO.GN2C   <- respair.all.newM[respair.all.newM$ensembl_gene_id%in%unique(selGO.GID2C), "external_gene_name"]  ## 146
sel.respair.all.newM <- subset(respair.all.newM, external_gene_name!="") ## 17297
Fig2C.mat    <- subset(sel.respair.all.newM, padj < 0.05 & abs(log2FoldChange) >= log2(1.5))
Fig2C.mat    <- Fig2C.mat[Fig2C.mat$external_gene_name%in%selGO.GN2C==T,] ## 66

write.table(Fig2C.mat,file="GO6486_FC1.5_N66_Summary_Table_Oct_2021.txt", row.names=F)

Fig2C.plt    <- Fig2C.mat[,10:19]
rownames(Fig2C.plt) <- Fig2C.mat$external_gene_name
Fig2C.l2fc   <- as.matrix(Fig2C.mat$log2FoldChange, ncol=1)
colnames(Fig2C.l2fc) <- "L2FC"
rownames(Fig2C.l2fc) <- Fig2C.mat$external_gene_name
breaksListv2 <- seq(3, 4, by = 0.5)
cols.fnv2    <- colorRampPalette(colors = c("lightcyan", "darkblue"))(length(breaksListv2))
cols.fnl2    <- colorRamp2(c(-4, 0, 4), c("purple4", "white", "darkgreen"))

## swipe the control and Tg order 
Fig2C.plt <- Fig2C.plt[,c(6:10,1:5)]
o2Top100.v2  <- seriate(dist(t(Fig2C.plt )), method = "GW")

### Heatmap similar as version 1, except no highlight genes rowannotation.
Heatmap_Top_fn_v2  <- function(heatmat, col.deg1, gwhclust, col.deg2, ord.col, basemeans.plt){
  pht_heat = Heatmap(as.matrix(heatmat), col = col.deg1,  
                     name = "Heatmap", show_row_names=T,
                     show_column_names = T, width = unit(6, "cm"),
                     heatmap_legend_param =list(title = "Log2(Norm Counts+1)",
                                                title_position = "leftcenter-rot",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     top_annotation = ha,
                     cluster_rows = T, 
                     cluster_columns = reorder(as.dendrogram(gwhclust[[1]]), ord.col), 
                     show_column_dend = T,
                     column_title=NULL,
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side="right",
                     show_row_dend=T,
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 8, col = rep(c("blue","red"), each=5)),
                     column_labels = unlist(lapply(colnames(heatmat), function(x) strsplit(x, split="-")[[1]][2])))+
    Heatmap(basemeans.plt, name = "l2FC", col=col.deg2, width = unit(5, "mm"),column_names_gp = gpar(fontsize = 8, fontface="bold"),
            heatmap_legend_param = list(title = "Log2FoldChange")) 
  
  return(pht_heat)
}
modelv2 <- "Fig2C-GO0006486"
selTop100.hplt2C   <- Heatmap_Top_fn_v2(Fig2C.plt, cols.fnv2, o2Top100.v2, cols.fnl2, 1:10, Fig2C.l2fc)

pdf("Fig2C-GO0006486_top66_FC1.5_Fig_Heatmap.pdf", width=6, height= 10)
draw(selTop100.hplt2C, heatmap_legend_side = "right", merge_legend=T)
dev.off()


message("+---                  Gene Ontology Analysis          ---------------+")

message("+-------------------------------------------------------------------------------+")
message("+    GeneOntology Analysis                                                      +")
message("+-------------------------------------------------------------------------------+")

suppressPackageStartupMessages({
  library("clusterProfiler")
  library("DOSE")
  library("GSEABase")
  library("AnnotationHub")
  library("org.Hs.eg.db")
  library("BiocParallel")
  library("gage")
  library("gageData")
  library("enrichplot")
  library("reshape2")
  library("ggraph")
  library("ggforce")
  library("reactome.db")
  library("ReactomePA")
  library("Biostrings")
})

hub <- AnnotationHub()
query(hub, "Homo_sapiens")
register(MulticoreParam(2))
ensEMBL2id.GO <- getBM(attributes=c('ensembl_gene_id', "entrezgene_id"), mart = ensembl) 

message("+---              Data sorting for clusterProfiler input format --------------+")
message("+---- The total number of input genes and bkgenes will slightly differ from the output GeneRatio and BKgeneRatio                  ------------------------------------+")
message("+-- That is due to 1. If your input gene id contains duplicated IDs, those duplicated will be removed. 2. Those genes that do not have GO annotation will be removed. +")

## significant DEGs
sel.column.names <-c("entrezgene_id", "ensembl_gene_id", "external_gene_name", "log2FoldChange")
res.pair.sig     <- subset(res.allmat, padj < 0.05 & abs(log2FoldChange)>=1)
resdf.ann        <- merge(res.pair.sig, ensEMBL2id.GO, by = "ensembl_gene_id", all.x=T)
resdf            <- resdf.ann[,colnames(resdf.ann)%in%sel.column.names==T]
resdf            <- resdf[order(-resdf$log2FoldChange),sel.column.names]
resdf            <- resdf[-which(duplicated(resdf$ensembl_gene_id)==T),]
colnames(resdf)  <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC")  ## 3194
print(dim(resdf))  ## 3194
print(length(unique(resdf$ENSEMBL))) ## 3194
print(length(unique(resdf$SYMBOL))) ## 2843, 
print(length(unique(resdf$ENTREZID))) ## 2636

## Background genes list
background_genes <- merge(res.allmat, ensEMBL2id.GO, by ="ensembl_gene_id", all.x=T) ## 20817
bkgdf            <- background_genes[,colnames(background_genes)%in%sel.column.names==T]
bkgdf            <- bkgdf[order(-bkgdf$log2FoldChange),sel.column.names]
bkgdf            <- bkgdf[-which(duplicated(bkgdf$ensembl_gene_id)==T),]
colnames(bkgdf)  <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC") 
print(dim(bkgdf))  ## 20701
print(length(unique(bkgdf$ENSEMBL))) ## 20701
print(length(unique(bkgdf$SYMBOL))) ## 17277
print(length(unique(bkgdf$ENTREZID))) ## 15077


message("+------------------ enrichR analysis  (revised 29/09/2021)          ----------------------------+")

#install.packages("enrichR")
library(enrichR)
dbs             <- c("KEGG_2021_Human", "Reactome_2016",  
                     "MSigDB_Hallmark_2020","GO_Biological_Process_2018",
                     "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
enrichR_res     <- enrichr(resdf$SYMBOL, databases = dbs)
sdims           <- lapply(enrichR_res, dim)

enrichR_res.sig <- lapply(enrichR_res, function(x) subset(x, Adjusted.P.value < 0.05))
sigdims         <- lapply(enrichR_res.sig, dim)
for(i in 1:length(dbs)){
  enrichR_res.sig[[i]]$pathways <- as.character(dbs[i])
  enrichR_res.sig
}
save(enrichR_res.sig, enrichR_res, file = paste0(Project, "-enrichR_all_padj05_summary_nofilter_2021.RData"))


library(plyr)
enrichR_res.sig.sel <- list(enrichR_res.sig[[1]], enrichR_res.sig[[2]], enrichR_res.sig[[3]],enrichR_res.sig[[4]],
                            enrichR_res.sig[[5]],enrichR_res.sig[[6]])
merRpathways        <- ldply(enrichR_res.sig.sel, data.frame)
merRpathways$paths  <- rep(c("KEGG","Reactome", "Hallmark", "BP", "CC", "MF"),
                           c(unlist(lapply(enrichR_res.sig.sel,function(x) dim(x)[1]))))
table(merRpathways$paths)

enrichR_res.sel     <- list(enrichR_res[[1]], enrichR_res[[2]], enrichR_res[[3]],enrichR_res[[4]],
                            enrichR_res[[5]],enrichR_res[[6]])
merRpathways.all    <- ldply(enrichR_res.sel, data.frame)
merRpathways.all$paths  <- rep(c("KEGG","Reactome", "Hallmark", "BP", "CC", "MF"),
                               c(unlist(lapply(enrichR_res.sel,function(x) dim(x)[1]))))
table(merRpathways.all$paths)

write.csv(merRpathways, file = paste0(Project, "-enrichr_kegg6_reactome16_hallmark26_BP6_CC2_MF3_sigPadj_summary_Sep_2021.csv"), row.names=F)
write.csv(merRpathways.all, file = paste0(Project, "-enrichr_kegg312_reactome1336_hallmark50_BP4426_CC374_MF976_all_summary_Sep_2021.csv"), row.names=F)

message("+----    WikiPathways analysis, entrezID -----------------------------------------------+")

library(magrittr)
geneList         <- resdf$L2FC
names(geneList)  <- resdf$ENTREZID
geneList         <- geneList[-which(is.na(names(geneList)))]
geneList         <- geneList[-which(duplicated(names(geneList))==T)]
genes            <- names(geneList) ## 2634

## background gene list
bkgeneList       <- bkgdf$L2FC
names(bkgeneList)<- bkgdf$ENTREZID
bkgeneList       <- bkgeneList[-which(is.na(names(bkgeneList)))]
bkgeneList       <- bkgeneList[-which(duplicated(names(bkgeneList))==T)]
bkgenes          <- names(bkgeneList) ## 15076


wpgmtfile        <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene          <- read.gmt(wpgmtfile)
wp2gene          <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene        <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name        <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp              <- enricher(genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, universe = bkgenes)
head(ewp) ## Gene856, Bg4436

ewp2             <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
head(ewp2)

library(org.Hs.eg.db)
ewp              <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
ewp2             <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
head(ewp); head(ewp2)

save(ewp, ewp2, file = paste0(Project, "-WikiPathways_enrichr_gsea.RData"))

write.csv(ewp, file = paste0(Project, "-WikiPathways_enrichr_N9_nofilter_OG2612_OBG15326_WG856_WBG4436.csv"), row.names=F)
write.csv(ewp2, file = paste0(Project, "-WikiPathways_gsea_N4_nofilter_OG2612_OBG15326_WG856_WBG4436.csv"), row.names=F)

message("+---           Cell Markers                                            -------+")

cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))

y           <- enricher(genes, TERM2GENE=cell_markers, minGSSize=1, universe = bkgenes)
y           <- setReadable(y, org.Hs.eg.db, keyType = "ENTREZID")

save(y, file = paste0(Project, "-cellMarkers_enrichr_nofilter.RData"))

write.csv(y, file = paste0(Project, "-cellMarkers_enrichr_N9_nofilter.csv"), row.names=F)

message("+----MSigDb analysis (Molecular Signatures Database)--------------------------+")
message("+H: hallmark gene sets C1: positional gene sets C2: curated gene sets --------+")
message("+C3: motif gene sets C4: computational gene sets C5: GO gene sets     --------+")
message("+        C6: oncogenic signatures C7: immunologic signatures                  +")

library(msigdbr)
msigdbr_species()
m_df        <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

m_t2g       <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em          <- enricher(genes, TERM2GENE=m_t2g)
em          <- setReadable(em, org.Hs.eg.db, keyType = "ENTREZID")
em2         <- GSEA(geneList, TERM2GENE = m_t2g)
em2         <- setReadable(em2, org.Hs.eg.db, keyType = "ENTREZID")

m_t2g.new   <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em3         <- GSEA(geneList, TERM2GENE = m_t2g.new)
em3         <- setReadable(em3, org.Hs.eg.db, keyType = "ENTREZID")

save(em,em2, em3, file = paste0(Project, "-MSigDb_enrihcer_gsea_nofilter.RData"))
write.csv(em2, file = paste0(Project, "-cMSigDb_gsea_N19_BP18_MF1_nofilter.csv"), row.names=F)

message("+--- GO over-representation test                             --------------+")
## "GO:0006487"| "GO:0018279"| "GO:0006486"
selIDs           <- c("GO:0006487", "GO:0018279", "GO:0006486")
ego              <- enrichGO( gene          = genes,
                              OrgDb         = org.Hs.eg.db,
                              universe      = bkgenes,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.05,
                              keyType       = "ENTREZID",
                              readable      = TRUE)
## "protein glycosylation" "GO:0006486"
ego.kegg         <- enrichKEGG( gene        = genes,
                                organism      = "hsa",
                                universe      = bkgenes,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                keyType       = "kegg")

write.csv(ego, file=paste0(Project, "-GeneOntology_BP158_CC19_MF30_summaryTable.csv"), row.names=F)
write.csv(ego.kegg, file=paste0(Project, "-GeneOntology_KEGG6_summaryTable.csv"), row.names=F)


message("+---           gseGO                                            --------------+")

ego.gse          <- gseGO( geneList     = geneList,
                           OrgDb        = org.Hs.eg.db,
                           ont          = "ALL",
                           keyType      = "ENTREZID",
                           minGSSize    = 10,
                           maxGSSize    = 500,
                           pvalueCutoff = 0.05,
                           verbose      = FALSE,
                           pAdjustMethod = "BH",
                           by           = "fgsea")
ego.gse         <- setReadable(ego.gse, org.Hs.eg.db, keyType = "ENTREZID")


save(ego.gse, file="CTR_gjb2_0002-clusterProfiler_gseGO_EntrezID_N83.RData")
write.csv(as.data.frame(ego.gse), file ="CTR_gjb2_0002-clusterProfiler_gseGO_EntrezID_N83.csv")

kegg.gse        <- gseKEGG(geneList    = geneList,
                           organism     = "hsa",
                           keyType      = "kegg",
                           # nPerm        = 1000,
                           minGSSize    = 10,
                           maxGSSize    = 500,
                           pvalueCutoff = 0.05,
                           verbose      = FALSE,
                           pAdjustMethod = "BH"
)
kegg.gse         <- setReadable(kegg.gse, org.Hs.eg.db, keyType = "ENTREZID")
save(ego, ego.kegg, ego.gse, kegg.gse, file=paste0(Project, "-GOenrich_GOGSE_KEGG.RData"))

write.csv(ego.gse, file=paste0(Project, "-GSE_BP47_CC2_MF1_summaryTable.csv"), row.names=F)
write.csv(kegg.gse, file=paste0(Project, "-GSE_KEGG4_summaryTable.csv"), row.names=F)



message("+--   semantic similarity among GO terms  rrvgo package (1.2.0), esp BP       ---------------------------+")
message("+Reduce and visualize lists of Gene Ontology terms by identifying redudance based on semantic similarity.+")

library("rrvgo")
library("GOSemSim")
ego           <- read.csv(paste0(Project, "-GeneOntology_BP158_CC19_MF30_summaryTable.csv"),header=T)
ego.dat       <- as.data.frame(ego)
ego.BPdat     <- subset(ego.dat, ONTOLOGY=="BP")

GOSeSimdata.BP   <- godata(OrgDb = "org.Hs.eg.db", keytype = "ENTREZID", ont="BP", computeIC = TRUE)
simMatrix.BP     <- calculateSimMatrix(ego.BPdat$ID, orgdb="org.Hs.eg.db", ont="BP",
                                       semdata=GOSeSimdata.BP, method="Rel")
scores.BP        <- setNames(-log10(ego.BPdat$qvalue), ego.BPdat$ID)
reducedTerms.BP  <- reduceSimMatrix(simMatrix.BP, scores.BP, threshold=0.9, orgdb="org.Hs.eg.db")

## Removed 1 terms that were not found in orgdb for BP
write.csv(reducedTerms.BP, file = paste0(Project, "-GOSeSim_BP157_cluster14_summaryTable.csv"), row.names=F)


options(ggrepel.max.overlaps = Inf)

pdf(paste0(Project,"-GOSeSim_BP_rrvgo_scatterplot.pdf"))
scatterPlot(simMatrix.BP, reducedTerms.BP, size="score", labelSize = 2)
dev.off()

pdf(paste0(Project,"-GOSeSim_BP_rrvgo_heatmap.pdf"), height= 16, width= 14)
heatmapPlot(simMatrix.BP,reducedTerms.BP, annotateParent=TRUE,
            annotationLabel="parentTerm", fontsize=6) 

dev.off()

message("+--MSigDb GSEA plot for selected pathways in version 10- Fig2A-------------------+")

Selgsea.em2       <- c("GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE",
                       "GOBP_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",
                       "GOBP_ERAD_PATHWAY",
                       "GOBP_ENDOPLASMIC_RETICULUM_MANNOSE_TRIMMING")

Selgsea.em2ID     <- unlist(lapply(Selgsea.em2, function(x) which(em2$ID==x)))

pdf(paste0(Project,"-MSigDb_GSEA_plot_sel2BP_Fig2A.pdf"), width=14, height=5)
gseaplot2(em2, geneSetID = c(3,6), pvalue_table = TRUE)
dev.off()

## check the number of genes direction in the above two selected pathways,
## the genes are all upregulated in Tg treatment group 

em2.selGList      <- em2[Selgsea.em2ID, "core_enrichment"]
em2.selGenes      <- lapply(em2.selGList, function(x) strsplit(x, split="/"))
## DERL2 is the only different genes for endoplasmic reticulum unfolded protein response (37 genes) to cellular (42 genes)
## "SDF2L1" "OPTN"   "HSPA13" "UGGT1"  "RNF185" are the extra genes for cellular compare to endoplasmic.

pdf(paste0(Project,"-MSigDb_GSEA_sel2BP_cnetPlot.pdf"))
cnetSelBP         <- cnetplot(em2, categorySize="pvalue", foldChange=geneList, showCategory = Selgsea.em2[1:2],
                              cex_label_category = 0.5, cex_label_gene = 0.6)
print(cnetSelBP)
dev.off()


message("+---------------Glycosylation Heatmap plot         ----------------------------+")

## "protein glycosylation" "GO:0006486";
## "protein N-linked glycosylation"GO:0006487"
## "protein N-linked glycosylation via asparagine"GO:0018279"
ensEMBL2id.pro    <- getBM(attributes=c('external_gene_name', 'go_id'),  mart = ensembl, useCache = FALSE)   
res.matMer.GO     <- merge(res.matMer.mind, ensEMBL2id.pro, by = "external_gene_name")
res.matMer.selGO  <- lapply(selIDs, function(x) res.matMer.GO[res.matMer.GO$go_id%in%x==T,])
res.matMer.selGO.genes <- lapply(res.matMer.selGO, function(x) x$external_gene_name)
names(res.matMer.selGO.genes) <- unique(c(res.matMer.selGO[[1]]$go_id, res.matMer.selGO[[2]]$go_id,res.matMer.selGO[[3]]$go_id))
## devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")
pdf(paste0(Project,"-GO_glycosylation_overlaps_vennDiagram.pdf"))

ggvenn(
  res.matMer.selGO.genes, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

## enrichGO we identified "protein glycosylation" "GO:0006486";
GlyHeat.genes         <- unlist(strsplit(ego[which(ego$ID=="GO:0006486"),"geneID"], split="/"))
GlyHeat.mat           <- res.pair.sig[res.pair.sig$external_gene_name%in%GlyHeat.genes, c(2,10:19)]
rownames(GlyHeat.mat) <- GlyHeat.mat$external_gene_name
GlyHeat.mat           <- GlyHeat.mat[,-1]
cname                 <- ifelse(colData(vsd)$condition=="control", "Control", "Tg")
colnames(GlyHeat.mat) <- paste0(cname, "_", colData(vsd)$Individual)
colnames(GlyHeat.mat) <- as.factor(colnames(GlyHeat.mat))
GlyHeat.mat           <- GlyHeat.mat[,sort(colnames(GlyHeat.mat))]

pdf(paste0(Project,"-enrichGO_proteinglycosylation_N47_Heatmap.pdf"))
o2.gly    <- seriate(dist(t(GlyHeat.mat)), method = "GW")
Tg.mpht   <- Heatmap_Top_fn(GlyHeat.mat, ScaleCols.deg, o2.gly, 10:1, ha)
draw(Tg.mpht, heatmap_legend_side = "right", merge_legend=T)
dev.off()

message("+---- GSA plus GSEA analysis and plots for supplementary Figures ---------------------------------------+")
## use fgsea and gage to get GSEA analysis output

GSEA = function(gene_list, GO_file, pval, GOName) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=10,
                        maxSize=600,
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(10,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  print(dim(rbind(ups,downs)))
  ## Define up / down pathways which are significant in both tests
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( keepups$pathway, keepdowns$pathway))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  #sel.num = ifelse(dim(fgRes)[1] < 10, dim(fgRes)[1], 10)
  #filtRes = rbind(head(fgRes, n = 10),
  #                tail(fgRes, n = 10 ))
  filtRes = fgRes
  
  
  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtRes$Index = as.factor(1:nrow(filtRes))
  
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0("GSEA - ", GOName)) + 
    theme_minimal() + theme_bw() +
    theme(legend.position = "none")
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}

gene_list <- resdf[,4]
names(gene_list) <- resdf[,3]

gseaHallmark <- GSEA(gene_list, "h.all.v7.4.symbols.gmt", 0.05, "HALLMARK")
gseaReactome <- GSEA(gene_list, "c2.cp.reactome.v7.4.symbols.gmt", 0.05, "REACTOME")
gseaKegg     <- GSEA(gene_list, "c2.cp.kegg.v7.4.symbols.gmt", 0.05, "KEGG")
gseaBP       <- GSEA(gene_list, "c5.go.bp.v7.4.symbols.gmt", 0.05, "GO_BP")
filt_p       <- c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
                  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
                  "GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE")
GO.Files     <- c("h.all.v7.4.symbols.gmt", "c2.cp.reactome.v7.4.symbols.gmt", "c5.go.bp.v7.4.symbols.gmt")
selGSEA      <- list()
for(glib in 1:length(GO.Files)){
  
  myGO  = fgsea::gmtPathways(GO.Files[glib])
  selGSEA[[glib]] <- myGO[names(myGO)%in%filt_p==T]
  selGSEA
}
pdf(paste0("GSEA-", filt_p[1],".pdf"),height=5,width=7) # change height and width parameter
mytbl <- t(gseaHallmark$Results[grep("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", gseaHallmark$Results$pathway), c("padj", "NES", "size")])
mytbl <- cbind(rownames(mytbl), mytbl)
mytbl <- data.frame(Stats=mytbl[,1], Values=as.numeric(mytbl[,2]))

mytbl[1,2] <- round(mytbl[1,2], digits = 4)
mytbl[2,2] <- round(mytbl[2,2], digits = 4)
mytbl[3,2] <- round(mytbl[3,2], digits = 1)

plt <- plotEnrichment(pathway = selGSEA[[1]]$HALLMARK_UNFOLDED_PROTEIN_RESPONSE, gene_list) + 
  labs(title=filt_p[1]) + theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  annotation_custom(tableGrob(mytbl, rows=NULL), 
                    xmin=2500, xmax=3000, ymin=0.3, ymax=0.4)
print(plt)
dev.off()

pdf(paste0("GSEA-", filt_p[2],".pdf"),height=5,width=7) # change height and width parameter

mytbl <- t(gseaReactome$Results[grep("REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR", gseaReactome$Results$pathway), c("padj", "NES", "size")])
mytbl <- cbind(rownames(mytbl), mytbl)
mytbl <- data.frame(Stats=mytbl[,1], Values=as.numeric(mytbl[,2]))

mytbl[1,2] <- round(mytbl[1,2], digits = 4)
mytbl[2,2] <- round(mytbl[2,2], digits = 4)
mytbl[3,2] <- round(mytbl[3,2], digits = 1)
plt <- plotEnrichment(pathway = selGSEA[[2]]$REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR, gene_list) + 
  labs(title=filt_p[2]) + theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  annotation_custom(tableGrob(mytbl, rows=NULL), 
                    xmin=2500, xmax=3000, ymin=0.3, ymax=0.4)
print(plt)
dev.off()

pdf(paste0("GSEA-", filt_p[3],".pdf"),height=5,width=7) # change height and width parameter
mytbl <- t(gseaBP$Results[grep("GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS", gseaBP$Results$pathway), c("padj", "NES", "size")])
mytbl <- cbind(rownames(mytbl), mytbl)
mytbl <- data.frame(Stats=mytbl[,1], Values=as.numeric(mytbl[,2]))

mytbl[1,2] <- round(mytbl[1,2], digits = 4)
mytbl[2,2] <- round(mytbl[2,2], digits = 4)
mytbl[3,2] <- round(mytbl[3,2], digits = 1)
plt <- plotEnrichment(pathway = selGSEA[[3]]$GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS, gene_list) + 
  labs(title=filt_p[3]) + theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  annotation_custom(tableGrob(mytbl, rows=NULL), 
                    xmin=2500, xmax=3000, ymin=0.3, ymax=0.4)
print(plt)
dev.off()

pdf(paste0("GSEA-", filt_p[4],".pdf"),height=5,width=7) # change height and width parameter
mytbl <- t(gseaBP$Results[grep("GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE", gseaBP$Results$pathway), c("padj", "NES", "size")])
mytbl <- cbind(rownames(mytbl), mytbl)
mytbl <- data.frame(Stats=mytbl[,1], Values=as.numeric(mytbl[,2]))

mytbl[1,2] <- round(mytbl[1,2], digits = 4)
mytbl[2,2] <- round(mytbl[2,2], digits = 4)
mytbl[3,2] <- round(mytbl[3,2], digits = 1)
plt <- plotEnrichment(pathway = selGSEA[[3]]$GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE, gene_list) + 
  labs(title=filt_p[4]) + theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  annotation_custom(tableGrob(mytbl, rows=NULL), 
                    xmin=2500, xmax=3000, ymin=0.3, ymax=0.4)
print(plt)
dev.off()

gseaHallmark_dat <- gseaHallmark$Results
gseaReactome_dat <- gseaReactome$Results
gseaKegg_dat     <- gseaKegg$Results
gseaBP_dat       <- gseaBP$Results

gsea_all         <- rbind(gseaHallmark_dat, gseaReactome_dat, gseaKegg_dat, gseaBP_dat)
gsea_all         <- apply(gsea_all,2,as.character)
write.csv(gsea_all, file= "GSEA_Hallmark_Reactome_Kegg_BP.csv", quote=T)

save(gseaHallmark, gseaReactome, gseaKegg, gseaBP, file = "GSEA_Hallmark_Reactome_Kegg_BP.RData")







message("+------ GeneOntology Figures, 3 Supplementary heatmap ---------------")

S456  <- read.xlsx("S4_S5_and_S6.xlsx", sheet=1)
S4L   <- S456[1:118,1]
S5L   <- S456[1:28,2]
S6L   <- S456[1:16,3]

S4L.mat <- sel.respair.all.newM[sel.respair.all.newM$external_gene_name%in%S4L,]
S4L.plt <- S4L.mat[,c(15:19,10:14)]
rownames(S4L.plt) <- S4L.mat$external_gene_name
S4L.l2fc <- as.matrix(S4L.mat$log2FoldChange, ncol=1)
colnames(S4L.l2fc) <- "L2FC"
rownames(S4L.l2fc) <- S4L.mat$external_gene_name

S5L.mat <- sel.respair.all.newM[sel.respair.all.newM$external_gene_name%in%S5L,]
S5L.plt <- S5L.mat[,c(15:19,10:14)]
rownames(S5L.plt) <- S5L.mat$external_gene_name
S5L.l2fc <- as.matrix(S5L.mat$log2FoldChange, ncol=1)
colnames(S5L.l2fc) <- "L2FC"
rownames(S5L.l2fc) <- S5L.mat$external_gene_name

S6L.mat <- sel.respair.all.newM[sel.respair.all.newM$external_gene_name%in%S6L,]
S6L.plt <- S6L.mat[,c(15:19,10:14)]
rownames(S6L.plt) <- S6L.mat$external_gene_name
S6L.l2fc <- as.matrix(S6L.mat$log2FoldChange, ncol=1)
colnames(S6L.l2fc) <- "L2FC"
rownames(S6L.l2fc) <- S6L.mat$external_gene_name

breaksListv3 <- seq(3, 5, by = 0.5)
cols.fnv3    <- colorRampPalette(colors = c("lightcyan", "darkblue"))(length(breaksListv3))
cols.fnl3    <- colorRamp2(c(-4, 0, 6), c("purple4", "white", "darkgreen"))
o2Top100.vS4  <- seriate(dist(t(S4L.plt)), method = "GW")
o2Top100.vS5  <- seriate(dist(t(S5L.plt)), method = "GW")
o2Top100.vS6  <- seriate(dist(t(S6L.plt)), method = "GW")


S4L.hplt   <- Heatmap_Top_fn_v2(S4L.plt, cols.fnv3, o2Top100.vS4, cols.fnl3, 1:10, S4L.l2fc)
S5L.hplt   <- Heatmap_Top_fn_v2(S5L.plt, cols.fnv3, o2Top100.vS5, cols.fnl3, 1:10, S5L.l2fc)
S6L.hplt   <- Heatmap_Top_fn_v2(S6L.plt, cols.fnv3, o2Top100.vS4, cols.fnl3, 1:10, S6L.l2fc)


pdf( "-SFig3-1_GSEA1_Fig_Heatmap.pdf", width=6, height= 18)
draw(S4L.hplt, heatmap_legend_side = "right", merge_legend=T)
dev.off()
pdf("-SFig3-2_GSEA2_Fig_Heatmap.pdf", width=6, height=6)
draw(S5L.hplt, heatmap_legend_side = "right", merge_legend=T)
dev.off()
pdf("-SFig3-3_GSEA3_Fig_Heatmap.pdf", width=6, height= 6)
draw(S6L.hplt, heatmap_legend_side = "right", merge_legend=T)
dev.off()





message("+------------ Check the top 100/200 how many genes with ER stress -----------------+")

Heatdat.fn.1  <- function(res.mat, sel.cols, selGO.GN, topN){
  rvmat              <- as.matrix(res.mat[, sel.cols])
  rownames(rvmat)    <- res.mat$external_gene_name  
  rvmat.sel          <- rvmat[-which(rownames(rvmat)==""),]
  print(dim(rvmat.sel))
  rvmat.sel          <- rvmat.sel[!is.na(rownames(rvmat.sel)), ]
  print(dim(rvmat.sel))
  rv                 <- rowVars(rvmat.sel)
  select             <- order(rv, decreasing = TRUE)[seq_len(min(topN, length(rv)))]
  heatvsd            <- rvmat.sel[select,]
  heatvsd.overGO     <- heatvsd[rownames(heatvsd)%in%selGO.GN, ]
  print(dim(heatvsd.overGO))
  colnames(heatvsd)  <- gsub("[.]", "_", colnames(heatvsd))
  heatvsd            <- heatvsd[,sort(colnames(heatvsd))]
  heatvsd
}
sel.cols     <- c(10:19)
res.mat      <- res.allmat
topNs        <- c(100,200)
selGO.GN     <- res.selGO.GN

top100.heat  <-  Heatdat.fn.1(res.mat, sel.cols, selGO.GN, topNs[1])
top200.heat  <-  Heatdat.fn.1(res.mat, sel.cols, selGO.GN, topNs[2])

## Generate the heat map with highlight the above  12/100 and 13/200 GO:0034976
library(seriation)
condition <- rep(c("Control", "Tg"), each=5)
ha  = HeatmapAnnotation(Condition = condition,  col = list(Condition = c("Control" = "blue", "Tg"= "red")), annotation_name_side = "left")
o2.100                   <- seriate(dist(t(top100.heat)), method = "GW")
o2.200                   <- seriate(dist(t(top200.heat)), method = "GW")

basemeans1               <- res.allmat[res.allmat$external_gene_name%in%rownames(top100.heat),c("external_gene_name", "log2FoldChange")]
basemeans1               <- basemeans1[order(match(basemeans1$external_gene_name, rownames(top100.heat))), ]
basemeans1.plt           <- as.matrix(basemeans1[,-1])
rownames(basemeans1.plt) <- basemeans1$external_gene_name

basemeans2               <- res.allmat[res.allmat$external_gene_name%in%rownames(top200.heat),c("external_gene_name", "log2FoldChange")]
basemeans2               <- basemeans1[order(match(basemeans2$external_gene_name, rownames(top200.heat))), ]
basemeans2.plt           <- as.matrix(basemeans2[,-1])
rownames(basemeans2.plt) <- basemeans2$external_gene_name

colnames(basemeans1.plt) <- "l2FC"
colnames(basemeans2.plt) <- "l2FC"

save(ha, o2.100,o2.200,top100.heat,top200.heat,basemeans1.plt, basemeans2.plt, selGO.GN, file = "testHeatmap_Sep_2021.RData")

Heatmap_Top_fn1  <- function(heatmat, col.deg1, gwhclust, ord.col, col.deg2, basemeans.plt, ha, selGO.GN){
  pht_heat = Heatmap(as.matrix(heatmat), col = col.deg1,  
                     name = "Heatmap", show_row_names=T,
                     show_column_names = T, width = unit(6, "cm"),
                     heatmap_legend_param =list(title = "Expression \n(Normalised Counts)",
                                                title_position = "topleft",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     top_annotation = ha,
                     cluster_rows = T, 
                     cluster_columns = reorder(as.dendrogram(gwhclust[[1]]), ord.col), 
                     show_column_dend = T,
                     column_title=NULL,
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side="left",
                     show_row_dend=T,
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 6, col = rep(c("blue","red"), each=5)),
                     column_labels = unlist(lapply(colnames(heatmat), function(x) strsplit(x, split="_")[[1]][2])))+
    Heatmap(basemeans.plt, name = "l2FC", col=col.deg2, width = unit(5, "mm"),
            heatmap_legend_param = list(title = "Log2FoldChange")) +
    rowAnnotation(link = anno_mark(at = which(rownames(heatmat)%in%selGO.GN==T), 
                                   labels = rownames(heatmat)[which(rownames(heatmat)%in%selGO.GN==T)], 
                                   labels_gp = gpar(fontsize = 8, col="red"), padding = unit(1, "mm")))
  
  
  pht_heat
}

breaksList.deg           <- seq(8, 20, by = 1)
ScaleCols.deg            <- colorRampPalette(colors = c("skyblue", "dodgerblue4"))(length(breaksList.deg))
breaksList.l2            <- seq(-10,10, by=1)
ScaleCols.l2             <- colorRampPalette(colors = c("purple4", "white", "darkgreen"))(length(breaksList.l2))

heatmap.T100            <- Heatmap_Top_fn1(top100.heat, ScaleCols.deg, o2.100, 1:10, ScaleCols.l2, basemeans1.plt, ha, selGO.GN)
heatmap.T200            <- Heatmap_Top_fn1(top200.heat, ScaleCols.deg, o2.200, 1:10, ScaleCols.l2, basemeans2.plt, ha, selGO.GN)
model                   <- "pair"



pdf(paste0(Project, "-", model, "_top100_Fig.newHeatmap.pdf"), width=6, height= 8)
draw(heatmap.T100, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_top200_Fig.newHeatmap.pdf"), width=6, height= 16)
draw(heatmap.T200, heatmap_legend_side = "right", merge_legend=T)
dev.off()

message("+--- DESeq2 heatmap newest version Oct,2021, with 5 paired columns log2FoldChange   ----------+")
res.sig        <- read.table(paste0(Project, "-DESeq2_condition_pair_res_l2fc1_padj05_N3194.txt"), header=T)
res.sig.sel    <- subset(res.sig, external_gene_name!="")
res.sig.ord    <- res.sig.sel[order(res.sig.sel$padj), ] 
heat100        <- res.sig.ord[1:100,c(2,6)]
heat200        <- res.sig.ord[1:200,c(2,6)]

res.matMer.100  <- res.matMer[res.matMer$external_gene_name%in%heat100$external_gene_name,]
pair1.l2fca     <- res.matMer.100[,"Tg-A003"]-res.matMer.100[,"Control-A001"]
pair2.l2fca     <- res.matMer.100[,"Tg-A010"]-res.matMer.100[,"Control-A009"]
pair3.l2fca     <- res.matMer.100[,"Tg-A018"]-res.matMer.100[,"Control-A016"]
pair4.l2fca     <- res.matMer.100[,"Tg-A021"]-res.matMer.100[,"Control-A020"]
pair5.l2fca     <- res.matMer.100[,"Tg-A025"]-res.matMer.100[,"Control-A023"]
## 
res.matMer.200  <- res.matMer[res.matMer$external_gene_name%in%heat200$external_gene_name,]
pair1.l2fcb     <- res.matMer.200[,"Tg-A003"]-res.matMer.200[,"Control-A001"]
pair2.l2fcb     <- res.matMer.200[,"Tg-A010"]-res.matMer.200[,"Control-A009"]
pair3.l2fcb     <- res.matMer.200[,"Tg-A018"]-res.matMer.200[,"Control-A016"]
pair4.l2fcb     <- res.matMer.200[,"Tg-A021"]-res.matMer.200[,"Control-A020"]
pair5.l2fcb     <- res.matMer.200[,"Tg-A025"]-res.matMer.200[,"Control-A023"]

## heatmat matrix with log2FoldChange for 5 pairs
hplt.100        <- cbind(pair1.l2fca, pair2.l2fca, pair3.l2fca, pair4.l2fca, pair5.l2fca)
hplt.200        <- cbind(pair1.l2fcb, pair2.l2fcb, pair3.l2fcb, pair4.l2fcb, pair5.l2fcb)
rownames(hplt.100) <- res.matMer.100$external_gene_name
rownames(hplt.200) <- res.matMer.200$external_gene_name
colnames(hplt.100) <- c("A003/A001", "A010/A009", "A018/A016", "A021/A020", "A025/A023")
colnames(hplt.200) <- c("A003/A001", "A010/A009", "A018/A016", "A021/A020", "A025/A023")

breaksList.l2            <- seq(-5,5, by=1)
ScaleCols.l2             <- colorRampPalette(colors = c("purple4", "white", "darkgreen"))(length(breaksList.l2))

load("CTR_gjb2_0002-Ensembl_merge_resall_GOID.RData")
selGOID        <- "GO:0034976"
res.selGOID    <- subset(res.merENID, go_id == selGOID)
res.selGO.GN   <- unique(res.selGOID$external_gene_name) ## 81 genes involved
res.selGO.GID  <- unique(res.selGOID$ensembl_gene_id) ## 81

Heatmap_Top_fn_newL2  <- function(heatmat, col.deg, ord.col, selGO.GN){
  pht_heat = Heatmap(as.matrix(heatmat), col = col.deg,  
                     name = "Heatmap", show_row_names=T,
                     show_column_names = T, width = unit(6, "cm"),
                     heatmap_legend_param =list(title = "log2FoldChange",
                                                title_position = "lefttop-rot",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     cluster_rows = T, 
                     cluster_columns =F,
                     column_order = ord.col,
                     show_column_dend = T,
                     column_title=NULL,
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side="left",
                     show_row_dend=T,
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 8, fontface="bold"))+
    rowAnnotation(link = anno_mark(at = which(rownames(heatmat)%in%selGO.GN==T), 
                                   labels = rownames(heatmat)[which(rownames(heatmat)%in%selGO.GN==T)], 
                                   labels_gp = gpar(fontsize = 8, col="red"), padding = unit(1, "mm")))
  
  
  pht_heat
}

model                        <- "pair"
ord.col                      <- c(1:5)
heatmap.T100.l2pair          <- Heatmap_Top_fn_newL2(hplt.100, ScaleCols.l2 , ord.col, res.selGO.GN)
heatmap.T200.l2pair          <- Heatmap_Top_fn_newL2(hplt.200, ScaleCols.l2 , ord.col, res.selGO.GN)


pdf(paste0(Project, "-", model, "_top100_Fig.l2fc.newPairHeatmap.pdf"), width=6, height= 10)
draw(heatmap.T100.l2pair, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_top200_Fig.l2fc.newPairHeatmap.pdf"), width=6, height= 18)
draw(heatmap.T200.l2pair, heatmap_legend_side = "right", merge_legend=T)
dev.off()




message("+ Heatmap for all genes relating to GO:0006486 pathways.    ------------------- +")

gly.heatall       <- unique(res.selGOGly.dat[,c(2,6,9)])
gly.1.3           <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(1.3),2]) ## 85
gly.1.4           <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(1.4),2]) ## 76
gly.1.5           <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(1.5),2]) ## 66
gly.2.0           <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(2.0),2]) ## 35

gly.1.3.up        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange >= log2(1.3),2]) 
gly.1.3.dw        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange <= -log2(1.3),2]) 
gly.1.4.up        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange >= log2(1.4),2]) 
gly.1.4.dw        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange <= -log2(1.4),2]) 
gly.1.5.up        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange >= log2(1.5),2]) 
gly.1.5.dw        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange <= -log2(1.5),2]) 
gly.2.0.up        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange >= log2(2.0),2]) 
gly.2.0.dw        <- unique(res.selGOGly.dat[res.selGOGly.dat$padj<0.05& res.selGOGly.dat$log2FoldChange <= -log2(2.0),2]) 

gly.heat.plt      <- as.matrix(gly.heatall[,2], ncol=1)
rownames(gly.heat.plt) <- gly.heatall[,1]
colnames(gly.heat.plt) <- "log2FoldChange"

breaksList.l2            <- seq(-5,5, by=1)
ScaleCols.l2             <- colorRampPalette(colors = c("purple4", "white", "darkgreen"))(length(breaksList.l2))

Heatmap_Top_fn3  <- function(heatmat, col.deg, selGO.GN.up, selGO.GN.dw){
  pht_heat = Heatmap(as.matrix(heatmat), col = col.deg,  
                     name = "Heatmap", show_row_names=T,
                     show_column_names = F, width = unit(2, "cm"),
                     heatmap_legend_param =list(title = "log2FoldChange",
                                                title_position = "lefttop-rot",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10),
                                                legend_width = unit(0.8,"cm"),
                                                legend_height= unit(4, "cm")),
                     cluster_rows = T, 
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side="left",
                     show_row_dend=T)+
    rowAnnotation(link = anno_mark(at = which(rownames(heatmat)%in%selGO.GN.up==T), 
                                   labels = rownames(heatmat)[which(rownames(heatmat)%in%selGO.GN.up==T)], 
                                   labels_gp = gpar(fontsize = 8, col="red"), padding = unit(1, "mm")))+
    rowAnnotation(link = anno_mark(at = which(rownames(heatmat)%in%selGO.GN.dw==T), 
                                   labels = rownames(heatmat)[which(rownames(heatmat)%in%selGO.GN.dw==T)], 
                                   labels_gp = gpar(fontsize = 8, col="blue"), padding = unit(1, "mm")))
  
  
  pht_heat
}
heatmap.allgly.l2         <- Heatmap(as.matrix(gly.heat.plt), col = ScaleCols.l2,  
                                     name = "Heatmap", show_row_names=T,
                                     show_column_names = F, width = unit(2, "cm"),
                                     heatmap_legend_param =list(title = "log2FoldChange",
                                                                title_position = "lefttop-rot",
                                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                                labels_gp = gpar(fontsize = 10),
                                                                legend_width = unit(0.8,"cm"),
                                                                legend_height= unit(4, "cm")),
                                     cluster_rows = T, 
                                     row_names_gp = gpar(fontsize = 6),
                                     row_names_side="left",
                                     show_row_dend=T)
heatmap.fc1.3.l2          <- Heatmap_Top_fn2(gly.heat.plt, ScaleCols.l2 , gly.1.3) #.up, gly.1.3.dw)
heatmap.fc1.4.l2          <- Heatmap_Top_fn2(gly.heat.plt, ScaleCols.l2 , gly.1.4)#.up, gly.1.4.dw)
heatmap.fc1.5.l2          <- Heatmap_Top_fn2(gly.heat.plt, ScaleCols.l2 , gly.1.5)#.up, gly.1.5.dw)
heatmap.fc2.0.l2          <- Heatmap_Top_fn2(gly.heat.plt, ScaleCols.l2 , gly.2.0)#.up, gly.2.0.dw)
model                     <- "pair"



pdf(paste0(Project, "-", model, "_GO0006486_allG_Fig.l2fc.newHeatmap.pdf"), width=4, height= 20)
draw(heatmap.allgly.l2, heatmap_legend_side = "right", merge_legend=T)
dev.off()


pdf(paste0(Project, "-", model, "_GO0006486_allG_fc1.3_Fig.l2fc.newHeatmap.pdf"), width=4, height= 20)
draw(heatmap.fc1.3.l2, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_GO0006486_allG_fc1.4_Fig.l2fc.newHeatmap.pdf"), width=4, height= 20)
draw(heatmap.fc1.4.l2, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_GO0006486_allG_fc1.5_Fig.l2fc.newHeatmap.pdf"), width=4, height= 20)
draw(heatmap.fc1.5.l2, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-", model, "_GO0006486_allG_fc2.0_Fig.l2fc.newHeatmap.pdf"), width=4, height= 20)
draw(heatmap.fc2.0.l2, heatmap_legend_side = "right", merge_legend=T)
dev.off()



message("+           summary table and test for lose the fold change  --------------+")

fc1.3  <- table(res.selGOGly.dat$padj<0.05, abs(res.selGOGly.dat$log2FoldChange) >= log2(1.3))
fc1.4  <- table(res.selGOGly.dat$padj<0.05, abs(res.selGOGly.dat$log2FoldChange) >= log2(1.4))
fc1.5  <- table(res.selGOGly.dat$padj<0.05, abs(res.selGOGly.dat$log2FoldChange) >= log2(1.5))
fc2    <- table(res.selGOGly.dat$padj<0.05, abs(res.selGOGly.dat$log2FoldChange) >= log2(2))

## fisher-exact test with two-tailed
p1.3 <- fisher.test(fc1.3, alternative = "two.sided")$p.value
p1.4 <- fisher.test(fc1.4, alternative = "two.sided")$p.value
p1.5 <- fisher.test(fc1.5, alternative = "two.sided")$p.value
p2.0 <- fisher.test(fc2, alternative = "two.sided")$p.value

table.data <- matrix(c(fc1.3[2,2], fc1.4[2,2],fc1.5[2,2],fc2[2,2],p1.3,p1.4,p1.5,p2.0), ncol=2, byrow=F)
colnames(table.data) <- c("No.Genes", "Fisher_Exact_test_2tail")
rownames(table.data) <- paste0("FC", c(1.3,1.4,1.5,2.0))

write.table(table.data, file = "CTR_gjb2_0002-GO0006486_sigP_l2fc_var_summary_Table.txt")

write.table(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(1.3),],
            file = "CTR_gjb2_0002-GO0006486_fc1.3_Table.txt")
write.table(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(1.4),],
            file = "CTR_gjb2_0002-GO0006486_fc1.4_Table.txt")
write.table(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(1.5),],
            file = "CTR_gjb2_0002-GO0006486_fc1.5_Table.txt")
write.table(res.selGOGly.dat[res.selGOGly.dat$padj<0.05&abs(res.selGOGly.dat$log2FoldChange) >= log2(2.0),],
            file = "CTR_gjb2_0002-GO0006486_fc2.0_Table.txt")

message("+-------------------------------------------------------------------------------+")
message("+Generate Heatmap for specific pathways GO:0006486 with sigGenes ---------------+")
message("+-------------------------------------------------------------------------------+")

Heatmap.1.3 <- subset(res.selGOGly.dat, padj < 0.05 & abs(res.selGOGly.dat$log2FoldChange) >= log2(1.3))[,c(2,10:19)]
Heatmap.1.4 <- subset(res.selGOGly.dat, padj < 0.05 & abs(res.selGOGly.dat$log2FoldChange) >= log2(1.4))[,c(2,10:19)]
Heatmap.1.5 <- subset(res.selGOGly.dat, padj < 0.05 & abs(res.selGOGly.dat$log2FoldChange) >= log2(1.5))[,c(2,10:19)]
Heatmap.2.0 <- subset(res.selGOGly.dat, padj < 0.05 & abs(res.selGOGly.dat$log2FoldChange) >= log2(2.0))[,c(2,10:19)]

Heatmap.1.3.plt           <- Heatmap.1.3[,-1]
rownames(Heatmap.1.3.plt) <- Heatmap.1.3[,1]
colnames(Heatmap.1.3.plt) <- gsub("[.]", "-", colnames(Heatmap.1.3.plt))

Heatmap.1.4.plt           <- Heatmap.1.4[,-1]
rownames(Heatmap.1.4.plt) <- Heatmap.1.4[,1]
colnames(Heatmap.1.4.plt) <- gsub("[.]", "-", colnames(Heatmap.1.4.plt))

Heatmap.1.5.plt           <- Heatmap.1.5[,-1]
rownames(Heatmap.1.5.plt) <- Heatmap.1.5[,1]
colnames(Heatmap.1.5.plt) <- gsub("[.]", "-", colnames(Heatmap.1.5.plt))

Heatmap.2.0.plt           <- Heatmap.2.0[,-1]
rownames(Heatmap.2.0.plt) <- Heatmap.2.0[,1]
colnames(Heatmap.2.0.plt) <- gsub("[.]", "-", colnames(Heatmap.2.0.plt))

o2.1.3     <- seriate(dist(t(Heatmap.1.3.plt)), method = "GW")
o2.1.4     <- seriate(dist(t(Heatmap.1.4.plt)), method = "GW")
o2.1.5     <- seriate(dist(t(Heatmap.1.5.plt)), method = "GW")
o2.2.0     <- seriate(dist(t(Heatmap.2.0.plt)), method = "GW")

pdf(paste0(Project, "-GO0006486_FC1.3_Fig.Heatmap.pdf"), width=6, height= 8)
hplt.1.3   <- Heatmap_Top_fn(Heatmap.1.3.plt, ScaleCols.deg, o2.1.3, 1:10, ha)
draw(hplt.1.3, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-GO0006486_FC1.4_Fig.Heatmap.pdf"), width=6, height= 8)
hplt.1.4   <- Heatmap_Top_fn(Heatmap.1.4.plt, ScaleCols.deg, o2.1.4, 1:10, ha)
draw(hplt.1.4, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-GO0006486_FC1.5_Fig.Heatmap.pdf"), width=6, height= 8)
hplt.1.5   <- Heatmap_Top_fn(Heatmap.1.5.plt, ScaleCols.deg, o2.1.5, 1:10, ha)
draw(hplt.1.5, heatmap_legend_side = "right", merge_legend=T)
dev.off()

pdf(paste0(Project, "-GO0006486_FC2.0_Fig.Heatmap.pdf"), width=6, height= 8)
hplt.2.0   <- Heatmap_Top_fn(Heatmap.2.0.plt, ScaleCols.deg, o2.2.0, 1:10, ha)
draw(hplt.2.0, heatmap_legend_side = "right", merge_legend=T)
dev.off()

message("+-------------------------------------------------------------------------------+")
message("+Check the number of genes relating to GO:0006487 pathways. ------------------- +")
message("+-------------------------------------------------------------------------------+")

res.selGONGly.dat  <- res.merENID[res.merENID$go_id=="GO:0006487",]
res.selGONGly.dat  <- res.selGONGly.dat[-which(is.na(res.selGONGly.dat$go_id)),]
allNGly.Glength    <- unique(res.selGONGly.dat$ensembl_gene_id) ## 45

write.table(res.selGONGly.dat, file = "CTR_gjb2_0002-GO0006487_summary_Table.txt")

nfc1.3  <- table(res.selGONGly.dat$padj<0.05, abs(res.selGONGly.dat$log2FoldChange) >= log2(1.3))
nfc1.4  <- table(res.selGONGly.dat$padj<0.05, abs(res.selGONGly.dat$log2FoldChange) >= log2(1.4))
nfc1.5  <- table(res.selGONGly.dat$padj<0.05, abs(res.selGONGly.dat$log2FoldChange) >= log2(1.5))
nfc2    <- table(res.selGONGly.dat$padj<0.05, abs(res.selGONGly.dat$log2FoldChange) >= log2(2))

## fisher-exact test with two-tailed
np1.3 <- fisher.test(nfc1.3, alternative = "two.sided")$p.value
np1.4 <- fisher.test(nfc1.4, alternative = "two.sided")$p.value
np1.5 <- fisher.test(nfc1.5, alternative = "two.sided")$p.value
np2.0 <- fisher.test(nfc2, alternative = "two.sided")$p.value

ntable.data <- matrix(c(nfc1.3[2,2], nfc1.4[2,2],nfc1.5[2,2],nfc2[2,2],np1.3,np1.4,np1.5,np2.0), ncol=2, byrow=F)
colnames(ntable.data) <- c("No.Genes", "Fisher_Exact_test_2tail")
rownames(ntable.data) <- paste0("FC", c(1.3,1.4,1.5,2.0))

write.table(ntable.data, file = "CTR_gjb2_0002-GO0006487_sigP_l2fc_var_summary_Table.txt")

save(fc1.3,fc1.4,fc1.5,fc2,nfc1.3,nfc1.4,nfc1.5,nfc2,file="CTR_gjb2_0002-GO486_487_FisherExacTest_table.RData")

write.table(res.selGONGly.dat[res.selGONGly.dat$padj<0.05&abs(res.selGONGly.dat$log2FoldChange) >= log2(1.3),],
            file = "CTR_gjb2_0002-GO0006487_fc1.3_Table.txt")
write.table(res.selGONGly.dat[res.selGONGly.dat$padj<0.05&abs(res.selGONGly.dat$log2FoldChange) >= log2(1.4),],
            file = "CTR_gjb2_0002-GO0006487_fc1.4_Table.txt")
write.table(res.selGONGly.dat[res.selGONGly.dat$padj<0.05&abs(res.selGONGly.dat$log2FoldChange) >= log2(1.5),],
            file = "CTR_gjb2_0002-GO0006487_fc1.5_Table.txt")
write.table(res.selGONGly.dat[res.selGONGly.dat$padj<0.05&abs(res.selGONGly.dat$log2FoldChange) >= log2(2.0),],
            file = "CTR_gjb2_0002-GO0006487_fc2.0_Table.txt")

message("+ ---------------------Finished----------------------------------------------------------------------+")


























##------------FINISH---------------------------------------------------------##

library(Seurat)
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)
# H9 and Lonza separately graphed----

h9.ecto <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/H9_Ecto_dense_expression_matrix.csv")
h9.meso <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/H9_Meso_dense_expression_matrix.csv")
h9.endo <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/H9_Endo_dense_expression_matrix.csv")
h9.es <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/H9_ES_dense_expression_matrix.csv")

ipsc.ecto <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/iPSC_Lonza_Ecto_dense_expression_matrix.csv")
ipsc.meso <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/iPSC_Lonza_Meso_dense_expression_matrix.csv")
ipsc.endo <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/iPSC_Lonza_Endo_dense_expression_matrix.csv")
ipsc.es <- fread("/data/NCATS_ifx/iPSC/IS006/gene_symbol_matrices/iPSC_Lonza_ES_dense_expression_matrix.csv")


#
reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  x <- CreateSeuratObject(raw.data = x, project = "NPC")
  x@meta.data$sample <- samplename
  x <- NormalizeData(x)
  x <- ScaleData(x, display.progress = F)
  return(x)
}

h9.ecto <- reformat_for_seurat(h9.ecto, "h9.ecto")
h9.endo <- reformat_for_seurat(h9.endo, "h9.endo")
h9.meso <- reformat_for_seurat(h9.meso, "h9.meso")
h9.es <- reformat_for_seurat(h9.es, "h9.es")

is006.h9 <- MergeSeurat(h9.ecto, h9.endo, project="NPC")
is006.h9 <- MergeSeurat(is006.h9, h9.meso, project="NPC")
is006.h9 <- MergeSeurat(is006.h9, h9.es, project="NPC")

is006.h9 <- SetAllIdent(object = is006.h9, id = "sample")
is006.h9 <- NormalizeData(object = is006.h9, normalization.method = "LogNormalize", scale.factor = 10000)
is006.h9 <- ScaleData(object = is006.h9, vars.to.regress = c("nUMI"))
is006.h9 <- FindVariableGenes(object = is006.h9, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

is006.h9 <- RunPCA(object = is006.h9)
PCAPlot(object=is006.h9)
is006.h9 <- RunTSNE(is006.h9)
TSNEPlot(object = is006.h9)


pdf("IS006.H9.TSNE.pdf",width=10,height=10) 
TSNEPlot(object = is006.h9)
dev.off()

## iPSC----
ipsc.ecto <- reformat_for_seurat(ipsc.ecto, "ipsc.ecto")
ipsc.endo <- reformat_for_seurat(ipsc.endo, "ipsc.endo")
ipsc.meso <- reformat_for_seurat(ipsc.meso, "ipsc.meso")
ipsc.es <- reformat_for_seurat(ipsc.es, "ipsc.es")

is006.ipsc <- MergeSeurat(ipsc.ecto, ipsc.endo, project="NPC")
is006.ipsc <- MergeSeurat(is006.ipsc, ipsc.meso, project="NPC")
is006.ipsc <- MergeSeurat(is006.ipsc, ipsc.es, project="NPC")

is006.ipsc <- SetAllIdent(object = is006.ipsc, id = "sample")
is006.ipsc <- NormalizeData(object = is006.ipsc, normalization.method = "LogNormalize", scale.factor = 10000)
is006.ipsc <- ScaleData(object = is006.ipsc, vars.to.regress = c("nUMI"))
is006.ipsc <- FindVariableGenes(object = is006.ipsc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

is006.ipsc <- RunPCA(object = is006.ipsc)
PCAPlot(object=is006.ipsc)
is006.ipsc <- RunTSNE(is006.ipsc)
TSNEPlot(object = is006.ipsc)



pdf("IS006.iPSC.TSNE.pdf",width=10,height=10) 
TSNEPlot(object = is006.ipsc)
dev.off()

# loading in from cluster data work -----
rm(list=ls())
setwd("/Volumes/ncatssctl/NGS_related/Chromium/IS006/expression_matrices/gene_symbol_matrices")
load("IS006.H9.Seurat.Rdata")

# DE analysis, H9.-----

es.markers <- FindMarkers(object = is006.h9, ident.1 = "h9.es", min.pct = 0.25)
es.markers <- row.names(es.markers[1:20,])
meso.markers <- FindMarkers(object = is006.h9, ident.1 = "h9.meso", min.pct = 0.25)
meso.markers <- row.names(meso.markers[1:20,])
endo.markers <- FindMarkers(object = is006.h9, ident.1 = "h9.endo", min.pct = 0.25)
endo.markers <- row.names(endo.markers[1:20,])
ecto.markers <- FindMarkers(object = is006.h9, ident.1 = "h9.ecto", min.pct = 0.25)
ecto.markers <- row.names(ecto.markers[1:20,])

# h9 violin plots----
data.h9 <- as.data.frame(as.matrix(is006.h9@data))

scale.data.h9 <- as.data.frame(as.matrix(is006.h9@scale.data))
scale.data.h9$GeneId <- row.names(scale.data.h9)
scale.data.h9 <- as.data.table(scale.data.h9)
scale.data.h9[GeneId =="T",]
scale.data.h9[GeneId =="TBXT",]
scale.data.h9[GeneId =="ESRG"]
scale.data.h9[GeneId =="HESRG"]

raw.data.h9 <- as.data.frame(as.matrix(is006.h9@raw.data))
raw.data.h9$GeneId <- row.names(raw.data.h9)
raw.data.h9 <- as.data.table(raw.data.h9)
raw.data.h9[GeneId =="ESRG"]
raw.data.h9[GeneId =="HESRG"]

data.h9 <- as.data.frame(as.matrix(is006.h9@data))
data.h9$GeneId <- row.names(data.h9)
data.h9 <- as.data.table(data.h9)
h9.subset <- data.h9[GeneId %in% c("TBXT", "SOX17", "POU5F1", "PAX6"),]
h9.subset <- as.data.table(t(h9.subset))
names(h9.subset)<- unlist(h9.subset[16583,])
h9.subset <- h9.subset[1:16582,]
h9.subset$barcode <- names(as.data.frame(as.matrix(is006.h9@data)))
h9.subset[,tissue:= tstrsplit(barcode, "_")[2]]
h9.subset[,PAX6:=as.numeric(PAX6)]
h9.subset[,TBXT:=as.numeric(TBXT)]
h9.subset[,SOX17:=as.numeric(SOX17)]
h9.subset[,POU5F1:=as.numeric(POU5F1)]

h9.subset[,tissue.reorder := factor(tissue, levels = c("H9.ES","H9.meso","H9.endo", "H9.ecto"))]
h9.subset[tissue.reorder=="H9.ES",tissue.reorder:="ES"]
h9.subset[tissue.reorder=="H9.endo",tissue.reorder:="Endo"]
h9.subset[tissue.reorder=="H9.ecto",tissue.reorder:="Ecto"]
h9.subset[tissue.reorder=="H9.meso",tissue.reorder:="Meso"]

pax6.plot <- ggplot(data=h9.subset, aes(x=tissue.reorder, y=PAX6)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="PAX6, normalized expression", title="PAX6")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto")) + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))

tbxt.plot <- ggplot(data=h9.subset, aes(x=tissue.reorder, y=TBXT)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="TBXT, normalized expression", title="TBXT")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto")) + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))

sox17.plot <- ggplot(data=h9.subset, aes(x=tissue.reorder, y=SOX17)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="SOX17, normalized expression", title="SOX17")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto")) + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))

POU5F1.plot <- ggplot(data=h9.subset, aes(x=tissue.reorder, y=POU5F1)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="POU5F1, normalized expression", title="POU5F1")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto")) + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))


plot.page <- cowplot::plot_grid(tbxt.plot, sox17.plot, POU5F1.plot, pax6.plot, ncol=2, nrow=2)
plot.page

#ipsc violin plots-----

data.ipsc <- as.data.frame(as.matrix(is006.ipsc@data))
data.ipsc$GeneId <- row.names(data.ipsc)
data.ipsc <- as.data.table(data.ipsc)
ipsc.subset <- data.ipsc[GeneId %in% c("TBXT", "SOX17", "POU5F1", "PAX6"),]
ipsc.subset <- as.data.table(t(ipsc.subset))
names(ipsc.subset)<- unlist(ipsc.subset[17184,])
ipsc.subset <- ipsc.subset[1:17183,]
ipsc.subset$barcode <- names(as.data.frame(as.matrix(is006.ipsc@data)))
ipsc.subset[,tissue:= tstrsplit(barcode, "_")[2]]
ipsc.subset[,PAX6:=as.numeric(PAX6)]
ipsc.subset[,TBXT:=as.numeric(TBXT)]
ipsc.subset[,SOX17:=as.numeric(SOX17)]
ipsc.subset[,POU5F1:=as.numeric(POU5F1)]

ipsc.subset[,tissue.reorder := factor(tissue, levels = c("iPSC.Lonza.ES","iPSC.Lonza.Meso","iPSC.Lonza.Endo", "iPSC.Lonza.Ecto"))]
ipsc.subset[tissue.reorder=="iPSC.Lonza.ES",tissue.reorder:="ES"]
ipsc.subset[tissue.reorder=="iPSC.Lonza.Endo",tissue.reorder:="Endo"]
ipsc.subset[tissue.reorder=="iPSC.Lonza.Ecto",tissue.reorder:="Ecto"]
ipsc.subset[tissue.reorder=="iPSC.Lonza.Meso",tissue.reorder:="Meso"]


pax6.plot <- ggplot(data=ipsc.subset, aes(x=tissue.reorder, y=PAX6)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="PAX6, normalized expression", title="PAX6")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto")) + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))

tbxt.plot <- ggplot(data=ipsc.subset, aes(x=tissue.reorder, y=TBXT)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="TBXT, normalized expression", title="TBXT")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto"))  + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))

sox17.plot <- ggplot(data=ipsc.subset, aes(x=tissue.reorder, y=SOX17)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="SOX17, normalized expression", title="SOX17")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto"))  + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))

POU5F1.plot <- ggplot(data=ipsc.subset, aes(x=tissue.reorder, y=POU5F1)) + theme_bw() + geom_jitter(aes(color=tissue.reorder)) + geom_boxplot(width=0.1) +
  labs(x="Tissue", y="POU5F1, normalized expression", title="POU5F1")+
  guides(color=FALSE)+
  scale_x_discrete(limits=c("ES","Meso","Endo", "Ecto"))  + 
  scale_colour_manual(values = c('#F8766D','#C77CFF','#00BFC4', '#7CAE00'))


plot.page <- cowplot::plot_grid(tbxt.plot, sox17.plot, POU5F1.plot, pax6.plot, ncol=2, nrow=2)


# heatmap-----
DoHeatmap(is006.h9, genes.use = c(es.markers, ecto.markers, endo.markers, meso.markers), col.low="green", col.high="red", slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))

DoHeatmap(is006.h9, genes.use = c("TERF1", "FOXD3-AS1", "POU5F1", "NANOG", "LCK",
                                "PAX6", "SOX11", "MAP2", "SOX21", "TPBG",
                                "LEFTY1", "SOX17", "FOXA2", "NODAL", "EOMES",
                                "TBXT", "DKK1", "TNFRSF11B", "MSX1", "VIM"),
          slim.col.label=TRUE, remove.key=TRUE, col.low="green", col.high="red",
          group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso")) ## note T = TBXT here

DoHeatmap(is006.h9, genes.use = c("TERF1", "FOXD3-AS1", "POU5F1", "NANOG", "RARRES2",
                                  "NNAT", "DLK1", "SOX11", "TFPI", "PTN",
                                  "LEFTY1", "SOX17", "FOXA2", "NODAL", "EOMES",
                                  "TBXT", "DKK1", "TNFRSF11B", "MSX1", "VIM"),
          slim.col.label=TRUE, remove.key=TRUE, col.low="green", col.high="red",
          group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso")) ## note T = TBXT here


# h9.subset <- data.h9[GeneId %in% c("TBXT", "SOX17", "POU5F1", "PAX6"),]
# h9.subset <- as.data.table(t(h9.subset))
# names(h9.subset)<- unlist(h9.subset[16583,])
# h9.subset <- h9.subset[1:16582,]
# h9.subset$barcode <- names(as.data.frame(as.matrix(is006.h9@data)))
# h9.subset[,tissue:= tstrsplit(barcode, "_")[2]]
# h9.subset[,PAX6:=as.numeric(PAX6)]
# h9.subset[,TBXT:=as.numeric(TBXT)]
# h9.subset[,SOX17:=as.numeric(SOX17)]
# h9.subset[,POU5F1:=as.numeric(POU5F1)]

# h9.subset[,tissue.reorder := factor(tissue, levels = c("H9.ES","H9.meso","H9.endo", "H9.ecto"))]



# DE analysis, iPSC-----
es.markers <- FindMarkers(object = is006.ipsc, ident.1 = "ipsc.es", min.pct = 0.25)
es.markers <- row.names(es.markers[1:20,])
meso.markers <- FindMarkers(object = is006.ipsc, ident.1 = "ipsc.meso", min.pct = 0.25)
meso.markers <- row.names(meso.markers[1:20,])
endo.markers <- FindMarkers(object = is006.ipsc, ident.1 = "ipsc.endo", min.pct = 0.25)
endo.markers <- row.names(endo.markers[1:20,])
ecto.markers <- FindMarkers(object = is006.ipsc, ident.1 = "ipsc.ecto", min.pct = 0.25)
ecto.markers <- row.names(ecto.markers[1:20,])

# data, ipsc---
data.ipsc <- as.data.frame(as.matrix(is006.ipsc@data))
data.ipsc$GeneId <- row.names(data.ipsc)
data.ipsc <- as.data.table(data.ipsc)

DoHeatmap(is006.ipsc, genes.use = c(es.markers, ecto.markers, endo.markers, meso.markers), col.low="green", col.high="red", slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))

DoHeatmap(is006.ipsc, genes.use = c("TERF1", "FOXD3-AS1", "POU5F1", "NANOG", "LCK",
                                  "PAX6", "SOX11", "MAP2", "SOX21", "TPBG",
                                  "LEFTY1", "SOX17", "FOXA2", "NODAL", "EOMES",
                                  "TBXT", "DKK1", "TNFRSF11B", "MSX1", "VIM"),
          slim.col.label=TRUE, remove.key=TRUE, col.low="green", col.high="red",
          group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso")) ## note T = TBXT here

DoHeatmap(is006.ipsc, genes.use = c("TERF1", "FOXD3-AS1", "POU5F1", "NANOG", "RARRES2",
                                  "NNAT", "DLK1", "SOX11", "TFPI", "PTN",
                                  "LEFTY1", "SOX17", "FOXA2", "NODAL", "EOMES",
                                  "TBXT", "DKK1", "TNFRSF11B", "MSX1", "VIM"),
          slim.col.label=TRUE, remove.key=TRUE, col.low="green", col.high="red",
          group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))


# convert ensg to gene_symbol-----

library(biomaRt)

ensg_to_symbol <- function(x){
  names(x)[1] <- "ENSG"
  x[,ENSG:=tstrsplit(ENSG, "\\.")[1]]
  x <- x[!duplicated(x$ENSG),]
  genes <- x$ENSG
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  out <- merge(x,G_list,by.x="ENSG",by.y="ensembl_gene_id", all=T)
  out <- out[!duplicated(out$external_gene_name),]
  out <- na.omit(out)
  return(out)
}

setwd("/Volumes/ncatssctl/NGS_related/Chromium/IS009/expression_matrices/")
files <- Sys.glob('*.csv')

for (i in 1:length(files)){
  sample <- substr(files[i], 7,(nchar(files[i])-28))
  file <- fread(files[i])
  file <- ensg_to_symbol(file)
  fwrite(file, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS009/expression_matrices/gene_symbol/", "IS009_", sample, "_dense_expression_matrix_genesymbol.csv"), sep="\t", col.names=T, row.names=F, quote=F)
}

# new interesting figure ----
load("/Volumes/ncatssctl/NGS_related/Chromium/IS006/Analysis/IS006.H9.Seurat.Rdata")
load("/Volumes/ncatssctl/NGS_related/Chromium/IS006/Analysis/IS006.iPSC.Seurat.Rdata")

is006 <- MergeSeurat(is006.h9, is006.ipsc, project="IS006")
rm(is006.h9, is006.ipsc)
is006 <- SetAllIdent(object = is006, id = "sample")
is006 <- ScaleData(is006)
is006 <- FindVariableGenes(is006, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is006 <- RunPCA(is006)
is006 <- RunTSNE(is006)

pdfPath <- getwd()
pdf(file=paste0(pdfPath, "/IS006.tSNE.pdf"), width=10, height=10)
TSNEPlot(is006)
dev.off()

# SC3 plots (run on cluster with 100g memory)----
source("https://bioconductor.org/biocLite.R")
biocLite("SC3")
library("SC3")
library("SingleCellExperiment")

load("IS006.H9.Seurat.Rdata")
load("IS006.iPSC.Seurat.Rdata")

is006.h9.data <- as.data.frame(as.matrix(is006.h9@raw.data))
is006.ipsc.data <- as.data.frame(as.matrix(is006.ipsc@raw.data))
rm(is006.h9, is006.ipsc)
gc()

ann.h9 <- data.table(cell_type1 = names(is006.h9.data))
ann.h9[, day := tstrsplit(cell_type1, "_")[2]]
ann.h9 <- as.data.frame(ann.h9)
row.names(ann.h9) <- ann.h9$cell_type1
ann.h9 <- ann.h9[-1]
names(ann.h9)[1]<-"cell_type1"

ann.ipsc <- data.table(cell_type1 = names(is006.ipsc.data))
ann.ipsc[, day := tstrsplit(cell_type1, "_")[2]]
ann.ipsc <- as.data.frame(ann.ipsc)
row.names(ann.ipsc) <- ann.ipsc$cell_type1
ann.ipsc <- ann.ipsc[-1]
names(ann.ipsc)[1] <- "cell_type1"

ann <- rbind(ann.h9, ann.ipsc)
#ann

yan <- cbind(is006.h9.data, is006.ipsc.data)# keeping the same object names as in http://127.0.0.1:30357/library/SC3/doc/SC3.html

# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(yan),
    logcounts = log2(as.matrix(yan) + 1)
  ), 
  colData = ann
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)

# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# pca plot
pdfPath<-getwd()
pdf(file=paste0(pdfPath, "/IS006.Scater.PCA.pdf"), width=10, height=10)
scater::plotPCA(sce, colour_by = "cell_type1")
dev.off()

#sce <- sc3(sce, ks=1:8, biology=TRUE) # takes a long while - 5 hours - doesn't finish
sce <- sc3(sce, ks=1:8, biology=FALSE) # does finish.


col_data <- colData(sce) # many NAs
row_data <- rowData(sce)

sc3_plot_consensus(
  sce, k = 8, 
  show_pdata = c(
    "cell_type1", 
    "log10_total_features",
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)

pdf(file=paste0(pdfPath, "/IS006.Scater-SCE.heatmap.4clusters.pdf"), width=20, height=20)
sc3_plot_consensus(
  sce, k = 8, 
  show_pdata = c(
    "cell_type1",
    "sc3_4_clusters"
  )
)
dev.off()
pdf(file=paste0(pdfPath, "/IS006.Scater-SCE.markers.pdf"), width=20, height=20)
sc3_plot_markers(
  sce, k = 3, 
  show_pdata = c(
    "cell_type1",
    "sc3_3_clusters"
  )
)
dev.off()

# IS006 heatmaps august 20----
glycolysis <- c("AKR1A1", "ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "GALM", "ADH7", "LDHAL6A", "DLAT", "DLD", "ENO1", "ENO2", "ENO3", "ALDH2", "ALDH3A1", "ALDH1B1", "ALDH1A3", "FBP1", "ALDH3B1", "ALDH3B2", "ALDH9A1", "ALDH3A2", "ALDOA", "ALDOB", "ALDOC", "G6PC", "GAPDH", "GCK", "GPI", "HK1", "HK2", "HK3", "LDHA", "LDHB", "LDHC", "PGAM4", "ALDH7A1", "PCK1", "PCK2", "PDHA1", "PDHA2", "PDHB", "PFKL", "PFKM", "PFKP", "PGAM1", "PGAM2", "PGK1", "PGK2", "PGM1", "PKLR", "PKM2", "PGM2", "ACSS2", "G6PC2", "BPGM", "TPI1", "ACSS1", "FBP2", "LDHAL6B")
glycolysis.keep <- c("AKR1A1", "ADH5", "GALM", "DLAT", "DLD", "ENO1", "ENO2", "ENO3", "ALDH2", "ALDH1B1", "ALDH1A3", "FBP1", "ALDH9A1", "ALDH3A2", "ALDOA", "ALDOC", "GAPDH", "GCK", "GPI", "HK1", "HK2", "LDHA", "LDHB", "ALDH7A1", "PCK2", "PDHA1", "PDHB", "PFKL", "PFKM", "PFKP", "PGAM1", "PGAM2", "PGK1", "PGM1", "PKM2", "PGM2", "ACSS2", "BPGM", "TPI1")
tca <- c("CS", "DLD", "DLST", "FH", "NNT", "LOC283398", "IDH2", "IDH3A", "IDH3B", "IDH3G", "MDH2", "OGDH", "ACO2", "SDHA", "SDHB", "SDHC", "SDHD", "SUCLA2P1", "LOC646675", "LOC646677", "LOC650667", "LOC650674", "LOC651820", "SUCLG2", "SUCLG1", "SUCLA2")
aerores <- c("ME3", "FAM54A", "FAM36A", "MDH1B", "COX10", "CS", "DLAT", "DLD", "DLST", "FH", "SIRT3", "NNT", "FXN", "BLOC1S1", "UQCR10", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "MDH1", "MDH2", "NDUFV1", "ACO1", "OGDH", "ACO2", "OXA1L", "PDHA1", "PDHA2", "PDHB", "DHTKD1", "OGDHL", "FAM54B", "SDHA", "SDHB", "SDHC", "SDHD", "SURF1", "UCN", "UQCRB", "UQCRC1", "UQCRC2", "UQCRH", "PANK2", "CHCHD5", "CAT", "SUCLG2", "SUCLG1", "SUCLA2", "SLC25A14", "COX19", "PMPCB", "MTFR1")
aerores.keep <- c("ME3", "FAM54A", "FAM36A", "MDH1B", "COX10", "CS", "DLAT", "DLD", "DLST", "FH", "SIRT3", "NNT", "FXN", "BLOC1S1", "UQCR10", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "MDH1", "MDH2", "NDUFV1", "ACO1", "OGDH", "ACO2", "OXA1L", "PDHA1", "PDHB", "DHTKD1", "OGDHL", "FAM54B", "SDHA", "SDHB", "SDHC", "SDHD", "SURF1", "UQCRB", "UQCRC1", "UQCRC2", "UQCRH", "PANK2", "CHCHD5", "CAT", "SUCLG2", "SUCLG1", "SUCLA2", "SLC25A14", "COX19", "PMPCB", "MTFR1")
ppp <- c("RPIA", "SHPK", "G6PD", "PGLS", "DERA", "PGD", "PGM2", "RPE", "TALDO1", "TKT", "TPI1", "LOC729020", "H6PD")
pluri <- c("DNMT3B", "HESX1", "LCK", "NANOG", "POU5F1", "SOX2", "CDH1", "CD30 ", "SSEA4", "SSEA3", "GDF3", "EPCAM", "REX1", "LIN28", "PODXL", "SALL4", "EHMT2", "APOE", "CDH3", "ERBB3", "LEFTY1", "GCTM2", "CD24", "DUSP6", "ZFP42", "ESRG")

DoHeatmap(is006.ipsc, genes.use=c(glycolysis.keep), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))
DoHeatmap(is006.ipsc, genes.use=c(tca), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))
DoHeatmap(is006.ipsc, genes.use=c(aerores.keep), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))
DoHeatmap(is006.ipsc, genes.use=c(ppp), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))
DoHeatmap(is006.ipsc, genes.use=c(pluri), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))

DoHeatmap(is006.h9, genes.use=c(glycolysis.keep), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))
DoHeatmap(is006.h9, genes.use=c(tca), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))
DoHeatmap(is006.h9, genes.use=c(aerores.keep), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))
DoHeatmap(is006.h9, genes.use=c(ppp), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))
DoHeatmap(is006.h9, genes.use=c(pluri), col.low="blue", col.mid="white", col.high="red", slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))


# IS006 heatmaps for peroxisomal genes, august 24---
# on the cluster
load("IS006.H9.Seurat.Rdata")
peroxi <- c("ABCD1", "ABCD2", "ABCD3", "ABCD4", "ABP1", "ACAA1", "ACOT4", "ACOT8", "ACOX1", "ACOX2", "ACOX3", "AGPS", "AGXT", "ALDH3A2", "AMACR", "ATAD1", "BAAT", "CAT", "CNOT1", "CRAT", "CROT", "DAO", "DDO", "DHRS4", "DNM1L", "ECH1", "EHHADH", "EPHX2", "FAR1", "FAR2", "FIS1", "GNPAT", "GSTK1", "HACL1", "HAO1", "HAO2", "HMGCL", "HMGCR", "HSD17B4", "HSDL2", "IDE", "IDH1", "IDI1", "IDI2", "LKAP", "LONP2", "MAP2K2", "MAVS", "MFF", "MLYCD", "MPV17L", "MUL1", "NOS2", "NUDT12", "NUDT7", "PECI", "PECR", "PEX1", "PEX10", "PEX11A", "PEX11B", "PEX11G", "PEX12", "PEX13", "PEX14", "PEX16", "PEX19", "PEX26", "PEX3", "PEX5", "PEX6", "PEX7", "PHYH", "PIPOX", "PMVK", "PNPLA8", "PRDX5", "PXMP2", "PXMP3", "PXMP4", "SCP2", "SLC25A17", "SLC27A2", "SOD1", "TRIM37", "TTC1", "TYSND1", "VIM")

pdfPath <- getwd()
getwd()
pdf(file=paste0(pdfPath, "/IS006.H9.peroxisome.pdf"), width=10, height=20)
DoHeatmap(is006.h9, genes.use=c(peroxi), col.low="green", col.high="red",
          slim.col.label = T, remove.key = T, group.order = c("h9.es", "h9.ecto", "h9.endo", "h9.meso"))
dev.off()

pdf(file=paste0(pdfPath, "/IS006.iPSC.peroxisome.pdf"), width=10, height=20)
DoHeatmap(is006.ipsc, genes.use=c(peroxi), col.low="green", col.high="red",
          slim.col.label = T, remove.key = T, group.order = c("ipsc.es", "ipsc.ecto", "ipsc.endo", "ipsc.meso"))
dev.off()


#


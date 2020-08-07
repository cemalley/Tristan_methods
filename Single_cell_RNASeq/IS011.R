# Single cell RNASeq analysis pipeline with Seurat and gene set enrichment analysis
# Claire Malley 2019
# NCATS NIH

library(Seurat)
library(data.table)
library(biomaRt)

setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS011/Countfiles_ENSG/')
files <- Sys.glob('*dense*csv')

for (file in files){
  dt <- as.data.frame(fread(file))
  names(dt)[1] <- 'ENSG'
  dt <- as.data.table(dt)
  ensg.genes <- data.table("ENSG" = dt$ENSG)
  genes <- ensg.genes$ENSG
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG %in% ensg.genes$ENSG)
  dt <- merge(dt, ensg.genes, by="ENSG")
  dt <- dt[,ENSG:=NULL]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "_dense_expression_matrix.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS011/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
  
}

setwd('../Countfiles_gene_symbol/')
files <- Sys.glob('*gene_symbol.csv')

for (file in files){
  dt <- fread(file, header=T)
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  sample_id <- str_split_fixed(file, "_gene_symbol.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS011/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
}

#

ipsc.meso <- fread("/Volumes/ncatssctl/NGS_related/Chromium/IS006/expression_matrices/gene_symbol_matrices/iPSC_Lonza_Meso_dense_expression_matrix.csv")
ipsc.endo <- fread("/Volumes/ncatssctl/NGS_related/Chromium/IS006/expression_matrices/gene_symbol_matrices/iPSC_Lonza_Endo_dense_expression_matrix.csv")
ipsc.es <- fread("/Volumes/ncatssctl/NGS_related/Chromium/IS006/expression_matrices/gene_symbol_matrices/iPSC_Lonza_ES_dense_expression_matrix.csv")


#
reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  x <- CreateSeuratObject(raw.data = x, project = "NPC")
  x@meta.data$sample <- samplename
  #x <- NormalizeData(x)
 # x <- ScaleData(x, display.progress = F)
  return(x)
}

ipsc.ecto <- reformat_for_seurat(ipsc.ecto, "ipsc.ecto")
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

# /data/NCATS_ifx/iPSC/IS006/IS006.H9.Seurat.Rdata
# /data/NCATS_ifx/iPSC/IS011/IS006.11.Seurat.RData
# merge.

x <- is006.11.ipsc.h9

#x <- NormalizeData(x)
x <- FindVariableFeatures(x)
x <- ScaleData(x)
x <- RunPCA(x)
x <- FindNeighbors(x)
x <- FindClusters(x)
x <- RunTSNE(x)

x@meta.data$sample_rename <- x@meta.data$sample
x@meta.data$sample_rename <- gsub('ipsc.ES', 'iPSC Pluripotent', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('ipsc.endo', 'iPSC Endoderm', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('ipsc.meso', 'iPSC Mesoderm', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('ipsc.ecto', 'iPSC Ectoderm', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('h9.ecto', 'ESC Ectoderm', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('h9.endo', 'ESC Endoderm', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('h9.meso', 'ESC Mesoderm', x@meta.data$sample_rename)
x@meta.data$sample_rename <- gsub('h9.es', 'ESC Pluripotent', x@meta.data$sample_rename)
x@meta.data$sample_rename <- factor(x@meta.data$sample_rename, levels=
              c('ESC Pluripotent','iPSC Pluripotent',
                'ESC Ectoderm', 'iPSC Ectoderm',
                'ESC Endoderm','iPSC Endoderm',
                'ESC Mesoderm','iPSC Mesoderm'))
Idents(x) <- 'sample_rename'

tsne.plot <- TSNEPlot(x, pt.size=1)
pdf("IS006-11.iPSC.H9.tSNE.pdf", height = 4.8, width= 7.6)
plot(tsne.plot)
dev.off()

pca.plot <- PCAPlot(x, pt.size=1)
pdf("IS006-11.iPSC.H9.PCA.pdf", height = 4.8, width= 7.6)
plot(pca.plot)
dev.off()

pca.plot2 <- PCAPlot(x,dims=2:3, pt.size=1)
pdf("IS006-11.iPSC.H9.PCA_dims2-3.pdf", height = 4.8, width= 7.6)
plot(pca.plot2)
dev.off()


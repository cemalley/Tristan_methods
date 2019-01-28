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


reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- x[external_gene_name != "NA",]
  x <- subset(x, select=c("external_gene_name", unlist(names(x))[2:(length(x)-1)] ))
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  barcodes <- names(x)[2:length(x)]
  barcodes <- paste0(barcodes, paste0(".", samplename))
  names(x)[2:length(x)] <- barcodes
  x <- CreateSeuratObject(raw.data = x, project = "IS011")
  x@meta.data$sample <- samplename
  #x <- NormalizeData(x)
  #x <- ScaleData(x, display.progress = F)
  return(x)
}


IS011_2D_endo <- reformat_for_seurat(fread(files[1]), "2D_endo")
IS011_iRPE <- reformat_for_seurat(fread(files[2]), "iRPE")
IS011_3D_endo<- reformat_for_seurat(fread(files[3]), "3D_endo")
IS011_Lonza_NSC_auto_2D<- reformat_for_seurat(fread(files[5]), "Lonza_NSC_auto_2D")


# merging IS011 and IS013, only partial samples merge-----

is011.13 <- MergeSeurat(IS011_2D_endo, IS011_iRPE, project = "IS011.13", do.normalize = F)
is011.13 <- MergeSeurat(is011.13, IS011_3D_endo, project="IS011.13", do.normalize = F)
is011.13 <- MergeSeurat(is011.13, IS013_2D_RPE, project='IS011.13', do.normalize=F)
is011.13 <- MergeSeurat(is011.13, IS013_tissue_RPE, project="IS011.13", do.normalize = F)


is011.13 <- NormalizeData(is011.13, normalization.method = "LogNormalize", scale.factor = 10000)
is011.13 <- ScaleData(is011.13)
is011.13 <- FindVariableGenes(object = is011.13, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
is011.13 <- SetAllIdent(is011.13, id="sample")
is011.13 <- RunPCA(object = is011.13, pcs.print = 1:5, genes.print = 5)
#PCAPlot(is011.13)

is011.13 <- RunTSNE(object = is011.13, dims.use = 1:8, do.fast = TRUE)
#tsneplotdata <- TSNEPlot(object = is011.13)

pdfPath <- '/Volumes/ncatssctl/NGS_related/Chromium/IS011/'
pdf(file=paste0(pdfPath, "/IS011.13.PCA.pdf"), width=10, height=10)
print(PCAPlot(is011.13))
dev.off()
pdf(file=paste0(pdfPath, "/IS011.13.TSNE.pdf"), width=10, height=10)
print(TSNEPlot(is011.13))
dev.off()
setwd(pdfPath)
save(is011.13, file="IS011.13.partialmerge.Seurat.RData")

# top genes for each group ----

IS011_2D_endo_markers <- FindMarkers(object = is011.13, ident.1 = "2D_endo", min.pct = 0.25)
IS011_iRPE_markers <- FindMarkers(object = is011.13, ident.1 = "iRPE", min.pct = 0.25)
IS011_3D_endo_markers <- FindMarkers(object = is011.13, ident.1 = "3D_endo", min.pct = 0.25)
IS013_2D_RPE <- FindMarkers(object = is011.13, ident.1 = "2D_RPE", min.pct = 0.25)
IS013_tissue_RPE <- FindMarkers(object = is011.13, ident.1 = "Tissue_RPE", min.pct = 0.25)

markers <- c(row.names(IS011_2D_endo_markers)[1:10], row.names(IS011_iRPE_markers )[1:10], row.names(IS011_3D_endo_markers)[1:10], row.names(IS013_2D_RPE)[1:10], row.names(IS013_tissue_RPE )[1:10])
#markers <- unique(markers)


fwrite(IS011_2D_endo_markers, 'IS011_2D_endo_markers.csv', sep=',', col.names = T, row.names = T, quote=F)
fwrite(IS011_iRPE_markers, 'IS011_iRPE_markers.csv', sep=',', col.names = T, row.names = T, quote=F)
fwrite(IS011_3D_endo_markers, 'IS011_3D_endo_markers.csv', sep=',', col.names = T, row.names = T, quote=F)
fwrite(IS013_2D_RPE, 'IS013_2D_RPE.csv', sep=',', col.names = T, row.names = T, quote=F)
fwrite(IS013_tissue_RPE, 'IS013_tissue_RPE.csv', sep=',', col.names = T, row.names = T, quote=F)

# heatmap----

pdf(file=paste0(pdfPath, "/IS011.13.heatmap.pdf"), width=11, height=10)
print(DoHeatmap(is011.13, genes.use=c(markers), col.low = "green", col.high="red", remove.key = F, slim.col.label = TRUE, group.order = c('2D_endo', 'iRPE', '3D_endo', '2D_RPE', 'Tissue_RPE'), group.label.rot = T))
dev.off()

#


#
cluster.averages.IS011.ecto <- AverageExpression(object = IS011_Lonza_NSC_auto_2D)

# try to merge this sample with older trilineage samples


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







library(Seurat)
library(data.table)
library(biomaRt)
library(stringr)
library(ggplot2)
#is018.counts <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS018/IS018_DNASG-CT-4283_S10_dense_expression_matrix.csv')
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS018/Countfiles_ENSG')

files <- Sys.glob('*csv')
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
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS018/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
  
}

setwd('../Countfiles_gene_symbol/')
files <- Sys.glob('*gene_symbol.csv')

for (file in files){
  dt <- fread(file, header=T)
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  sample_id <- str_split_fixed(file, "_gene_symbol.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS018/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
}


reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- x[external_gene_name != "NA",]
  x <- subset(x, select=c("external_gene_name", unlist(names(x))[2:(length(x)-1)] ))
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  barcodes <- names(x)
  barcodes <- paste0(barcodes, paste0(".", samplename))
  names(x) <- barcodes
  x <- CreateSeuratObject(counts = x, project = "IS018")
  x@meta.data$sample <- samplename
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x, features.use=c(as.data.frame(as.matrix(x@assays$RNA@data))))
  x <- RunPCA(x)
  x <- FindNeighbors(x)
  x <- FindClusters(x)
  x <- RunTSNE(x)
  return(x)
}


is018 <- reformat_for_seurat(fread(files[1]), "neuron")

DimPlot(object = is018, reduction = "tsne", do.return=T) + labs(title='IS018 neuron cells tSNE')

VlnPlot(is018,c("nCount_RNA","nFeature_RNA"))
VariableFeaturePlot(is018)

PCAPlot(is018)

RidgePlot(object = is018, feature = c("ZIC2", "PAX6", "POU4F1"))

DimHeatmap(object = is018, dims = 1:5, balanced = TRUE)

is018 <- JackStraw(object = is018, num.replicate = 100)
is018 <- ScoreJackStraw(object = is018, dims = 1:20)
JackStrawPlot(object = is018, dims = 1:20)
ElbowPlot(object = is018) # up to the max 20


is018 <- FindNeighbors(object = is018, dims = 1:20)
is018 <- FindClusters(object = is018, resolution = 0.5)
is018 <- RunUMAP(object = is018, dims = 1:20)
DimPlot(object = is018, reduction = "umap")

save(is018, file='../IS018.neurons.RData')
load('/Volumes/ncatssctl/NGS_related/Chromium/IS018/IS018.neurons.RData')


# markers for each cluster----

TSNEPlot(is018, label=T)

cluster0.2.10.vs.3.4.7.markers <- FindMarkers(object = is018, ident.1 = c(0,2,10), ident.2 = c(3,4,7), min.pct = 0.25)
cluster0.2.10.vs.3.4.7.markers

cluster0.2.10.vs.1.6.8.9.markers <- FindMarkers(object = is018, ident.1 = c(0,2,10), ident.2 = c(1,6,8,9), min.pct = 0.25)
cluster0.2.10.vs.1.6.8.9.markers

# ok, regressing out cell cycle.----

g2.m.genes <- fread('/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/Seurat/cell.cycle/G2M.genes.txt', header= F)
g2.m.genes <- g2.m.genes$V1
s.genes <- fread('/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/Seurat/cell.cycle/S.genes.txt', header=F)
s.genes <- s.genes$V1

is018 <- CellCycleScoring(object = is018, s.features = s.genes, g2m.features = g2.m.genes, set.ident = TRUE)

is018.nocycle <- ScaleData(object = is018, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = TRUE)
Idents(is018.nocycle) <- "Phase"
is018.nocycle <- RunPCA(is018.nocycle)
is018.nocycle <- RunTSNE(is018.nocycle)
TSNEPlot(object = is018.nocycle, label=T)
TSNEPlot(is018, label=T)

Idents(is018.nocycle) <- 'sample'
is018.nocycle <- FindNeighbors(object = is018.nocycle, dims = 1:20)
is018.nocycle <- FindClusters(object = is018.nocycle, resolution = 0.5)
is018.nocycle <- FindClusters(object = is018.nocycle, resolution = 0.8)
Idents(is018.nocycle) <- 'RNA_snn_res.0.5'
TSNEPlot(is018.nocycle, label=T)

Idents(is018.nocycle) <- 'RNA_snn_res.0.8'
TSNEPlot(is018.nocycle, label=T)

Idents(is018) <- 'RNA_snn_res.0.8'
TSNEPlot(is018, label=T)


save(is018.nocycle, file='../IS018.nocyle.Seurat.RData')

cluster.1.2.4.9.vs.0 <- FindMarkers(is018.nocycle, ident.1 = c(1,2,4,9), ident.2 = c(0), min.pct = 0.25)
cluster.1.2.4.9.vs.3.5 <- FindMarkers(is018.nocycle, ident.1 = c(1,2,4,9), ident.2 = c(3,5), min.pct = 0.25)

cluster.3.5.vs.0 <- FindMarkers(is018.nocycle, ident.1 = c(3,5), ident.2 = c(0), min.pct = 0.25)
cluster.3.5.vs.7.8 <- FindMarkers(is018.nocycle, ident.1 = c(3,5), ident.2 = c(7,8), min.pct = 0.25)

cluster.10.vs.rest <- FindMarkers(is018.nocycle, ident.1 = c(10), min.pct = 0.25)
cluster.6.vs.rest <- FindMarkers(is018.nocycle, ident.1 = c(6), min.pct = 0.25)

cluster.genes <- rbind(cluster.1.2.4.9.vs.0, cluster.1.2.4.9.vs.3.5, cluster.3.5.vs.0, cluster.3.5.vs.7.8, cluster.10.vs.rest, cluster.6.vs.rest)
cluster.genes <- subset(cluster.genes, p_val_adj <= 0.001 & abs(avg_logFC) >= 1)
cluster.genes$GeneId <- row.names(cluster.genes)
cluster.genes <- as.data.table(cluster.genes)
cluster.genes <- cluster.genes[order(p_val_adj, -abs(avg_logFC))]

DoHeatmap(is018.nocycle, features = c(cluster.genes$GeneId[1:50]) )

library(dplyr)
#top10 <- is018.nocycle %>% dplyr::group_by(seurat_clusters) %>% top_n(n = 10, wt = avg_logFC) doesnt work.
#DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()


cluster0 <- FindMarkers(is018.nocycle, ident.1=0, min.pct=0.25)
cluster1 <- FindMarkers(is018.nocycle, ident.1=1, min.pct=0.25)
cluster2 <- FindMarkers(is018.nocycle, ident.1=2, min.pct=0.25)
cluster3 <- FindMarkers(is018.nocycle, ident.1=3, min.pct=0.25)
cluster4 <- FindMarkers(is018.nocycle, ident.1=4, min.pct=0.25)
cluster5 <- FindMarkers(is018.nocycle, ident.1=5, min.pct=0.25)
cluster6 <- FindMarkers(is018.nocycle, ident.1=6, min.pct=0.25)
cluster7 <- FindMarkers(is018.nocycle, ident.1=7, min.pct=0.25)
cluster8 <- FindMarkers(is018.nocycle, ident.1=8, min.pct=0.25)
cluster9 <- FindMarkers(is018.nocycle, ident.1=9, min.pct=0.25)
cluster10 <- FindMarkers(is018.nocycle, ident.1=10, min.pct=0.25)

cluster0$cluster <- 0
cluster1$cluster <- 1
cluster2$cluster <- 2
cluster3$cluster <- 3
cluster4$cluster <- 4
cluster5$cluster <- 5
cluster6$cluster <- 6
cluster7$cluster <- 7
cluster8$cluster <- 8
cluster9$cluster <- 9
cluster10$cluster <- 10

cluster.genes <- rbind(cluster0, cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,cluster8,cluster9,cluster10)
cluster.genes$GeneId <- row.names(cluster.genes)

library(dplyr)
cluster.genes.top10 <- cluster.genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cluster.genes.top10
DoHeatmap(is018.nocycle, features = c(cluster.genes.top10$GeneId) , raster=F) + labs(title='IS018 neurons: top 10 genes per cluster')

plot <- DimPlot(object=is018, reduction="tsne")

topcluster <- CellSelector(plot = plot)
leftcluster <- CellSelector(plot=plot)
rightcluster <- CellSelector(plot=plot)
bottom2 <- CellSelector(plot=plot)
bottom1 <- CellSelector(plot=plot)

new_grouping <- data.table(cluster=c(rep(1,length(topcluster)), rep(2, length(leftcluster)), rep(3, length(rightcluster)), rep(4, length(bottom1)), rep(5, length(bottom2))), c(topcluster, leftcluster, rightcluster, bottom1, bottom2))

new_grouping <- new_grouping[!duplicated(V2)]
names(new_grouping)[2] <- 'barcode'

new_grouping

is018.grouped <- SubsetData(is018.nocycle, cells=new_grouping$barcode)
is018.grouped@meta.data

data.order <- data.table('meta.cells'= row.names(is018.grouped@meta.data), 'order'=1:3875)
data.order <- merge(data.order, new_grouping, by.x='meta.cells', by.y='barcode')
data.order <- data.order[order(order)]

all(row.names(is018.grouped@meta.data) == data.order$meta.cells)

is018.grouped@meta.data$new_grouping <- data.order$cluster
Idents(is018.grouped) <- 'new_grouping'

TSNEPlot(is018.grouped)

cells.add <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.add

unique(data.order[meta.cells %in% cells.add,cluster])
data.order[meta.cells %in% cells.add, cluster := 1]

is018.grouped@meta.data$new_grouping <- data.order$cluster

Idents(is018.grouped) <- 'new_grouping'

DimPlot(object=is018.grouped, reduction='tsne', label=T)

cells.add.7.8 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))

data.order[meta.cells %in% cells.add.7.8, cluster := 6]
is018.grouped@meta.data$new_grouping <- data.order$cluster

Idents(is018.grouped) <- 'new_grouping'

DimPlot(object=is018.grouped, reduction='tsne', label=T)

cells.add.5 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.add.3 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.add.2 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))


data.order[meta.cells %in% cells.add.5, cluster := 5]
is018.grouped@meta.data$new_grouping <- data.order$cluster

data.order[meta.cells %in% cells.add.3, cluster := 3]
is018.grouped@meta.data$new_grouping <- data.order$cluster

data.order[meta.cells %in% cells.add.2, cluster := 2]
is018.grouped@meta.data$new_grouping <- data.order$cluster


Idents(is018.grouped) <- 'new_grouping'

DimPlot(object=is018.grouped, reduction='tsne', label=T)

cluster.markers.1 <- FindMarkers(is018.grouped,1,thresh.use = 2,test.use = "roc")
cluster.markers.2 <- FindMarkers(is018.grouped,2,thresh.use = 2,test.use = "roc")
cluster.markers.3 <- FindMarkers(is018.grouped,3,thresh.use = 2,test.use = "roc")
cluster.markers.4 <- FindMarkers(is018.grouped,4,thresh.use = 2,test.use = "roc")
cluster.markers.5 <- FindMarkers(is018.grouped,5,thresh.use = 2,test.use = "roc")
cluster.markers.6 <- FindMarkers(is018.grouped,6,thresh.use = 2,test.use = "roc")

cluster.markers.1$cluster <- 1
cluster.markers.2$cluster <- 2
cluster.markers.3$cluster <- 3
cluster.markers.4$cluster <- 4
cluster.markers.5$cluster <- 5
cluster.markers.6$cluster <- 6

cluster.markers.1 <- as.data.frame(cluster.markers.1)
cluster.markers.1$GeneId <- row.names(cluster.markers.1)

cluster.markers.2 <- as.data.frame(cluster.markers.2)
cluster.markers.2$GeneId <- row.names(cluster.markers.2)

cluster.markers.3 <- as.data.frame(cluster.markers.3)
cluster.markers.3$GeneId <- row.names(cluster.markers.3)

cluster.markers.4 <- as.data.frame(cluster.markers.4)
cluster.markers.4$GeneId <- row.names(cluster.markers.4)

cluster.markers.5 <- as.data.frame(cluster.markers.5)
cluster.markers.5$GeneId <- row.names(cluster.markers.5)

cluster.markers.6 <- as.data.frame(cluster.markers.6)
cluster.markers.6$GeneId <- row.names(cluster.markers.6)

cluster.markers.1 <- as.data.table(cluster.markers.1)
cluster.markers.2 <- as.data.table(cluster.markers.2)
cluster.markers.3 <- as.data.table(cluster.markers.3)
cluster.markers.4 <- as.data.table(cluster.markers.4)
cluster.markers.5 <- as.data.table(cluster.markers.5)
cluster.markers.6 <- as.data.table(cluster.markers.6)


cluster.genes <- rbind(cluster.markers.1,cluster.markers.2,cluster.markers.3,cluster.markers.4,cluster.markers.5,cluster.markers.6)

fwrite(cluster.genes, 'Cluster.genes.afterCC.consolidated.csv')


View(cluster.genes)

is018.grouped <- BuildClusterTree(is018.grouped)

PlotClusterTree(is018.grouped)

print(head(cluster.markers.1,10))

DotPlot(is018.grouped,features = rownames(cluster.markers.1)[1:10],cex.use=4)

cluster.genes.top10 <- cluster.genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

fwrite(cluster.genes.top10, 'Cluster.genes.top10.afterCC.consolidated.fixed.csv')


is018.grouped@meta.data$new_grouping_reordered <- factor(is018.grouped@meta.data$new_grouping, levels = c(1:6))

Idents(is018.grouped) <- 'new_grouping_reordered'
DoHeatmap(is018.grouped, features = c(cluster.genes.top10$GeneId) , raster=F) + labs(title='IS018 neurons: top 10 genes per cluster after consolidating')

cortical_genes <- fread('../Neuron_subpops/cortical_neuron_genes.txt', header=F)
cortical_genes <- unique(cortical_genes)

FindMarkers(is018.nocycle, ident.1 = c(1,2,4,9), ident.2 = c(0), min.pct = 0.25)

cortical.markers.1.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(1), min.pct = 0.25)
cortical.markers.1.vs.rest  <- subset(cortical.markers.1.vs.rest, row.names(cortical.markers.1.vs.rest) %in% cortical_genes$V1)
cortical.markers.1.vs.rest
DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.1.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 1 vs. rest')

cortical.markers.2.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(2), min.pct = 0.25)
cortical.markers.2.vs.rest  <- subset(cortical.markers.2.vs.rest, row.names(cortical.markers.2.vs.rest) %in% cortical_genes$V1)
cortical.markers.2.vs.rest
DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.2.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 2 vs. rest')

cortical.markers.3.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(3), min.pct = 0.25)
cortical.markers.3.vs.rest  <- subset(cortical.markers.3.vs.rest, row.names(cortical.markers.3.vs.rest) %in% cortical_genes$V1)
cortical.markers.3.vs.rest
DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.3.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 3 vs. rest')

cortical.markers.4.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(4), min.pct = 0.25)
cortical.markers.4.vs.rest  <- subset(cortical.markers.4.vs.rest, row.names(cortical.markers.4.vs.rest) %in% cortical_genes$V1)
cortical.markers.4.vs.rest
DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.4.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 4 vs. rest')

cortical.markers.5.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(5), min.pct = 0.25)
cortical.markers.5.vs.rest  <- subset(cortical.markers.5.vs.rest, row.names(cortical.markers.5.vs.rest) %in% cortical_genes$V1)
cortical.markers.5.vs.rest
DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.5.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 5 vs. rest')

cortical.markers.6.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(6), min.pct = 0.25)
cortical.markers.6.vs.rest  <- subset(cortical.markers.6.vs.rest, row.names(cortical.markers.6.vs.rest) %in% cortical_genes$V1)
cortical.markers.6.vs.rest
DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.6.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 6 vs. rest')

cortical_genes_in_scale <- cortical_genes$V1[c(cortical_genes$V1 %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts))))]
var.cortical <- cortical_genes_in_scale[cortical_genes_in_scale %in% is018.grouped@assays$RNA@var.features]


DoHeatmap(is018.grouped, features= cortical_genes_in_scale, raster=F)+ labs(title='IS018 neurons: all cortical neuron genes')
DoHeatmap(is018.grouped, features= var.cortical, raster=F)+ labs(title='IS018 neurons: variant cortical neuron genes')

library(ComplexHeatmap)

mat <- as.data.frame(as.matrix(is018.grouped@assays$RNA@data))
mat <- subset(mat, row.names(mat) %in% var.cortical)

mat.names <- data.table(mat.names = names(mat))
mat.names$order <- 1:3875

cluster.metadata <- is018.grouped@meta.data
cluster.metadata$barcode <- row.names(cluster.metadata)
cluster.metadata <- subset(cluster.metadata , select=c('barcode', 'new_grouping_reordered'))
cluster.metadata <- as.data.table(cluster.metadata)

mat.names <- merge(mat.names, cluster.metadata, by.x='mat.names', by.y='barcode')

mat.names <- mat.names[order(order)]

all(names(mat) == mat.names$mat.names)

#names(mat) <- unlist(mat.names$new_grouping_reordered)

mat.names <- mat.names[order(new_grouping_reordered)]

mat <- subset(mat, select=mat.names$mat.names)
library(circlize)

ha = HeatmapAnnotation(cluster = mat.names$new_grouping_reordered,
                       col = list(cluster=c('1'='red','2'='orange','3'='gold','4'='green','5'='blue', '6'='purple')))
#ha = HeatmapAnnotation(df = data.frame(cluster = mat.names$new_grouping_reordered))
draw(ha, 1:3000)



draw(HeatmapAnnotation(type = c(rep("a", 5), rep("b", 5)),
                       age = sample(1:20, 10),
                       col = list(type = c("a" = "red", "b" = "blue"),
                                  age = colorRamp2(c(0, 20), c("white", "red")))))

ht <- Heatmap(as.matrix(mat), cluster_rows= TRUE, cluster_columns=FALSE, top_annotation = ha, show_column_names = F)

save(is018.grouped, file='../IS018.grouped.RData')

# subsetting to more variable genes.
var.cortical.2 <- c('EOMES', 'NEUROD1', 'NHLH1', 'INSM1', 'LMO1', 'ARID5B', 'E2F1', 'MYBL2', 'TGIF1', 'ZFP36L2', 'NFIB', 'XBP1', 'CTNNB1', 'NEUROG2', 'NEUROG1', 'TFAP2C', 'POU3F2', 'CREM', 'CREB5', 'GLI3', 'ID1', 'JUN', 'DLX2', 'DLX1', 'DLX5', 'POU3F4', 'SP8', 'ASCL1', 'TSHZ2', 'ATF4', 'ZFP36L1', 'HES1', 'SOX3', 'NR2F1', 'EMX2', 'SIX3', 'MAFB', 'ZFHX4', 'RUNX1T1', 'HMGA1', 'HES6', 'HES4', 'ID3', 'HMGB2', 'PTTG1', 'SOX2', 'PAX6', 'LMO4', 'LHX9', 'PBX3', 'SOX11', 'SOX4')

mat2 <- subset(mat, row.names(mat) %in% var.cortical.2)


ht2 <- Heatmap( as.matrix(mat2), cluster_rows= TRUE, cluster_columns=FALSE, top_annotation = ha, show_column_names = F, row_dend_reorder = c(8, 7, 42, 28, 26, 6, 30, 29, 38, 39, 34, 27, 36, 44, 46, 3, 48, 5, 25, 10, 15, 43, 12, 23, 9, 51, 33, 20, 47, 16, 49, 11, 17, 40, 37, 19, 52, 4, 2, 18, 24, 50, 13, 31, 32, 45, 1, 21, 22, 35, 41, 14 ))
ht2

#

# make heatmap of top 10 genes per cluster after grouping some clusters together---

cluster.genes.after <- rbind(cluster.markers.1, cluster.markers.2,cluster.markers.3,cluster.markers.4,cluster.markers.5,cluster.markers.6)
cluster.genes.after$GeneId <- row.names(cluster.genes.after)

library(dplyr)
cluster.genes.top10.after <- cluster.genes.after %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cluster.genes.top10.after
DoHeatmap(is018.grouped, features = c(cluster.genes.top10.after$GeneId) , raster=F) + labs(title='IS018 neurons: top 10 genes per cluster after consolidation')

# some of the gene names are missing from the scale.data and I think they were corrupted with an extra digit here and there...it was because of row.names.

genes.in.scaledata <- row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@scale.data)))
genes.in.scaledata

genes.in.scaledata[grep('SIX31',genes.in.scaledata)]

genes.in.data <- row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@data)))
genes.in.data[grep('NEUROG22', genes.in.data)]


# prepare for GSEA----
res <- fread('Cluster.genes.afterCC.consolidated.fixed.csv')
res[,abslogFC := abs(avg_logFC)]
res <- res[order(-abslogFC, -myAUC)]
res <- res[abslogFC >= 1 & myAUC >= 0.6,]
res.lists <- res %>% group_by(cluster) %>% top_n(n=100, wt = -abslogFC)

res.lists %>% group_by(cluster)

res.lists <- as.data.table(res.lists)
cat(res.lists[cluster==1,GeneId], sep='\n')
cat(res.lists[cluster==2,GeneId], sep='\n')
cat(res.lists[cluster==3,GeneId], sep='\n')
cat(res.lists[cluster==4,GeneId], sep='\n')
cat(res.lists[cluster==5,GeneId], sep='\n')
cat(res.lists[cluster==6,GeneId], sep='\n')

#add extra markers from carlos-----
extra.markers <- c('RELN', 'RASGRF2', 'RORB', 'PCP4', 'BCL11B', 'FOXP2', 'GABA', 'CNR1', 'CALB1', 'NECAB1', 'CNTC6', 'B3BGALT2', 'TLE4', 'vGLUT', 'CHRNA7', 'CUX1', 'CACNG5', 'TRIB2', 'KCNK2', 'CTGF', 'O4', 'C4orf31', 'IGSF11', 'CHRNA3', 'CPNE7', 'PCP4', 'CDH2', 'GFAP', 'INPP4B', 'KCNIP2', 'GRIK4', 'ETV1', 'PDE1A', 'CYR61', 'MAP2', 'CXCL14', 'PVRL3', 'KCNIP1', 'FAM3C', 'RPRM', 'NTNG2', 'TUBB3', 'SYT17', 'PDYN', 'TOX', 'RXFP1', 'SYT10', 'TH', 'WFS1', 'VAT1L', 'GABRA5', 'SYT6', 'TBR1', 'C1QL2', 'KIAA1456', 'KCNA1', 'TH', 'C20orf103', 'HTR2C', 'TMEM163', 'CARTPT', 'AKR1C2', 'CCK', 'AKR1C3', 'FXYD6', 'ANXA1', 'PENK', 'NPY2R', 'CACNA1E', 'OPRK1', 'KCNH4', 'PCDH17', 'SCN3B', 'SEMA3C', 'COL24A1', 'SYNPR', 'CRYM', 'ADRA2A', 'TPBG', 'CTGF', 'BEND5', 'NR4A2', 'COL6A1', 'PRSS12', 'SCN4B', 'SYT2', 'LGALS1', 'MFGE8', 'SV2C', 'SNCG')

extra.markers %in% res.lists$GeneId
# need to add in these markers in case they are good distinguishers.--------
cortical_genes_extra <- unique(c(cortical_genes$V1, extra.markers))

cortical.markers.1.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes_extra %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(1), min.pct = 0.25)
cortical.markers.1.vs.rest  <- subset(cortical.markers.1.vs.rest, row.names(cortical.markers.1.vs.rest) %in% cortical_genes_extra)
cortical.markers.1.vs.rest
#DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.1.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 1 vs. rest')

cortical.markers.2.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes_extra %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(2), min.pct = 0.25)
cortical.markers.2.vs.rest  <- subset(cortical.markers.2.vs.rest, row.names(cortical.markers.2.vs.rest) %in% cortical_genes_extra)
cortical.markers.2.vs.rest
#DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.2.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 2 vs. rest')

cortical.markers.3.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes_extra %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(3), min.pct = 0.25)
cortical.markers.3.vs.rest  <- subset(cortical.markers.3.vs.rest, row.names(cortical.markers.3.vs.rest) %in% cortical_genes_extra)
cortical.markers.3.vs.rest
#DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.3.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 3 vs. rest')

cortical.markers.4.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes_extra %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(4), min.pct = 0.25)
cortical.markers.4.vs.rest  <- subset(cortical.markers.4.vs.rest, row.names(cortical.markers.4.vs.rest) %in% cortical_genes_extra)
cortical.markers.4.vs.rest
#DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.4.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 4 vs. rest')

cortical.markers.5.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes_extra %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(5), min.pct = 0.25)
cortical.markers.5.vs.rest  <- subset(cortical.markers.5.vs.rest, row.names(cortical.markers.5.vs.rest) %in% cortical_genes_extra)
cortical.markers.5.vs.rest
#DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.5.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 5 vs. rest')

cortical.markers.6.vs.rest <- FindMarkers(is018.grouped, features = c(cortical_genes_extra %in% row.names(as.data.frame(as.matrix(is018.grouped@assays$RNA@counts)))), ident.1 = c(6), min.pct = 0.25)
cortical.markers.6.vs.rest  <- subset(cortical.markers.6.vs.rest, row.names(cortical.markers.6.vs.rest) %in% cortical_genes_extra)
cortical.markers.6.vs.rest
#DoHeatmap(is018.grouped, features = c(row.names(cortical.markers.6.vs.rest)) , raster=F) + labs(title='IS018 neurons: cortical neuron genes cluster 6 vs. rest')


# combine new output lists-----
cortical.markers.1.vs.rest$cluster <- 1
cortical.markers.2.vs.rest$cluster <- 2
cortical.markers.3.vs.rest$cluster <- 3
cortical.markers.4.vs.rest$cluster <- 4
cortical.markers.5.vs.rest$cluster <- 5
cortical.markers.6.vs.rest$cluster <- 6

cortical.markers.1.vs.rest <- as.data.frame(cortical.markers.1.vs.rest)
cortical.markers.1.vs.rest$GeneId <- row.names(cortical.markers.1.vs.rest)

cortical.markers.2.vs.rest <- as.data.frame(cortical.markers.2.vs.rest)
cortical.markers.2.vs.rest$GeneId <- row.names(cortical.markers.2.vs.rest)

cortical.markers.3.vs.rest <- as.data.frame(cortical.markers.3.vs.rest)
cortical.markers.3.vs.rest$GeneId <- row.names(cortical.markers.3.vs.rest)

cortical.markers.4.vs.rest <- as.data.frame(cortical.markers.4.vs.rest)
cortical.markers.4.vs.rest$GeneId <- row.names(cortical.markers.4.vs.rest)

cortical.markers.5.vs.rest <- as.data.frame(cortical.markers.5.vs.rest)
cortical.markers.5.vs.rest$GeneId <- row.names(cortical.markers.5.vs.rest)

cortical.markers.6.vs.rest <- as.data.frame(cortical.markers.6.vs.rest)
cortical.markers.6.vs.rest$GeneId <- row.names(cortical.markers.6.vs.rest)

cortical.markers.1.vs.rest <- as.data.table(cortical.markers.1.vs.rest)
cortical.markers.2.vs.rest <- as.data.table(cortical.markers.2.vs.rest)
cortical.markers.3.vs.rest <- as.data.table(cortical.markers.3.vs.rest)
cortical.markers.4.vs.rest <- as.data.table(cortical.markers.4.vs.rest)
cortical.markers.5.vs.rest <- as.data.table(cortical.markers.5.vs.rest)
cortical.markers.6.vs.rest <- as.data.table(cortical.markers.6.vs.rest)


cortical.genes.combined <- rbind(cortical.markers.1.vs.rest,cortical.markers.2.vs.rest,cortical.markers.3.vs.rest,cortical.markers.4.vs.rest,cortical.markers.5.vs.rest,cortical.markers.6.vs.rest)

cortical.genes.combined.subset <- cortical.genes.combined[abs(avg_logFC) > 1,]


DoHeatmap(is018.grouped, features= cortical.genes.combined.subset$GeneId , raster=F)+ labs(title='IS018 neurons: DE significant cortical neuron genes')

DoHeatmap(is018.grouped, features= extra.markers , raster=F)+ labs(title='IS018 neurons: cortical layer genes')

#
# april 29 correct cortical layer gene names-------
cortex.genes <- unique(c('RELN', 'RASGRF2', 'RORB', 'PCP4', 'BCL11B', 'FOXP2', 'ABAT', 'CNR1', 'CALB1', 'NECAB1', 'TLE4', 'SLC17A6', 'CHRNA7', 'CUX1', 'CACNG5', 'TRIB2', 'KCNK2', 'CCN2', 'FOXO4', 'NDNF', 'IGSF11', 'CHRNA3', 'CPNE7', 'PCP4', 'CDH2', 'GFAP', 'INPP4B', 'KCNIP2', 'GRIK4', 'ETV1', 'PDE1A', 'CCN1', 'MAP2', 'CXCL14', 'NECTIN3', 'KCNIP1', 'FAM3C', 'RPRM', 'NTNG2', 'TUBB3', 'SYT17', 'PDYN', 'TOX', 'RXFP1', 'SYT10', 'TH', 'WFS1', 'VAT1L', 'GABRA5', 'SYT6', 'TBR1', 'C1QL2', 'TRMT9B', 'KCNA1', 'TH', 'LAMP5', 'HTR2C', 'TMEM163', 'CARTPT', 'AKR1C2', 'CCK', 'AKR1C3', 'FXYD6', 'ANXA1', 'PENK', 'NPY2R', 'CACNA1E', 'OPRK1', 'KCNH4', 'PCDH17', 'SCN3B', 'SEMA3C', 'COL24A1', 'SYNPR', 'CRYM', 'ADRA2A', 'TPBG', 'BEND5', 'NR4A2', 'COL6A1', 'PRSS12', 'SCN4B', 'SYT2', 'LGALS1', 'MFGE8', 'SV2C', 'SNCG'))

library(sctransform)


#all.genes <- rownames(x = pbmc)
#pbmc <- ScaleData(object = pbmc, features = all.genes)

#groups <- is018.grouped


is018 <- CellCycleScoring(object = is018, s.features = s.genes, g2m.features = g2.m.genes, set.ident = TRUE)

is018.nocycle <- ScaleData(object = is018, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = TRUE, features=rownames(is018))

is018.nocycle <- RunTSNE(is018.nocycle)
TSNEPlot(is018.nocycle)

Idents(is018.nocycle) <- 'sample'

TSNEPlot(is018.nocycle, label=T)

setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS018/')
save(is018.nocycle, file='IS018.nocyle.Seurat.RData')

is018.nocycle <- FindClusters(object = is018.nocycle, resolution = 0.8, force.recalc=T)
Idents(is018.nocycle) <- 'RNA_snn_res.0.8'
TSNEPlot(is018.nocycle, label=T)

DoHeatmap(is018.nocycle, features=cortex.genes, raster=F)


layer1 <- c('RELN', 'CNR1', 'CHRNA7', 'NDNF', 'INPP4B', 'CXCL14')
layer2.3 <-unique( c('RASGRF2', 'CALB1', 'CUX1', 'IGSF11', 'KCNIP2', 'NECTIN3', 'SYT17', 'WFS1', 'C1QL2', 'LAMP5', 'CARTPT', 'CCK', 'FXYD6', 'PENK', 'CACNA1E', 'KCNH4', 'SCN3B', 'COL24A1', 'CRYM', 'TPBG', 'BEND5', 'COL6A1', 'PRSS12', 'SCN4B', 'SYT2', 'LGALS1', 'MFGE8', 'SV2C', 'SNCG'))
layer4 <-c('RORB', 'NECAB1', 'CACNG5', 'CHRNA3', 'GRIK4', 'KCNIP1', 'PDYN')
layer5 <- unique(c('PCP4', 'TRIB2', 'CPNE7', 'ETV1', 'FAM3C', 'TOX', 'VAT1L', 'TRMT9B', 'HTR2C'))
layer5.6 <- unique(c('BCL11B', 'KCNK2', 'PCP4', 'PDE1A', 'RPRM', 'RXFP1', 'GABRA5', 'KCNA1'))
layer6 <- unique(c('FOXP2', 'TLE4', 'CCN2', 'CDH2', 'CCN1', 'NTNG2', 'SYT10', 'SYT6', 'TH', 'TMEM163', 'AKR1C2', 'AKR1C3', 'ANXA1', 'NPY2R', 'OPRK1', 'PCDH17', 'SEMA3C', 'SYNPR', 'ADRA2A', 'NR4A2'))
misc <- c('ABAT', 'SLC17A6', 'FOXO4', 'GFAP', 'MAP2', 'TUBB3', 'TH', 'TBR1')

DoHeatmap(is018.nocycle, features=layer1, raster=F) + labs(title='Cortical layer 1 genes')
DoHeatmap(is018.nocycle, features=layer2.3, raster=F) + labs(title='Cortical layer 2-3 genes')
DoHeatmap(is018.nocycle, features=layer4, raster=F) + labs(title='Cortical layer 4 genes')
DoHeatmap(is018.nocycle, features=layer5, raster=F) + labs(title='Cortical layer 5 genes')
DoHeatmap(is018.nocycle, features=layer5.6, raster=F) + labs(title='Cortical layer 5-6 genes')
DoHeatmap(is018.nocycle, features=layer6, raster=F) + labs(title='Cortical layer 6 genes')
DoHeatmap(is018.nocycle, features=misc, raster=F) + labs(title='Cortical misc genes')

FeaturePlot(is018.nocycle, features=c('TBR1'), pt.size = 2)

sublist <- gtools::mixedsort(c('CXCL14', 'CALB1', 'NECTIN3', 'LAMP5', 'FXYD6', 'SCN3B', 'TPBG', 'PCP4', 'TRIB2', 'ETV1', 'TOX', 'BCL11B', 'PCP4', 'PDE1A', 'TLE4', 'CCN2', 'CDH2', 'CCN1', 'TMEM163', 'ANXA1', 'NR4A2', 'SLC17A6', 'TBR1', 'VAT1L'))

DoHeatmap(is018.nocycle, features=sublist, raster=F) + labs(title='Cortical genes')

# regrouping.


cortical.genes.AUC <- FindMarkers(is018.nocycle, features = cortex.genes, test.use = 'roc')

is018.grouped <- is018.nocycle

cells.0.1.9 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.2.5 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.3.4.6.8 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))


data.order <- data.table('meta.cells'= row.names(is018.grouped@meta.data), 'order'=1:3875, cluster.old=as.character(is018.grouped@meta.data$RNA_snn_res.0.8))
data.order[,cluster.old := as.numeric(cluster.old)]
data.order[,cluster.new := as.numeric(cluster.old)]

data.order[meta.cells %in% cells.0.1.9, cluster.new := 1]
#is018.grouped@meta.data$new_grouping <- data.order$cluster

data.order[meta.cells %in% cells.2.5, cluster.new := 2]
#is018.grouped@meta.data$new_grouping <- data.order$cluster

data.order[meta.cells %in% cells.3.4.6.8, cluster.new := 3]

data.order[cluster.old ==7,cluster.new := 4]
data.order[cluster.old ==10,cluster.new := 5]

unique(data.order$cluster.new)

all(data.order$meta.cells == row.names(is018.grouped@meta.data))
is018.grouped@meta.data$new_clusters <- data.order$cluster.new

Idents(is018.grouped) <- 'new_clusters'
TSNEPlot(is018.grouped, label=T)

cells.5 <- CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.7 <-  CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))
cells.1 <-  CellSelector(plot=DimPlot(object=is018.grouped, reduction="tsne"))

data.order[meta.cells %in% cells.5, cluster.new := 5]
data.order[meta.cells %in% cells.7, cluster.new := 4]
data.order[meta.cells %in% cells.1, cluster.new := 1]

is018.grouped@meta.data$new_clusters <- data.order$cluster.new

Idents(is018.grouped) <- 'new_clusters'
TSNEPlot(is018.grouped, label=T)

save(is018.grouped, file='IS018.grouped.allScaled.Seurat.RData')

FeaturePlot(is018.grouped, features = c('TPBG', 'CALB1', 'LAMP5','SCN3B'), pt.size = 1)
FeaturePlot(is018.grouped, features = c('CXCL14'), pt.size = 1) # strong in bottom right cluster
FeaturePlot(is018.grouped, features = c('BCL11B', 'PCP4', 'PDE1A', 'TLE4'), pt.size = 1)
FeaturePlot(is018.grouped, features = c('CCN2', 'CDH2','CCN1','TMEM163'), pt.size = 1)
FeaturePlot(is018.grouped, features = c('ANXA1','NR4A2','SLC17A6','TBR1'), pt.size = 1)
FeaturePlot(is018.grouped, features = c('MAP2'), pt.size = 1)

cluster.1.markers <- FindMarkers(is018.grouped, ident.1 = c(1), ident.2 = c(2,3,4,5))
cluster.2.markers <- FindMarkers(is018.grouped, ident.1 = c(2), ident.2 = c(1,3,4,5))
cluster.3.markers <- FindMarkers(is018.grouped, ident.1 = c(3), ident.2 = c(1,2,4,5))
cluster.4.markers <- FindMarkers(is018.grouped, ident.1 = c(4), ident.2 = c(1,2,3,5))
cluster.5.markers <- FindMarkers(is018.grouped, ident.1 = c(5), ident.2 = c(1,2,3,4))

cluster.1.markers$GeneId <- row.names(cluster.1.markers)
cluster.2.markers$GeneId <- row.names(cluster.2.markers)
cluster.3.markers$GeneId <- row.names(cluster.3.markers)
cluster.4.markers$GeneId <- row.names(cluster.4.markers)
cluster.5.markers$GeneId <- row.names(cluster.5.markers)

setwd('Neuron_subpops/Scaling_all_genes_and_CC/')
fwrite(cluster.1.markers, 'Cluster1.DE.csv')
fwrite(cluster.2.markers, 'Cluster2.DE.csv')
fwrite(cluster.3.markers, 'Cluster3.DE.csv')
fwrite(cluster.4.markers, 'Cluster4.DE.csv')
fwrite(cluster.5.markers, 'Cluster5.DE.csv')

# doing enrichr for all genes with >= 1 log2FC
# ABAT, SLC17A6

FeaturePlot(is018.grouped, features = c('ABAT'), pt.size = 1)
FeaturePlot(is018.grouped, features = c('SLC17A6'), pt.size = 1)

# may 8, categorize by glial neuron gene list----
genelist <- as.data.table(read_xlsx('/Volumes/ncatssctl/NGS_related/marker_sets/glia neuron gene list.xlsx'))
genelist
length(unique(genelist$gene))
genelist$gene[duplicated(genelist$gene)]

uniqued <- genelist[!duplicated(gene),]

is018.grouped@meta.data$new_clusters <- factor(is018.grouped@meta.data$new_clusters, levels = c(1,2,3,4,5))
Idents(is018.grouped) <- 'new_clusters'

DoHeatmap(is018.grouped, features = c(uniqued$gene) , raster=F) + labs(title='IS018 neurons: glia neuron gene list')

DoHeatmap(is018.grouped, features = c(uniqued[category=='Cortical neurons',gene]) , raster=F) + labs(title='IS018 neurons: cortical neurons')
DoHeatmap(is018.grouped, features = c(uniqued[category=='Progenitors',gene]) , raster=F) + labs(title='IS018 neurons: progenitor neurons')
DoHeatmap(is018.grouped, features = c(uniqued[category=='Glia',gene]) , raster=F) + labs(title='IS018 neurons: glial neurons')
DoHeatmap(is018.grouped, features = c(uniqued[category=='Lower cortex',gene]) , raster=F) + labs(title='IS018 neurons: lower cortex')
DoHeatmap(is018.grouped, features = c(uniqued[category=='Upper cortex',gene]) , raster=F) + labs(title='IS018 neurons: upper cortex')
DoHeatmap(is018.grouped, features = c(uniqued[category=='Other',gene]) , raster=F) + labs(title='IS018 neurons: other genes')
DoHeatmap(is018.grouped, features = c(uniqued[category=='Neural crest',gene]) , raster=F) + labs(title='IS018 neurons: neural crest')

TSNEPlot(is018.grouped)

is018.grouped.cluster1 <- SubsetData(is018.grouped, ident.use = '1')
is018.grouped.cluster1@meta.data

lower.cortical.genes <- intersect(c('SNAP25', 'GRIA2', 'CNTNAP2', 'CELF4', 'NSG2', 'SYT1', 'YWHAH', 'SNCA', 'BASP1', 'DOK6', 'RTN1', 'RUNX1T1', 'FAM49A', 'MAP1B', 'SYT4', 'B3GALT2', 'GABRB2', 'LMO3', 'SCG3', 'UCHL1', 'VAMP2', 'TMEM161B-AS1', 'LY6H', 'MAPT', 'CDKN2D', 'RAB3A'), row.names(is018.grouped.cluster1@assays$RNA@data))

upper.cortical.genes <- intersect(c('MEF2C', 'STMN2', 'NSG2', 'ARPP21', 'STMN4', 'MAPT', 'GRIN2B', 'CALM1', 'NELL2', 'SCD5', 'SATB2', 'PKIA', 'MAP1B', 'INA', 'STMN1', 'NEUROD6', 'VAMP2', 'DOK5', 'RASL11B', 'SNCA', 'R3HDM1', 'TTC9B', 'RAC3', 'CXADR', 'HN1', 'CAMK2B', 'RTN1', 'CHL1', 'NSG1', 'TUBB2A', 'GABBR2', 'RBFOX2', 'CRMP1', 'GAP43', 'UCHL1', 'CDKN2D', 'NCAM1', 'MSRA', 'GPR85', 'DAAM1'), row.names(is018.grouped.cluster1@assays$RNA@data))

is018.grouped.cluster1 <- FindNeighbors(object = is018.grouped.cluster1, dims = 1:20)
is018.grouped.cluster1 <- FindClusters(object = is018.grouped.cluster1, resolution = 0.5)
Idents(is018.grouped.cluster1)
is018.grouped.cluster1 <- RunTSNE(is018.grouped.cluster1)
TSNEPlot(is018.grouped.cluster1)

DoHeatmap(is018.grouped.cluster1, features = c(lower.cortical.genes) , raster=F) + labs(title='IS018 neurons in cluster 1: lower cortical genes')
DoHeatmap(is018.grouped.cluster1, features = c(upper.cortical.genes) , raster=F) + labs(title='IS018 neurons in cluster 1: upper cortical genes')

lower.DE <- FindMarkers(is018.grouped.cluster1, features = c(lower.cortical.genes), ident.1 = '1')
lower.DE <- FindMarkers(is018.grouped.cluster1, features = c(lower.cortical.genes), ident.1 = '0')

lower.avgs <- AverageExpression(is018.grouped.cluster1, features = lower.cortical.genes)
lower.avgs.sums <- colSums(lower.avgs$RNA)

upper.avgs <- AverageExpression(is018.grouped.cluster1, features = upper.cortical.genes)
upper.avgs.sums <- colSums(upper.avgs$RNA)

cortical.genes <- intersect(c('SOX11', 'NEUROD2', 'SOX2', 'GPM6A', 'SOX4', 'MLLT11', 'CCNI', 'SLA', 'MARCKSL1', 'DCX', 'NES'), row.names(is018.grouped.cluster1@assays$RNA@data))

cortical.genes

cortical.avgs <- AverageExpression(is018.grouped.cluster1, features = cortical.genes)
cortical.avgs.sums <- colSums(cortical.avgs$RNA)

upper.avgs.sums
lower.avgs.sums

TSNEPlot(is018.grouped)

is018.grouped@meta.data
Idents(is018.grouped) <- 'Phase'
TSNEPlot(is018.grouped)

Idents(is018.grouped) <- 'RNA_snn_res.0.8'
TSNEPlot(is018.grouped)

is018.grouped <- FindClusters(object = is018.grouped, resolution = 1)
Idents(is018.grouped) <- 'RNA_snn_res.1'
TSNEPlot(is018.grouped)

is018.grouped <- FindClusters(object = is018.grouped, resolution = 0.3)
Idents(is018.grouped) <- 'RNA_snn_res.0.3'
TSNEPlot(is018.grouped)

is018.grouped <- FindClusters(object = is018.grouped, resolution = 0.1)
Idents(is018.grouped) <- 'RNA_snn_res.0.1'
TSNEPlot(is018.grouped)

VlnPlot(is018.grouped, features = 'SOX11')

progenitor.genes <- unique(intersect(c('ANXA2', 'GYPC', 'SPARC', 'SDC2', 'CRABP2', 'NTRK2', 'CCND1', 'LGALS1', 'SERF2', 'MDK', 'VGLL3', 'S100A13', 'PDLIM7', 'ANXA5', 'PRSS23', 'RPL41', 'NPC2', 'SEC11A', 'PRDX6', 'TPM1', 'RHOC', 'NEAT1', 'RPL12', 'RPL7A', 'EEF1A1', 'RPL28', 'RPS6', 'RPL23A', 'TIMP1', 'RPL8', 'METRN', 'WLS', 'RPL27A', 'CTGF', 'RCN1', 'PFN1', 'PMP22', 'ITGB8', 'SERPINH1', 'VIM', 'NME4', 'RPS7', 'MYL12A', 'RPS20', 'RPS2', 'RPLP1', 'RAB13', 'TUBB6', 'CRNDE', 'TTYH1', 'RPL23', 'RPS19', 'RPL29', 'RPS14', 'RPL3', 'SLC25A6', 'SPATS2L', 'QPRT', 'RPL35', 'RPS18', 'CLIC1', 'RPS3', 'RPL10A', 'RPS28', 'CD63', 'PDPN', 'ACTG1', 'CCNG1', 'CD99', 'B2M', 'CHCHD10', 'RPLP0', 'RPS27L', 'COL1A2', 'PFN2', 'UBB', 'RPL37', 'CRABP1', 'RPL7', 'FSTL1', 'RPL36', 'RPL19', 'FGFR1', 'ENO1', 'RPS15', 'MYL6', 'GSTP1', 'PODXL', 'CNN3', 'GNG11', 'RPS4Y1', 'AHNAK', 'CST3', 'RPS23', 'RPL13A'), row.names(is018.grouped.cluster1@assays$RNA@data)))

progenitor.avgs <- AverageExpression(is018.grouped.cluster1, features = progenitor.genes)
progenitor.avgs.sums <- colSums(progenitor.avgs$RNA)

str(is018.grouped@reductions$tsne)
is018.grouped@reductions$tsne@cell.embeddings

tsne.blend <- as.data.frame(is018.grouped@reductions$tsne@cell.embeddings)
tsne.blend$barcode <- row.names(tsne.blend)
tsne.blend <- as.data.table(tsne.blend)

tsne.blend

glia.genes <- intersect(c('SFRP1', 'SOX2', 'C1orf61', 'FABP7', 'SLC1A3', 'SYNE2', 'PAX6', 'HMGN3', 'ID4', 'MYO10', 'DBI', 'PTN', 'QKI', 'LINC01158', 'ZFHX4', 'HES1', 'HMGB2', 'LHX2'), row.names(is018.grouped.cluster1@assays$RNA@data))

#all.gene.sets <- data.table('gene'= c(lower.cortical.genes, upper.cortical.genes, cortical.genes, progenitor.genes, glia.genes), 'category'=c(rep('Lower cortical', length(lower.cortical.genes)), rep('Upper cortical', length(upper.cortical.genes)), rep('Cortical', length(cortical.genes)), rep('Progenitors', length(progenitor.genes)), rep('Glial', length(glia.genes)) ))

is018.grouped.data <- as.data.frame(is018.grouped@assays$RNA@data)
is018.grouped.data$GeneId <- row.names(is018.grouped.data)
is018.grouped.data <- as.data.table(is018.grouped.data)
glia.avgs <- colMeans(is018.grouped.data[GeneId %in% glia.genes,-c('GeneId')])

all(names(glia.avgs) == tsne.blend$barcode)

tsne.blend[,glia:=colMeans(is018.grouped.data[GeneId %in% glia.genes,-c('GeneId')])]
tsne.blend[,lower_cortical:=colMeans(is018.grouped.data[GeneId %in% lower.cortical.genes,-c('GeneId')])]
tsne.blend[,upper_cortical:=colMeans(is018.grouped.data[GeneId %in% upper.cortical.genes,-c('GeneId')])]
tsne.blend[,cortical:=colMeans(is018.grouped.data[GeneId %in% cortical.genes,-c('GeneId')])]
tsne.blend[,progenitor:=colMeans(is018.grouped.data[GeneId %in% progenitor.genes,-c('GeneId')])]
tsne.blend[,other:=colMeans(is018.grouped.data[GeneId %in% uniqued[category=='Other',gene],-c('GeneId')])]

ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=glia)) + geom_point() + theme_bw()

Idents(is018.grouped) <- 'new_clusters'
DimPlot(is018.grouped, label=TRUE, label.size = 10)

ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=lower_cortical)) + geom_point() + theme_bw()
ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=upper_cortical)) + geom_point() + theme_bw()
ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=cortical)) + geom_point() + theme_bw()
ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=progenitor)) + geom_point() + theme_bw()
ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=other)) + geom_point() + theme_bw()

VlnPlot(is018.grouped, features = c('SLC1A3','GRIA2','BCL11B','SATB2','CRABP1','NEUROD2','NES','PAX6','SOX2','FOXG1'), pt.size = 0)

cluster1.DE <- FindMarkers(is018.grouped, ident.1 = '1')
cluster2.DE <- FindMarkers(is018.grouped, ident.1 = '2')
cluster3.DE <- FindMarkers(is018.grouped, ident.1 = '3')
cluster4.DE <- FindMarkers(is018.grouped, ident.1 = '4')
cluster5.DE <- FindMarkers(is018.grouped, ident.1 = '5')

cluster1.DE$gene <- row.names(cluster1.DE)
cluster2.DE$gene <- row.names(cluster2.DE)
cluster3.DE$gene <- row.names(cluster3.DE)
cluster4.DE$gene <- row.names(cluster4.DE)
cluster5.DE$gene <- row.names(cluster5.DE)

Idents(is018.grouped) <- 'new_clusters'

fwrite(cluster1.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/IS018.grouped.cluster1.DE.csv')
fwrite(cluster2.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/IS018.grouped.cluster2.DE.csv')
fwrite(cluster3.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/IS018.grouped.cluster3.DE.csv')
fwrite(cluster4.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/IS018.grouped.cluster4.DE.csv')
fwrite(cluster5.DE, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/IS018.grouped.cluster5.DE.csv')

cluster1.DE.sub <- subset(cluster1.DE, abs(cluster1.DE$avg_logFC) > 1)
cluster2.DE.sub <- subset(cluster2.DE, abs(cluster2.DE$avg_logFC) > 1)
cluster3.DE.sub <- subset(cluster3.DE, abs(cluster3.DE$avg_logFC) > 1)
cluster4.DE.sub <- subset(cluster4.DE, abs(cluster4.DE$avg_logFC) > 1)
cluster5.DE.sub <- subset(cluster5.DE, abs(cluster5.DE$avg_logFC) > 1)

cluster1.DE.sub$gene <- row.names(cluster1.DE.sub)
cluster2.DE.sub$gene <- row.names(cluster2.DE.sub)
cluster3.DE.sub$gene <- row.names(cluster3.DE.sub)
cluster4.DE.sub$gene <- row.names(cluster4.DE.sub)
cluster5.DE.sub$gene <- row.names(cluster5.DE.sub)

cluster1.DE.sub <- as.data.table(cluster1.DE.sub)
cluster1.DE.sub <- cluster1.DE.sub[order(-avg_logFC)]

cluster2.DE.sub <- as.data.table(cluster2.DE.sub)
cluster2.DE.sub <- cluster2.DE.sub[order(-avg_logFC)]

cluster3.DE.sub <- as.data.table(cluster3.DE.sub)
cluster3.DE.sub <- cluster3.DE.sub[order(-avg_logFC)]

cluster4.DE.sub <- as.data.table(cluster4.DE.sub)
cluster4.DE.sub <- cluster4.DE.sub[order(-avg_logFC)]

cluster5.DE.sub <- as.data.table(cluster5.DE.sub)
cluster5.DE.sub <- cluster5.DE.sub[order(-avg_logFC)]

is018.grouped@meta.data$new_clusters_refactored <- factor(is018.grouped@meta.data$new_clusters, levels=c(1,2,3,4,5))
Idents(is018.grouped) <- 'new_clusters_refactored'
save(is018.grouped, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018.grouped.allScaled.Seurat.RData')

DoHeatmap(is018.grouped, features = unique(c(cluster1.DE.sub$gene[1:5], cluster2.DE.sub$gene[1:5], cluster3.DE.sub$gene[1:5], cluster4.DE.sub$gene[1:5], cluster5.DE.sub$gene[1:5] )) , raster=F) + labs(title='IS018 neurons: top DE genes by cluster') + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


all.markers <- FindAllMarkers(is018.grouped, logfc.threshold = 1)
fwrite(all.markers, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/IS018.grouped.allmarkers.DE.csv')

all.markers.logFC2 <- subset(all.markers, abs(all.markers$avg_logFC) > 2)
all.markers.logFC2


# 12/10/19: find genes expressed in a subpopulation of IS018 cells.----
# these are unbiased clusters 1 and 2.

load('/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Scaling_all_genes_and_CC/after_regrouping/Glia_neuron_categories/tSNEs.RData')

tsne.blend

load('/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018.grouped.allScaled.Seurat.RData')
is018.grouped <- UpdateSeuratObject(is018.grouped)

is018.grouped@meta.data

is018.1_2.cells <- subset(is018.grouped@meta.data, is018.grouped@meta.data$new_clusters == 1 | is018.grouped@meta.data$new_clusters == 2)
is018.1_2.cells

is018.raw <- as.data.frame(as.matrix(is018.grouped@assays$RNA@counts))
is018.raw <- subset(is018.raw, select=c(row.names(is018.1_2.cells)))
is018.raw

all(names(is018.raw) == row.names(is018.1_2.cells))

is018.1_2 <- CreateSeuratObject(counts=is018.raw, project='IS018', assay='RNA')

finish_seurat <- function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x, features=c(row.names(as.data.frame(as.matrix(x@assays$RNA@data)))))
  x <- RunPCA(x)
  x <- FindNeighbors(x)
  x <- FindClusters(x)
  x <- RunTSNE(x)
  x <- RunUMAP(x, reduction='pca', dims=1:20, verbose = F)
  return(x)
}


is018.1_2 <- finish_seurat(is018.1_2) 
is018.1_2@meta.data
TSNEPlot(is018.1_2)

is018.1_2 <- FindClusters(is018.1_2, resolution = 0.2)

cortex.genes <- unique(c('RELN', 'RASGRF2', 'RORB', 'PCP4', 'BCL11B', 'FOXP2', 'ABAT', 'CNR1', 'CALB1', 'NECAB1', 'TLE4', 'SLC17A6', 'CHRNA7', 'CUX1', 'CACNG5', 'TRIB2', 'KCNK2', 'CCN2', 'FOXO4', 'NDNF', 'IGSF11', 'CHRNA3', 'CPNE7', 'PCP4', 'CDH2', 'GFAP', 'INPP4B', 'KCNIP2', 'GRIK4', 'ETV1', 'PDE1A', 'CCN1', 'MAP2', 'CXCL14', 'NECTIN3', 'KCNIP1', 'FAM3C', 'RPRM', 'NTNG2', 'TUBB3', 'SYT17', 'PDYN', 'TOX', 'RXFP1', 'SYT10', 'TH', 'WFS1', 'VAT1L', 'GABRA5', 'SYT6', 'TBR1', 'C1QL2', 'TRMT9B', 'KCNA1', 'TH', 'LAMP5', 'HTR2C', 'TMEM163', 'CARTPT', 'AKR1C2', 'CCK', 'AKR1C3', 'FXYD6', 'ANXA1', 'PENK', 'NPY2R', 'CACNA1E', 'OPRK1', 'KCNH4', 'PCDH17', 'SCN3B', 'SEMA3C', 'COL24A1', 'SYNPR', 'CRYM', 'ADRA2A', 'TPBG', 'BEND5', 'NR4A2', 'COL6A1', 'PRSS12', 'SCN4B', 'SYT2', 'LGALS1', 'MFGE8', 'SV2C', 'SNCG'))

Idents(is018.1_2) <- 'RNA_snn_res.0.1'

is018.1_2.markers <- FindAllMarkers(is018.1_2, features=c(cortex.genes))
View(is018.1_2.markers)
fwrite(is018.1_2.markers, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Cortical_genes_DE_clusters1-2.csv')
save(is018.1_2, file='/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018.clusters1-2.RData')

load('/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018.clusters1-2.RData')

is018.1_2.markers <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuron_subpops/Cortical_genes_DE_clusters1-2.csv')

DoHeatmap(is018.grouped, features = c(cortex.genes), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


DoHeatmap(is018.grouped, features = c('SOX11', 'NEUROD2', 'SOX2', 'GPM6A', 'SOX4', 'MLLT11', 'CCNI', 'SLA', 'MARCKSL1', 'DCX', 'NES'), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)



DoHeatmap(is018.grouped, features = c('SOX11', 'GPM6A', 'SOX4', 'MLLT11', 'CCNI', 'MARCKSL1', 'DCX', 'FXYD6', 'MAP2', 'TBR1', 'CDH2', 'SLC17A6', 'TRIB2', 'SCB3B', 'BCL11B'), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

VlnPlot(is018.1_2, features=c('SOX11', 'GPM6A', 'SOX4', 'MLLT11', 'CCNI', 'MARCKSL1', 'DCX', 'FXYD6', 'MAP2', 'TBR1', 'CDH2', 'SLC17A6', 'TRIB2', 'SCB3B', 'BCL11B'))

Idents(is018.grouped) <- 'new_clusters'
TSNEPlot(is018.grouped)

VlnPlot(is018.grouped, features=c('SOX11', 'GPM6A', 'SOX4', 'MLLT11', 'CCNI', 'MARCKSL1', 'DCX', 'FXYD6', 'MAP2', 'TBR1', 'CDH2', 'SLC17A6', 'TRIB2', 'SCB3B', 'BCL11B'))

Idents(is018.grouped)
genes.test <- unique(c(cortex.genes, 'SOX11', 'GPM6A', 'SOX4', 'MLLT11', 'CCNI', 'MARCKSL1', 'DCX', 'FXYD6', 'MAP2', 'TBR1', 'CDH2', 'SLC17A6', 'TRIB2', 'SCB3B', 'BCL11B'))

is018.cortical.DE <- FindAllMarkers(is018.grouped, features=c(genes.test %in% row.names(as.data.frame(is018.grouped@assays$RNA@data))))

is018.cortical.DE.1_2_vs_rest <- FindMarkers(is018.grouped, features=c(genes.test %in% row.names(as.data.frame(is018.grouped@assays$RNA@data))), ident.1=c('1','2'))
is018.cortical.DE.1_2_vs_rest.subset <- subset(is018.cortical.DE.1_2_vs_rest, row.names(is018.cortical.DE.1_2_vs_rest) %in% c(genes.test))
View(is018.cortical.DE.1_2_vs_rest.subset)


DoHeatmap(is018.grouped, features = c(row.names(is018.cortical.DE.1_2_vs_rest.subset)), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


is018.cortical.DE.1_2_vs_rest.dt <- is018.cortical.DE.1_2_vs_rest
is018.cortical.DE.1_2_vs_rest.dt$Gene <- row.names(is018.cortical.DE.1_2_vs_rest.dt)
is018.cortical.DE.1_2_vs_rest.dt<- as.data.table(is018.cortical.DE.1_2_vs_rest.dt )
is018.cortical.DE.1_2_vs_rest.dt <- is018.cortical.DE.1_2_vs_rest.dt[order(-avg_logFC)]

is018.cortical.DE.1_2_vs_rest.dt.low <- is018.cortical.DE.1_2_vs_rest.dt[order(avg_logFC)]

is018.grouped@meta.data$new_clusters_ref <- factor(is018.grouped@meta.data$new_clusters, levels=c(1,2,3,4,5))
Idents(is018.grouped) <- 'new_clusters'

DoHeatmap(is018.grouped, features = c(is018.cortical.DE.1_2_vs_rest.dt$Gene[1:50], is018.cortical.DE.1_2_vs_rest.dt.low$Gene[1:50]), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)
fwrite(is018.cortical.DE.1_2_vs_rest.dt, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018_clusters1_2_vs_rest_DE.csv')

# 12/16/19: make heatmap with cortical and glial genes----
library(Seurat)
library(data.table)
library(biomaRt)
library(stringr)
library(ggplot2)
load("/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018.grouped.allScaled.Seurat.RData")

is018 <- is018.grouped
rm(is018.grouped)

is018
cortex.genes <- intersect(unique(c('RELN', 'RASGRF2', 'RORB', 'PCP4', 'BCL11B', 'FOXP2', 'ABAT', 'CNR1', 'CALB1', 'NECAB1', 'TLE4', 'SLC17A6', 'CHRNA7', 'CUX1', 'CACNG5', 'TRIB2', 'KCNK2', 'CCN2', 'FOXO4', 'NDNF', 'IGSF11', 'CHRNA3', 'CPNE7', 'PCP4', 'CDH2', 'GFAP', 'INPP4B', 'KCNIP2', 'GRIK4', 'ETV1', 'PDE1A', 'CCN1', 'MAP2', 'CXCL14', 'NECTIN3', 'KCNIP1', 'FAM3C', 'RPRM', 'NTNG2', 'TUBB3', 'SYT17', 'PDYN', 'TOX', 'RXFP1', 'SYT10', 'TH', 'WFS1', 'VAT1L', 'GABRA5', 'SYT6', 'TBR1', 'C1QL2', 'TRMT9B', 'KCNA1', 'TH', 'LAMP5', 'HTR2C', 'TMEM163', 'CARTPT', 'AKR1C2', 'CCK', 'AKR1C3', 'FXYD6', 'ANXA1', 'PENK', 'NPY2R', 'CACNA1E', 'OPRK1', 'KCNH4', 'PCDH17', 'SCN3B', 'SEMA3C', 'COL24A1', 'SYNPR', 'CRYM', 'ADRA2A', 'TPBG', 'BEND5', 'NR4A2', 'COL6A1', 'PRSS12', 'SCN4B', 'SYT2', 'LGALS1', 'MFGE8', 'SV2C', 'SNCG')), row.names(is018@assays$RNA@data))


glia.genes <- intersect(c('SFRP1', 'C1orf61', 'FABP7', 'SLC1A3', 'SYNE2', 'PAX6', 'HMGN3', 'ID4', 'MYO10', 'DBI', 'PTN', 'QKI', 'LINC01158', 'ZFHX4', 'HES1', 'HMGB2', 'LHX2'), row.names(is018@assays$RNA@data))

cortical.genes <- intersect(c("SOX11","NEUROD2","GPM6A","SOX4","MLLT11","CCNI","SLA","MARCKSL1","DCX","NES"), row.names(is018@assays$RNA@data))

sox2.genes <- c('SOX2')

genelist <- data.table(GeneId =c(sox2.genes,glia.genes, cortical.genes), Category=c(rep('Neuron marker', 1),
                                                                                    rep('Glia markers',length(glia.genes)),
                                                                          rep('Cortical markers',length(cortical.genes))))


# DoHeatmap(is018, features = c(cortex.genes, glia.genes), raster=F) + 
#   scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
#   guides(color=FALSE)
# 
# 
# 
# DoHeatmap(is018, features = c(glia.genes), raster=F) + 
#   scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
#   guides(color=FALSE)
# 
# is018@meta.data$new_clusters <- factor(is018@meta.data$new_clusters, levels=c(1,2,3,4,5))
# 
# 
# 
# DotPlot(is018, features=c(glia.genes), scale.max=100) + theme(axis.text.x=element_text(angle=90, hjust=0))
# 

# avg.expr <- AverageExpression(is018, features=c(glia.genes, cortical.genes))
# avg.expr <- as.data.frame(avg.expr$RNA)
# avg.expr$GeneId <- row.names(avg.expr)
# avg.expr <- as.data.table(avg.expr)
# avg.expr 
# 
# avg.expr <- melt(avg.expr, id.vars = 'GeneId')
# names(avg.expr)[2] <- 'Cluster'
# avg.expr <- as.data.table(avg.expr)
# avg.expr
# names(avg.expr)[3] <- 'Expression'

exp <- FetchData(is018, c(genelist$GeneId))
exp$barcode <- row.names(exp)
exp <- as.data.table(exp)
all(exp$barcode == is018@meta.data$new_clusters)

meta.mini <- is018@meta.data
meta.mini <- subset(meta.mini, select=c('new_clusters'))
meta.mini$barcode <- row.names(meta.mini)
meta.mini <- as.data.table(meta.mini)
meta.mini

exp <- merge(exp, meta.mini, by='barcode')
exp

exp.melt <- melt(exp, id.vars=c('barcode','new_clusters'))

exp.melt <- exp.melt[order(new_clusters)]
perc.expr <- exp.melt %>% group_by(new_clusters, variable) %>% summarize(n_cells = n(),
                                                  n_expressed = sum(value > 0),
                                                  p_expressed = n_expressed / n_cells)

perc.expr

exp.final <- merge(exp.melt, perc.expr, by=c('new_clusters','variable'))
exp.final

range(exp.final$p_expressed)

exp.final[,new_clusters := factor(new_clusters, levels=c(5,4,3,2,1))]
exp.final <- merge(exp.final, genelist, by.x='variable', by.y='GeneId')
# fancy dot plot 01/13/20-----

ggplot(exp.final, aes(x=variable, y=new_clusters, size=p_expressed)) + geom_point(aes(color=value))+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(size='Percent expressed', fill='Average expression', color='Average expression',x='Gene', y='Cluster')+
  scale_color_viridis_c(option = 'A', direction = -1, begin=0.1, end= 0.8)+
  scale_size_continuous(range = c(5,15))+
  facet_wrap(~Category, shrink = T,scales = "free_x",switch = 'x')+
  theme(axis.text=element_text(size=15), axis.text.x = element_text(angle=90, vjust=0.1), axis.title = element_text(size=15), title=element_text(size=16),strip.text.x = element_text(size=15))

# rerunning later 04/16/20------
cortical.genelist <- c('CCNI', 'DCX', 'GPM6A', 'MARCKSL1', 'MLLT11', 'NES', 'NEUROD2', 'SLA', 'SOX11', 'SOX4')
glia.genelist <- c('C1orf61', 'DBI', 'FABP7', 'HES1', 'HMGN3', 'ID4', 'LHX2', 'MYO10', 'PAX6', 'PTN', 'QKI', 'SFRP1', 'SLC1A3', 'SYNE2', 'ZFHX4')

is018 <- is018.grouped
exp <- FetchData(is018, c(cortical.genelist, glia.genelist))
exp$barcode <- row.names(exp)
exp <- as.data.table(exp)
all(exp$barcode == row.names(is018@meta.data))

meta.mini <- is018@meta.data
meta.mini <- subset(meta.mini, select=c('new_clusters'))
meta.mini$barcode <- row.names(meta.mini)
meta.mini <- as.data.table(meta.mini)
meta.mini

exp <- merge(exp, meta.mini, by='barcode')
exp

exp.melt <- melt(exp, id.vars=c('barcode','new_clusters'))

exp.melt <- exp.melt[order(new_clusters)]
perc.expr <- exp.melt %>% group_by(new_clusters, variable) %>% summarize(n_cells = n(),
                                                                         n_expressed = sum(value > 0),
                                                                         p_expressed = n_expressed / n_cells)

perc.expr

exp.final <- merge(exp.melt, perc.expr, by=c('new_clusters','variable'))
exp.final

range(exp.final$p_expressed)

exp.final[,new_clusters := factor(new_clusters, levels=c(5,4,3,2,1))]

genelist <- data.table(GeneId = c(cortical.genelist, glia.genelist), 
                       Category=c(rep('Cortical', length(cortical.genelist)),
                                  rep('Glia', length(glia.genelist))))

exp.final <- merge(exp.final, genelist, by.x='variable', by.y='GeneId')
exp.final
fwrite(exp.final, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018_dotplot_expfinal.csv')


ggplot(exp.final, aes(x=variable, y=new_clusters, size=p_expressed)) + geom_point(aes(color=value))+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(size='Percent expressed', fill='Average expression', color='Average expression',x='Gene', y='Cluster')+
  scale_color_viridis_c(option = 'A', direction = -1, begin=0.1, end= 0.8)+
  scale_size_continuous(range = c(5,15))+
  facet_wrap(~Category, shrink = T,scales = "free_x",switch = 'x')+
  theme(axis.text=element_text(size=15), axis.text.x = element_text(angle=90, vjust=0.1), axis.title = element_text(size=15), title=element_text(size=16),strip.text.x = element_text(size=15))



# 04/17/20 tsne resizing -----
TSNEPlot(is018.grouped, pt.size=1, label=TRUE,
         label.size = 10) + scale_color_manual(values=c('1'='#F8766D', '2'='#A3A500',
                                        '3'='#00BF7D','4'='#00B0F6','5'='#E76BF3'))+
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20))+
  guides(color=FALSE)

# rerunning dotplot again, april 23, 2020-----

tele <- c('BARHL2', 'MAB21L1', 'ENC1', 'MAPT', 'LHX9')
vent <- c('SST', 'GAD2', 'DLX1', 'DLX2', 'DLX5')
npc <- c('ID3', 'FABP7', 'VIM', 'PCLAF', 'SOX2', 'BIRC5', 'MDK')
pons <- c('SPARC', 'MYL9', 'SPATS2L', 'TPBG', 'SAT1', 'FAM107A', 'HSPB1', 'ZNF503', 'PDLIM2')
medu <- c('ANXA2', 'EIF5', 'AEN', 'PTGR1', 'TNFRSF12A')

is018 <- is018.grouped
exp <- FetchData(is018, c(tele, vent, npc, pons, medu))
exp$barcode <- row.names(exp)
exp <- as.data.table(exp)
all(exp$barcode == row.names(is018@meta.data))

meta.mini <- is018@meta.data
meta.mini <- subset(meta.mini, select=c('new_clusters'))
meta.mini$barcode <- row.names(meta.mini)
meta.mini <- as.data.table(meta.mini)
meta.mini

exp <- merge(exp, meta.mini, by='barcode')
exp

exp.melt <- melt(exp, id.vars=c('barcode','new_clusters'))

exp.melt <- exp.melt[order(new_clusters)]
perc.expr <- exp.melt %>% group_by(new_clusters, variable) %>% summarize(n_cells = n(),
                                                                         n_expressed = sum(value > 0),
                                                                         p_expressed = n_expressed / n_cells)

perc.expr

exp.final <- merge(exp.melt, perc.expr, by=c('new_clusters','variable'))
exp.final

range(exp.final$p_expressed)

exp.final[,new_clusters := factor(new_clusters, levels=c(5,4,3,2,1))]

genelist <- data.table(GeneId = c(tele, vent, npc, pons, medu), 
Category=c(rep('Telencephalic or dorsal forebrain neuron',
               length(tele)),
    rep('Ventral forebrain neuron', length(vent)),
    rep('Neural progenitors', length(npc)),
    rep('Pons', length(pons)),
    rep('Medulla', length(medu))
    ))

exp.final <- merge(exp.final, genelist, by.x='variable', by.y='GeneId')
exp.final
#fwrite(exp.final, '/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018_dotplot_expfinal.csv')


ggplot(exp.final, aes(x=variable, y=new_clusters, size=p_expressed)) +
  geom_point(aes(color=value))+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(size='Percent expressed', fill='Average expression',
       color='Average expression',x='Gene', y='Cluster')+
  scale_color_viridis_c(option = 'A', direction = -1, begin=0.1, end= 0.8)+
  scale_size_continuous(range = c(5,10))+
  facet_wrap(~Category, shrink = T,scales = "free_x",switch = 'x')+
  theme(axis.text=element_text(size=15), axis.text.x = element_text(angle=90,
        vjust=0.1), axis.title = element_text(size=15),
        title=element_text(size=16),strip.text.x = element_text(size=15))

is018.grouped@meta.data
is018.grouped@meta.data$new_clusters_reorder_rev <- factor(is018.grouped@meta.data$new_clusters, levels=c(1,2,3,4,5))
is018.grouped@meta.data$new_clusters_reorder_rev <- is018.grouped@meta.data$new_clusters_reorder
is018.grouped@meta.data$new_clusters_reorder_rev <- gsub(5, 'five', is018.grouped@meta.data$new_clusters_reorder_rev)
is018.grouped@meta.data$new_clusters_reorder_rev <- gsub(4, 'four', is018.grouped@meta.data$new_clusters_reorder_rev)
is018.grouped@meta.data$new_clusters_reorder_rev <- gsub(3, 'three', is018.grouped@meta.data$new_clusters_reorder_rev)
is018.grouped@meta.data$new_clusters_reorder_rev <- gsub(2, 'two', is018.grouped@meta.data$new_clusters_reorder_rev)
is018.grouped@meta.data$new_clusters_reorder_rev <- gsub(1, 'one', is018.grouped@meta.data$new_clusters_reorder_rev)

Idents(is018.grouped) <- 'new_clusters_reorder_rev'

is018.grouped@meta.data$new_clusters_reorder_rev <- factor(is018.grouped@meta.data$new_clusters_reorder_rev , levels=c('five','four','three','two','one'))

DotPlot(is018.grouped, features=c(tele, vent, npc, pons, medu))+theme(axis.text.x = element_text(angle=90))


# glutamate receptors plot may 4, 2020-----
load("/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/IS018.grouped.allScaled.Seurat.RData")

library(Seurat)

Idents(is018.grouped)

genelist <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Chromium/IS018/Analysis/Neuroreceptors.xlsx'))
genelist

DotPlot(is018.grouped, features=c(genelist$GeneId))+theme(axis.text.x = element_text(angle=90))

DoHeatmap(is018.grouped, features = c(genelist$GeneId), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

Idents(is018.grouped)

is018.grouped@meta.data$new_clusters
is018.grouped@active.ident

is018.grouped@meta.data$new_clusters_reordered <- factor(is018.grouped@meta.data$new_clusters, levels=c(1,2,3,4,5))
Idents(is018.grouped) <- 'new_clusters_reordered'

DotPlot(is018.grouped, features=c(genelist$GeneId))+theme(axis.text.x = element_text(angle=90))


avgs <- AverageExpression(is018.grouped, features=c(genelist$GeneId[genelist$GeneId %in% row.names(is018.grouped@assays$RNA@counts)] ))$RNA
avgs$rowSum <- rowSums(avgs)
avgs <- subset(avgs, rowSum >0.5)
avgs$GeneId <- row.names(avgs)
avgs <- as.data.table(avgs)
#avgs <- avgs[order(-rowSum)]
avgs <- merge(avgs, genelist, by='GeneId', all.x=T, all.y=F)
avgs

DotPlot(is018.grouped, features=c(rev(avgs$GeneId) ))+
  theme(axis.text.y = element_text(angle=0), axis.text.x=
          element_text(angle=90))





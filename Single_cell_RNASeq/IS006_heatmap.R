library(Seurat)

setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS006/')
load('Analysis/IS006.H9.Seurat.Rdata')

H9.ecto.markers <- FindMarkers(is006.h9, ident.1='h9.ecto', logfc.threshold = 1)
H9.ES.markers <- FindMarkers(is006.h9, ident.1='h9.es', logfc.threshold = 1)
H9.meso.markers <- FindMarkers(is006.h9, ident.1='h9.meso', logfc.threshold = 1)
H9.endo.markers <- FindMarkers(is006.h9, ident.1='h9.endo', logfc.threshold = 1)

setwd('Analysis/')
H9.ecto.markers <- fread('H9.ecto.DE.markers.csv')
H9.ES.markers <- fread('H9.ES.DE.markers.csv')
H9.endo.markers <- fread('H9.endo.DE.markers.csv')
H9.meso.markers <- fread('H9.meso.DE.markers.csv')

H9.ES.markers$ident1 <- 'ES'
H9.ecto.markers$ident1 <- 'ecto'
H9.endo.markers$ident1 <- 'endo'
H9.meso.markers$ident1 <- 'meso'

merged_markers <- rbind(H9.ES.markers, H9.ecto.markers, H9.endo.markers, H9.meso.markers)
merged_markers <- merged_markers[order(GeneId, -avg_logFC)]
merged_markers <- merged_markers[!duplicated(GeneId)]
merged_markers <- merged_markers[order(ident1, -avg_logFC)]
merged_markers

markers_up <- setorder(setDT(merged_markers), -avg_logFC)[, head(.SD, 10), keyby='ident1']
markers_up


is006.h9@meta.data$sample_reorder <- factor(is006.h9@meta.data$sample,
                                            levels=c('h9.es','h9.ecto','h9.endo','h9.meso'))

Idents(is006.h9) <- 'sample_reorder'

genelist <- c('POU5F1', 'TERF1', 'FOXD3-AS1', 'RARRES2', 'CD24', 'SCGB3A2', 'DPPA4', 'CLDN6', 'SFRP2', 'THY1', 'TFPI', 'NNAT', 'DLK1', 'TMSB15A', 'CLU', 'PTN', 'HES4', 'TLE4', 'SAT1', 'MALAT1', 'LEFTY2', 'LEFTY1', 'PLSCR2', 'CER1', 'CYP26A1', 'APOA2', 'PTGR1', 'NPPB', 'RHOBTB3', 'FGF17', 'DKK1', 'RBP1', 'VIM', 'S100A11', 'NKD1', 'HOXB-AS3', 'ID2', 'MUCL1', 'TNFRSF11B', 'TRDC')

plot <- DoHeatmap(is006.h9, features = c(genelist), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

pdf("IS006_H9_topDE_heatmap.pdf")
plot(plot)
dev.off()

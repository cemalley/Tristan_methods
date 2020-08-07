# on biowulf via nomachine----
library(Seurat)
library(data.table)
library(biomaRt)
library(stringr)
library(ggplot2)

setwd('/data/NCATS_ifx/iPSC/IS020')

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
  x <- CreateSeuratObject(counts = x, project = "IS020.h9")
  x@meta.data$sample <- samplename
  return(x)
}

setwd('/data/NCATS_ifx/iPSC/IS020/Countfiles_ENSG')

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
  
  fwrite(dt, paste0("/data/NCATS_ifx/iPSC/IS020/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
  
}



setwd('../Countfiles_gene_symbol/')
files <- Sys.glob('*gene_symbol.csv')

for (file in files){
  dt <- fread(file, header=T)
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  sample_id <- str_split_fixed(file, "_gene_symbol.csv",2)[1]
  
  fwrite(dt, paste0("/data/NCATS_ifx/iPSC/IS020/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
}


is020.C1 <- reformat_for_seurat(fread(files[1]), "Lonza_iPSC_manual")
is020.C2 <- reformat_for_seurat(fread(files[2]), "H9_ES_manual")
is020.C3 <- reformat_for_seurat(fread(files[3]), "H9_ES_auto")
is020.C4 <- reformat_for_seurat(fread(files[4]), "Lonza_iPSC_auto")

is020.lonza <- merge(is020.C1, is020.C4)
is020.h9 <- merge(is020.C2, is020.C3)

finish_seurat <- function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)
  Idents(x) <- 'sample'
  x <- RunPCA(x)
  x <- FindNeighbors(x)
  x <- FindClusters(x)
  x <- RunTSNE(x)
  return(x)
}

is020.lonza <- finish_seurat(is020.lonza)
is020.h9 <- finish_seurat(is020.h9)


Idents(is020.lonza) <- 'sample'
Idents(is020.h9) <- 'sample'

setwd('..')
save(is020.lonza, file='IS020.Lonza.Seurat.Rdata')
save(is020.h9, file='IS020.H9.Seurat.Rdata')


# locally, h9------
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020')
load('IS020.H9.Seurat.Rdata')

Idents(is020.h9)
PCAPlot(is020.h9)
TSNEPlot(is020.h9)

# I found no real effect of cell cycle on tsne coords.

is020.h9.markers <- FindMarkers(is020.h9, ident.1 = c('H9_ES_manual'), ident.2 = c('H9_ES_auto'), min.pct = 0.25)
is020.h9.markers
is020.h9.markers$GeneID <- row.names(is020.h9.markers)
is020.h9.markers <- is020.h9.markers[,c(6,1:5)]
fwrite(is020.h9.markers, 'IS020.H9.DE.csv', sep=',', col.names = T, row.names = F, quote=T)
is020.h9.markers <- fread('IS020.H9.DE.csv')
is020.h9.markers <- is020.h9.markers[order(-avg_logFC)]

# lonza-----
load('IS020.Lonza.Seurat.Rdata')

Idents(is020.lonza)
PCAPlot(is020.lonza)
TSNEPlot(is020.lonza)

is020.lonza.markers <- FindMarkers(is020.lonza, ident.1 = c('Lonza_iPSC_manual'), ident.2 = c('Lonza_iPSC_auto'), min.pct = 0.25)
is020.lonza.markers
is020.lonza.markers$GeneID <- row.names(is020.lonza.markers)
is020.lonza.markers <- is020.lonza.markers[,c(6,1:5)]
fwrite(is020.lonza.markers, 'IS020.Lonza.DE.csv', sep=',', col.names = T, row.names = F, quote=T)
is020.lonza.markers <- fread('IS020.Lonza.DE.csv')
is020.lonza.markers <- is020.lonza.markers[order(-avg_logFC)]

# heatmaps-----
genelist <- c('MYC', 'CD44', 'TP53', 'BAD', 'NANOG', 'RPS6', 'TXN', 'NFKBIA', 'LCK', 'H2AFX', 'CCNA2', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK3', 'MAPK1', 'STAT3', 'POU5F1', 'SOX2', 'FUT4', 'CD9', 'CD24', 'CD81', 'CASP7', 'CASP3', 'ICOS', 'MKI67', 'EPCAM', 'MAPKAPK2')

DoHeatmap(is020.lonza, features = c(genelist),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

is020.lonza.markers.targeted <- FindMarkers(is020.lonza, features = c(genelist),ident.1 = c('Lonza_iPSC_manual'), ident.2 = c('Lonza_iPSC_auto'), min.pct = 0.25) # CD9 and SOX2.

is020.h9.markers.targeted <- FindMarkers(is020.h9, features = c(genelist),ident.1 = c('H9_ES_manual'), ident.2 = c('H9_ES_auto'), min.pct = 0.25) # none pass logFC threshold.

DoHeatmap(is020.lonza, features = c(is020.lonza.markers$GeneID),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

genelist_filter <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/35473_Gene_To_Filter_PF.csv')

is020.h9.markers <- is020.h9.markers[GeneID %nin% genelist_filter$GeneSymbol,]
is020.lonza.markers <- is020.lonza.markers[GeneID %nin% genelist_filter$GeneSymbol,]

DoHeatmap(is020.lonza, features = c(is020.lonza.markers$GeneID[c(1:20,24:44)]),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

DoHeatmap(is020.h9, features = c(is020.h9.markers$GeneID),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


# cell cycle phase by cell line and by auto/manual----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020')
h9.meta <- read.table('H9.afterCC.metadata.csv')
h9.meta

lonza.meta <- read.table('Lonza.afterCC.metadata.csv')


library(ggplot2)
library(dtplyr)
library(ggthemes)

h9.meta.summary <- as.data.table(h9.meta %>%
                                   group_by(sample, Phase) %>%
                                   summarize(Count=n())%>% 
                                   mutate(perc=Count/sum(Count)))

ggplot(h9.meta.summary, aes(x = sample, y = perc, fill=Phase)) +
  geom_bar(stat='identity', position='stack') + scale_fill_manual(values = c('G1'='#ece7f2', 'G2M'='#a6bddb', 'S'='#2b8cbe'))+
  theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=13), legend.text = element_text(size=13), legend.title=element_text(size=15))+
  labs(x='Sample',y='Percentage',title='Percent cells per phase and culturing method')


meta <- rbind(h9.meta, lonza.meta)
meta.summary <- as.data.table(meta %>%
                                   group_by(sample, Phase) %>%
                                   summarize(Count=n())%>% 
                                   mutate(perc=Count/sum(Count)))
ggplot(meta.summary, aes(x = sample, y = perc, fill=Phase)) +
  geom_bar(stat='identity', position='stack') + scale_fill_manual(values = c('G1'='#ece7f2', 'G2M'='#a6bddb', 'S'='#2b8cbe'))+
  theme_bw()+theme(axis.title = element_text(size=15), axis.text.x = element_text(size=13, angle = 45, hjust=1), legend.text = element_text(size=13), legend.title=element_text(size=15), axis.text.y=element_text(size=13))+
  labs(x='Sample',y='Percentage',title='Percent cells per phase and culturing method')

# glial and cortical gene expression in H9 and Lonza separately-----
markers <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/marker_sets/glia neuron gene list.xlsx'))
markers

cortical.genes <- unlist(markers[category=='Cortical neurons',c('gene')], use.names = F)
glia.genes <- unlist(markers[category=='Glia',c('gene')], use.names = F)

FeaturePlot(is020.h9, features = c(cortical.genes[1]))
VlnPlot(is020.h9, features=c(cortical.genes))

DoHeatmap(is020.h9, features = c(cortical.genes),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

VlnPlot(is020.h9, features=c(glia.genes))

DoHeatmap(is020.h9, features = c(glia.genes),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


DotPlot(is020.h9, features=c(cortical.genes, glia.genes))



# make venn diagram for DE genes similarity----
load("/Volumes/ncatssctl/NGS_related/Chromium/IS020/IS020.H9.Seurat.Rdata")
load("/Volumes/ncatssctl/NGS_related/Chromium/IS020/IS020.Lonza.Seurat.Rdata")
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020')
h9.de <- fread('IS020.H9.DE.csv')
ipsc.de <- fread('IS020.Lonza.DE.csv') #

h9.de
ipsc.de


nrow(is020.h9@assays$RNA@counts) # 33002
nrow(is020.lonza@assays$RNA@counts) #  33002
nrow(h9.de) #15
nrow(ipsc.de) #98

library(VennDiagram)
library(Hmisc)


length(row.names(is020.h9@assays$RNA@counts)[row.names(is020.h9@assays$RNA@counts) %nin% unique(c(h9.de$GeneID, ipsc.de$GeneID))])

myvenn <- venn.diagram(
  x = list( row.names(is020.h9@assays$RNA@counts), h9.de$GeneID, ipsc.de$GeneID),
  category.names = c("All genes" , "WA09" , "LiPSC GR1.1"),
  filename = NULL,
  fill=c("skyblue","darkorchid","darkorange")
)

ggsave(myvenn, file='iPSC_vs_H9_DE_venn.svg', device='svg')

library(Seurat)
is020.merged <- merge(is020.h9, is020.lonza)
is020.merged <- NormalizeData(is020.merged)
is020.merged <- FindVariableFeatures(is020.merged)
is020.merged <- ScaleData(is020.merged, features = row.names(as.data.frame(is020.merged@assays$RNA@data)))
Idents(is020.merged) <- 'sample'
is020.merged <- RunPCA(is020.merged)
is020.merged <- FindNeighbors(is020.merged)
is020.merged <- FindClusters(is020.merged)
is020.merged <- RunTSNE(is020.merged)


is020.merged@meta.data$sample_rename <- is020.merged@meta.data$sample
unique(is020.merged@meta.data$sample_rename)
is020.merged@meta.data$sample_rename <- gsub('H9_ES_manual', 'WA09 manual', is020.merged@meta.data$sample_rename)
is020.merged@meta.data$sample_rename <- gsub('H9_ES_auto', 'WA09 auto', is020.merged@meta.data$sample_rename)
is020.merged@meta.data$sample_rename <- gsub('Lonza_iPSC_manual', 'LiPSC GR1.1 manual', is020.merged@meta.data$sample_rename)
is020.merged@meta.data$sample_rename <- gsub('Lonza_iPSC_auto', 'LiPSC GR1.1 auto', is020.merged@meta.data$sample_rename)

Idents(is020.merged) <- 'sample_rename'
TSNEPlot(is020.merged, pt.size=1)


# new diagrams----

h9.de <- h9.de[order(avg_logFC)]
ipsc.de <- ipsc.de[order(avg_logFC)]

h9.avgs <- AverageExpression(is020.h9, features=c(h9.de$GeneID))
h9.avgs$RNA$GeneID <- row.names(h9.avgs$RNA)
ipsc.avgs <- AverageExpression(is020.lonza, features=c(ipsc.de$GeneID))
ipsc.avgs$RNA$GeneID <- row.names(ipsc.avgs$RNA)

h9.de.merged <- merge(h9.de, h9.avgs$RNA, by='GeneID', all=T)
ipsc.de.merged <- merge(ipsc.de, ipsc.avgs$RNA, by='GeneID', all=T)

# volcano plots for IS020-------
load("/Volumes/ncatssctl/NGS_related/Chromium/IS020/IS020.H9.Seurat.Rdata")
load("/Volumes/ncatssctl/NGS_related/Chromium/IS020/IS020.Lonza.Seurat.Rdata")
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020')
h9.de <- fread('IS020.H9.DE.csv')
ipsc.de <- fread('IS020.Lonza.DE.csv') #

fig3.h9.up <- as.data.table(readxl::read_xlsx('Figure3_H9_DE.xlsx'), sheet='up')
fig3.h9.down <- as.data.table(readxl::read_xlsx('Figure3_H9_DE.xlsx'), sheet='down')

fig3.ipsc.up <- as.data.table(readxl::read_xlsx('Figure3_iPSC_DE.xlsx'), sheet='up')
fig3.ipsc.down <- as.data.table(readxl::read_xlsx('Figure3_iPSC_DE.xlsx'), sheet='down')

# volcano
library(scales)
library(ggthemes)
library(ggrepel)


restoplot <- na.omit(h9.de)

restoplot[,threshold:=ifelse(p_val_adj <=(1*(10^(-2))) & avg_logFC > 2, '#007F00', 
                             ifelse( p_val_adj <=(1*(10^(-2))) & avg_logFC < -2, 'red', 'gray' ))]

restoplot[,padj_reduce := ifelse(p_val_adj <= (1*(10^(-10))), (1*(10^(-10))), p_val_adj)]

right <- restoplot[GeneID %in% c('RPL7', 'MALAT1', 'BNIP3', 'DDIT4'),]
left <- restoplot[GeneID %in% c('SFRP1', 'ATP5F1E', 'SLIRP', 'RPL22L1', 'HIST1H4C', 'BEX1', 'POLR2L','DNAJC15', 'HNRNPAB', 'APOE', 'COPS9'),]

condition1 <- 'WA09 manual'
condition2 <- 'WA09 automated'

# no gene labels:
attach(restoplot)
class(restoplot$threshold)
factor(restoplot$threshold)
restoplot <- restoplot[order(-avg_logFC)]

ggplot(data=restoplot) + geom_point(aes(x=avg_logFC, y=-log10(p_val_adj)), 
                            size=2,)+
  labs(title=paste0('Volcano plot: WA09 manual vs automated'),
       x='Effect size: log2(fold-change)', y='-log10(adjusted p-value)')+
  theme_hc()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=12))+
  xlim(-.37,.37)+
  geom_text_repel(data=right, aes(x=right$avg_logFC, y=-log10(right$p_val_adj), label=right$GeneID),
          segment.size = 0.5, nudge_x = 5 + right$avg_logFC)+
  geom_text_repel(data=left, aes(x=left$avg_logFC, y=-log10(left$p_val_adj), label=left$GeneID),
                  segment.size = 0.5, nudge_x = -20 - left$avg_logFC)
                                 

restoplot[threshold =='red',] #5 down
restoplot[threshold =='#007F00',] #56 up

h9.de.redo <- FindMarkers(is020.h9, ident.1='H9_ES_manual', ident.2='H9_ES_auto', logfc.threshold = 0, min.pct=0, min.diff.pct = 0)
write.csv(h9.de.redo, 'H9.DE.allgenes.csv')

h9.de.redo <- as.data.table(readxl::read_xlsx('IS020_H9_DE.xlsx'))

fig3.h9.up <- as.data.table(readxl::read_xlsx('Figure3_H9_DE.xlsx',sheet='up'))
fig3.h9.down <- as.data.table(readxl::read_xlsx('Figure3_H9_DE.xlsx',sheet='down'))
restoplot <- na.omit(h9.de.redo)
names(restoplot)[1] <- 'GeneId'
restoplot <- as.data.table(restoplot)
restoplot[,p_val_adj_limits := ifelse(p_val_adj==0, 1*(10^(-300)), p_val_adj)]
right <- restoplot[GeneId %in% c(fig3.h9.up$GeneId),]
left <- restoplot[GeneId %in% c(fig3.h9.down$GeneId),]
restoplot[,color:=ifelse(GeneId %in% left$GeneId, 'red',
                         ifelse(GeneId %in% right$GeneId, '#007F00','gray'))]


geom_text_data <- data.table(x=c(-1,.6), y=c(0,0), label=c('Up in automated', 'Up in manual'))

plot <- ggplot(data=restoplot) + geom_point(aes(x=avg_logFC, y=-log10(p_val_adj_limits)), 
                                            size=2,color=restoplot$color)+
  labs(title=paste0('Volcano plot: WA09 manual vs automated'),
       x='Effect size: log2(fold-change)', y='-log10(adjusted p-value)')+
  theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
        plot.title = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "gray"))+
  xlim(-1,1)+ylim(0, 320)+
  geom_text_repel(data=right, aes(x=right$avg_logFC, y=-log10(right$p_val_adj_limits),
                                  label=right$GeneId))+
  geom_text_repel(data=left, aes(x=left$avg_logFC, y=-log10(left$p_val_adj_limits),
                                 label=left$GeneId))+
  geom_text(data=geom_text_data, aes(x=x, y=y, label=label), hjust=0, size=5)

plot
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020/')
pdf("IS020_volcano_H9.pdf")
plot(plot)
dev.off()

h9.de.redo <- fread('H9.DE.allgenes.csv')
h9.de.redo

avgs.h9 <- AverageExpression(is020.h9, features=c(h9.de.redo$V1))
avgs.h9 <- avgs.h9$RNA
avgs.h9$GeneId <- row.names(avgs.h9)
avgs.h9

h9.de.redo <- merge(h9.de.redo, avgs.h9, by.x='V1', by.y='GeneId', all=T)
h9.de.redo[,avg_logFC := -1*(avg_logFC)]
h9.de.redo <- h9.de.redo[order(-avg_logFC)]
h9.de.redo <- h9.de.redo[order(p_val_adj)]
fwrite(h9.de.redo,'IS020_H9_DE.csv')

#volcano plot for iPSC ----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020/')
fig3.ipsc.up <- data.table(GeneId=c('DRAP1', 'MFGE8', 'PKM', 'LITAF', 'AL591030.1', 'FTL', 'NDUFS2', 'BST2', 'STC1', 'AC093866.1'))
fig3.ipsc.down <- data.table(GeneId=c('NUTF2', 'MT1X', 'ERH', 'PMAIP1', 'SEC61G', 'PABPN1', 'PCBP2', 'C4orf48', 'THOC6', 'RPS21'))


#ipsc.de.redo <- FindMarkers(is020.lonza, ident.1='Lonza_iPSC_manual', ident.2='Lonza_iPSC_auto', logfc.threshold = 0, min.pct=0, min.diff.pct = 0)
#write.csv(ipsc.de.redo, '/data/NCATS_ifx/iPSC/IS020/Lonza.DE.allgenes.csv')

ipsc.de.redo <- read.csv('Lonza.DE.allgenes.csv', row.names = 1)
ipsc.de.redo$GeneId <- row.names(ipsc.de.redo)

avgs.ipsc <- AverageExpression(is020.lonza, features=c(ipsc.de.redo$GeneId))
avgs.ipsc <- avgs.ipsc$RNA
avgs.ipsc$GeneId <- row.names(avgs.ipsc)

ipsc.de.redo <- merge(ipsc.de.redo, avgs.ipsc, by='GeneId')
ipsc.de.redo <- as.data.table(ipsc.de.redo)
ipsc.de.redo <- ipsc.de.redo[order(-avg_logFC)]
ipsc.de.redo[,avg_logFC := -1*(avg_logFC)]
ipsc.de.redo
ipsc.de.redo <- ipsc.de.redo[order(-avg_logFC)]
ipsc.de.redo
ipsc.de.redo[,avg_logFC := -1*(avg_logFC)]
ipsc.de.redo <- ipsc.de.redo[,c('GeneId','p_val', 'avg_logFC','pct.1','pct.2','p_val_adj',
                                'Lonza_iPSC_auto','Lonza_iPSC_manual')]
ipsc.de.redo

fwrite(ipsc.de.redo, 'IS020_iPSC_DE.csv')

restoplot <- na.omit(ipsc.de.redo)
restoplot <- as.data.table(restoplot)
restoplot[,p_val_adj_limits := ifelse(p_val_adj==0, 1*(10^(-300)), p_val_adj)]
restoplot[,color:=ifelse(GeneId %in% right$GeneId, 'red',
                         ifelse(GeneId %in% left$GeneId, '#007F00','gray'))]
right <- subset(restoplot, restoplot$GeneId %in% fig3.ipsc.up$GeneId)
right <- restoplot[GeneId %in% c(fig3.ipsc.up$GeneId),]
left <- restoplot[GeneId %in% c(fig3.ipsc.down$GeneId),]

geom_text_data <- data.table(x=c(-1,.6), y=c(0,0), label=c('Up in automated', 'Up in manual'))

plot <- ggplot(data=restoplot) + geom_point(aes(x=avg_logFC, y=-log10(p_val_adj_limits)), 
                                            size=2,color=restoplot$color)+
  labs(title=paste0('Volcano plot: LiPSC GR1.1 manual vs automated'),
       x='Effect size: log2(fold-change)', y='-log10(adjusted p-value)')+
  theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
        plot.title = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "gray"))+
  xlim(-1,1)+ylim(0, 320)+
  geom_text_repel(data=right, aes(x=right$avg_logFC, y=-log10(right$p_val_adj_limits),
                                  label=right$GeneId))+
  geom_text_repel(data=left, aes(x=left$avg_logFC, y=-log10(left$p_val_adj_limits),
                                 label=left$GeneId))+
  geom_text(data=geom_text_data, aes(x=x, y=y, label=label), hjust=0, size=5)

plot


setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS020/')
pdf("IS020_volcano_iPSC.pdf")
plot(plot)
dev.off()



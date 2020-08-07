
# Helios/CyTOF data anlysis R script 04/16/2020
# Script authors: Jaroslav Slamecka, Claire Malley
# Related to manuscript:
# Title: Robotic Cell Culture Enables Standardized Cell Expansion and Multi-Lineage Functional Differentiation of Human Pluripotent Stem Cells
# Authors: Carlos A. Tristan, Pinar Ormanoglu, Pei-Hsuan Chu, Claire Malley, Jaroslav Slamecka, Vukasin M. Jovanovic, Yeliz Gedik, John Braisted, Sunil K. Mallanna1, Yu Chen2, Dorjbal Dorjsuren, Michael Iannotti, Sam Michael, Anton Simeonov, Ilyas Singeç
# National Center for Advancing Translational Sciences (NCATS), Division of Preclinical Innovation, Stem Cell Translation Laboratory (SCTL), National Institutes of Health (NIH), Rockville, MD 20850, USA
# Correspondence: ilyas.singec@nih.gov



# load libraries ----

library(flowCore)
library(openCyto)
library(ggcyto)
library(cowplot)
library(diffcyt)



# BEAD NORMALIZATION ----

# using a Shiny app
library(premessa)
normalizer_GUI()



# READ and PRE-PROCESS data ----

# read in data
# creates a flowSet which is a list of flowFrames
Helios.data = read.flowSet(path="premessa normalized", pattern="normalized.fcs")

# flowSet[[i]]@description$FILENAME contains preprended path supplied in read.flowSet to the file name
# flowSet[[i]]@description$GUID only contains the file name
# copy GUID into FILENAME, otherwise later in the workflow, "prepData" function in CATALYST will throw an error
# alternatively, make sure the FCS files are in the working directory and are declared explicitly by name to "read.flowSet" function
for (i in 1:length(Helios.data)) {Helios.data[[i]]@description$FILENAME = Helios.data[[i]]@description$GUID}



# load PHENO data ----

# pre-process and append to pData
pheno = read.delim("phenoData.txt", header=TRUE, stringsAsFactors=FALSE)
pheno$condition = factor(pheno$culture, levels=unique(pheno$culture))



# create the antibody PANEL ----

# get the panel from the flowSet
# column "antigen" contains a convenient marker name
# column "fcs_channel" contains the corresponding metal conjugate
# "marker_class" contains one of three possible values: "type", "state", "none"
# for this particular dataset, this distinction was not relevant, therefore all markers were assigned the value "type"
# the column containing the marker class is required by CATALYST version 1.10.3 for several functions to work properly
# (in CATALYST version 1.10.1, this column was optional and disabling it in "prepData" function by argument class=NULL did not break the functions)
panel = Helios.data[[1]]@parameters@data[,1:2]
panel$antigen = panel$desc
panel$fcs_channel = panel$name
panel$name = NULL
panel$marker_class = "type"
# write into file to manually modify
write.table(panel, file="panel.extracted.txt", row.names = FALSE, sep="\t")
# modify column "antigen" to give the antigens suitable names based on "desc" column
# delete the metals in the antigen names to make the names shorter
# then delete all unused channels (rows) with no antigen assigned
# keep metals with antigens/markers (proteins, DNA, live-dead) assigned, these will be used in building of the SCE object
# read the panel back in
panel = read.delim("panel.extracted.modified.txt", header=TRUE, stringsAsFactors=FALSE)



# check that all panel columns are in the flowSet object
all(panel$fcs_channel %in% colnames(Helios.data))



# GATING and visualization ----


colnames(Helios.data) # get channel names

# explore the data first
ch1 = "Event_length"
ch2 = "Ir191Di"
p = ggcyto(Helios.data, aes_string(x = ch1, y = ch2))
p = p + geom_hex(bins = 128)
p


# transformation - using flowCore ----

# create transformation function for visualization purposes
# include parameters that will be used for gating
# this is just for visualiztion, the data will be transformed later when constructing the SCE object
asinhTrans = arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
myTrans = transformList(c("Ir191Di","Ir193Di","Ce140Di","Pt195Di"), asinhTrans)

# transform the channels that will be used for gating and write into a separate object
hel = transform(Helios.data, myTrans)



# 191Ir_DNA1 vs. Ce140Di (140Ce_EQ4_beads) ----

ch1 = "Ir191Di"
ch2 = "Ce140Di"

# display the data first
p = ggcyto(hel, aes_string(x=ch1, y=ch2))
p = p + geom_hex(bins = 256)
p

# define the gate
sqrcut = matrix(c(9,11.7,11.7,9,
                  1.8,1.8,4,4),ncol=2,nrow=4)
colnames(sqrcut) = c(ch1,ch2)
pg = polygonGate(filterId="Cells", .gate=sqrcut)
pg

# display at the gate
p = ggcyto(hel, aes_string(x = ch1, y = ch2))
p = p + geom_hex(bins = 128)
p = p + geom_gate(pg)
p
ggsave(p, file="DNA.vs.EQ4.beads.pdf", width=6, height=6)

# filter using the gate
fres = filter(hel, pg) # if CATALYST package is loaded, this will not work!!!
fres
summary(fres)

# the result of polygon filtering is a logical subset
hel = Subset(hel, fres)

# display the result
p = ggcyto(hel, aes_string(x=ch1, y=ch2))
p = p + geom_hex(bins = 128)
p + geom_gate(pg)




# cells >> Ir191Di (DNA) vs. Event_length ----

ch1 = "Ir191Di"
ch2 = "Event_length"

# display the data first
p = ggcyto(hel, aes_string(x=ch1, y=ch2))
p = p + geom_hex(bins = 128)
p
p + coord_cartesian(ylim=c(10,150))

# define the gate
sqrcut = matrix(c(9.2,11.1,11.1,9.2,
                  18,18,50,50),ncol=2,nrow=4)
colnames(sqrcut) = c(ch1,ch2)
pg = polygonGate(filterId="Singlets", .gate=sqrcut)
pg

# display the gate
p = ggcyto(hel, aes_string(x=ch1, y=ch2))
p = p + geom_hex(bins = 128)
p = p + geom_gate(pg)
p
ggsave(p, file="DNA.vs.event.length.pdf", width=6, height=6)

# filter using the gate
fres = filter(hel, pg)
fres
summary(fres)

# the result of polygon filtering is a logical subset
hel = Subset(hel, fres)

# display the result
p = ggcyto(hel, aes_string(x = ch1, y = ch2))
p = p + geom_hex(bins = 128)
p + geom_gate(pg)



# singlets >> Ir191Di (DNA) vs. Pt195Di ----

ch1 = "Ir191Di"
ch2 = "Pt195Di"

# display the data first
p = ggcyto(hel, aes_string(x=ch1, y=ch2))
p = p + geom_hex(bins = 128)
p

# define function to draw ellipsoid gates for each sample separately
draw.ellipsoid.gate = function(sidx, # index of sample
                               par1, par2, par3, par4, # parameters for the covariate matrix
                               mean.x, # x and y axis means
                               mean.y) {
  # define the covariate matrix
  cov = matrix(c(par1, par2, par3, par4), ncol=2,
               dimnames=list(c(ch1, ch2), c(ch1, ch2)))
  mean = c(mean.x, mean.y)
  eg = ellipsoidGate(filterId="Live.Cells", .gate=cov, mean=mean)
  # draw gate
  p = ggcyto(hel[[sidx]], aes_string(x=ch1, y=ch2))
  p = p + geom_hex(bins = 128)
  p = p + geom_gate(eg)
  p
  ggsave(p, file=paste0("DNA.vs.Pt.sample",sidx,".pdf"), width=4, height=4)
  # filter using the gate
  fres = filter(hel[[sidx]], eg)
  print(summary(fres))
  # return the filter to use in subsetting
  return(fres)
}

# the result of polygon filtering is a logical subset
# the function will draw pdfs with the gates
fres1 = draw.ellipsoid.gate(sidx=1, par1=0.18, par2=0.18, par3=0.18, par4=0.8, mean.x=10.4, mean.y=5.5)
fres2 = draw.ellipsoid.gate(sidx=2, par1=0.18, par2=0.18, par3=0.18, par4=0.6, mean.x=10.24, mean.y=5.4)
fres3 = draw.ellipsoid.gate(sidx=3, par1=0.18, par2=0.18, par3=0.18, par4=0.8, mean.x=10.2, mean.y=5.3)
fres4 = draw.ellipsoid.gate(sidx=4, par1=0.18, par2=0.18, par3=0.18, par4=0.8, mean.x=10.3, mean.y=5.5)

# subset
hel = Subset(hel, list(`20190618_H9_auto_01_normalized.fcs`=fres1,
                       `20190618_H9_manual_01_normalized.fcs`=fres2,
                       `20190618_Lonza_auto_01_normalized.fcs`=fres3,
                       `20190618_Lonza_manual_01_normalized.fcs`=fres4))

# display the result
p = ggcyto(hel, aes_string(x=ch1, y=ch2))
p = p + geom_hex(bins = 128)
p



# CyTOF workflow - version 4 ----

# based on Nowicka et al., 2019
# works with updated CATALYST that no longer creates daFrame objects but SCE objects instead

library(CATALYST)

# construct SCE object
sce = prepData(hel,
               panel,
               pheno,
               panel_cols=list(channel="fcs_channel",
                               antigen="antigen",
                               class="marker_class"))

# plot numbers of cells
n_cells(sce)
pdf(file="n.cells.pdf", width=6, height=6)
plotCounts(sce, color_by="condition")
dev.off()

# smoothened histograms
pdf(file="expr.plot.QC.pdf", width=12, height=8)
p = plotExprs(sce, color_by = "sample_id")
p$facet$params$ncol = 6
p
dev.off()

# SUBSET SCE to remove QC parameters, non-variable or problematic parameters/markers
sce = sce[!(rownames(sce) %in% c("DNA1", "DNA2", "Live_Dead", "CD278_ICOS", "TRA_1_60")), ]

# MDS plot
pdf(file="MDS.pdf", width=6, height=6)
CATALYST::plotMDS(sce, color_by = "condition")
dev.off()

# heatmap
pdf(file="heatmap.pdf", width=10, height=2.8)
plotExprHeatmap(sce, bin_anno = TRUE, row_anno = TRUE)
dev.off()

# plot NRS - non-redundancy score
pdf(file="NRS.condition.pdf", width=10, height=6)
plotNRS(sce, color_by = "condition")
dev.off()
pdf(file="NRS.sample_id.pdf", width=10, height=6)
plotNRS(sce, color_by = "sample_id")
dev.off()



# define analysis parameters ----
# parameters for clustering, dimensionality reduction and plotting

params = list()
params$metan = 8 # for maxK argument, maxK is a number of metaclusters
params$xdim = 8 # xdim,ydim are dimensions for number of clusters, e.g. x,ydim=8 (64 clusters) or 10 (100 clusters)
params$ydim = 8
params$seed = 3579
params$dr.ncells = 8000 # number of cells for dimensionality reduction algorithms
params$font.base.size = 28

# write parameters to file for documenting
library(data.table)
fwrite(as.data.frame(params), file="params.txt", sep="\t")



# CLUSTERING ----
# FlowSOM clustering and ConsensusClusterPlus metaclustering

# cols_to_use - a character vector, specifies which antigens to use for clustering
# the default for cols_to_use (NULL) uses type_markers(x)
# must be provided if colData(x)$marker_class has not been specified

# FlowSOM clustering and ConsensusClusterPlus metaclustering
set.seed(params$seed)
sce = cluster(sce,
              features=rownames(sce),
              xdim=params$xdim, ydim=params$ydim,
              maxK=params$metan,
              seed=params$seed)


# draw global cluster heatmap
pdf(file=paste0("FlowSOM.clusters.heatmap.meta",params$metan,".pdf"), width=10, height=3.6)
plotClusterHeatmap(sce, hm2=NULL, k=paste0("meta",params$metan), m=NULL, cluster_anno=TRUE, draw_freqs=TRUE)         
dev.off()

# split_by="sample_id" is used to draw one heatmap per sample, will create a multipage pdf
pdf(file="FlowSOM.clusters.heatmap.split.sample_id.pdf", width=8, height=2.8)
plotClusterHeatmap(sce, hm2=NULL, k=paste0("meta",params$metan), m=NULL, cluster_anno=TRUE, draw_freqs=TRUE, split_by="sample_id")
dev.off()

# add a side panel with CD24
pdf(file="FlowSOM.clusters.heatmap.CD24.pdf", width=12, height=4)
plotClusterHeatmap(sce, hm2="CD24", k=paste0("meta",params$metan), draw_freqs=TRUE)
dev.off()

# add a side panel with sample abundances
pdf(file="FlowSOM.clusters.heatmap.sample-abundances.pdf", width=12, height=4)
plotClusterHeatmap(sce, hm2="abundances", k=paste0("meta",params$metan), draw_freqs=TRUE)
dev.off()

# draw expression histograms in each cluster
pdf(file="FlowSOM.clusters.exprs.pdf", width=24, height=12)
plotClusterExprs(sce, k=paste0("meta",params$metan))
dev.off()



# run t-SNE & UMAP ----
# UMAP needs package "uwot"

set.seed(params$seed)
sce = runDR(sce, dr="TSNE", cells=params$dr.ncells)
sce = runDR(sce, dr="UMAP", cells=params$dr.ncells)

# draw DR plots - UMAP and tSNE ----

# metaclusters and CD24
pdf(file=paste0("tSNE.FlowSOM.clusters.meta",params$metan,".pdf"), width=8, height=8)
plotDR(sce, "TSNE", color_by=paste0("meta",params$metan)) + theme(text=element_text(size=params$font.base.size))
dev.off()
pdf(file="tSNE.FlowSOM.clusters.CD24.pdf", width=8, height=8)
plotDR(sce, "TSNE", color_by="CD24") + theme(text=element_text(size=params$font.base.size))
dev.off()
pdf(file=paste0("UMAP.FlowSOM.clusters.meta",params$metan,".pdf"), width=8, height=8)
plotDR(sce, "UMAP", color_by=paste0("meta",params$metan)) + theme(text=element_text(size=params$font.base.size))
dev.off()
pdf(file="UMAP.FlowSOM.clusters.CD24.pdf", width=8, height=8)
plotDR(sce, "UMAP", color_by="CD24") + theme(text=element_text(size=params$font.base.size))
dev.off()

# facet per sample
pdf(file=paste0("UMAP.FlowSOM.clusters.sample_id.meta",params$metan,".pdf"), width=16, height=8)
plotDR(sce, "UMAP", color_by=paste0("meta",params$metan)) + theme(text=element_text(size=params$font.base.size)) + facet_wrap("sample_id", nrow=1)
dev.off()


# color by sample - t-SNE
pdf(file="tSNE.sample_id.pdf", width=16, height=8)
plotDR(sce, "TSNE", color_by="sample_id") + theme(text=element_text(size=params$font.base.size))
dev.off()
# color by line - t-SNE
pdf(file="tSNE.line.pdf", width=16, height=8)
plotDR(sce, "TSNE", color_by="condition") + theme(text=element_text(size=params$font.base.size))
dev.off()
# color by sample - UMAP
pdf(file="UMAP.sample_id.pdf", width=16, height=8)
plotDR(sce, "UMAP", color_by="sample_id") + theme(text=element_text(size=params$font.base.size))
dev.off()
# color by line - UMAP
pdf(file="UMAP.condition.pdf", width=16, height=8)
plotDR(sce, "UMAP", color_by="condition") + theme(text=element_text(size=params$font.base.size))
dev.off()



# FINAL plots ----

dim.red="TSNE"
dim.red="UMAP"

markers = c("Oct4","Nanog","Sox2","CD24")

# make a list of plots and use cowplot to put them together
# changing the direction of the legend somehow removes grid lines, this is likely a bug
# panel.drid.major=element_line() puts them back

# first plot - DR plot colored by metacluster
# should be a list, even though it's of length 1, so that it's easy to combine with the list of marker plots
plot1 = list()
plot1[[1]] = plotDR(sce, dim.red, color_by=paste0("meta",params$metan)) +
  facet_wrap("sample_id", ncol = 1) +
  theme(text=element_text(size=params$font.base.size),
        #strip.text=element_blank(), # comment this to reveal sample names on the plots
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.direction="horizontal", legend.position="bottom", legend.key.size=unit(0.10,"in")) +
  guides(color=guide_legend(label.position="bottom", nrow=1, override.aes=list(size=6)))

# plot DR plots colored by marker expression
plots = lapply(markers, function(x) {
  plotDR(sce, dim.red, color_by=x) +
    facet_wrap("sample_id", ncol=1) +
    theme(text=element_text(size=params$font.base.size),
          #strip.text=element_blank(), # comment to reveal sample names on the plots
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          #legend.title=element_blank(), # comment this to reveal marker names on legends
          legend.position="bottom", legend.direction="horizontal", legend.key.size=unit(0.4,"in"),
          legend.title=element_text(size=22), # comment to keep size of the text identical to base size of the theme, or use params$font.base.size
          panel.grid.major=element_line())
  })

# combine plots into a single list
plots = c(plot1,plots)

# cowplot ----
library(cowplot)
pdf(file=paste0(dim.red,".FlowSOM.clusters.sample_id.markers.",paste0(markers,collapse="."),".transposed.pdf"), width=length(plots)*4, height=15)
plot_grid(plotlist=plots, ncol=length(plots), align="hv", axis="tblr")
dev.off()

rm(dim.red,markers,plot1,plots)



# plot abundances of individual clusters
pdf(file=paste0("abundances.sample_id.meta",params$metan,".pdf"), width=3.4, height=6)
plotAbundances(sce, k=paste0("meta",params$metan), by="sample_id")
dev.off()



# capture SESSION INFO ----
writeLines(capture.output(sessionInfo()), con="session.info.txt")


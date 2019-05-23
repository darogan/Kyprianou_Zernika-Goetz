#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# TITLE
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/darogan/CTR_mzg205_0007
#
#
# Analysis Performed by Russell S. Hamilton
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
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


library("Rtsne")
library("useful")
library("dplyr")
library("Seurat")
library("reshape2")
library("ggplot2")
library("cowplot")
library("viridis")
library("matrixStats")

library(umap)

#pip install umap-learn
reticulate::py_install(packages = 'umap-learn')

baseDir <- "~/Documents/CTR-Groups/Magda_Zernica-Goetz/CTR_mz205_0007_Christos_and_Neophytos/scRNA"
setwd(baseDir)

Project <- "CTR_mz205_0007"

sample2age           <- read.table("sample2age.tab.csv", header=F, sep=",")
sample2age           <- as.data.frame(t(sample2age))
names(sample2age)    <- c("Sample", "Age")
rownames(sample2age) <- sample2age$Sample
corner(sample2age)
unique(sample2age$Age)

rpkm_tab <- read.table("GSE109071_rpkm.txt", header=T)
corner(rpkm_tab)
dim(rpkm_tab)

#rpkm_mat <- as.matrix(rpkm_tab)
#set.seed(42) # Set a seed if you want reproducible results
#tsne_out <- Rtsne(rpkm_mat, check_duplicates = FALSE) # Run TSNE
#plot(tsne_out$Y)


read_tab <- read.table("GSE109071_read.txt", header=T, row.names=1, stringsAsFactors = TRUE)
corner(read_tab)
dim(read_tab)

#read_tab$variance <- rowVars(as.matrix(read_tab))
#read_tab          <- read_tab[1:3000,]
#read_tab          <- select(read_tab,-c(variance))
#dim(read_tab)


message("Read in read file in to Seurat and run tSNE clustering to identify EPI cells")

matrix.su <- CreateSeuratObject(rpkm_tab, project = "cheng_et_al", min.cells = 3, min.features=3)
matrix.su <- NormalizeData(matrix.su, normalization.method = "LogNormalize", scale.factor = 10000)
matrix.su <- FindVariableFeatures(matrix.su, selection.method = "vst", nfeatures = 3000)
matrix.su <- ScaleData(matrix.su, features = rownames(matrix.su))
matrix.su <- RunPCA(matrix.su, features = VariableFeatures(object = matrix.su))

matrix.su <- JackStraw(matrix.su, num.replicate = 100)
matrix.su <- ScoreJackStraw(matrix.su, dims = 1:20)
#JackStrawPlot(matrix.su, dims = 1:20)
#ElbowPlot(matrix.su)

matrix.su <- FindNeighbors(matrix.su, dims = 1:5)
matrix.su <- FindClusters(matrix.su, resolution = 0.2)
matrix.su <- RunUMAP(matrix.su, dims = 1:5)

plt.dim   <- DimPlot(matrix.su, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + 
                     theme(text=element_text(size=12,  family="sans"),
                           axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), 
                           axis.title.x=element_text(size=8), axis.title.y=element_text(size=8))

plt.ftrs <- FeaturePlot(matrix.su, cols=c("blue","lightgrey","red"), features = c("Pou5f1","Bmp4","Amn"), 
                        pt.size = 0.25, reduction="umap", combine=F) 
plt.ftrs <- lapply(X = plt.ftrs, FUN = function(x) x + 
                       theme(plot.title=element_text(size = 10), legend.text=element_text(size=4),
                             axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), 
                             axis.title.x=element_text(size=8), axis.title.y=element_text(size=8),
                             legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.2, "lines") ))
plt.ftr  <- CombinePlots(plots = plt.ftrs, ncol=1)
                        
new.cluster.ids        <- c("EPI", "ExE", "EPI", "VE", "VE", "VE", "EPI", "EPI")
names(new.cluster.ids) <- levels(matrix.su)
matrix.su              <- RenameIdents(matrix.su, new.cluster.ids)
plt.dim.lab            <- DimPlot(matrix.su, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + 
                                  theme(text=element_text(size=12,  family="sans"),
                                        axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), 
                                        axis.title.x=element_text(size=8), axis.title.y=element_text(size=8))

plt.dims               <- plot_grid(plt.dim, plt.dim.lab, ncol=1, labels = c("A", "C"))
plt.dim.red            <- plot_grid(plt.dims, plt.ftr, ncol=2, labels = c("", "B"))

pdf(paste0(Project, "_scUMAPs", ".pdf"), width=5,height=5)
par(bg=NA)
plt.dim.red
dev.off()

png(paste0(Project, "_scUMAPs", ".png"), units="cm", width=15, height=15, res=300)
par(bg=NA)
plt.dim.red 
dev.off()


# TSNE
#matrix.tsne         <- as.data.frame(Embeddings(object=matrix.su, reduction.type="tsne", dims.use=1:2))
#res <- 0.2
#reso                <- paste0("res.", res)
#clust               <- matrix.su@meta.data %>% dplyr::select(contains(reso))
#colnames(clust)     <- c("Cluster")
#matrix.tsne$Cluster <- clust$Cluster
#head(matrix.tsne)


EPI.cells.use <- WhichCells(object = matrix.su, idents = 'EPI')


EPI.cells.rpkm <- rpkm_tab[, colnames(rpkm_tab) %in% EPI.cells.use]
corner(EPI.cells.rpkm)
dim(EPI.cells.rpkm)


test <- rownames(EPI.cells.rpkm[grep(pattern = "Mmp", x = rownames(EPI.cells.rpkm)),])
test

MMPs <- c("Mmp24", "Mmp9", "Mmp16", "Mmp23", "Mmp17", 
          "Mmp21","Mmp2","Mmp15", "Mmp1a","Mmp1b",
          "Mmp12", "Mmp7", "Mmp11", "Mmp19", "Mmp28",
          "Mmp14", "Mmp25", "T", "Otx2", "Nodal")


EPI.cells.rpkm.mmp <- EPI.cells.rpkm[rownames(EPI.cells.rpkm) %in% MMPs, ]
dim(EPI.cells.rpkm.mmp)
corner(EPI.cells.rpkm.mmp)


EPI.cells.rpkm.mmp.m <- melt(t(EPI.cells.rpkm.mmp))
head(EPI.cells.rpkm.mmp.m)



test.data <- as.data.frame(t(EPI.cells.rpkm.mmp))
test.data$Sample <- rownames(test.data)
head(test.data)

head(sample2age)

test.data.ann <- merge(test.data, sample2age, by="Sample")
rownames(test.data.ann) <- test.data.ann$Sample
head(test.data.ann)


get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix   <- findInterval(x, dens$x)
  iy   <- findInterval(y, dens$y)
  ii   <- cbind(ix, iy)
  return(dens$z[ii]) }


functionPlotMmpVsMarkerCorrelations <- function(test.data.ann, mmpGene, markerGene, varXmax, varYmax, varScalemax) {
  
cor5.25 <- cor.test(subset(test.data.ann[, c(markerGene, "Age")], Age==5.25)[,c(markerGene)], 
                    subset(test.data.ann[, c(mmpGene,    "Age")], Age==5.25)[,c(mmpGene)], method="pearson")
cor5.50 <- cor.test(subset(test.data.ann[, c(markerGene, "Age")], Age==5.5 )[,c(markerGene)], 
                    subset(test.data.ann[, c(mmpGene,    "Age")], Age==5.5 )[,c(mmpGene)], method="pearson")
cor6.25 <- cor.test(subset(test.data.ann[, c(markerGene, "Age")], Age==6.25)[,c(markerGene)], 
                    subset(test.data.ann[, c(mmpGene,    "Age")], Age==6.25)[,c(mmpGene)], method="pearson")
cor6.50 <- cor.test(subset(test.data.ann[, c(markerGene, "Age")], Age==6.5 )[,c(markerGene)], 
                    subset(test.data.ann[, c(mmpGene,    "Age")], Age==6.5 )[,c(mmpGene)], method="pearson")

labels <- c( "5.25"=paste("5.25 \n[ r=", signif(cor5.25$estimate,2), " p=", signif(cor5.25$p.value,2),"]"), 
             "5.5"=paste("5.5  \n[ r=",  signif(cor5.50$estimate,2), " p=", signif(cor5.50$p.value,2),"]"),
             "6.25"=paste("6.25 \n[ r=", signif(cor6.25$estimate,2), " p=", signif(cor6.25$p.value,2),"]"),
             "6.5"=paste("6.5  \n[ r=",  signif(cor6.50$estimate,2), " p=", signif(cor6.50$p.value,2),"]") )

#test.data.ann$density <- get_density(log2(test.data.ann[,c(markerGene)]+1), log2(test.data.ann[,c(mmpGene)]+1))

test.data.ann$density <- tryCatch( { get_density(log2(test.data.ann[,c("T")]+1), log2(test.data.ann[,c( MMPs[i])]+1)) },
                            error = function(error_condition) 
                                   { rep(0.01, nrow(test.data.ann)) } )

plot <-   ggplot(test.data.ann, aes(x=log2(test.data.ann[,c(markerGene)]+1), y=log2(test.data.ann[,c(mmpGene)]+1), 
                                    group=as.factor(Age), colour=density) ) +
          geom_abline(slope=1, intercept=0, colour="black", linetype = "dashed", alpha=0.5) +
          #scale_colour_viridis() +
          scale_colour_gradientn(colors=viridis_pal()(9), limits=c(0,varScalemax), 
                                 breaks=seq(0,varScalemax,0.05), na.value = "#FDE725FF") + 
          xlim(0,varXmax) + ylim(0,varYmax) +
          xlab(paste0(markerGene, " log2(RPKM+1)")) +
          ylab(paste0(mmpGene, " log2(RPKM+1)")) +
          geom_point(alpha=0.75, size=2) +
          facet_grid( ~Age, labeller = labeller(Age = labels) ) +
          coord_equal() +
          theme(text=element_text(size=10,  family="sans"), 
                axis.text.x = element_text(size=8), #angle = 45, hjust = 1),
                axis.text.y = element_text(size=8),
                strip.text.x = element_text(size=8, face="bold"),
                legend.text=element_text(size=6),
                legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.2, "lines")
            )

return(plot)
}


p1 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "T",     varXmax=10, varYmax=10, varScalemax=0.10)  
p2 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "Otx2",  varXmax=10, varYmax=10, varScalemax=0.10)  
p3 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "Nodal", varXmax=10, varYmax=10, varScalemax=0.10)  
p4 <- plot_grid(p1, p2, p3, ncol=1)
p4


for (i in seq(1, (length(MMPs)-3), 1))
{
  message(paste0("Running ", i, " and ", MMPs[i]))
  
  p1 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "T",     varXmax=10, varYmax=10, varScalemax=0.125)  
  p2 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "Otx2",  varXmax=10, varYmax=10, varScalemax=0.125)  
  p3 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "Nodal", varXmax=10, varYmax=10, varScalemax=0.125)  
  p4 <- plot_grid(p1, p2, p3, ncol=1)
  
  pdf(paste0(Project, "_", MMPs[i],".pdf"), width=7,height=6)
  par(bg=NA)
  print(p4)
  dev.off()
  
  png(paste0(Project, "_", MMPs[i],".png"), units="cm", width=16, height=15, res=180)
  par(bg=NA)
  print(p4)
  dev.off()
}





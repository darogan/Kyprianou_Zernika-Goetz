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
library("plyr")
library("Seurat")
library("reshape2")
library("ggplot2")
library("cowplot")
library("viridis")
library("matrixStats")
library("data.table")

library(umap)

#pip install umap-learn
#reticulate::py_install(packages = 'umap-learn')

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




#
# message("Read in read file in to Seurat and run tSNE clustering to identify EPI cells")
#
matrix.su <- CreateSeuratObject(rpkm_tab, project = "cheng_et_al", min.cells = 3, min.features=3)
matrix.su <- NormalizeData(matrix.su, normalization.method = "LogNormalize", scale.factor = 10000)
matrix.su <- FindVariableFeatures(matrix.su, selection.method = "vst", nfeatures = 3000)

matrix.su <- ScaleData(matrix.su, features = rownames(matrix.su) )#, do.scale=F)

matrix.su <- RunPCA(matrix.su, features = VariableFeatures(object = matrix.su))
matrix.su <- JackStraw(matrix.su, num.replicate = 100)
matrix.su <- ScoreJackStraw(matrix.su, dims = 1:20)
#JackStrawPlot(matrix.su, dims = 1:20)
#ElbowPlot(matrix.su)

matrix.su <- FindNeighbors(matrix.su, dims = 1:5)
matrix.su <- FindClusters(matrix.su, resolution = 0.2)
matrix.su <- RunUMAP(matrix.su, dims = 1:5)

plt.dim   <- DimPlot(matrix.su, reduction = "umap", label = TRUE, pt.size = 0.5) + #NoLegend() + 
                     xlab("UMAP 1") + ylab("UMAP 2") +
                     theme(aspect.ratio=1,
                           text=element_text(size=12,  family="sans"),
                           axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), 
                           axis.title.x=element_text(size=8), axis.title.y=element_text(size=8),
                           legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.2, "lines") )

plt.ftrs <- FeaturePlot(matrix.su, cols=c("blue","lightgrey","red"), features = c("Pou5f1","Bmp4","Amn"), 
                        pt.size = 0.25, reduction="umap", combine=F) 
plt.ftrs <- lapply(X = plt.ftrs, FUN = function(x) x + 
                       xlab("UMAP 1") + ylab("UMAP 2") +
                       theme(aspect.ratio=1,
                             plot.title=element_text(size = 10), legend.text=element_text(size=4),
                             axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), 
                             axis.title.x=element_text(size=8), axis.title.y=element_text(size=8),
                             legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.2, "lines") ))
plt.ftr  <- CombinePlots(plots = plt.ftrs, ncol=3)
                        
new.cluster.ids        <- c("EPI", "ExE", "EPI", "VE", "VE", "VE", "EPI", "EPI")
names(new.cluster.ids) <- levels(matrix.su)
matrix.su              <- RenameIdents(matrix.su, new.cluster.ids)
plt.dim.lab            <- DimPlot(matrix.su, reduction = "umap", label = TRUE, pt.size = 0.5) + #NoLegend() + 
                                  xlab("UMAP 1") + ylab("UMAP 2") +
                                  theme(aspect.ratio=1,
                                        text=element_text(size=12,  family="sans"),
                                        axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), 
                                        axis.title.x=element_text(size=8), axis.title.y=element_text(size=8))

plt.dims               <- plot_grid(plt.dim, plt.dim.lab, NULL, ncol=3, labels = c("A", "B", ""))
plt.dim.red            <- plot_grid(plt.dims, plt.ftr, nrow=2, labels = c("", "C"))

pdf(paste0(Project, "_scUMAPs", ".pdf"), width=10,height=5)
par(bg=NA)
plt.dim.red
dev.off()


#Anotate the Cells with Ages
#
cell.idents           <- as.data.frame(Idents(matrix.su))
cell.idents$CellNames <- rownames(cell.idents)
colnames(cell.idents) <- c("CellType", "Sample")
cell.idents           <- merge(cell.idents,sample2age, by='Sample', sort=F)
rownames(cell.idents) <- cell.idents$Sample
matrix.su             <- AddMetaData( object = matrix.su, metadata=as.character(cell.idents$Age),       col.name='Age')
matrix.su             <- AddMetaData( object = matrix.su, metadata=as.character(cell.idents$CellType), col.name='CellType')

#CTR_mz205_0007.UMAP.Mmp11.pdf

meta.data.selection <- c("UMAP_1", "UMAP_2", "Age", "CellType", "seurat_clusters", 
                         "Mmp24", "Mmp9", "Mmp16", "Mmp23", "Mmp17", "Mmp21","Mmp2","Mmp15", "Mmp1a","Mmp1b",
                         "Mmp12", "Mmp7", "Mmp11", "Mmp19", "Mmp28", "Mmp14", "Mmp25", 
                         "T", "Otx2", "Nodal", "Fam25c", "Tdgf1",
                         "Col18a1", "Col4a1", "Col9a3", "Lama1", "Lama5", "Lamb1", "Lamc1" )

matrix.umap  <- FetchData(object = matrix.su, vars = meta.data.selection )
head(matrix.umap)

table(matrix.umap$Age)
table(matrix.umap$CellType)
table(matrix.umap$Age, matrix.umap$CellType)


# Distribution of cell ages and Types

plt.age.all <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Age)) +
               geom_point(alpha=0.5, size=0.5) +
               xlab("UMAP 1") + ylab("UMAP 2") +
               scale_colour_manual("", values = c("5.25"="blue", "5.5"="red", "6.25"="purple", "6.5"="green")) +
               theme(aspect.ratio=1, text = element_text(size=8), 
                     axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                     legend.position = "bottom") +
               guides(colour = guide_legend(override.aes = list(size=2))) 

plt.age.625 <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Age)) +
               geom_point(alpha=0.5, size=0.5) +
               xlab("UMAP 1") + ylab("UMAP 2") +
               scale_colour_manual("", values = c("5.25"="grey", "5.5"="grey", "6.25"="purple", "6.5"="grey")) +
               theme(aspect.ratio=1, text = element_text(size=8), 
                     axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), 
                     legend.position = "bottom") +
               guides(colour = guide_legend(override.aes = list(size=2)))

plt.age.combo <- plot_grid(plt.age.all, plt.age.625, ncol=2)

pdf(paste0(Project, "_scUMAPs.Ages", ".pdf"), width=7.5,height=4)
par(bg=NA)
plt.age.combo
dev.off()


head(matrix.umap)
matrix.umap.m <- melt(matrix.umap, id=c("UMAP_1","UMAP_2","Age","CellType","seurat_clusters"))
head(matrix.umap.m)


plt.all.Mmps <- ggplot(matrix.umap.m[ grep("Mmp", matrix.umap.m$variable), ], aes(x=UMAP_1, y=UMAP_2, fill=value, colour=value, group=variable)) +
                geom_point(alpha=1, size=0.5) +
                xlab("UMAP 1") + ylab("UMAP 2") +
                #scale_colour_gradient(name="RPKM", low="lightgrey", high="blue") +
                #cale_colour_gradientn(colours = c("lightgrey", "grey", "red"), values = scales::rescale(c(-0.5, 0.0, 0.01, 0.025, 3.25))) +
  scale_color_gradient2(midpoint=1.0, low="grey", mid="lightgrey", high="red", space ="Lab" ) +
                facet_wrap( ~variable, ncol=6) +
                theme(aspect.ratio=1, text = element_text(size=8), 
                axis.text.x = element_text(size=6), axis.text.y = element_text(size=6) )
plt.all.Mmps

pdf(paste0(Project, "_scUMAPs.All.Mmp", ".pdf"), width=6,height=3)
par(bg=NA)
plt.all.Mmps
dev.off()

plt.main.Mmps <- ggplot( subset(matrix.umap.m, variable == "Mmp2" | variable == "Mmp14" | variable == "Mmp25"), 
                         aes(x=UMAP_1, y=UMAP_2, colour=value, group=variable)) +
                 geom_point(alpha=1, size=0.25) +
                 scale_colour_gradient(name="RPKM", low="lightgrey", high="red") +
                 facet_wrap( ~variable, nrow=1) +
                 theme(aspect.ratio=1, text = element_text(size=8), 
                       axis.text.x = element_text(size=6), axis.text.y = element_text(size=6) )
pdf(paste0(Project, "_scUMAPs.Main.Mmp", ".pdf"), width=5,height=2)
par(bg=NA)
plt.main.Mmps
dev.off()

plt.main.Mmps.625 <- ggplot() +
                     geom_point(data=subset(matrix.umap.m, CellType== "EPI" & (variable == "Mmp2" | variable == "Mmp14")), 
                                aes(x=UMAP_1, y=UMAP_2), colour="lightgrey", alpha=1, size=0.5, show.legend=F) +
                     geom_point( data=subset(matrix.umap.m, 
                                             CellType== "EPI" & Age>=6 & (variable == "Mmp2" | variable == "Mmp14")), 
                                 aes(x=UMAP_1, y=UMAP_2, colour=value, group=variable), alpha=1, size=0.5) +
                     scale_colour_gradient2(name="RPKM", low="blue", mid="blue", high="red") +
                     xlim(-2.5,NA) +  ylim(NA,2.5) + xlab("UMAP 1") + ylab("UMAP 2") +
                     facet_wrap( ~variable, nrow=2) +
                     #ggtitle("Epiblast cells Age 6.25 & 6.5") +
                     theme(aspect.ratio=1, text = element_text(size=8), plot.title =  element_text(size=8),
                           axis.text.x = element_text(size=6), axis.text.y = element_text(size=6),
                           legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.3, "lines") )
plt.main.Mmps.625

pdf(paste0(Project, "_scUMAPs.Main.Mmp.EPI.6.plus", ".pdf"), width=6,height=2.5)
par(bg=NA)
plt.main.Mmps.625
dev.off()


plt.main.Mark <- ggplot( subset(matrix.umap.m, variable == "T" | variable == "Nodal" | variable == "Otx2" | variable == "Fam25c"), 
                         aes(x=UMAP_1, y=UMAP_2, colour=value, group=variable)) +
  geom_point(alpha=1, size=0.25) +
  scale_colour_gradient(name="RPKM", low="lightgrey", high="red") +
  facet_wrap( ~variable, nrow=1) +
  theme(aspect.ratio=1, text = element_text(size=8), 
        axis.text.x = element_text(size=6), axis.text.y = element_text(size=6) )
pdf(paste0(Project, "_scUMAPs.Markers", ".pdf"), width=5,height=2)
par(bg=NA)
plt.main.Mark
dev.off()

plt.main.Mark.625 <- ggplot() +
  geom_point(data=subset(matrix.umap.m, CellType== "EPI" & (variable == "T" | variable == "Nodal" | variable == "Tdgf1")), 
             aes(x=UMAP_1, y=UMAP_2), colour="lightgrey", alpha=1, size=0.1, show.legend=F) +
  geom_point( data=subset(matrix.umap.m, 
                          CellType== "EPI" & Age>=6 & (variable == "T" | variable == "Nodal" | variable == "Tdgf1")), 
              aes(x=UMAP_1, y=UMAP_2, colour=value, group=variable), alpha=1, size=0.1) +
  scale_colour_gradient2(name="RPKM", low="blue", mid="blue", high="red") +
  xlim(-2.5,NA) +  ylim(NA,2.5) +
  xlab("") + ylab("") +
  facet_wrap( ~variable, nrow=3) +
  theme(aspect.ratio=1, text = element_text(size=8), plot.title =  element_text(size=8), line = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        legend.key.height=unit(0.5, "lines"), legend.key.width=unit(0.3, "lines") )
plt.main.Mark.625


pdf(paste0(Project, "_scUMAPs.Figure", ".pdf"), width=5,height=4)
par(bg=NA)
plot_grid(plt.main.Mmps.625, plt.main.Mark.625, ncol=2, rel_widths=c(1,.75))
dev.off()



plt.all.Lams <- ggplot( subset(matrix.umap.m, grepl("Lam", variable)), aes(x=UMAP_1, y=UMAP_2, colour=value, group=variable)) +
  geom_point(alpha=1, size=0.25) +
  scale_colour_gradient(name="RPKM", low="lightgrey", high="red") +
  facet_wrap( ~variable, ncol=6) +
  theme(aspect.ratio=1, text = element_text(size=8), 
        axis.text.x = element_text(size=6), axis.text.y = element_text(size=6) )
pdf(paste0(Project, "_scUMAPs.Laminins", ".pdf"), width=5,height=2)
par(bg=NA)
plt.all.Lams
dev.off()


plt.all.Cols <- ggplot( subset(matrix.umap.m, grepl("Col", variable)), aes(x=UMAP_1, y=UMAP_2, colour=value, group=variable)) +
  geom_point(alpha=1, size=0.25) +
  scale_colour_gradient(name="RPKM", low="lightgrey", high="red") +
  facet_wrap( ~variable, ncol=6) +
  theme(aspect.ratio=1, text = element_text(size=8), 
        axis.text.x = element_text(size=6), axis.text.y = element_text(size=6) )
pdf(paste0(Project, "_scUMAPs.Collagens", ".pdf"), width=5,height=2)
par(bg=NA)
plt.all.Cols
dev.off()









EPI.cells.use <- WhichCells(object = matrix.su, idents = 'EPI')


EPI.cells.rpkm <- rpkm_tab[, colnames(rpkm_tab) %in% EPI.cells.use]
corner(EPI.cells.rpkm)
dim(EPI.cells.rpkm)


test <- rownames(EPI.cells.rpkm[grep(pattern = "Mmp", x = rownames(EPI.cells.rpkm)),])
test

MMPs <- c("Mmp24", "Mmp9", "Mmp16", "Mmp23", "Mmp17", 
          "Mmp21","Mmp2","Mmp15", "Mmp1a","Mmp1b",
          "Mmp12", "Mmp7", "Mmp11", "Mmp19", "Mmp28",
          "Mmp14", "Mmp25", "T", "Otx2", "Nodal", "Tdgf1")


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


message(paste("", mmpGene, markerGene,  
               signif(cor5.25$estimate,2), signif(cor5.25$p.value,2),
               signif(cor5.50$estimate,2), signif(cor5.50$p.value,2), 
               signif(cor6.25$estimate,2), signif(cor6.25$p.value,2), 
               signif(cor6.50$estimate,2), signif(cor6.50$p.value,2), "", sep=" | "))

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


p1 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "T",     varXmax=12, varYmax=12, varScalemax=0.10)  
p2 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "Otx2",  varXmax=12, varYmax=12, varScalemax=0.10)  
p3 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "Nodal", varXmax=12, varYmax=12, varScalemax=0.10)  
p4 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, "Mmp14", "Tdgf1", varXmax=12, varYmax=12, varScalemax=0.10)  

pX <- plot_grid(p1, p2, p3, p4, ncol=1)
pX


for (i in seq(1, (length(MMPs)-4), 1))
{
  message(paste0("Running ", i, " and ", MMPs[i]))
  
  p1 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "T",     varXmax=12, varYmax=12, varScalemax=0.125)  
  p2 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "Otx2",  varXmax=12, varYmax=12, varScalemax=0.125)  
  p3 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "Nodal", varXmax=12, varYmax=12, varScalemax=0.125) 
  p4 <- functionPlotMmpVsMarkerCorrelations(test.data.ann, MMPs[i], "Tdgf1", varXmax=12, varYmax=12, varScalemax=0.125)  
  pX <- plot_grid(p1, p2, p3, p4, ncol=1)
  
  pdf(paste0(Project, "_", MMPs[i],".pdf"), width=7,height=6)
  par(bg=NA)
  print(pX)
  dev.off()
  
  png(paste0(Project, "_", MMPs[i],".png"), units="cm", width=16, height=18, res=180)
  par(bg=NA)
  print(pX)
  dev.off()
}




test.data.ann.m <- melt(test.data.ann)
test.data.ann.m <- test.data.ann.m[grep(pattern = "Mmp", x = test.data.ann.m$variable),]
head(test.data.ann.m)
plt.mmp.age <- ggplot(test.data.ann.m, aes(x=Age, y=log2(value+1), group=as.factor(Age), colour=Age) ) +
               geom_boxplot(outlier.alpha = 0) +
               geom_jitter(alpha=0.75, size=0.1) +
               facet_wrap(variable ~ ., ncol=3) +
               ylab("log2(RPKM+1)") +
               ggtitle("Mmp expression vs Epiblast Cell Age") +
               theme(text=element_text(size=12,  family="sans"), legend.position = 'none')

pdf(paste0(Project, "_mmp.age",".pdf"), width=5,height=10)
par(bg=NA)
print(plt.mmp.age)
dev.off()

png(paste0(Project, "_mmp.age",".png"), units="cm", width=12, height=20, res=180)
par(bg=NA)
print(plt.mmp.age)
dev.off()


# Tabulate means
t             <- ddply(test.data.ann.m, c("Age", "variable"), summarise, mean = mean(value))
t.d           <- dcast(t, Age ~ variable)
rownames(t.d) <- t.d$Age
t.d.x         <-  t(t.d[ , !(names(t.d) %in% c("Age"))])
print(round(t.d.x, digits=2) ) 


test.data.ann.m.above0 <- test.data.ann.m
test.data.ann.m.above0[test.data.ann.m.above0 == 0] <- NA
t.above0             <- ddply(test.data.ann.m.above0, c("Age", "variable"), summarise, mean = mean(value, na.rm=TRUE))
t.d.above0           <- dcast(t.above0, Age ~ variable)
rownames(t.d.above0) <- t.d.above0$Age
t.d.x.above0         <-  t(t.d.above0[ , !(names(t.d.above0) %in% c("Age"))])
print(round(t.d.x.above0, digits=2) )


message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT ")
message("+-------------------------------------------------------------------------------")


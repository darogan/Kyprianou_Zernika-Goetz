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

#source("http://bioconductor.org/biocLite.R")
#biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
#biocLite("karyoploteR")
#biocLite("BSgenome.Mmusculus.UCSC.mm9")

library("useful")
library("dplyr")
library("reshape2")
library("ggplot2")
library("cowplot")
library("viridis")
library("matrixStats")
library("xml2")

library("karyoploteR")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("biomaRt")
library("regioneR")
library("BSgenome")
library("BSgenome.Mmusculus.UCSC.mm9")
library("org.Mm.eg.db")
library("GenomicRanges")
library("GenomicFeatures")

baseDir <- "/storage/CTR-Projects/CTR_mz205/CTR_mz205_0007/ChIP/"
setwd(baseDir)
Project   <- "CTR_mz205_0007"

txdb      <- TxDb.Mmusculus.UCSC.mm9.knownGene
all.genes <- genes(txdb)

ensembl   <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host='may2012.archive.ensembl.org')
#listAttributes(ensembl)
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_id', 'chromosome_name', 
                                 'start_position', 'end_position', 'strand', 'entrezgene'), mart = ensembl)  

histone.marks <- rev(c(Smad2_3_D3EB_AC_E14="GSM1782914_Smad2_3_D3EB_AC_E14.bigwig",
                       Smad2_3_D3EB_SB_E14="GSM1782915_Smad2_3_D3EB_SB_E14.bigwig",
                       #Tcf3_D3EB_AC_E14="GSM1782918_Tcf3_D3EB_AC_E14.bigwig",
                       #Input_D3EB_AC_E14="GSM1782919_Input_D3EB_AC_E14.bigwig",
                       #Smad2_3_D0ES_E14="GSM1782920_Smad2_3_D0ES_E14.bigwig",
                       #Smad2_3_D0ES_E14="GSM1782920_Smad2_3_D0ES_E14.bigwig",
                       #Tcf3_D0ES_E14="GSM1782923_Tcf3_D0ES_E14.bigwig",
                       #Input_D0ES_E14="GSM1782924_Input_D0ES_E14.bigwig",
                       #p53_D3EB_AC_E14="GSM1782925_p53_D3EB_AC_E14.bigwig",
                       #p53_D3EB_SB_E14="GSM1782926_p53_D3EB_SB_E14.bigwig",
                       #Input_for_p53_chip_D3EB_AC="GSM1782927_Input_for_p53_chip_D3EB_AC.bigwig",
                       Smad2_3_D3EB_AC_p53WT="GSM1782928_Smad2_3_D3EB_AC_p53WT.bigwig",
                       Smad2_3_D3EB_SB_p53WT="GSM1782929_Smad2_3_D3EB_SB_p53WT.bigwig"
                       #Input_D3EB_AC_p53WT="GSM1782930_Input_D3EB_AC_p53WT.bigwig",
                       #Smad2_3_D3EB_AC_p53KO_p73sh="GSM1782931_Smad2_3_D3EB_AC_p53KO_p73sh.bigwig",
                       #Smad2_3_D3EB_SB_p53KO_p73sh="GSM1782932_Smad2_3_D3EB_SB_p53KO_p73sh.bigwig",
                       #Input_D3EB_AC_p53KO_p73sh="GSM1782933_Input_D3EB_AC_p53KO_p73sh.bigwig" 
                       ))
histone.colours <- c("red", "blue", "red", "blue")
histone.ymax    <- c( 200,   200,    200,   200)

MMPs <- c("Mmp24", "Mmp9", "Mmp16", "Mmp23", "Mmp17", 
          "Mmp21","Mmp2","Mmp15", "Mmp1a","Mmp1b",
          "Mmp12", "Mmp7", "Mmp11", "Mmp19", "Mmp28",
          "Mmp14", "Mmp25" )#, "T", "Otx2", "Nodal")

ensEMBL.mmp           <- ensEMBL2id[ensEMBL2id$external_gene_id %in% MMPs, ]
rownames(ensEMBL.mmp) <- ensEMBL.mmp$external_gene_id

mm9.entrez2Symbol.x         <- org.Mm.egSYMBOL
mm9.mapped_genes            <- mappedkeys(mm9.entrez2Symbol.x)
mm9.entrez2Symbol           <- as.data.frame(mm9.entrez2Symbol.x[mm9.mapped_genes])
rownames(mm9.entrez2Symbol) <- mm9.entrez2Symbol$gene_id


left_win        <- 10000
right_win       <- 10000
tick.dist       <- 10000
mtick.dist      <- 2000

pp                <- getDefaultPlotParams(plot.type=1)
pp$leftmargin     <- 0.25
pp$topmargin      <- 15
pp$bottommargin   <- 15
pp$ideogramheight <- 5
pp$data1inmargin  <- 10

# If a subset is to be run ...
#MMPs <- c("Mmp24" , "Mmp14","Mmp2", "Mmp25"  )
for(gene in MMPs) 
  {
    gene.info <- ensEMBL.mmp[gene,]
    message( paste0(gene, " ", "chr", gene.info$chromosome_name, ",", 
                    (gene.info$start_position-left_win), ",", (gene.info$end_position+right_win) ))

    test.region    <- toGRanges( paste0("chr", gene.info$chromosome_name),
                                        (gene.info$start_position-left_win),(gene.info$end_position+right_win), 
                                        genome="mm9") 

    pdf(paste0(Project, ".ChIP.", gene, ".pdf"),width=7,height=4.75)
    par(bg=NA)

      kp <- plotKaryotype(zoom = test.region, cex=1.0, plot.params = pp, genome="mm9", cytobands = NULL)
      kpAddBaseNumbers(kp, tick.dist = tick.dist, minor.tick.dist = tick.dist, add.units = TRUE, cex=0.75, digits = 6)
      kpAddMainTitle(kp, gene, cex=1.0)

      genes.data                        <- makeGenesDataFromTxDb(txdb=TxDb.Mmusculus.UCSC.mm9.knownGene, karyoplot=kp,
                                                                 plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
      
      genes.data$genes$external_gene_id <- mm9.entrez2Symbol[genes.data$genes$gene_id,]$symbol
      gn2tn                             <- genes.data$genes$external_gene_id
      names(gn2tn)                      <- genes.data$genes$gene_id
      
      kpPlotGenes(kp, data=genes.data, gene.names=gn2tn, r0=0, r1=0.075, gene.name.cex = 0.8, plot.transcripts=T, 
                  add.transcript.names=F, add.strand.marks=TRUE, mark.height=0.4, col="darkblue")

      #ChIP data tracks
      total.tracks <- length(histone.marks)
      out.at       <- autotrack(1:length(histone.marks), total.tracks, margin = 0.3, r0=0.14, r1=1.175)

      for(i in seq_len(length(histone.marks))) 
         {
            bg.file <- paste0(baseDir, "/", histone.marks[i])
            message(paste0("+ Adding ", i, "/", length(histone.marks), ": ", histone.marks[i]))
            at <- autotrack(i, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.2)
            gr <- toGRanges(histone.marks[[i]])
            kp <- kpPlotBigWig(kp, data=bg.file, ymax=histone.ymax[i], r0=at$r0, r1=at$r1, col = histone.colours[i], border=NULL) # rainbow(total.tracks)[i]
            kpAxis(kp, ymin=0.0, ymax=histone.ymax[i], numticks = 2, r0=at$r0, r1=at$r1, cex=0.5)
            kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, cex=0.6, label.margin = 0.035, srt=0)
         }
    dev.off()
  }


message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT ")
message("+-------------------------------------------------------------------------------")

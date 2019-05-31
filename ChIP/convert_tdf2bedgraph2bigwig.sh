#!/bin/bash

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

# Need a chromosome size file
# cat Mus_musculus.NCBIM37.67.dna.chromosome.all.fa.fai | awk '{ print "chr"$1"\t"$2 } ' | sed 's/^chrMT/chrM/g' > /storage/Genomes/Mus_musculus/GRCm37/mm9.chrom.sizes
# OR
# fetchChromSizes mm9 > /storage/Genomes/Mus_musculus/GRCm37/mm9.chrom.sizes

module load bedtools2/2.25.0
module load hgdownload.cse.ucsc.edu/20190415

for i in *.tdf;
do
  echo "${i}"
  cd /storage/Software/packages/IGV_2.5.2/
  ./igvtools tdftobedgraph /storage/CTR-Projects/CTR_mz205/CTR_mz205_0007/ChIP/${i} /storage/CTR-Projects/CTR_mz205/CTR_mz205_0007/ChIP/${i/.tdf/.bedgraph} 

  cd /storage/CTR-Projects/CTR_mz205/CTR_mz205_0007/ChIP
  # https://gist.github.com/taoliu/2469050
  ./bdg2bw.sh ${i/.tdf/.bedgraph} /storage/Genomes/Mus_musculus/GRCm37/mm9.chrom.sizes

  rm ${i/.tdf/.bedgraph}

done

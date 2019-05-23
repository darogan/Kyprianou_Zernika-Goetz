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



grep "Sample_title" GSE109071_series_matrix.txt | sed 's/Mouse embryo cell //g' | sed 's/"//g' | tr '\t' ', ' | sed 's/.Sample_title//g' | sed 's/^,//g' > sample2age.tab.csv


grep "Sample_characteristics_ch1" GSE109071_series_matrix.txt | grep age | sed 's/age: embryonic day //g' | sed 's/"//g' | tr '\t' ', ' | sed 's/.Sample_characteristics_ch1//g' | sed 's/^,//g' >> sample2age.tab.csv

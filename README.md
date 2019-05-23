# Paper Title
### Paper Details


## Investigating Mmp expression in scRNA-Seq

### Step 1: Get expression for Mmp genes in published scRNA-seq data

* Cheng, S., Pei, Y., He, L., Peng, G., Reinius, B., Tam, P.P.L, Jing, N. and Deng, Q. (2019) Single-Cell {RNA}-Seq Reveals Cellular Heterogeneity of Pluripotency Transition and X Chromosome Dynamics during Early Mouse Development. Cell Reports. 26:10, 2593--2607.e3. [DOI](https://doi.org/10.1016/j.celrep.2019.02.031)

Download supplimental files from [GSE109071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109071)

* [GSE109071_read.txt](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109071&format=file&file=GSE109071%5Fread%2Etxt%2Egz)
* [GSE109071_rpkm.txt](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109071&format=file&file=GSE109071%5Frpkm%2Etxt%2Egz)
* [GSE109071_series_matrix.txt](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109071/matrix/GSE109071_series_matrix.txt.gz)

Use bash script [make_sample2Age_table.sh](scRNA/make_sample2Age_table.sh) to make a look up table of cell idetifiers to cell ages [sample2age.tab.csv](scRNA/sample2age.tab.csv)

Using the R script [CTR_mz205_0007.scRNA.R](scRNA/CTR_mz205_0007.scRNA.R) to first cluster and plot a dimensional reduction of the scRNA-seq data. To identify the Epiblast cells (exclude the ExE (marker = _Bmp4_) and VE (marker = _Amn_) from the analysis) _Pou5f1_ was used as a marker. The figure below shows: A. Finding clusters, B. cluster identification and C. Identified clusters.

![Dimensionality Reduction](scRNA/CTR_mz205_0007_scUMAPs.png?raw=true=100x)

__Figure__  [[PDF](scRNA/CTR_mz205_0007_scUMAPs.pdf)] [[PNG](scRNA/CTR_mz205_0007_scUMAPs.png)]

### Step 2: Analyse expression for Mmp genes in published scRNA-seq data

After extraction of the Epiblast cells from the scRNA-seq the next step is to plot the Mmp expression values against a selection of markers: T, Otx2 and Nodal. r values are Pearson's correlatins with p.value.

Example for Mmp14:
![Mp14 Correlations](scRNA/CTR_mz205_0007_Mmp14.png?raw=true=100x)

| Mmp   | PDF | PNG  |
| ----- | --- | ---- |
| Mmp1a | [CTR_mz205_0007_Mmp1a.png](scRNA/CTR_mz205_0007_Mmp1a.png) | [CTR_mz205_0007_Mmp1a.png](scRNA/CTR_mz205_0007_Mmp1a.pdf) |
| Mmp1b | [CTR_mz205_0007_Mmp1b.png](scRNA/CTR_mz205_0007_Mmp1b.png) | [CTR_mz205_0007_Mmp1b.png](scRNA/CTR_mz205_0007_Mmp1b.pdf) |
| __Mmp2__  | [CTR_mz205_0007_Mmp2.png](scRNA/CTR_mz205_0007_Mmp2.png)   | [CTR_mz205_0007_Mmp2.png](scRNA/CTR_mz205_0007_Mmp2.pdf) |
| Mmp3 | _missing_ | _missing_ |
| Mmp7  | [CTR_mz205_0007_Mmp7.png](scRNA/CTR_mz205_0007_Mmp7.png)   | [CTR_mz205_0007_Mmp7.png](scRNA/CTR_mz205_0007_Mmp7.pdf) |
| Mmp8 | _missing_ | _missing_ |
| Mmp9  | [CTR_mz205_0007_Mmp9.png](scRNA/CTR_mz205_0007_Mmp9.png)   | [CTR_mz205_0007_Mmp9.png](scRNA/CTR_mz205_0007_Mmp9.pdf) |
| Mmp10 | _missing_ | _missing_ |
| Mmp11 | [CTR_mz205_0007_Mmp11.png](scRNA/CTR_mz205_0007_Mmp11.png) | [CTR_mz205_0007_Mmp11.png](scRNA/CTR_mz205_0007_Mmp11.pdf) |
| Mmp12 | [CTR_mz205_0007_Mmp12.png](scRNA/CTR_mz205_0007_Mmp12.png) | [CTR_mz205_0007_Mmp12.png](scRNA/CTR_mz205_0007_Mmp12.pdf) |
| Mmp13 | _missing_ | _missing_ |
| __Mmp14__ | [CTR_mz205_0007_Mmp14.png](scRNA/CTR_mz205_0007_Mmp14.png) | [CTR_mz205_0007_Mmp14.png](scRNA/CTR_mz205_0007_Mmp14.pdf) |
| Mmp15 | [CTR_mz205_0007_Mmp15.png](scRNA/CTR_mz205_0007_Mmp15.png) | [CTR_mz205_0007_Mmp15.png](scRNA/CTR_mz205_0007_Mmp15.pdf) |
| Mmp16 | [CTR_mz205_0007_Mmp16.png](scRNA/CTR_mz205_0007_Mmp16.png) | [CTR_mz205_0007_Mmp16.png](scRNA/CTR_mz205_0007_Mmp16.pdf) |
| Mmp17 | [CTR_mz205_0007_Mmp17.png](scRNA/CTR_mz205_0007_Mmp17.png) | [CTR_mz205_0007_Mmp17.png](scRNA/CTR_mz205_0007_Mmp17.pdf) |
| Mmp19 | [CTR_mz205_0007_Mmp19.png](scRNA/CTR_mz205_0007_Mmp19.png) | [CTR_mz205_0007_Mmp19.png](scRNA/CTR_mz205_0007_Mmp19.pdf) |
| Mmp21 | [CTR_mz205_0007_Mmp21.png](scRNA/CTR_mz205_0007_Mmp21.png) | [CTR_mz205_0007_Mmp21.png](scRNA/CTR_mz205_0007_Mmp21.pdf) |
| Mmp23 | [CTR_mz205_0007_Mmp23.png](scRNA/CTR_mz205_0007_Mmp23.png) | [CTR_mz205_0007_Mmp23.png](scRNA/CTR_mz205_0007_Mmp23.pdf) |
| Mmp24 | [CTR_mz205_0007_Mmp24.png](scRNA/CTR_mz205_0007_Mmp24.png) | [CTR_mz205_0007_Mmp24.png](scRNA/CTR_mz205_0007_Mmp24.pdf) |
| __Mmp25__ | [CTR_mz205_0007_Mmp25.png](scRNA/CTR_mz205_0007_Mmp25.png) | [CTR_mz205_0007_Mmp25.png](scRNA/CTR_mz205_0007_Mmp25.pdf) |
| Mmp26 | _missing_ | _missing_ |
| Mmp27 | _missing_ | _missing_ |
| Mmp28 | [CTR_mz205_0007_Mmp28.png](CTR_mz205_0007_Mmp28.png) | [CTR_mz205_0007_Mmp28.png](CTR_mz205_0007_Mmp28.png) |

Missing Mmps: Mmp3, Mmp8, Mmp10, Mmp13, Mmp26, Mmp27


## Investigating Mmp ChIP-Seq associations

To be completed soon

## Contact

Contact Russell S. Hamilton (rsh46 -at- cam.ac.uk)

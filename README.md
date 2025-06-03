# wSIR

wSIR: ***W***eighted ***S***liced ***I***nverse ***R***egression for supervised dimension reduction of spatial transcriptomics and single cell gene expression data

This is an R package for computation of the wSIR low-dimensional embedding of gene expression data, using its spatial coordinates. It includes functions to analyse the resulting low-dimensional space, as well as the mapping from high-dimensional gene expression data to low-dimensional space. There is a function to project new gene expression data (where the spatial coordinates are unknown) into a low-dimensional space which has the ability to predict each cell's spatial location. 

For an overview of the method and examples, see the vignette at this [website](https://sydneybiox.github.io/wSIR/articles/wSIR_vignette.html).

## Installation

To install wSIR via Bioconductor: (note this won't work yet, only once wSIR is actually on Bioconductor)
```{r}
library(BiocManager)
BiocManager::install("wSIR")
```
To install wSIR from [GitHub] (https://github.com/SydneyBioX/wSIR):

```{r}
library(devtools)
install_github("SydneyBioX/wSIR")
```

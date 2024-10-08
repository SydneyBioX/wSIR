% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MouseData}
\alias{MouseData}
\alias{sample1_exprs}
\alias{sample1_coords}
\alias{sample1_cell_types}
\alias{sample2_exprs}
\alias{sample2_coords}
\alias{sample2_cell_types}
\alias{sample3_exprs}
\alias{sample3_coords}
\alias{sample3_cell_types}
\title{MouseGastrulationData}
\format{
\code{sample1_exprs} has a row for each cell in sample 1 and a
column for the
expression level of each gene. \code{sample1_coords} has a row for each
cell in sample 1
and a column for its position in each of the two spatial axes.
\code{sample1_cell_types}
a vector whose i'th entry contains the cell type of the i'th cell of sample
1.
\code{sample2_exprs} has a row for each cell in sample 2 and a column
for the expression
level of each gene. \code{sample2_coords} has a row for each cell in
sample 2 and a column
for its position in each of the two spatial axes. \code{sample2_cell_types}
a vector whose
i'th entry contains the cell type of the i'th cell of sample 2.
\code{sample3_exprs} has
a row for each cell in sample 3 and a column for the expression level of
each gene.
\code{sample3_coords} has a row for each cell in sample 3 and a column
for its position
in each of the two spatial axes. \code{sample3_cell_types} a vector,
whose i'th entry
contains the cell type of the i'th cell of sample 3.
}
\source{
Integration of spatial and single-cell transcriptomic data
elucidates mouse
organogenesis, \emph{Nature Biotechnology}, 2022.  Webpage:
\url{https://www.nature.com/articles/s41587-021-01006-2}
}
\description{
Data set consists of spatial transcriptomics data from a mouse embryo.
There are three samples, for each we have gene expression data (351 genes),
spatial
coordinates on the two-dimensional spatial plane, and cell type labels.
Sample 1 contains 19451 cells, sample 2 contains 14891 cells and sample 3
contains 23194 cells. We have randomly sampled 20\% of the cells from
each of these datasets for the vignette and examples to stop the data being
too large. After the random sampling, we have 3848 cells in sample 1, 2986
cells in sample 2 and 4607 cells in sample 3.
}
\keyword{datasets}

#' Example list with three count data sets in different formats
#'
#' A named list with three elements, each corresponding to a (real or simulated)
#' count data set. One of them is provided as a DESeqDataset, one as a count
#' data frame and one as a count matrix.
#'
#' The \code{Original} data set represents a subset of 10,000 genes and 11 cells
#' from the GSE74596 single-cell RNA-seq data set, obtained from the conquer
#' repository (http://imlspenticton.uzh.ch:3838/conquer/). The \code{Sim1} and
#' \code{Sim2} data sets similarly represent subsets of scRNA-seq data sets
#' simulated with two different simulation methods, using the real GSE74596 data
#' set as the basis for parameter estimation.

#' @format A named list with three elements, each corresponding to a (real or
#'   simulated) count data set.
#'
#' @return A named list with three elements, each corresponding to a (real or
#'   simulated) count data set.
"countsimExample_dfmat"

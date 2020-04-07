# countsimQC
[![Travis-CI Build Status](https://travis-ci.org/csoneson/countsimQC.svg?branch=master)](https://travis-ci.org/csoneson/countsimQC)
[![R build status](https://github.com/csoneson/countsimQC/workflows/R-CMD-check/badge.svg)](https://github.com/csoneson/countsimQC/actions)

`countsimQC` is an R package that provides functionality to create a 
comprehensive report comparing many different characteristics across multiple 
count data sets. One important use case is comparing one or more 
synthetic (e.g., RNA-seq) count matrices to a real count matrix, possibly the 
one based on which the synthetic data sets were generated. However, any 
collection of one or more count matrices can be visualized and compared.

If you use `countsimQC` for your work, we appreciate if you cite the 
accompanying paper:

- Soneson C and Robinson MD: [Towards unified quality verification of synthetic count data with countsimQC](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btx631/4345646/Towards-unified-quality-verification-of-synthetic). Bioinformatics 34(4):691-692 (2018).

## Installation
`countsimQC` depends on a number of other R packages. The following commands 
check whether the dependencies are available and installs them otherwise 
(note that R version >= 3.5 and Bioconductor version >= 3.8 are required in 
order to use the `BiocManager` package). If you have an older version of R
(3.4), you can still install `countsimQC` v0.5.4 (see "Releases"). Please see
the `NEWS` file for differences between versions.

```
## Install `BiocManager` if needed
if (!("BiocManager" %in% installed.packages()[, "Package"])) {
  install.packages("BiocManager")
}

## List dependencies
pkg <- c("rmarkdown", "edgeR", "DESeq2", "dplyr", "tidyr", "ggplot2", 
         "SummarizedExperiment", "genefilter", "DT", "GenomeInfoDbData",
         "caTools", "randtests", "stats", "utils", "methods")

## Check if dependencies are already installed
pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]

## If some dependency is missing, install it
if (length(pkg) > 0) {
	BiocManager::install(pkg, dependencies = TRUE, ask = FALSE)
}
```

Once all dependencies are available, `countsimQC` can be installed using 
the `BiocManager` package:

```
## Install countsimQC
BiocManager::install("csoneson/countsimQC")
```

## Getting started
To run `countsimQC` and generate a report, you simply need to call the
function `countsimQCReport()`, with an input consisting of a named list of
`DESeqDataSets` (see the
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
package for a description of this class). Each `DESeqDataSet` should
correspond to one data set and contain a count matrix, a data frame with sample
information and a design formula, which is needed for proper dispersion 
calculations. To generate a `DESeqDataSet` from a count matrix `counts`, a 
sample information data frame `sample_df` and a design formula `formula` 
(of the form `~ predictors`), you can do as follows:

```
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_df,
                              design = formula)
```
There are many other ways of generating valid `DESeqDataSets`, depending on in 
what form your counts are (e.g., reading directly from 
[HTSeq](http://htseq.readthedocs.io/en/release_0.9.1/) output, or from a [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html) 
output object (see the 
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) [vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)). 
 
`countsimQC` contains an small example list with subsets of three data sets: 
two synthetic ones and the real data set that was used to generate them. The 
following code generates a comparative report for these three data sets:

```
library(countsimQC)
data(countsimExample)
countsimQCReport(ddsList = countsimExample, 
                 outputFile = "countsimReport.html", 
                 outputDir = "./", 
                 description = "This is a comparison of three count data sets.")
```

For more detailed information about how to use the package, we refer to the vignette:

```
browseVignettes("countsimQC")
```

## Example reports

- [Comparison of 16S microbiome species count matrices for four body subsites from the Human Microbiome Project](http://imlspenticton.uzh.ch/robinson_lab/countsimQC_example_reports/HMP_sampled_datasets_countsimQC.html)
- [Comparison of three real bulk RNA-seq data sets](http://imlspenticton.uzh.ch/robinson_lab/countsimQC_example_reports/bulkrnaseq_crossdataset_countsimQC.html)
- [Comparison of gene- and transcript-level count matrices for a single-cell RNA-seq data set](http://imlspenticton.uzh.ch/robinson_lab/countsimQC_example_reports/GSE74596_genevstx_countsimQC.html) 
- [Comparison of four real scRNA-seq data sets](http://imlspenticton.uzh.ch/robinson_lab/countsimQC_example_reports/scrnaseq_crossdataset_countsimQC.html)
- [Comparison of two simulated scRNA-seq data sets to the underlying real data set](http://imlspenticton.uzh.ch/robinson_lab/countsimQC_example_reports/GSE48968-GPL13112_2simulations_countsimQC.html)
- [Comparison of six simulated bulk RNA-seq data set with different number of genes](http://imlspenticton.uzh.ch/robinson_lab/countsimQC_example_reports/countsimQC_compcodeR_simulations.html)

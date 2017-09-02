# countsimQC
[![Travis-CI Build Status](https://travis-ci.org/csoneson/countsimQC.svg?branch=master)](https://travis-ci.org/csoneson/countsimQC)

`countsimQC` is an R package that provides functionality to create a 
comprehensive report comparing many different characteristics of multiple 
count data sets. One of the main use cases is comparing one or more 
synthetic (e.g., RNA-seq) count matrices to a real count matrix, possibly the 
one based on which the synthetic data sets were generated. However, any collection of one or more count matrices can be visualized and compared.

## Installation
`countsimQC` depends on a number of other R packages. The following commands check whether the dependencies are available and installs them otherwise:

```
pkg <- c("rmarkdown", "edgeR", "DESeq2", "dplyr", "tidyr", "ggplot2", "grDevices", "tools", "SummarizedExperiment", "genefilter", "DT", "GenomeInfoDbData", "caTools", "randtests")
pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(pkg) > 0) {
	source("https://bioconductor.org/biocLite.R")
	biocLite(pkg, dependencies = TRUE, ask = FALSE)
}
```

Once all dependencies are available, `countsimQC` can be installed using the `devtools` package:

```
## Install devtools if needed
install.packages(devtools)

## Install countsimQC
devtools::install_github("csoneson/countsimQC")
```

## Getting started
To run `countsimQC` and generate the report, you simply need to call the
function `countsimQCReport()`, with an input consisting of a named list of
`DESeqDataSets` (see the
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
package for the description of this class). Each `DESeqDataSet` should
correspond to one data set and contain a count matrix, a data frame with sample
information and a design formula, which is needed for dispersion calculations. To generate a `DESeqDataSet` from a count matrix `counts`, a sample information data frame `sample_df` and a design formula `formula` (of the form `~ predictors`), you can do as follows:

```
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_df,
                              design = formula)
```
There are many other ways of generating valid `DESeqDataSets`, depending on in what form your counts are (e.g., reading directly from HTSeq output, or from a [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html) output object (see the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) [vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)). 
 
`countsimQC` contains an small example list with subsets of three data sets: two synthetic ones and the real data set that was used to generate them. To generate a comparative report from this list, simply run the following code:

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

# countsimQC
[![Travis-CI Build Status](https://travis-ci.org/csoneson/countsimQC.svg?branch=master)](https://travis-ci.org/csoneson/countsimQC)

`countsimQC` is an R package that provides functionality to create a 
comprehensive report comparing many different characteristics of multiple 
RNA-seq count matrices. One of the main use cases is comparing one or more 
synthetic (e.g., RNA-seq) count matrices to a real count matrix, possibly the 
one underlying the simulation. However, any collection of one or more count
matrices can be visualized and compared.

`countsimQC` can be installed using the `devtools` package:

```
## Install devtools if needed
install.packages(devtools)

## Install countsimQC
devtools::install_github("csoneson/countsimQC")
```

To run `countsimQC` and generate the report, you simply need to call the
function `countsimQCReport()`, with an input consisting of a named list of
`DESeqDataSets` (see the
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
package for the description of this class). Each `DESeqDataSet` should
correspond to one data set and contain a count matrix, a data frame with sample
information and a design formula, which is needed for dispersion calculations.
The package contains an example data set with subsets of three data sets:
two synthetic ones and the real data set that was used to generate them.

```
data(countsimExample)
countsimQCReport(ddsList = countsimExample, outputFile = "countsimReport.html", outputDir = "./", description = "This is a comparison of three count data sets.")
```

For more detailed information about how to use the package, see the vignette:

```
browseVignettes("countsimQC")
```

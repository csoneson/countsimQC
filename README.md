# countsimQC

`countsimQC` is an R package that provides functionality to create a comprehensive report 
comparing many different characteristics of multiple RNA-seq count matrices. 
It is mainly intended as a simple way to compare one or more synthetic RNA-seq count 
matrices to a real count matrix, possibly the one underlying the simulation. 
However, any collection of count matrices can be compared 

`countsimQC` can be installed using the `devtools` package:

```
## Install devtools if needed
install.packages(devtools)

## Install countsimQC
devtools::install_github("csoneson/countsimQC")
```

To run `countsimQC`, you simply need to call the function `countsimQC_report()`, with an input consisting of a named list of DESeqDataSets (see the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package for the description of this class). Each DESeqDataSet should correspond to one data set and contain a count matrix, a data frame with sample information and a design formula, which is needed for dispersion calculations. The package contains an example data set, containing three data sets. 

```
data(countsim_example)
countsimQC_report(countsim_example, output_file = "countsim_report.html", output_dir = "./")
```
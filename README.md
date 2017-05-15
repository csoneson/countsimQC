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

#' Calculate Spearman correlation between sample pairs
#'
#' This function is typically not called directly by the user, but is used when
#' generating the countsimQC report.
#'
#' @param ddsList List of lists, with one element per data set. Each element is
#'   a list containing a DGEList and a DESeqDataSet, with calculated
#'   dispersions (e.g., output from \code{\link{calculateDispersionsddsList}}).
#' @param maxNForCorr Maximal number of samples to use for correlation
#'   calculation. If the number of samples in a data set exceeds
#'   \code{maxNForCorr}, \code{maxNForCorr} samples will be randomly selected
#'   for calculation of correlations.
#' @param seed The random seed used for downsampling/selection of samples.
#'
#' @return A data frame with pairwise sample correlations for each data set
#' @author Charlotte Soneson
#' @importFrom stats cor
#' @importFrom edgeR cpm
#'
#' @examples
#' data(countsimExample)
#' ## Subset genes to speed up calculations
#' sub <- lapply(countsimExample, function(x) x[1:100, ])
#' ## Calculate dispersions
#' obj <- calculateDispersionsddsList(ddsList = sub,
#'                                    maxNForDisp = 5, seed = 123)
#' ## Calculate sample correlations
#' corrs <- calculateSampleCorrs(obj, maxNForCorr = 5, seed = 123)
#'
calculateSampleCorrs <- function(ddsList, maxNForCorr, seed = 123) {
  sampleCorrDF <- lapply(ddsList, function(x) {
    cpms <- edgeR::cpm(x$dge, prior.count = 2, log = TRUE)
    if (ncol(cpms) > maxNForCorr) {
      set.seed(seed)
      cpms <- cpms[, sample(seq_len(ncol(cpms)), maxNForCorr, replace = FALSE)]
    }
    corrs <- stats::cor(cpms, use = "pairwise.complete.obs", method = "spearman")
    data.frame(
      Correlation = corrs[upper.tri(corrs)]
    )
  })
  ns <- sapply(sampleCorrDF, nrow)
  do.call(rbind, sampleCorrDF) %>%
    dplyr::mutate(dataset = rep(names(sampleCorrDF), ns))
}

#' Calculate Spearman correlation between feature pairs
#'
#' This function is typically not called directly by the user, but is used when
#' generating the countsimQC report.
#'
#' @param ddsList List of lists, with one element per data set. Each element is
#'   a list containing a DGEList and a DESeqDataSet, with calculated
#'   dispersions (e.g., output from \code{\link{calculateDispersionsddsList}}).
#' @param maxNForCorr Maximal number of features to use for calculation of
#'   correlations. If the number of features in a data set exceeds
#'   \code{maxNForCorr}, \code{maxNForCorr} features will be randomly selected
#'   for calculation of correlations.
#' @param seed The random seed used for downsampling/selection of features.
#'
#' @return A data frame with pairwise feature correlations for each data set
#' @author Charlotte Soneson
#' @importFrom stats cor
#' @importFrom edgeR cpm
#' @importFrom genefilter rowVars
#'
#' @examples
#' data(countsimExample)
#' ## Subset genes to speed up calculations
#' sub <- lapply(countsimExample, function(x) x[1:100, ])
#' ## Calculate dispersions
#' obj <- calculateDispersionsddsList(ddsList = sub,
#'                                    maxNForDisp = 5, seed = 123)
#' ## Calculate feature correlations
#' corrs <- calculateFeatureCorrs(obj, maxNForCorr = 5, seed = 123)
#'
calculateFeatureCorrs <- function(ddsList, maxNForCorr, seed = 123) {
  featureCorrDF <- lapply(ddsList, function(x) {
    cpms <- edgeR::cpm(x$dge, prior.count = 2, log = TRUE)
    cpms <- cpms[genefilter::rowVars(cpms) > 0, ]
    if (nrow(cpms) > maxNForCorr) {
      set.seed(seed)
      cpms <- cpms[sample(seq_len(nrow(cpms)), maxNForCorr, replace = FALSE), ]
    }
    corrs <- stats::cor(t(cpms), use = "pairwise.complete.obs", method = "spearman")
    data.frame(
      Correlation = corrs[upper.tri(corrs)]
    )
  })
  ns <- sapply(featureCorrDF, nrow)
  do.call(rbind, featureCorrDF) %>%
    dplyr::mutate(dataset = rep(names(featureCorrDF), ns))

}

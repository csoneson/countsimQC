#' Calculate Spearman correlation between sample pairs
#'
#' @param ddsList List of lists, with one element per data set. Each element is
#'   a list containing a DGEList and a DESeqDataSet, with calculated
#'   dispersions.
#' @param maxNForCorr If the number of samples in a data set exceeds
#'   \code{maxNForCorr}, \code{maxNForCorr} samples will be randomly selected
#'   for correlation calculation.
#' @param seed The seed for downsampling
#'
#' @return A data frame with pairwise sample correlations for each data set
#' @author Charlotte Soneson
#' @importFrom stats cor
#' @importFrom edgeR cpm
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

#' Calculate Spearman correlation between variable pairs
#'
#' @param ddsList List of lists, with one element per data set. Each element is
#'   a list containing a DGEList and a DESeqDataSet, with calculated
#'   dispersions.
#' @param maxNForCorr If the number of variables in a data set exceeds
#'   \code{maxNForCorr}, \code{maxNForCorr} variables will be randomly selected
#'   for correlation calculation.
#' @param seed The seed for downsampling
#'
#' @return A data frame with pairwise variable correlations for each data set
#' @author Charlotte Soneson
#' @importFrom stats cor
#' @importFrom edgeR cpm
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

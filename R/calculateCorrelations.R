#' Calculate Spearman correlation between sample pairs
#'
#' @param ddsList List of lists, with one element per data set. Each element is
#'   a list containing a DGEList and a DESeqDataSet, with calculated
#'   dispersions (e.g., output from \code{\link{calculateDispersionsddsList}}).
#' @param maxNForCorr Maximal number of samples to use for correlation
#'   calculation. If the number of samples in a data set exceeds
#'   \code{maxNForCorr}, \code{maxNForCorr} samples will be randomly selected
#'   for calculation of correlations.
#'
#' @return A data frame with pairwise sample correlations for each data set
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
#' @importFrom stats cor
#' @importFrom edgeR cpm
#'
calculateSampleCorrs <- function(ddsList, maxNForCorr) {
  sampleCorrDF <- lapply(ddsList, function(x) {
    ## Calculate logCPMs
    cpms <- edgeR::cpm(x$dge, prior.count = 2, log = TRUE)

    ## Subsample columns if required
    if (ncol(cpms) > maxNForCorr) {
      cpms <- cpms[, sample(seq_len(ncol(cpms)), maxNForCorr, replace = FALSE)]
    }

    ## Calculate Spearman correlations
    corrs <- stats::cor(cpms, use = "pairwise.complete.obs",
                        method = "spearman")
    data.frame(
      Correlation = corrs[upper.tri(corrs)]
    )
  })
  ## Merge correlations from all data sets
  ns <- vapply(sampleCorrDF, nrow, 0)
  do.call(rbind, sampleCorrDF) %>%
    dplyr::mutate(dataset = rep(names(sampleCorrDF), ns))
}

#' Calculate Spearman correlation between feature pairs
#'
#' @param ddsList List of lists, with one element per data set. Each element is
#'   a list containing a DGEList and a DESeqDataSet, with calculated
#'   dispersions (e.g., output from \code{\link{calculateDispersionsddsList}}).
#' @param maxNForCorr Maximal number of features to use for calculation of
#'   correlations. If the number of features in a data set exceeds
#'   \code{maxNForCorr}, \code{maxNForCorr} features will be randomly selected
#'   for calculation of correlations.
#'
#' @return A data frame with pairwise feature correlations for each data set
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
#' @importFrom stats cor
#' @importFrom edgeR cpm
#' @importFrom genefilter rowVars
#'
calculateFeatureCorrs <- function(ddsList, maxNForCorr) {
  featureCorrDF <- lapply(ddsList, function(x) {
    ## Calculate logCPMs, keep only non-constant features
    cpms <- edgeR::cpm(x$dge, prior.count = 2, log = TRUE)
    cpms <- cpms[genefilter::rowVars(cpms) > 0, ]

    ## Subsample rows if required
    if (nrow(cpms) > maxNForCorr) {
      cpms <- cpms[sample(seq_len(nrow(cpms)), maxNForCorr, replace = FALSE), ]
    }

    ## Calculate Spearman correlations
    corrs <- stats::cor(t(cpms), use = "pairwise.complete.obs",
                        method = "spearman")
    data.frame(
      Correlation = corrs[upper.tri(corrs)]
    )
  })
  ## Merge correlations from all data sets
  ns <- vapply(featureCorrDF, nrow, 0)
  do.call(rbind, featureCorrDF) %>%
    dplyr::mutate(dataset = rep(names(featureCorrDF), ns))

}

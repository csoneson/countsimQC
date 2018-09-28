#' Calculate dispersions
#'
#' Calculate the dispersions for each data set in a list of DESeqDataSets, using
#' both edgeR and DESeq2.
#'
#' @param ddsList A list of DESeqDataSets
#' @param maxNForDisp If any data set contains more than \code{maxNForDisp}
#'   samples, \code{maxNForDisp} of them will be randomly sampled before the
#'   dispersions are calculated, in order to speed up calculations
#' @param seed The seed set for sampling
#'
#' @return A list of the same length as the input list. Each element in the list
#'   is itself a list, containing a DGEList and a DESeqDataSet with calculated
#'   dispersions.
#'
#' @author Charlotte Soneson
#' @import edgeR DESeq2 SummarizedExperiment
#' @importFrom stats model.matrix
#'
calculateDispersionsddsList <- function(ddsList, maxNForDisp, seed = 123) {
  lapply(ddsList, function(ds) {
    ## --------------------------- edgeR --------------------------- ##
    ## Define DGEList
    dge <- edgeR::DGEList(counts = DESeq2::counts(ds))
    ## Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge)
    ## Subset DGEList if the number of samples exceeds maxNForDisp and estimate
    ## dispersion
    if (ncol(dge) > maxNForDisp) {
      set.seed(seed)
      keepSamples <- sample(seq_len(ncol(dge)), maxNForDisp, replace = FALSE)
      dgetmp <- dge[, keepSamples]
      destmp <- stats::model.matrix(
        DESeq2::design(ds),
        data = droplevels(SummarizedExperiment::colData(ds)[keepSamples, ,
                                                            drop = FALSE])
      )
      dgetmp <- edgeR::estimateDisp(dgetmp, design = destmp)
      stopifnot(all(rownames(dge) == rownames(dgetmp)))
      dge$tagwise.dispersion <- dgetmp$tagwise.dispersion
      dge$common.dispersion <- dgetmp$common.dispersion
      dge$trended.dispersion <- dgetmp$trended.dispersion
      dge$prior.df <- dgetmp$prior.df
      dge$AveLogCPM <- edgeR::aveLogCPM(dge)
      dge$AveLogCPMDisp <- dgetmp$AveLogCPM
    } else {
      des <- stats::model.matrix(DESeq2::design(ds),
                                 data = SummarizedExperiment::colData(ds))
      dge <- edgeR::estimateDisp(dge, design = des)
      dge$AveLogCPMDisp <- dge$AveLogCPM
    }
    ## --------------------------- DESeq2 -------------------------- ##
    ## Calculate size factors
    dds <- DESeq2::estimateSizeFactors(ds, type = "poscounts")
    ## Subset DESeqDataSet if the number of samples exceeds maxNForDisp and
    ## estimate dispersion
    if (ncol(dds) > maxNForDisp) {
      set.seed(seed)
      keepSamples <- sample(seq_len(ncol(dds)), maxNForDisp, replace = FALSE)
      ddstmp <- dds[, keepSamples]
      colData(ddstmp) <- droplevels(colData(ddstmp))
      ddstmp <- DESeq2::estimateDispersions(ddstmp, quiet = TRUE)
      stopifnot(all(rownames(dds) == rownames(ddstmp)))
      rowData(dds)$dispGeneEst <- rowData(ddstmp)$dispGeneEst
      rowData(dds)$dispFit <- rowData(ddstmp)$dispFit
      rowData(dds)$dispersion <- rowData(ddstmp)$dispersion
      rowData(dds)$baseMean <- rowMeans(counts(dds, normalized = TRUE))
      rowData(dds)$baseMeanDisp <- rowData(ddstmp)$baseMean
    } else {
      dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
      rowData(dds)$baseMeanDisp <- rowData(dds)$baseMean
    }
    ## Check that the genes are in the same order in the DGEList and
    ## DESeqDataSet
    stopifnot(all(rownames(dge) == rownames(dds)))
    list(dge = dge, dds = dds)
  })
}

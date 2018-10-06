#' Calculate dispersions
#'
#' Calculate the dispersions for each data set in a list of DESeqDataSets, using
#' both edgeR and DESeq2.
#'
#' @param ddsList A list of DESeqDataSets
#' @param maxNForDisp If any data set contains more than \code{maxNForDisp}
#'   samples, \code{maxNForDisp} of them will be randomly sampled before the
#'   dispersions are calculated, in order to speed up calculations
#'
#' @return A list of the same length as the input list. Each element in the list
#'   is itself a list, containing a DGEList and a DESeqDataSet with calculated
#'   dispersions.
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
#' @import edgeR DESeq2 SummarizedExperiment
#' @importFrom stats model.matrix
#'
calculateDispersionsddsList <- function(ddsList, maxNForDisp) {
  lapply(ddsList, function(ds) {
    ## If the number of samples exceeds maxNForDisp, select random subset of
    ## samples to use for dispersion calculations (do it here so that the same
    ## set of samples is used for edgeR and DESeq2)
    if (ncol(ds) > maxNForDisp) {
      keepSamples <- sample(seq_len(ncol(ds)), maxNForDisp, replace = FALSE)
    }

    ## --------------------------- edgeR --------------------------- ##
    ## Define DGEList
    dge <- edgeR::DGEList(counts = DESeq2::counts(ds))
    ## Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge)
    if (ncol(dge) > maxNForDisp) {
      ## Subset DGEList if the number of samples exceeds maxNForDisp and
      ## estimate dispersions
      dgetmp <- dge[, keepSamples]
      destmp <- stats::model.matrix(
        DESeq2::design(ds),
        data = droplevels(SummarizedExperiment::colData(ds)[keepSamples, ,
                                                            drop = FALSE])
      )
      dgetmp <- edgeR::estimateDisp(dgetmp, design = destmp)

      ## Add dispersion estimates to original DGEList
      stopifnot(all(rownames(dge) == rownames(dgetmp)))
      dge$tagwise.dispersion <- dgetmp$tagwise.dispersion
      dge$common.dispersion <- dgetmp$common.dispersion
      dge$trended.dispersion <- dgetmp$trended.dispersion
      dge$prior.df <- dgetmp$prior.df
      dge$AveLogCPM <- edgeR::aveLogCPM(dge)
      dge$AveLogCPMDisp <- dgetmp$AveLogCPM
    } else {
      ## No subsampling, estimate dispersion from full DGEList
      des <- stats::model.matrix(DESeq2::design(ds),
                                 data = SummarizedExperiment::colData(ds))
      dge <- edgeR::estimateDisp(dge, design = des)
      dge$AveLogCPMDisp <- dge$AveLogCPM
    }
    ## --------------------------- DESeq2 -------------------------- ##
    ## Calculate size factors
    dds <- DESeq2::estimateSizeFactors(ds, type = "poscounts")
    if (ncol(dds) > maxNForDisp) {
      ## Subset DESeqDataSet if the number of samples exceeds maxNForDisp and
      ## estimate dispersions
      ddstmp <- dds[, keepSamples]
      colData(ddstmp) <- droplevels(colData(ddstmp))
      ddstmp <- DESeq2::estimateDispersions(ddstmp, quiet = TRUE)

      ## Add dispersion estimates to original DESeqDataSet
      stopifnot(all(rownames(dds) == rownames(ddstmp)))
      rowData(dds)$dispGeneEst <- rowData(ddstmp)$dispGeneEst
      rowData(dds)$dispFit <- rowData(ddstmp)$dispFit
      rowData(dds)$dispersion <- rowData(ddstmp)$dispersion
      rowData(dds)$baseMean <- rowMeans(counts(dds, normalized = TRUE))
      rowData(dds)$baseMeanDisp <- rowData(ddstmp)$baseMean
    } else {
      ## No subsampling, estimate dispersion from full DESeqDataSet
      dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
      rowData(dds)$baseMeanDisp <- rowData(dds)$baseMean
    }
    ## Check that the genes are in the same order in the DGEList and
    ## DESeqDataSet
    stopifnot(all(rownames(dge) == rownames(dds)))
    list(dge = dge, dds = dds)
  })
}

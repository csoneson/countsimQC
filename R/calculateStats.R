#' Calculate statistics for pairwise comparison of data sets
#'
#' Calculate a range of statistics and p-values for comparison of two data sets.
#'
#' @param df The input data frame. Must contain at least a column named
#'   'dataset' and an additional column with values to use as the basis for the
#'   comparison.
#' @param ds1,ds2 The names of the two data sets to be compared.
#' @param column The name of the column(s) of \code{df} to be used as the basis
#'   for the comparison.
#' @param subsampleSize The number of observations for which certain
#'   time-consuming statistics will be calculated. The observations will be
#'   selected randomly among the rows of \code{df}.
#' @param permute Whether to permute the dataset column of \code{df} before
#'   calculating the statistics.
#' @param kmin,kfrac For statistics that require the extraction of k nearest
#'   neighbors of a given point, the number of neighbors will be max(kmin, kfrac
#'   * nrow(df)).
#' @param xmin,xmax Smallest and largest value of \code{column}, used to
#'   normalize the x-axis when calculating the area between the eCDFs.
#'
#' @return A vector with statistics and p-values
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
#' @importFrom stats ks.test ecdf chisq.test
#'
calculateStats <- function(df, ds1, ds2, column, subsampleSize,
                           permute = FALSE, kmin, kfrac, xmin, xmax) {
  ## Remove rows with NA values in column(s) of interest
  df <- df[rowSums(is.na(df[, column, drop = FALSE])) == 0, ]

  ## Permute dataset column if requested
  if (permute) {
    df$dataset <- sample(df$dataset, nrow(df))
  }

  if (length(column) == 1) {
    ## Kolmogorov-Smirnov test
    kst <- stats::ks.test(df[, column][df$dataset == ds1],
                          df[, column][df$dataset == ds2])

    ## Area between ecdfs
    e1 <- stats::ecdf(df[, column][df$dataset == ds1])
    e2 <- stats::ecdf(df[, column][df$dataset == ds2])
    xv <- sort(unique(df[, column]))
    ediff <- abs(e1(xv) - e2(xv))
    earea <- caTools::trapz((xv - xmin)/(xmax - xmin), ediff)
  }

  ## Silhouette width and nearest neighbor distribution for subsample of
  ## observations
  ## Calculate overall proportion of values from the two data sets
  overallprop <- table(factor(df$dataset, levels = c(ds1, ds2)))
  overallprop <- overallprop/sum(overallprop)
  nrep <- 1
  chisilh <- t(vapply(seq_len(nrep), function(m) {
    ## Select subset of observations for calculation of silhouette widths and
    ## nearest neighbor distributions
    idx <- sample(seq_len(nrow(df)), min(nrow(df), subsampleSize))

    ## For each selected observation, calculate statistics
    tmp <- t(vapply(idx, function(j) {

      ## Distances to all observations
      dists <- sqrt(rowSums((sweep(df[, column, drop = FALSE],
                                   2, unlist(df[j, column,
                                                drop = FALSE]), "-")) ^ 2))
      dists_this <- dists[df$dataset == df$dataset[j]]
      dists_other <- dists[df$dataset != df$dataset[j]]

      ## Silhouette width
      silh <- (mean(dists_other, na.rm = TRUE) -
                 mean(dists_this, na.rm = TRUE)) /
        max(mean(dists_other, na.rm = TRUE), mean(dists_this, na.rm = TRUE))

      ## Local silhouette width
      dists_this_local <-
        sort(dists_this)[seq_len(max(kmin, kfrac * nrow(df)) *
                                   overallprop[df$dataset[j]])]
      dists_other_local <-
        sort(dists_other)[seq_len(max(kmin, kfrac * nrow(df)) *
                                    (1 - overallprop[df$dataset[j]]))]
      silh_local <-
        (mean(dists_other_local, na.rm = TRUE) -
           mean(dists_this_local, na.rm = TRUE)) /
        max(mean(dists_other_local, na.rm = TRUE),
            mean(dists_this_local, na.rm = TRUE))
      if (all(c(dists_this_local, dists_other_local) == 0)) {
        silh_local <- 0
      }

      ## Chi-square test comparing distribution of data set labels among
      ## nearest neighbors to the overall distribution
      x <- df$dataset[order(dists)][seq_len(max(kmin, kfrac * nrow(df)))]
      suppressWarnings({
        chisqp <- stats::chisq.test(x = table(factor(x, levels = c(ds1, ds2))),
                                    p = overallprop)$p.value
      })
      ## Return values
      c(silh = silh, chisqp = chisqp, silh_local = silh_local)
    }, c(silh = NA_real_, chisqp = NA_real_, silh_local = NA_real_)))
    c(silh = mean(tmp[, "silh"], na.rm = TRUE),
      silhlocal = mean(tmp[, "silh_local"], na.rm = TRUE),
      chisqp = mean(tmp[, "chisqp"] <= 0.05, na.rm = TRUE))
  }, c(silh = NA_real_, silhlocal = NA_real_, chisqp = NA_real_)))

  ## Runs test
  if (length(column) == 1) {
    df <- df %>% arrange_(column)
    runs_res <- randtests::runs.test(as.numeric(df$dataset == ds1),
                                     threshold = 0.5,
                                     alternative = "left.sided", plot = FALSE)
  }

  if (length(column) == 1)
    c(ksstatistic = signif(kst$statistic, 3),
      kspvalue = signif(kst$p.value, 3),
      diffarea = signif(earea, 3),
      runsstatistic = signif(runs_res$statistic, 3),
      runspvalue = signif(runs_res$p.value, 3),
      NNmismatch = signif(mean(chisilh[, "chisqp"]), 3),
      avesilh = signif(mean(chisilh[, "silh"]), 3),
      avesilhlocal = signif(mean(chisilh[, "silhlocal"]), 3))
  else
    c(NNmismatch = signif(mean(chisilh[, "chisqp"]), 3),
      avesilh = signif(mean(chisilh[, "silh"]), 3),
      avesilhlocal = signif(mean(chisilh[, "silhlocal"]), 3))
}

#' Return a vector of NA scores
#'
#' @param n Number of columns to use for the comparison
#' @param withP Whether or not to include p-value columns
#'
#' @return A vector with NA values for all applicable statistics
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
defaultStats <- function(n, withP = FALSE) {
  if (n == 1) {
    if (withP) {
      c(ksstatistic = NA_real_, kspvalue = NA_real_, diffarea = NA_real_,
        runsstatistic = NA_real_, runspvalue = NA_real_, NNmismatch = NA_real_,
        avesilh = NA_real_, avesilhlocal = NA_real_, NNmismatchP = NA_real_,
        avesilhP = NA_real_, avesilhlocalP = NA_real_, diffareaP = NA_real_)
    } else {
      c(ksstatistic = NA_real_, kspvalue = NA_real_, diffarea = NA_real_,
        runsstatistic = NA_real_, runspvalue = NA_real_, NNmismatch = NA_real_,
        avesilh = NA_real_, avesilhlocal = NA_real_)
    }
  } else {
    if (withP) {
      c(NNmismatch = NA_real_, avesilh = NA_real_, avesilhlocal = NA_real_,
        NNmismatchP = NA_real_, avesilhP = NA_real_, avesilhlocalP = NA_real_)
    } else {
      c(NNmismatch = NA_real_, avesilh = NA_real_, avesilhlocal = NA_real_)
    }
  }
}

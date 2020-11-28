library(countsimQC)
context("Check individual functions")

local({
  data(countsimExample)
  x <- lapply(countsimExample, function(w) w[1:50, ])
  set.seed(1)
  xsub <- calculateDispersionsddsList(x, maxNForDisp = 4)
  x <- calculateDispersionsddsList(x, maxNForDisp = Inf)

  test_that("Dimensions are correct", {
    expect_equivalent(unlist(lapply(x, function(w) c(nrow(w$dge), nrow(w$dds)))),
                      rep(50, 6))
    expect_equivalent(unlist(lapply(x, function(w) c(ncol(w$dge), ncol(w$dds)))),
                      rep(11, 6))
    expect_equivalent(unlist(lapply(xsub, function(w) c(nrow(w$dge), nrow(w$dds)))),
                      rep(50, 6))
    expect_equivalent(unlist(lapply(xsub, function(w) c(ncol(w$dge), ncol(w$dds)))),
                      rep(11, 6))
  })

  ## Sample correlations
  ## With calculateSampleCorrs
  set.seed(1)
  sample_corrs <- calculateSampleCorrs(x, maxNForCorr = 3)
  ## Manually
  set.seed(1)
  cpms <- edgeR::cpm(x$Original$dge, prior.count = 2, log = TRUE)
  samples_to_keep <- sample(seq_len(ncol(cpms)), 3, replace = FALSE)
  test_that("Sample correlations are correct", {
    expect_equal(subset(sample_corrs, dataset == "Original")$Correlation[1],
                 stats::cor(cpms[, samples_to_keep[1]], cpms[, samples_to_keep[2]],
                            method = "spearman", use = "pairwise.complete.obs"))
    expect_equal(subset(sample_corrs, dataset == "Original")$Correlation[2],
                 stats::cor(cpms[, samples_to_keep[1]], cpms[, samples_to_keep[3]],
                            method = "spearman", use = "pairwise.complete.obs"))
  })

  ## Feature correlations
  ## With calculateFeatureCorrs
  set.seed(1)
  feature_corrs <- calculateFeatureCorrs(x, maxNForCorr = 3)
  ## Manually
  set.seed(1)
  cpms <- edgeR::cpm(x$Original$dge, prior.count = 2, log = TRUE)
  cpms <- cpms[genefilter::rowVars(cpms) > 0, ]
  features_to_keep <- sample(seq_len(nrow(cpms)), 3, replace = FALSE)
  test_that("Feature correlations are correct", {
    expect_equal(subset(feature_corrs, dataset == "Original")$Correlation[1],
                 stats::cor(cpms[features_to_keep[1], ], cpms[features_to_keep[2], ],
                            method = "spearman", use = "pairwise.complete.obs"))
    expect_equal(subset(feature_corrs, dataset == "Original")$Correlation[2],
                 stats::cor(cpms[features_to_keep[1], ], cpms[features_to_keep[3], ],
                            method = "spearman", use = "pairwise.complete.obs"))
  })

  ## Stats
  sampleDF <- lapply(x, function(w) {
    data.frame(
      Libsize = colSums(w$dge$counts),
      Fraczero = colMeans(w$dge$counts == 0),
      TMM = w$dge$samples$norm.factors
    ) %>% dplyr::mutate(EffLibsize = Libsize * TMM)
  })
  ns <- sapply(sampleDF, nrow)
  sampleDF <- do.call(rbind, sampleDF) %>%
    dplyr::mutate(dataset = rep(names(sampleDF), ns))

  test_that("Summary stats are correct", {
    expect_equal(subset(sampleDF, dataset == "Sim1")$Fraczero[1],
                 mean(x$Sim1$dge$counts[, 1] == 0))
    expect_equal(subset(sampleDF, dataset == "Sim2")$Libsize[1],
                 sum(x$Sim2$dge$counts[, 1]))
    expect_equivalent(subset(sampleDF, dataset == "Sim2")$TMM[1],
                      edgeR::calcNormFactors(x$Sim2$dge$counts)[1])
    expect_equivalent(subset(sampleDF, dataset == "Sim2")$EffLibsize[1],
                      edgeR::calcNormFactors(x$Sim2$dge$counts)[1] *
                        sum(x$Sim2$dge$counts[, 1]))
  })

  statDF <- makeDF(df = sampleDF, column = "TMM",
                   permutationPvalues = FALSE, nPermutations = 0,
                   subsampleSize = 11, kmin = 3, kfrac = 0.1)$x$data
  set.seed(1)
  statDFperm <- makeDF(df = sampleDF, column = "TMM",
                       permutationPvalues = TRUE, nPermutations = 100,
                       subsampleSize = 11, kmin = 3, kfrac = 0.1)$x$data

  ## Average silhouette width is different since we're taking a subsample
  compare_cols <- setdiff(colnames(statDF), "Average silhouette width")
  test_that("Stats are the same whether permutation p-values are calculated or not", {
    expect_equal(statDF[, compare_cols], statDFperm[, compare_cols])
  })

  test_that("Derived stats are correct", {
    expect_equivalent(signif(stats::ks.test(sampleDF$TMM[sampleDF$dataset == "Original"],
                                            sampleDF$TMM[sampleDF$dataset == "Sim1"])$stat, 3),
                      statDF$`K-S statistic`[statDF$dataset1 == "Original" & statDF$dataset2 == "Sim1"])
    expect_equivalent(signif(stats::ks.test(sampleDF$TMM[sampleDF$dataset == "Original"],
                                            sampleDF$TMM[sampleDF$dataset == "Sim1"])$p.value, 3),
                      statDF$`K-S p-value`[statDF$dataset1 == "Original" & statDF$dataset2 == "Sim1"])
  })

  expect_equal(length(countsimQC:::defaultStats(n = 1, withP = FALSE)), 8L)
  expect_equal(length(countsimQC:::defaultStats(n = 1, withP = TRUE)), 12L)
  expect_equal(length(countsimQC:::defaultStats(n = 2, withP = FALSE)), 3L)
  expect_equal(length(countsimQC:::defaultStats(n = 2, withP = TRUE)), 6L)

  dtd <- countsimQC:::defineTableDesc(
    calculateStatistics = TRUE, subsampleSize = 5, kfrac = 0.25, kmin = 4,
    obstype = "sample", aspect = "aspect", minvalue = 0, maxvalue = 1,
    permutationPvalues = TRUE, nPermutations = 10, nDatasets = 3)
  expect_is(dtd, "list")
  expect_length(dtd, 2)
  expect_named(dtd, c("tabledesc", "tabledesc2d"), ignore.order = TRUE)

  dtdnoperm <- countsimQC:::defineTableDesc(
    calculateStatistics = TRUE, subsampleSize = 5, kfrac = 0.25, kmin = 4,
    obstype = "sample", aspect = "aspect", minvalue = 0, maxvalue = 1,
    permutationPvalues = FALSE, nPermutations = 10, nDatasets = 3)
  expect_is(dtdnoperm, "list")
  expect_length(dtdnoperm, 2)
  expect_named(dtdnoperm, c("tabledesc", "tabledesc2d"), ignore.order = TRUE)

  dtd <- countsimQC:::defineTableDesc(
    calculateStatistics = TRUE, subsampleSize = 5, kfrac = 0.25, kmin = 4,
    obstype = "sample", aspect = "aspect", minvalue = 0, maxvalue = 1,
    permutationPvalues = TRUE, nPermutations = 10, nDatasets = 1)
  expect_is(dtd, "list")
  expect_length(dtd, 2)
  expect_named(dtd, c("tabledesc", "tabledesc2d"), ignore.order = TRUE)
  expect_equal(dtd[[1]],
               "No statistics were calculated, since there is only one data set.")
})

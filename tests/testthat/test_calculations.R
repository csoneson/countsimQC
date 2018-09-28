library(countsimQC)
context("Check individual functions")

local({
  x <- lapply(countsimExample, function(w) w[1:50, ])
  x <- calculateDispersionsddsList(x, maxNForDisp = Inf)

  test_that("Dimensions are correct", {
    expect_equivalent(unlist(lapply(x, function(w) c(nrow(w$dge), nrow(w$dds)))),
                      rep(50, 6))
    expect_equivalent(unlist(lapply(x, function(w) c(ncol(w$dge), ncol(w$dds)))),
                      rep(11, 6))
  })

  ## Correlations
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

})

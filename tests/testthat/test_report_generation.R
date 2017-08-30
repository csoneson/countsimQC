library(countsimQC)
context("Check report generation")

local({
  x <- lapply(countsimExample, function(w) w[1:50, ])
  x2 <- lapply(countsimExample_dfmat, function(w) w[1:50, ])
  if (file.exists("report_from_test.Rmd")) {
    file.remove("report_from_test.Rmd")
  }
  if (file.exists("report_from_test.html")) {
    file.remove("report_from_test.html")
  }
  if (file.exists("report_from_test.pdf")) {
    file.remove("report_from_test.pdf")
  }
  if (file.exists("report_from_test_ggplots.rds")) {
    file.remove("report_from_test_ggplots.rds")
  }

  test_that("Input argument errors are recognized", {
    ## outputFormat
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html",
                                  outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE))
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "report_from_test.pdf",
                                  outputDir = "./", savePlots = TRUE))
    ## ddsList
    expect_error(countsimQCReport(ddsList = x[[1]], outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE))
    expect_error(countsimQCReport(ddsList = list(A = 1:3),
                                  outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE))
    y <- x
    names(y) <- NULL
    expect_error(countsimQCReport(ddsList = y, outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE))
    y <- x
    names(y) <- c("A", "A", "B")
    expect_error(countsimQCReport(ddsList = y, outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE))

    ## matrix/data frame input
    expect_is(countsimQCReport(ddsList = x2, outputFormat = "html_document",
                               outputFile = "report_from_test.html",
                               outputDir = "./", savePlots = FALSE,
                               calculateStatistics = FALSE,
                               forceOverwrite = TRUE), "character")

    ## rmd_template
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE,
                                  rmdTemplate = "nonexistent.Rmd"))

    ## existing rmd output
    file.create("test_generated.Rmd")
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "test_generated.html",
                                  outputDir = "./", savePlots = TRUE))
    file.remove("test_generated.Rmd")
  })

  countsimQCReport(ddsList = x, outputFile = "report_from_test.html",
                   outputDir = "./", savePlots = TRUE, forceOverwrite = TRUE,
                   permutationPvalues = TRUE, nPermutations = 3)

  test_that("Report is not generated if forceOverwrite = FALSE", {
    ## Report should not be generated if it already exists and forceOverwrite = FALSE
    expect_error(countsimQCReport(ddsList = x, outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE,
                                  forceOverwrite = FALSE))

    ## Report should be generated if it already exists and forceOverwrite = TRUE
    expect_is(countsimQCReport(ddsList = x, outputFile = "report_from_test.html",
                               outputDir = "./", savePlots = TRUE,
                               calculateStatistics = FALSE,
                               forceOverwrite = TRUE), "character")

    ## Report should not be generated if the ggplots object exists and forceOverwrite = FALSE
    file.remove("./report_from_test.html")
    expect_error(countsimQCReport(ddsList = x, outputFile = "report_from_test.html",
                                  outputDir = "./", savePlots = TRUE,
                                  forceOverwrite = FALSE))

    ## Report should be generated even if the ggplots object exists, if savePlots = FALSE
    expect_is(countsimQCReport(ddsList = x, outputFile = "report_from_test.html",
                               outputDir = "./", savePlots = FALSE,
                               calculateStatistics = FALSE,
                               forceOverwrite = FALSE), "character")

    ## pdf report should be generated even if the html report exists
    expect_is(countsimQCReport(ddsList = x, outputFormat = "pdf_document",
                               outputFile = "report_from_test.pdf",
                               outputDir = "./", savePlots = FALSE,
                               calculateStatistics = FALSE,
                               forceOverwrite = FALSE), "character")
  })
})

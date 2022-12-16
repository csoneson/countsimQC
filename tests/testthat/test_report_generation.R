library(countsimQC)
context("Check report generation")

local({
  tempDir <- tempdir()
  data(countsimExample)
  x <- lapply(countsimExample, function(w) w[1:25, ])
  if (file.exists(file.path(tempDir, "report_from_test.Rmd"))) {
    file.remove(file.path(tempDir, "report_from_test.Rmd"))
  }
  if (file.exists(file.path(tempDir, "report_from_test.html"))) {
    file.remove(file.path(tempDir, "report_from_test.html"))
  }
  if (file.exists(file.path(tempDir, "report_from_test.pdf"))) {
    file.remove(file.path(tempDir, "report_from_test.pdf"))
  }
  if (file.exists(file.path(tempDir, "report_from_test_ggplots.rds"))) {
    file.remove(file.path(tempDir, "report_from_test_ggplots.rds"))
  }

  test_that("Input argument errors are recognized", {
    ## outputFormat
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html",
                                  outputFile = "report_from_test.html",
                                  outputDir = tempDir, savePlots = TRUE))
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "report_from_test.pdf",
                                  outputDir = tempDir, savePlots = TRUE))
    ## ddsList
    expect_error(countsimQCReport(ddsList = x[[1]], outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = tempDir, savePlots = TRUE))
    expect_error(countsimQCReport(ddsList = list(A = 1:3),
                                  outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = tempDir, savePlots = TRUE))
    y <- x
    names(y) <- NULL
    expect_error(countsimQCReport(ddsList = y, outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = tempDir, savePlots = TRUE))
    y <- x
    names(y) <- c("A", "A", "B")
    expect_error(countsimQCReport(ddsList = y, outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = tempDir, savePlots = TRUE))

    ## rmd_template
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "report_from_test.html",
                                  outputDir = tempDir, savePlots = TRUE,
                                  rmdTemplate = "nonexistent.Rmd"))

    ## existing rmd output
    file.create(file.path(tempDir, "test_generated.Rmd"))
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "test_generated.html",
                                  outputDir = tempDir, savePlots = TRUE))
    file.remove(file.path(tempDir, "test_generated.Rmd"))

    ## existing plot output
    file.create(file.path(tempDir, "test_generated2_ggplots.rds"))
    expect_error(countsimQCReport(ddsList = x, outputFormat = "html_document",
                                  outputFile = "test_generated2.html",
                                  outputDir = tempDir, savePlots = TRUE,
                                  forceOverwrite = FALSE))
    file.remove(file.path(tempDir, "test_generated2_ggplots.rds"))
  })

  ## This doesn't generate a warning ("Chi-squared approximation may be incorrect")
  ## anymore (starting from 2022-12-15).
  ## Also the progress messages are not displayed in the console anymore.
  ## Maybe due to changes in rmarkdown 2.19? They are still displayed with 2.18.
  wns <- capture_warnings({
    cqr <- countsimQCReport(ddsList = x, outputFormat = "html_document",
                            outputFile = "test_generated.html",
                            outputDir = tempDir, savePlots = TRUE,
                            description = NULL)
    })
  expect_is(cqr, "character")
  expect_error(countsimQCReport(ddsList = x, outputFormat = NULL,
                                outputFile = "test_generated.html",
                                outputDir = tempDir, savePlots = TRUE,
                                forceOverwrite = FALSE))
  ## Test that existing report will be overwritten if requested
  # expect_warning(
  #   cqr <- countsimQCReport(ddsList = x, outputFormat = NULL,
  #                           outputFile = "test_generated.html",
  #                           outputDir = tempDir, savePlots = TRUE,
  #                           forceOverwrite = TRUE))
  # expect_is(cqr, "character")

  # expect_warning(generateIndividualPlots(
  #   ggplotsRds = file.path(tempDir, "test_generated_ggplots.rds"),
  #   device = "png",
  #   outputDir = tempDir, nDatasets = 3
  # ))
  expect_warning(generateIndividualPlots(
    ggplotsRds = file.path(tempDir, "test_generated_ggplots.rds"),
    device = "png",
    outputDir = file.path(tempDir, "subdir"), nDatasets = 3
  ))

  expect_error(generateIndividualPlots(
    ggplotsRds = 1, device = "png",
    outputDir = file.path(tempDir, "subdir"), nDatasets = 3
  ), "The provided ggplotsRds object")
  expect_error(generateIndividualPlots(
    ggplotsRds = list(1, 2), device = "png",
    outputDir = file.path(tempDir, "subdir"), nDatasets = 3
  ), "The elements of the provided ggplotsRds object")

  expect_warning(countsimQCReport(ddsList = x, outputFormat = NULL,
                                  outputFile = "test_generated.html",
                                  outputDir = tempDir, savePlots = FALSE,
                                  forceOverwrite = TRUE),
                 "already exists, but will not be overwritten")

  ## Run with matrix input
  expect_warning(
    cqr <- countsimQCReport(ddsList = lapply(x, function(y) DESeq2::counts(y)),
                            outputFormat = NULL,
                            outputFile = "test_generated.html",
                            outputDir = tempDir, savePlots = FALSE,
                            forceOverwrite = TRUE))
  expect_is(cqr, "character")

  if (file.exists(file.path(tempDir, "report_from_test.html"))) {
    file.remove(file.path(tempDir, "report_from_test.html"))
  }
  if (file.exists(file.path(tempDir, "report_from_test.pdf"))) {
    file.remove(file.path(tempDir, "report_from_test.pdf"))
  }
  if (file.exists(file.path(tempDir, "report_from_test_ggplots.rds"))) {
    file.remove(file.path(tempDir, "report_from_test_ggplots.rds"))
  }
})

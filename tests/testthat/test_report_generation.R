library(countsimQC)
context("Check report generation")

local({
  x <- lapply(countsimExample, function(w) w[1:50, ])
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

  if (file.exists("report_from_test.html")) {
    file.remove("report_from_test.html")
  }
  if (file.exists("report_from_test.pdf")) {
    file.remove("report_from_test.pdf")
  }
  if (file.exists("report_from_test_ggplots.rds")) {
    file.remove("report_from_test_ggplots.rds")
  }
})

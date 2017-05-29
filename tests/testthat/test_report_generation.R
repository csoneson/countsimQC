library(countsimQC)
context("Check report generation")

local({
  x <- lapply(countsim_example, function(w) w[1:100, ])
  if (file.exists("report_from_test.rmd")) {
    file.remove("report_from_test.rmd")
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
    ## output_format
    expect_error(countsimQC_report(dds_list = x, output_format = "html",
                                   output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE))
    expect_error(countsimQC_report(dds_list = x, output_format = "html_document",
                                   output_file = "report_from_test.pdf",
                                   output_dir = "./", save_ggplots = TRUE))
    ## dds_list
    expect_error(countsimQC_report(dds_list = x[[1]], output_format = "html_document",
                                   output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE))
    expect_error(countsimQC_report(dds_list = list(A = 1:3),
                                   output_format = "html_document",
                                   output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE))
    y <- x
    names(y) <- NULL
    expect_error(countsimQC_report(dds_list = y, output_format = "html_document",
                                   output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE))
    y <- x
    names(y) <- c("A", "A", "B")
    expect_error(countsimQC_report(dds_list = y, output_format = "html_document",
                                   output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE))
    ## rmd_template
    expect_error(countsimQC_report(dds_list = x, output_format = "html_document",
                                   output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE,
                                   rmd_template = "nonexistent.Rmd"))

    ## existing rmd output
    file.create("test_generated.Rmd")
    expect_error(countsimQC_report(dds_list = x, output_format = "html_document",
                                   output_file = "test_generated.html",
                                   output_dir = "./", save_ggplots = TRUE))
    file.remove("test_generated.Rmd")
  })

  countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                    output_dir = "./", save_ggplots = TRUE)

  test_that("Report is not generated if force_overwrite = FALSE", {
    ## Report should not be generated if it already exists and force_overwrite = FALSE
    expect_error(countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE,
                                   force_overwrite = FALSE))

    ## Report should be generated if it already exists and force_overwrite = TRUE
    expect_is(countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                                output_dir = "./", save_ggplots = TRUE,
                                force_overwrite = TRUE), "character")

    ## Report should not be generated if the ggplots object exists and force_overwrite = FALSE
    file.remove("./report_from_test.html")
    expect_error(countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE,
                                   force_overwrite = FALSE))

    ## Report should be generated even if the ggplots object exists, if save_ggplots = FALSE
    expect_is(countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                                output_dir = "./", save_ggplots = FALSE,
                                force_overwrite = FALSE), "character")

    ## pdf report should be generated even if the html report exists
    expect_is(countsimQC_report(dds_list = x, output_format = "pdf_document",
                                output_file = "report_from_test.pdf",
                                output_dir = "./", save_ggplots = FALSE,
                                force_overwrite = FALSE), "character")
  })
})

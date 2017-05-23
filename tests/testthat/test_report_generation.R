library(countsimQC)
context("Check report generation")

local({
  x <- lapply(countsim_example, function(w) w[1:500, ])
  if (file.exists("report_from_test.html")) {
    file.remove("report_from_test.html")
  }
  if (file.exists("report_from_test.pdf")) {
    file.remove("report_from_test.pdf")
  }
  if (file.exists("report_from_test_ggplots.rds")) {
    file.remove("report_from_test_ggplots.rds")
  }

  countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                    output_dir = "./", save_ggplots = TRUE)

  test_that("Report is not generated if force_overwrite = FALSE", {
    ## Report should not be generated if it already exists and force_overwrite = FALSE
    expect_error(countsimQC_report(dds_list = x, output_file = "report_from_test.html",
                                   output_dir = "./", save_ggplots = TRUE,
                                   force_overwrite = FALSE))

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

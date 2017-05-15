#' Generate report
#'
#' @param dds_list Named list of DESeqDataSets to compare
#' @param show_code Whether or not to include code in report
#' @param output_format The output format of the report. If set to NULL, an html
#'   report will be generated. Other options are "html_document" or
#'   "pdf_document", or "all" to generate both an html report and a pdf report.
#' @param output_file The name of the report
#' @param output_dir The directory where the report will be saved
#' @param ... Other arguments that will be passed to \code{rmarkdown::render}
#'
#' @author Charlotte Soneson, using a template from Nicholas Hamilton
#'   (http://stackoverflow.com/questions/37097535/generate-report-in-r)
#'
#' @export
#' @importFrom rmarkdown render
#' @return No value is returned, but a report is generated in the
#'   \code{output_dir} directory.
#'
generate_report <- function(dds_list, show_code = FALSE, output_format = NULL,
                            output_file = NULL, output_dir = "./", ...){
  ## This function was inspired by code from Nicholas Hamilton, obtained from
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r

  ## Path to the template file
  template_file <- system.file("extdata", "simulation_comparison_template_list.Rmd",
                               package = "countsimQC")

  ## Process the arguments
  args <- list(...)
  args$input <- template_file
  args$output_dir <- output_dir
  args$output_format <- output_format
  args$output_file <- output_file

  ## Render the report
  output_file <- do.call("render", args = args)
  invisible(output_file)
}

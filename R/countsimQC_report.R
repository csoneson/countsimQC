#' Generate countsimQC report
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
countsimQC_report <- function(dds_list, show_code = FALSE, output_format = NULL,
                              output_file = NULL, output_dir = "./", ...){
  ## This function was inspired by code from Nicholas Hamilton, provided at
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r

  ## Raise an error if the input list is not named
  if (length(setdiff(unique(names(dds_list)), c("", NA, NULL))) != length(dds_list))
    stop("The dds_list must be a named list, with unique names for all elements.")

  ## Raise an error if the output directory already exists
  if (file.exists(paste0(output_dir, "/", gsub(".html", "_figures", basename(output_file)))))
    stop("The folder ",
         paste0(output_dir, "/", gsub(".html", "_figures", basename(output_file))),
         ", where the figures will be generated, already exists. ",
         "Please remove or rename the folder, or provide another value of output_file.")

  ## Path to the template file
  template_file <- system.file("extdata", "simulation_comparison_template_list.Rmd",
                               package = "countsimQC")

  ## Process the arguments
  args <- list(...)
  args$input <- template_file
  args$output_format <- output_format
  args$output_file <- output_file

  ## Render the report
  output_file <- do.call("render", args = args)

  ## Move output files
  file.rename(output_file, paste0(output_dir, "/", basename(output_file)))
  file.rename(gsub(".html", "_figures", output_file),
              paste0(output_dir, "/", gsub(".html", "_figures", basename(output_file))))

  invisible(output_file)
}

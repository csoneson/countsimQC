#' Generate countsimQC report
#'
#' Generate a report comparing various characteristics of a collection of count
#' data sets.
#'
#' @param dds_list Named list of \code{DESeqDataSets} to compare. See the
#'   \code{DESeq2} Bioconductor package
#'   (http://bioconductor.org/packages/release/bioc/html/DESeq2.html) for more
#'   information about the \code{DESeqDataSet} class. Each \code{DESeqDataSet}
#'   object in the list should contain a count matrix, a data frame with sample
#'   information and a design formula. The sample information and design formula
#'   will be used to calculate dispersions appropriately.
#' @param show_code Whether or not to include the code in report.
#' @param output_format The output format of the report. If set to NULL, an html
#'   report will be generated. Other options are "html_document" or
#'   "pdf_document" to generate an html or a pdf report, respectively.
#' @param output_file The file name of the report.
#' @param output_dir The directory where the report should be saved.
#' @param rmd_template The Rmarkdown (.Rmd) file that will be used as the
#'   template for generating the report. If set to NULL (default), the template
#'   provided with the countsimQC package will be used. See Details for more
#'   information.
#' @param force_overwrite Whether to force overwrite existing output files and
#'   directories when saving the generated report and figures.
#' @param save_ggplots Whether to save the ggplot objects for all the output
#'   figures, to allow additional fine-tuning and generation of individual
#'   plots. Note, however, that the resulting file can be quite large,
#'   especially when large or many data sets are compared.
#' @param ... Other arguments that will be passed to \code{rmarkdown::render}.
#'
#' @author Charlotte Soneson, using a template from Nicholas Hamilton
#'   (http://stackoverflow.com/questions/37097535/generate-report-in-r)
#'
#' @details When the function is called, the template file (specified by
#'   rmd_template) will be copied into the output folder, and
#'   \code{rmarkdown::render} will be called to generate the final reports. If
#'   there is already a .Rmd file with the same name in the output folder, the
#'   function will raise an error and stop, to avoid overwriting the existing
#'   file. The reason for this behaviour is that the copied template in the
#'   output folder will be deleted once the report is generated.
#'
#' @export
#' @importFrom rmarkdown render
#' @importFrom tools file_ext
#' @import DESeq2 edgeR dplyr tidyr ggplot2
#' @return No value is returned, but a report is generated in the
#'   \code{output_dir} directory.
#'
countsimQC_report <- function(dds_list, show_code = FALSE, output_format = NULL,
                              output_file = NULL, output_dir = "./",
                              rmd_template = NULL, force_overwrite = FALSE,
                              save_ggplots = FALSE, ...){
  ## This function was inspired by code from Nicholas Hamilton, provided at
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r

  if (is.null(output_format)) output_format <- "html_document"

  ## ------------------------------------------------------------------------ ##
  ## -------------------------- Check input arguments ----------------------- ##
  ## ------------------------------------------------------------------------ ##

  ## Raise an error if output_format is not one of the allowed
  if (!(output_format %in% c("pdf_document", "html_document")))
    stop("The provided output_format is currently not supported. Please use ",
         "either 'html_document' (or NULL) or 'pdf_document'.", call. = FALSE)

  ## Raise an error if the output format and file name extension don't match
  if (output_format == "html_document" & !(file_ext(output_file) == "html"))
    stop("File name extension of output_file doesn't agree with the ",
         "output_format, should be .html.", call. = FALSE)
  if (output_format == "pdf_document" & !(file_ext(output_file) == "pdf"))
    stop("File name extension of output_file doesn't agree with the ",
         "output_format, should be .pdf.", call. = FALSE)

  ## Raise an error if dds_list is not a list of DESeqDataSets
  if (class(dds_list) != "list")
    stop("dds_list must be a list.", call. = FALSE)
  if (!all(sapply(dds_list, class) == "DESeqDataSet"))
    stop("All elements of dds_list must be DESeqDataSet objects. See the ",
         "DESeq2 Bioconductor package ",
         "(http://bioconductor.org/packages/release/bioc/html/DESeq2.html) ",
         "for more information about the DESeqDataSet class.", call. = FALSE)

  ## Raise an error if dds_list is not named
  if (length(setdiff(unique(names(dds_list)),
                     c("", NA, NULL))) != length(dds_list))
    stop("The dds_list must be a named list, ",
         "with a unique name for each element.", call. = FALSE)

  ## Raise an error or warnings if any output files/directories already exist
  if (output_format == "pdf_document") exts <- ".pdf"
  if (output_format == "html_document") exts <- ".html"

  ## Report
  fe <- paste0(output_dir, "/", basename(output_file))
  if (file.exists(fe)) {
    if (!force_overwrite) {
      stop("The file ", fe,
           " already exists. Please remove or rename the file, provide ",
           "another value of output_file, or set force_overwrite = TRUE.",
           call. = FALSE)
    } else {
      warning("The file ", fe,
              " already exists and will be overwritten, since ",
              "force_overwrite = TRUE.", immediate. = TRUE, call. = FALSE)
    }
  }

  ## ggplot2 objects
  fe <- paste0(output_dir, "/", gsub(exts, "_ggplots.rds",
                                     basename(output_file)))
  if (file.exists(fe)) {
    if (save_ggplots & !force_overwrite) {
      stop("The file ", fe,
           ", where the ggplot objects will be generated, already exists. ",
           "Please remove or rename the file, provide another value ",
           "of output_file, or set force_overwrite = TRUE.", call. = FALSE)
    } else if (save_ggplots & force_overwrite) {
      warning("The file ", fe,
              ", where the ggplot objects will be generated, already exists ",
              "and will be overwritten, since force_overwrite = TRUE.",
              immediate. = TRUE, call. = FALSE)
    } else if (!save_ggplots) {
      warning("The file ", fe, " already exists, but will not be overwritten ",
              "since save_ggplots = FALSE. Note that the ggplots in the ",
              "existing file will not correspond to the newly ",
              "generated report.", immediate. = TRUE, call. = FALSE)
    }
  }

  ## Path to the template file
  if (is.null(rmd_template)) {
    template_file <- system.file("extdata",
                                 "simulation_comparison_template_list.Rmd",
                                 package = "countsimQC")
  } else {
    template_file <- rmd_template
  }
  if (file.exists(template_file)) {
    if (file.exists(paste0(output_dir, "/", gsub(exts, ".Rmd", output_file)))) {
      stop("There is already an .Rmd file ",
           paste0(output_dir, "/", gsub(exts, ".Rmd", output_file)),
           ". Please remove or rename this file, or choose another ",
           "output_file name.", call. = FALSE)
    } else {
      file.copy(from = template_file,
                to = paste0(output_dir, "/", gsub(exts, ".Rmd", output_file)),
                overwrite = FALSE)
      template_file <- paste0(output_dir, "/", gsub(exts, ".Rmd", output_file))
    }
  } else {
    stop("The intended Rmd template file ", template_file, " does not exist.",
         call. = FALSE)
  }

  ## ------------------------------------------------------------------------ ##
  ## ----------------------- Process the arguments -------------------------- ##
  ## ------------------------------------------------------------------------ ##

  args <- list(...)
  args$input <- template_file
  args$output_format <- output_format
  args$output_file <- output_file

  ## ------------------------------------------------------------------------ ##
  ## ------------------------ Render the report ----------------------------- ##
  ## ------------------------------------------------------------------------ ##

  output_file <- do.call("render", args = args)

  ## ------------------------------------------------------------------------ ##
  ## --------------------- Remove temporary file ---------------------------- ##
  ## ------------------------------------------------------------------------ ##

  file.remove(template_file)

  invisible(output_file)
}

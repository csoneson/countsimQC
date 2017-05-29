#' Generate countsimQC report
#'
#' Generate a report comparing a range of characteristics across a collection of
#' count data sets.
#'
#' @param dds_list Named list of \code{DESeqDataSets} to compare. See the
#'   \code{DESeq2} Bioconductor package
#'   (http://bioconductor.org/packages/release/bioc/html/DESeq2.html) for more
#'   information about the \code{DESeqDataSet} class. Each \code{DESeqDataSet}
#'   object in the list should contain a count matrix, a data frame with sample
#'   information and a design formula. The sample information and design formula
#'   will be used to calculate dispersions appropriately.
#' @param show_code Whether or not to include the code in the final report.
#' @param output_format The output format of the report. If set to NULL or
#'   "html_document", an html report will be generated. If set to
#'   "pdf_document", a pdf report will be generated.
#' @param output_file The file name of the final report.
#' @param output_dir The directory where the final report should be saved.
#' @param rmd_template The Rmarkdown (.Rmd) file that will be used as the
#'   template for generating the report. If set to NULL (default), the template
#'   provided with the countsimQC package will be used. See Details for more
#'   information.
#' @param force_overwrite Whether to force overwrite existing output files when
#'   saving the generated report and figures.
#' @param save_ggplots Whether to save the ggplot objects for all the output
#'   figures, to allow additional fine-tuning and generation of individual
#'   plots. Note, that the resulting file can be quite large, especially when
#'   many or large data sets are compared.
#' @param ... Other arguments that will be passed to \code{rmarkdown::render}.
#'
#' @author Charlotte Soneson, using a template from Nicholas Hamilton
#'   (http://stackoverflow.com/questions/37097535/generate-report-in-r)
#'
#' @details When the function is called, the template file (specified by
#'   \code{rmd_template}) will be copied into the output folder, and
#'   \code{rmarkdown::render} will be called to generate the final report. If
#'   there is already a .Rmd file with the same name in the output folder, the
#'   function will raise an error and stop, to avoid overwriting the existing
#'   file. The reason for this behaviour is that the copied template in the
#'   output folder will be deleted once the report is generated.
#'
#' @export
#' @importFrom rmarkdown render
#' @importFrom tools file_ext file_path_sans_ext
#' @import DESeq2 edgeR dplyr tidyr ggplot2
#' @return No value is returned, but a report is generated in the
#'   \code{output_dir} directory.
#'
countsimQC_report <- function(dds_list, show_code = FALSE, output_format = NULL,
                              output_file, output_dir = "./",
                              rmd_template = NULL, force_overwrite = FALSE,
                              save_ggplots = FALSE, ...){
  ## This function was inspired by code from Nicholas Hamilton, provided at
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r

  if (is.null(output_format)) output_format <- "html_document"

  ## ------------------------------------------------------------------------ ##
  ## --------------------- Check input arguments ---------------------------- ##
  ## ------------------------------------------------------------------------ ##

  ## ------------------------ output_format --------------------------------- ##
  ## Raise an error if output_format is not one of the allowed
  if (!(output_format %in% c("pdf_document", "html_document")))
    stop("The provided output_format is currently not supported. Please use ",
         "either 'html_document' (or NULL) or 'pdf_document'.", call. = FALSE)

  ## Raise an error if the output format and file name extension don't match
  if (output_format != paste0(file_ext(output_file), "_document"))
    stop(paste0("File name extension of output_file doesn't agree with the ",
                "output_format, should be .",
                gsub("_document", "", output_format)), call. = FALSE)

  ## --------------------------- dds_list ----------------------------------- ##
  ## Raise an error if dds_list is not a named list of DESeqDataSets
  if (class(dds_list) != "list")
    stop("dds_list must be a list.", call. = FALSE)
  if (!all(sapply(dds_list, class) == "DESeqDataSet"))
    stop("All elements of dds_list must be DESeqDataSet objects. See the ",
         "DESeq2 Bioconductor package ",
         "(http://bioconductor.org/packages/release/bioc/html/DESeq2.html) ",
         "for more information about the DESeqDataSet class.", call. = FALSE)
  if (length(setdiff(unique(names(dds_list)),
                     c("", NA, NULL))) != length(dds_list))
    stop("The dds_list must be a named list, ",
         "with a unique name for each element.", call. = FALSE)

  ## ------------------------- output files --------------------------------- ##
  output_report <- paste0(output_dir, "/", basename(output_file))
  output_ggplots <- paste0(output_dir, "/",
                           file_path_sans_ext(basename(output_file)), "_ggplots.rds")
  output_rmd <- paste0(output_dir, "/",
                       file_path_sans_ext(basename(output_file)), ".Rmd")

  ## Report
  if (file.exists(output_report)) {
    if (!force_overwrite) {
      stop("The file ", output_report,
           " already exists. Please remove or rename the file, provide ",
           "another value of output_file, or set force_overwrite = TRUE.",
           call. = FALSE)
    } else {
      warning("The file ", output_report,
              " already exists and will be overwritten, since ",
              "force_overwrite = TRUE.", immediate. = TRUE, call. = FALSE)
    }
  }

  ## ggplot2 objects
  if (file.exists(output_ggplots)) {
    if (save_ggplots & !force_overwrite) {
      stop("The file ", output_ggplots,
           ", where the ggplot objects will be generated, already exists. ",
           "Please remove or rename the file, provide another value ",
           "of output_file, or set force_overwrite = TRUE.", call. = FALSE)
    } else if (save_ggplots & force_overwrite) {
      warning("The file ", output_ggplots,
              ", where the ggplot objects will be generated, already exists ",
              "and will be overwritten, since force_overwrite = TRUE.",
              immediate. = TRUE, call. = FALSE)
    } else if (!save_ggplots) {
      warning("The file ", output_ggplots,
              " already exists, but will not be overwritten ",
              "since save_ggplots = FALSE. Note that the ggplots in the ",
              "existing file will not correspond to the newly ",
              "generated report.", immediate. = TRUE, call. = FALSE)
    }
  }

  ## ------------------------- Rmd template --------------------------------- ##
  ## Path to the template file
  if (is.null(rmd_template)) {
    template_file <- system.file("extdata",
                                 "simulation_comparison_template_list.Rmd",
                                 package = "countsimQC")
  } else {
    template_file <- rmd_template
  }
  if (file.exists(template_file)) {
    if (file.exists(output_rmd)) {
      stop("There is already an .Rmd file ", output_rmd,
           ". Please remove or rename this file, or choose another ",
           "output_file name.", call. = FALSE)
    } else {
      file.copy(from = template_file, to = output_rmd, overwrite = FALSE)
    }
  } else {
    stop("The intended Rmd template file ", template_file, " does not exist.",
         call. = FALSE)
  }

  ## ------------------------------------------------------------------------ ##
  ## ----------------------- Process the arguments -------------------------- ##
  ## ------------------------------------------------------------------------ ##

  args <- list(...)
  args$input <- output_rmd
  args$output_format <- output_format
  args$output_file <- output_file

  ## ------------------------------------------------------------------------ ##
  ## ------------------------ Render the report ----------------------------- ##
  ## ------------------------------------------------------------------------ ##

  output_file <- do.call("render", args = args)

  ## ------------------------------------------------------------------------ ##
  ## --------------------- Remove temporary file ---------------------------- ##
  ## ------------------------------------------------------------------------ ##

  file.remove(output_rmd)

  invisible(output_file)
}

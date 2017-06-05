#' Generate countsimQC report
#'
#' Generate a report comparing a range of characteristics across a collection of
#' one or more count data sets.
#'
#' @param ddsList Named list of \code{DESeqDataSets} to compare. See the
#'   \code{DESeq2} Bioconductor package
#'   (http://bioconductor.org/packages/release/bioc/html/DESeq2.html) for more
#'   information about the \code{DESeqDataSet} class. Each \code{DESeqDataSet}
#'   object in the list should contain a count matrix, a data frame with sample
#'   information and a design formula. The sample information and design formula
#'   will be used to calculate dispersions appropriately.
#' @param outputFile The file name of the final report. The extension must match
#'   the selected \code{outputFormat} (i.e., either .html or .pdf).
#' @param outputDir The directory where the final report should be saved.
#' @param outputFormat The output format of the report. If set to NULL or
#'   "html_document", an html report will be generated. If set to
#'   "pdf_document", a pdf report will be generated.
#' @param showCode Whether or not to include the code in the final report.
#' @param rmdTemplate The Rmarkdown (.Rmd) file that will be used as the
#'   template for generating the report. If set to NULL (default), the template
#'   provided with the countsimQC package will be used. See Details for more
#'   information.
#' @param forceOverwrite Whether to force overwrite existing output files when
#'   saving the generated report and figures.
#' @param savePlots Whether to save the ggplot objects for all the output
#'   figures, to allow additional fine-tuning and generation of individual
#'   plots. Note that the resulting file can be quite large, especially when
#'   many and/or large data sets are compared.
#' @param description A string (of arbitrary length) describing the content of
#'   the generated report. This will be included in the beginning of the report.
#'   If set to NULL, a default description listing the number and names of the
#'   included data sets will be used.
#' @param ... Other arguments that will be passed to \code{rmarkdown::render}.
#'
#' @author Charlotte Soneson
#'
#' @details When the function is called, the template file (specified by
#'   \code{rmdTemplate}) will be copied into the output folder, and
#'   \code{rmarkdown::render} will be called to generate the final report. If
#'   there is already a .Rmd file with the same name in the output folder, the
#'   function will raise an error and stop, to avoid overwriting the existing
#'   file. The reason for this behaviour is that the copied template in the
#'   output folder will be deleted once the report is generated.
#'
#' @export
#' @importFrom rmarkdown render
#' @importFrom tools file_ext file_path_sans_ext
#' @import SummarizedExperiment DESeq2 edgeR dplyr tidyr ggplot2
#' @return No value is returned, but a report is generated in the
#'   \code{outputDir} directory.
#'
#' @examples
#' \dontrun{
#' data(countsimExample)
#' countsimQCReport(countsimExample, outputDir = "./",
#'                  outputFile = "example.html")
#' }
#'
countsimQCReport <- function(ddsList, outputFile, outputDir = "./",
                             outputFormat = NULL, showCode = FALSE,
                             rmdTemplate = NULL, forceOverwrite = FALSE,
                             savePlots = FALSE, description = NULL, ...){
  ## This function was inspired by code from Nicholas Hamilton, provided at
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r

  if (is.null(outputFormat)) outputFormat <- "html_document"

  ## ------------------------------------------------------------------------ ##
  ## --------------------- Check input arguments ---------------------------- ##
  ## ------------------------------------------------------------------------ ##

  ## ------------------------ outputFormat --------------------------------- ##
  ## Raise an error if outputFormat is not one of the allowed
  if (!(outputFormat %in% c("pdf_document", "html_document")))
    stop("The provided outputFormat is currently not supported. Please use ",
         "either 'html_document' (or NULL) or 'pdf_document'.", call. = FALSE)

  ## Raise an error if the output format and file name extension don't match
  if (outputFormat != paste0(tools::file_ext(outputFile), "_document"))
    stop(paste0("File name extension of outputFile doesn't agree with the ",
                "outputFormat, should be .",
                gsub("_document$", "", outputFormat)), call. = FALSE)

  ## --------------------------- ddsList ----------------------------------- ##
  ## Raise an error if ddsList is not a named list of DESeqDataSets
  if (class(ddsList) != "list")
    stop("ddsList must be a list.", call. = FALSE)
  if (!all(sapply(ddsList, class) == "DESeqDataSet"))
    stop("All elements of ddsList must be DESeqDataSet objects. See the ",
         "DESeq2 Bioconductor package ",
         "(http://bioconductor.org/packages/release/bioc/html/DESeq2.html) ",
         "for more information about the DESeqDataSet class.", call. = FALSE)
  if (length(setdiff(unique(names(ddsList)),
                     c("", NA, NULL))) != length(ddsList))
    stop("The ddsList must be a named list, ",
         "with a unique name for each element.", call. = FALSE)

  ## ------------------------- output files --------------------------------- ##
  outputReport <- paste0(outputDir, "/", basename(outputFile))
  outputPlots <- paste0(outputDir, "/",
                           tools::file_path_sans_ext(basename(outputFile)),
                        "_ggplots.rds")
  outputRmd <- paste0(outputDir, "/",
                       tools::file_path_sans_ext(basename(outputFile)), ".Rmd")

  ## Report
  if (file.exists(outputReport)) {
    if (!forceOverwrite) {
      stop("The file ", outputReport,
           " already exists. Please remove or rename the file, provide ",
           "another value of outputFile, or set forceOverwrite = TRUE.",
           call. = FALSE)
    } else {
      warning("The file ", outputReport,
              " already exists and will be overwritten, since ",
              "forceOverwrite = TRUE.", immediate. = TRUE, call. = FALSE)
    }
  }

  ## ggplot2 objects
  if (file.exists(outputPlots)) {
    if (savePlots & !forceOverwrite) {
      stop("The file ", outputPlots,
           ", where the ggplot objects will be generated, already exists. ",
           "Please remove or rename the file, provide another value ",
           "of outputFile, or set forceOverwrite = TRUE.", call. = FALSE)
    } else if (savePlots & forceOverwrite) {
      warning("The file ", outputPlots,
              ", where the ggplot objects will be generated, already exists ",
              "and will be overwritten, since forceOverwrite = TRUE.",
              immediate. = TRUE, call. = FALSE)
    } else if (!savePlots) {
      warning("The file ", outputPlots,
              " already exists, but will not be overwritten ",
              "since savePlots = FALSE. Note that the ggplots in the ",
              "existing file may not correspond to the newly ",
              "generated report.", immediate. = TRUE, call. = FALSE)
    }
  }

  ## ------------------------- Rmd template --------------------------------- ##
  ## Path to the template file
  if (is.null(rmdTemplate)) {
    templateFile <- system.file("extdata",
                                 "simulation_comparison_template_list.Rmd",
                                 package = "countsimQC")
  } else {
    templateFile <- rmdTemplate
  }
  if (file.exists(templateFile)) {
    if (file.exists(outputRmd)) {
      stop("There is already an .Rmd file ", outputRmd,
           ". Please remove or rename this file, or choose another ",
           "outputFile name.", call. = FALSE)
    } else {
      file.copy(from = templateFile, to = outputRmd, overwrite = FALSE)
    }
  } else {
    stop("The intended Rmd template file ", templateFile, " does not exist.",
         call. = FALSE)
  }

  ## ------------------------------------------------------------------------ ##
  ## ----------------- Add description if not provided ---------------------- ##
  ## ------------------------------------------------------------------------ ##

  if (is.null(description)) {
    description <- sprintf("This report shows the characteristics of %s %s,
                           named %s.",
                           length(ddsList),
                           ifelse(length(ddsList) == 1, "data set", "data sets"),
                           paste(names(ddsList), collapse = ", "))
  }

  ## ------------------------------------------------------------------------ ##
  ## ----------------------- Process the arguments -------------------------- ##
  ## ------------------------------------------------------------------------ ##

  args <- list(...)
  args$input <- outputRmd
  args$output_format <- outputFormat
  args$output_file <- outputFile

  ## ------------------------------------------------------------------------ ##
  ## ------------------------ Render the report ----------------------------- ##
  ## ------------------------------------------------------------------------ ##

  outputFile <- do.call("render", args = args)

  ## ------------------------------------------------------------------------ ##
  ## --------------------- Remove temporary file ---------------------------- ##
  ## ------------------------------------------------------------------------ ##

  file.remove(outputRmd)

  invisible(outputFile)
}

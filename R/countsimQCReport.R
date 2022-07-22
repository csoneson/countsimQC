#' Check whether pandoc and pandoc-citeproc is available
#'
#' @author Charlotte Soneson
#'
#' @param ignorePandoc logical. If TRUE, just give a warning if pandoc or
#'   pandoc-citeproc is not available. If FALSE, stop.
#'
#' @keywords internal
#'
#' @return A logical(1), indicating whether pandoc can be run or not.
#'   In addition, raises either a warning or an error (depending on the
#'   value of \code{ignorePandoc}) if pandoc or pandoc-citeproc is not
#'   available.
#'
#' @importFrom rmarkdown pandoc_available pandoc_exec
#'
.checkPandoc <- function(ignorePandoc) {
  ## Initialize output to TRUE
  doRender <- TRUE

  ## First check whether pandoc is available
  if (!rmarkdown::pandoc_available()) {
    doRender <- FALSE
    ## If pandoc is not available, either give a warning or an error,
    ## depending on the value of ignorePandoc
    if (ignorePandoc) {
      ## If ignorePandoc is TRUE, just give a warning
      warning("pandoc is not available! ",
              "The final report will not be generated.",
              immediate. = TRUE)
    } else {
      ## If ignorePandoc is FALSE, stop
      stop("pandoc is not available!")
    }
  } else {
    ## If pandoc is available, check for pandoc-citeproc
    ## Only do this if the pandoc version is <2.11, since
    ## pandoc-citeproc is not included (or needed) in v2.11 and later.
    if (!rmarkdown::pandoc_available(version = "2.11")) {
      ## TRUE if the available pandoc version is not 2.11 or newer
      ## pandoc-citeproc should be found in the path, or in the
      ## same folder as the pandoc executable
      if (Sys.which("pandoc-citeproc") == "" &&
          !file.exists(file.path(dirname(rmarkdown::pandoc_exec()),
                                 "pandoc-citeproc"))) {
        doRender <- FALSE
        ## pandoc-citeproc is required, but not found
        if (ignorePandoc) {
          ## If ignorePandoc is TRUE, just give a warning
          warning("pandoc-citeproc is not available! ",
                  "The final report will not be generated.",
                  immediate. = TRUE)
        } else {
          ## If ignorePandoc is FALSE, stop
          stop("pandoc-citeproc is not available!")
        }
      }
    }
  }
  return(doRender)
}

#' Generate countsimQC report
#'
#' Generate a report comparing a range of characteristics across a collection of
#' one or more count data sets.
#'
#' @param ddsList Named list of \code{DESeqDataSets} or count matrices to
#'   compare. See the \code{DESeq2} Bioconductor package
#'   (http://bioconductor.org/packages/release/bioc/html/DESeq2.html) for more
#'   information about the \code{DESeqDataSet} class. Each \code{DESeqDataSet}
#'   object in the list should contain a count matrix, a data frame with sample
#'   information and a design formula. The sample information and design formula
#'   will be used to calculate dispersions appropriately. If count matrices are
#'   provided, it is assumed that all columns represent replicate samples, and
#'   the design formula ~1 will be used.
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
#' @param maxNForCorr The maximal number of samples (features) for which
#'   pairwise correlation coefficients will be calculated. If the number of
#'   samples (features) exceeds this number, they will be randomly subsampled.
#' @param maxNForDisp The maximal number of samples that will be used to
#'   estimate dispersions. By default, all samples are used. This can be lowered
#'   to speed up calculations (and obtain approximate results) for large data
#'   sets.
#' @param calculateStatistics Whether to calculate quantitative pairwise
#'   statistics for comparing data sets in addition to generating the plots.
#' @param subsampleSize The number of randomly selected observations (samples,
#'   features or pairs of samples or features) for which certain
#'   (time-consuming) statistics will be calculated. Only used if
#'   \code{calculateStatistics} = TRUE.
#' @param kmin,kfrac For statistics that require the extraction of the k nearest
#'   neighbors of a given point, the number of neighbors will be max(kmin, kfrac
#'   * nrow(df))
#' @param permutationPvalues Whether to calculate permutation p-values for
#'   selected pairwise data set comparison statistics.
#' @param nPermutations The number of permutations to perform when calculating
#'   permutation p-values for data set comparison statistics. Only used if
#'   \code{permutationPvalues} = TRUE.
#' @param knitrProgress Whether to show the progress bar when the report is
#'   generated.
#' @param quiet Whether to suppress warnings and progress messages when the
#'   report is generated.
#' @param ignorePandoc Determines what to do if pandoc or pandoc-citeproc is
#'   missing (if Sys.which("pandoc") or Sys.which("pandoc-citeproc") is ""). If
#'   ignorePandoc is TRUE, only a warning is given. The figures will be
#'   generated, but not the final report. If ignorePandoc is FALSE (default),
#'   the execution stops immediately.
#' @param useRAGG Logical scalar, indicating whether to use ragg_png as the
#'   graphics device in the report rather than the default png.
#' @param dpi Numeric scalar, setting the dpi of the generated plots. Only
#'   used if \code{useRAGG} is \code{TRUE}.
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
#' @importFrom caTools trapz
#' @importFrom randtests runs.test
#' @importFrom methods is
#' @import SummarizedExperiment DESeq2 edgeR dplyr tidyr ggplot2
#' @import genefilter DT GenomeInfoDbData ragg
#' @return No value is returned, but a report is generated in the
#'   \code{outputDir} directory.
#'
#' @examples
#' ## Load example data
#' data(countsimExample)
#' \dontrun{
#' ## Generate report
#' countsimQCReport(countsimExample, outputDir = "./",
#'                  outputFile = "example.html")
#' }
#'
countsimQCReport <- function(ddsList, outputFile, outputDir = "./",
                             outputFormat = NULL, showCode = FALSE,
                             rmdTemplate = NULL, forceOverwrite = FALSE,
                             savePlots = FALSE, description = NULL,
                             maxNForCorr = 500, maxNForDisp = Inf,
                             calculateStatistics = TRUE, subsampleSize = 500,
                             kfrac = 0.01, kmin = 5,
                             permutationPvalues = FALSE, nPermutations = NULL,
                             knitrProgress = FALSE, quiet = FALSE,
                             ignorePandoc = FALSE, useRAGG = FALSE,
                             dpi = 96, ...){
  ## This function was inspired by code from Nicholas Hamilton, provided at
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r

  if (is.null(outputFormat)) outputFormat <- "html_document"

  doRender <- .checkPandoc(ignorePandoc)

  ## ------------------------------------------------------------------------ ##
  ## --------------------- Check input arguments ---------------------------- ##
  ## ------------------------------------------------------------------------ ##

  ## ------------------------ outputFormat --------------------------------- ##
  ## Raise an error if outputFormat is not one of the allowed
  if (!(outputFormat %in% c("pdf_document", "html_document"))) {
    stop("The provided outputFormat is currently not supported. Please use ",
         "either 'html_document' (or NULL) or 'pdf_document'.", call. = FALSE)
  }

  ## Raise an error if the output format and file name extension don't match
  if (outputFormat != paste0(tools::file_ext(outputFile), "_document")) {
    stop(paste0("File name extension of outputFile doesn't agree with the ",
                "outputFormat, should be .",
                gsub("_document$", "", outputFormat)), call. = FALSE)
  }

  ## --------------------------- ddsList ----------------------------------- ##
  ## Raise an error if ddsList is not a named list of DESeqDataSets/data
  ## frames/matrices
  if (!is(ddsList, "list")) {
    stop("ddsList must be a list.", call. = FALSE)
  }
  if (length(setdiff(unique(names(ddsList)),
                     c("", NA, NULL))) != length(ddsList)) {
    stop("ddsList must be a named list, ",
         "with a unique name for each element.", call. = FALSE)
  }
  if (!all(vapply(ddsList, function(w) {
    is(w, "DESeqDataSet") | is(w, "data.frame") | is(w, "matrix")
  }, FALSE))) {
    stop("All elements of ddsList must be DESeqDataSet objects, data.frames ",
         "or matrices. See the DESeq2 Bioconductor package ",
         "(http://bioconductor.org/packages/release/bioc/html/DESeq2.html) ",
         "for more information about the DESeqDataSet class.", call. = FALSE)
  }

  ## If some objects are data frames or matrices, convert them into
  ## DESeqDataSets
  ddsList <- lapply(ddsList, function(ds) {
    if (is(ds, "DESeqDataSet")) {
      ds
    } else {
      DESeq2::DESeqDataSetFromMatrix(
        countData = round(as.matrix(ds)),
        colData = data.frame(sample = seq_len(ncol(ds))),
        design = ~ 1)
    }
  })
  stopifnot(all(vapply(ddsList, function(w) is(w, "DESeqDataSet"), FALSE)))

  ## ------------------------- output files --------------------------------- ##
  outputReport <- file.path(outputDir, basename(outputFile))
  outputPlots <- file.path(
    outputDir,
    paste0(tools::file_path_sans_ext(basename(outputFile)), "_ggplots.rds"))
  outputRmd <- file.path(
    outputDir,
    paste0(tools::file_path_sans_ext(basename(outputFile)), ".Rmd"))

  ## Report
  if (file.exists(outputReport)) {
    if (!forceOverwrite) {
      stop("The file ", outputReport,
           " already exists. Please remove or rename the file, provide ",
           "another value of outputFile, or set forceOverwrite = TRUE.",
           call. = FALSE)
    } else {
      if (!quiet) {
        warning("The file ", outputReport,
                " already exists and will be overwritten, since ",
                "forceOverwrite = TRUE.", immediate. = TRUE, call. = FALSE)
      }
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
      if (!quiet) {
        warning("The file ", outputPlots,
                ", where the ggplot objects will be generated, already exists ",
                "and will be overwritten, since forceOverwrite = TRUE.",
                immediate. = TRUE, call. = FALSE)
      }
    } else if (!savePlots) {
      if (!quiet) {
        warning("The file ", outputPlots,
                " already exists, but will not be overwritten ",
                "since savePlots = FALSE. Note that the ggplots in the ",
                "existing file may not correspond to the newly ",
                "generated report.", immediate. = TRUE, call. = FALSE)
      }
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
    description <- sprintf(
      "This report shows the characteristics of %s %s, named %s.",
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
  args$quiet <- !knitrProgress
  args$run_pandoc <- doRender

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

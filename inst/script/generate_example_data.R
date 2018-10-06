suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(splatter)
  library(DESeq2)
})

mae <- readRDS(url("http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE74596.rds"))
config <- list(groupid = "source_name_ch1",
               keepgroups = c("Single_cell_RNA-seq_NKT0", "Single_cell_RNA-seq_NKT17")
)
mae <- updateObject(mae)
pdata <- colData(mae)
groupid <- config$groupid
keepgroups <- config$keepgroups

counts <- assays(experiments(mae)[["gene"]])[["count_lstpm"]]
stopifnot(all(colnames(counts) == rownames(pdata)))
keep <- which(pdata[, groupid] %in% keepgroups)
counts <- round(counts[, keep])
countsnz <- counts[rowSums(counts > 0) > 1, ]
countsz <- counts[rowSums(counts > 0) <= 1, ]
group <- as.character(pdata[keep, groupid])

dds_orig <- DESeqDataSetFromMatrix(
  countData = countsnz,
  colData = data.frame(group = group, sample = colnames(countsnz),
                       row.names = colnames(countsnz),
                       stringsAsFactors = FALSE),
  design = ~ group)

## splat
params <- splatEstimate(countsnz)
sim <- splatSimulate(params, method = "groups", group.prob = c(table(group))/length(group))
dds_splat <- DESeqDataSetFromMatrix(countData = counts(sim),
                                    colData = colData(sim),
                                    design = ~ Group)

## Lun
params <- lunEstimate(countsnz)
sim <- lunSimulate(params, groupCells = c(table(group)))
dds_lun <- DESeqDataSetFromMatrix(countData = counts(sim),
                                  colData = colData(sim),
                                  design = ~ Group)

countsimExample <- list(Original = dds_orig[1:10000, c(1:5, 84:89)],
                        Sim1 = dds_splat[1:10000, c(1:5, 84:89)],
                        Sim2 = dds_lun[1:10000, c(1:5, 84:89)])
devtools::use_data(countsimExample, pkg = "../..", overwrite = TRUE)

countsimExample_dfmat <- list(Original = dds_orig[1:10000, c(1:5, 84:89)],
                              Sim1 = as.matrix(counts(dds_splat[1:10000, c(1:5, 84:89)])),
                              Sim2 = as.data.frame(counts(dds_lun[1:10000, c(1:5, 84:89)])))
devtools::use_data(countsimExample_dfmat, pkg = "../..", overwrite = TRUE)

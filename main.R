rm(list = ls())
remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"
test.folder <- c("GSE171330")
opt.method <- "manual"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(dplyr)
library(harmony)
library(ggplot2)
library(SingleR)
library(BiocParallel)

for (folder in test.folder) {
  fname <-
    paste(saved.dir,
      paste(folder, opt.method, "macrophage.RDdata", sep = "_"),
      sep = "/"
    )
  obj <- readRDS(file = fname)
  # normalize/scale data
  # gene expression between different conditions
}

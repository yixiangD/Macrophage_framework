remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"

library(Biobase)
library(GEOquery)
gset <- getGEO("GSE178209", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

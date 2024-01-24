remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package

data.10x <- readxl::read_excel(paste(local.path, "secondary/10x.xlsx", sep = "/"))
data.10x$type <- gsub(".mtx.gz", "", data.10x$type)
data.10x$type <- gsub(".tsv.gz", "", data.10x$type)

all10x.list <- list()
for (f in unique(data.10x$folder)) {
  data.subset <- data.10x[data.10x$folder == f,]
  seurat.list <- list()
  for (samp in unique(data.subset$sample)) {
    data.samp <- data.subset[data.subset$sample == samp, ]
    data.samp$full.path <- paste(data.samp$folder, data.samp$file, sep = "/")
    mtx <- paste(local.path, data.samp$full.path[grepl("matrix", data.samp$full.path)], sep = "/")
    cells <- paste(local.path, data.samp$full.path[grepl("barcode", data.samp$full.path)], sep = "/")
    features <- paste(local.path, data.samp$full.path[grepl("feature|gene", data.samp$full.path)], sep = "/")
    # reading
    counts <- ReadMtx(mtx = mtx, cells = cells, features = features)
    obj <- CreateSeuratObject(counts = counts, project = samp)
    seurat.list <- append(seurat.list, obj)
  }
  names(seurat.list) <- unique(data.subset$sample)
  all10x.list <- append(all10x.list, seurat.list)
}
names(all10x.list) <- unique(data.10x$folder)
# Retrieve datasets
saveRDS(all10x.list, file = paste(local.path, "secondary/all10x.rds", sep = "/"))
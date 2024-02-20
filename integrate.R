rm(list = ls())
remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(dplyr)
library(harmony)

data.10x <- readxl::read_excel(paste(local.path, "secondary/10x.xlsx", sep = "/"))
data.10x$type <- gsub(".mtx.gz", "", data.10x$type)
data.10x$type <- gsub(".tsv.gz", "", data.10x$type)

df <- data.10x %>%
  mutate(
    splits = strsplit(sample, "_")
  ) %>%
  rowwise() %>%
  mutate(
    code = splits[1],
    case = paste(splits[2:length(splits)], collapse = "_")
  ) %>%
  select(-splits)

df <- df %>%
  group_by(folder, type) %>%
  mutate(id = row_number())

df <- df[order(df$folder, df$id), ]

test.folder <- c("GSE171330")

for (folder in test.folder) {
  df.sub <- df[df$folder == folder, ]
  example.path <- paste(local.path, folder, sep = "/")
  obj.list <- c()
  for (samp in unique(df.sub$sample)) {
    print(samp)
    mtx <- paste(example.path, df.sub$file[grepl("matrix", df.sub$file) & grepl(samp, df.sub$file)], sep = "/")
    cells <- paste(example.path, df.sub$file[grepl("barcode", df.sub$file) & grepl(samp, df.sub$file)], sep = "/")
    features <- paste(example.path, df.sub$file[grepl("feature|gene", df.sub$file) & grepl(samp, df.sub$file)], sep = "/")
    counts <- ReadMtx(mtx = mtx, cells = cells, features = features)
    seurat_object <- CreateSeuratObject(counts = counts, project = samp)
    seurat_object$condition <- unique(df.sub$case[grepl(samp, df.sub$file)])
    obj.list <- c(obj.list, seurat_object)
  }
  obj.final <- merge(obj.list[1][[1]], obj.list[2][[1]])
  obj.final <- NormalizeData(obj.final) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)
  obj.final <- RunHarmony(obj.final, group.by.vars = "condition")
  obj.final <- RunUMAP(obj.final, reduction = "harmony", dims = 1:30)
  obj.final <- FindNeighbors(obj.final, reduction = "harmony", dims = 1:30) %>% FindClusters()
  DimPlot(obj.final, group.by = c("condition", "ident", "seurat_annotations"), ncol = 3)
}

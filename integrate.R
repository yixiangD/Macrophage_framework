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
    seurat_object <- NormalizeData(seurat_object) %>%
      FindVariableFeatures()

    obj.list <- c(obj.list, seurat_object)
  }
  obj.final <- merge(obj.list[1][[1]], obj.list[2][[1]])
  obj.final <- obj.final %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)
  obj.final <- RunHarmony(obj.final, group.by.vars = "condition")
  obj.final <- RunUMAP(obj.final, reduction = "harmony", dims = 1:30)
  obj.final <- FindNeighbors(obj.final, reduction = "harmony", dims = 1:30) %>% FindClusters()
  # DimPlot(obj.final, group.by = c("condition", "ident", "seurat_annotations"), ncol = 3)
  DimPlot(obj.final, reduction = "umap", split.by = "condition")
  obj.final <- JoinLayers(obj.final)
  obj.final <- FindClusters(obj.final, resolution = 0.5)
}

# find markers----
# manual
# markers <- FindConservedMarkers(obj.final, ident.1 = 6, grouping.var = "condition", verbose = FALSE)

library(celldex)
library(SingleR)
library(BiocParallel)

ref <- celldex::MouseRNAseqData()
sceObject <- as.SingleCellExperiment(obj.final)
singleRResults <- SingleR(test = sceObject, ref = ref, labels = ref$label.main, clusters = obj.final$seurat_clusters)
# Viewing the first few annotations
head(singleRResults$labels)
# Summary of predicted cell types
table(singleRResults$labels)

scores <- as.data.frame(singleRResults@listData$scores)

plotScoreHeatmap(singleRResults, clusters = singleRResults$labels, fontsize.row = 9, show_colnames = T)

new.cluster.ids <- singleRResults$pruned.labels
names(new.cluster.ids) <- levels(obj.final)
levels(obj.final)
obj.final <- RenameIdents(obj.final, new.cluster.ids)

UMAPPlot(object = obj.final, pt.size = 0.5, label = TRUE)

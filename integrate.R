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
lower_bound <- 200
upper_bound <- 2500
mito_threshold <- 5

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
    seurat_object[["percent.mito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

    fig <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
    print(fig)

    seurat_object <- subset(seurat_object, subset = nFeature_RNA > lower_bound & nFeature_RNA < upper_bound & percent.mito < mito_threshold)

    seurat_object$condition <- unique(df.sub$case[grepl(samp, df.sub$file)])

    seurat_object <- NormalizeData(seurat_object) %>%
      FindVariableFeatures()

    obj.list <- c(obj.list, seurat_object)
    # dim.plot <- DimPlot(seurat_object, reduction = "umap")
    # print(dim.plot)
  }
  obj.final <- merge(obj.list[1][[1]], obj.list[2][[1]])
  obj.final <- RunHarmony(obj.final, group.by.vars = "condition")

  obj.final <- obj.final %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

  obj.final <- RunUMAP(obj.final, reduction = "harmony", dims = 1:30)
  obj.final <- FindNeighbors(obj.final, reduction = "harmony", dims = 1:30) %>% FindClusters()
  dim.plot <- DimPlot(obj.final, group.by = c("condition", "ident"), ncol = 2)
  # # dim.plot <- DimPlot(obj.final, reduction = "umap", split.by = "condition")
  # # print(dim.plot)
  # obj.final <- JoinLayers(obj.final)
  # obj.final <- FindClusters(obj.final, resolution = 0.5)
}

dim.plot <- DimPlot(obj.final)
print(dim.plot)

# find markers----
# manual
# markers <- FindConservedMarkers(obj.final, ident.1 = 6, grouping.var = "condition", verbose = FALSE)
# macropage features
features <- c("Adgre1", "Csf1r", "H2-Ab1", "Cd68", "Lyz2", "Itgam", "Mertk")
FeaturePlot(obj.final, features = features)

obj.small <- subset(obj.final, Adgre1 > 1 | Csf1r > 1 | `H2-Ab1` > 1 | Cd68 > 1 | Lyz2 > 1 | Itgam > 1 | Mertk > 1)
obj.small <- RunUMAP(obj.small, reduction = "harmony", dims = 1:30)
obj.small <- FindNeighbors(obj.small, reduction = "harmony", dims = 1:30) %>% FindClusters()
dim.plot.small <- DimPlot(obj.small, group.by = c("condition", "ident"), ncol = 2)

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

score.heatmap <- plotScoreHeatmap(singleRResults, clusters = singleRResults$labels, fontsize.row = 9, show_colnames = T)

new.cluster.ids <- singleRResults$pruned.labels
names(new.cluster.ids) <- levels(obj.final)
levels(obj.final)
obj.final <- RenameIdents(obj.final, new.cluster.ids)

umap.integrated <- UMAPPlot(object = obj.final, pt.size = 0.5, label = TRUE)
print(umap.integrated)

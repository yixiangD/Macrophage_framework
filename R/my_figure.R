remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"

# load Data provide by David, i.e., reference data
c0 <- readRDS(paste(remote.path, "data.rds", sep = "/"))
df.10x.loc <- "secondary/10x.xlsx"
data.10x <- readxl::read_excel(paste(local.path, df.10x.loc, sep = "/"))
data.10x$type <- gsub(".mtx.gz", "", data.10x$type)
data.10x$type <- gsub(".tsv.gz", "", data.10x$type)

file <- data.10x$folder[1]
data.subset <- data.10x[data.10x$folder == file, ]

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
  print(samp)
  print(obj)
  # seuratObject <- NormalizeData(obj)
  # seuratObject <- SCTransform(seuratObject, verbose = FALSE)  # Or use NormalizeData for normalization
  # seuratObject <- ScaleData(seuratObject, verbose = FALSE)
  # seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject), verbose = FALSE)
  # seuratObject <- RunUMAP(seuratObject, dims = 1:10)
  # seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
  # seuratObject <- FindClusters(seuratObject, resolution = 0.5)
  #
  # seurat.list <- append(seurat.list, seuratObject)
}


# pval # genes annotated with p values
# fully processed reference dataset
library(Seurat)
library(ggplot2)
seurat.list <- PrepSCTIntegration(seurat.list, verbose = FALSE)
# Assuming you have two Seurat objects: seuratObject1 and seuratObject2
features <- intersect(rownames(seurat.list[[1]]), rownames(seurat.list[[2]]))
# If you have more than two objects, continue finding the intersection with other datasets
seurat.list
max.f <- 500
n <- 1
for (i in seurat.list) {
  DefaultAssay(i) <- "RNA"
  i <- DietSeurat(i, assays = "RNA")
  if ("stage.score" %in% colnames(i@meta.data)) {
    t <- subset(i, subset = stage.score > 0.8, slot = "counts") # retain high score cells from query datasets
  } else {
    t <- subset(i, cells = sample(colnames(i), 400)) # retain 500 cells from reference data
  }
  t <- SCTransform(object = t, assay = "RNA") # SCTransform for normalization
  if (nrow(t) < max.f) {
    max.f <- nrow(t)
  }
  seurat.list[[n]] <- t
  n <- n + 1
}

anchors <- FindIntegrationAnchors(
  object.list = seurat.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:50
)

# integrate dataset
d <- IntegrateData(
  anchorset = anchors,
  new.assay.name = "integrated",
  normalization.method = "SCT",
  dims = 1:50
)

# Batch correct with Harmony
d <- d %>%
  RunHarmony(
    group.by.vars = "tissue",
    assay.use = "SCT", reduction = "pca"
  ) %>%
  RunUMAP(
    reduction = "harmony", dims = 1:20,
    reduction.name = "harmony.umap",
    assay = "SCT", reduction.key = "hUMAP_"
  )

# Fig. 1B cartoon of experimental set up then UMAP####
f1b <- DimPlot(
  object = c0, reduction = "umap",
  pt.size = 0.1, label = T,
  label.size = 5
) +
  theme(
    plot.title = element_text(face = "plain", size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank()
  ) +
  NoLegend()
print(f1b)
# ggsave("~/Downloads/Fig1B.pdf", f1b)

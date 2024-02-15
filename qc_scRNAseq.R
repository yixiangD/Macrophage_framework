remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(ggplot2)

data.10x <- readxl::read_excel(paste(local.path, "secondary/10x.xlsx", sep = "/"))
data.10x$type <- gsub(".mtx.gz", "", data.10x$type)
data.10x$type <- gsub(".tsv.gz", "", data.10x$type)

folder <- data.10x$folder[1]
# samp <- data.10x$sample[1] # dark
samp <- data.10x$sample[1] # light
example.path <- paste(local.path, folder, sep = "/")

mtx <- paste(example.path, data.10x$file[grepl("matrix", data.10x$file) & grepl(samp, data.10x$file)], sep = "/")
cells <- paste(example.path, data.10x$file[grepl("barcode", data.10x$file) & grepl(samp, data.10x$file)], sep = "/")
features <- paste(example.path, data.10x$file[grepl("feature|gene", data.10x$file) & grepl(samp, data.10x$file)], sep = "/")
# reading
counts <- ReadMtx(mtx = mtx, cells = cells, features = features)
seurat_object <- CreateSeuratObject(counts = counts, project = samp)

# QC-------
# Filter out low-quality cells
# This will involve setting thresholds for metrics like number of detected genes, percentage of mitochondrial genes, etc.

seurat_object[["percent.mito"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

fig <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(fig)

plt1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mito")
plt2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plt1 + plt2)

# Quality control - filtering cells
# lower_bound: 200 to 500 genes
# upper_bound: 2,500 to 6,000 genes

# Cells with too few detected genes might be dying or of poor quality,
# while cells with too many detected genes could be doublets (two cells captured as one).

# mito_threshold: 5% to 20%
lower_bound <- 200
upper_bound <- 2500
mito_threshold <- 5

qc_object <- subset(seurat_object, subset = nFeature_RNA > lower_bound & nFeature_RNA < upper_bound & percent.mito < mito_threshold)

# Normalize and identify highly variable genes, if not already done
qc_object <- NormalizeData(qc_object)
qc_object <- FindVariableFeatures(qc_object)
top10 <- head(VariableFeatures(qc_object), 10)
plot1 <- VariableFeaturePlot(qc_object)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(qc_object)

# Cluster the cells
qc_object <- ScaleData(qc_object)

qc_object <- RunPCA(qc_object)
DimPlot(qc_object, reduction = "pca")
ElbowPlot(qc_object, ndims = 50)

qc_object <- FindNeighbors(qc_object)
qc_object <- FindClusters(qc_object, resolution = 0.5) # resolution can vary

# Project the data (e.g., UMAP or t-SNE)
qc_object <- RunUMAP(qc_object, dims = 1:10) # or RunTSNE
FeaturePlot(qc_object, features = c("Cd36", "C1qc"))

# cluster
qc_object <- FindNeighbors(qc_object, dims = 1:10) # louvain claster, graph based
qc_object <- FindClusters(qc_object, resolution = 0.5) # higher resolution leads to more clusters?
DimPlot(qc_object, reduction = "umap", group.by = "seurat_clusters", label = T)
FeaturePlot(qc_object, features = top10, order = T, pt.size = 0.01)

# co-expression if needed; overlay more than 2 features/genes
cells <- WhichCells(qc_object, expression = Spp1 > 1 & Rho > 1)
DimPlot(qc_object, reduction = "umap", cells.highlight = cells)

# below are just place holders need to annotate using biological info
new.cluster.id <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T")
names(new.cluster.id) <- levels(qc_object)
qc_object <- RenameIdents(qc_object, new.cluster.id)

## DE analysis----
markers <- FindAllMarkers(qc_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
CD4.mem.DEGS <- FindMarkers(qc_object, ident.1 = "Memory CD4 T", ident.2 = "Naive CD4 T", min.pct = 0.25)

## gene signature analysis----
exhaustion.genes <- list(top10)
qc_object <- AddModuleScore(qc_object, features = exhaustion.genes, name = "exhaustion.score")
FeaturePlot(qc_object, features = "exhaustion.score1", reduction = "umap")
VlnPlot(qc_object, features = "exhaustion.score1", group.by = "seurat_clusters")

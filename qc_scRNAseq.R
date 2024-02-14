remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package

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
library(ggplot2)
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

# Cluster the cells
qc_object <- ScaleData(qc_object)
qc_object <- RunPCA(qc_object)
qc_object <- FindNeighbors(qc_object)
qc_object <- FindClusters(qc_object, resolution = 0.5) # resolution can vary

# Project the data (e.g., UMAP or t-SNE)
qc_object <- RunUMAP(qc_object, dims = 1:10) # or RunTSNE

# Assuming you have a Seurat object named 'seurat_obj'
# And you want to sample 500 cells from it

# Set seed for reproducibility
set.seed(42)

# Get a vector of cell names (cell barcodes)
cell_names <- colnames(qc_object)

# Randomly sample cell names
sampled_cells <- sample(cell_names, size = 500, replace = FALSE)

# Subset the Seurat object to keep only the sampled cells
qc_obj_sampled <- subset(qc_object, cells = sampled_cells)

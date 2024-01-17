remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"
lm.file <- "GSE171328/GSM5223078_LM"
# how to create n1 objects
# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package

# Retrieve datasets
raw_data <- Read10X(data.dir = paste(local.path, lm.file, sep = "/"))

# Preprocess the data
# Filter out low-quality cells
# This will involve setting thresholds for metrics like number of detected genes, percentage of mitochondrial genes, etc.
seurat_object <- CreateSeuratObject(counts = raw_data)
mito_genes <- grep(pattern = "^MT-", rownames(seurat_object@assays$RNA@counts), value = TRUE)
seurat_object[["percent.mito"]] <- PercentageFeatureSet(seurat_object, features = mito_genes)

# Quality control - filtering cells
# lower_bound: 200 to 500 genes
# upper_bound: 2,500 to 6,000 genes

# Cells with too few detected genes might be dying or of poor quality,
# while cells with too many detected genes could be doublets (two cells captured as one).

# mito_threshold: 5% to 20%
lower_bound <- 200
upper_bound <- 2500
mito_threshold <- 0.05

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

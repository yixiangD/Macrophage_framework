# Load necessary libraries
library(Seurat)
library(clustree)
rm(list = ls())
remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"
david.data <- "data.rds"
data <- readRDS(paste(remote.path, david.data, sep = "/"))
print("reference data: ")
print(data)

# GSE157313 and LM GSE171328
# Load data (replace with your actual data loading code)
primary.file <- "GSE157313/GSM4761285_PRIMARY"
lm.file <- "GSE171328/GSM5223078_LM"
lm <- Read10X(data.dir = paste(local.path, lm.file, sep = "/"))
lm_obj <- CreateSeuratObject(counts = lm, project = "lm")
print("lm data: ")
print(lm_obj)
print(paste0("number of features: ", length(rownames(lm_obj))))
print(paste0("number of cells: ", length(colnames(lm_obj))))

mito_genes <- grep(pattern = "^MT-", rownames(lm_obj@assays$RNA@counts), value = TRUE)

# Calculate percent of mitochondrial gene expression
lm_obj[["percent.mito"]] <- PercentageFeatureSet(lm_obj, features = mito_genes)


# Quality control - filtering cells

# lower_bound: 200 to 500 genes
# upper_bound: 2,500 to 6,000 genes

# Cells with too few detected genes might be dying or of poor quality,
# while cells with too many detected genes could be doublets (two cells captured as one).

# mito_threshold: 5% to 20%
lower_bound <- 200
upper_bound <- 2500
mito_threshold <- 0.05

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > lower_bound & nFeature_RNA < upper_bound & percent.mito < mito_threshold)

# Normalization using SCTransform
seurat_obj <- SCTransform(seurat_obj)

# Integration with rpca
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50)

# Clustering resolution analysis with Clustree
# Loop over a range of resolutions
for (res in seq(0.2, 1.2, by = 0.2)) {
  # Find clusters at this resolution
  seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = 3)

  # Store the results in the metadata with a specific naming pattern
  colname <- paste0("RNA_snn_res.", res)
  seurat_obj@meta.data[[colname]] <- seurat_obj@active.ident
}

clustree(seurat_obj@meta.data, prefix = "RNA_snn_res.")

# List all genes in your Seurat object
all_genes <- rownames(seurat_obj@assays$RNA@counts)
# Check if your marker genes are in the list
markers <- c("Adgre1", "Csf1r", "H2-Ab1", "Cd68", "Lyz2", "Itgam", "Mertk")
missing_markers <- markers[!markers %in% all_genes]
for (m in markers) {
  if (m %in% all_genes) print(m)
}

raw_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
selected_genes_counts <- raw_counts[markers, ]
summary(selected_genes_counts)
barplot(colSums(selected_genes_counts))


expression_threshold <- 1e-5

subset_seurat <- subset(seurat_obj, subset = Adgre1 >= expression_threshold &
  Csf1r >= expression_threshold &
  `H2-Ab1` >= expression_threshold &
  Cd68 >= expression_threshold &
  Lyz2 >= expression_threshold &
  Itgam >= expression_threshold &
  Mertk >= expression_threshold)

num_cells <- ncol(subset_seurat[["RNA"]])
print(num_cells)
# Randomly select 500 cells per condition
rd_seurat_obj <- subset(subset_seurat, cells = sample(Cells(subset_seurat), 500))

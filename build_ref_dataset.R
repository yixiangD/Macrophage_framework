# Load necessary libraries
library(Seurat)
library(clustree)

remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"

data <- readRDS(paste(remote.path, "data.RDS", sep = "/"))
print(data)
break
# Read and preprocess data
seurat_objects <- lapply(count_matrices, function(count_matrix) {
    # Create a Seurat object
    seurat_object <- CreateSeuratObject(counts = count_matrix, project = "YourProject", min.cells = 3, min.features = 200)

    # Filter out poor quality cells
    seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 5)

    # Normalize data using SCTransform
    seurat_object <- SCTransform(seurat_object, verbose = FALSE)

    return(seurat_object)
})

# Integration
features_to_integrate <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)
seurat_objects <- lapply(seurat_objects, function(x) PrepSCTIntegration(x, anchor.features = features_to_integrate))
anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT", anchor.features = features_to_integrate)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)

# Run PCA and UMAP
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "rpca", dims = 1:50)

# Clustering
for (res in seq(0, 1.2, by = 0.1)) {
    integrated <- FindClusters(integrated, resolution = res)
    # Plot the clustree here if needed
}

# Subset for macrophages based on markers and downsample to 500 cells per condition
markers <- c("Adgre1", "Csf1r", "H2-Ab1", "Cd68", "Lyz2", "Itgam", "Mertk")
macrophages <- subset(integrated, subset = rowSums(integrated@assays$RNA@counts[markers, ] > 1) >= length(markers))

# Downsample cells - Assuming 'condition' is a metadata column
# Note: Implement a function or method for downsampling

# Split and reprocess each condition
conditions <- unique(macrophages$condition)
macrophage_datasets <- list()
for (cond in conditions) {
    subset_macrophages <- subset(macrophages, subset = condition == cond)
    # Repeat normalization, integration, PCA, and UMAP for each condition
    # Add the processed data for each condition to the macrophage_datasets list
}


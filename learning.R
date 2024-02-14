# data source: https://support.10xgenomics.com/single-cell-gene-expression/datasets
data <- Seurat::Read10X_h5("./data/SC3pv3_GEX_Human_PBMC_raw_feature_bc_matrix.h5", use.names = T)
obj <- Seurat::CreateSeuratObject(data, project = "test")
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
library(ggplot2)

fig <- Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(fig)

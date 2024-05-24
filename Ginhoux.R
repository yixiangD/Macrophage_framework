remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"
library(Seurat)
# library(Biobase)
# library(GEOquery)
# gset <- getGEO("GSE178209", GSEMatrix = TRUE, AnnotGPL = TRUE)

fname <- "~/Downloads/2021_MNP_Verse.RDS"
g.data <- readRDS(fname)

feat1 <- c("S100A8", "S100A9")
feat2 <- c("C1QA", "C1QB", "CD68")
features <- c(feat1, feat2)
# genes <- rownames(g.data)
# Seurat::DotPlot(g.data, features = features)

health.data <- subset(g.data, subset = Status == "Healthy")

# count_matrix <- health.data@assays$RNA@counts
# head(count_matrix)
# dim(count_matrix)
# count_df <- as.data.frame(as.matrix(count_matrix))

unique(health.data$Tissue)
sel.tissue <- c("Liver", "Lung")
# sel.tissue <- c("Breast","Liver","Lung","Spleen","Kidney","Pancreas")
obj.sel <- subset(health.data, subset = Tissue %in% sel.tissue)
rm(g.data)
rm(health.data)
rm(genes)
# normalize/scale data
memory.limit(size = 32000) # Increase to 16GB, adjust as necessary

obj.sel <- ScaleData(object = obj.sel, features = rownames(obj.sel))
obj.sel <-
  RunPCA(health.data, features = VariableFeatures(object = obj.sel))

obj.sel <- FindNeighbors(obj.sel, dims = 1:10)
obj.sel <- FindClusters(obj.sel)
obj.sel <- RunUMAP(obj.sel, dims = 1:10)
head(obj.sel@meta.data)

# identify differentially expressed genes----
table(obj.sel$Tissue)
table(obj.sel$Clusters)
# Set 'condition' as the active identity
Idents(obj.sel) <- obj.sel$condition

# Run differential expression analysis
de_results <-
  FindMarkers(
    obj.sel,
    ident.1 = "CHOW",
    ident.2 = "HFD",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
# View top differentially expressed genes
head(de_results)

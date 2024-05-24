rm(list = ls())
remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"
test.folder <- c("GSE171330")
# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(dplyr)
library(harmony)
library(ggplot2)
library(SingleR)
library(BiocParallel)

opt.method <- "manual"
# features/genes to plot
feat1 <- c("Csf1r", "Cd68")
feat2 <- c("Cd86", "Mertk", "Marco", "Aif1")
features <- c(feat1, feat2)
# features <- c("Adgre1", "Csf1r", "H2-Ab1", "Cd68", "Lyz2", "Itgam", "Mertk")
# Csf1r > 1 or Cd68 > 1
# All the other genes (include Cd86, mertk, csf1r (cd115), marco, aif1)

if (opt.method != "manual") {
  ref <- celldex::MouseRNAseqData()
}

for (folder in test.folder) {
  fname.integrated <-
    paste(saved.dir,
      paste(folder, "integrated.RDdata", sep = "_"),
      sep = "/"
    )
  obj.final <- readRDS(file = fname.integrated)
  obj.final <- JoinLayers(obj.final)
  obj.final <- NormalizeData(obj.final, verbose = TRUE)
  # Retrieve normalized data
  normalized_data <-
    GetAssayData(obj.final, assay = "RNA", layer = "data")

  if (opt.method == "manual") {
    # opt1. manually find macrophages----
    # Check expression of typical macrophage function genes
    DotPlot(obj.final, features = features)

    # Normalize data if not already normalized
    # Step 1----
    expression_above_threshold <- normalized_data[feat1, ] > 0
    cells_to_keep <- apply(expression_above_threshold, 2, any)

    # Subset the Seurat object to include only cells where all specified genes are above the threshold
    obj.med <- obj.final[, cells_to_keep]

    # Step 2----
    normalized_data <-
      GetAssayData(obj.med, assay = "RNA", layer = "data")
    expression_above_threshold <- normalized_data[feat2, ] > 0
    cells_to_keep <- apply(expression_above_threshold, 2, any)

    # Subset the Seurat object to include only cells where all specified genes are above the threshold
    obj.manual <- obj.med[, cells_to_keep]
    VlnPlot(obj.manual, features = features)
    DotPlot(obj.manual, features = features)
    # export to temporary data
    fname.macro <-
      paste(saved.dir,
        paste(folder, "manual_macrophage.RDdata", sep = "_"),
        sep = "/"
      )
    saveRDS(obj.manual, file = fname.macro)
  } else {
    # opt2. automate finding macrophages----

    sceObject <- as.SingleCellExperiment(obj.final)
    singleRResults <-
      SingleR(
        test = sceObject,
        ref = ref,
        labels = ref$label.main,
        clusters = obj.final$seurat_clusters
      )
    # Viewing the first few annotations
    head(singleRResults$labels)
    # Summary of predicted cell types
    table(singleRResults$labels)

    scores <- as.data.frame(singleRResults@listData$scores)
    # plot
    score.heatmap <-
      plotScoreHeatmap(singleRResults,
        # clusters = singleRResults$labels,
        fontsize.row = 9,
        show_colnames = T
      )

    # extract macrophages only
    # Subset for a specific cluster, e.g., cluster 0
    head(obj.final@meta.data)

    # Subset for a specific annotated cell type
    obj.final$SingleR.labels <- singleRResults$labels
    obj.auto <-
      subset(obj.final, subset = SingleR.labels == "Monocytes")

    # export to temporary data
    fname.macro <-
      paste(saved.dir,
        paste(folder, "auto_macrophage.RDdata", sep = "_"),
        sep = "/"
      )
    saveRDS(obj.auto, file = fname.macro)
  }
}

rm(list = ls())
remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"
# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(dplyr)
library(harmony)
library(ggplot2)

library(celldex)
library(SingleR)
library(BiocParallel)
# find markers----
# manual
# markers <- FindConservedMarkers(obj.final, ident.1 = 6, grouping.var = "condition", verbose = FALSE)
# macropage features
test.folder <- c("GSE171330")
for (folder in test.folder) {
  fname.integrated <- paste(saved.dir, paste(folder, "integrated.RDdata", sep = "_"), sep = "/")
  obj.final <- readRDS(file = fname.integrated)
  features <-
    c("Adgre1", "Csf1r", "H2-Ab1", "Cd68", "Lyz2", "Itgam", "Mertk")
  feat.plot <- FeaturePlot(obj.final, features = features)
  
  ggsave(
    paste(saved.dir, paste0(
      paste(folder, "feat.plot", sep = "_"), ".pdf"
    ), sep = "/"),
    feat.plot,
    width = 10,
    height = 5
  )
  
  obj.small <-
    subset(
      obj.final,
      Adgre1 > 1 |
        Csf1r > 1 |
        `H2-Ab1` > 1 | Cd68 > 1 | Lyz2 > 1 | Itgam > 1 | Mertk > 1
    )
  obj.small <- RunUMAP(obj.small, reduction = "harmony", dims = 1:30)
  obj.small <-
    FindNeighbors(obj.small, reduction = "harmony", dims = 1:30) %>% FindClusters()
  
  dim.plot.small <-
    DimPlot(obj.small,
            group.by = c("condition", "ident"),
            ncol = 2
    )
  ggsave(
    paste(saved.dir, paste0(
      paste(folder, "dim.plot.small", sep = "_"), ".pdf"
    ), sep = "/"),
    dim.plot.small,
    width = 10,
    height = 5
  )
  
  ref <- celldex::MouseRNAseqData()
  obj.final <- JoinLayers(obj.final)
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
  
  score.heatmap <-
    plotScoreHeatmap(
      singleRResults,
      clusters = singleRResults$labels,
      fontsize.row = 9,
      show_colnames = T
    )
  ggsave(
    paste(saved.dir, paste0(
      paste(folder, "score.heatmap", sep = "_"), ".pdf"
    ), sep = "/"),
    score.heatmap,
    width = 8,
    height = 6
  )
  
  new.cluster.ids <- singleRResults$pruned.labels
  names(new.cluster.ids) <- levels(obj.final)
  levels(obj.final)
  obj.final <- RenameIdents(obj.final, new.cluster.ids)
  
  umap.integrated <-
    UMAPPlot(
      object = obj.final,
      pt.size = 0.5,
      label = TRUE
    )
  print(umap.integrated)
  
  ggsave(
    paste(saved.dir, paste0(
      paste(folder, "clusterUMAP", sep = "_"), ".pdf"
    ), sep = "/"),
    umap.integrated,
    width = 6,
    height = 5
  )
}

library(harmony)
library(Seurat)
library(SeuratData)
# InstallData("pbmcsca")
# data("pbmcsca")
#
# pbmcsca <- UpdateSeuratObject(pbmcsca)
#
# pbmcsca <- NormalizeData(pbmcsca) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
# pbmcsca <- RunHarmony(pbmcsca, group.by.vars = "Method")
# pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:30)
# pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:30) %>% FindClusters()
# DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)

InstallData("ifnb")
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
ifnb <- NormalizeData(ifnb) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)
ifnb <- RunHarmony(ifnb, group.by.vars = "stim")
ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:30)
ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:30) %>% FindClusters()
DimPlot(ifnb, group.by = c("stim", "ident", "seurat_annotations"), ncol = 3)

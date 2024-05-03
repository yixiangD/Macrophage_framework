rm(list = ls())
remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"
test.folder <- c("GSE171330")
opt.method <- "manual"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(dplyr)
library(harmony)
library(ggplot2)
library(SingleR)
library(slingshot)
library(SingleCellExperiment)
library(BiocParallel)

for (folder in test.folder) {
  fname <-
    paste(saved.dir,
      paste(folder, opt.method, "macrophage.RDdata", sep = "_"),
      sep = "/"
    )
  obj <- readRDS(file = fname)
  # normalize/scale data
  obj.sel <- ScaleData(object = obj, features = rownames(obj))
  obj.sel <-
    RunPCA(obj.sel, features = VariableFeatures(object = obj.sel))
  obj.sel <- FindNeighbors(obj.sel, dims = 1:10)
  obj.sel <- FindClusters(obj.sel)
  obj.sel <- RunUMAP(obj.sel, dims = 1:10)
  head(obj.sel@meta.data)

  # identify differentially expressed genes----
  table(obj.sel$condition)

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
  # Dot plot for multiple DE genes
  DotPlot(obj.sel, features = rownames(de_results)[1:10])

  # Set up and run a pseudotime analysis using Slingshot----
  sce <- as.SingleCellExperiment(obj.sel)

  sce <-
    slingshot(sce,
      clusterLabels = colData(sce)$seurat_clusters,
      reducedDim = "UMAP"
    )

  slingshot_data <- metadata(sce)$SlingshotDataSet

  # Extract pseudotime and store it
  colData(sce)$slingshot_pseudotime <-
    slingPseudotime(SlingshotDataSet(sce), na = FALSE)
  colData(sce)$slingshot_lineages <-
    slingLineages(SlingshotDataSet(sce))

  # Assuming slingshot_lineages is incorrectly stored as a list of lists or similar
  # Convert it to a usable format, such as a factor or numeric vector

  lineage_vector <-
    sapply(colData(sce)$slingshot_lineages, function(x) {
      x[[1]]
    })
  colData(sce)$slingshot_lineages <- as.factor(lineage_vector)

  # Assign colors based on the number of unique lineages
  num_lineages <- length(unique(colData(sce)$slingshot_lineages))
  lineage_colors <- rainbow(num_lineages)

  # Plot using the corrected lineages and assigned colors
  plot(
    reducedDims(sce)$UMAP[, 1],
    reducedDims(sce)$UMAP[, 2],
    col = lineage_colors[as.integer(colData(sce)$slingshot_lineages)],
    pch = 16,
    xlab = "UMAP 1",
    ylab = "UMAP 2",
    main = "UMAP with Slingshot Trajectories"
  )
  lines(SlingshotDataSet(sce),
    lwd = 2,
    col = "black"
  ) # Add trajectory lines

  # Extract pseudotime values
  pseudotime <- slingPseudotime(sce, na = FALSE)

  # Examine gene expression along pseudotime
  gene_of_interest <-
    rownames(de_results)[1:10] # Change to your gene of interest
  plot(
    pseudotime,
    assay(sce)[gene_of_interest[1], ],
    # fix dimension
    pch = 16,
    cex = 0.6,
    xlab = "Pseudotime",
    ylab = "Expression level of gene"
  )
}

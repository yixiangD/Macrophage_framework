library(Seurat)
memory.limit(size = 8000)
fname <- "./data/2021_MoMac_VERSE.RDS"
print("Reading...")
g.data <- readRDS(fname)
print("End Reading...")

print("Inspect whole data...")
print(unique(g.data$Tissue))
feat1 <- c("S100A8", "S100A9")
feat2 <- c("C1QA", "C1QB", "CD68")
features <- c(feat1, feat2)
# genes <- rownames(g.data)
# Seurat::DotPlot(g.data, features = features)

print("Inspect healthy data...")
health.data <- subset(g.data, subset = Status == "Healthy")

# count_matrix <- health.data@assays$RNA@counts
# head(count_matrix)
# dim(count_matrix)
# count_df <- as.data.frame(as.matrix(count_matrix))

print(unique(health.data$Tissue))
sel.tissue <- c("Lung", "Kidney", "Skin", "Breast")
# sel.tissue <- c("Breast","Liver","Lung","Spleen","Kidney","Pancreas")
obj.tissue <- subset(health.data, subset = Tissue %in% sel.tissue)
print(unique(obj.tissue$Tissue))
# rm(g.data)
# rm(health.data)
# normalize/scale data

cells.keep <- grepl("Macrophage", obj.tissue$Clusters) & grepl("-2|-3|-6|-16", obj.tissue$Clusters)
obj.sel <- subset(obj.tissue, cells = Cells(obj.tissue)[cells.keep])
print(unique(obj.sel$Clusters))
exit()
obj.sel <- Seurat::ScaleData(object = obj.sel, features = rownames(obj.sel))
obj.sel <-
  Seurat::RunPCA(obj.sel, features = Seurat::VariableFeatures(object = obj.sel))

obj.sel <- Seurat::FindNeighbors(obj.sel, dims = 1:10)
obj.sel <- Seurat::FindClusters(obj.sel)
obj.sel <- Seurat::RunUMAP(obj.sel, dims = 1:10)
head(obj.sel@meta.data)

# identify differentially expressed genes----
table(obj.sel$Tissue)
table(obj.sel$Clusters)
# Set 'condition' as the active identity
Idents(obj.sel) <- obj.sel$Tissue

# Run differential expression analysis
id1 <- "Lung"
id2 <- "Liver"
for (tis in sel.tissue) {
  if (tis != "Lung") {
    de_results <-
      Seurat::FindMarkers(
        obj.sel,
        ident.1 = "Lung",
        ident.2 = id2,
        min.pct = 0.1,
        logfc.threshold = 0.25,
        test.use = "wilcox"
      )
    # View top differentially expressed genes
    head(de_results)

    df <- de_results
    df$gene <- rownames(df)
    df$log10_p_val <- -log10(df$p_val)

    # Create the volcano plot
    library(ggplot2)
    log2FC_threshold <- 1
    p_val_threshold <- 0.05

    # Add a column to indicate significance
    df$significant <- with(df, log10_p_val > -log10(p_val_threshold) & abs(avg_log2FC) > log2FC_threshold)

    # Create the volcano plot with highlights
    fig <- ggplot(df, aes(x = avg_log2FC, y = log10_p_val)) +
      geom_point(aes(color = significant), size = 2) +
      scale_color_manual(values = c("grey", "red")) +
      geom_text(aes(label = ifelse(significant, as.character(gene), "")), hjust = 0.5, vjust = -0.5, size = 3) +
      geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "blue") +
      theme_minimal() +
      ggtitle("Volcano Plot of Differential Gene Expression") +
      xlab("log2 Fold Change") +
      ylab("-log10 p-value")

    ggsave(paste0("./output/", id1, "_", id2, ".pdf"), fig)
  }
}

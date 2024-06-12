remote.path <-
  "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <-
  "/Users/yd973/Documents/research/Macrophage_framework/data"
saved.dir <- "~/Downloads/macrophage"
library(Seurat)
# library(Biobase)
# library(GEOquery)
# gset <- getGEO("GSE178209", GSEMatrix = TRUE, AnnotGPL = TRUE)
memory.limit(size=8000)
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
sel.tissue <- c("Liver", "Lung", "Kidney", "Skin", "Breast")
# sel.tissue <- c("Breast","Liver","Lung","Spleen","Kidney","Pancreas")
obj.sel <- subset(health.data, subset = Tissue %in% sel.tissue)
print(unique(obj.sel$Tissue))
#rm(g.data)
#rm(health.data)
# normalize/scale data

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
de_results <-
  Seurat::FindMarkers(
    obj.sel,
    ident.1 = "Lung",
    ident.2 = "Liver",
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

print("plotting")
# Create the volcano plot with highlights
fig <- ggplot(df, aes(x = avg_log2FC, y = log10_p_val)) +
  geom_point(aes(color = significant), size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text(aes(label = ifelse(significant, as.character(gene), '')), hjust = 0.5, vjust = -0.5, size = 3) +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "blue") +
  theme_minimal() +
  ggtitle("Volcano Plot of Differential Gene Expression") +
  xlab("log2 Fold Change") +
  ylab("-log10 p-value")

ggsave("./output/lung_liver.pdf", fig)

# read sc RNA seq data

sc.dir <- list.dirs(path = "./data", full.names = TRUE, recursive = FALSE)

# scRNA_seq_data now contains your data frame

for (d in sc.dir) {
  files <- list.files(d)
  # print(d)
  grps <- c()
  for (f in files) {
    example <- f[grepl("barcode", f)]
    grp <- gsub("_barcodes.tsv.gz", "", example)
    grps <- c(grps, grp)
  }
  # different platform: Illumina Nextseq 500
  if (length(grps) == 0) {
    print(d)
  }
}
# Path to your .txt.gz file
file_path <- "/Users/yd973/Documents/research/Macrophage_framework/data/GSE144707/GSM4294086_countTable_nerveStSt.txt.gz"

# Reading the gzipped file directly
scRNA_seq_data <- read.table(gzfile(file_path), header = TRUE, sep = "\t")
seurat_obj <- Seurat::CreateSeuratObject(counts = scRNA_seq_data, project = "test")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)

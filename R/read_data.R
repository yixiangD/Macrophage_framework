# read sc RNA seq data

sc.dir <- list.dirs(path = "./data", full.names = TRUE, recursive = FALSE)
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

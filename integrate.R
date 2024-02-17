rm(list = ls())
remote.path <- "/Users/yd973/Dropbox (Partners HealthCare)/macrophage/data"
local.path <- "/Users/yd973/Documents/research/Macrophage_framework/data"

# Load necessary libraries
library(Seurat) # or another scRNA-seq analysis package
library(dplyr)

data.10x <- readxl::read_excel(paste(local.path, "secondary/10x.xlsx", sep = "/"))
data.10x$type <- gsub(".mtx.gz", "", data.10x$type)
data.10x$type <- gsub(".tsv.gz", "", data.10x$type)

df <- data.10x %>%
  mutate(
    splits = strsplit(sample, "_")
  ) %>%
  rowwise() %>%
  mutate(
    code = splits[1],
    case = paste(splits[2:length(splits)], collapse = "_")
  ) %>%
  select(-splits)

df <- df %>%
  group_by(folder, type) %>%
  mutate(id = row_number())

df <- df[order(df$folder, df$id), ]

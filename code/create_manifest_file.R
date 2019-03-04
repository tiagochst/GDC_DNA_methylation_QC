# Export manifest files
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readr)
setwd("/mnt/home/tiagochst/paper_elmer/GDAN/")
files <- dir(pattern = "correlation_final.csv",recursive = T,full.names = T)
tumors <- unique(basename(dirname(dirname(files))))

all.samples.hg38 <- c()
for(tumor in tumors){
  print(tumor)
  # This will keep samples with the both DNA methylation and gene expression
  # Take the log2(exp + 1)
  # Put samples in the right order
  # Keep only probes in 2KB close to TSS
  # Remove probes that have NA for more than 50% of samples
  mae.file <- paste0(paste0(tumor),"/",tumor,"_mae_hg38.rda")
  if(file.exists(mae.file)) {
    mae <- get(load(mae.file))
    all.samples.hg38 <- c(all.samples.hg38,mae$sample)
  } 
}

all.samples.hg19 <- c()
for(tumor in tumors){
  print(tumor)
  mae.file <- paste0(paste0(tumor,"/hg19/"),"/",tumor,"_mae_hg19.rda")
  if(file.exists(mae.file)) {
    mae <- get(load(mae.file))
    all.samples.hg19 <- c(all.samples.hg19,mae$sample)
  } 
}
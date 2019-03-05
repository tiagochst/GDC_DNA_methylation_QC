library(ELMER)

#
# Plot DNA methylation vs Gene expression
#
ntop <- 100 # plot the top correlations

tumors <- TCGAbiolinks:::getGDCprojects()$project_id
tumors <- tumors[grepl('^TCGA', tumors, perl = TRUE)]

# hg38
files <- dir(pattern = "correlation.csv", recursive = T, full.names = T)
for(tumor in tumors){
  mae.file <- grep(tumor,dir(pattern = "mae_hg38", recursive = T, full.names = T), value = TRUE)
  if(length(mae.file) == 0) next
  load(mae.file)
  f <- sort(grep(tumor, files,value = T))
  hg38_correlation <- read_csv(f[2])
  colnames(hg38_correlation)[grep("ENSG", hg38_correlation[1,])] <- "ensembl_gene_id"
  for(i in 1:ntop){
    library(ELMER)
    dir.create(paste0(tumor,"/hg38/plots/"), recursive = T)
    scatter.plot(data = mae,
                 dir.out = paste0(tumor,"/hg38/plots/"),
                 save = TRUE,
                 byPair = list(probe = as.character(hg38_correlation$probe[i]),
                               gene = as.character(hg38_correlation$ensembl_gene_id[i])),
                 category = "definition")
  }
}


# hg19
files <- dir(pattern = "correlation.csv", recursive = TRUE, full.names = TRUE)
for(tumor in tumors){
  f <- sort(grep(tumor, files,value = T))
  hg19_correlation <- read_csv(f[1])
  colnames(hg19_correlation)[grep("cg",hg19_correlation[1,])] <- "probe"
  mae.file <- grep(tumor,dir(pattern = "mae_hg19", recursive = TRUE, full.names = TRUE), value = TRUE)
  if(length(mae.file) == 0) next
  load(mae.file)
  for(i in 1:ntop){
    dir.create(paste0(tumor, "/hg19/plots/"), recursive = T)
    scatter.plot(data = mae,
                 dir.out = paste0(tumor, "/hg19/plots/"),
                 save = TRUE,
                 byPair = list(probe = as.character(hg19_correlation$probe[i]),
                               gene = as.character(hg19_correlation$ensembl_gene_id[i])),
                 category = "definition")
  }
}

library(readr)
library(dplyr)

load("/mnt/home/tiagochst/paper_elmer/GDAN/associated.genes_id_hg19.rda")
load("/mnt/home/tiagochst/paper_elmer/GDAN/associated.genes_id.rda")
colnames(associated.genes_id.hg19)[2] <- "probe"
associated.genes_id$distance <- NULL
associated.genes_id <- unique(associated.genes_id)

files <- dir(pattern = "correlation.csv",recursive = T,full.names = T)
tumors <- basename(dirname(dirname(files)))
snc.in.hg38.not.in.hg19 <- NULL
for(tumor in tumors){
  print(tumor)
  f <- sort(grep(tumor, files,value = T))
  hg19_correlation <- read_csv(f[1])
  colnames(hg19_correlation)[grep("cg",hg19_correlation[1,])] <- "probe"
  hg38_correlation <- read_csv(f[2])
  colnames(hg38_correlation)[grep("ENSG",hg38_correlation[1,])] <- "ensembl_gene_id"
  hg19_correlation <- left_join(as_tibble(hg19_correlation),unique(as_tibble(associated.genes_id.hg19)))
  hg38_correlation <- left_join(as_tibble(hg38_correlation),unique(as_tibble(associated.genes_id)))
  save(hg19_correlation,hg38_correlation, file = paste0(tumor,"_correlation_final.rda"))
  
  if(any(hg19_correlation$status == "Strongly negatively correlated (SNC)")){
    y <- hg19_correlation[hg19_correlation$status == "Strongly negatively correlated (SNC)",]
    z <- hg38_correlation[hg38_correlation$status == "Strongly negatively correlated (SNC)",]
    z <- unique(z[!z$ensembl_gene_id %in% y$ensembl_gene_id,c(8:9)])
    z$tumor <- tumor
    snc.in.hg38.not.in.hg19 <- rbind(snc.in.hg38.not.in.hg19,z)
  }
}
readr::write_csv(unique(as.data.frame(snc.in.hg38.not.in.hg19)),path = "snc_in_hg38_not_in_hg19.csv")

# Evaluating  only negative
library(readr)
load("/mnt/home/tiagochst/paper_elmer/GDAN/associated.genes_id_hg19.rda")
load("/mnt/home/tiagochst/paper_elmer/GDAN/associated.genes_id.rda")
colnames(associated.genes_id.hg19)[2] <- "probe"
associated.genes_id <- associated.genes_id[as.numeric(as.character(associated.genes_id$distance)) < 0,]
associated.genes_id$distance <- NULL
associated.genes_id <- unique(associated.genes_id)

files <- dir(pattern = "correlation.csv",recursive = T,full.names = T)
tumors <- basename(dirname(dirname(files)))
snc.in.hg38.not.in.hg19 <- NULL
wnc.in.hg38.not.in.hg19 <- NULL
wnc.in.hg38.not.in.hg19.all <- NULL
snc.in.hg38.not.in.hg19.all <- NULL
for(tumor in tumors){
  print(tumor)
  f <- sort(grep(tumor, files,value = T))
  hg19_correlation <- read_csv(f[1])
  colnames(hg19_correlation)[grep("cg",hg19_correlation[1,])] <- "probe"
  hg38_correlation <- read_csv(f[2])
  colnames(hg38_correlation)[grep("ENSG",hg38_correlation[1,])] <- "ensembl_gene_id"
  hg19_correlation <- left_join(as_tibble(hg19_correlation),unique(as_tibble(associated.genes_id.hg19)))
  hg38_correlation <- left_join(as_tibble(hg38_correlation),unique(as_tibble(associated.genes_id)))
  save(hg19_correlation,hg38_correlation, file = paste0(tumor,"_correlation_only_negative.rda"))
  
  if(any(hg19_correlation$status == "Strongly negatively correlated (SNC)")){
    y <- hg19_correlation[hg19_correlation$status ==  "Strongly negatively correlated (SNC)",]
    z <- hg38_correlation[hg38_correlation$status ==  "Strongly negatively correlated (SNC)",]
    z$tumor <- tumor
    aux <- unique(z[!z$ensembl_gene_id %in% y$ensembl_gene_id,c(8:9)])
    snc.in.hg38.not.in.hg19 <- rbind(snc.in.hg38.not.in.hg19,aux)
    aux <- unique(z[!z$ensembl_gene_id %in% y$ensembl_gene_id,])
    snc.in.hg38.not.in.hg19.all <- rbind(snc.in.hg38.not.in.hg19.all,aux)
  }
  if(any(hg19_correlation$status == "Weakly negatively correlated (WNC)")){
    y <- hg19_correlation[hg19_correlation$status == "Weakly negatively correlated (WNC)",]
    z <- hg38_correlation[hg38_correlation$status == "Weakly negatively correlated (WNC)",]
    z$tumor <- tumor
    aux <- unique(z[!z$ensembl_gene_id %in% y$ensembl_gene_id,c(8:9)])
    wnc.in.hg38.not.in.hg19 <- rbind(wnc.in.hg38.not.in.hg19,aux)
    aux <- unique(z[!z$ensembl_gene_id %in% y$ensembl_gene_id,])
    wnc.in.hg38.not.in.hg19.all <- rbind(wnc.in.hg38.not.in.hg19.all,aux)
  }
}
readr::write_csv(unique(as.data.frame(snc.in.hg38.not.in.hg19.all)),
                 path = "snc_in_hg38_negative_1500bp_not_in_hg19_negative_1500bp_all_info.csv")
readr::write_csv(unique(as.data.frame(wnc.in.hg38.not.in.hg19.all)),
                 path = "wnc_in_hg38_negative_1500bp_not_in_hg19_negative_1500bp_all_info.csv")

readr::write_csv(unique(as.data.frame(snc.in.hg38.not.in.hg19)),
                 path = "snc_in_hg38_negative_1500bp_not_in_hg19_negative_1500bp.csv")
readr::write_csv(unique(as.data.frame(wnc.in.hg38.not.in.hg19)),
                 path = "wnc_in_hg38_negative_1500bp_not_in_hg19_negative_1500bp.csv")

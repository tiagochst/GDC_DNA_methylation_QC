library(doParallel)
library(parallel)
library(dplyr)
library(ELMER)
library(SummarizedExperiment)
library(TCGAbiolinks)
root <- "/mnt/home/tiagochst/paper_elmer/GDAN/"
setwd(root)

# Illumina manifest
# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534
metadata.hg19 <- readr::read_tsv(file = paste0(root,"/450k-manifest-GPL13534-11288.txt"),skip = 37)
metadata.hg19 <- metadata.hg19[!is.na(metadata.hg19$UCSC_RefGene_Name),]
metadata.hg19 <- metadata.hg19[grep("TSS1500|TSS200",metadata.hg19$UCSC_RefGene_Group),]

file <- "associated.genes_id_hg19.rda"
associated.genes_id.hg19 <- NULL
if(file.exists(file)) {
  associated.genes_id.hg19 <- get(load(file))
  colnames(associated.genes_id.hg19)[3] <- "probe"
}


if(is.null(associated.genes_id.hg19)){
  associated.genes_id.hg19 <- plyr::adply(.data = metadata.hg19,
                                          .margins = 1,
                                          .fun = function(x){
                                            genes <- x["UCSC_RefGene_Name"] %>% stringr::str_split(";") %>% unlist
                                            groups <- x["UCSC_RefGene_Group"] %>% stringr::str_split(";") %>% unlist
                                            idx <- grep("TSS",groups)
                                            data.frame(gene = genes[idx],group_type = groups[idx])
                                          },.id = NULL,.expand = T,.progress = T,.parallel = TRUE)

  # get all LOC gene names and update them there are not in ENSEMBL database with LOC
  symbols <- gsub("LOC","",unique(grep("^LOC",associated.genes_id.hg19$gene,value = T)))
  mapping <- select(org.Hs.eg.db, symbols, c("ENTREZID","GENENAME","ENSEMBL","SYMBOL"))

  # Replace LOC by gene names. LOC should be removed as this gene names are LOC + entrez ID
  aux <- associated.genes_id.hg19
  idx <- grep("^LOC",aux$gene)
  aux$gene[idx] <- as.character(mapping$SYMBOL[match(gsub("LOC","",aux$gene[idx]),mapping$ENTREZID)])
  idx <- !aux$gene %in% tss$external_gene_name
  aux$gene[idx] <- HGNChelper::checkGeneSymbols(aux$gene[idx])$Suggested.Symbol
  associated.genes_id.hg19 <- aux

  associated.genes_id.hg19 <- merge(associated.genes_id.hg19,
                                    values(rna)[,1:2],
                                    by.x = "gene",
                                    by.y = "external_gene_name")

  save(associated.genes_id.hg19,file = file)
}

library(doParallel)
library(parallel)
library(dplyr)
library(ELMER)
cores <- 2
registerDoParallel(cores)


tumors <- TCGAbiolinks:::getGDCprojects()$project_id
tumors <- tumors[grepl('^TCGA', tumors, perl = TRUE)]

for(tumor in tumors){
  print(tumor)
  if(file.exists(paste0(paste0(tumor,"/hg19/"),"/",tumor,"_hg19_correlation.csv"))) next
  dir.create(paste0(tumor,"/hg19/"),showWarnings = FALSE,recursive = TRUE)
  met <- get(load(paste0(root,"/Data/",tumor,"/",tumor,"_meth_hg38_no_filter.rda")))
  met <- met[unique(associated.genes_id.hg19$ID),]
  rna <- get(load(paste0(root,"/Data/",tumor,"/",tumor,"_RNA_hg38.rda")))

  # This will keep samples with the both DNA methylation and gene expression
  # Take the log2(exp + 1)
  # Put samples in the right order
  # Keep only probes in 2KB close to TSS
  # Remove probes that have NA for more than 50% of samples
  dir.create(paste0(paste0(tumor,"/hg19/")),showWarnings = F,recursive = T)
  mae.file <- paste0(paste0(tumor,"/hg19/"),"/",tumor,"_mae_hg19.rda")
  if(file.exists(mae.file)) {
    mae <- get(load(mae.file))
  } else {
    mae <- createMAE(exp = rna,
                     met = met,
                     TCGA = TRUE,
                     linearize.exp = T,
                     save = T,
                     save.filename = mae.file,
                     met.platform = "450K",
                     genome = "hg19",
                     met.na.cut = 0.5)
    save(mae, file = mae.file)

  }
  correlations <- plyr::adply(associated.genes_id.hg19,
                              .margins = 1,
                              .fun =  function(x){
                                GeneID <- as.character(x$ensembl_gene_id)
                                probe <- as.character(x$probe)
                                if(!probe %in% names(getMet(mae))) {
                                  return(tibble::tibble(rho = NA, pval = NA, ensembl_gene_id = GeneID, probe))
                                }
                                cor <- cor.test(x = as.numeric(assay(getMet(mae)[probe])),
                                                y = as.numeric(assay(getExp(mae)[GeneID])),
                                                method = c("spearman"))
                                corval <- as.numeric(cor$estimate)
                                pvalue <- as.numeric(cor$p.value)
                                tibble::tibble(rho = corval,pval = pvalue, ensembl_gene_id = GeneID, probe)
                              },.id = NULL,
                              .expand = T,
                              .progress = "text",
                              .parallel = T)
  #correlations <- correlations[!is.na(correlations$rho),] # Remove associations not found
  correlations$rho <- as.numeric(as.character(correlations$rho))
  correlations$pval <- as.numeric(as.character(correlations$pval))
  correlations$FDR <- p.adjust(as.numeric(as.character(correlations$pval)),method = "BH")

  # If FDR > 0.05 - Insignificant
  # If FDR <= 0.05:
  # - Strongly negatively correlated (SNC) when the rho value was less than -0.5;
  # - Weakly negatively correlated (WNC) when the rho value was between -0.5 and -0.25;
  # - No negative correlation (NNC) when the rho value was greater than -0.25.
  correlations <- correlations[!is.na( as.numeric(correlations$rho)),]
  correlations$status <- "Insignificant"
  correlations[correlations$FDR <= 0.05 & as.numeric(correlations$rho) >= -0.25,]$status <- "No negative correlation"
  correlations[correlations$FDR <= 0.05 & correlations$rho < -0.25 & correlations$rho > -0.5,]$status <- "Weakly negatively correlated (WNC)"
  correlations[correlations$FDR <= 0.05 & correlations$rho <= -0.5,]$status <- "Strongly negatively correlated (SNC)"
  correlations <- correlations[with(correlations,order(rho,pval)),]

  # save
  readr::write_csv(unique(as.data.frame(correlations)),
                   path = paste0(paste0(tumor,"/hg19/"),"/",tumor,"_hg19_correlation.csv"))
}

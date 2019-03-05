library(doParallel)
library(parallel)
library(dplyr)
library(ELMER)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readr)
root <- "/mnt/home/tiagochst/paper_elmer/GDAN/"
setwd(root)

# Parallelizing code
cores <- 2
registerDoParallel(cores)

file <- "associated.genes_id.rda"
associated.genes_id <- NULL
if(file.exists(file)) {
  associated.genes_id <- unique(get(load(file)))
  associated.genes_id$distance <- NULL
  associated.genes_id <- unique(associated.genes_id)
}

tumors <- TCGAbiolinks:::getGDCprojects()$project_id
tumors <- tumors[grepl('^TCGA', tumors, perl = TRUE)]
for(tumor in tumors){
  print(tumor)
  if(file.exists(paste0(paste0(tumor,"/hg38/"),tumor,"_hg38_correlation.csv"))) next
  dir.create(paste0(tumor,"/hg38/"),showWarnings = FALSE,recursive = T)
  met <- get(load(paste0(root,"/Data/",tumor,"/",tumor,"_meth_hg38_no_filter.rda")))
  rna <- get(load(paste0(root,"/Data/",tumor,"/",tumor,"_RNA_hg38.rda")))
  gene.metadata <- values(rna)[,1:2]
  met <- met[rowRanges(met)$Gene_Symbol != ".",]

  if(is.null(associated.genes_id)){
    associated.genes_id <- plyr::adply(.data = values(met),
                                       .margins = 1,
                                       .fun = function(x){
                                         genes <- x["Gene_Symbol"] %>% stringr::str_split(";") %>% unlist
                                         groups <- x["Gene_Type"] %>% stringr::str_split(";") %>% unlist
                                         distance <- x["Position_to_TSS"] %>% stringr::str_split(";") %>% unlist
                                         idx <- which(abs(as.numeric(distance)) <= 1500)
                                         data.frame(gene = genes[idx],
                                                    group_type = groups[idx],
                                                    distance = distance[idx],
                                                    probe = rep(x$Composite.Element.REF,
                                                                length(idx)))
                                       },.id = NULL,.expand = T,.progress = T,.parallel = T)
    save(associated.genes_id, file = "associated.genes_id_not_merged_max_dist_1500.rda")

    load("Human_genes__GRCh38_p12__tss.rda")
    idx <- !associated.genes_id$gene %in% tss$external_gene_name
    associated.genes_id$gene[idx] <- as.character(HGNChelper::checkGeneSymbols(associated.genes_id$gene[idx])$Suggested.Symbol)

    associated.genes_id <- unique(
      merge(associated.genes_id,
            unique(tss[,c("external_gene_name","ensembl_gene_id")]),
            by.x = "gene",
            by.y = "external_gene_name")
    )
    save(associated.genes_id,file = file)
  }

  # This will keep samples with the both DNA methylation and gene expression
  # Take the log2(exp + 1)
  # Put samples in the right order
  # Keep only probes in 2KB close to TSS
  # Remove probes that have NA for more than 50% of samples
  mae.file <- paste0(paste0(tumor),"/",tumor,"_mae_hg38.rda")
  if(file.exists(mae.file)) {
    mae <- get(load(mae.file))
  } else {
    mae <- createMAE(exp = rna,
                     met = met[rownames(met) %in% associated.genes_id$probe,],
                     TCGA = TRUE,
                     linearize.exp = T,
                     save = T,
                     save.filename = paste0(tumor,"/",tumor,"_mae_hg38.rda"),
                     met.platform = "450K",
                     genome = "hg38",
                     met.na.cut = 0.5)
  }
  correlations <- plyr::adply(unique(associated.genes_id[,3:4]),
                              .margins = 1,
                              .fun =  function(x){
                                GeneID <- as.character(x$ensembl_gene_id)
                                probe <- as.character(x$probe)
                                if(!probe %in% names(getMet(mae))) {
                                  return(tibble::tibble(rho = NA,
                                                        pval = NA,
                                                        probe,
                                                        ensembl_gene_id =  GeneID))
                                }
                                cor <- cor.test(x = as.numeric(assay(getMet(mae)[probe])),
                                                y = as.numeric(assay(getExp(mae)[GeneID])),
                                                method = c("spearman"))
                                corval <- as.numeric(cor$estimate)
                                pvalue <- as.numeric(cor$p.value)

                                df <- tibble::tibble(rho = corval,
                                                     pval = pvalue,
                                                     probe,
                                                     ensembl_gene_id = GeneID)
                                return(df)
                              },.id = NULL,.expand = T,.progress = "text",.parallel = F)
  #correlations <- correlations[!is.na(correlations$rho),] # Remove associations not found
  correlations$rho <- as.numeric(as.character(correlations$rho))
  correlations$pval <- as.numeric(as.character(correlations$pval))
  correlations$FDR <- p.adjust(as.numeric(as.character(correlations$pval)),method = "BH")

  # If FDR > 0.05 - Insignificant
  # If FDR <= 0.05:
  # - Strongly negatively correlated (SNC) when the rho value was less than -0.5;
  # - Weakly negatively correlated (WNC) when the rho value was between -0.5 and -0.25;
  # - No negative correlation (NNC) when the rho value was greater than -0.25.

  correlations$status <- "Insignificant"
  correlations$status[correlations$FDR <= 0.05 & as.numeric(correlations$rho) >= -0.25] <- "No negative correlation"
  correlations$status[correlations$FDR <= 0.05 & correlations$rho < -0.25 & correlations$rho > -0.5] <- "Weakly negatively correlated (WNC)"
  correlations$status[correlations$FDR <= 0.05 & correlations$rho <= -0.5] <- "Strongly negatively correlated (SNC)"
  correlations <- correlations[with(correlations,order(rho,pval)),]

  # save
  readr::write_csv(as.data.frame(correlations),
                   path = paste0(paste0(tumor,"/hg38/"),"/",tumor,"_hg38_correlation.csv"))

  aux <- left_join(as_tibble(correlations),
                   unique(as_tibble(associated.genes_id)))
  readr::write_csv(aux,
                   path = paste0(paste0(tumor,"/hg38/"),"/",tumor,"_hg38_correlation_all_info.csv"))

}

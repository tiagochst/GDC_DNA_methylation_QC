library(ELMER)
library(TCGAbiolinks)
root <- "/mnt/home/tiagochst/paper_elmer/GDAN/"
setwd(root)

# get TCGA projects
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]

# Download RNA-seq and DNA methylation for hg19 and hg38
for(proj in gsub("TCGA-", "", projects)){
  for(genome in c("hg38","hg19")){
    getTCGA(disease = proj, # TCGA disease abbreviation (BRCA, BLCA, GBM, LGG, etc)
            basedir = "Data", # Where data will be downloaded
            genome  = genome) # Genome of reference "hg38" or "hg19"
  }
}

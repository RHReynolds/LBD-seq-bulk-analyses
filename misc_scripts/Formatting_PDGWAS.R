# Description: Script to format PD and PD AOO GWAS for use with LDSC

#----Load libraries and data----
library(data.table)
library(LDSCforRyten)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(tidyverse)
library(stringr)

dbsnp_151 <- SNPlocs.Hsapiens.dbSNP151.GRCh38


#----PD GWAS------

print("Formatting PD GWAS")

GWAS <- fread("/data/LDScore/GWAS/PD2019_meta5_ex23andMe/nallsEtAl2019_excluding23andMe_allVariants.tab")

hg38 <- GWAS %>% 
  tidyr::separate(col = SNP, into = c("CHR", "BP"), sep = ":") %>% 
  LDSCforRyten::liftover_hg19_to_hg38(., path_to_chain = "/data/liftover/hg19/hg19ToHg38.over.chain") %>% 
  LDSCforRyten::add_RS_to_GWAS(., dbsnp_151)

fwrite(hg38, "/data/LDScore/GWAS/PD2019_meta5_ex23andMe/PD2019_ex23andMe_hg38.txt.gz")


#----PD AOO-------

print("Formatting PD AOO GWAS")

GWAS <- fread("/data/LDScore/GWAS/PD2018_AOO/sorted_AAO_april3_18_final_discovery.txt")

hg38 <- GWAS %>% 
  dplyr::rename(SNP = MarkerName) %>% 
  tidyr::separate(col = SNP, into = c("CHR", "BP"), sep = ":") %>% 
  LDSCforRyten::liftover_hg19_to_hg38(., path_to_chain = "/data/liftover/hg19/hg19ToHg38.over.chain") %>% 
  LDSCforRyten::add_RS_to_GWAS(., dbsnp_151)

fwrite(hg38, "/data/LDScore/GWAS/PD2018_AOO/PD2018_AOO_hg38.txt.gz")
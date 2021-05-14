# Description: Script to format progression GWAS for use with LDSC

#---Function----

format_PD_progression <- function(path_to_GWAS, ref_file_path, path_to_chain){
  
  library(data.table)
  library(LDSCforRyten)
  library(tidyverse)
  library(stringr)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  
  dbsnp <- SNPlocs.Hsapiens.dbSNP151.GRCh38
  
  ref <- fread(ref_file_path)
  
  GWAS_df <- data.frame(paths = list.files(path = path_to_GWAS, pattern = ".txt.gz", full.names = T),
                        base_or_surival = list.files(path = path_to_GWAS, pattern = ".txt.gz", full.names = F) %>% 
                          str_replace("_.*", ""),
                        phenotype = list.files(path = path_to_GWAS, pattern = ".txt.gz", full.names = F) %>% 
                          str_replace(".txt.gz", "") %>% 
                          str_replace(".*_", "")
  ) %>% 
    dplyr::filter(!phenotype == "reference")
  
  for(i in 1:nrow(GWAS_df)){
    
    print(str_c("Formatting: ", GWAS_df$base_or_surival[i], "_", GWAS_df$phenotype[i]))
    
    GWAS <- fread(GWAS_df$paths[i] %>%
                    as.character())
    
    GWAS <- GWAS %>%
      dplyr::inner_join(ref %>%
                          dplyr::select(SNP, CHR, START, REF, ALT, MAF)) %>%
      dplyr::rename(BP = START,
                    A1 = ALT,
                    A2 = REF) %>%
      dplyr::select(CHR, BP, A1, A2, MAF, BETA, SE, P, N, NSTUDY)
    
    hg38 <- LDSCforRyten::liftover_hg19_to_hg38(GWAS, path_to_chain) %>% 
      LDSCforRyten::add_RS_to_GWAS(., dbSNPref = dbsnp) %>% 
      dplyr::select(SNP, everything())
    
    fwrite(hg38, file = str_c("/data/LDScore/GWAS/PD2019_Progression/",
                              GWAS_df$base_or_surival[i], "_",
                              GWAS_df$phenotype[i], "_hg38.txt.gz"))
    
  }
  
}

#---Main----

format_PD_progression(path_to_GWAS = "/data/LDScore/GWAS/PD2019_Progression/",
                      ref_file_path = "/data/LDScore/GWAS/PD2019_Progression/reference.txt.gz",
                      path_to_chain = "/data/liftover/hg19/hg19ToHg38.over.chain")
# Description: Script to run analysis of motif enrichment

#---Load Libraries and data----------------------------------------------------------------------------
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(doParallel)
library(here)
library(foreach)
library(tidyverse)
library(rtracklayer)

# Set wd
here::i_am("misc_scripts/rbp_ame_analysis.R")

# Source fixed file paths
source(here::here("R", "file_paths.R"))

#---Functions------------------------------------------------------------------------------------------

source(here::here("R", "rbp_functions.R"))
source(here::here("R", "get_regions.R"))

#---Main-----------------------------------------------------------------------------------------------

print(str_c(Sys.time(), " - starting script..."))

# Load query list
query_list <- 
  readRDS(
    file.path(
      path_to_raw_data,
      "rbp_analysis/rbp_query_list.Rds"
    )
  )

# Filter to include only proximal intronic regions
query_list <- 
  query_list %>% 
  lapply(., function(list){
    
    list[str_detect(names(list), "prox")]
    
  })


print(str_c(Sys.time(), " - creating fasta files"))

# Export fasta files from query list
query_list %>% 
  lapply(., function(list){
    
    for(i in 1:length(list)){
      
      file_name <- str_c(unique(list[[i]]$comparison), ":", names(list)[i])
      
      gr <- list[[i]]
      
      query <- gr[c(abs(gr$deltapsi) >= 0.1 & gr$p.adjust < 0.05)]
      control <- gr[c(gr$p.adjust > 0.05)]
      
      generate_fasta_from_gr(gr = query,
                             genome = BSgenome.Hsapiens.UCSC.hg38.masked,
                             style = "UCSC",
                             remove_star_strand = T,
                             output_path = file.path(path_to_raw_data, "rbp_analysis/rbp_enrichment_fasta/query/"),
                             file_name = str_c(file_name, ":query"),
                             reduce = T)
      
      generate_fasta_from_gr(gr = control,
                             genome = BSgenome.Hsapiens.UCSC.hg38.masked,
                             style = "UCSC",
                             remove_star_strand = T,
                             output_path = file.path(path_to_raw_data, "rbp_analysis/rbp_enrichment_fasta/control/"),
                             file_name = str_c(file_name, ":control"),
                             reduce = T)
      
    }
    
  })

print(str_c(Sys.time(), " - running ame analysis"))

# Run MEME analysis
run_ame_analysis(path_to_ame ="/home/ssethi/meme/bin/ame",
                 output_dir = file.path(path_to_results, "rbp_analysis/rbp_enrichment_fasta/"),
                 path_to_query_fasta_dir = file.path(path_to_raw_data, "rbp_analysis/rbp_enrichment_fasta/"),
                 path_to_pwm_dir = file.path(path_to_raw_data, "rbp_analysis/attract_db"), 
                 query_name = "query",
                 control_name = "control",
                 cores = 9)


print(str_c(Sys.time(), " - done!"))

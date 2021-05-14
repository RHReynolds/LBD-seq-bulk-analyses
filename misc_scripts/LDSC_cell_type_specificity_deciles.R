# Description: Script to create cell-type specificity decile annotations for LDSC and run them

#---Load Libraries and data----------------------------------------------------------------------------

library(devtools)
library(doParallel)
library(foreach)
library(tidyverse)
library(qdapTools)
library(rtracklayer)
devtools::load_all(path = "/home/rreynolds/packages/LDSCforRyten/")

# Source fixed file paths
source(here::here("R", "file_paths.R"))

ctd_files <- 
  list.files(
    file.path(
      path_to_results,
      "snRNA/specificity_matrices/2020_Feb/"
    ), 
    pattern = "ctd", 
    full.names = T)
ctd_list <- vector(mode = "list", length = length(ctd_files))

# Loop to load data 
for(i in 1:length(ctd_files)){
  
  # Load ctd
  ctd_list[[i]] <- readRDS(ctd_files[i])
  
  # Extract file name
  title <- 
    ctd_files[i] %>% 
    str_replace(".*/", "") %>% 
    str_replace("\\..*", "") %>% 
    str_replace(".*_", "")
  
  # Name list
  names(ctd_list)[i] <- title
 
}

#---Creating annotations ------------------------------------------------------------------------------

# Create deciles dataframe
decile_df <- ctd_list %>% 
  lapply(., function(ctd){
    
    # Only using level 1 annotations
    ctd[[1]]$specificity %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "gene") %>% 
      dplyr::mutate(Vascular = Endo + Per) %>% 
      tidyr::gather(key = "cell_type", value = "specificity", -gene) %>% 
      dplyr::group_by(cell_type) %>% 
      dplyr::mutate(quantile = ntile(specificity, 10))
  }) %>% 
  qdapTools::list_df2df(col = "Disease_Group")

# As we have so many genes (many of which repeat), use ensembl version 97 .gtf to get gene start/end
ref <- import(path_to_ref_gtf)
# Keep only chr 1-22
ref <- ref %>% keepSeqlevels(c(1:22), pruning.mode = "coarse") 
ref <- ref[ref$type == "gene"]

# Extract positional information and add +/- 100kb window around gene start/end within reference
all_genes <-  ref %>% 
  as.data.frame() %>% 
  dplyr::select(seqnames, start, end, source, type, gene_id, gene_name, gene_biotype) %>% 
  # +/- 100 kb to start and end
  dplyr::mutate(start = ifelse(start - 100000 < 0, 0, start - 100000),
                end = end + 100000) %>% 
  # Keep only distinct gene names
  dplyr::distinct(gene_name, .keep_all = TRUE)

# Join positional info to decile dataframe
decile_df <- decile_df %>% 
  inner_join(all_genes %>% 
               dplyr::select(seqnames, start, end, gene_name), by = c("gene" = "gene_name")) %>% 
  dplyr::mutate(joint_name = str_c(Disease_Group, cell_type, quantile, sep = ":")) %>%
  dplyr::rename(start_pos = start,
                end_pos = end,
                CHR = seqnames) %>% 
  dplyr::arrange(joint_name) %>% 
  dplyr::select(joint_name, everything())

# Split dataframe by Disease_Group/cell_type/quantile into individual dataframes and store within a named list
ldsc_list <- setNames(decile_df %>% 
                        group_split(joint_name),
                      decile_df %>% 
                        .[["joint_name"]] %>% 
                        unique() %>% 
                        sort())

# # Loading baseline and creating Granges object from baseline model
# BM <- LDSCforRyten::creating_baseline_df(baseline_model = "53")
# BM_GR <- LDSCforRyten::df2GRanges(BM, 
#                                   seqname_col = "CHR", 
#                                   start_col = "BP",
#                                   end_col = "BP")


# As we have such a large list divide work across 15 iterations
subset_ids <- c(1:length(ldsc_list)) %>% 
  split(., sort(rep_len(1:15, length(.))))

# for(i in 1:length(subset_ids)){
for(i in 1:length(subset_ids)){
  
  print(str_c("Running loop: ", i, " in ", length(subset_ids)))
  
  # # Find overlapping regions between genes in each annotation and BM (baseline model)
  # print(Sys.time())
  # print("Find overlapping regions between genes in each annotation and BM...")
  # 
  ids <- subset_ids[[i]]
  # 
  # hits <- ldsc_list[ids] %>%
  #   LDSCforRyten::overlap_annot_list(query_GR = BM_GR, 
  #                                    seqname_col = "CHR", 
  #                                    start_col = "start_pos", 
  #                                    end_col = "end_pos",
  #                                    cores = 22)
  # 
  # # Annotate BM with 1s where annotation overlaps with the baseline model
  # print(Sys.time())
  # print("Annotate BM with 1s where annotation overlaps with the baseline model...")
  # 
  # list_BM <- hits %>%
  #   LDSCforRyten::overlap_annot_hits_w_baseline(BM = BM, 
  #                                               cores = 22)
  # 
  # # Exporting files
  # annot_dir <- file.path(path_to_raw_data, "ldsc_annotations/celltype.deciles/")
  # LDSCforRyten::create_annot_file_and_export(list_BM, annot_basedir = annot_dir)
  
  #---Running sLDSC -------------------------------------------------------------------------------------
  
  #---Defining arguments----
  annot_basedir <- file.path(path_to_raw_data, "ldsc_annotations/")
  annot_name <- "celltype.deciles"
  annot_subcategories <- names(ldsc_list[ids])
  baseline_model <- "53"
  gwas_df <- Create_GWAS_df() %>% 
    # dplyr::filter(str_detect(Original.name, "hg38") | Original.name == "AD2019") %>% 
    dplyr::filter(Original.name == "LBD2020.hg38")
  
  #---Fixed arguments---
  fixed_args <- get_LDSC_fixed_args(Baseline_model = baseline_model)
  
  #---Running LDSC---
  print(Sys.time())
  print("Start running LDSC...")
  
  # LDSCforRyten::Calculate_LDscore(Annotation_Basedir = annot_basedir, 
  #                                 Annot_name = annot_name, 
  #                                 Annotation_Subcategories = annot_subcategories, 
  #                                 Fixed_Arguments = fixed_args,
  #                                 cores = 2)
  
  LDSCforRyten::Calculate_H2(Annotation_Basedir = annot_basedir, 
                             Annot_name = annot_name, 
                             Annotation_Subcategories = annot_subcategories, 
                             Fixed_Arguments = fixed_args, 
                             GWAS_df = gwas_df,
                             cores = 8)
  
}

print(Sys.time())
print("Done!")


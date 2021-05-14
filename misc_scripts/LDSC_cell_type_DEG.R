# Description: Script to create cell-type DEG annotations for LDSC and run them

#---Load Libraries and data----------------------------------------------------------------------------

library(devtools)
library(doParallel)
library(foreach)
library(tidyverse)
library(qdapTools)
library(readxl)
library(rtracklayer)
devtools::load_all(path = "/home/rreynolds/packages/LDSCforRyten/")

# Source fixed file paths
source(here::here("R", "file_paths.R"))

DEG_files <- 
  list.files(
    file.path(
      path_to_results,
      "snRNA/differential_gene_expr/2021_Mar"
    ), 
    pattern = "_Allfiles2",
    full.names = T)

# Load across files
DEG_list <- setNames(DEG_files %>%
                       lapply(., function(file){
                         
                         file %>% 
                           readRDS() %>%
                           lapply(., function(df){
                             df %>%
                               dplyr::rename(gene = primerid)
                             
                           })
                         
                       }),
                     DEG_files %>% 
                       str_replace(".*/", "") %>% 
                       str_replace("_Allfiles2.*", ""))

#---Creating annotations ------------------------------------------------------------------------------

# Thresholds: FDR < 0.05 and |logFC| > log2(1.5)
fdr_threshold <- 0.05
fc_threshold <- log2(1.5)

# Create dataframe with FDR < 0.05 and |logFC| > log2(1.5)
DEG_df <- DEG_list %>% 
  lapply(., function(DEG){
    
    DEG %>% 
      qdapTools::list_df2df(col1 = "cell_type")
    
  }) %>% 
  qdapTools::list_df2df(col1 = "comparison") %>% 
  # Change naming convention to reflect our convention
  tidyr::separate(comparison, into = c("comparison_2", "comparison_1"), sep = "_") %>% 
  dplyr::mutate(comparison_1 = str_replace(comparison_1, pattern = "C", "Control"),
                comparison = str_c(comparison_1, ".vs.", comparison_2),
                # Add direction of effect
                direction_of_effect = ifelse(coef > 0, "up", "down"),
                cell_type = 
                  case_when(cell_type == "Astrocytes" ~ "Astro",
                            cell_type == "Oligodendrocytes" ~ "Oligo",
                            TRUE ~ cell_type)) %>% 
  dplyr::select(comparison, everything(), -comparison_1, -comparison_2) %>% 
  dplyr::filter(fdr < fdr_threshold, abs(coef) > fc_threshold) %>% 
  as_tibble()

DEG_df <- DEG_df %>% 
  dplyr::select(comparison, cell_type, gene, direction_of_effect)

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

# Join positional info to DEG dataframe and filter for lists (i.e. list with up + down DE genes) with more than 20 genes
DEG_df <- DEG_df %>% 
  left_join(all_genes %>% 
               dplyr::select(seqnames, start, end, gene_name), by = c("gene" = "gene_name")) %>% 
  dplyr::mutate(joint_name = str_c(comparison, cell_type, direction_of_effect, sep = ":")) %>%
  dplyr::rename(start_pos = start,
                end_pos = end,
                CHR = seqnames) %>% 
  dplyr::arrange(joint_name) %>% 
  dplyr::select(joint_name, everything()) %>% 
  dplyr::group_by(comparison, cell_type) %>% 
  dplyr::filter(n() >= 20) %>% 
  dplyr::ungroup()

# n_genes without positional information
print(str_c("Number of rows without positional information: ", 
            DEG_df %>% 
              dplyr::filter(is.na(start_pos)) %>% 
              nrow()))
print(str_c("Total rows: ", nrow(DEG_df)))

# Remove genes without positional information ~ 
DEG_df <-  DEG_df %>% 
  dplyr::filter(!is.na(start_pos))

# Split dataframe by into individual dataframes and store within a named list
# Run lists split by up/down-regulated genes as well as joint lists.
ldsc_list <- 
  c(
    # Creating joint lists
    setNames(DEG_df %>% 
               dplyr::mutate(joint_name = str_c(comparison, cell_type, sep = ":")) %>%
               group_split(joint_name),
             DEG_df %>% 
               dplyr::mutate(joint_name = str_c(comparison, cell_type, sep = ":")) %>%
               .[["joint_name"]] %>% 
               unique() %>% 
               sort()),
    setNames(DEG_df %>% 
               group_split(joint_name),
             DEG_df %>% 
               .[["joint_name"]] %>% 
               unique() %>% 
               sort())
    )

# Loading baseline and creating Granges object from baseline model
BM <- LDSCforRyten::creating_baseline_df(baseline_model = "53")
BM_GR <- LDSCforRyten::df2GRanges(BM,
                                  seqname_col = "CHR",
                                  start_col = "BP",
                                  end_col = "BP")


# As we have such a large list divide work across 10 iterations
subset_ids <- c(1:length(ldsc_list)) %>% 
  split(., sort(rep_len(1:10, length(.))))

for(i in 1:length(subset_ids)){

  print(str_c("Running loop: ", i, " in ", length(subset_ids)))
  
  # Find overlapping regions between genes in each annotation and BM (baseline model)
  print(Sys.time())
  print("Find overlapping regions between genes in each annotation and BM...")
  
  ids <- subset_ids[[i]]
  
  hits <- ldsc_list[ids] %>%
    LDSCforRyten::overlap_annot_list(query_GR = BM_GR,
                                     seqname_col = "CHR",
                                     start_col = "start_pos",
                                     end_col = "end_pos",
                                     cores = 10)

  # Annotate BM with 1s where annotation overlaps with the baseline model
  print(Sys.time())
  print("Annotate BM with 1s where annotation overlaps with the baseline model...")

  list_BM <- hits %>%
    LDSCforRyten::overlap_annot_hits_w_baseline(BM = BM,
                                                cores = 10)

  # Exporting files
  annot_dir <- file.path(path_to_raw_data, "ldsc_annotations/celltype.DEG/")
  LDSCforRyten::create_annot_file_and_export(list_BM, annot_basedir = annot_dir)
  
  #---Running sLDSC -------------------------------------------------------------------------------------
  
  #---Defining arguments----
  annot_basedir <- file.path(path_to_raw_data, "ldsc_annotations/")
  annot_name <- "celltype.DEG"
  annot_subcategories <- names(ldsc_list[ids])
  baseline_model <- "53"
  gwas_df <- Create_GWAS_df() %>% 
    dplyr::filter(str_detect(Original.name, "hg38") | Original.name == "AD2019")
    # dplyr::filter(Original.name == "LBD2020.hg38")
  
  #---Fixed arguments---
  fixed_args <- get_LDSC_fixed_args(Baseline_model = baseline_model)
  
  #---Running LDSC---
  print(Sys.time())
  print("Start running LDSC...")
  
  LDSCforRyten::Calculate_LDscore(Annotation_Basedir = annot_basedir,
                                  Annot_name = annot_name,
                                  Annotation_Subcategories = annot_subcategories,
                                  Fixed_Arguments = fixed_args,
                                  cores = 2)
  
  LDSCforRyten::Calculate_H2(Annotation_Basedir = annot_basedir, 
                             Annot_name = annot_name, 
                             Annotation_Subcategories = annot_subcategories, 
                             Fixed_Arguments = fixed_args, 
                             GWAS_df = gwas_df,
                             cores = 8)
  
}

print(Sys.time())
print("Done!")


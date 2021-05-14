# Description: Script to create new "_celltypes.txt" files for Scaden usage.

#---Load Libraries and data----------------------------------------------------------------------------

library(data.table)
library(readxl)
library(MASS)
library(qdapTools)
library(tidyverse)
library(stringr)

# Source fixed file paths
source(here::here("R", "file_paths.R"))

#---Notes----------------------------------------------------------------------------------------------

# 18/02/2019: new cell type labels provided, which effectively merged pericytes and endothelial cells into one cluster.

#---Main ----------------------------------------------------------------------------------------------

file_path <- file.path(path_to_raw_data, "deconvolution/scaden/2019Oct/")
output_path <- file.path(path_to_raw_data, "deconvolution/scaden/2020Feb/")

# Create dataframe of scaden-formatted files
sn_df <- data.frame(
  file_path = list.files(path = file_path, full.names = T, include.dirs = F),
  files = list.files(path = file_path, include.dirs = F),
  sample_name = list.files(path = file_path, include.dirs = F) %>% 
    str_replace("_.*", ""),
  file_type = list.files(path = file_path, include.dirs = F) %>% 
    str_replace(".txt", "") %>% 
    str_remove("^.*?_")) %>% 
  dplyr::filter(file_type %in% c("norm_counts_all", "celltypes"))
  
# Read in new cell type labels and format barcodes to suit norm_count_all files
ct <- read_delim(file = file.path(path_to_raw_data, "Cell.Type.Annotations.For.All.Samples.csv"), 
                 delim = " ", 
                 col_names = c("ID", "Celltype"), 
                 skip = 1) %>% 
  tidyr::separate(ID, into = c("sample_name", "disease_group", "barcode"), sep = "_") %>% 
  dplyr::mutate(ID = str_c(sample_name, "_", barcode)) 
  
# Extract sample names for loop
samples <- sn_df %>% .[["sample_name"]] %>% as.character() %>% unique()

# Creating new celltype files, with same order as count matrices
for (i in 1:length(samples)) {
  
  print(str_c("Loading '_celltypes.txt' file: ", samples[i]))
  
  # Read in count file
  norm_counts <- fread(sn_df %>% 
                         dplyr::filter(sample_name == samples[i], file_type == "norm_counts_all") %>% 
                         .[["file_path"]] %>% 
                         as.character()) %>% 
    as_tibble()
  
  # Filter for correct sample and merge endo and per into vascular
  celltype <- ct %>% 
    dplyr::filter(sample_name == samples[i]) %>% 
    dplyr::mutate(Celltype = case_when(Celltype %in% c("Endo", "Per") ~ "Vascular", 
                                       TRUE ~ Celltype))
  
  # Perform check that number of nuclei pre- and post-merge are the same
  merge_check <- c(ct %>% 
    dplyr::filter(sample_name == samples[i]) %>% 
    dplyr::filter(Celltype %in% c("Endo", "Per")) %>% nrow()) == 
    c(celltype %>% 
        dplyr::filter(Celltype == "Vascular") %>% 
        nrow())
  
  if(merge_check){
    print("Number of endo/per nuclei match the number of vascular nuclei after merge.")
  } else{
    stop("Number of endo/per did not match the number of vascular nuclei after merge.")
    }
  
  # Save cell type assignment in scaden format
  # Scaden needs _celltypes.txt (n_cells x 2)
  # Columns needed: Celltype and ID

  # Arrange order of IDs to match order of cell IDs in count data
  celltype <- celltype[match(norm_counts %>%
                               .[["cell_id"]],
                             celltype %>%
                               .[["ID"]]), ] %>% 
    dplyr::select(Celltype, ID)
  
  
  write_delim(celltype,
              path = str_c(output_path, "/",
                           samples[i],
                           "_celltypes.txt"), delim = "\t")
  
}

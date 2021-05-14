# Description: Script to transform snRNA-seq data into non-logarithmic space and to format data for Scaden usage.

#---Load Libraries and data----------------------------------------------------------------------------

library(readxl)
library(MASS)
library(qdapTools)
library(tidyverse)
library(stringr)

# Source fixed file paths
source(here::here("R", "file_paths.R"))

#---Notes----------------------------------------------------------------------------------------------

# 16/10/2019: new snRNA-seq count matrices provided due to change in filtering, thus following changes made:
# a. New file_path variable
# b. No library sizes provided, thus for now all operations re. library size have been commented out
# c. As there are no library sizes, cannot calculate raw counts
# d. Celltype item with new matrices no longer provided as a dataframe, thus commented out old code and made new

# 31/10/2019: library size file now provided, thus calculated raw counts. Only saved raw counts.

#---Main ----------------------------------------------------------------------------------------------

# file_path <- file.path(path_to_raw_data, "scRNAseq_counts")
# output_path <- file.path(path_to_raw_data, "deconvolution/scaden/")

file_path <- file.path(path_to_raw_data, "scRNAseq_counts_NEWfilteredmatrices")
output_path <- file.path(path_to_raw_data, "deconvolution/scaden/2019Oct/")

# Create dataframe of single nuclear files
sn_df <- data.frame(
  file_path = list.files(path = file_path, full.names = T),
  files = list.files(path = file_path),
  sample_name = list.files(path = file_path) %>% 
    str_replace("\\..*", ""),
  file_type = list.files(path = file_path) %>% 
    str_replace(".rds", "") %>% 
    str_replace(".*\\.", ""))

# Filter for only files paired with bulk data
bulk_samples <- read_excel(
  path = 
    file.path(path_to_raw_data,
              "sample_details/20190719_MasterFile_SampleInfo.xlsx"),
  sheet = "20190719_MasterFile_SampleInfo", 
  skip = 1
  ) %>% 
  dplyr::filter(Sample_Type == "Tissue section" & sent_to_bulk_seq == "yes") %>% 
  .[["CaseNo"]]

paired_df <- sn_df %>% 
  dplyr::filter(sample_name %in% bulk_samples)

# Extract sample names for loop
samples <- paired_df %>% .[["sample_name"]] %>% as.character() %>% unique()

# # Load library size file
library_sizes <- readRDS(file = str_c(file_path, "/PD_project.Sample.Library.sizes.rds"))

for (i in 1:length(samples)) {
  
  print(str_c("Transforming into non-logarithmic space: ", samples[i]))
  
  # Filter paired_df for relevant sample in loop
  sample_df <- paired_df %>% 
    dplyr::filter(sample_name == samples[i])
  
  # # Filter library sizes list for relevant sample
  sample_lib_size <- library_sizes[[samples[i]]]
  
  # Load sparse matrix
  count_matrix <- readRDS(file = c(sample_df %>% 
                                     dplyr::filter(file_type == "dgCMatrix") %>% 
                                     .[["file_path"]] %>% 
                                     as.character()))
  
  # Rename colnames so check can be performed to see they're in same order
  colnames(count_matrix) <- colnames(count_matrix) %>%
    str_replace("_.*_", "_")
  names(sample_lib_size) <- names(sample_lib_size) %>%
    str_replace("_.*_", "_")

  # Order count matrix columns by library size names to ensure same order when multiplying matrices
  count_matrix <- count_matrix[,match(names(sample_lib_size), colnames(count_matrix))]

  # Create a matrix of dimensions 1 x n_cells of the sample library sizes
  lib_matrix <- sample_lib_size %>%
    as.matrix() %>%
    t()
  
  # Repeat the 1st row in lib_matrix such that lib_matrix has the same dimensions as the count_matrix
  lib_matrix <- lib_matrix[rep(1, times = nrow(count_matrix)),]
  
  # Counts were originally log normalised using NormalizeData(..., normalization.method = "LogNormalize", scale.factor = 10000, ...) 
  # https://rdrr.io/cran/Seurat/man/NormalizeData.html
  # LogNormalize: Feature counts for each cell are divided by the total counts for that cell and 
  # multiplied by the scale.factor. This is then natural-log transformed using log1p, which computes
  # the natural logarithm of (value + 1). I.e. log counts (gene x, cell y) = ln(((gene_x_counts_cell_y/total_counts_cell_y) * 10000) + 1)
  # To transform back to raw counts need to perform:
  # raw_counts (gene x, cell y) = total_counts_cell_y * ((e^(log_counts) - 1)/10000)
  # Note that (e^(log_counts) - 1) = (raw_counts(gene x, cell y)/total_counts_cell_y) * 10000
  # In other words, it is the normalised count (normalised by library size) multipled by the scale factor
  norm_counts <- (exp(count_matrix)-1)
  raw_counts <- (norm_counts/10000) * lib_matrix
  
  print(str_c("Saving norm_counts and raw_counts: ", samples[i]))
  
  # Save both as dgCMatrix
  # saveRDS(norm_counts %>% as("dgCMatrix"), 
  #         file = str_c(file_path, "/norm_counts/", samples[i], "_norm_counts.dgCMatrix.rds"))
  
  saveRDS(raw_counts %>% as("dgCMatrix"),
          file = str_c(file_path, "/raw_counts/", samples[i], "_raw_counts.dgCMatrix.rds"))
  
  # Remove raw_counts from memory
  rm(raw_counts)
  
  print(str_c("Saving scaden-formatted data: ", samples[i]))

  # Save for use with Scaden
  norm_counts <- norm_counts %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "cell_id")

  write_delim(norm_counts,
              path = str_c(output_path, "/",
                           samples[i],
                           "_norm_counts_all.txt"), delim = "\t")

  # # Save cell type assignment in scaden format
  # # Scaden needs _celltypes.txt (n_cells x 2)
  # # Columns needed: Celltype and ID
  # celltype <- readRDS(file = c(sample_df %>%
  #                                dplyr::filter(file_type == "CellID_Type") %>%
  #                                .[["file_path"]] %>%
  #                                as.character())) %>%
  #   dplyr::mutate(Celltype = CellType,
  #                 ID = str_replace(ID, "_.*_", "_")) %>%
  #   dplyr::select(Celltype, ID)

  # New matrix celltype assignments provided as list of named vectors - turn into dataframe
  celltype <- readRDS(file = c(sample_df %>%
                                 dplyr::filter(file_type == "CellID_Type") %>%
                                 .[["file_path"]] %>%
                                 as.character())) %>%
    qdapTools::list2df(col1 = "Celltype", col2 = "ID") %>%
    dplyr::mutate(ID = str_replace(ID, "_.*_", "_")) %>%
    dplyr::select(Celltype, ID)


  # Arrange order of IDs to match order of cell IDs in count data
  celltype <- celltype[match(norm_counts %>%
                               .[["cell_id"]],
                             celltype %>%
                               .[["ID"]]), ]


  write_delim(celltype,
              path = str_c(output_path, "/",
                           samples[i],
                           "_celltypes.txt"), delim = "\t")
  
}

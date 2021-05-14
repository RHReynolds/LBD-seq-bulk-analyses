#' Calculate snRNA proportions
#'
#' Function to calculate snRNA cell-type proportions from cell-type assignment
#' files.
#'
#' @param path_to_snRNAseq_celltype_files File path to directory containing
#'   cell-type assignment files. These should have 2 columns named 'ID' and
#'   'Celltype'. (1) ID: cell/nuclei ID (2) Celltype: cell type.
#'
#' @return
#' @export
#' 

calculate_snRNA_proportions <- function(path_to_snRNAseq_celltype_files){
  
  library(tidyverse)
  library(stringr)
  
  # Create dataframe of single nuclear files
  sn_df <- data.frame(
    file_path = list.files(path = path_to_snRNAseq_celltype_files, full.names = T, include.dirs = F),
    files = list.files(path = path_to_snRNAseq_celltype_files, include.dirs = F),
    sample_name = list.files(path = path_to_snRNAseq_celltype_files, include.dirs = F) %>% 
      str_replace("_.*", ""),
    file_type = list.files(path = path_to_snRNAseq_celltype_files, include.dirs = F) %>% 
      str_replace(".*_", "") %>% 
      str_replace(".txt.", "")) %>% 
    dplyr::filter(sample_name != "bulkRNAseq")
  
  
  # Extract sample names for loop
  samples <- sn_df %>% .[["sample_name"]] %>% as.character() %>% unique()
  
  for (i in 1:length(samples)) {
    
    # Filter paired_df for relevant sample in loop
    sample_df <- sn_df %>% 
      dplyr::filter(sample_name == samples[i]) %>% 
      dplyr::filter(file_type == "celltypes.txt")
    
    celltype <- read_delim(file = sample_df$file_path[1] %>% as.character(), delim = "\t")
    
    celltype <- celltype %>% 
      dplyr::group_by(Celltype) %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::mutate(sample_id = samples[i],
                    snRNA_proportions = n/sum(n)) %>% 
      dplyr::select(-n) %>% 
      tidyr::spread(key = Celltype, value = snRNA_proportions)
    
    if(i == 1){
      
      master_df <- celltype
      
    } else{
      
      master_df <- master_df %>% 
        dplyr::bind_rows(celltype)
      
    }
    
    
  }
  
  return(master_df)
  
}
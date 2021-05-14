#' md5 checksum
#'
#' Function to check that md5s from original and transferred files match up.
#'
#' @param file_paths Path to transferred files that require an md5 check for
#'   corrupted file transfer.
#' @param original_md5 A dataframe with md5s from the files prior to transfer.
#'   Important that column with file names is named "file_name".
#' @param column_to_join_by Name of column in original_md5 dataframe containing
#'   the file names, which will be used to join original and transfer md5
#'   dataframes.
#'
#' @return Dataframe with column 'transfer_corrupted' with values TRUE/FALSE;
#'   FALSE if file is not corrupted, and TRUE if corrupted.
#' @export

md5_checksum <- function(file_paths, original_md5, column_to_join_by){
  
  library(tidyverse)
  
  # Create vectors of transferred md5s and file names.
  md5 <- md5sum(file_paths)  
  file_name <- names(md5) %>% 
    str_remove("/.*/")
  
  # Create dataframe
  transfer_md5 <- tibble::enframe(file_name, name = NULL, value = "file_name") %>% 
    bind_cols(tibble::enframe(md5, name = NULL, value = "transfer_md5"))
  
  # Create filtering vector for inner join
  filtering_vector <- "file_name"
  names(filtering_vector) <- column_to_join_by
  
  # Join dataframe
  check_df <- original_md5 %>% 
    inner_join(transfer_md5, by = filtering_vector) %>% 
    dplyr::mutate(transfer_corrupted = ifelse(original_md5 == transfer_md5, F, T))
  
  return(check_df)
  
}
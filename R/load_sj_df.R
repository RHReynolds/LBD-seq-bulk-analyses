#' Load splice junctions from first pass STAR alignment.
#'
#' @param sj_dir_path chr. Path to directory containing splice junctions files
#'   (SJ.out.tab).
#'
#' @return
#' @export
#' 

load_sj_df <- function(sj_dir_path){
  
  library(tidyverse)
  
  paths <- list.files(path = sj_dir_path, pattern = "SJ.out.tab", full.names = TRUE)
  
  for(i in 1:length(paths)){
    
    sample_name <- paths[i] %>% 
      str_replace("/.*/", "") %>%
      str_replace("_SJ.out.tab", "")
    
    cat("Loading splice junctions from:", paths[i],"\n")
    sj_df <- read_delim(file = paths[i], 
                        delim = "\t", 
                        col_names = c("chr", "intron_start", "intron_end", "strand", "intron_motif", "in_annotation", "unique_reads_junction", "multi_map_reads_junction", "max_splice_alignment_overhang"),
                        col_types = "cdddddddd") %>% 
      dplyr::mutate(Sample = sample_name)
    
    if(i == 1){
      master_df <- sj_df
    } else {
      master_df <- master_df %>% 
        bind_rows(sj_df)
    }
    
  }
  
  return(master_df)
  
}
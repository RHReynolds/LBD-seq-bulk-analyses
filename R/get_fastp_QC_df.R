#' Get fastp QC data.
#'
#' Allows users to extract a number of QC parameters outputted by fastp in a
#' json format.
#'
#' @param fastp_json_dir_paths Path in which the .json files are stored.
#'
#' @return Dataframe with QC parameters
#' @export

get_fastp_QC_df <- function(fastp_json_dir_paths){
  
  fastp_json_paths <- list.files(path = fastp_json_dir_paths, pattern = ".json", full.names = T)
  
  for (i in 1:length(fastp_json_paths)) {
    
    QC <- jsonlite::fromJSON(fastp_json_paths[i], flatten=TRUE)
    
    fastp_QC_df <- 
      tibble(SampleID = fastp_json_paths[i] %>% 
               str_replace("/.*/", "") %>% 
               str_replace("_fastp.json", ""),
             before_filtering_total_reads = QC$summary$before_filtering$total_reads,
             after_filtering_total_reads = QC$summary$after_filtering$total_reads,
             before_filtering_read1_mean_length = QC$summary$before_filtering$read1_mean_length,
             before_filtering_read2_mean_length = QC$summary$before_filtering$read2_mean_length,
             after_filtering_read1_mean_length = QC$summary$after_filtering$read1_mean_length,
             after_filtering_read2_mean_length = QC$summary$after_filtering$read2_mean_length,
             reads_unknown_insert_size_percent = (QC$insert_size$unknown/(sum(QC$insert_size$histogram)+QC$insert_size$unknown)) * 100)
    
    if (i == 1) {
      
      master_df <- fastp_QC_df
      
    } else{
      
      master_df <- master_df %>% 
        bind_rows(fastp_QC_df)
    }
    
  }
  
  return(master_df)
  
}

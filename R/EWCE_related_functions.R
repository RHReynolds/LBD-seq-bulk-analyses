#' Return labelled dendrogram from cell-type data outputted by \code{EWCE::generate.celltype.data()}.
#'
#' @param ctd Cell-type data in the format outputted by \code{EWCE::generate.celltype.data()}.
#' @param level integer. Level 1 or 2 annotations from ctd.
#' @param rotate logical. If TRUE dendrogram returned rotated. Default = FALSE.
#'
#' @return
#' @export
#'
plot_dendrogram <- function(ctd, level = c(1,2), rotate = FALSE){
  
  binned_file_dist <- ctd[[level]]$specificity_quantiles %>% 
    t() %>% 
    dist() # euclidean distances between the rows
  binned_file_dist_hclust <- binned_file_dist %>% hclust()
  plot <- ggdendro::ggdendrogram(binned_file_dist_hclust, 
                                 segments = TRUE,
                                 leaf_labels = TRUE,
                                 rotate = rotate, 
                                 type = "rectangle")
  
  return(plot)
  
}